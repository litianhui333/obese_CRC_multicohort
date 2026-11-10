from pathlib import Path
import subprocess
import csv
import hashlib
import json

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
META = ROOT / '02.processed_data/02_aligned_data/metadata_aligned_20260208.csv'
SOURCE_HASH = hashlib.sha256(META.read_bytes()).hexdigest()
R_CODE = r'''
args <- commandArgs(trailingOnly=TRUE)
m <- read.csv(args[1],check.names=FALSE,na.strings=c('NA','','NaN'))
out <- args[2]
stopifnot(nrow(m)==4385,!anyDuplicated(m$Sample.Source.ID.processed))
m$disease <- factor(m$Disease,levels=c('CTRL','AA','CRC'))
m$study <- factor(m$Study.name)
m$age <- m$Age
m$male <- ifelse(is.na(m$Sex),NA,as.integer(m$Sex=='Male'))
stopifnot(all(na.omit(m$Sex) %in% c('Male','Female')))
omni <- list(); pairs <- list(); coverage <- list(); diagnostics <- list()
append_row <- function(x, row) {x[[length(x)+1]] <- row; x}
robust_cov <- function(fit) {
 X <- model.matrix(fit)
 stopifnot(qr(X)$rank==ncol(X))
 h <- hatvalues(fit)
 stopifnot(all(1-h>1e-8))
 bread <- solve(crossprod(X))
 scores <- X * as.numeric(residuals(fit)/(1-h))
 bread %*% crossprod(scores) %*% bread
}
for(b in c('Normal','Ow','Ob')) {
 for(v in c('age','male')) {
  allrows <- droplevels(m[m$BMI_Group1==b & !is.na(m[[v]]),,drop=FALSE])
  study_split <- split(allrows,allrows$study,drop=TRUE)
  eligible <- vapply(study_split,function(d) length(unique(d$disease))>=2 && (v=='age' || length(unique(d$male))==2),logical(1))
  retained <- names(eligible)[eligible]
  dat <- droplevels(allrows[as.character(allrows$study) %in% retained,,drop=FALSE])
  stopifnot(all(table(dat$disease)>0))
  for(st in names(study_split)) {
   d <- study_split[[st]]
   coverage <- append_row(coverage,data.frame(BMI=b,Variable=v,Study=st,N=nrow(d),Disease_groups=paste(sort(unique(as.character(d$disease))),collapse='|'),
     Outcome_levels=length(unique(d[[v]])),Included=eligible[[st]],Reason=if(eligible[[st]]) 'retained' else if(length(unique(d$disease))<2) 'single_disease_group_no_within_cohort_contrast' else 'single_sex_no_within_cohort_outcome_variation'))
  }
  for(mode in c('pooled_all','unadjusted_same_samples','cohort_adjusted')) {
   d <- if(mode=='pooled_all') allrows else dat
   adjusted <- mode=='cohort_adjusted'
   f <- reformulate(c('disease',if(adjusted) 'study'),response=v)
   f0 <- reformulate(if(adjusted) 'study' else '1',response=v)
   warnings <- character()
   if(v=='age') {
    fit <- lm(f,data=d)
    V <- robust_cov(fit)
    ids <- grep('^disease',names(coef(fit)))
    beta <- coef(fit)[ids]
    stat <- as.numeric(t(beta)%*%solve(V[ids,ids,drop=FALSE],beta))
    p <- pchisq(stat,length(ids),lower.tail=FALSE)
    test <- 'OLS joint Wald chi-square; HC3 covariance'
    converged <- TRUE
   } else {
    fit <- withCallingHandlers(glm(f,data=d,family=binomial(),control=glm.control(maxit=100)),warning=function(w) {warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    fit0 <- withCallingHandlers(glm(f0,data=d,family=binomial(),control=glm.control(maxit=100)),warning=function(w) {warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    stopifnot(fit$converged,fit0$converged,all(is.finite(coef(fit))),fit$rank==length(coef(fit)),abs(coef(fit)['diseaseAA'])<15,abs(coef(fit)['diseaseCRC'])<15)
    stat <- deviance(fit0)-deviance(fit)
    ids <- grep('^disease',names(coef(fit)))
    p <- pchisq(stat,fit$rank-fit0$rank,lower.tail=FALSE)
    test <- 'Binomial logistic regression likelihood-ratio chi-square'
    V <- vcov(fit)
    converged <- TRUE
   }
   stopifnot(length(ids)==2,all(is.finite(V)))
   ns <- table(d$disease)
   omni <- append_row(omni,data.frame(BMI=b,Variable=v,Model=mode,N=nrow(d),Cohorts=nlevels(d$study),
     N_CTRL=as.integer(ns['CTRL']),N_AA=as.integer(ns['AA']),N_CRC=as.integer(ns['CRC']),
     N_available_all=nrow(allrows),Excluded_for_cohort_overlap=nrow(allrows)-nrow(dat),Test=test,Statistic=stat,DF=2,P=p))
   diagnostics <- append_row(diagnostics,data.frame(BMI=b,Variable=v,Model=mode,Converged=converged,
     Design_rank=fit$rank,Parameters=length(coef(fit)),Warnings=paste(unique(warnings),collapse=';'),
     Max_leverage=if(v=='age') max(hatvalues(fit)) else NA))
   if(mode=='pooled_all') next
   for(cmp in list(c('AA','CTRL'),c('CRC','CTRL'),c('CRC','AA'))) {
    L <- setNames(rep(0,length(coef(fit))),names(coef(fit)))
    if(cmp[1]!='CTRL') L[paste0('disease',cmp[1])] <- 1
    if(cmp[2]!='CTRL') L[paste0('disease',cmp[2])] <- -1
    estimate <- sum(L*coef(fit))
    se <- sqrt(as.numeric(t(L)%*%V%*%L))
    z <- estimate/se
    pp <- 2*pnorm(-abs(z))
    ci <- estimate+c(-1,1)*qnorm(.975)*se
    pairs <- append_row(pairs,data.frame(BMI=b,Variable=v,Model=mode,Group_1=cmp[1],Group_2=cmp[2],
      Effect=if(v=='age') 'mean_age_difference_years' else 'male_odds_ratio',
      Estimate=if(v=='age') estimate else exp(estimate),CI95_low=if(v=='age') ci[1] else exp(ci[1]),
      CI95_high=if(v=='age') ci[2] else exp(ci[2]),Z=z,P=pp,N=nrow(d)))
   }
  }
 }
}
omni <- do.call(rbind,omni); pairs <- do.call(rbind,pairs)
for(mode in unique(omni$Model)) {i <- omni$Model==mode;omni$BH_q[i] <- p.adjust(omni$P[i],'BH');omni$BH_family_N[i] <- sum(i)}
for(mode in unique(pairs$Model)) {i <- pairs$Model==mode;pairs$BH_q[i] <- p.adjust(pairs$P[i],'BH');pairs$BH_family_N[i] <- sum(i)}
stopifnot(nrow(omni)==18,nrow(pairs)==36,all(is.finite(omni$P)),all(is.finite(pairs$P)))
write_tsv <- function(x,name) write.table(x,file.path(out,'tables',name),sep='\t',row.names=FALSE,quote=FALSE,na='')
write_tsv(omni,'cohort_adjusted_omnibus.tsv')
write_tsv(pairs,'cohort_adjusted_pairwise.tsv')
write_tsv(do.call(rbind,coverage),'cohort_adjusted_cohort_eligibility.tsv')
write_tsv(do.call(rbind,diagnostics),'cohort_adjusted_model_diagnostics.tsv')
cat(R.version.string,'; cohort-adjusted analyses completed\n')
'''
result = subprocess.run([r'D:\Program Files\R-4.5.1\bin\x64\Rscript.exe','-',str(META),str(OUT)],input=R_CODE,text=True,capture_output=True)
if result.returncode:
    raise RuntimeError(result.stderr)
print(result.stdout.strip())

def read(name):
    with (OUT/'tables'/name).open(encoding='utf-8-sig',newline='') as handle:
        return list(csv.DictReader(handle,delimiter='\t'))
rows = read('cohort_adjusted_omnibus.tsv')
pairrows = read('cohort_adjusted_pairwise.tsv')
assert all(r['Converged']=='TRUE' and not r['Warnings'] for r in read('cohort_adjusted_model_diagnostics.tsv'))
assert hashlib.sha256(META.read_bytes()).hexdigest()==SOURCE_HASH
labels={'Normal':'Nw','Ow':'Ow','Ob':'Ob','age':'年龄','male':'性别（男性）'}
report = '''# BMI 层内疾病比较：队列调整分析

Status: exploratory

本次仅回答：在 Nw、Ow、Ob 各自内部，CTRL、adenoma（AA）、CRC 的年龄与性别是否存在差异，并调整 Study.name 队列信息。不以同一疾病内的 BMI 比较作为本次结论。

## 方法和解释

- 年龄：线性模型 Age ~ Disease + Study，队列固定效应，HC3 异方差稳健协方差；疾病总体效应使用 2 自由度 Wald 卡方检验。两两效应为调整队列后的平均年龄差（岁），并给出未做同时覆盖校正的 95% CI。这检验均值差异，与前次 Kruskal–Wallis 的秩/分布检验不是同一假设。
- 性别：二项 logistic 模型 Male ~ Disease + Study，队列固定效应；疾病总体效应使用完整与去除 Disease 模型的似然比检验。两两效应为男性 OR 和未做同时覆盖校正的 Wald 95% CI。OR>1 表示 Group_1 的男性 odds 高于 Group_2，不等于男性比例之比。
- 不额外互相调整年龄/性别，因为本次目的是检查这两个临床特征本身的组间分布；没有把研究问题改成年龄调整后的性别效应。
- 主分析仅保留该 BMI 层内有至少两种疾病、能够提供队列内疾病对比的队列；性别模型还要求同一队列中同时有男、女。单一疾病队列不能识别队列内疾病差异，单一性别队列不能为固定效应 logistic 的疾病 OR 提供可识别信息。排除原因逐队列列于 eligibility 表。三个疾病组通过共同队列网络比较；共同疾病效应假设未建模为疾病×队列交互。
- 输出 pooled_all（全部有效数据、不调整）、unadjusted_same_samples（与调整模型相同样本、不调整）和 cohort_adjusted（相同样本、调整队列）三种结果。判断队列调整的影响时应比较后两种，避免把样本集改变误当作调整效应。
- 每种模型的 6 项总体检验（3 BMI×2 特征）各自共同 BH 校正；相同样本未调整、队列调整两种模型的 18 项两两比较分别共同 BH 校正。不使用上一轮 14/42 项家族的 q 值替代本次 q 值。总体检验与两两比较是假设不同的检验，可能得出不同显著性结论。
- 每个变量按有效数据分析，不插补。保留原始样本、疾病和队列标签；描述表仍展示全部 4,385 例，模型的 N 另列。唯一处理后样本 ID 作为分析单位，不能据此排除同一参与者多次采样。
- 当前文件中肿瘤分期与部位仅 CRC 有数据，因此不能在 CTRL、AA、CRC 三种疾病之间作有意义的比较。腺瘤 AJCC 分期不适用，腺瘤部位信息亦未提供。

## 总体结果

| BMI | 特征 | 全数据未调整 P | 相同样本未调整 P | 队列调整 P | 队列调整 BH q | 调整模型 N（CTRL/AA/CRC） | 队列数 |
|---|---|---:|---:|---:|---:|---|---:|
'''
for bmi in ('Normal','Ow','Ob'):
    for var in ('age','male'):
        sel={r['Model']:r for r in rows if r['BMI']==bmi and r['Variable']==var}
        a=sel['cohort_adjusted']
        report += f"| {labels[bmi]} | {labels[var]} | {float(sel['pooled_all']['P']):.5g} | {float(sel['unadjusted_same_samples']['P']):.5g} | {float(a['P']):.5g} | {float(a['BH_q']):.5g} | {a['N']}（{a['N_CTRL']}/{a['N_AA']}/{a['N_CRC']}） | {a['Cohorts']} |\n"
report += '\n## 队列调整后的两两比较\n\n| BMI | 特征 | 比较（Group_1 vs Group_2） | 效应［95% CI］ | P | BH q |\n|---|---|---|---|---:|---:|\n'
for r in pairrows:
    if r['Model']!='cohort_adjusted': continue
    report += f"| {labels[r['BMI']]} | {labels[r['Variable']]} | {r['Group_1']} vs {r['Group_2']} | {float(r['Estimate']):.3f} [{float(r['CI95_low']):.3f}, {float(r['CI95_high']):.3f}] | {float(r['P']):.5g} | {float(r['BH_q']):.5g} |\n"
report += '\n以上为基线组间关联，不作因果解释，也不直接证明微生物结果受混杂。队列固定效应控制所纳入队列的截距差异，不能保证消除所有队列相关偏倚。年龄是均值差的稳健推断，性别是共同队列条件 OR；没有检验各队列效果异质性。没有重跑微生物关联模型。\n'
report += '\n复现：使用原项目记录的 bundled Python 执行本目录 run_cohort_adjusted.py；它仅依赖 Python 标准库并调用 base R 4.5.1。模型满秩、HC3 杠杆值、logistic 收敛及警告、样本数、输入文件 SHA256 均已检查。\n'
(OUT/'cohort_adjusted_report_zh.md').write_text(report,encoding='utf-8')
(OUT/'tables'/'cohort_adjusted_QA.json').write_text(json.dumps(dict(metadata_sha256=SOURCE_HASH,sources_unchanged=True,
    R=result.stdout.strip(),omnibus_rows=len(rows),pairwise_rows=len(pairrows),all_models_converged=True,model_warnings=[]),indent=2),encoding='utf-8')
for r in rows:
    if r['Model']=='cohort_adjusted': print(r)
