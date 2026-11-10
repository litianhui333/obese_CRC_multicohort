from pathlib import Path
import csv
from collections import Counter
import statistics

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
PROCESSED = ROOT / '02.processed_data/01_metadata_processed/metadata_processed_20260208_restored.csv'
ALIGNED = ROOT / '02.processed_data/02_aligned_data/metadata_aligned_20260208.csv'
SUPPLEMENT = Path('E:/proj_ob_AA_CRC/supplementary_tables/Supplementary_Table_1_metadata/Supplementary_Table_1_metadata_filtered.csv')

def read(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))

def write(name, rows):
    path = OUT / 'tables' / name
    with path.open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

def present(value):
    return str(value).strip().lower() not in ('', 'na', 'nan', 'none', 'null')

processed, aligned, supplement = map(read, (PROCESSED, ALIGNED, SUPPLEMENT))
idcol = 'Sample.Source.ID.processed'
keep = {r[idcol] for r in aligned}
excluded = [r for r in processed if r[idcol] not in keep]
assert len(processed) == 4550 and len(aligned) == 4385 and len(excluded) == 165
assert len({r[idcol] for r in aligned}) == len(aligned)
expected = {('Normal','CTRL'):997, ('Normal','AA'):404, ('Normal','CRC'):652,
            ('Ow','CTRL'):648, ('Ow','AA'):367, ('Ow','CRC'):511,
            ('Ob','CTRL'):310, ('Ob','AA'):226, ('Ob','CRC'):270}
counts = Counter((r['BMI_Group1'],r['Disease']) for r in aligned)
assert counts == expected
assert Counter((r['BMI_Group1'],r['Disease']) for r in supplement) == expected
# Name may be a processed ID in the supplement, so check both source Name and processed ID.
lookup = {(r['Study.name'], r['Name']) for r in aligned} | {(r['Study.name'], r[idcol]) for r in aligned}
assert all((r['Study.name'], r['Name']) in lookup for r in supplement)
write('excluded_from_aligned_165_records.tsv', excluded)
write('excluded_by_cohort_disease.tsv', [dict(Study_name=k[0], Disease=k[1], N=v)
    for k,v in sorted(Counter((r['Study.name'],r['Disease']) for r in excluded).items())])
summary, categorical, missing = [], [], []
for bmi in ('Normal','Ow','Ob'):
    for disease in ('CTRL','AA','CRC'):
        group = [r for r in aligned if r['BMI_Group1']==bmi and r['Disease']==disease]
        ages = sorted(float(r['Age']) for r in group if present(r['Age']))
        quartiles = statistics.quantiles(ages, n=4, method='inclusive')
        summary.append(dict(BMI_Group1=bmi,Disease=disease,N=len(group),Age_available=len(ages),
            Age_missing=len(group)-len(ages),Age_median=statistics.median(ages),Age_Q1=quartiles[0],
            Age_Q3=quartiles[2],Age_mean=statistics.mean(ages),Age_SD=statistics.stdev(ages)))
        for field in ('Age','Sex','Tumor.Staging.AJCC','Primary.Tumor.Location'):
            available = sum(present(r[field]) for r in group)
            applicable = field in ('Age','Sex') or disease=='CRC'
            missing.append(dict(BMI_Group1=bmi,Disease=disease,Variable=field,N=len(group),
                Applicable=applicable,Available=available,Missing_if_applicable=len(group)-available if applicable else 'not applicable'))
            if field=='Age':
                continue
            valid_n = available
            for level,n in sorted(Counter(r[field] if present(r[field]) else 'Missing' for r in group).items()):
                categorical.append(dict(BMI_Group1=bmi,Disease=disease,Variable=field,Level=level,
                    Applicable=applicable,N=n,Group_N=len(group),Available_N=valid_n,
                    Percent_of_group=100*n/len(group),Percent_of_available=100*n/valid_n if valid_n and level!='Missing' else ''))
write('nine_group_age_summary.tsv', summary)
write('nine_group_categorical_summary.tsv', categorical)
write('nine_group_missingness.tsv', missing)
write('sample_count_reconciliation.tsv', [dict(BMI_Group1=bmi,Disease=disease,
    Processed_N=sum(r['BMI_Group1']==bmi and r['Disease']==disease for r in processed),
    Aligned_N=n,Not_aligned_N=sum(r['BMI_Group1']==bmi and r['Disease']==disease for r in excluded))
    for (bmi,disease),n in expected.items()])
print('Verified: processed 4550; aligned 4385; excluded records 165; aligned IDs unique; supplement groups and membership match.')
print('Metadata reconciliation complete; baseline analyses follow. Source files unchanged.')

# Exploratory baseline analyses requested on 2026-11-02.
# Uses an existing Python installation with scipy/numpy; no environment installation.
import numpy as np
import scipy
from scipy import stats
from itertools import combinations
import json
import hashlib
import sys

BMI_LEVELS = ('Normal', 'Ow', 'Ob')
DISEASE_LEVELS = ('CTRL', 'AA', 'CRC')
FIELDS = ('Age', 'Sex', 'Tumor.Staging.AJCC', 'Primary.Tumor.Location')
DISPLAY = {'Normal':'Nw', 'Ow':'Ow', 'Ob':'Ob', 'CTRL':'CTRL', 'AA':'Adenoma', 'CRC':'CRC'}
LEVELS = {'Sex':('Female','Male'), 'Tumor.Staging.AJCC':('I','II','III','IV'),
          'Primary.Tumor.Location':('left_sided','right_sided')}
source_hashes = {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (PROCESSED,ALIGNED,SUPPLEMENT)}
for r in aligned:
    if present(r['Age']):
        assert np.isfinite(float(r['Age']))
    for field, levels in LEVELS.items():
        assert not present(r[field]) or r[field] in levels, (field,r[field])

omnibus, pairwise = [], []

def bh(rows):
    order = sorted(range(len(rows)), key=lambda i:rows[i]['P'])
    previous = 1.0
    for rank_index in range(len(order)-1,-1,-1):
        i = order[rank_index]
        previous = min(previous, rows[i]['P']*len(rows)/(rank_index+1))
        rows[i]['BH_q'] = previous
        rows[i]['BH_family_N'] = len(rows)

# Explicitly include six stratified comparison sets; no unstructured nine-group test.
scopes = [('within_BMI',b,'Disease',DISEASE_LEVELS,[r for r in aligned if r['BMI_Group1']==b]) for b in BMI_LEVELS]
scopes += [('within_disease',d,'BMI_Group1',BMI_LEVELS,[r for r in aligned if r['Disease']==d]) for d in DISEASE_LEVELS]
for scope, stratum, group_field, groups, rows in scopes:
    fields = list(FIELDS[:2])
    if scope=='within_disease' and stratum=='CRC':
        fields += list(FIELDS[2:])
    for field in fields:
        valid = [r for r in rows if present(r[field])]
        base = dict(Scope=scope,Stratum=stratum,Variable=field,Grouping=group_field,
                    Groups='|'.join(groups),N_total=len(rows),N_available=len(valid),N_missing=len(rows)-len(valid))
        group_n = {g:sum(r[group_field]==g for r in valid) for g in groups}
        base['N_by_group'] = json.dumps(group_n)
        if field=='Age':
            arrays = [np.array([float(r[field]) for r in valid if r[group_field]==g]) for g in groups]
            test = stats.kruskal(*arrays)
            effect = max(0.,(test.statistic-len(groups)+1)/(len(valid)-len(groups)))
            omnibus.append(dict(base,Test='Kruskal-Wallis (tie corrected)',Statistic=float(test.statistic),
                DF=len(groups)-1,P=float(test.pvalue),Effect_name='rank_epsilon_squared',Effect=effect,Min_expected=''))
        else:
            levels = LEVELS[field]
            observed = np.array([[sum(r[group_field]==g and r[field]==level for r in valid) for level in levels] for g in groups])
            chi,p,df,expected_cells = stats.chi2_contingency(observed,correction=False)
            # All observed analysis tables must meet this prespecified approximation check.
            assert expected_cells.min()>=5, (scope,stratum,field,expected_cells.min())
            effect = np.sqrt(chi/(observed.sum()*min(observed.shape[0]-1,observed.shape[1]-1)))
            omnibus.append(dict(base,Test='Pearson chi-square',Statistic=float(chi),DF=int(df),P=float(p),
                Effect_name='Cramers_V',Effect=float(effect),Min_expected=float(expected_cells.min())))
        for first,second in combinations(groups,2):
            pair = dict(Scope=scope,Stratum=stratum,Variable=field,Grouping=group_field,
                        Group_1=first,Group_2=second,N_1=group_n[first],N_2=group_n[second])
            if field=='Age':
                x,y = [np.array([float(r[field]) for r in valid if r[group_field]==g]) for g in (first,second)]
                result = stats.mannwhitneyu(x,y,alternative='two-sided',method='asymptotic',use_continuity=True)
                pairwise.append(dict(pair,Test='Mann-Whitney U (two-sided; tie corrected)',Statistic=float(result.statistic),
                    P=float(result.pvalue),Effect_name='rank_biserial_positive_if_group1_older',
                    Effect=float(2*result.statistic/(len(x)*len(y))-1),Median_difference_years=float(np.median(x)-np.median(y)),
                    Male_percentage_point_difference='',Min_expected=''))
            else:
                levels = LEVELS[field]
                observed = np.array([[sum(r[group_field]==g and r[field]==level for r in valid) for level in levels] for g in (first,second)])
                chi,p,df,expected_cells = stats.chi2_contingency(observed,correction=False)
                assert expected_cells.min()>=5
                effect = np.sqrt(chi/(observed.sum()*min(observed.shape[0]-1,observed.shape[1]-1)))
                difference = 100*(observed[0,1]/observed[0].sum()-observed[1,1]/observed[1].sum()) if field=='Sex' else ''
                pairwise.append(dict(pair,Test='Pearson chi-square',Statistic=float(chi),P=float(p),Effect_name='Cramers_V',
                    Effect=float(effect),Median_difference_years='',Male_percentage_point_difference=difference,
                    Min_expected=float(expected_cells.min())))
assert len(omnibus)==14 and len(pairwise)==42
bh(omnibus)
bh(pairwise)
write('baseline_omnibus_tests.tsv',omnibus)
write('baseline_pairwise_tests.tsv',pairwise)

# Readable nine-group baseline table. Missing percentages use total group N;
# nonmissing category percentages use field-specific available N.
wide = []
columns = [(b,d,DISPLAY[b]+'_'+DISPLAY[d]) for b in BMI_LEVELS for d in DISEASE_LEVELS]
for label,field,level in [('N','N',None),('Age, median [Q1, Q3]','Age',None),('Age missing, n (%)','Age','missing'),
    ('Female, n/N available (%)','Sex','Female'),('Male, n/N available (%)','Sex','Male'),('Sex missing, n (%)','Sex','missing')]+[
    ('AJCC '+level+', n/N available (%)','Tumor.Staging.AJCC',level) for level in LEVELS['Tumor.Staging.AJCC']]+[
    ('AJCC missing, n (%)','Tumor.Staging.AJCC','missing')]+[
    ('Location '+level+', n/N available (%)','Primary.Tumor.Location',level) for level in LEVELS['Primary.Tumor.Location']]+[
    ('Location missing, n (%)','Primary.Tumor.Location','missing')]:
    row = {'Characteristic':label}
    for bmi,disease,col in columns:
        group = [r for r in aligned if r['BMI_Group1']==bmi and r['Disease']==disease]
        if field=='N':
            value = str(len(group))
        elif field in FIELDS[2:] and disease!='CRC':
            value = 'Not applicable'
        else:
            available = [r for r in group if present(r[field])]
            if level=='missing':
                n = len(group)-len(available)
                value = f'{n} ({100*n/len(group):.1f}%)'
            elif field=='Age':
                ages = [float(r[field]) for r in available]
                q1,median,q3 = np.quantile(ages,[.25,.5,.75])
                value = f'{median:.1f} [{q1:.1f}, {q3:.1f}]'
            else:
                n = sum(r[field]==level for r in available)
                value = f'{n}/{len(available)} ({100*n/len(available):.1f}%)'
        row[col] = value
    wide.append(row)
write('baseline_nine_groups_publication_table.tsv',wide)
# Cohort availability and composition support interpretation of pooled tests.
cohort_rows = []
for study in sorted({r['Study.name'] for r in aligned}):
    for bmi,disease,col in columns:
        group = [r for r in aligned if r['Study.name']==study and r['BMI_Group1']==bmi and r['Disease']==disease]
        if not group:
            continue
        cohort_rows.append(dict(Study_name=study,BMI_Group1=bmi,Disease=disease,N=len(group),
            **{field+'_available':sum(present(r[field]) for r in group) for field in FIELDS}))
write('baseline_cohort_composition_and_availability.tsv',cohort_rows)

notes = '''# Nine-group baseline analysis (including adenoma)

Status: exploratory

Analysis population: 4,385 aligned records with unique processed sample IDs, from 21 cohorts. Each record is treated as an independent sample. Participant-level repeated-sample independence cannot be established beyond these metadata IDs.

All nine disease-by-BMI groups include adenoma (AA). Age is reported as median [Q1, Q3]. Categorical percentages use available data for the specific characteristic; missingness percentages use total group N. Tumor AJCC stage and location apply to CRC only in this dataset. No raw values or labels were renamed; display labels use Nw and Adenoma.

The 14 omnibus tests comprise age and sex across diseases within each BMI stratum (6), age and sex across BMI strata within each disease including adenoma (6), and CRC stage and location across BMI strata (2). Kruskal-Wallis tests use tie correction; Pearson chi-square tests use no continuity correction. All categorical expected-cell counts are checked to be at least five. BH correction covers all 14 omnibus tests together. All 42 pairwise contrasts are reported in a separate, secondary BH family (42 tests), irrespective of omnibus significance; they are not selectively chosen. Mann-Whitney tests are two-sided, asymptotic with tie and continuity correction.

Omnibus effect sizes are rank epsilon-squared for age (bounded below by zero) and Cramer's V for categorical characteristics. Pairwise age effects use signed rank-biserial correlation, positive when group 1 tends older; differences in age medians are descriptive, not the specific null tested by Mann-Whitney. Sex pairwise tables include group-1 minus group-2 male percentage points. Kruskal-Wallis tests distribution/rank differences, not necessarily medians alone. No confidence intervals are supplied for these exploratory rank/association effect sizes.

These are pooled, unadjusted descriptive comparisons. Differences may reflect cohort composition, region, sampling or missingness. They do not establish an independent BMI effect or demonstrate confounding of microbiome findings. Non-significance does not establish equivalence. No microbiome models or covariate-adjusted sensitivity analyses were rerun. Stage/location missing values are excluded variable-wise; no imputation was performed. The cohort availability table helps identify uneven coverage.

'''
notes += '## Omnibus results\n\n| Scope | Stratum | Variable | Available N | P | BH q | Effect |\n|---|---|---|---:|---:|---:|---:|\n'
for r in omnibus:
    notes += f"| {r['Scope']} | {DISPLAY.get(r['Stratum'],r['Stratum'])} | {r['Variable']} | {r['N_available']} | {r['P']:.6g} | {r['BH_q']:.6g} | {r['Effect']:.4f} |\n"
notes += '\n## Descriptive table\n\n'
notes += '| '+' | '.join(wide[0])+' |\n|'+ '|'.join(['---']*len(wide[0]))+'|\n'
for row in wide:
    notes += '| '+' | '.join(row.values())+' |\n'
notes += '\n## Reproduction and checks\n\nExisting Python: '+sys.executable+'; Python '+sys.version.split()[0]+', NumPy '+np.__version__+', SciPy '+scipy.__version__+'.\n'
notes += 'Source hashes, nine-group counts, unique aligned IDs, clinical category levels, expected-cell checks, 14/42 test counts and BH adjustment were verified. Sources are read-only.\n'
(OUT/'baseline_analysis_report.md').write_text(notes,encoding='utf-8')
assert source_hashes == {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (PROCESSED,ALIGNED,SUPPLEMENT)}
(OUT/'tables'/'baseline_run_QA.json').write_text(json.dumps(dict(source_hashes=source_hashes,
    Python=sys.version,Numpy=np.__version__,Scipy=scipy.__version__,N=4385,
    omnibus_tests=len(omnibus),pairwise_tests=len(pairwise),sources_unchanged=True,
    minimum_categorical_expected_count=min(float(r['Min_expected']) for r in pairwise if r['Min_expected']!='')),
    ensure_ascii=False,indent=2),encoding='utf-8')
print('Baseline complete: 14 omnibus tests, 42 pairwise contrasts, adenoma included; BH families 14 and 42.')
for r in omnibus:
    print(r['Scope'],r['Stratum'],r['Variable'],f"P={r['P']:.6g}",f"q={r['BH_q']:.6g}")


# Independent verification in base R, including every P value and both BH families.
# Required because this installed SciPy warns about its NumPy support range.
import subprocess
R_SCRIPT = r'''
args <- commandArgs(trailingOnly=TRUE)
m <- read.csv(args[1], check.names=FALSE, na.strings=c('NA','','NaN'))
o <- read.delim(file.path(args[2],'tables','baseline_omnibus_tests.tsv'),check.names=FALSE,fileEncoding='UTF-8-BOM')
p <- read.delim(file.path(args[2],'tables','baseline_pairwise_tests.tsv'),check.names=FALSE,fileEncoding='UTF-8-BOM')
get_rows <- function(row) {
  if(row$Scope=='within_BMI') m[m$BMI_Group1==row$Stratum,,drop=FALSE] else m[m$Disease==row$Stratum,,drop=FALSE]
}
validate <- function(tab, pair=FALSE) {
  rp <- numeric(nrow(tab))
  for(i in seq_len(nrow(tab))) {
    row <- tab[i,]; dat <- get_rows(row)
    if(pair) dat <- dat[dat[[row$Grouping]] %in% c(row$Group_1,row$Group_2),,drop=FALSE]
    dat <- dat[!is.na(dat[[row$Variable]]),,drop=FALSE]
    g <- factor(dat[[row$Grouping]])
    if(row$Variable=='Age') {
      if(pair) rp[i] <- wilcox.test(dat[[row$Variable]]~g,exact=FALSE,correct=TRUE)$p.value else rp[i] <- kruskal.test(dat[[row$Variable]]~g)$p.value
    } else rp[i] <- chisq.test(table(g,dat[[row$Variable]]),correct=FALSE)$p.value
  }
  rq <- p.adjust(rp,method='BH')
  stopifnot(all(abs(rp-tab$P)<=1e-12+1e-7*abs(rp)),all(abs(rq-tab$BH_q)<=1e-12+1e-7*abs(rq)))
  data.frame(Family=if(pair) 'pairwise42' else 'omnibus14',Scope=tab$Scope,Stratum=tab$Stratum,
    Variable=tab$Variable,Python_P=tab$P,R_P=rp,Absolute_P_difference=abs(rp-tab$P),Python_BH_q=tab$BH_q,R_BH_q=rq)
}
checks <- rbind(validate(o),validate(p,TRUE))
write.table(checks,file.path(args[2],'tables','baseline_base_R_validation.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cat(R.version.string,'; independently verified all 56 P values and BH q values\n')
'''
check = subprocess.run([r'D:\Program Files\R-4.5.1\bin\x64\Rscript.exe','-',str(ALIGNED),str(OUT)],
    input=R_SCRIPT,text=True,capture_output=True,check=False)
print(check.stdout.strip())
if check.returncode:
    raise RuntimeError(check.stderr)
qa_path = OUT/'tables'/'baseline_run_QA.json'
qa = json.loads(qa_path.read_text(encoding='utf-8'))
qa['independent_base_R_verification'] = check.stdout.strip()
qa['runtime_note'] = 'Installed SciPy warns NumPy 2.5.3 is outside supported range; all 56 P values and BH values independently verified in base R.'
qa_path.write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
with (OUT/'baseline_analysis_report.md').open('a',encoding='utf-8') as handle:
    handle.write('\nIndependent verification: '+check.stdout.strip()+'. The local SciPy/NumPy compatibility warning was addressed by independently recalculating every omnibus and pairwise P value and BH q value in base R (absolute tolerance 1e-12 plus relative tolerance 1e-7); all passed.\n')



zh = '''# 九组临床基线分析（含 adenoma）

Status: exploratory

使用与正文一致的 metadata_aligned_20260208.csv，共 4,385 条、21 个队列；其中 CTRL 1,955、adenoma（AA）997、CRC 1,433。保留九组原始人数和样本标签。

## 主要结果

- AA 在三个 BMI 组之间，年龄总体检验 BH q=0.00243，性别总体检验 BH q=0.000396。
- AA 年龄中位数［IQR］：Nw 63［57–68］岁；Ow 64［58–68］岁；Ob 65［61–69］岁。
- AA 已知性别者中的男性比例：Nw 113/358（31.6%）；Ow 137/334（41.0%）；Ob 54/219（24.7%）。
- AA 两两比较（42 项两两检验共同 BH 校正）：年龄 Nw–Ob q=0.00160、Ow–Ob q=0.0135、Nw–Ow q=0.457；性别 Nw–Ow q=0.0214、Ow–Ob q=0.000265、Nw–Ob q=0.139。
- CTRL 在三个 BMI 组之间的年龄、性别分布也有差异；CRC 内的年龄、性别、AJCC 分期、部位均未检出差异。未检出差异不代表等同。
- 在每个 BMI 层内比较 CTRL/AA/CRC：Nw 和 Ow 的年龄分布有差异，Ob 年龄未检出差异；三个 BMI 层的性别分布均有差异。

## 方法与边界

年龄报告中位数［Q1, Q3］，总体检验使用含并列值校正的 Kruskal–Wallis；两两比较使用双侧 Mann–Whitney U，含并列值和连续性校正。分类变量使用 Pearson 卡方检验；全部列联表最小期望频数不少于 5。14 项总体检验共同 BH 校正；42 项两两比较作为单独的探索性检验家族共同 BH 校正，完整报告全部结果。

分类百分比按各变量非缺失人数计算；缺失率按各组总人数计算，不插补。AA 的年龄和性别缺失人数分别均为 Nw 46、Ow 33、Ob 7。AJCC 分期不适用于腺瘤；当前文件也没有可用于分析的腺瘤部位信息，因此本次肿瘤部位、分期仅比较 CRC。CRC 分期缺失 266/1433（18.6%），部位缺失 404/1433（28.2%）。

以上为合并队列、未经协变量调整的探索性分布比较。年龄总体检验并非只检验中位数；两两年龄中位数差仅为描述值。独立单位暂按唯一处理后样本 ID 计；现有 metadata 不能进一步验证是否存在同一参与者多次采样。结果可能受队列构成、地区及信息缺失影响，不能据此认定独立 BMI 效应或微生物结果已被混杂。本次不重跑微生物关联模型。建议结合研究设计，在相关微生物模型中考虑年龄、性别、队列/地区；CRC 部位与分期适合作为 CRC 限定的敏感性分析，需考虑缺失覆盖。

## 输出与验证

- tables/baseline_nine_groups_publication_table.tsv：可用 Excel 打开的九组描述表，含 AA、临床数据分母和缺失情况。
- tables/baseline_omnibus_tests.tsv：14 项总体检验，含原始 P、BH q、效应量、有效样本数。
- tables/baseline_pairwise_tests.tsv：42 项两两比较，含 AA。
- tables/baseline_cohort_composition_and_availability.tsv：按队列、BMI 和疾病的组成与临床信息覆盖。
- baseline_analysis_report.md：完整方法、全部总体结果及九组宽表。

Python/SciPy 产生本地 NumPy 兼容范围提示；已使用 base R 4.5.1 独立复算全部 56 个 P 值及两组 BH 校正值，结果一致。输入文件 SHA256 未改变，核验见 tables/baseline_run_QA.json 与 tables/baseline_base_R_validation.tsv。未生成科研图，因此 figures/ 和 data_for_figures/ 保持为空。
'''
(OUT/'baseline_summary_zh.md').write_text(zh,encoding='utf-8')
