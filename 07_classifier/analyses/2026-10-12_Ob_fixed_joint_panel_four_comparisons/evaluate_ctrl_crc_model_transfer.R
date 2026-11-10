#!/usr/bin/env Rscript

# Direct transfer of the originally locked Ob CTRL-vs-CRC classifiers.
# Marker membership, training subjects, transformations and tuning are fixed.
.libPaths(c('D:/R/library', .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(ranger); library(glmnet); library(e1071)
  library(xgboost); library(pROC); library(ggplot2)
})

ROOT <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-10-12_Ob_fixed_joint_panel_four_comparisons'
SOURCE <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-10-11_Ob_fixed_panel_xMarkerFinder_like'
STEM <- 'Ob22_CTRLCRC_model_transfer'
P <- c(
  dsv='H:/My Drive/proj_ob_AA_CRC/02.processed_data/05_SV_result/all_BMI_crossBMI/all_BMI_dsgv.csv',
  vsv='H:/My Drive/proj_ob_AA_CRC/02.processed_data/05_SV_result/all_BMI_crossBMI/all_BMI_vsgv.csv',
  meta='H:/My Drive/proj_ob_AA_CRC/02.processed_data/05_SV_result/all_BMI_crossBMI/metadata_allBMI_SGV_aligned.rds',
  abundance='H:/My Drive/proj_ob_AA_CRC/03.results/04_differential_abundance/heatmap_final/cache/abd_corrected_meta_clean.rds',
  feature_map='H:/My Drive/proj_ob_AA_CRC/analyses/2026-09-28_allBMI_BMIstratified_SV_landscape/tables/01_allBMI_full_SV_feature_to_SGB_map.tsv',
  callability='H:/My Drive/proj_ob_AA_CRC/analyses/2026-09-28_allBMI_BMIstratified_SV_landscape/tables/02_BMIstratified_SGB_callability.tsv',
  locked=file.path(SOURCE,'tables/restricted_ObAssoc_L1L2_locked_marker_panels.tsv'),
  split=file.path(SOURCE,'tables/restricted_ObAssoc_L1L2_fixed_development_test_split.tsv'),
  old_scores=file.path(SOURCE,'tables/restricted_ObAssoc_L1L2_algo_multi_algorithm_test_predictions.tsv'),
  tuning=file.path(SOURCE,'tables/restricted_ObAssoc_L1L2_algo_multi_algorithm_development_CV_tuning.tsv'),
  checkpoint=file.path(SOURCE,'checkpoints/Ob__restricted_ObAssoc_L1L2_joint.rds')
)
stopifnot(all(file.exists(P)))
for (d in c('figures','tables','data_for_figures'))
  dir.create(file.path(ROOT,d), showWarnings=FALSE, recursive=TRUE)

# Recreate the exact feature IDs used for the locked panel, reading only its SVs.
locked <- fread(P['locked'])[Panel=='joint'][order(Rank)]
selected <- readRDS(P['checkpoint'])$selected
stopifnot(nrow(locked)==22L, identical(locked$Feature,selected),
          sum(locked$Type=='species')==10L, sum(locked$Type=='SV')==12L)
sp_features <- locked[Type=='species',Feature]
sv_features <- locked[Type=='SV',Feature]
fm <- fread(P['feature_map'])
ca <- fread(P['callability'])
ca[,Retained:=as.character(Retained)%in%c('TRUE','True','true','1')]
sgb_union <- unique(ca[BMI_Group%in%c('Ob','Ow','Nw') & Retained,SGB])
fm <- fm[SGB%in%sgb_union & SV_type%in%c('dSGV','vSGV')]
fm[,Feature:=sprintf('%s%05d',ifelse(SV_type=='dSGV','DSV','VSV'),rowid(SV_type))]
stopifnot(all(sv_features%in%fm$Feature),
          all(locked[Type=='SV',Original_feature] ==
              fm[match(sv_features,Feature),SV_feature]))
read_sv <- function(path,type) {
  mp <- fm[SV_type==type & Feature%in%sv_features]
  z <- fread(path, select=c('V1',mp$SV_feature),
             check.names=FALSE, showProgress=FALSE)
  setnames(z,1L,'sample_id')
  old <- setdiff(names(z),'sample_id')
  setnames(z,old,mp[match(old,SV_feature),Feature])
  z
}
dd <- read_sv(P['dsv'],'dSGV')
vv <- read_sv(P['vsv'],'vSGV')
stopifnot(identical(dd$sample_id,vv$sample_id))
meta <- as.data.table(readRDS(P['meta']))
stopifnot(!anyDuplicated(meta$SGV_ID))
ii <- which(dd$sample_id%in%meta$SGV_ID)
stopifnot(length(ii)==nrow(meta))
dd <- dd[ii]; vv <- vv[ii]
sid <- dd$sample_id
dmat <- 1-as.matrix(dd[,-1]); vmat <- as.matrix(vv[,-1])
storage.mode(dmat) <- 'double'; storage.mode(vmat) <- 'double'
svmat <- cbind(dmat,vmat)
rownames(svmat) <- sid
meta <- meta[match(sid,SGV_ID)]
rm(dd,vv,dmat,vmat); invisible(gc())

ab <- readRDS(P['abundance'])
abd <- as.matrix(ab$abd_corrected)
storage.mode(abd) <- 'double'
abd[!is.finite(abd) | abd<0] <- 0
rownames(abd) <- sprintf('SP%05d',seq_len(nrow(abd)))
species_raw <- t(abd)
if (max(species_raw)>1.5) species_raw <- species_raw/100
ameta <- as.data.table(ab$meta_clean)
ameta[,species_id:=colnames(ab$abd_corrected)]
smap <- merge(meta[,.(SGV_ID,Study.name,Name)],
              ameta[,.(Study.name,Name,species_id)],
              by=c('Study.name','Name'),all.x=TRUE,sort=FALSE)
stopifnot(!anyDuplicated(smap$SGV_ID))
mi <- match(meta$SGV_ID,smap$SGV_ID)
ok <- !is.na(smap$species_id[mi])
meta <- meta[ok]
svmat <- svmat[ok,,drop=FALSE]
species_raw <- species_raw[smap$species_id[mi[ok]],,drop=FALSE]
rownames(species_raw) <- meta$SGV_ID
rm(ab,abd); invisible(gc())

keep <- meta$BMI_Group1=='Ob' & meta$Disease%in%c('CTRL','CRC','AA') &
        is.finite(meta$Age) & !is.na(meta$Sex) & meta$Sex!=''
meta <- meta[keep]
svmat <- svmat[keep,,drop=FALSE]
species_raw <- species_raw[keep,,drop=FALSE]
stopifnot(identical(rownames(svmat),meta$SGV_ID),
          identical(rownames(species_raw),meta$SGV_ID),
          all(sp_features%in%colnames(species_raw)),
          all(sv_features%in%colnames(svmat)))

split <- fread(P['split'])[order(row_id)]
stopifnot(nrow(split)==477L,
          identical(split$sample_id,meta[Disease%in%c('CTRL','CRC'),SGV_ID]),
          all(split$split%in%c('development','test')))
dev_ids <- split[split=='development',sample_id]
test_ids <- split[split=='test',sample_id]
aa_ids <- meta[Disease=='AA',SGV_ID]
score_ids <- c(test_ids,aa_ids)
stopifnot(length(intersect(dev_ids,score_ids))==0L,
          length(dev_ids)+length(test_ids)==477L,
          length(aa_ids)==195L,
          all(meta[match(test_ids,SGV_ID),Disease]==split[split=='test',Disease]))

# Original CLR pseudocount and SV imputation are learned only from CTRL/CRC development.
dev_sp <- species_raw[dev_ids,,drop=FALSE]
pos <- dev_sp[dev_sp>0]
pseudo <- if(length(pos)) min(pos)/2 else 1e-10
ids <- c(dev_ids,score_ids)
clr <- log(species_raw[ids,,drop=FALSE]+pseudo)
clr <- sweep(clr,1,rowMeans(clr),'-')
x <- matrix(NA_real_,nrow=length(ids),ncol=length(selected),
            dimnames=list(ids,selected))
x[,sp_features] <- clr[,sp_features,drop=FALSE]
x[,sv_features] <- svmat[ids,sv_features,drop=FALSE]
n_dev <- length(dev_ids)
for (j in seq_len(ncol(x))) {
  m <- median(x[seq_len(n_dev),j][is.finite(x[seq_len(n_dev),j])])
  stopifnot(is.finite(m))
  x[!is.finite(x[,j]),j] <- m
}
meta_aligned <- meta[match(ids,SGV_ID)]
x <- cbind(x,Age=meta_aligned$Age,
           SexMale=as.integer(meta_aligned$Sex=='Male'))
xtr <- x[seq_len(n_dev),,drop=FALSE]
xscore <- x[-seq_len(n_dev),,drop=FALSE]
ytr <- factor(meta_aligned$Disease[seq_len(n_dev)],levels=c('CTRL','CRC'))
wgt <- function(y) length(y)/(2*as.numeric(table(y)[as.character(y)]))

old_scores <- fread(P['old_scores'])[Panel=='joint']
tuning <- fread(P['tuning'])[Panel=='joint' & Selected==TRUE]
stopifnot(setequal(old_scores$sample_id,test_ids),
          uniqueN(tuning$Model)==4L, nrow(tuning)==4L)
methods <- c('Random forest','LASSO','Elastic Net','Radial SVM','XGBoost')
model_scores <- lapply(methods,function(method) {
  if (method=='Random forest') {
    set.seed(61004L)
    fit <- ranger(x=as.data.frame(xtr),y=ytr,probability=TRUE,
                  num.trees=500L,mtry=max(1L,floor(sqrt(ncol(xtr)))),
                  min.node.size=5L,case.weights=wgt(ytr),
                  importance='permutation',seed=61004L,num.threads=1L)
    p <- as.numeric(predict(fit,as.data.frame(xscore))$predictions[,'CRC'])
  } else {
    hp <- tuning[Model==method]
    if (method%in%c('LASSO','Elastic Net')) {
      fit <- glmnet(as.matrix(xtr),as.integer(ytr=='CRC'),family='binomial',
                    alpha=if(method=='LASSO') 1 else .5,
                    lambda=hp$lambda,weights=as.numeric(wgt(ytr)),standardize=TRUE)
      p <- as.numeric(predict(fit,newx=as.matrix(xscore),
                              type='response',s=hp$lambda))
    } else if (method=='Radial SVM') {
      set.seed(74001L)
      fit <- svm(x=as.data.frame(xtr),y=ytr,kernel='radial',
                 cost=hp$cost,gamma=hp$gamma_scale/ncol(xtr),
                 scale=TRUE,probability=TRUE,
                 class.weights=setNames(length(ytr)/(2*as.numeric(table(ytr))),levels(ytr)))
      pp <- predict(fit,as.data.frame(xscore),probability=TRUE)
      p <- as.numeric(attr(pp,'probabilities')[,'CRC'])
    } else {
      dtr <- xgb.DMatrix(as.matrix(xtr),label=as.integer(ytr=='CRC'),
                         weight=as.numeric(wgt(ytr)))
      set.seed(74002L)
      fit <- xgb.train(params=list(objective='binary:logistic',eval_metric='auc',
                   max_depth=hp$max_depth,eta=hp$eta,subsample=.8,
                   colsample_bytree=.8,min_child_weight=1,nthread=1,
                   tree_method='hist'),data=dtr,nrounds=hp$nrounds,verbose=0)
      p <- as.numeric(predict(fit,xgb.DMatrix(as.matrix(xscore))))
    }
  }
  stopifnot(length(p)==length(score_ids),all(is.finite(p)))
  data.table(Model=method,sample_id=score_ids,score=p)
})
pred <- rbindlist(model_scores)
pred <- merge(pred,meta[,.(sample_id=SGV_ID,Study=Study.name,Disease,
                           AJCC_stage=as.character(Tumor.Staging.AJCC))],
              by='sample_id',all.x=TRUE,sort=FALSE)
pred[,Source:=ifelse(sample_id%in%test_ids,'Original CTRL/CRC test','Untrained adenoma')]
stopifnot(nrow(pred)==length(methods)*length(score_ids),!anyNA(pred$Disease))

# Reject transfer estimates unless the refit reproduces every saved original test score.
check <- merge(pred[Source=='Original CTRL/CRC test',.(Model,sample_id,score)],
               old_scores[,.(Model,sample_id,old_score=score)],
               by=c('Model','sample_id'))
check[,Absolute_error:=abs(score-old_score)]
qc <- check[,.(N=.N,Max_absolute_error=max(Absolute_error),
                Mean_absolute_error=mean(Absolute_error)),by=Model]
stopifnot(nrow(check)==nrow(old_scores),all(qc$N==length(test_ids)),
          all(qc$Max_absolute_error<1e-6))

comp <- data.table(
  Comparison=c('CTRL vs CRC (original test)','Adenoma vs HC',
               'Adenoma vs CRC','Early CRC vs HC','Late CRC vs HC'),
  Control=c('CTRL','CTRL','AA','CTRL','CTRL'),
  Case=c('CRC','AA','CRC','CRC_I_II','CRC_III_IV'))
pred[,Group:=as.character(Disease)]
pred[,StageGroup:=fifelse(Disease=='CRC' & AJCC_stage%in%c('I','II'),
                         'CRC_I_II',fifelse(Disease=='CRC' &
                         AJCC_stage%in%c('III','IV'),'CRC_III_IV',
                         NA_character_))]
matches_group <- function(d,target) {
  if (target%in%c('CRC_I_II','CRC_III_IV'))
    !is.na(d$StageGroup) & d$StageGroup==target
  else d$Group==target
}
results <- rbindlist(lapply(seq_len(nrow(comp)),function(i) {
  a <- comp[i]
  rbindlist(lapply(methods,function(method) {
    d <- pred[Model==method]
    use <- (matches_group(d,a$Control) | matches_group(d,a$Case)) &
           (d$Disease=='AA' | d$Source=='Original CTRL/CRC test')
    z <- d[which(use)]
    z[,Outcome:=factor(ifelse(matches_group(z,a$Case),'Case','Control'),
                       levels=c('Control','Case'))]
    stopifnot(all(c('Case','Control')%in%as.character(z$Outcome)))
    ro <- pROC::roc(z$Outcome,z$score,
                     levels=c('Control','Case'),direction='<',quiet=TRUE)
    ci <- as.numeric(pROC::ci.auc(ro,method='delong'))
    data.table(Comparison=a$Comparison,Model=method,
      N_control=sum(z$Outcome=='Control'),N_case=sum(z$Outcome=='Case'),
      AUC=as.numeric(pROC::auc(ro)),CI_low=ci[1],CI_high=ci[3])
  }))
}))
expected_counts <- data.table(Comparison=comp$Comparison,
  N_control=c(47L,47L,195L,47L,47L),
  N_case=c(50L,195L,50L,18L,23L))
stopifnot(nrow(results)==25L,
          all(results[,.(Comparison,N_control,N_case)] ==
              expected_counts[match(results$Comparison,Comparison),
                              .(Comparison,N_control,N_case)]))

fwrite(pred,file.path(ROOT,'tables',paste0(STEM,'_predictions.tsv')),sep='\t')
fwrite(results,file.path(ROOT,'tables',paste0(STEM,'_AUC.tsv')),sep='\t')
fwrite(qc,file.path(ROOT,'tables',paste0(STEM,'_refit_QC.tsv')),sep='\t')

# XGBoost figure uses empirical ROC thresholds for the same frozen classifier.
roc_points <- function(z) {
  setorder(z,-score)
  n_case <- sum(z$Outcome=='Case'); n_control <- sum(z$Outcome=='Control')
  rr <- z[,.(score,FPR=cumsum(Outcome=='Control')/n_control,
             TPR=cumsum(Outcome=='Case')/n_case)]
  rbind(data.table(FPR=0,TPR=0),
        rr[!duplicated(score,fromLast=TRUE),.(FPR,TPR)])
}
figure_points <- rbindlist(lapply(seq_len(nrow(comp)),function(i) {
  a <- comp[i]
  d <- pred[Model=='XGBoost']
  use <- (matches_group(d,a$Control) | matches_group(d,a$Case)) &
         (d$Disease=='AA' | d$Source=='Original CTRL/CRC test')
  z <- copy(d[which(use)])
  z[,Outcome:=ifelse(matches_group(z,a$Case),'Case','Control')]
  q <- roc_points(z)
  q[,`:=`(Comparison=a$Comparison,Point_order=seq_len(.N),
          N_control=sum(z$Outcome=='Control'),N_case=sum(z$Outcome=='Case'))]
  q
}))
figure_points[,Comparison:=factor(Comparison,levels=comp$Comparison)]
labels <- results[Model=='XGBoost',.(Comparison,AUC,N_control,N_case)]
labels[,Comparison:=factor(Comparison,levels=comp$Comparison)]
labels[,`:=`(Label=sprintf('AUC = %.3f',AUC),Label_x=.97,Label_y=.08)]
plot_data <- rbindlist(list(
  figure_points[,.(Record_type='ROC point',Comparison,N_control,N_case,
                  FPR,TPR,Point_order,Label=NA_character_,
                  Label_x=NA_real_,Label_y=NA_real_)],
  labels[,.(Record_type='AUC annotation',Comparison,N_control,N_case,
             FPR=NA_real_,TPR=NA_real_,Point_order=NA_integer_,
             Label,Label_x,Label_y)]))
fwrite(plot_data,file.path(ROOT,'data_for_figures',paste0(STEM,'_XGBoost_ROC.tsv')),
       sep='\t',na='NA')
g <- ggplot(figure_points,aes(FPR,TPR)) +
  geom_abline(slope=1,intercept=0,linetype='dashed',
              colour='#BDBDBD',linewidth=.35) +
  geom_step(direction='hv',colour='#C44E52',linewidth=.9) +
  geom_text(data=labels,aes(x=Label_x,y=Label_y,label=Label),
            inherit.aes=FALSE,hjust=1,size=3.1,fontface='bold',
            colour='#C44E52') +
  facet_wrap(~Comparison,ncol=3,axes='all',axis.labels='all') +
  scale_x_continuous(breaks=c(0,.5,1),limits=c(0,1),expand=c(0,0)) +
  scale_y_continuous(breaks=c(0,.5,1),limits=c(0,1),expand=c(0,0)) +
  labs(title='Locked Ob CTRL vs CRC model: cross-stage transfer',
       x='False-positive rate',y='True-positive rate') +
  theme_classic(base_size=10,base_family='Arial') +
  theme(strip.background=element_blank(),strip.text=element_text(face='bold',size=10),
        plot.title=element_text(face='bold',size=12),
        axis.text=element_text(size=8,colour='#333333'),
        panel.spacing=grid::unit(.8,'cm'))
stem <- file.path(ROOT,'figures',paste0(STEM,'_XGBoost_ROC'))
grDevices::cairo_pdf(paste0(stem,'.pdf'),width=190/25.4,height=125/25.4,
                     family='Arial')
print(g); dev.off()
svglite::svglite(paste0(stem,'.svg'),width=190/25.4,height=125/25.4)
print(g); dev.off()
ragg::agg_png(paste0(stem,'.png'),width=190/25.4,height=125/25.4,
              units='in',res=300)
print(g); dev.off()
cat('TRANSFER_COMPLETE\n')
