#!/usr/bin/env Rscript
.libPaths(c('D:/R/library', .libPaths()))
suppressPackageStartupMessages({library(data.table);library(ranger);library(pROC);library(ggplot2)})

# Exploratory cross-BMI transfer sensitivity analysis. Marker panels and tuning
# parameters are locked from the prior Ob development workflow; Nw/Ow are never
# used for feature selection, tuning, preprocessing, or fitting.
ROOT <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-10-11_Ob_fixed_panel_xMarkerFinder_like'
PREFIX <- 'restricted_ObAssoc_L1L2_crossBMI_otherAlgorithms_'
PATHS <- c(
 dsv='H:/My Drive/proj_ob_AA_CRC/02.processed_data/05_SV_result/all_BMI_crossBMI/all_BMI_dsgv.csv',
 vsv='H:/My Drive/proj_ob_AA_CRC/02.processed_data/05_SV_result/all_BMI_crossBMI/all_BMI_vsgv.csv',
 sv_meta='H:/My Drive/proj_ob_AA_CRC/02.processed_data/05_SV_result/all_BMI_crossBMI/metadata_allBMI_SGV_aligned.rds',
 abundance='H:/My Drive/proj_ob_AA_CRC/03.results/04_differential_abundance/heatmap_final/cache/abd_corrected_meta_clean.rds',
 feature_map='H:/My Drive/proj_ob_AA_CRC/analyses/2026-09-28_allBMI_BMIstratified_SV_landscape/tables/01_allBMI_full_SV_feature_to_SGB_map.tsv',
 callability='H:/My Drive/proj_ob_AA_CRC/analyses/2026-09-28_allBMI_BMIstratified_SV_landscape/tables/02_BMIstratified_SGB_callability.tsv',
 markers=file.path(ROOT,'tables','restricted_ObAssoc_L1L2_locked_marker_panels.tsv'),
 split=file.path(ROOT,'tables','restricted_ObAssoc_L1L2_fixed_development_test_split.tsv'),
 tuning=file.path(ROOT,'tables','restricted_ObAssoc_L1L2_algo_multi_algorithm_development_CV_tuning.tsv'))
stopifnot(all(file.exists(PATHS)))
dir.create(file.path(ROOT,'tables'),showWarnings=FALSE,recursive=TRUE)
dir.create(file.path(ROOT,'figures'),showWarnings=FALSE,recursive=TRUE)
dir.create(file.path(ROOT,'data_for_figures'),showWarnings=FALSE,recursive=TRUE)

PANELS <- c(species='Species only',sv='SV only',joint='Species + SV')
PANEL_ORDER <- names(PANELS)
MODELS <- c('Random forest','LASSO','Elastic Net','Radial SVM','XGBoost')
BMI_ORDER <- c('Ob_test','Nw','Ow')
N_REPEATS <- 5L
SPLIT_SEEDS <- 78101L + seq_len(N_REPEATS)
FIT_SEED_BASE <- 79100L
N_TREES <- 500L
wgt <- function(y) length(y)/(2*as.numeric(table(y)[as.character(y)]))
write_tsv <- function(x,name,where='tables') fwrite(x,file.path(ROOT,where,paste0(PREFIX,name,'.tsv')),sep='\t',na='NA')

markers <- fread(PATHS['markers'])
split <- fread(PATHS['split'])
tuning <- fread(PATHS['tuning'])
stopifnot(setequal(unique(markers$Panel),PANEL_ORDER),
          all(c('Feature','Panel','Rank')%in%names(markers)),
          all(c('Panel','Model','Selected')%in%names(tuning)))
selected_tuning <- tuning[Selected==TRUE & Model!='Random forest']
stopifnot(nrow(selected_tuning)==length(PANEL_ORDER)*4L,
          uniqueN(selected_tuning[,.(Panel,Model)])==length(PANEL_ORDER)*4L)

message('Loading matched SV and species data...')
fm <- fread(PATHS['feature_map']); ca <- fread(PATHS['callability'])
ca[,Retained:=as.character(Retained)%in%c('TRUE','True','true','1')]
sgb_union <- unique(ca[BMI_Group%in%c('Ob','Ow','Nw')&Retained==TRUE,SGB])
fm <- fm[SGB%in%sgb_union&SV_type%in%c('dSGV','vSGV')]
fm[,Feature:=sprintf('%s%05d',ifelse(SV_type=='dSGV','DSV','VSV'),rowid(SV_type))]
sv_need <- unique(markers[Type=='SV',Feature])
read_sv_selected <- function(path,type) {
 mp <- fm[SV_type==type&Feature%in%sv_need]
 z <- fread(path,select=c('V1',mp$SV_feature),check.names=FALSE,showProgress=FALSE)
 setnames(z,1L,'sample_id'); old <- setdiff(names(z),'sample_id')
 setnames(z,old,mp[match(old,SV_feature),Feature]); z
}
dd <- read_sv_selected(PATHS['dsv'],'dSGV'); vv <- read_sv_selected(PATHS['vsv'],'vSGV')
stopifnot(identical(dd$sample_id,vv$sample_id))
meta <- as.data.table(readRDS(PATHS['sv_meta']))
ii <- match(meta$SGV_ID,dd$sample_id); stopifnot(!anyNA(ii))
dd <- dd[ii]; vv <- vv[ii]
dmat <- 1-as.matrix(dd[,-1]); storage.mode(dmat)<-'double'
vmat <- as.matrix(vv[,-1]); storage.mode(vmat)<-'double'
svmat <- cbind(dmat,vmat); rownames(svmat)<-meta$SGV_ID
rm(dd,vv,dmat,vmat); invisible(gc())

aa <- readRDS(PATHS['abundance']); abd <- as.matrix(aa$abd_corrected)
storage.mode(abd)<-'double'; abd[!is.finite(abd)|abd<0]<-0
spmap <- data.table(Feature=sprintf('SP%05d',seq_len(nrow(abd))),Original_feature=rownames(abd))
rownames(abd)<-spmap$Feature; species_raw<-t(abd)
if(max(species_raw)>1.5) species_raw<-species_raw/100
ameta<-as.data.table(aa$meta_clean);ameta[,species_id:=colnames(aa$abd_corrected)]
smap<-merge(meta[,.(SGV_ID,Study.name,Name)],ameta[,.(Study.name,Name,species_id)],
            by=c('Study.name','Name'),all.x=TRUE,sort=FALSE)
mi<-match(meta$SGV_ID,smap$SGV_ID); keep_map<-!is.na(smap$species_id[mi])
meta<-meta[keep_map];svmat<-svmat[keep_map,,drop=FALSE]
species_raw<-species_raw[smap$species_id[mi[keep_map]],,drop=FALSE]
rownames(species_raw)<-meta$SGV_ID
eligible<-meta$BMI_Group1%in%c('Ob','Ow','Normal') & meta$Disease%in%c('CTRL','CRC') &
          is.finite(meta$Age)&!is.na(meta$Sex)&meta$Sex!=''
meta<-meta[eligible];svmat<-svmat[eligible,,drop=FALSE];species_raw<-species_raw[eligible,,drop=FALSE]
rownames(svmat)<-meta$SGV_ID;rownames(species_raw)<-meta$SGV_ID

dev_ids<-split[split=='development',sample_id]
ob_dev<-match(dev_ids,meta$SGV_ID);ob_test<-match(split[split=='test',sample_id],meta$SGV_ID)
stopifnot(!anyNA(ob_dev),!anyNA(ob_test),all(meta$BMI_Group1[ob_dev]=='Ob'),all(meta$BMI_Group1[ob_test]=='Ob'))
eval_ix<-list(Ob_test=ob_test,Nw=which(meta$BMI_Group1=='Normal'),Ow=which(meta$BMI_Group1=='Ow'))
for(g in names(eval_ix)) stopifnot(length(eval_ix[[g]])>0L,length(unique(meta$Disease[eval_ix[[g]]]))==2L)

# Common stratified resamples within the original Ob development pool.
split_index <- function(seed) {
 z<-meta[ob_dev,.(idx=seq_along(ob_dev),Study=Study.name,Disease,u=NA_real_)]
 set.seed(seed);z[,u:=runif(.N)]
 z[,n_train:=if(.N>1L)max(1L,min(.N-1L,round(.8*.N)))else .N,by=.(Study,Disease)]
 z[,Train:=rank(u,ties.method='first')<=n_train,by=.(Study,Disease)]
 z[,.(idx,Train)]
}

fit_score <- function(model,xtr,ytr,xev,params,seed) {
 ytr<-factor(ytr,levels=c('CTRL','CRC'))
 if(model=='Random forest') {
  set.seed(seed)
  fit<-ranger(x=as.data.frame(xtr),y=ytr,probability=TRUE,num.trees=N_TREES,
   mtry=max(1L,floor(sqrt(ncol(xtr)))),min.node.size=5L,
   case.weights=wgt(ytr),importance='none',seed=seed,num.threads=1L)
  return(as.numeric(predict(fit,as.data.frame(xev))$predictions[,'CRC']))
 }
 if(model%in%c('LASSO','Elastic Net')) {
  alpha<-if(model=='LASSO')1 else .5;lam<-as.numeric(params$lambda[1])
  fit<-glmnet::glmnet(as.matrix(xtr),as.integer(ytr=='CRC'),family='binomial',alpha=alpha,
   lambda=lam,weights=as.numeric(wgt(ytr)),standardize=TRUE)
  return(as.numeric(predict(fit,newx=as.matrix(xev),type='response',s=lam)))
 }
 if(model=='Radial SVM') {
  set.seed(seed)
  fit<-e1071::svm(x=as.data.frame(xtr),y=ytr,kernel='radial',
   cost=as.numeric(params$cost[1]),gamma=as.numeric(params$gamma_scale[1])/ncol(xtr),
   scale=TRUE,probability=TRUE,
   class.weights=setNames(length(ytr)/(2*as.numeric(table(ytr))),levels(ytr)))
  pp<-predict(fit,as.data.frame(xev),probability=TRUE)
  return(as.numeric(attr(pp,'probabilities')[,'CRC']))
 }
 if(model=='XGBoost') {
  set.seed(seed)
  dtr<-xgboost::xgb.DMatrix(as.matrix(xtr),label=as.integer(ytr=='CRC'),weight=as.numeric(wgt(ytr)))
  dev<-xgboost::xgb.DMatrix(as.matrix(xev))
  fit<-xgboost::xgb.train(params=list(objective='binary:logistic',eval_metric='auc',
   max_depth=as.integer(params$max_depth[1]),eta=as.numeric(params$eta[1]),
   subsample=.8,colsample_bytree=.8,min_child_weight=1,nthread=1,tree_method='hist'),
   data=dtr,nrounds=as.integer(params$nrounds[1]),verbose=0)
  return(as.numeric(predict(fit,dev)))
 }
 stop('Unsupported model: ',model)
}

predictions<-list(); split_manifest<-list()
for(rep_i in seq_len(N_REPEATS)) {
 message('Training resample ',rep_i,' / ',N_REPEATS)
 ss<-split_index(SPLIT_SEEDS[rep_i]);tr<-which(ss$Train);stopifnot(length(unique(meta$Disease[ob_dev[tr]]))==2L)
 train_rows<-ob_dev[tr]
 pos<-species_raw[train_rows,,drop=FALSE];pos<-pos[pos>0]
 pseudo<-if(length(pos))min(pos)/2 else 1e-10
 clr<-log(species_raw+pseudo);clr<-sweep(clr,1,rowMeans(clr),'-')
 split_manifest[[rep_i]]<-data.table(Repeat=rep_i,Seed=SPLIT_SEEDS[rep_i],
  sample_id=meta$SGV_ID[ob_dev],Study=meta$Study.name[ob_dev],Disease=meta$Disease[ob_dev],
  Split=ifelse(ss$Train,'training','not_training'))
 for(panel in PANEL_ORDER) {
  ff<-markers[Panel==panel][order(Rank),Feature]
  sp<-intersect(ff,colnames(clr));sv<-intersect(ff,colnames(svmat))
  stopifnot(length(c(sp,sv))==length(ff))
  x<-cbind(clr[,sp,drop=FALSE],svmat[,sv,drop=FALSE])[,ff,drop=FALSE]
  med<-vapply(seq_len(ncol(x)),function(j)median(x[train_rows,j][is.finite(x[train_rows,j])],na.rm=TRUE),numeric(1))
  stopifnot(all(is.finite(med)))
  for(j in seq_len(ncol(x))) x[!is.finite(x[,j]),j]<-med[j]
  xt<-cbind(x,Age=meta$Age,SexMale=as.integer(meta$Sex=='Male'))
  ytr<-meta$Disease[train_rows]
  for(model in MODELS) {
   params<-if(model=='Random forest')data.table() else selected_tuning[Panel==panel&Model==model]
   if(model!='Random forest')stopifnot(nrow(params)==1L)
   scores<-fit_score(model,xt[train_rows,,drop=FALSE],ytr,xt,params,
                     FIT_SEED_BASE+rep_i*100L+match(panel,PANEL_ORDER)*10L+match(model,MODELS))
   for(g in names(eval_ix)) {
    ix<-eval_ix[[g]];yy<-factor(meta$Disease[ix],levels=c('CTRL','CRC'))
    rr<-roc(yy,scores[ix],levels=c('CTRL','CRC'),direction='<',quiet=TRUE)
    predictions[[length(predictions)+1L]]<-data.table(Panel=panel,Model=model,Repeat=rep_i,
      BMI_evaluation=g,sample_id=meta$SGV_ID[ix],Study=meta$Study.name[ix],Disease=as.character(yy),
      score=scores[ix],AUC=as.numeric(auc(rr)),N=length(ix))
   }
  }
 }
}

pred<-rbindlist(predictions);splits<-rbindlist(split_manifest)
auc_rep<-unique(pred[,.(Panel,Model,Repeat,BMI_evaluation,AUC,N)])
auc_summary<-auc_rep[,.(Mean_AUC=mean(AUC),SD_AUC=sd(AUC),Min_AUC=min(AUC),Max_AUC=max(AUC),
                        N_repeats=.N,N=.SD$N[1]),by=.(Panel,Model,BMI_evaluation)]
# Display one observed empirical curve nearest the five-repeat mean AUC, while
# retaining the mean AUC itself as the annotated summary statistic.
representative<-auc_rep[auc_summary,on=.(Panel,Model,BMI_evaluation),nomatch=0L]
representative[,Distance:=abs(AUC-Mean_AUC)]
representative<-representative[order(Distance,Repeat),.SD[1],by=.(Panel,Model,BMI_evaluation)]
rep_ids<-representative[,.(Panel,Model,BMI_evaluation,Representative_repeat=Repeat,
                           Mean_AUC,SD_AUC,N)]
setkey(rep_ids,Panel,Model,BMI_evaluation)
roc_parts<-list()
for(i in seq_len(nrow(rep_ids))) {
 q<-rep_ids[i];z<-pred[Panel==q$Panel&Model==q$Model&BMI_evaluation==q$BMI_evaluation&Repeat==q$Representative_repeat]
 rr<-roc(factor(z$Disease,levels=c('CTRL','CRC')),z$score,levels=c('CTRL','CRC'),direction='<',quiet=TRUE)
 roc_parts[[i]]<-data.table(Panel=q$Panel,Model=q$Model,BMI_evaluation=q$BMI_evaluation,
  Representative_repeat=q$Representative_repeat,FPR=1-rr$specificities,TPR=rr$sensitivities,
  Mean_AUC=q$Mean_AUC,SD_AUC=q$SD_AUC,N=q$N,
  Curve_label=sprintf('%s (n=%s; mean AUC %.3f)',q$BMI_evaluation,q$N,q$Mean_AUC))
}
roc_data<-rbindlist(roc_parts)
setorder(roc_data,Panel,Model,BMI_evaluation,FPR,TPR)
write_tsv(pred,'repeat_predictions');write_tsv(auc_rep,'repeat_AUC');write_tsv(auc_summary,'AUC_summary')
write_tsv(rep_ids,'representative_curve_manifest');write_tsv(splits,'Ob_development_split_manifest')
write_tsv(roc_data,'Fig_transfer_ROC_representative_curves','data_for_figures')

panel_labels<-c(species='Species only (12)',sv='SV only (20)',joint='Species + SV (22)')
model_levels<-MODELS
roc_data[,Panel_label:=factor(unname(panel_labels[Panel]),levels=unname(panel_labels[PANEL_ORDER]))]
roc_data[,Model:=factor(Model,levels=model_levels)]
roc_data[,BMI_label:=factor(BMI_evaluation,levels=BMI_ORDER,labels=c('Ob test','Nw','Ow'))]
ann<-unique(roc_data[,.(Panel_label,Model,BMI_label,Mean_AUC,N)])
ann[,y:=c('Ob test'=.20,'Nw'=.13,'Ow'=.06)[as.character(BMI_label)]]
ann[,Label:=sprintf('%s mean AUC %.3f',BMI_label,Mean_AUC)]
bmi_cols<-c('Ob test'='#C0523D','Nw'='#477DB3','Ow'='#E5A052')
g<-ggplot(roc_data,aes(FPR,TPR,colour=BMI_label,group=interaction(Panel_label,Model,BMI_label)))+
 geom_abline(slope=1,intercept=0,linetype='dashed',colour='#B8B8B8',linewidth=.4)+
 geom_step(direction='vh',linewidth=.55,alpha=.82)+
 geom_text(data=ann,aes(x=.98,y=y,label=Label,colour=BMI_label),inherit.aes=FALSE,
  hjust=1,size=2.05,fontface='bold',show.legend=FALSE)+
 facet_grid(Model~Panel_label,switch='y')+scale_colour_manual(values=bmi_cols,drop=FALSE)+
 coord_equal(xlim=c(0,1),ylim=c(0,1),expand=FALSE)+
 labs(x='False-positive rate',y='True-positive rate',
  title='Transfer of fixed Ob marker panels across BMI groups',
  subtitle='Five Ob-development resamples; curves show the repeat nearest mean AUC',colour=NULL)+
 theme_classic(base_size=8.5,base_family='Arial')+
 theme(axis.line=element_line(linewidth=.35,colour='#333333'),axis.ticks=element_line(linewidth=.3),
  axis.text=element_text(size=6.5,colour='#333333'),axis.title=element_text(size=7.8),
  strip.background=element_rect(fill='#F1F1F1',colour=NA),strip.text=element_text(face='bold',size=7.2),
  strip.text.y.left=element_text(angle=0),strip.placement='outside',plot.title=element_text(face='bold',size=10),
  plot.subtitle=element_text(size=7.4,colour='#555555'),legend.position='bottom',
  legend.text=element_text(size=7),legend.key.width=grid::unit(.4,'cm'),
  panel.spacing=grid::unit(.3,'cm'),plot.margin=margin(5,5,5,5))
stem<-file.path(ROOT,'figures',paste0(PREFIX,'Fig_ROC'))
ggsave(paste0(stem,'.pdf'),g,width=155/25.4,height=170/25.4,device=cairo_pdf,family='Arial')
svglite::svglite(paste0(stem,'.svg'),width=155/25.4,height=170/25.4);print(g);dev.off()
ragg::agg_png(paste0(stem,'.png'),width=155/25.4,height=170/25.4,units='in',res=400);print(g);dev.off()

qc<-data.table(Check=c('Five Ob resamples per panel/model',
 'One representative ROC per panel/model/evaluation group',
 'Representative repeat is closest to mean AUC',
 'AUC means summarize all five fitted repeats',
 'ROC points are empirical and ordered from (0,0) to (1,1)',
 'AUC and plot files are complete'),Passed=c(
  all(auc_rep[, .N,by=.(Panel,Model,BMI_evaluation)]$N==N_REPEATS),
  nrow(rep_ids)==length(PANEL_ORDER)*length(MODELS)*length(BMI_ORDER),
  all(rep_ids$Representative_repeat%in%1:N_REPEATS),
  all(auc_summary$N_repeats==N_REPEATS),
  all(roc_data[,.(ok=FPR[1]==0&&TPR[1]==0&&FPR[.N]==1&&TPR[.N]==1),by=.(Panel,Model,BMI_evaluation)]$ok),
  all(file.exists(c(paste0(stem,c('.pdf','.svg','.png')),
    file.path(ROOT,'data_for_figures',paste0(PREFIX,'Fig_transfer_ROC_representative_curves.tsv')))))))
write_tsv(qc,'QC');stopifnot(all(qc$Passed));print(auc_summary);print(qc)
cat('CROSS_BMI_OTHER_ALGORITHMS_COMPLETE\n')
