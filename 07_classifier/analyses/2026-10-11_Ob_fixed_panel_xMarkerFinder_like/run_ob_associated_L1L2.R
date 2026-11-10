#!/usr/bin/env Rscript
.libPaths(c('D:/R/library', .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(Boruta); library(ranger); library(pROC); library(ggplot2)
})
ROOT <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-10-11_Ob_fixed_panel_xMarkerFinder_like'
OUT_PREFIX <- 'restricted_ObAssoc_L1L2_'
arg <- commandArgs(TRUE)
MODE <- if(length(arg)>=1) arg[1] else 'audit'
BMI <- if(length(arg)>=2) arg[2] else 'Ob'
PANEL_ARG <- if(length(arg)>=3) arg[3] else 'all'
stopifnot(MODE %in% c('audit','smoke','full','summarize'),BMI%in%c('Ob','Ow','Nw'))
stopifnot(PANEL_ARG%in%c('species','sv','joint','all'))
for(d in c('figures','tables','data_for_figures','checkpoints','logs'))
 dir.create(file.path(ROOT,d),recursive=TRUE,showWarnings=FALSE)

PATHS <- c(
 dsv='H:/My Drive/proj_ob_AA_CRC/02.processed_data/05_SV_result/all_BMI_crossBMI/all_BMI_dsgv.csv',
 vsv='H:/My Drive/proj_ob_AA_CRC/02.processed_data/05_SV_result/all_BMI_crossBMI/all_BMI_vsgv.csv',
 sv_meta='H:/My Drive/proj_ob_AA_CRC/02.processed_data/05_SV_result/all_BMI_crossBMI/metadata_allBMI_SGV_aligned.rds',
 abundance='H:/My Drive/proj_ob_AA_CRC/03.results/04_differential_abundance/heatmap_final/cache/abd_corrected_meta_clean.rds',
 feature_map='H:/My Drive/proj_ob_AA_CRC/analyses/2026-09-28_allBMI_BMIstratified_SV_landscape/tables/01_allBMI_full_SV_feature_to_SGB_map.tsv',
 callability='H:/My Drive/proj_ob_AA_CRC/analyses/2026-09-28_allBMI_BMIstratified_SV_landscape/tables/02_BMIstratified_SGB_callability.tsv',
 ob_species='H:/My Drive/proj_ob_AA_CRC/analyses/2026-07-13_bmi_union_overlap_heatmaps/tables/Feature_List_Ob_Refined.csv',
 ob_l1_l2='H:/My Drive/proj_ob_AA_CRC/analyses/2026-09-29_BMI_enriched_SV_function_modelQC/2026-09-30_BMI_specificity_layers/tables/L1_L2_candidate_master.tsv')
stopifnot(all(file.exists(PATHS)))

SPLIT_SEED<-20261011L;TEST_FRACTION<-.20
INNER_FOLDS<-5L;INNER_REPEATS<-5L;INNER_SEEDS<-31001L+seq_len(INNER_REPEATS)
RANK_REPEATS<-20L;RANK_SEEDS<-41001L+seq_len(RANK_REPEATS)
BORUTA_SEED<-51001L;BORUTA_MAX_RUNS<-100L;BORUTA_P<-.01
N_TREES<-500L;COR_THRESHOLD<-.80
MIN_SPECIES_PREV<-.10;MIN_SV_OBS<-.20;MIN_PER_CLASS<-20L;MIN_STUDIES<-3L;MIN_DSV_MINOR<-5L
PANELS<-c(species='Species only',sv='SV only',joint='Species + SV')
PANELS_RUN<-if(PANEL_ARG=='all')names(PANELS)else PANEL_ARG
wt<-function(x,n,w='tables')fwrite(x,file.path(ROOT,w,paste0(OUT_PREFIX,n,'.tsv')),sep='\t',na='NA')
wgt<-function(y)length(y)/(2*as.numeric(table(y)[as.character(y)]))
aucx<-function(y,p)as.numeric(auc(roc(y,p,levels=c('CTRL','CRC'),direction='<',quiet=TRUE)))

# Load only mapped SV columns, then align the species and SV sample spaces.
fm<-fread(PATHS['feature_map']);ca<-fread(PATHS['callability'])
ca[,Retained:=as.character(Retained)%in%c('TRUE','True','true','1')]
sgb_union<-unique(ca[BMI_Group%in%c('Ob','Ow','Nw')&Retained==TRUE,SGB])
fm<-fm[SGB%in%sgb_union&SV_type%in%c('dSGV','vSGV')]
stopifnot(!anyDuplicated(fm[,.(SV_feature,SV_type)]))
fm[,Feature:=sprintf('%s%05d',ifelse(SV_type=='dSGV','DSV','VSV'),rowid(SV_type))]
read_sv<-function(path,type){
 mp<-fm[SV_type==type]
 z<-fread(path,select=c('V1',mp$SV_feature),check.names=FALSE,showProgress=FALSE)
 setnames(z,1L,'sample_id');old<-setdiff(names(z),'sample_id')
 setnames(z,old,mp[match(old,SV_feature),Feature]);z
}
message('Reading SV matrices...')
dd<-read_sv(PATHS['dsv'],'dSGV');vv<-read_sv(PATHS['vsv'],'vSGV')
stopifnot(identical(dd$sample_id,vv$sample_id))
meta<-as.data.table(readRDS(PATHS['sv_meta']));stopifnot(!anyDuplicated(meta$SGV_ID))
ii<-which(dd$sample_id%in%meta$SGV_ID);stopifnot(length(ii)==nrow(meta))
dd<-dd[ii];vv<-vv[ii];sid<-dd$sample_id
dmat<-1-as.matrix(dd[,-1]);storage.mode(dmat)<-'double'
vmat<-as.matrix(vv[,-1]);storage.mode(vmat)<-'double'
svmat<-cbind(dmat,vmat);rownames(svmat)<-sid;rm(dd,vv,dmat,vmat);gc()
meta<-meta[match(sid,SGV_ID)];stopifnot(identical(meta$SGV_ID,sid))

aa<-readRDS(PATHS['abundance']);abd<-as.matrix(aa$abd_corrected)
storage.mode(abd)<-'double';abd[!is.finite(abd)|abd<0]<-0
spmap<-data.table(Feature=sprintf('SP%05d',seq_len(nrow(abd))),Original_feature=rownames(abd))
rownames(abd)<-spmap$Feature;species_raw<-t(abd)
if(max(species_raw)>1.5)species_raw<-species_raw/100
ameta<-as.data.table(aa$meta_clean);ameta[,species_id:=colnames(aa$abd_corrected)]
smap<-merge(meta[,.(SGV_ID,Study.name,Name)],ameta[,.(Study.name,Name,species_id)],
 by=c('Study.name','Name'),all.x=TRUE,sort=FALSE)
stopifnot(!anyDuplicated(smap$SGV_ID));mi<-match(meta$SGV_ID,smap$SGV_ID)
ok<-!is.na(smap$species_id[mi]);meta<-meta[ok];svmat<-svmat[ok,,drop=FALSE]
species_raw<-species_raw[smap$species_id[mi[ok]],,drop=FALSE];rownames(species_raw)<-meta$SGV_ID
rm(abd,aa);gc()
keep<-meta$BMI_Group1==BMI&meta$Disease%in%c('CTRL','CRC')&
 is.finite(meta$Age)&!is.na(meta$Sex)&meta$Sex!=''
meta<-meta[keep];svmat<-svmat[keep,,drop=FALSE];species_raw<-species_raw[keep,,drop=FALSE]
meta[,row_id:=seq_len(.N)]
stopifnot(identical(rownames(svmat),meta$SGV_ID),identical(rownames(species_raw),meta$SGV_ID))
manifest<-rbindlist(list(
 spmap[,.(Feature,Type='species',Original_feature,SGB=NA_character_,SV_type=NA_character_)],
 fm[,.(Feature,Type='SV',Original_feature=SV_feature,SGB,SV_type)]),fill=TRUE)
obsp<-fread(PATHS['ob_species'])
stopifnot(all(c('feature','overlap')%in%names(obsp)),!anyDuplicated(obsp$feature))
obsv<-fread(PATHS['ob_l1_l2'])[BMI_enriched_group=='Ob']
stopifnot(all(c('SV_feature','L2_nominal','L2_selected_BH')%in%names(obsv)),!anyDuplicated(obsv$SV_feature))
sp_ann<-obsp[,.(Original_feature=feature,Ob_species_class=overlap)]
sv_ann<-obsv[,.(Original_feature=SV_feature,Ob_L1=TRUE,L2_nominal,L2_selected_BH)]
manifest<-merge(manifest,sp_ann,by='Original_feature',all.x=TRUE,sort=FALSE)
manifest<-merge(manifest,sv_ann,by='Original_feature',all.x=TRUE,sort=FALSE)
species_candidate_features<-manifest[Type=='species'&!is.na(Ob_species_class),Feature]
sv_candidate_features<-manifest[Type=='SV'&Ob_L1==TRUE,Feature]
candidate_summary<-data.table(
 Candidate_set=c('Ob-associated species','Ob L1 SV','Ob L2 nominal SV','Ob L2 BH-selected SV'),
 N=c(length(species_candidate_features),length(sv_candidate_features),
     obsv[L2_nominal==TRUE,.N],obsv[L2_selected_BH==TRUE,.N]))
wt(manifest,'feature_manifest')
wt(candidate_summary,'prespecified_candidate_pool')
wt(meta[,.(N=.N),by=.(BMI=BMI_Group1,Study=Study.name,Disease)],'population_counts')

# One fixed, prespecified test split; marker selection is development-only.
set.seed(SPLIT_SEED)
split_dt<-meta[,.(row_id,sample_id=SGV_ID,Study=Study.name,Disease,u=runif(.N))]
split_dt[,n_test:=if(.N>=2L)max(1L,round(TEST_FRACTION*.N))else 0L,by=.(Study,Disease)]
split_dt[,split:=ifelse(rank(u,ties.method='first')<=n_test[1L],'test','development'),by=.(Study,Disease)]
split_dt[,c('u','n_test'):=NULL];wt(split_dt,'fixed_development_test_split')
dev<-split_dt[split=='development',row_id];test<-split_dt[split=='test',row_id]
stopifnot(all(c('CTRL','CRC')%in%meta$Disease[dev]),all(c('CTRL','CRC')%in%meta$Disease[test]))

# Fit label-independent transformations and adequacy filters on development only.
prepare_data<-function(dev,test){
 y<-factor(meta$Disease[dev],levels=c('CTRL','CRC'))
 pos<-species_raw[dev,][species_raw[dev,]>0];pseudo<-if(length(pos))min(pos)/2 else 1e-10
 zz<-log(species_raw[c(dev,test),,drop=FALSE]+pseudo);zz<-sweep(zz,1,rowMeans(zz),'-')
 sxdev<-zz[seq_along(dev),,drop=FALSE];sxtest<-zz[length(dev)+seq_along(test),,drop=FALSE]
 prev<-colMeans(species_raw[dev,,drop=FALSE]>0)
 sp_ok<-intersect(names(prev)[prev>=MIN_SPECIES_PREV&apply(sxdev,2,function(v)is.finite(sd(v))&&sd(v)>0)],
                  species_candidate_features)
 retained<-character();sgb_audit<-list()
 for(sgb in unique(fm$SGB)){
  ff<-fm[SGB==sgb,Feature];comp<-rowMeans(is.finite(svmat[dev,ff,drop=FALSE]));pr<-mean(comp>=.5)
  if(pr>=.3)retained<-c(retained,sgb)
  sgb_audit[[sgb]]<-data.table(SGB=sgb,N_features=length(ff),Callable_prevalence=pr,Retained=pr>=.3)
 }
 sm<-fm[SGB%in%retained&Feature%in%sv_candidate_features]
 sv_audit<-rbindlist(lapply(seq_len(nrow(sm)),function(i){
  f<-sm$Feature[i];v<-svmat[dev,f];obs<-is.finite(v)
  inf<-data.table(Study=meta$Study.name[dev],Disease=meta$Disease[dev],Observed=obs)[,
   .(hasCTRL=any(Observed&Disease=='CTRL'),hasCRC=any(Observed&Disease=='CRC')),by=Study]
  minor<-if(sm$SV_type[i]=='dSGV'&&any(obs))min(table(factor(v[obs],levels=c(0,1))))else NA_integer_
  good<-mean(obs)>=MIN_SV_OBS&&sum(obs&y=='CTRL')>=MIN_PER_CLASS&&sum(obs&y=='CRC')>=MIN_PER_CLASS&&
   sum(inf$hasCTRL&inf$hasCRC)>=MIN_STUDIES&&length(unique(v[obs]))>1&&
   (sm$SV_type[i]!='dSGV'||minor>=MIN_DSV_MINOR)
  data.table(Feature=f,SGB=sm$SGB[i],SV_type=sm$SV_type[i],Observed_fraction=mean(obs),
   N_CTRL=sum(obs&y=='CTRL'),N_CRC=sum(obs&y=='CRC'),N_informative_studies=sum(inf$hasCTRL&inf$hasCRC),
   Minor_state=minor,Eligible=good)
 }))
 sv_ok<-sv_audit[Eligible==TRUE,Feature]
 list(xdev=cbind(sxdev[,sp_ok,drop=FALSE],svmat[dev,sv_ok,drop=FALSE]),
      xtest=cbind(sxtest[,sp_ok,drop=FALSE],svmat[test,sv_ok,drop=FALSE]),
      y=y,sp_ok=sp_ok,sv_ok=sv_ok,pseudo=pseudo,sgb_audit=rbindlist(sgb_audit),sv_audit=sv_audit)
}
impute_pair<-function(a,b){
 a<-as.matrix(a);b<-as.matrix(b)
 for(j in seq_len(ncol(a))){m<-median(a[is.finite(a[,j]),j]);stopifnot(is.finite(m))
  a[!is.finite(a[,j]),j]<-m;b[!is.finite(b[,j]),j]<-m}
 list(a=a,b=b)
}
prep<-prepare_data(dev,test);wt(prep$sgb_audit,'development_SGB_callability_audit')
wt(prep$sv_audit,'development_SV_eligibility_audit')
audit<-data.table(BMI=BMI,Split_seed=SPLIT_SEED,N_development=length(dev),N_test=length(test),
 Development_CTRL=sum(meta$Disease[dev]=='CTRL'),Development_CRC=sum(meta$Disease[dev]=='CRC'),
 Test_CTRL=sum(meta$Disease[test]=='CTRL'),Test_CRC=sum(meta$Disease[test]=='CRC'),
 Eligible_species=length(prep$sp_ok),Eligible_SV=length(prep$sv_ok),Test_used_for_Boruta_or_model_fitting=FALSE,
 Upstream_candidate_lists_used_full_dataset=TRUE,
 Validation_status='Exploratory: upstream outcome-based candidate lists were derived in the full Ob dataset',
 Panel_size_rule='smallest k within one SE of maximum development CV AUC')
wt(audit,'design_audit')
if(MODE=='audit'){print(audit);cat('AUDIT_COMPLETE\n');quit(save='no')}

make_folds<-function(seed,nfolds=INNER_FOLDS){
 set.seed(seed);z<-meta[dev,.(row_id,Study=Study.name,Disease,u=runif(.N))]
 z[,Fold:=((rank(u,ties.method='first')-1L)%%nfolds)+1L,by=.(Study,Disease)]
 z[,u:=NULL];z
}
fit_rf<-function(xtr,ytr,xte,seed,importance='none'){
 set.seed(seed);fit<-ranger(x=as.data.frame(xtr),y=ytr,probability=TRUE,num.trees=N_TREES,
  mtry=max(1L,floor(sqrt(ncol(xtr)))),min.node.size=5L,case.weights=wgt(ytr),
  importance=importance,seed=seed,num.threads=1L)
 list(fit=fit,pred=as.numeric(predict(fit,as.data.frame(xte))$predictions[,'CRC']))
}
run_boruta<-function(x,y,seed){
 set.seed(seed)
 fit<-Boruta(as.data.frame(x),y,pValue=BORUTA_P,mcAdj=TRUE,maxRuns=BORUTA_MAX_RUNS,
  holdHistory=FALSE,doTrace=0,getImp=getImpRfZ,num.trees=N_TREES,num.threads=1,
  case.weights=wgt(y))
 data.table(Feature=names(fit$finalDecision),Decision=as.character(fit$finalDecision))
}
rank_features<-function(x,y,features){
 out<-lapply(seq_along(RANK_SEEDS),function(i){
  fit<-fit_rf(x[,features,drop=FALSE],y,x[,features,drop=FALSE],RANK_SEEDS[i],'permutation')$fit
  imp<-fit$variable.importance
  data.table(Feature=names(imp),Repeat=i,Seed=RANK_SEEDS[i],Importance=as.numeric(imp),
   Rank=frank(-as.numeric(imp),ties.method='average'))
 })
 long<-rbindlist(out)
 sm<-long[,.(Mean_importance=mean(Importance),SD_importance=sd(Importance),
  Mean_rank=mean(Rank),Median_rank=median(Rank)),by=Feature]
 setorder(sm,Mean_rank,-Mean_importance,Feature);sm[,Initial_rank:=seq_len(.N)]
 list(long=long,summary=sm)
}
prune_cor<-function(x,ranking){
 ord<-ranking$Feature
 if(length(ord)==1)return(data.table(Feature=ord,Correlation_keep=TRUE,Correlated_with=NA_character_,Abs_rho=NA_real_))
 cm<-suppressWarnings(cor(x[,ord,drop=FALSE],method='spearman',use='pairwise.complete.obs'))
 kept<-character();ans<-vector('list',length(ord))
 for(i in seq_along(ord)){
  f<-ord[i]
  if(!length(kept)){ans[[i]]<-data.table(Feature=f,Correlation_keep=TRUE,Correlated_with=NA_character_,Abs_rho=NA_real_);kept<-f;next}
  rr<-abs(cm[f,kept]);rr[!is.finite(rr)]<-0;j<-which.max(rr);keep<-rr[j]<COR_THRESHOLD
  ans[[i]]<-data.table(Feature=f,Correlation_keep=keep,
   Correlated_with=if(keep)NA_character_ else kept[j],Abs_rho=rr[j])
  if(keep)kept<-c(kept,f)
 }
 rbindlist(ans)
}
size_curve<-function(xraw,y,ordered,panel,smoke=FALSE){
 seeds<-if(smoke)INNER_SEEDS[1] else INNER_SEEDS
 nf<-if(smoke)3L else INNER_FOLDS
 fold_manifest<-rbindlist(lapply(seq_along(seeds),function(r){
  z<-make_folds(seeds[r],nf);z[,':='(Repeat=r,Seed=seeds[r],Panel=panel)];z
 }))
 curves<-vector('list',length(seeds));ks<-seq_along(ordered)
 for(r in seq_along(seeds)){
  fd<-fold_manifest[Repeat==r];pm<-matrix(NA_real_,length(y),length(ks))
  for(kfold in seq_len(nf)){
   va<-match(fd[Fold==kfold,row_id],dev);tr<-setdiff(seq_along(dev),va)
   stopifnot(!anyNA(va),length(unique(y[tr]))==2,length(unique(y[va]))==2)
   q<-impute_pair(xraw[tr,ordered,drop=FALSE],xraw[va,ordered,drop=FALSE])
   a<-cbind(q$a,Age=meta$Age[dev[tr]],SexMale=as.integer(meta$Sex[dev[tr]]=='Male'))
   b<-cbind(q$b,Age=meta$Age[dev[va]],SexMale=as.integer(meta$Sex[dev[va]]=='Male'))
   for(k in ks){ff<-c(ordered[seq_len(k)],'Age','SexMale')
    pm[va,k]<-fit_rf(a[,ff,drop=FALSE],y[tr],b[,ff,drop=FALSE],seeds[r]+1000L*kfold+k)$pred}
  }
  curves[[r]]<-data.table(Panel=panel,Repeat=r,Seed=seeds[r],k=ks,
   AUC=vapply(ks,function(k)aucx(y,pm[,k]),numeric(1)))
 }
 long<-rbindlist(curves);sm<-long[,.(Mean_AUC=mean(AUC),SD_AUC=sd(AUC),
  SE_AUC=if(.N>1)sd(AUC)/sqrt(.N) else NA_real_,N_repeats=.N),by=.(Panel,k)]
 best<-sm[which.max(Mean_AUC)];se<-if(is.finite(best$SE_AUC))best$SE_AUC else 0
 threshold<-best$Mean_AUC-se;selected_k<-min(sm[Mean_AUC>=threshold,k])
 sm[,':='(Best_mean_AUC=best$Mean_AUC,One_SE_threshold=threshold,Selected_k=selected_k,
  Selected_by_one_SE=k==selected_k)]
 list(long=long,summary=sm,folds=fold_manifest,selected_k=selected_k)
}

run_panel<-function(panel,smoke=FALSE){
  cp<-file.path(ROOT,'checkpoints',paste0(BMI,'__',OUT_PREFIX,panel,'.rds'))
 if(!smoke&&file.exists(cp)){message('Reuse ',basename(cp));return(readRDS(cp))}
 pool<-switch(panel,species=prep$sp_ok,sv=prep$sv_ok,joint=c(prep$sp_ok,prep$sv_ok))
 xd<-prep$xdev[,pool,drop=FALSE];xt<-prep$xtest[,pool,drop=FALSE];allimp<-impute_pair(xd,xt)
 if(smoke){
  confirmed<-head(pool,min(8L,length(pool)));dec<-data.table(Feature=confirmed,Decision='Smoke_only')
  ranks<-list(long=data.table(Feature=confirmed,Repeat=1L,Seed=1L,Importance=rev(seq_along(confirmed)),Rank=seq_along(confirmed)),
   summary=data.table(Feature=confirmed,Mean_importance=rev(seq_along(confirmed)),SD_importance=NA_real_,
    Mean_rank=seq_along(confirmed),Median_rank=seq_along(confirmed),Initial_rank=seq_along(confirmed)))
 }else{
  dec<-run_boruta(allimp$a,prep$y,BORUTA_SEED+1000L*match(panel,names(PANELS)))
  confirmed<-dec[Decision=='Confirmed',Feature];if(!length(confirmed))stop('No confirmed features: ',panel)
  ranks<-rank_features(allimp$a,prep$y,confirmed)
 }
 ca<-prune_cor(allimp$a,ranks$summary);ordered<-ca[Correlation_keep==TRUE,Feature]
 if(smoke)ordered<-head(ordered,min(4L,length(ordered)))
 curve<-size_curve(xd,prep$y,ordered,panel,smoke);selected<-ordered[seq_len(curve$selected_k)]
 q<-impute_pair(xd[,selected,drop=FALSE],xt[,selected,drop=FALSE])
 a<-cbind(q$a,Age=meta$Age[dev],SexMale=as.integer(meta$Sex[dev]=='Male'))
 b<-cbind(q$b,Age=meta$Age[test],SexMale=as.integer(meta$Sex[test]=='Male'))
 fit<-fit_rf(a,prep$y,b,61001L+match(panel,names(PANELS)),'permutation')
 yt<-factor(meta$Disease[test],levels=c('CTRL','CRC'))
 ro<-roc(yt,fit$pred,levels=c('CTRL','CRC'),direction='<',quiet=TRUE);ci<-as.numeric(ci.auc(ro,method='delong'))
 out<-list(panel=panel,pool=pool,decisions=dec,confirmed=confirmed,ranking_long=ranks$long,
  ranking=ranks$summary,correlation=ca,ordered=ordered,curve=curve,selected=selected,
  predictions=data.table(Panel=panel,sample_id=meta$SGV_ID[test],Study=meta$Study.name[test],
   Disease=meta$Disease[test],score=fit$pred),
  summary=data.table(Panel=panel,Panel_label=unname(PANELS[panel]),N_markers=length(selected),
   N_species=sum(selected%in%prep$sp_ok),N_SV=sum(selected%in%prep$sv_ok),Test_AUC=as.numeric(auc(ro)),
   CI_low=ci[1],CI_high=ci[3],N_test=length(test),N_CTRL=sum(yt=='CTRL'),N_CRC=sum(yt=='CRC')),
  roc=data.table(Panel=panel,FPR=1-ro$specificities,TPR=ro$sensitivities,Threshold=ro$thresholds),
  importance=data.table(Panel=panel,Feature=names(fit$fit$variable.importance),
   Importance=as.numeric(fit$fit$variable.importance)))
 if(!smoke)saveRDS(out,cp);out
}

if(MODE=='smoke'){
 z<-run_panel('joint',TRUE)
 stopifnot(all(is.finite(z$predictions$score)),length(z$selected)==z$curve$selected_k)
 wt(z$curve$summary,'smoke_panel_size_curve');wt(z$summary,'smoke_test_summary')
 cat('SMOKE_COMPLETE\n');quit(save='no')
}
if(MODE=='full'){
  lock<-file.path(ROOT,paste0('run_',OUT_PREFIX,BMI,'_',PANEL_ARG,'.lock'))
 if(!dir.create(lock,showWarnings=FALSE))stop('Run lock exists: ',lock)
 writeLines(as.character(Sys.getpid()),file.path(lock,'pid.txt'))
 tryCatch(invisible(lapply(PANELS_RUN,run_panel)),finally={
  pf<-file.path(lock,'pid.txt');if(file.exists(pf))unlink(pf)
  stopifnot(identical(normalizePath(dirname(lock)),normalizePath(ROOT)))
  unlink(lock,recursive=TRUE)
 })
}

available<-names(PANELS)[file.exists(file.path(ROOT,'checkpoints',paste0(BMI,'__',OUT_PREFIX,names(PANELS),'.rds')))]
if(MODE=='summarize'&&PANEL_ARG!='all')available<-intersect(available,PANEL_ARG)
if(MODE=='full'&&PANEL_ARG=='all'&&length(available)<3L)stop('Not all checkpoints completed')
if(!length(available))stop('No completed checkpoints')
res<-lapply(available,function(p)readRDS(file.path(ROOT,'checkpoints',paste0(BMI,'__',OUT_PREFIX,p,'.rds'))))

dec<-rbindlist(lapply(res,function(z){q<-copy(z$decisions);q[,Panel:=z$panel];q}),fill=TRUE)
rank_long<-rbindlist(lapply(res,function(z){q<-copy(z$ranking_long);q[,Panel:=z$panel];q}))
ranked<-rbindlist(lapply(res,function(z){q<-merge(copy(z$ranking),z$correlation,by='Feature',all.x=TRUE)
 q[,':='(Panel=z$panel,Selected=Feature%in%z$selected)];q}))
curve_long<-rbindlist(lapply(res,function(z)z$curve$long))
curve_summary<-rbindlist(lapply(res,function(z)z$curve$summary))
folds<-rbindlist(lapply(res,function(z)z$curve$folds))
pred<-rbindlist(lapply(res,function(z)z$predictions))
test_summary<-rbindlist(lapply(res,function(z)z$summary))
roc_points<-rbindlist(lapply(res,function(z)z$roc))
importance<-rbindlist(lapply(res,function(z)z$importance))
markers<-rbindlist(lapply(res,function(z)data.table(Panel=z$panel,Rank=seq_along(z$selected),Feature=z$selected)))
markers<-merge(markers,manifest,by='Feature',all.x=TRUE,sort=FALSE);setorder(markers,Panel,Rank)

wt(dec,'Boruta_decisions');wt(rank_long,'RF_importance_repeats')
wt(ranked,'feature_ranking_and_correlation_pruning')
wt(curve_long,'development_CV_panel_size_repeats');wt(curve_summary,'development_CV_panel_size_summary')
wt(folds,'development_CV_fold_manifest');wt(markers,'locked_marker_panels')
wt(pred,'locked_panel_test_predictions');wt(test_summary,'locked_panel_test_AUC')

# Secondary paired DeLong comparison; the one-sided direction is prespecified as joint > species.
cmp<-data.table()
if(all(c('species','joint')%in%available)){
 wide<-dcast(pred[Panel%in%c('species','joint')],sample_id+Disease~Panel,value.var='score')
 yy<-factor(wide$Disease,levels=c('CTRL','CRC'))
 rs<-roc(yy,wide$species,levels=c('CTRL','CRC'),direction='<',quiet=TRUE)
 rj<-roc(yy,wide$joint,levels=c('CTRL','CRC'),direction='<',quiet=TRUE)
 cmp<-data.table(Contrast='Joint vs species',AUC_species=as.numeric(auc(rs)),AUC_joint=as.numeric(auc(rj)),
  Delta_AUC=as.numeric(auc(rj)-auc(rs)),
  P_two_sided=as.numeric(roc.test(rs,rj,paired=TRUE,method='delong',alternative='two.sided')$p.value),
  P_one_sided_joint_greater=as.numeric(roc.test(rs,rj,paired=TRUE,method='delong',alternative='less')$p.value))
}
wt(cmp,'paired_test_AUC_comparisons')

# Exact plotting data and compact vector/raster exports.
theme_pub<-function()theme_classic(base_size=12,base_family='Arial')+theme(
 axis.line=element_line(linewidth=.45,colour='#333333'),
 strip.background=element_rect(fill='#F3F3F3',colour=NA),strip.text=element_text(face='bold'),
 plot.title=element_text(face='bold'),legend.title=element_blank())
save_plot<-function(g,stem,w,h){
 stem<-paste0(OUT_PREFIX,stem)
 grDevices::cairo_pdf(file.path(ROOT,'figures',paste0(stem,'.pdf')),
  width=w,height=h,family='Arial');print(g);grDevices::dev.off()
 ggsave(file.path(ROOT,'figures',paste0(stem,'.png')),g,width=w,height=h,dpi=300)
 if(requireNamespace('svglite',quietly=TRUE)){
  svglite::svglite(file.path(ROOT,'figures',paste0(stem,'.svg')),width=w,height=h)
  print(g);grDevices::dev.off()
 }
}
pd<-copy(curve_summary);pd[,Panel_label:=factor(PANELS[Panel],levels=PANELS)]
wt(pd,'Fig_panel_size_selection','data_for_figures')
g1<-ggplot(pd,aes(k,Mean_AUC))+geom_ribbon(aes(ymin=Mean_AUC-SE_AUC,ymax=Mean_AUC+SE_AUC),
 fill='#A9C8C3',alpha=.45)+geom_line(linewidth=.8,colour='#356E70')+
 geom_point(data=pd[Selected_by_one_SE==TRUE],size=2.8,shape=21,stroke=.7,fill='#D8614C',colour='#333333')+
 geom_hline(aes(yintercept=One_SE_threshold),linetype='dashed',linewidth=.45,colour='#777777')+
 facet_wrap(~Panel_label,scales='free_x',nrow=1)+
 labs(x='Number of markers',y='Development CV AUC',title=paste0(BMI,' fixed-panel size selection'))+theme_pub()
save_plot(g1,'Fig_panel_size_selection',7.2,3.4)

rd<-merge(roc_points,test_summary[,.(Panel,Panel_label,Test_AUC,N_markers)],by='Panel')
setorder(rd,Panel,FPR,TPR)
rd[,Curve_label:=sprintf('%s (%d markers; AUC %.3f)',Panel_label,N_markers,Test_AUC)]
wt(rd,'Fig_test_ROC','data_for_figures')
base_cols<-c('Species only'='#777777','SV only'='#4E8C8A','Species + SV'='#C0523D')
labs_now<-test_summary[,sprintf('%s (%d markers; AUC %.3f)',Panel_label,N_markers,Test_AUC)]
cols<-setNames(base_cols[test_summary$Panel_label],labs_now)
g2<-ggplot(rd,aes(FPR,TPR,colour=Curve_label))+geom_abline(slope=1,intercept=0,linetype='dashed',
 colour='#B8B8B8',linewidth=.5)+geom_step(linewidth=.9,direction='vh')+
 scale_colour_manual(values=cols)+coord_equal(xlim=c(0,1),ylim=c(0,1),expand=FALSE)+
 labs(x='False-positive rate',y='True-positive rate',title=paste0(BMI,' CRC classification: fixed test set'))+
 theme_pub()+theme(legend.position=c(.98,.05),legend.justification=c(1,0),
 legend.background=element_rect(fill='white',colour=NA))
save_plot(g2,'Fig_test_ROC',5.2,4.6)

imp<-importance[Panel=='joint'&Feature%in%markers[Panel=='joint',Feature]]
imp<-merge(imp,markers[,.(Panel,Feature,Rank,Type,Original_feature,SGB,SV_type)],
 by=c('Panel','Feature'),all.x=TRUE)
imp[,Display:=fifelse(
 Type=='species',
 gsub('_',' ',sub('\\|t__',' ',sub('^s__','',Original_feature))),
 paste0(sub('^t__','',SGB),' ',SV_type,': ',gsub('_','-',sub('.*:','',Original_feature)))
)]
imp[,Panel_label:=factor(PANELS[Panel],levels=PANELS)]
imp[,Display:=factor(Display,levels=unique(Display[order(Importance)]))]
wt(imp,'Fig_locked_marker_importance','data_for_figures')
g3<-ggplot(imp,aes(Importance,Display,colour=Type))+
 geom_segment(aes(x=0,xend=Importance,yend=Display),linewidth=.55,colour='#C9C9C9')+
 geom_point(size=2.2)+
 scale_colour_manual(values=c(species='#356E70',SV='#C0523D'))+
 labs(x='Permutation importance',y=NULL,title=paste0(BMI,' locked joint marker panel'))+
 theme_pub()+theme(axis.text.y=element_text(size=9),legend.position='top',
  plot.margin=margin(6,10,6,6))
save_plot(g3,'Fig_locked_marker_importance',7.2,max(4.8,1.8+.20*nrow(imp)))

qc<-data.table(Check=c('No development-test overlap','Finite test predictions',
 'Locked markers came from eligible pools','Figure plotting data exist'),
 Passed=c(!length(intersect(dev,test)),all(is.finite(pred$score)),
  all(vapply(res,function(z)all(z$selected%in%z$pool),logical(1))),
  all(file.exists(file.path(ROOT,'data_for_figures',paste0(OUT_PREFIX,c('Fig_panel_size_selection.tsv',
   'Fig_test_ROC.tsv','Fig_locked_marker_importance.tsv')))))))
wt(qc,'QC');stopifnot(all(qc$Passed))
capture.output(sessionInfo(),file=file.path(ROOT,'tables',paste0(OUT_PREFIX,'sessionInfo.txt')))
writeLines(list.files(ROOT,recursive=TRUE,full.names=TRUE),file.path(ROOT,paste0('output_files_',OUT_PREFIX,'.txt')))
cat('COMPLETED\n')




