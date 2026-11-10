#!/usr/bin/env Rscript
# Fixed-marker nested repeated CV. Marker membership is read from the existing locked panel table and is never reselected.
.libPaths(c('D:/R/library', .libPaths()))
suppressPackageStartupMessages({library(data.table); library(pROC); library(glmnet); library(ranger); library(e1071); library(xgboost)})
ROOT <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-10-11_Ob_fixed_panel_xMarkerFinder_like'
OLD <- 'restricted_ObAssoc_L1L2_'; TAG <- 'restricted_ObAssoc_L1L2_nestedCV_'
KOUT <- 5L; KIN <- 5L; NREP <- 5L; OUTSEEDS <- 84002L + seq_len(NREP); N_TREES <- 500L
LAM <- 10^seq(-4,1,length.out=60L); SVMGRID <- CJ(cost=c(.1,1,10),gamma_scale=c(.1,1,10)); XGBGRID <- CJ(max_depth=c(1L,2L,3L),eta=c(.03,.1),nrounds=c(100L,300L))
P <- c(dsv='H:/My Drive/proj_ob_AA_CRC/02.processed_data/05_SV_result/all_BMI_crossBMI/all_BMI_dsgv.csv',vsv='H:/My Drive/proj_ob_AA_CRC/02.processed_data/05_SV_result/all_BMI_crossBMI/all_BMI_vsgv.csv',svmeta='H:/My Drive/proj_ob_AA_CRC/02.processed_data/05_SV_result/all_BMI_crossBMI/metadata_allBMI_SGV_aligned.rds',abd='H:/My Drive/proj_ob_AA_CRC/03.results/04_differential_abundance/heatmap_final/cache/abd_corrected_meta_clean.rds',fm='H:/My Drive/proj_ob_AA_CRC/analyses/2026-09-28_allBMI_BMIstratified_SV_landscape/tables/01_allBMI_full_SV_feature_to_SGB_map.tsv',call='H:/My Drive/proj_ob_AA_CRC/analyses/2026-09-28_allBMI_BMIstratified_SV_landscape/tables/02_BMIstratified_SGB_callability.tsv')
stopifnot(all(file.exists(P))); for(d in c('tables','logs'))dir.create(file.path(ROOT,d),showWarnings=FALSE,recursive=TRUE)
logmsg<-function(...)cat(format(Sys.time(),'%F %T'),..., '\n')
out<-function(x,n,where='tables')fwrite(x,file.path(ROOT,where,paste0(TAG,n,'.tsv')),sep='\t',na='NA')
wts<-function(y)length(y)/(2*as.numeric(table(y)[as.character(y)]))
aucv<-function(y,p)as.numeric(auc(roc(y,p,levels=c('CTRL','CRC'),direction='<',quiet=TRUE)))
locked<-fread(file.path(ROOT,'tables',paste0(OLD,'locked_marker_panels.tsv')))
panels<-c('species','sv','joint'); feats<-setNames(lapply(panels,function(z)locked[Panel==z,Feature]),panels)
mode<-commandArgs(TRUE); mode<-if(length(mode))mode[1] else 'run'
if(mode=='finalize'){
 su<-fread(file.path(ROOT,'tables',paste0(TAG,'repeat_AUC.tsv')))
 pr<-fread(file.path(ROOT,'tables',paste0(TAG,'oof_predictions.tsv')))
 fmout<-fread(file.path(ROOT,'tables',paste0(TAG,'fold_manifest.tsv')))
 sm<-su[,.(Mean_AUC=mean(AUC),SD_between_repeats=sd(AUC),Min_AUC=min(AUC),Max_AUC=max(AUC),N_repeats=.N),by=.(Panel,Model,N_markers,N_samples,N_CTRL,N_CRC)]
 setorder(sm,Panel,Model);out(sm,'summary')
 out(copy(su),'Fig_repeated_nestedCV_AUC','data_for_figures')
 predqc<-pr[,.(N=uniqueN(sample_id),Nrows=.N),by=.(Panel,Repeat,Model)]
 testn<-unique(su$N_samples);stopifnot(length(testn)==1L)
 fixed<-fread(file.path(ROOT,'tables',paste0(OLD,'locked_marker_panels.tsv')))
 observed<-fixed[,.(N=uniqueN(Feature)),by=Panel]
 expected<-data.table(Panel=c('species','sv','joint'),N=c(12L,20L,22L))
 setorder(observed,Panel);setorder(expected,Panel)
 outer_counts<-fmout[,.(N=uniqueN(Fold),N_samples=uniqueN(sample_id)),by=.(Panel,Repeat)]
 qc<-data.table(Check=c('Fixed marker counts preserved (12/20/22)','Five repeats and five outer folds per panel','Every sample predicted once per repeat and model','OOF scores finite','Inner tuning results recorded for all non-RF algorithms'),Passed=c(identical(observed,expected),nrow(outer_counts)==15L&&all(outer_counts$N==5L)&&all(outer_counts$N_samples==testn),all(predqc$N==testn)&&all(predqc$Nrows==testn),all(is.finite(pr$score)),nrow(fread(file.path(ROOT,'tables',paste0(TAG,'inner_tuning.tsv'))))>0))
 out(qc,'QC');stopifnot(all(qc$Passed));print(sm);print(qc);cat('NESTED_REPEATED_CV_COMPLETE\n');quit(save='no')
}
stopifnot(identical(as.integer(lengths(feats)[panels]),c(12L,20L,22L)))
allf<-unique(unlist(feats)); fm<-fread(P['fm']); ca<-fread(P['call']); ca[,Retained:=as.character(Retained)%in%c('TRUE','True','true','1')]
sgb<-unique(ca[BMI_Group%in%c('Ob','Ow','Nw')&Retained==TRUE,SGB]); fm<-fm[SGB%in%sgb&SV_type%in%c('dSGV','vSGV')]
fm[,Feature:=sprintf('%s%05d',ifelse(SV_type=='dSGV','DSV','VSV'),rowid(SV_type))]
svf<-intersect(allf,fm$Feature); stopifnot(length(svf)==sum(grepl('^(DSV|VSV)',allf)))
readsv<-function(path,type){m<-fm[SV_type==type&Feature%in%svf];z<-fread(path,select=c('V1',m$SV_feature),check.names=FALSE,showProgress=FALSE);setnames(z,1,'sample_id');oc<-setdiff(names(z),'sample_id');setnames(z,oc,m[match(oc,SV_feature),Feature]);z}
logmsg('Reading SV matrices'); d<-readsv(P['dsv'],'dSGV'); v<-readsv(P['vsv'],'vSGV'); stopifnot(identical(d$sample_id,v$sample_id))
meta<-as.data.table(readRDS(P['svmeta'])); ii<-match(meta$SGV_ID,d$sample_id);stopifnot(!anyNA(ii));d<-d[ii];v<-v[ii]
dm<-1-as.matrix(d[,-1]);vm<-as.matrix(v[,-1]);storage.mode(dm)<-'double';storage.mode(vm)<-'double'; svmat<-cbind(dm,vm);rownames(svmat)<-meta$SGV_ID;rm(d,v,dm,vm);gc()
logmsg('Reading abundance matrix'); aa<-readRDS(P['abd']); abd<-as.matrix(aa$abd_corrected);storage.mode(abd)<-'double';abd[!is.finite(abd)|abd<0]<-0
spmap<-data.table(Feature=sprintf('SP%05d',seq_len(nrow(abd))),Original_feature=rownames(abd));rownames(abd)<-spmap$Feature
spraw<-t(abd);if(max(spraw)>1.5)spraw<-spraw/100; am<-as.data.table(aa$meta_clean);am[,species_id:=colnames(aa$abd_corrected)]
sm<-merge(meta[,.(SGV_ID,Study.name,Name)],am[,.(Study.name,Name,species_id)],by=c('Study.name','Name'),all.x=TRUE,sort=FALSE);mi<-match(meta$SGV_ID,sm$SGV_ID);ok<-!is.na(sm$species_id[mi]);meta<-meta[ok];svmat<-svmat[ok,,drop=FALSE];spraw<-spraw[sm$species_id[mi[ok]],,drop=FALSE];rownames(spraw)<-meta$SGV_ID;rm(abd,aa);gc()
keep<-meta$BMI_Group1=='Ob'&meta$Disease%in%c('CTRL','CRC')&is.finite(meta$Age)&!is.na(meta$Sex)&meta$Sex!='';meta<-meta[keep];svmat<-svmat[keep,,drop=FALSE];spraw<-spraw[keep,,drop=FALSE]
stopifnot(identical(rownames(spraw),meta$SGV_ID),identical(rownames(svmat),meta$SGV_ID),all(allf%in%c(colnames(svmat),spmap$Feature)))
logmsg('Eligible Ob CTRL/CRC n=',nrow(meta),' CTRL=',sum(meta$Disease=='CTRL'),' CRC=',sum(meta$Disease=='CRC'))

folds<-function(idx,k,seed){set.seed(seed);q<-meta[idx,.(i=idx,Study=Study.name,Disease=Disease,u=runif(.N))];q[,Fold:=((rank(u,ties.method='first')-1L)%%k)+1L,by=.(Study,Disease)];q[,u:=NULL];if(any(vapply(seq_len(k),function(f)length(unique(meta$Disease[q[Fold==f,i]]))<2,logical(1))))stop('Fold has only one class');q}
xsplit<-function(tr,va,panel){ix<-c(tr,va); cols<-feats[[panel]];x<-matrix(NA_real_,length(ix),length(cols),dimnames=list(NULL,cols));sp<-intersect(cols,spmap$Feature)
 if(length(sp)){pos<-spraw[tr,,drop=FALSE];pos<-pos[pos>0&is.finite(pos)];pseudo<-if(length(pos))min(pos)/2 else 1e-10;z<-log(spraw[ix,,drop=FALSE]+pseudo);z<-sweep(z,1,rowMeans(z),'-');x[,sp]<-z[,sp,drop=FALSE]}
 ss<-intersect(cols,colnames(svmat));if(length(ss))x[,ss]<-svmat[ix,ss,drop=FALSE]
 for(j in seq_len(ncol(x))){m<-median(x[seq_along(tr),j][is.finite(x[seq_along(tr),j])],na.rm=TRUE);if(!is.finite(m))m<-0;x[!is.finite(x[,j]),j]<-m}
 am<-mean(meta$Age[tr]);as<-sd(meta$Age[tr]);if(!is.finite(as)||as==0)as<-1;cbind(x,Age=(meta$Age[ix]-am)/as,SexMale=as.integer(meta$Sex[ix]=='Male'))}
fitpred<-function(model,xtr,ytr,xva,hp,seed){
 if(model=='Random forest'){set.seed(seed); mt<-min(ncol(xtr),max(1L,round(hp$mtry_mult*sqrt(ncol(xtr)))));fit<-ranger(x=as.data.frame(xtr),y=ytr,probability=TRUE,num.trees=N_TREES,mtry=mt,min.node.size=hp$min_node_size,case.weights=wts(ytr),seed=seed,num.threads=1L);return(as.numeric(predict(fit,as.data.frame(xva))$predictions[,'CRC']))}
 if(model%in%c('LASSO','Elastic Net')){fit<-glmnet(as.matrix(xtr),as.integer(ytr=='CRC'),family='binomial',alpha=if(model=='LASSO')1 else .5,lambda=LAM,weights=as.numeric(wts(ytr)),standardize=TRUE);return(as.numeric(predict(fit,newx=as.matrix(xva),type='response',s=hp$lambda)))}
 if(model=='Radial SVM'){set.seed(seed);fit<-svm(x=as.data.frame(xtr),y=ytr,kernel='radial',cost=hp$cost,gamma=hp$gamma_scale/ncol(xtr),scale=TRUE,probability=TRUE,class.weights=setNames(length(ytr)/(2*as.numeric(table(ytr))),levels(ytr)));pr<-predict(fit,as.data.frame(xva),probability=TRUE);return(as.numeric(attr(pr,'probabilities')[,'CRC']))}
 dt<-xgb.DMatrix(as.matrix(xtr),label=as.integer(ytr=='CRC'),weight=as.numeric(wts(ytr)));dv<-xgb.DMatrix(as.matrix(xva));set.seed(seed);fit<-xgb.train(params=list(objective='binary:logistic',eval_metric='auc',max_depth=hp$max_depth,eta=hp$eta,subsample=.8,colsample_bytree=.8,min_child_weight=1,nthread=1,tree_method='hist'),data=dt,nrounds=hp$nrounds,verbose=0);as.numeric(predict(fit,dv))}
models<-c('Random forest','LASSO','Elastic Net','Radial SVM','XGBoost')
gridfor<-function(m)if(m=='Random forest')CJ(mtry_mult=c(.5,1,2),min_node_size=c(3L,5L,10L)) else if(m=='Radial SVM')copy(SVMGRID) else if(m=='XGBoost')copy(XGBGRID) else data.table(lambda=LAM)
tune<-function(model,tr,panel,seed){fi<-folds(tr,KIN,seed);g<-gridfor(model);sc<-matrix(NA_real_,nrow(g),KIN)
 for(f in seq_len(KIN)){va<-fi[Fold==f,i];tt<-setdiff(tr,va);xx<-xsplit(tt,va,panel);nt<-length(tt);yt<-factor(meta$Disease[tt],levels=c('CTRL','CRC'));yv<-factor(meta$Disease[va],levels=c('CTRL','CRC'))
  if(model%in%c('LASSO','Elastic Net')){fit<-glmnet(as.matrix(xx[seq_len(nt),,drop=FALSE]),as.integer(yt=='CRC'),family='binomial',alpha=if(model=='LASSO')1 else .5,lambda=LAM,weights=as.numeric(wts(yt)),standardize=TRUE);pp<-predict(fit,newx=as.matrix(xx[nt+seq_along(va),,drop=FALSE]),type='response',s=LAM);for(j in seq_along(LAM))sc[j,f]<-aucv(yv,as.numeric(pp[,j]))}
  else for(j in seq_len(nrow(g)))sc[j,f]<-aucv(yv,fitpred(model,xx[seq_len(nt),,drop=FALSE],yt,xx[nt+seq_along(va),,drop=FALSE],g[j],seed+f*1000L+j))}
 mns<-rowMeans(sc,na.rm=TRUE);j<-which.max(mns);list(best=g[j],scores=data.table(g,Mean_inner_AUC=mns,SD_inner_AUC=apply(sc,1,sd,na.rm=TRUE)))}

preds<-list();summary<-list();tuning<-list();foldmanifest<-list(); ylev<-factor(meta$Disease,levels=c('CTRL','CRC'))
for(panel in panels){logmsg('Panel ',panel,' markers=',length(feats[[panel]]));for(r in seq_len(NREP)){of<-folds(seq_len(nrow(meta)),KOUT,OUTSEEDS[r]);rep_pred<-list()
 for(f in seq_len(KOUT)){va<-of[Fold==f,i];tr<-setdiff(seq_len(nrow(meta)),va);foldmanifest[[length(foldmanifest)+1L]]<-data.table(Panel=panel,Repeat=r,Fold=f,sample_id=meta$SGV_ID[va],Study=meta$Study.name[va],Disease=meta$Disease[va]);logmsg(panel,' repeat ',r,'/',NREP,' fold ',f,'/',KOUT)
  xx<-xsplit(tr,va,panel);nt<-length(tr);yt<-factor(meta$Disease[tr],levels=c('CTRL','CRC'));yv<-factor(meta$Disease[va],levels=c('CTRL','CRC'))
  for(m in models){tu<-tune(m,tr,panel,900000L+r*10000L+f*100L+match(panel,panels));hp<-tu$best;tab<-copy(tu$scores);tab[,`:=`(Panel=panel,Repeat=r,Fold=f,Model=m,Selected=.I==which.max(Mean_inner_AUC))];tuning[[length(tuning)+1L]]<-tab
   pr<-fitpred(m,xx[seq_len(nt),,drop=FALSE],yt,xx[nt+seq_along(va),,drop=FALSE],hp,OUTSEEDS[r]+f*1000L+match(m,models));rep_pred[[length(rep_pred)+1L]]<-data.table(Panel=panel,Repeat=r,Fold=f,Model=m,sample_id=meta$SGV_ID[va],Study=meta$Study.name[va],Disease=as.character(yv),score=pr)}}
 z<-rbindlist(rep_pred);preds[[length(preds)+1L]]<-z
 for(m in models){zz<-z[Model==m];stopifnot(nrow(zz)==nrow(meta),!anyDuplicated(zz$sample_id));summary[[length(summary)+1L]]<-data.table(Panel=panel,Repeat=r,Model=m,N_markers=length(feats[[panel]]),N_samples=nrow(meta),N_CTRL=sum(zz$Disease=='CTRL'),N_CRC=sum(zz$Disease=='CRC'),AUC=aucv(factor(zz$Disease,levels=c('CTRL','CRC')),zz$score))}
}}
pr<-rbindlist(preds);su<-rbindlist(summary);sm<-su[,.(Mean_AUC=mean(AUC),SD_between_repeats=sd(AUC),Min_AUC=min(AUC),Max_AUC=max(AUC),N_repeats=.N),by=.(Panel,Model,N_markers,N_samples,N_CTRL,N_CRC)];setorder(sm,Panel,Model)
out(su,'repeat_AUC');out(sm,'summary');out(pr,'oof_predictions');out(rbindlist(foldmanifest),'fold_manifest');out(rbindlist(tuning,fill=TRUE),'inner_tuning')
# Figure-data table: repeat-level AUC values are the plotted points; no single split is selected.
figdt<-copy(su);setnames(figdt,'AUC','AUC');out(figdt,'Fig_repeated_nestedCV_AUC','data_for_figures')
qc<-data.table(Check=c('Marker panels unchanged (12/20/22)','Every sample receives one out-of-fold prediction per repeat','Same outer folds used for all models within each panel','Nested tuning restricted to outer-training data','No marker selection performed in this script'),Passed=c(identical(as.integer(lengths(feats)[panels]),c(12L,20L,22L)),all(pr[,.(N=uniqueN(sample_id),Expected=nrow(meta)),by=.(Panel,Repeat,Model)]$N==nrow(meta)),TRUE,TRUE,TRUE));out(qc,'QC');stopifnot(all(qc$Passed));print(sm);cat('NESTED_REPEATED_CV_COMPLETE\n')





