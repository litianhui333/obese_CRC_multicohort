#!/usr/bin/env Rscript
# Fixed 22-marker Ob-associated joint panel; four Obese case-control comparisons.
# Marker membership is locked and never reselected. The pipeline mirrors the
# existing 5-repeat nested 5-fold workflow for the fixed panel.
.libPaths(c('D:/R/library', .libPaths()))
suppressPackageStartupMessages({library(data.table); library(pROC); library(glmnet); library(ranger); library(e1071); library(xgboost); library(ggplot2)})

ROOT <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-10-12_Ob_fixed_joint_panel_four_comparisons'
SOURCE <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-10-11_Ob_fixed_panel_xMarkerFinder_like'
TAG <- 'Ob22_fourComp_nested5x5_'
KOUT <- 5L; KIN <- 5L; NREP <- 5L; OUTSEEDS <- 84002L + seq_len(NREP); N_TREES <- 500L
LAM <- 10^seq(-4, 1, length.out=60L)
SVMGRID <- CJ(cost=c(.1,1,10), gamma_scale=c(.1,1,10))
XGBGRID <- CJ(max_depth=c(1L,2L,3L), eta=c(.03,.1), nrounds=c(100L,300L))
P <- c(dsv='H:/My Drive/proj_ob_AA_CRC/02.processed_data/05_SV_result/all_BMI_crossBMI/all_BMI_dsgv.csv',
       vsv='H:/My Drive/proj_ob_AA_CRC/02.processed_data/05_SV_result/all_BMI_crossBMI/all_BMI_vsgv.csv',
       svmeta='H:/My Drive/proj_ob_AA_CRC/02.processed_data/05_SV_result/all_BMI_crossBMI/metadata_allBMI_SGV_aligned.rds',
       abd='H:/My Drive/proj_ob_AA_CRC/03.results/04_differential_abundance/heatmap_final/cache/abd_corrected_meta_clean.rds',
       fm='H:/My Drive/proj_ob_AA_CRC/analyses/2026-09-28_allBMI_BMIstratified_SV_landscape/tables/01_allBMI_full_SV_feature_to_SGB_map.tsv',
       call='H:/My Drive/proj_ob_AA_CRC/analyses/2026-09-28_allBMI_BMIstratified_SV_landscape/tables/02_BMIstratified_SGB_callability.tsv',
       panel=file.path(SOURCE,'tables/restricted_ObAssoc_L1L2_locked_marker_panels.tsv'))
for(d in c('figures','tables','data_for_figures','logs')) dir.create(file.path(ROOT,d), showWarnings=FALSE, recursive=TRUE)
logmsg <- function(...) cat(format(Sys.time(),'%F %T'), ..., '\n')
writeout <- function(x,n,where='tables') fwrite(x,file.path(ROOT,where,paste0(TAG,n,'.tsv')),sep='\t',na='NA')
wts <- function(y) length(y)/(2*as.numeric(table(y)[as.character(y)]))
aucv <- function(y,p) as.numeric(auc(roc(y,p,levels=c('Control','Case'),direction='<',quiet=TRUE)))
models <- c('Random forest','LASSO','Elastic Net','Radial SVM','XGBoost')
comparisons <- data.table(
  Comparison=c('Adenoma vs HC','Adenoma vs CRC','Early CRC vs HC','Late CRC vs HC'),
  Control=c('CTRL','AA','CTRL','CTRL'), Case=c('AA','CRC','CRC_I_II','CRC_III_IV'),
  Definition=c('AA cases; CTRL controls','CRC cases; AA controls','AJCC I/II CRC cases; CTRL controls','AJCC III/IV CRC cases; CTRL controls'))

locked <- fread(P['panel'])
feats <- locked[Panel=='joint', Feature]
stopifnot(length(feats)==22L, uniqueN(feats)==22L)
spfeatures <- feats[grepl('^SP[0-9]+$',feats)]
svfeatures <- feats[grepl('^(DSV|VSV)[0-9]+$',feats)]
stopifnot(length(spfeatures)+length(svfeatures)==22L)
writeout(locked[Panel=='joint'],'fixed_22_marker_panel')

logmsg('Read input matrices and metadata')
fm <- fread(P['fm']); ca <- fread(P['call']); ca[,Retained:=as.character(Retained)%in%c('TRUE','True','true','1')]
sgb <- unique(ca[BMI_Group%in%c('Ob','Ow','Nw') & Retained==TRUE,SGB])
fm <- fm[SGB%in%sgb & SV_type%in%c('dSGV','vSGV')]
fm[,Feature:=sprintf('%s%05d',ifelse(SV_type=='dSGV','DSV','VSV'),rowid(SV_type))]
svfeatures_all <- intersect(feats, fm$Feature)
stopifnot(length(svfeatures_all)==length(svfeatures))
readsv <- function(path,type){
  map <- fm[SV_type==type & Feature%in%svfeatures_all]
  z <- fread(path,select=c('V1',map$SV_feature),check.names=FALSE,showProgress=FALSE)
  setnames(z,1,'sample_id'); old <- setdiff(names(z),'sample_id')
  setnames(z,old,map[match(old,SV_feature),Feature]); z
}
d <- readsv(P['dsv'],'dSGV'); v <- readsv(P['vsv'],'vSGV'); stopifnot(identical(d$sample_id,v$sample_id))
meta <- as.data.table(readRDS(P['svmeta'])); ii <- match(meta$SGV_ID,d$sample_id); stopifnot(!anyNA(ii))
d <- d[ii]; v <- v[ii]
dm <- 1-as.matrix(d[,-1]); vm <- as.matrix(v[,-1]); storage.mode(dm)<-'double'; storage.mode(vm)<-'double'
svmat <- cbind(dm,vm); rownames(svmat)<-meta$SGV_ID; rm(d,v,dm,vm); invisible(gc())
aa <- readRDS(P['abd']); abd <- as.matrix(aa$abd_corrected); storage.mode(abd)<-'double'; abd[!is.finite(abd)|abd<0]<-0
spmap <- data.table(Feature=sprintf('SP%05d',seq_len(nrow(abd))), Original_feature=rownames(abd)); rownames(abd)<-spmap$Feature
spraw <- t(abd); if(max(spraw)>1.5) spraw<-spraw/100
am <- as.data.table(aa$meta_clean); am[,species_id:=colnames(aa$abd_corrected)]
sm <- merge(meta[,.(SGV_ID,Study.name,Name)],am[,.(Study.name,Name,species_id)],by=c('Study.name','Name'),all.x=TRUE,sort=FALSE)
mi <- match(meta$SGV_ID,sm$SGV_ID); ok <- !is.na(sm$species_id[mi])
meta<-meta[ok]; svmat<-svmat[ok,,drop=FALSE]; spraw<-spraw[sm$species_id[mi[ok]],,drop=FALSE]; rownames(spraw)<-meta$SGV_ID
rm(abd,aa); invisible(gc())

stage <- as.character(meta$Tumor.Staging.AJCC)
meta[,OutcomeGroup:=fifelse(Disease=='CTRL','CTRL',fifelse(Disease=='AA','AA',
  fifelse(Disease=='CRC' & stage%in%c('I','II'),'CRC_I_II',
  fifelse(Disease=='CRC' & stage%in%c('III','IV'),'CRC_III_IV',NA_character_))))]
meta_full <- copy(meta); spraw_full <- spraw; svmat_full <- svmat
sample_manifest <- list(); repeat_auc <- list(); oof_all <- list(); fold_manifest <- list(); tuning_all <- list(); roc_all <- list()

for(ci in seq_len(nrow(comparisons))){
  meta <- meta_full; spraw <- spraw_full; svmat <- svmat_full
  comp <- comparisons[ci]; selected <- c(comp$Control,comp$Case)
  if(comp$Comparison=='Adenoma vs HC') selected <- c('CTRL','AA')
  if(comp$Comparison=='Adenoma vs CRC') selected <- c('AA','CRC')
  if(comp$Comparison=='Early CRC vs HC') selected <- c('CTRL','CRC_I_II')
  if(comp$Comparison=='Late CRC vs HC') selected <- c('CTRL','CRC_III_IV')
  if(comp$Comparison=='Adenoma vs CRC'){
    keep <- meta$BMI_Group1=='Ob' & meta$Disease%in%c('AA','CRC') & is.finite(meta$Age) & !is.na(meta$Sex) & meta$Sex!=''
  } else {
    keep <- meta$BMI_Group1=='Ob' & meta$OutcomeGroup%in%selected & is.finite(meta$Age) & !is.na(meta$Sex) & meta$Sex!=''
  }
  m0 <- copy(meta[keep]); m0[,Outcome:=fifelse(OutcomeGroup==selected[2L],'Case','Control')]
  if(comp$Comparison=='Adenoma vs CRC') m0[,Outcome:=fifelse(Disease=='CRC','Case','Control')]
  # Explicitly verify the direction and counts for each contrast.
  stopifnot(setequal(unique(m0$Outcome),c('Control','Case')))
  sample_manifest[[ci]] <- data.table(Comparison=comp$Comparison,Definition=comp$Definition,
    Control_group=selected[1L],Case_group=selected[2L],N_control=sum(m0$Outcome=='Control'),N_case=sum(m0$Outcome=='Case'),
    N_studies=uniqueN(m0$Study.name),N_complete=nrow(m0))
  logmsg(comp$Comparison,': control=',selected[1L],' n=',sum(m0$Outcome=='Control'),'; case=',selected[2L],' n=',sum(m0$Outcome=='Case'))

  keep_idx <- which(keep); sp_i <- spraw[keep_idx,,drop=FALSE]; sv_i <- svmat[keep_idx,,drop=FALSE]
  # Functions below see only this contrast's eligible Ob samples.
  meta <- m0; spraw <- sp_i; svmat <- sv_i
  folds <- function(idx,k,seed){
    set.seed(seed); q<-meta[idx,.(i=idx,Study=Study.name,Outcome=Outcome,u=runif(.N))]
    q[,Fold:=((rank(u,ties.method='first')-1L)%%k)+1L,by=.(Study,Outcome)]; q[,u:=NULL]
    if(any(vapply(seq_len(k),function(f)length(unique(meta$Outcome[q[Fold==f,i]]))<2,logical(1)))) stop('Outer/inner fold lacks one class for ',comp$Comparison)
    q
  }
  xsplit <- function(tr,va){
    ix<-c(tr,va); x<-matrix(NA_real_,length(ix),length(feats),dimnames=list(NULL,feats)); sp<-intersect(feats,spmap$Feature)
    if(length(sp)){pos<-spraw[tr,,drop=FALSE];pos<-pos[pos>0&is.finite(pos)];pseudo<-if(length(pos))min(pos)/2 else 1e-10;z<-log(spraw[ix,,drop=FALSE]+pseudo);z<-sweep(z,1,rowMeans(z),'-');x[,sp]<-z[,sp,drop=FALSE]}
    ss<-intersect(feats,colnames(svmat)); if(length(ss)) x[,ss]<-svmat[ix,ss,drop=FALSE]
    for(j in seq_len(ncol(x))){m<-median(x[seq_along(tr),j][is.finite(x[seq_along(tr),j])],na.rm=TRUE);if(!is.finite(m))m<-0;x[!is.finite(x[,j]),j]<-m}
    am<-mean(meta$Age[tr]);as<-sd(meta$Age[tr]);if(!is.finite(as)||as==0)as<-1
    cbind(x,Age=(meta$Age[ix]-am)/as,SexMale=as.integer(meta$Sex[ix]=='Male'))
  }
  fitpred <- function(model,xtr,ytr,xva,hp,seed){
    if(model=='Random forest'){set.seed(seed);mt<-min(ncol(xtr),max(1L,round(hp$mtry_mult*sqrt(ncol(xtr)))));fit<-ranger(x=as.data.frame(xtr),y=ytr,probability=TRUE,num.trees=N_TREES,mtry=mt,min.node.size=hp$min_node_size,case.weights=wts(ytr),seed=seed,num.threads=1L);return(as.numeric(predict(fit,as.data.frame(xva))$predictions[,'Case']))}
    if(model%in%c('LASSO','Elastic Net')){fit<-glmnet(as.matrix(xtr),as.integer(ytr=='Case'),family='binomial',alpha=if(model=='LASSO')1 else .5,lambda=LAM,weights=as.numeric(wts(ytr)),standardize=TRUE);return(as.numeric(predict(fit,newx=as.matrix(xva),type='response',s=hp$lambda)))}
    if(model=='Radial SVM'){set.seed(seed);fit<-svm(x=as.data.frame(xtr),y=ytr,kernel='radial',cost=hp$cost,gamma=hp$gamma_scale/ncol(xtr),scale=TRUE,probability=TRUE,class.weights=setNames(length(ytr)/(2*as.numeric(table(ytr))),levels(ytr)));pr<-predict(fit,as.data.frame(xva),probability=TRUE);return(as.numeric(attr(pr,'probabilities')[,'Case']))}
    dt<-xgb.DMatrix(as.matrix(xtr),label=as.integer(ytr=='Case'),weight=as.numeric(wts(ytr)));dv<-xgb.DMatrix(as.matrix(xva));set.seed(seed);fit<-xgb.train(params=list(objective='binary:logistic',eval_metric='auc',max_depth=hp$max_depth,eta=hp$eta,subsample=.8,colsample_bytree=.8,min_child_weight=1,nthread=1,tree_method='hist'),data=dt,nrounds=hp$nrounds,verbose=0);as.numeric(predict(fit,dv))
  }
  gridfor<-function(m)if(m=='Random forest')CJ(mtry_mult=c(.5,1,2),min_node_size=c(3L,5L,10L)) else if(m=='Radial SVM')copy(SVMGRID) else if(m=='XGBoost')copy(XGBGRID) else data.table(lambda=LAM)
  tune<-function(model,tr,seed){
    fi<-folds(tr,KIN,seed);g<-gridfor(model);sc<-matrix(NA_real_,nrow(g),KIN)
    for(f in seq_len(KIN)){va<-fi[Fold==f,i];tt<-setdiff(tr,va);xx<-xsplit(tt,va);nt<-length(tt);yt<-factor(meta$Outcome[tt],levels=c('Control','Case'));yv<-factor(meta$Outcome[va],levels=c('Control','Case'))
      if(model%in%c('LASSO','Elastic Net')){fit<-glmnet(as.matrix(xx[seq_len(nt),,drop=FALSE]),as.integer(yt=='Case'),family='binomial',alpha=if(model=='LASSO')1 else .5,lambda=LAM,weights=as.numeric(wts(yt)),standardize=TRUE);pp<-predict(fit,newx=as.matrix(xx[nt+seq_along(va),,drop=FALSE]),type='response',s=LAM);for(j in seq_along(LAM))sc[j,f]<-aucv(yv,as.numeric(pp[,j]))}
      else for(j in seq_len(nrow(g)))sc[j,f]<-aucv(yv,fitpred(model,xx[seq_len(nt),,drop=FALSE],yt,xx[nt+seq_along(va),,drop=FALSE],g[j],seed+f*1000L+j))
    }
    means<-rowMeans(sc,na.rm=TRUE);j<-which.max(means);list(best=g[j],scores=data.table(g,Mean_inner_AUC=means,SD_inner_AUC=apply(sc,1,sd,na.rm=TRUE)))
  }
  n<-nrow(meta); preds_comp<-list()
  for(r in seq_len(NREP)){
    of<-folds(seq_len(n),KOUT,OUTSEEDS[r]); pred_rep<-list()
    for(f in seq_len(KOUT)){
      va<-of[Fold==f,i];tr<-setdiff(seq_len(n),va); fold_manifest[[length(fold_manifest)+1L]]<-data.table(Comparison=comp$Comparison,Repeat=r,Fold=f,sample_id=meta$SGV_ID[va],Outcome=meta$Outcome[va],Study=meta$Study.name[va])
      logmsg(comp$Comparison,' repeat ',r,'/',NREP,' fold ',f,'/',KOUT)
      xx<-xsplit(tr,va);nt<-length(tr);yt<-factor(meta$Outcome[tr],levels=c('Control','Case'))
      for(mod in models){tu<-tune(mod,tr,900000L+ci*100000L+r*10000L+f*100L);hp<-tu$best
        tab<-copy(tu$scores);tab[,`:=`(Comparison=comp$Comparison,Repeat=r,Fold=f,Model=mod,Selected=.I==which.max(Mean_inner_AUC))];tuning_all[[length(tuning_all)+1L]]<-tab
        pr<-fitpred(mod,xx[seq_len(nt),,drop=FALSE],yt,xx[nt+seq_along(va),,drop=FALSE],hp,OUTSEEDS[r]+f*1000L+match(mod,models))
        pred_rep[[length(pred_rep)+1L]]<-data.table(Comparison=comp$Comparison,Repeat=r,Fold=f,Model=mod,sample_id=meta$SGV_ID[va],Study=meta$Study.name[va],Outcome=as.character(meta$Outcome[va]),score=pr)
      }
    }
    z<-rbindlist(pred_rep);preds_comp[[r]]<-z
    for(mod in models){zz<-z[Model==mod];stopifnot(nrow(zz)==n,uniqueN(zz$sample_id)==n,all(is.finite(zz$score)))
      repeat_auc[[length(repeat_auc)+1L]]<-data.table(Comparison=comp$Comparison,Repeat=r,Model=mod,N_markers=22L,N_samples=n,N_control=sum(zz$Outcome=='Control'),N_case=sum(zz$Outcome=='Case'),AUC=aucv(factor(zz$Outcome,levels=c('Control','Case')),zz$score))
    }
  }
  pp<-rbindlist(preds_comp);oof_all[[ci]]<-pp
  auccomp<-rbindlist(repeat_auc)[Comparison==comp$Comparison]
  for(mod in models){
    aa<-auccomp[Model==mod]; rep_pick<-aa[which.min(abs(AUC-mean(AUC))),Repeat]
    zz<-pp[Model==mod & Repeat==rep_pick]
    rocobj<-roc(factor(zz$Outcome,levels=c('Control','Case')),zz$score,levels=c('Control','Case'),direction='<',quiet=TRUE)
    rr<-data.table(FPR=1-as.numeric(rocobj$specificities),TPR=as.numeric(rocobj$sensitivities))
    rr[,`:=`(Comparison=comp$Comparison,Model=mod,Representative_repeat=rep_pick,Representative_AUC=aa[Repeat==rep_pick,AUC],Mean_AUC=mean(aa$AUC),SD_AUC=sd(aa$AUC))]
    roc_all[[length(roc_all)+1L]]<-rr
  }
}

repeat_dt<-rbindlist(repeat_auc); pred_dt<-rbindlist(oof_all); roc_dt<-rbindlist(roc_all); folds_dt<-rbindlist(fold_manifest)
summary_dt<-repeat_dt[,.(Mean_AUC=mean(AUC),SD_between_repeats=sd(AUC),Min_AUC=min(AUC),Max_AUC=max(AUC),N_repeats=.N),by=.(Comparison,Model,N_markers,N_samples,N_control,N_case)]
setorder(summary_dt,Comparison,Model)
writeout(rbindlist(sample_manifest),'sample_counts');writeout(repeat_dt,'repeat_AUC');writeout(summary_dt,'summary');writeout(pred_dt,'oof_predictions');writeout(folds_dt,'fold_manifest');writeout(rbindlist(tuning_all,fill=TRUE),'inner_tuning');writeout(roc_dt,'Fig_representative_ROC_curves','data_for_figures')

for(comp in comparisons$Comparison){
  dd<-roc_dt[Comparison==comp]; labels<-summary_dt[Comparison==comp]
  p<-ggplot(dd,aes(FPR,TPR,color=Model,group=Model))+geom_step(linewidth=.8,direction='vh')+
    geom_abline(slope=1,intercept=0,linetype='dashed',color='grey75')+
    coord_equal(xlim=c(0,1),ylim=c(0,1),expand=FALSE)+labs(title=comp,x='False-positive rate',y='True-positive rate',color=NULL)+
    theme_classic(base_size=12)+theme(legend.position='bottom',plot.title=element_text(face='bold',hjust=.5))+
    geom_text(data=labels,aes(x=.98,y=seq(.10,.10+0.07*(nrow(labels)-1),length.out=nrow(labels)),label=sprintf('%s  %.3f ± %.3f',Model,Mean_AUC,SD_between_repeats)),inherit.aes=FALSE,hjust=1,size=3.0,color='#333333')
  stem<-gsub('[^A-Za-z0-9]+','_',comp)
  ggsave(file.path(ROOT,'figures',paste0(TAG,'ROC_',stem,'.pdf')),p,width=5.7,height=4.8,device=cairo_pdf)
  ggsave(file.path(ROOT,'figures',paste0(TAG,'ROC_',stem,'.png')),p,width=5.7,height=4.8,dpi=300)
}

qc<-data.table(Check=c('Exactly 22 unique locked joint markers','Five repeats per comparison/model','Five outer folds per repeat','Five inner folds used for tuning','One OOF prediction per eligible sample/repeat/model','No LODO performed','Early/late CRC use AJCC I/II and III/IV'),Passed=c(length(feats)==22L&&uniqueN(feats)==22L,all(repeat_dt[,.(n=.N),by=.(Comparison,Model)]$n==NREP),all(folds_dt[,.(n=uniqueN(Fold)),by=.(Comparison,Repeat)]$n==KOUT),TRUE,all(pred_dt[,.(n=.N,u=uniqueN(sample_id)),by=.(Comparison,Repeat,Model)]$n==pred_dt[,.(n=.N,u=uniqueN(sample_id)),by=.(Comparison,Repeat,Model)]$u),TRUE,TRUE))
pred_check<-pred_dt[,.(n=.N,u=uniqueN(sample_id)),by=.(Comparison,Repeat,Model)]
qc[Check=='One OOF prediction per eligible sample/repeat/model',Passed:=all(pred_check$n==pred_check$u)]
writeout(qc,'QC');stopifnot(all(qc$Passed));logmsg('COMPLETE; outputs written to ',ROOT)
