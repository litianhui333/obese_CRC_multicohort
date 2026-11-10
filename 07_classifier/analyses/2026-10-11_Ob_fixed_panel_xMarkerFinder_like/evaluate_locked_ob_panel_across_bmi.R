#!/usr/bin/env Rscript
.libPaths(c('D:/R/library', .libPaths()))
suppressPackageStartupMessages({
  library(data.table); library(ranger); library(pROC); library(ggplot2)
})

ROOT <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-10-11_Ob_fixed_panel_xMarkerFinder_like'
PREFIX <- 'restricted_ObAssoc_L1L2_'
PATHS <- c(
 dsv='H:/My Drive/proj_ob_AA_CRC/02.processed_data/05_SV_result/all_BMI_crossBMI/all_BMI_dsgv.csv',
 vsv='H:/My Drive/proj_ob_AA_CRC/02.processed_data/05_SV_result/all_BMI_crossBMI/all_BMI_vsgv.csv',
 sv_meta='H:/My Drive/proj_ob_AA_CRC/02.processed_data/05_SV_result/all_BMI_crossBMI/metadata_allBMI_SGV_aligned.rds',
 abundance='H:/My Drive/proj_ob_AA_CRC/03.results/04_differential_abundance/heatmap_final/cache/abd_corrected_meta_clean.rds',
 feature_map='H:/My Drive/proj_ob_AA_CRC/analyses/2026-09-28_allBMI_BMIstratified_SV_landscape/tables/01_allBMI_full_SV_feature_to_SGB_map.tsv',
 callability='H:/My Drive/proj_ob_AA_CRC/analyses/2026-09-28_allBMI_BMIstratified_SV_landscape/tables/02_BMIstratified_SGB_callability.tsv',
 markers=file.path(ROOT,'tables',paste0(PREFIX,'locked_marker_panels.tsv')),
 split=file.path(ROOT,'tables',paste0(PREFIX,'fixed_development_test_split.tsv')),
 expected=file.path(ROOT,'tables',paste0(PREFIX,'locked_panel_test_AUC.tsv')))
stopifnot(all(file.exists(PATHS)))

N_TREES <- 500L
PANELS <- c(species='Species only',sv='SV only',joint='Species + SV')
FINAL_SEEDS <- setNames(61001L + seq_along(PANELS),names(PANELS))
wt <- function(x,n,w='tables') fwrite(x,file.path(ROOT,w,paste0(PREFIX,n,'.tsv')),sep='\t',na='NA')
wgt <- function(y) length(y)/(2*as.numeric(table(y)[as.character(y)]))

markers <- fread(PATHS['markers'])
split <- fread(PATHS['split'])
expected <- fread(PATHS['expected'])
stopifnot(all(names(PANELS)%in%unique(markers$Panel)))

# Reconstruct the exact internal feature identifiers used by the locked analysis.
fm <- fread(PATHS['feature_map']); ca <- fread(PATHS['callability'])
ca[,Retained:=as.character(Retained)%in%c('TRUE','True','true','1')]
sgb_union <- unique(ca[BMI_Group%in%c('Ob','Ow','Nw')&Retained==TRUE,SGB])
fm <- fm[SGB%in%sgb_union&SV_type%in%c('dSGV','vSGV')]
fm[,Feature:=sprintf('%s%05d',ifelse(SV_type=='dSGV','DSV','VSV'),rowid(SV_type))]
sv_need <- unique(markers[Type=='SV',Feature])

read_sv_selected <- function(path,type){
 mp <- fm[SV_type==type&Feature%in%sv_need]
 z <- fread(path,select=c('V1',mp$SV_feature),check.names=FALSE,showProgress=FALSE)
 setnames(z,1L,'sample_id'); old <- setdiff(names(z),'sample_id')
 setnames(z,old,mp[match(old,SV_feature),Feature]); z
}
dd <- read_sv_selected(PATHS['dsv'],'dSGV')
vv <- read_sv_selected(PATHS['vsv'],'vSGV')
stopifnot(identical(dd$sample_id,vv$sample_id))

meta <- as.data.table(readRDS(PATHS['sv_meta']))
ii <- which(dd$sample_id%in%meta$SGV_ID)
stopifnot(length(ii)==nrow(meta))
dd <- dd[ii]; vv <- vv[ii]; sid <- dd$sample_id
dmat <- 1-as.matrix(dd[,-1]); storage.mode(dmat)<-'double'
vmat <- as.matrix(vv[,-1]); storage.mode(vmat)<-'double'
svmat <- cbind(dmat,vmat); rownames(svmat)<-sid
meta <- meta[match(sid,SGV_ID)]

# Recreate the original all-species CLR transformation; the pseudocount is learned
# only from Ob development samples below.
aa <- readRDS(PATHS['abundance']); abd <- as.matrix(aa$abd_corrected)
storage.mode(abd)<-'double'; abd[!is.finite(abd)|abd<0]<-0
spmap <- data.table(Feature=sprintf('SP%05d',seq_len(nrow(abd))),Original_feature=rownames(abd))
rownames(abd)<-spmap$Feature; species_raw<-t(abd)
if(max(species_raw)>1.5) species_raw<-species_raw/100
ameta<-as.data.table(aa$meta_clean); ameta[,species_id:=colnames(aa$abd_corrected)]
smap<-merge(meta[,.(SGV_ID,Study.name,Name)],ameta[,.(Study.name,Name,species_id)],
 by=c('Study.name','Name'),all.x=TRUE,sort=FALSE)
mi<-match(meta$SGV_ID,smap$SGV_ID); ok<-!is.na(smap$species_id[mi])
meta<-meta[ok]; svmat<-svmat[ok,,drop=FALSE]
species_raw<-species_raw[smap$species_id[mi[ok]],,drop=FALSE]
rownames(species_raw)<-meta$SGV_ID

eligible <- meta$BMI_Group1%in%c('Ob','Ow','Normal') & meta$Disease%in%c('CTRL','CRC') &
 is.finite(meta$Age)&!is.na(meta$Sex)&meta$Sex!=''
meta<-meta[eligible]; svmat<-svmat[eligible,,drop=FALSE]; species_raw<-species_raw[eligible,,drop=FALSE]
rownames(svmat)<-meta$SGV_ID; rownames(species_raw)<-meta$SGV_ID

dev_ids <- split[split=='development',sample_id]
test_ids <- split[split=='test',sample_id]
dev <- match(dev_ids,meta$SGV_ID); ob_test <- match(test_ids,meta$SGV_ID)
stopifnot(!anyNA(dev),!anyNA(ob_test),all(meta$BMI_Group1[dev]=='Ob'),all(meta$BMI_Group1[ob_test]=='Ob'))

pos <- species_raw[dev,][species_raw[dev,]>0]
pseudo <- min(pos)/2
zz <- log(species_raw+pseudo); zz <- sweep(zz,1,rowMeans(zz),'-')

fit_rf <- function(xtr,ytr,xte,seed){
 set.seed(seed)
 fit<-ranger(x=as.data.frame(xtr),y=ytr,probability=TRUE,num.trees=N_TREES,
  mtry=max(1L,floor(sqrt(ncol(xtr)))),min.node.size=5L,case.weights=wgt(ytr),
  importance='permutation',seed=seed,num.threads=1L)
 as.numeric(predict(fit,as.data.frame(xte))$predictions[,'CRC'])
}

group_indices <- list(Ob_test=ob_test,
 Nw=which(meta$BMI_Group1=='Normal'),Ow=which(meta$BMI_Group1=='Ow'))
pred_out<-list(); auc_out<-list(); avail_out<-list(); roc_out<-list(); tests<-list()

for(panel in names(PANELS)){
 ff <- markers[Panel==panel,Feature]
 sp <- intersect(ff,colnames(zz)); sv <- intersect(ff,colnames(svmat))
 stopifnot(length(c(sp,sv))==length(ff))
 x <- cbind(zz[,sp,drop=FALSE],svmat[,sv,drop=FALSE])[,ff,drop=FALSE]
 med <- vapply(seq_len(ncol(x)),function(j)median(x[dev,j],na.rm=TRUE),numeric(1))
 stopifnot(all(is.finite(med)))
 for(j in seq_len(ncol(x))) x[!is.finite(x[,j]),j] <- med[j]
 x <- cbind(x,Age=meta$Age,SexMale=as.integer(meta$Sex=='Male'))
 ydev <- factor(meta$Disease[dev],levels=c('CTRL','CRC'))
 fitscore <- fit_rf(x[dev,,drop=FALSE],ydev,x,FINAL_SEEDS[panel])
 panel_rocs<-list()
 for(g in names(group_indices)){
  ix<-group_indices[[g]]; yy<-factor(meta$Disease[ix],levels=c('CTRL','CRC'))
  stopifnot(length(unique(yy))==2)
  ro<-roc(yy,fitscore[ix],levels=c('CTRL','CRC'),direction='<',quiet=TRUE)
  ci<-as.numeric(ci.auc(ro,method='delong'))
  auc_out[[length(auc_out)+1L]]<-data.table(Panel=panel,Panel_label=PANELS[panel],BMI_evaluation=g,
   N=length(ix),N_CTRL=sum(yy=='CTRL'),N_CRC=sum(yy=='CRC'),AUC=as.numeric(auc(ro)),CI_low=ci[1],CI_high=ci[3])
  pred_out[[length(pred_out)+1L]]<-data.table(Panel=panel,BMI_evaluation=g,sample_id=meta$SGV_ID[ix],
   Study=meta$Study.name[ix],Disease=meta$Disease[ix],score=fitscore[ix])
  rp<-data.table(Panel=panel,BMI_evaluation=g,FPR=1-ro$specificities,TPR=ro$sensitivities,
   Threshold=ro$thresholds); setorder(rp,FPR,TPR); roc_out[[length(roc_out)+1L]]<-rp
  panel_rocs[[g]]<-ro
  if(length(sv)) avail_out[[length(avail_out)+1L]]<-data.table(Panel=panel,BMI_evaluation=g,
   Feature=sv,Observed_fraction=vapply(seq_along(sv),function(j)mean(is.finite(svmat[ix,sv[j]])),numeric(1)))
 }
 for(g in c('Nw','Ow')){
  tt<-roc.test(panel_rocs[['Ob_test']],panel_rocs[[g]],paired=FALSE,method='delong',alternative='greater')
  tests[[length(tests)+1L]]<-data.table(Panel=panel,Contrast=paste0('Ob_test > ',g),
   P_one_sided=as.numeric(tt$p.value))
 }
}

auc_dt<-rbindlist(auc_out); pred_dt<-rbindlist(pred_out); roc_dt<-rbindlist(roc_out)
avail_dt<-rbindlist(avail_out); test_dt<-rbindlist(tests)

# Exact reproducibility check against the previously saved Ob test AUCs.
chk<-merge(auc_dt[BMI_evaluation=='Ob_test',.(Panel,AUC)],expected[,.(Panel,Expected=Test_AUC)],by='Panel')
stopifnot(all(abs(chk$AUC-chk$Expected)<1e-12))

wt(auc_dt,'crossBMI_transfer_AUC')
wt(pred_dt,'crossBMI_transfer_predictions')
wt(avail_dt,'crossBMI_transfer_SV_availability')
wt(test_dt,'crossBMI_transfer_AUC_tests')

plot_dt<-merge(roc_dt,auc_dt[,.(Panel,BMI_evaluation,AUC,N)],by=c('Panel','BMI_evaluation'))
plot_dt[,BMI_label:=factor(BMI_evaluation,levels=c('Ob_test','Nw','Ow'),
 labels=c('Ob test','Nw','Ow'))]
plot_dt[,Panel_label:=factor(PANELS[Panel],levels=PANELS)]
plot_dt[,Curve_label:=sprintf('%s (n=%d; AUC %.3f)',BMI_label,N,AUC)]
wt(plot_dt,'Fig_crossBMI_transfer_ROC','data_for_figures')
cols<-c('Ob test'='#C0523D','Nw'='#477DB3','Ow'='#E5A052')
ann<-copy(auc_dt)
ann[,BMI_label:=factor(BMI_evaluation,levels=c('Ob_test','Nw','Ow'),labels=c('Ob test','Nw','Ow'))]
ann[,Panel_label:=factor(PANELS[Panel],levels=PANELS)]
ann[,y:=c('Ob test'=.20,'Nw'=.13,'Ow'=.06)[as.character(BMI_label)]]
ann[,Label:=sprintf('%s: %.3f',BMI_label,AUC)]
g<-ggplot(plot_dt,aes(FPR,TPR,colour=BMI_label))+geom_abline(slope=1,intercept=0,
 linetype='dashed',colour='#BBBBBB',linewidth=.45)+geom_step(linewidth=.8,direction='vh')+
 geom_text(data=ann,aes(x=.98,y=y,label=Label,colour=BMI_label),inherit.aes=FALSE,
  hjust=1,size=3.2,fontface='bold')+
 scale_colour_manual(values=cols,drop=FALSE)+facet_wrap(~Panel_label,nrow=1)+
 coord_equal(xlim=c(0,1),ylim=c(0,1),expand=FALSE)+
 labs(x='False-positive rate',y='True-positive rate',title='Transfer of Ob-trained fixed panels')+
 theme_classic(base_size=12,base_family='Arial')+theme(strip.background=element_rect(fill='#F3F3F3',colour=NA),
 strip.text=element_text(face='bold'),plot.title=element_text(face='bold'),legend.position='none')
grDevices::cairo_pdf(file.path(ROOT,'figures',paste0(PREFIX,'Fig_crossBMI_transfer_ROC.pdf')),
 width=8.4,height=3.8,family='Arial');print(g);dev.off()
ggsave(file.path(ROOT,'figures',paste0(PREFIX,'Fig_crossBMI_transfer_ROC.png')),g,width=8.4,height=3.8,dpi=300)
if(requireNamespace('svglite',quietly=TRUE)){
 svglite::svglite(file.path(ROOT,'figures',paste0(PREFIX,'Fig_crossBMI_transfer_ROC.svg')),width=8.4,height=3.8)
 print(g);dev.off()
}
cat('CROSS_BMI_TRANSFER_COMPLETED\n')
