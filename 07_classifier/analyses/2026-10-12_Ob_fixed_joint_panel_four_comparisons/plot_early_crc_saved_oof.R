#!/usr/bin/env Rscript
# Target: Early CRC only, saved CRC-CTRL OOF predictions, no fitting or tuning.
.libPaths(c('D:/R/library', .libPaths()))
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
ROOT <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-10-12_Ob_fixed_joint_panel_four_comparisons'
SOURCE_ROOT <- file.path(dirname(ROOT),'2026-10-11_Ob_fixed_panel_xMarkerFinder_like')
STEM <- 'Fig_Ob_EarlyCRC_originalModel_OOF_withinStudy_noLASSO'
MODELS <- c('Random forest','XGBoost','Radial SVM','Elastic Net')
COLORS <- c('Random forest'='#C34E40','XGBoost'='#D99A30','Radial SVM'='#8062A4','Elastic Net'='#72A090')
INPUTS <- c(predictions=file.path(SOURCE_ROOT,'tables/restricted_ObAssoc_L1L2_nestedCV_oof_predictions.tsv'),
 metadata='H:/My Drive/proj_ob_AA_CRC/02.processed_data/05_SV_result/all_BMI_crossBMI/metadata_allBMI_SGV_aligned.rds',
 audit=file.path(ROOT,'tables/Ob22_CTRLCRC_stage_transfer_nested5x5_repeat_AUC.tsv'))
existing <- list.files(file.path(ROOT,'figures'),full.names=TRUE)
before <- file.info(existing[!startsWith(basename(existing),STEM)])
pred <- fread(INPUTS['predictions'],colClasses=list(character='sample_id'))[Panel=='joint' & Model %in% MODELS]
meta <- as.data.table(readRDS(INPUTS['metadata']))
stopifnot(!anyDuplicated(meta$SGV_ID),!anyDuplicated(pred[,.(Repeat,Model,sample_id)]),setequal(pred$Repeat,1:5))
m <- match(pred$sample_id,as.character(meta$SGV_ID))
stopifnot(!anyNA(m),all(pred$Disease==meta$Disease[m]),all(pred$Study==meta$Study.name[m]))
pred[,AJCC_stage:=as.character(meta$Tumor.Staging.AJCC[m])]
# Keep every CTRL and only CRC with recorded AJCC I/II; other CRC stages are out of scope.
z <- pred[Disease=='CTRL' | (Disease=='CRC' & AJCC_stage %in% c('I','II'))]
z[,Role:=ifelse(Disease=='CTRL','negative','positive')]
stopifnot(z[,all(.N==336L),by=.(Model,Repeat)]$V1,
          z[,all(sum(Role=='positive')==103L & sum(Role=='negative')==233L),by=.(Model,Repeat)]$V1)
audit <- fread(INPUTS['audit'])[Comparison=='Early CRC vs HC' & Model %in% MODELS]
# Copy only pure ROC geometry functions from the targeted plotting script.
helper_file <- file.path(ROOT,'plot_within_study_transfer_four_panels.R')
helper_ast <- parse(helper_file)
for(n in c('local_roc','trap_auc','roc_at')) {
  selected <- Filter(function(x) is.call(x) && identical(x[[1L]],as.name('<-')) && identical(x[[2L]],as.name(n)),as.list(helper_ast))
  stopifnot(length(selected)==1L); eval(selected[[1L]])
}
z[,Plot_score:=score]
counts <- z[,.(N_negative=sum(Role=='negative'),N_positive=sum(Role=='positive')),by=.(Model,Repeat,Study)]
counts[,`:=`(Pairs=N_negative*N_positive,Eligible=N_negative>0 & N_positive>0)]
counts[,Weight:=Pairs/sum(Pairs),by=.(Model,Repeat)]
z <- merge(z,counts[,.(Model,Repeat,Study,Eligible)],by=c('Model','Repeat','Study'),all.x=TRUE,sort=FALSE)
local_list <- list(); repeat_list <- list(); output_curves <- list(); labels <- list()
for(method in MODELS) {
  local_curves <- list(); local_weights <- numeric()
  method_ct <- counts[Model==method & Eligible==TRUE]
  for(repeat_id in 1:5) {
    ct <- method_ct[Repeat==repeat_id]
    local_auc <- numeric(nrow(ct))
    for(k in seq_len(nrow(ct))) {
      d <- z[Model==method & Repeat==repeat_id & Study==ct$Study[k]]
      q <- local_roc(d)
      local_auc[k] <- trap_auc(q)
      local_curves[[length(local_curves)+1L]] <- q
      local_weights <- c(local_weights,ct$Weight[k]/5)
      loc <- copy(q)
      loc[,`:=`(Model=method,Repeat=repeat_id,Study=ct$Study[k],Weight=ct$Weight[k],Point_order=.I)]
      local_list[[length(local_list)+1L]] <- loc
    }
    auc <- sum(ct$Weight*local_auc)
    ar <- audit[Model==method & Repeat==repeat_id]
    stopifnot(nrow(ar)==1L,abs(auc-ar$AUC_same_study_pairs)<1e-12,sum(ct$Pairs)==ar$N_same_study_pairs)
    repeat_list[[length(repeat_list)+1L]] <- data.table(Model=method,Repeat=repeat_id,AUC=auc,
      Audit_AUC=ar$AUC_same_study_pairs,Absolute_error=abs(auc-ar$AUC_same_study_pairs))
  }
  grid <- sort(unique(unlist(lapply(local_curves,function(q) q$FPR))))
  mean_y <- matrix(0,nrow=length(grid),ncol=2)
  for(k in seq_along(local_curves)) mean_y <- mean_y+local_weights[k]*roc_at(local_curves[[k]],grid)
  curve <- data.table(FPR=rep(grid,each=2L),TPR=as.vector(t(mean_y)))
  aucs <- rbindlist(repeat_list)[Model==method,AUC]
  stopifnot(length(aucs)==5L,abs(trap_auc(curve)-mean(aucs))<1e-12)
  curve[,`:=`(Model=method,Point_order=.I,Mean_AUC=mean(aucs),SD_AUC=sd(aucs))]
  output_curves[[length(output_curves)+1L]] <- curve
  used <- z[Model==method & Repeat==1 & Eligible==TRUE]
  labels[[length(labels)+1L]] <- data.table(Model=method,Mean_AUC=mean(aucs),SD_AUC=sd(aucs),
    N_early_all=103L,N_CTRL_all=233L,N_early_used=sum(used$Role=='positive'),N_CTRL_used=sum(used$Role=='negative'),
    Studies=method_ct[Repeat==1,.N],Pairs=sum(method_ct[Repeat==1,Pairs]),
    Label=sprintf('%s %.3f',method,mean(aucs)),Label_x=.98,Label_y=.30-(match(method,MODELS)-1)*.065)
}
roc <- rbindlist(output_curves); ann <- rbindlist(labels); repeats <- rbindlist(repeat_list)
roc[,Model:=factor(Model,levels=MODELS)]; ann[,Model:=factor(Model,levels=MODELS)]
setorder(roc,Model,Point_order); setorder(ann,Model)
stopifnot(nrow(ann)==4L,nrow(repeats)==20L,!any(as.character(roc$Model)=='LASSO'))
for(folder in c('figures','tables','data_for_figures')) dir.create(file.path(ROOT,folder),showWarnings=FALSE)
fwrite(roc,file.path(ROOT,'data_for_figures',paste0(STEM,'_curves.tsv')),sep='\t')
fwrite(ann,file.path(ROOT,'data_for_figures',paste0(STEM,'_labels.tsv')),sep='\t')
fwrite(z,file.path(ROOT,'data_for_figures',paste0(STEM,'_samples.tsv')),sep='\t',na='NA')
fwrite(rbindlist(local_list),file.path(ROOT,'data_for_figures',paste0(STEM,'_study_repeat_curves.tsv')),sep='\t')
fwrite(repeats,file.path(ROOT,'tables',paste0(STEM,'_repeat_AUC_QC.tsv')),sep='\t')
fwrite(counts,file.path(ROOT,'tables',paste0(STEM,'_study_weights.tsv')),sep='\t')
fwrite(ann,file.path(ROOT,'tables',paste0(STEM,'_summary.tsv')),sep='\t')
nn <- ann[1L]
fig <- ggplot(roc,aes(FPR,TPR,colour=Model,group=Model))+
 geom_abline(slope=1,intercept=0,linetype='dashed',colour='#C5C5C5',linewidth=.4)+
 geom_path(linewidth=.75)+
 geom_text(data=ann,aes(x=Label_x,y=Label_y,label=Label),hjust=1,fontface='bold',size=3.2,show.legend=FALSE)+
 scale_colour_manual(values=COLORS)+
 scale_x_continuous(breaks=c(0,.25,.5,.75,1),expand=expansion(mult=.005))+
 scale_y_continuous(breaks=c(0,.25,.5,.75,1),expand=expansion(mult=.005))+
 coord_fixed(xlim=c(0,1),ylim=c(0,1),clip='off')+
 labs(title='Obese group\nEarly CRC (Stage I-II) vs CTRL',
 subtitle=sprintf('Original CRC-CTRL out-of-fold predictions\n%d early CRC, %d CTRL; %d within-study pairs',nn$N_early_used,nn$N_CTRL_used,nn$Pairs),
 x='Within-study false-positive rate',y='Within-study true-positive rate',
 caption='Mean ROC and AUC over five CV repeats; no stage-specific refitting.')+
 theme_classic(base_size=12,base_family='Arial')+
 theme(legend.position='none',plot.title=element_text(face='bold',size=12,hjust=.5),
 plot.subtitle=element_text(size=9,hjust=.5,margin=margin(b=6)),axis.text=element_text(size=10),
 axis.title=element_text(size=11),plot.caption=element_text(size=8,hjust=.5),plot.margin=margin(10,12,8,10))
base <- file.path(ROOT,'figures',STEM)
grDevices::cairo_pdf(paste0(base,'.pdf'),width=105/25.4,height=121/25.4,family='Arial'); print(fig); dev.off()
svglite::svglite(paste0(base,'.svg'),width=105/25.4,height=121/25.4); print(fig); dev.off()
ggsave(paste0(base,'.png'),fig,width=105,height=121,units='mm',device=ragg::agg_png,dpi=300,bg='white')
after <- file.info(rownames(before)); stopifnot(identical(before$size,after$size),identical(before$mtime,after$mtime))
legend <- c('Early CRC subgroup evaluation using saved out-of-fold predictions of the original CRC-CTRL classification workflow.',
'Four classifiers use the fixed 22-marker species-SV panel with age and sex. Predictions originate from five repeats of nested five-fold CRC-CTRL cross-validation; this targeted script performs no fitting, tuning or stage-specific model optimization. This is not application of a single frozen estimator.',
'Of 103 early CRC (AJCC I/II) and 233 CTRL participants available in the OOF subset, only studies containing both groups contribute to the within-study ROC. The figure states eligible participant counts. Case-control pair counts are not independent sample sizes.',
'Each study ROC is weighted by its number of within-study early-CRC/CTRL pairs; the curves are then averaged over all five repeats at a common false-positive rate. Vertical segments and linear tied-score segments are retained; no smoothing is applied. Curve areas equal the five-repeat mean audited within-study pairwise AUC. Labels show mean AUC; between-repeat SD is in the summary table and is not a confidence interval.',
'The fixed marker panel was outcome-informed using the full Ob cohort. Estimates remain exploratory and do not independently validate marker discovery. No confidence intervals or significance tests are displayed.')
writeLines(legend,file.path(ROOT,'tables',paste0(STEM,'_legend.txt')))
section <- c('', '## Early CRC from original CRC-CTRL OOF scores (2026-10-21)', 'Status: exploratory',
'Purpose: assess early CRC without early-stage-specific fitting or tuning, using every existing eligible OOF prediction; LASSO omitted as requested.',
'Archetype: single quantitative ROC panel. R-only drawing; inherited colors and typography; original four-panel fixed-test figure is unchanged.',
'ROC: study ROC weights = n_CTRL * n_earlyCRC, averaged across the original five CV repeats. No smoothing or best-repeat selection. Integrated AUC agrees with the saved stage OOF audit within 1e-12.',
'Filters: joint panel; four requested models; CTRL plus AJCC I/II CRC. Studies without both groups are excluded from curve only; complete before/after counts and participants are exported.',
'CV OOF evaluation differs from application of a single frozen model. Original outcome-informed feature discovery remains a limitation. SD across repeats is descriptive, not a CI.',
'Run: Rscript plot_early_crc_saved_oof.R',paste0('Output stem: ',STEM),
'Exports: 105 x 121 mm PDF/SVG and 300 dpi PNG preview; Arial. Source curves, labels, per-study/repeat curves, participants and numerical QC accompany the figure.',
'Helper dependency: three pure ROC geometry functions are loaded from plot_within_study_transfer_four_panels.R via parsed function assignments only; the full script is never sourced or executed.')
readme <- file.path(ROOT,'README.md')
if(!any(grepl(section[2],readLines(readme,warn=FALSE),fixed=TRUE))) write(section,file=readme,append=TRUE)
outs <- c(paste0('figures/',STEM,c('.pdf','.svg','.png')),
 paste0('data_for_figures/',STEM,c('_curves.tsv','_labels.tsv','_samples.tsv','_study_repeat_curves.tsv')),
 paste0('tables/',STEM,c('_repeat_AUC_QC.tsv','_study_weights.tsv','_summary.tsv','_legend.txt')),'plot_early_crc_saved_oof.R')
for(manifest in c('input_files.txt','output_files.txt')) {
 target <- file.path(ROOT,manifest); old <- readLines(target,warn=FALSE)
 values <- if(manifest=='input_files.txt') c(unname(INPUTS),helper_file) else outs
 write(setdiff(values,old),file=target,append=TRUE)
}
print(ann); cat('EARLY_CRC_OOF_PLOT_COMPLETE\n')
