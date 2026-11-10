#!/usr/bin/env Rscript
# Targeted plotting only: saved fixed-model predictions; no fitting or tuning.
.libPaths(c('D:/R/library', .libPaths()))
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
ROOT <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-10-12_Ob_fixed_joint_panel_four_comparisons'
STEM <- 'Fig_Ob_fixedModel_withinStudy_fourComparisons_AAreversed_noLASSO'
REVERSE_ADENOMA <- TRUE
MODEL_ORDER <- c('Random forest','XGBoost','Radial SVM','Elastic Net')
COLORS <- c('Random forest'='#C34E40','XGBoost'='#D99A30','Radial SVM'='#8062A4','Elastic Net'='#72A090')
SOURCE <- file.path(ROOT,'tables/Ob22_CTRLCRC_model_transfer_predictions.tsv')
AUDIT <- file.path(ROOT,'tables/Ob22_CTRLCRC_model_transfer_audit_pairwise_AUC.tsv')
existing_figures <- list.files(file.path(ROOT,'figures'),full.names=TRUE)
BEFORE <- file.info(existing_figures[!startsWith(basename(existing_figures),STEM)])
raw <- fread(SOURCE, colClasses=list(character='sample_id'))[Model %in% MODEL_ORDER]
audit <- fread(AUDIT)
stopifnot(!anyDuplicated(raw[,.(Model,sample_id)]),all(is.finite(raw$score)),!anyNA(raw$Study))
raw[,Group:=fcase(Disease=='CTRL','CTRL',Disease=='AA','Adenoma',Disease=='CRC','CRC')]
raw[,Stage:=fcase(Disease=='CRC' & AJCC_stage %in% c('I','II'),'Early CRC',
                  Disease=='CRC' & AJCC_stage %in% c('III','IV'),'Late CRC',default=Group)]
comparisons <- data.table(Panel=LETTERS[1:4],
  Comparison=c('Early CRC vs CTRL','Late CRC vs CTRL','Adenoma vs CTRL','CRC vs Adenoma'),
  Audit_comparison=c('Early CRC vs HC','Late CRC vs HC','Adenoma vs HC','Adenoma vs CRC'),
  Negative=c('CTRL','CTRL','CTRL','Adenoma'), Positive=c('Early CRC','Late CRC','Adenoma','CRC'),
  Title=c('Obese group\nEarly CRC (Stage I-II) vs CTRL','Obese group\nLate CRC (Stage III-IV) vs CTRL',
          'Obese group\nAdenoma vs CTRL','Obese group\nCRC vs Adenoma'))
# Local ROC uses grouped score ties: segments are linear for tied case/control scores.
local_roc <- function(z) {
  s <- z[,.(n0=sum(Role=='negative'),n1=sum(Role=='positive')),by=Plot_score]
  setorder(s,-Plot_score)
  rbind(data.table(FPR=0,TPR=0),s[,.(FPR=cumsum(n0)/sum(n0),TPR=cumsum(n1)/sum(n1))])
}
trap_auc <- function(z) sum(diff(z$FPR)*(head(z$TPR,-1)+tail(z$TPR,-1))/2)
# Preserve both ends of vertical segments; interpolate only across nonzero FPR spans.
roc_at <- function(curve,grid) {
  v <- curve[,.(lo=min(TPR),hi=max(TPR)),by=FPR][order(FPR)]
  ans <- matrix(NA_real_,nrow=length(grid),ncol=2)
  for(j in seq_along(grid)) {
    x <- grid[j]; exact <- match(x,v$FPR)
    if(!is.na(exact)) ans[j,] <- c(v$lo[exact],v$hi[exact]) else {
      k <- findInterval(x,v$FPR)
      y <- v$hi[k]+(x-v$FPR[k])/(v$FPR[k+1]-v$FPR[k])*(v$lo[k+1]-v$hi[k])
      ans[j,] <- c(y,y)
    }
  }
  ans
}
curves <- list(); summaries <- list(); study_tables <- list(); plot_samples <- list(); study_rocs <- list()
for(i in seq_len(nrow(comparisons))) for(method in MODEL_ORDER) {
  cp <- comparisons[i]; d <- copy(raw[Model==method])
  d[,Role:=fcase((if(cp$Negative=='CRC') Group else Stage)==cp$Negative,'negative',
                (if(cp$Positive=='CRC') Group else Stage)==cp$Positive,'positive',default='excluded')]
  z <- d[Role!='excluded']
  reverse <- REVERSE_ADENOMA && cp$Positive=='Adenoma'
  z[,Plot_score:=if(reverse) 1-score else score]
  ct <- z[,.(N_negative=sum(Role=='negative'),N_positive=sum(Role=='positive')),by=Study]
  ct[,`:=`(Pair_count=N_negative*N_positive,Eligible=N_negative>0 & N_positive>0)]
  ct[,Weight:=Pair_count/sum(Pair_count)]
  ct[,`:=`(Model=method,Panel=cp$Panel,Comparison=cp$Comparison,Reversed=reverse)]
  study_tables[[length(study_tables)+1L]] <- ct
  eligible <- ct[Eligible==TRUE]
  stopifnot(nrow(eligible)>0L)
  cs <- lapply(eligible$Study,function(study_id) local_roc(z[Study==study_id]))
  grid <- sort(unique(unlist(lapply(cs,function(v) v$FPR))))
  yy <- matrix(0,nrow=length(grid),ncol=2)
  for(k in seq_along(cs)) {
    yy <- yy+eligible$Weight[k]*roc_at(cs[[k]],grid)
    sr <- copy(cs[[k]])
    sr[,`:=`(Model=method,Panel=cp$Panel,Comparison=cp$Comparison,Study=eligible$Study[k],
             Weight=eligible$Weight[k],Point_order=.I,Reversed=reverse)]
    study_rocs[[length(study_rocs)+1L]] <- sr
  }
  curve <- data.table(FPR=rep(grid,each=2L),TPR=as.vector(t(yy)))
  ar <- audit[Comparison==cp$Audit_comparison & Model==method]
  stopifnot(nrow(ar)==1L,sum(ct$Pair_count)==ar$N_same_study_pairs)
  target <- if(reverse) 1-ar$AUC_same_study_pairs else ar$AUC_same_study_pairs
  auc <- trap_auc(curve)
  stopifnot(abs(auc-target)<1e-12,all(diff(curve$FPR)>=0),all(diff(curve$TPR)> -1e-12))
  used <- z[Study %in% eligible$Study]
  ss <- data.table(Panel=cp$Panel,Comparison=cp$Comparison,Model=method,AUC=auc,
    Audit_AUC=target,Absolute_AUC_error=abs(auc-target),Reversed=reverse,
    N_negative_all=sum(z$Role=='negative'),N_positive_all=sum(z$Role=='positive'),
    N_negative_eligible=sum(used$Role=='negative'),N_positive_eligible=sum(used$Role=='positive'),
    N_eligible_studies=nrow(eligible),N_same_study_pairs=sum(ct$Pair_count),
    Model_order=match(method,MODEL_ORDER),Label=sprintf('%s %.3f',method,auc),
    Label_x=.98,Label_y=.30-(match(method,MODEL_ORDER)-1L)*.065)
  summaries[[length(summaries)+1L]] <- ss
  curve[,`:=`(Panel=cp$Panel,Comparison=cp$Comparison,Model=method,AUC=auc,
               Reversed=reverse,Point_order=.I,Model_order=match(method,MODEL_ORDER))]
  curves[[length(curves)+1L]] <- curve
  z[,`:=`(Panel=cp$Panel,Comparison=cp$Comparison,Included_in_curve=Study %in% eligible$Study,Reversed=reverse)]
  plot_samples[[length(plot_samples)+1L]] <- z
}
roc <- rbindlist(curves); summary <- rbindlist(summaries)
counts <- rbindlist(study_tables); sample_data <- rbindlist(plot_samples)
roc[,Model:=factor(Model,levels=MODEL_ORDER)]; summary[,Model:=factor(Model,levels=MODEL_ORDER)]
setorder(roc,Panel,Model,Point_order); setorder(summary,Panel,Model)
stopifnot(nrow(summary)==16L,!any(as.character(roc$Model)=='LASSO'))
for(folder in c('figures','tables','data_for_figures')) dir.create(file.path(ROOT,folder),showWarnings=FALSE)
fwrite(roc,file.path(ROOT,'data_for_figures',paste0(STEM,'_curves.tsv')),sep='\t')
fwrite(summary,file.path(ROOT,'data_for_figures',paste0(STEM,'_labels.tsv')),sep='\t')
fwrite(sample_data,file.path(ROOT,'data_for_figures',paste0(STEM,'_samples.tsv')),sep='\t',na='NA')
fwrite(rbindlist(study_rocs),file.path(ROOT,'data_for_figures',paste0(STEM,'_study_curves.tsv')),sep='\t')
fwrite(counts,file.path(ROOT,'tables',paste0(STEM,'_study_weights.tsv')),sep='\t')
fwrite(summary,file.path(ROOT,'tables',paste0(STEM,'_AUC_QC.tsv')),sep='\t')
panels <- lapply(seq_len(nrow(comparisons)),function(i) {
  cp <- comparisons[i]; ss <- summary[Panel==cp$Panel]; nn <- ss[1L]
  note <- sprintf('Within-study pairs: %s; CTRL/AA/CRC participants retained by study',nn$N_same_study_pairs)
  subtitle <- sprintf('n = %d vs %d; %d within-study pairs',nn$N_positive_eligible,nn$N_negative_eligible,nn$N_same_study_pairs)
  if(cp$Panel=='C') subtitle <- paste0('Reversed CRC score (1 - score)\n',subtitle)
  ggplot(roc[Panel==cp$Panel],aes(FPR,TPR,colour=Model,group=Model))+
    geom_abline(slope=1,intercept=0,linetype='dashed',colour='#C5C5C5',linewidth=.4)+
    geom_path(linewidth=.7)+
    geom_text(data=ss,aes(x=Label_x,y=Label_y,label=Label),hjust=1,fontface='bold',size=3.2,show.legend=FALSE)+
    scale_colour_manual(values=COLORS,drop=FALSE)+
    scale_x_continuous(breaks=c(0,.25,.5,.75,1),expand=expansion(mult=.005))+
    scale_y_continuous(breaks=c(0,.25,.5,.75,1),expand=expansion(mult=.005))+
    coord_fixed(xlim=c(0,1),ylim=c(0,1),clip='off')+
    labs(title=cp$Title,subtitle=subtitle,x='Within-study false-positive rate',y='Within-study true-positive rate')+
    theme_classic(base_size=12,base_family='Arial')+
    theme(legend.position='none',plot.title=element_text(face='bold',size=12,hjust=.5),
      plot.subtitle=element_text(size=9,hjust=.5,margin=margin(b=6)),axis.text=element_text(size=10),
      axis.title=element_text(size=11),plot.margin=margin(9,12,8,10))
})
fig <- wrap_plots(panels,ncol=2)+plot_annotation(tag_levels='A',
  caption='Fixed CRC-CTRL model; no refitting or tuning. Curves weighted by within-study case-control pairs.',
  theme=theme(plot.tag=element_text(face='bold',size=16),plot.caption=element_text(size=9,hjust=0)))
base <- file.path(ROOT,'figures',STEM)
grDevices::cairo_pdf(paste0(base,'.pdf'),width=183/25.4,height=192/25.4,family='Arial')
print(fig); grDevices::dev.off()
svglite::svglite(paste0(base,'.svg'),width=183/25.4,height=192/25.4)
print(fig); grDevices::dev.off()
ggsave(paste0(base,'.png'),fig,width=183,height=192,units='mm',device=ragg::agg_png,dpi=300,bg='white')
# Guard against accidental changes to unrelated pre-existing figures.
AFTER <- file.info(rownames(BEFORE))
stopifnot(identical(BEFORE$size,AFTER$size),identical(BEFORE$mtime,AFTER$mtime))
legend <- c('Fixed-model within-study ROC evaluation across colorectal neoplasia comparisons in obese individuals.',
'(A) Early CRC versus CTRL. (B) Late CRC versus CTRL. (C) Adenoma versus CTRL, using reversed CRC scores (1 - score), selected after inspecting the original score direction and therefore exploratory. (D) CRC versus adenoma.',
'Predictions were taken from the original fixed 22-marker CRC-CTRL classifier (10 species, 12 SVs, age and sex); no classifier was fitted or tuned in this plotting workflow. Original CTRL/CRC test subjects and previously untrained adenoma subjects were used.',
'Only studies containing both groups in a given comparison contributed to its curve. Within each study, empirical ROC coordinates used linear segments for tied scores; local curves were averaged at a common false-positive rate with weights proportional to the number of within-study case-control pairs. Integrated curve areas equal the previously audited within-study pairwise AUCs. These are weighted within-study curves, not ROC curves at a single common raw-score threshold across studies.',
'Panel n values indicate eligible participants (positive group followed by negative group); pair counts are not independent sample sizes. Colored labels indicate AUCs; dashed diagonal lines indicate chance discrimination. LASSO is omitted as requested. No confidence intervals or significance tests are displayed.',
'The upstream marker list was outcome-informed in the full Ob cohort. These results are exploratory fixed-panel evaluations, not fully independent marker-discovery validation.')
writeLines(legend,file.path(ROOT,'tables',paste0(STEM,'_legend.txt')))
readme <- file.path(ROOT,'README.md')
section <- c('', '## Targeted within-study fixed-model four-panel ROC (2026-10-21)',
'Status: exploratory',
'Conclusion: the fixed CRC classifier retains within-study CRC discrimination; adenoma versus CTRL remains weak even when its score is reversed.',
'Archetype: quantitative 2-by-2 grid. Style-only inheritance from the earlier ROC layout; the aggregation statistic is newly implemented.',
'Scope: four comparisons in Ob; Random forest, XGBoost, Radial SVM, Elastic Net only. No model training, tuning, marker selection or original-script edits.',
'Input: saved direct-transfer predictions and the audited within-study pairwise AUC table. Studies lacking one comparison group are excluded only for that comparison; all counts, exclusions and weights are exported.',
'ROC definition: weighted mean of study-specific empirical ROC curves with weights n_negative * n_positive. Both limits of vertical segments and linear tied-score segments are retained. Numerical integration must match audited AUC within 1e-12.',
'Adenoma versus CTRL uses 1 - CRC score by explicit user request after the initial results were inspected; this reversal is exploratory. Other comparisons retain original score direction.',
'Participant counts are those in eligible studies. Pair counts are dependent observations. No CIs or hypothesis tests were added. Very few early-stage within-study pairs limit interpretation.',
'Exports: 183 x 192 mm PDF/SVG, PNG at 300 dpi; Arial, readable labels; matching curve, annotation, participant and per-study plotting data.',
'Run: Rscript plot_within_study_transfer_four_panels.R',paste0('Output stem: ',STEM),
'QC: no LASSO; 16 model-comparison estimates; exact AUC agreement; unrelated existing figure sizes/timestamps checked unchanged.')
if (!any(grepl(section[2L],readLines(readme,warn=FALSE),fixed=TRUE))) write(section,file=readme,append=TRUE)
out_files <- c(paste0('figures/',STEM,c('.pdf','.svg','.png')),
  paste0('data_for_figures/',STEM,c('_curves.tsv','_labels.tsv','_samples.tsv','_study_curves.tsv')),
  paste0('tables/',STEM,c('_study_weights.tsv','_AUC_QC.tsv','_legend.txt')),
  'plot_within_study_transfer_four_panels.R')
existing_outputs <- readLines(file.path(ROOT,'output_files.txt'),warn=FALSE)
write(setdiff(out_files,existing_outputs),file=file.path(ROOT,'output_files.txt'),append=TRUE)
print(summary[,.(Panel,Model,AUC,N_positive_eligible,N_negative_eligible,N_same_study_pairs,Absolute_AUC_error)])
cat('TARGETED_WITHIN_STUDY_ROC_COMPLETE\n')


