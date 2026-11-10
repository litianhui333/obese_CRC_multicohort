#!/usr/bin/env Rscript
# Targeted update: original Ob L1 heatmap display, abundance-adjusted specificity.
# No model fits or residualization; all legacy outputs and scripts are read-only.
.libPaths(c('D:/R/library',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(ComplexHeatmap);library(circlize);library(grid)})
SOURCE_ROOT <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-09-29_BMI_enriched_SV_function_modelQC/2026-09-30_BMI_specificity_layers'
OUTPUT_DIR <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-10-13_SV_abundance_sensitivity'
SUFFIX <- '_abundanceAdjusted_20261026'
FORMAL <- file.path(OUTPUT_DIR,'tables/formal_effect_specific_candidate_reporting.tsv')
MODELS <- file.path(OUTPUT_DIR,'tables/run_model_results_long.tsv')
DOTPLOT <- file.path(OUTPUT_DIR,'data_for_figures/Fig_Ob_specific_abundanceAdjusted_crossBMI.tsv')
formal <- fread(FORMAL)[BMI_enriched_group=='Ob']
adjusted <- fread(MODELS)[Model=='abundance_adjusted' & BMI=='Ob']
dot <- fread(DOTPLOT)
stopifnot(nrow(formal)==31L,formal[SV_type=='dSGV',.N]==21L,formal[SV_type=='vSGV',.N]==10L,
 setequal(formal$BMI_SV_ID,unique(dot$BMI_SV_ID)),!anyDuplicated(adjusted$BMI_SV_ID))
SOURCE_SCRIPTS <- c(dSGV=file.path(SOURCE_ROOT,'plot_ob_L1_dSV_allSV_labels_final.R'),
 vSGV=file.path(SOURCE_ROOT,'plot_ob_L1_vSV_clustered_taxon_SVspecificity.R'))
BASE_STEMS <- c(dSGV='Heatmap_L1_Ob_dSGV_cohort_retention_allSV_labels_final',
 vSGV='Heatmap_L1_Ob_vSGV_multicolor_clustered_taxon_SVspecificity')
NEW_STEMS <- paste0(BASE_STEMS,SUFFIX);names(NEW_STEMS)<-names(BASE_STEMS)
for(d in c('figures','tables','data_for_figures'))dir.create(file.path(OUTPUT_DIR,d),showWarnings=FALSE)
# Protect every existing legacy figure, plus all unrelated sensitivity figures and sources.
legacy <- list.files(file.path(SOURCE_ROOT,'figures'),full.names=TRUE)
current <- list.files(file.path(OUTPUT_DIR,'figures'),full.names=TRUE)
current <- current[!vapply(basename(current),function(n)any(startsWith(n,NEW_STEMS)),logical(1))]
protected <- unique(c(legacy,current,SOURCE_SCRIPTS,FORMAL,MODELS,DOTPLOT,file.path(SOURCE_ROOT,'prepared_data.rds')))
before <- tools::md5sum(protected)
assignment_name <- function(e) if(is.call(e)&&identical(e[[1L]],as.name('<-'))&&is.symbol(e[[2L]]))as.character(e[[2L]])else ''
# Redirect only paths explicitly built using target STEM. All source reads stay in SOURCE_ROOT.
redirect_outputs <- function(e) {
 if(!is.call(e))return(e)
 if(identical(e[[1L]],as.name('file.path')) && length(e)>=3L && identical(e[[2L]],as.name('ROOT')) && 'STEM' %in% all.names(e))e[[2L]]<-as.name('OUTPUT_DIR')
 for(j in seq_along(e))if(j>1L && is.call(e[[j]]))e[[j]]<-redirect_outputs(e[[j]])
 e
}
upgrade_rows <- function(rows,type) {
 rows<-copy(rows)
 old_beta<-rows$Target_Beta
 rows[,`:=`(L2_nominal_before_abundance=L2_nominal,Target_Beta_before_abundance=Target_Beta)]
 mm<-adjusted[match(rows$BMI_SV_ID,BMI_SV_ID)]
 stopifnot(!anyNA(mm$BMI_SV_ID),all(is.finite(mm$Disease_Beta)),all(is.finite(mm$Disease_SE)),
  all(is.na(mm$Error)|!nzchar(mm$Error)),identical(sign(old_beta),sign(mm$Disease_Beta)))
 rows[,`:=`(L2_nominal=BMI_SV_ID %in% formal$BMI_SV_ID,L2_selected_BH=FALSE,
  Target_Beta=mm$Disease_Beta,Target_SE=mm$Disease_SE,Target_P=mm$Disease_P,
  Adjusted_model_N=mm$N,Adjusted_model_N_CTRL=mm$N_CTRL,Adjusted_model_N_CRC=mm$N_CRC,
  Specificity_basis='Both abundance-adjusted Ob-versus-Nw/Ow Wald P < 0.05; exploratory')]
 expected<-formal[SV_type==type,BMI_SV_ID]
 stopifnot(setequal(rows[L2_nominal==TRUE,BMI_SV_ID],expected),sum(rows$L2_nominal)==if(type=='dSGV')21L else 10L)
 rows
}
# Export with explicit editable vector devices; inherit the original dimensions.
export_original_dimensions <- function(render,path,width,height) {
 grDevices::cairo_pdf(paste0(path,'.pdf'),width=width,height=height,family='Arial')
 render();grDevices::dev.off()
 grDevices::png(paste0(path,'.png'),width=width,height=height,units='in',res=300,type='cairo')
 render();grDevices::dev.off()
 svglite::svglite(paste0(path,'.svg'),width=width,height=height)
 render();grDevices::dev.off()
}
outputs<-character();all_rows<-list();qa<-list()
for(type in names(SOURCE_SCRIPTS)) {
 env<-new.env(parent=globalenv());env$OUTPUT_DIR<-OUTPUT_DIR
 ast<-as.list(parse(SOURCE_SCRIPTS[type]));updated<-FALSE
 # Reuse original Arial typography, palettes, dimensions, original row/sample order,
 # cohort retention aggregation, missingness, and raw-z vSV display.
 for(e in ast) {
  n<-assignment_name(e)
  head_name <- if(is.call(e))paste(deparse(e[[1L]]),collapse='')else ''
  if(head_name %in% c('grDevices::cairo_pdf','png','svglite::svglite','render','dev.off'))next
  if(n=='append_unique')break # stop before original manifest-writing footer
  if(n=='ROOT'){env$ROOT<-SOURCE_ROOT;next}
  if(n=='STEM'){env$STEM<-unname(NEW_STEMS[type]);next}
  if(type=='vSGV' && n=='STYLE'){env$STYLE<-'raw_z';next}
  # Old specificity counts are the only obsolete assertions in original scripts.
  txt<-paste(deparse(e,width.cutoff=500L),collapse='\n')
  if(grepl('sum(rows$L2_nominal)',txt,fixed=TRUE)) {
   txt<-gsub(if(type=='dSGV')'== 24L'else'== 11L',if(type=='dSGV')'== 21L'else'== 10L',txt,fixed=TRUE)
   e<-parse(text=txt)[[1L]]
  }
  e<-redirect_outputs(e)
  eval(e,envir=env)
  if(n %in% c('figure_path','path'))export_original_dimensions(env$render,get(n,envir=env),env$width,env$height)
  if(!updated && n=='rows' && exists('rows',envir=env,inherits=FALSE)) {
   env$rows<-upgrade_rows(env$rows,type);updated<-TRUE
  }
 }
 stopifnot(updated)
 rows<-copy(env$rows)
 stopifnot(nrow(rows)==if(type=='dSGV')68L else 32L,
  setequal(rows[L2_nominal==TRUE,BMI_SV_ID],formal[SV_type==type,BMI_SV_ID]))
 # Verify source plotting quantities themselves were not residualized or altered.
 if(type=='dSGV') {
  old_long<-fread(file.path(SOURCE_ROOT,'data_for_figures',paste0(BASE_STEMS[type],'_long.tsv')))
  now_long<-copy(env$long);setorder(old_long,BMI_SV_ID,Stage,Cohort);setorder(now_long,BMI_SV_ID,Stage,Cohort)
  stopifnot(isTRUE(all.equal(old_long$Retention_rate,now_long$Retention_rate,tolerance=1e-12)),
   isTRUE(all.equal(old_long$Observed_N,now_long$Observed_N,tolerance=0,check.attributes=FALSE)))
 } else {
  old_mat<-fread(file.path(SOURCE_ROOT,'data_for_figures',paste0(BASE_STEMS[type],'.tsv')))
  now_mat<-data.table(BMI_SV_ID=rownames(env$plotmat),as.data.table(env$plotmat))
  stopifnot(identical(names(old_mat),names(now_mat)),identical(old_mat$BMI_SV_ID,now_mat$BMI_SV_ID),
   isTRUE(all.equal(as.matrix(old_mat[,-1]),as.matrix(now_mat[,-1]),tolerance=1e-12,check.attributes=FALSE)))
 }
 rows[,SV_panel:=if(type=='dSGV')'dSV'else'vSV']
 all_rows[[type]]<-rows[,.(BMI_SV_ID,SGB,SV_panel,L2_nominal_before_abundance,L2_nominal,Target_Beta_before_abundance,Target_Beta,Target_SE,Target_P,Specificity_basis)]
 qa[[type]]<-data.table(SV_type=type,N_L1=nrow(rows),N_old_specific=sum(rows$L2_nominal_before_abundance),
  N_adjusted_specific=sum(rows$L2_nominal),N_lost=sum(rows$L2_nominal_before_abundance & !rows$L2_nominal),
  N_gained=sum(!rows$L2_nominal_before_abundance & rows$L2_nominal),Matches_dotplot=TRUE,Underlying_heatmap_values_unchanged=TRUE)
 generated<-unlist(lapply(c('figures','data_for_figures','tables'),function(folder)
  list.files(file.path(OUTPUT_DIR,folder),pattern=unname(NEW_STEMS[type]),full.names=TRUE)))
 outputs<-c(outputs,generated)
}
annotation<-rbindlist(all_rows)
stopifnot(setequal(annotation[L2_nominal==TRUE,BMI_SV_ID],formal$BMI_SV_ID),sum(annotation$L2_nominal)==31L)
matched<-merge(annotation[L2_nominal==TRUE],dot[BMI=='Ob',.(BMI_SV_ID,Dotplot_Beta=Disease_Beta,Dotplot_SE=Disease_SE,Dotplot_P=Disease_P)],by='BMI_SV_ID')
stopifnot(nrow(matched)==31L,max(abs(matched$Target_Beta-matched$Dotplot_Beta))<1e-12,
 max(abs(matched$Target_SE-matched$Dotplot_SE))<1e-12,max(abs(matched$Target_P-matched$Dotplot_P))<1e-12)
qa_path<-file.path(OUTPUT_DIR,'tables/Ob_heatmap_abundanceAdjusted_20261026_QA.tsv')
changes_path<-file.path(OUTPUT_DIR,'tables/Ob_heatmap_abundanceAdjusted_20261026_specificity_crosswalk.tsv')
fwrite(rbindlist(qa),qa_path,sep='\t');fwrite(annotation,changes_path,sep='\t')
stopifnot(identical(unname(before),unname(tools::md5sum(protected))))
legend_path<-file.path(OUTPUT_DIR,'tables/Ob_heatmap_abundanceAdjusted_20261026_legend.txt')
writeLines(c('Ob L1 SV landscape with abundance-adjusted putative obesity-specific annotations.',
'The original Ob L1 display universe is retained: 68 dSVs and 32 vSVs. SV specificity annotation is updated using the same 31 candidates as the abundance-adjusted cross-BMI dotplot (21 dSVs, 10 vSVs). Six old candidates lose the annotation and two gain it.',
'Specificity requires nominal two-sided Ob-versus-Nw and Ob-versus-Ow effect-difference P < 0.05 in abundance-adjusted mixed models; none of the 31 survives the reported BH correction. Age, sex, study random intercept and parent-SGB abundance are included as in the dotplot. This annotation does not indicate confirmed BMI specificity.',
'Heatmap colors remain descriptive measured SV values: dSV panels show cohort-level retention fractions among observed samples by CTRL/Adenoma/CRC; vSV panels retain the original row-z display, clipping, sample filtering and ordering. Colors are NOT abundance-residualized SV values. Adjusted CRC effect estimates determine direction labels; all 100 signs agree with the prior version.',
'Taxon abundance status is a separate, unchanged annotation of SGB-level abundance associations; it must not be confused with adjustment for abundance in SV models. Missing values remain missing.',
'Original display dimensions, palette, labels and stage/sample ordering are inherited. New file suffix abundanceAdjusted_20261026; all old figures, original scripts, source model tables and dotplot hashes were checked unchanged.'),legend_path)
section<-c('', '## Ob heatmap specificity synchronized with abundance-adjusted dotplot (2026-10-26)', 'Status: exploratory',
'User scope: regenerate the two Ob overview heatmaps under different filenames while preserving edited legacy figures. No model rerun.',
'Figure contract: 68-dSV cohort retention heatmap plus 32-vSV sample-level row-z heatmap; the specific annotation is the exact dotplot candidate set (21+10=31). Comparative quantitative grids; same R-only rendering and prior style.',
'Adaptation: original source scripts are parsed and rendered in isolated environments. Source ROOT remains the original folder. Only output paths explicitly containing target STEM are redirected to this analysis folder, with suffix abundanceAdjusted_20261026; original manifest-writing footer is never executed.',
'Adjusted model estimates replace the target effect/SE/P used for annotation and direction. Row order, original display universe, taxon status annotation, colors, dimensions, sample ordering and SV matrix values remain unchanged. All 100 directions are unchanged.',
'Candidate transitions versus old 35: 29 retained, 6 removed (5 dSV,1 vSV), 2 added (both dSV). Updated annotation counts: 21 dSV,10 vSV. This is abundance-adjusted selection/annotation, not residualization of heatmap measurements.',
'QC: exact candidate-set equality with dotplot and formal table; adjusted coefficient/SE/P equality; underlying plotted values equal previous exports; original figures/scripts/model tables/dotplot verified unchanged by MD5.',
'Raw display labels and typography inherited by request, including compact row fonts. SVG/PDF retain editable text; PNG at 300 dpi is a preview. Source plotting frames, adjusted annotations, transition crosswalk and QA are exported.',
'Run: Rscript plot_Ob_heatmaps_abundanceAdjusted_20261026.R')
readme<-file.path(OUTPUT_DIR,'README.md');if(!any(grepl(section[2L],readLines(readme,warn=FALSE),fixed=TRUE)))write(section,file=readme,append=TRUE)
inputs<-unique(c(FORMAL,MODELS,DOTPLOT,SOURCE_SCRIPTS,file.path(SOURCE_ROOT,'prepared_data.rds'),
 file.path(SOURCE_ROOT,'tables/Ob_SV_taxon_abundance_crosswalk.tsv'),
 file.path(SOURCE_ROOT,'data_for_figures/Heatmap_L1_dSGV_allBMI_rows.tsv'),
 file.path(SOURCE_ROOT,'data_for_figures',paste0(BASE_STEMS['dSGV'],'_long.tsv')),
 file.path(SOURCE_ROOT,'data_for_figures',paste0(BASE_STEMS['vSGV'],'.tsv')),
 file.path(SOURCE_ROOT,'data_for_figures',paste0('Heatmap_L1_Ob_vSGV_multicolor_clustered_compact',c('.tsv','_rows.tsv','_columns.tsv','_palette.tsv')))))
outputs<-c(outputs,qa_path,changes_path,legend_path,file.path(OUTPUT_DIR,'plot_Ob_heatmaps_abundanceAdjusted_20261026.R'))
for(manifest in c('input_files.txt','output_files.txt')) {
 f<-file.path(OUTPUT_DIR,manifest);old<-readLines(f,warn=FALSE)
 write(setdiff(if(manifest=='input_files.txt')inputs else outputs,old),file=f,append=TRUE)
}
print(rbindlist(qa));cat('ABUNDANCE_ADJUSTED_HEATMAPS_COMPLETE\n')


