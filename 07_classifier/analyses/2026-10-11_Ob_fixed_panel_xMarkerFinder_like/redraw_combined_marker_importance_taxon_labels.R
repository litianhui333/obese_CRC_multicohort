.libPaths(c('D:/R/library', .libPaths()))
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT_OB <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-10-11_Ob_fixed_panel_xMarkerFinder_like'
ROOT_NW <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-10-13_Nw_Ow_locked_L1_Ob_method'
MAP_FILE <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-09-29_BMI_enriched_SV_function_modelQC/2026-09-30_BMI_specificity_layers/tables/L1_L2_candidate_master.tsv'
OB_TAXON_FILE <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-07-13_bmi_union_overlap_heatmaps/tables/Feature_List_Ob_Refined.csv'
OB_DSV_SPECIFICITY_ROWS <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-09-29_BMI_enriched_SV_function_modelQC/2026-09-30_BMI_specificity_layers/data_for_figures/Heatmap_L1_Ob_dSGV_cohort_retention_allSV_labels_final_rows.tsv'
OB_VSV_SPECIFICITY_ROWS <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-09-29_BMI_enriched_SV_function_modelQC/2026-09-30_BMI_specificity_layers/data_for_figures/Heatmap_L1_Ob_vSGV_multicolor_clustered_taxon_SVspecificity_rows.tsv'
INPUTS <- c(
  Nw=file.path(ROOT_NW,'data_for_figures/seedScreen_Nw_20261012_Fig_locked_marker_importance.tsv'),
  Ob=file.path(ROOT_OB,'data_for_figures/restricted_ObAssoc_L1L2_Fig_locked_marker_importance.tsv'))
stopifnot(all(file.exists(INPUTS)),file.exists(MAP_FILE),
  file.exists(OB_DSV_SPECIFICITY_ROWS),file.exists(OB_VSV_SPECIFICITY_ROWS))
master <- fread(MAP_FILE)
stopifnot(all(c('BMI_enriched_group','SGB','Species','SGB_label') %in% names(master)))
ob_taxa <- fread(OB_TAXON_FILE)
stopifnot(all(c('feature','overlap') %in% names(ob_taxa)))
ob_taxa[,SGB_key:=sub('.*\\|t__','',feature)]
ob_taxa_map <- unique(ob_taxa[,.(SGB_key,Ob_taxon_status=overlap)])
if(anyDuplicated(ob_taxa_map$SGB_key)) stop('Ob taxon-specificity map has conflicting SGB labels')
ob_sv_specificity <- unique(rbindlist(list(
  fread(OB_DSV_SPECIFICITY_ROWS)[,.(Original_feature=SV_feature,SV_specificity)],
  fread(OB_VSV_SPECIFICITY_ROWS)[,.(Original_feature=SV_feature,SV_specificity)]
),use.names=TRUE))
if(anyDuplicated(ob_sv_specificity$Original_feature)) stop('Ob SV-specificity map has conflicting feature labels')

make_one <- function(group, input, root, stem, title) {
  d <- fread(input)[Panel=='joint']
  cleanup_cols <- intersect(names(d),c('SGB_key','Taxon_name','Taxon_short','Region',
    'Ob_taxon_status','SV_specificity','Label_red','Display_html'))
  if(length(cleanup_cols)) d[,(cleanup_cols):=NULL]
  stopifnot(nrow(d)>0L,all(d$Type %in% c('species','SV')))
  d[,SGB_key:=sub('^t__','',SGB)]
  d[is.na(SGB_key),SGB_key:=sub('.*\\|t__','',Original_feature)]
  tax <- unique(master[BMI_enriched_group==group,
    .(SGB_key=sub('^t__','',SGB),Taxon_name=Species,Taxon_short=SGB_label)])
  tax <- unique(tax)
  if(anyDuplicated(tax$SGB_key)) stop('Taxon mapping is not unique for ',group)
  d <- merge(d,tax,by='SGB_key',all.x=TRUE,sort=FALSE)
  if(group=='Ob') d <- merge(d,ob_taxa_map,by='SGB_key',all.x=TRUE,sort=FALSE)
  else d[,Ob_taxon_status:=NA_character_]
  if(group=='Ob') {
    d <- merge(d,ob_sv_specificity,by='Original_feature',all.x=TRUE,sort=FALSE)
    if(any(is.na(d[Type=='SV',SV_specificity]))) stop('An Ob joint-panel SV is missing its SV-specificity mapping')
  } else d[,SV_specificity:=NA_character_]
  d[,Label_red:=group=='Ob' & ((Type=='species' & !is.na(Ob_taxon_status) & Ob_taxon_status=='Ob-specific') |
    (Type=='SV' & !is.na(SV_specificity) & SV_specificity=='Ob-specific SV'))]
  # Keep one or two intervals explicit. For longer composite regions, show the
  # genomic span and segment count; Original_feature retains every exact interval.
  compact_region <- function(feature) {
    region <- gsub('_','-',sub('.*:','',feature))
    pieces <- trimws(strsplit(region,';',fixed=TRUE)[[1]])
    if(length(pieces)>2L) {
      coords <- suppressWarnings(as.numeric(unlist(regmatches(pieces,
        gregexpr('[0-9]+',pieces,perl=TRUE)))))
      coords <- coords[is.finite(coords)]
      if(length(coords)>=2L) return(paste0(min(coords),'-',max(coords),
        ' (',length(pieces),' segments)'))
    }
    paste(pieces,collapse='; ')
  }
  d[,Region:=vapply(Original_feature,compact_region,character(1))]
  d[,Display:=fifelse(Type=='SV',
    paste0(fifelse(!is.na(Taxon_short)&nzchar(Taxon_short),Taxon_short,Taxon_name),
      ' (',SGB_key,'): ',Region,' (',SV_type,')'),
    fifelse(!is.na(Display)&nzchar(Display),Display,
      gsub('_',' ',sub('\\|t__',' ',sub('^s__','',Original_feature)))))]
  if(any(d[Type=='SV',is.na(Taxon_name)])) stop('Missing taxon name for an SV in ',group)
  # Put a single line break at the taxon/region boundary. This keeps labels
  # readable without allowing composite regions to collide with adjacent rows.
  wrap_label <- function(z,type) {
    vapply(seq_along(z),function(i) {
      s <- z[[i]]
      if(type[[i]]!='SV') s <- sub(': ', ':\n', s, fixed=TRUE)
      s
    },character(1))
  }
  d[,Display_plot:=wrap_label(Display,Type)]
  d[,Display_plot:=factor(Display_plot,levels=unique(Display_plot[order(Importance)]))]
  d[,Display_html:=gsub('\n','<br>',as.character(Display_plot),fixed=TRUE)]
  d[,Display_html:=ifelse(Label_red,
    paste0("<span style='color:#B33A3A'>",Display_html,'</span>'),Display_html)]
  setorder(d,Type,Importance)
  out_data <- file.path(root,'data_for_figures',paste0(stem,'.tsv'))
  out_base <- file.path(root,'figures',stem)
  fwrite(d,out_data,sep='\t',na='NA')
  g <- ggplot(d,aes(Importance,Display_plot,colour=Type))+
    geom_segment(aes(x=0,xend=Importance,yend=Display_plot),linewidth=.55,colour='#C9C9C9')+
    geom_point(size=2.2)+
    scale_colour_manual(values=c(species='#356E70',SV='#C0523D'))+
    scale_y_discrete(labels=function(x) d$Display_html[match(x,as.character(d$Display_plot))])+
    labs(x='Permutation importance',y=NULL,title=title)+
    theme_classic(base_size=12,base_family='Arial')+
    theme(axis.line=element_line(linewidth=.45,colour='#333333'),
      axis.text.y=ggtext::element_markdown(size=if(group=='Ob') 9 else 8.5,
        hjust=1,halign=1,lineheight=.95,
        margin=margin(r=if(group=='Ob') 7 else 4)),
      legend.title=element_blank(),legend.position='top',
      plot.title=element_text(face='bold'),
      plot.margin=if(group=='Ob') margin(5,6,5,20) else margin(6,12,6,6))
  # Narrow the canvas while ggplot keeps the y-axis label grobs at the same
  # physical size; the point-line panel therefore gives up most of the width.
  w <- if(group=='Nw') 5.2 else 5.9
  wrapped <- as.character(d$Display_plot)
  if(any(grepl('b$',wrapped))) stop('Unexpected truncated label detected')
  extra_lines <- sum(pmax(0L,lengths(strsplit(wrapped,'\n',fixed=TRUE))-1L))
  h <- if(group=='Ob') max(4.8,1.30+.17*nrow(d)+.08*extra_lines) else
    max(4.8,1.8+.22*nrow(d)+.17*extra_lines)
  # The axis-label grobs are sized in points, so narrowing the device preserves
  # their readable width while compressing the point-line plotting area.
  plot_width <- w
  grDevices::cairo_pdf(paste0(out_base,'.pdf'),width=plot_width,height=h,family='Arial')
  print(g);grDevices::dev.off()
  ggsave(paste0(out_base,'.png'),g,width=plot_width,height=h,dpi=300)
  if(requireNamespace('svglite',quietly=TRUE)){
    svglite::svglite(paste0(out_base,'.svg'),width=plot_width,height=h);print(g);grDevices::dev.off()
  }
  qc <- data.table(Group=group,N_markers=nrow(d),N_SV=sum(d$Type=='SV'),
    N_red_specific_taxon_labels=sum(d$Label_red & d$Type=='species'),
    N_red_ob_specific_SV_labels=sum(d$Label_red & d$Type=='SV'),
    Point_line_panel_width_reduction=1/3,Export_width_in=plot_width,Export_height_in=h,
    Ob_SV_labels_single_line=group!='Ob'||all(!grepl('\\n',d[Type=='SV',Display_plot],fixed=TRUE)),
    All_SV_taxa_mapped=all(!is.na(d[Type=='SV',Species_taxon])),
    Only_joint_panel=uniqueN(d$Panel)==1L&&unique(d$Panel)=='joint',
    Plot_data_exists=file.exists(out_data))
  if(!all(qc[,.(All_SV_taxa_mapped,Only_joint_panel,Plot_data_exists)])) stop('Plot QC failed: ',group)
  fwrite(qc,file.path(root,'tables',paste0(stem,'_QC.tsv')),sep='\t')
  message(group,': ',nrow(d),' combined markers, ',sum(d$Type=='SV'),' SV labels mapped; wrote ',out_base)
}
targets <- commandArgs(trailingOnly=TRUE)
if(!length(targets)) targets <- c('Nw','Ob')
if(any(!targets %in% c('Nw','Ob'))) stop('Use no argument, Nw, or Ob')
if('Nw' %in% targets) make_one('Nw',INPUTS[['Nw']],ROOT_NW,
  'seedScreen_Nw_20261012_Fig_locked_marker_importance',
  'Nw combined markers')
if('Ob' %in% targets) make_one('Ob',INPUTS[['Ob']],ROOT_OB,
  'restricted_ObAssoc_L1L2_Fig_locked_marker_importance',
  'Ob combined markers')
