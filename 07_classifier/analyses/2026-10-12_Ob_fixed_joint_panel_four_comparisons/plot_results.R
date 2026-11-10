#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(svglite); library(ragg)})
root <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-10-12_Ob_fixed_joint_panel_four_comparisons'
tag <- 'Ob22_fourComp_nested5x5_'
width_mm = 145
height_mm = 122
roc_data <- fread(file.path(root,'data_for_figures',paste0(tag,'Fig_representative_ROC_curves.tsv')))
summary_data <- fread(file.path(root,'tables',paste0(tag,'summary.tsv')))
comparisons <- unique(roc_data$Comparison)
model_colors <- c('Random forest'='#B64E42','LASSO'='#547A9E','Elastic Net'='#729B8B','Radial SVM'='#8067A1','XGBoost'='#D09A45')
for (comp in comparisons) {
  dd <- roc_data[Comparison==comp]
  labels <- summary_data[Comparison==comp]
  labels <- labels[match(names(model_colors),Model)]
  labels <- labels[!is.na(Model)]
  labels[,label_y:=seq(0.10,0.10+0.07*(.N-1),length.out=.N)]
  p <- ggplot(dd,aes(FPR,TPR,color=Model,group=Model)) +
    geom_step(linewidth=0.8,direction='vh') +
    geom_abline(slope=1,intercept=0,linetype='dashed',color='grey75') +
    geom_text(data=labels,aes(x=0.98,y=label_y,label=sprintf('%s  %.3f (SD %.3f)',Model,Mean_AUC,SD_between_repeats)),inherit.aes=FALSE,hjust=1,size=3,color='#333333') +
    scale_color_manual(values=model_colors,drop=FALSE) +
    coord_equal(xlim=c(0,1),ylim=c(0,1),expand=FALSE) +
    labs(title=comp,x='False-positive rate',y='True-positive rate',color=NULL) +
    theme_classic(base_size=12,base_family='Arial') +
    theme(legend.position='bottom',plot.title=element_text(face='bold',hjust=0.5),axis.line=element_line(linewidth=0.4),axis.ticks=element_line(linewidth=0.4))
  stem <- gsub('[^A-Za-z0-9]+','_',comp)
  svglite::svglite(file.path(root,'figures',paste0(tag,'ROC_',stem,'.svg')),width=width_mm/25.4,height=height_mm/25.4)
  print(p)
  dev.off()
  grDevices::cairo_pdf(file.path(root,'figures',paste0(tag,'ROC_',stem,'.pdf')),width=width_mm/25.4,height=height_mm/25.4,family='Arial')
  print(p)
  dev.off()
  ragg::agg_tiff(file.path(root,'figures',paste0(tag,'ROC_',stem,'.tiff')),width=width_mm/25.4,height=height_mm/25.4,units='in',res=600,compression='lzw')
  print(p)
  dev.off()
  ragg::agg_png(file.path(root,'figures',paste0(tag,'ROC_',stem,'.png')),width=width_mm/25.4,height=height_mm/25.4,units='in',res=300)
  print(p)
  dev.off()
}
cat('FIGURES_REGENERATED_FROM_SAVED_ROC_DATA\n')
