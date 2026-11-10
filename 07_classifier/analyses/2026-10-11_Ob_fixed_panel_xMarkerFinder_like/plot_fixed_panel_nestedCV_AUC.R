#!/usr/bin/env Rscript
.libPaths(c('D:/R/library', .libPaths()))
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

ROOT <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-10-11_Ob_fixed_panel_xMarkerFinder_like'
TAG <- 'restricted_ObAssoc_L1L2_nestedCV_'
summary_path <- file.path(ROOT,'tables',paste0(TAG,'summary.tsv'))
repeat_path <- file.path(ROOT,'tables',paste0(TAG,'repeat_AUC.tsv'))
stopifnot(file.exists(summary_path),file.exists(repeat_path))
fig_dir <- file.path(ROOT,'figures')
data_dir <- file.path(ROOT,'data_for_figures')
dir.create(fig_dir,showWarnings=FALSE,recursive=TRUE)
dir.create(data_dir,showWarnings=FALSE,recursive=TRUE)

panel_levels <- c('species','sv','joint')
panel_labels <- c(species='Species (12 markers)',sv='SV (20 markers)',joint='Species + SV (22 markers)')
model_levels <- c('Random forest','LASSO','Elastic Net','Radial SVM','XGBoost')
model_cols <- c('Random forest'='#777777','LASSO'='#C0523D','Elastic Net'='#356E70',
                'Radial SVM'='#D99A32','XGBoost'='#4778A8')

repdat <- fread(repeat_path)
summary <- fread(summary_path)
stopifnot(setequal(unique(repdat$Panel),panel_levels),
          setequal(unique(repdat$Model),model_levels),
          all(repdat[, .N, by=.(Panel,Model)]$N==5L),
          all(is.finite(repdat$AUC)), all(summary$N_repeats==5L))
repdat[,Panel:=factor(Panel,levels=panel_levels,labels=panel_labels[panel_levels])]
repdat[,Model:=factor(Model,levels=model_levels)]
summary[,Panel:=factor(Panel,levels=panel_levels,labels=panel_labels[panel_levels])]
summary[,Model:=factor(Model,levels=model_levels)]
summary[,`:=`(Lower=Mean_AUC-SD_between_repeats,Upper=Mean_AUC+SD_between_repeats)]

# Save exact inputs to each plotted geom, including facet/model ordering and error limits.
fwrite(repdat,file.path(data_dir,paste0(TAG,'Fig_AUC_repeat_points.tsv')),sep='\t',na='NA')
fwrite(summary,file.path(data_dir,paste0(TAG,'Fig_AUC_mean_SD.tsv')),sep='\t',na='NA')

p <- ggplot() +
  geom_point(data=repdat,aes(x=Model,y=AUC,colour=Model),
             position=position_jitter(width=.075,height=0,seed=136),
             size=1.9,alpha=.72,show.legend=TRUE) +
  geom_errorbar(data=summary,aes(x=Model,ymin=Lower,ymax=Upper),
                width=.18,linewidth=.65,colour='#333333',show.legend=FALSE) +
  geom_point(data=summary,aes(x=Model,y=Mean_AUC),
             shape=23,size=3.2,fill='white',colour='#222222',stroke=.75,show.legend=FALSE) +
  facet_grid(.~Panel,scales='fixed',space='fixed') +
  scale_colour_manual(values=model_cols,breaks=model_levels,drop=FALSE) +
  scale_y_continuous(limits=c(.78,.91),breaks=seq(.78,.90,.03),
                     labels=function(x)sprintf('%.2f',x),expand=expansion(mult=0)) +
  scale_x_discrete(drop=FALSE) +
  labs(title='Obese CRC classification with fixed marker panels',
       subtitle='Five repeats of nested 5-fold cross-validation',
       x=NULL,y='Out-of-fold AUC',colour=NULL,
       caption='Dots: repeat-level AUC (pooled out-of-fold predictions); diamonds: mean; bars: mean plus/minus 1 SD across five repeats.') +
  theme_classic(base_size=8.5,base_family='Arial') +
  theme(axis.line=element_line(linewidth=.4,colour='#333333'),
        axis.ticks=element_line(linewidth=.35,colour='#333333'),
        axis.text.x=element_text(angle=35,hjust=1,vjust=1,size=7.4,colour='#333333'),
        axis.text.y=element_text(size=7.5,colour='#333333'),
        axis.title.y=element_text(size=8.3,margin=margin(r=5)),
        strip.background=element_rect(fill='#F1F1F1',colour=NA),
        strip.text=element_text(face='bold',size=8.2,colour='#222222',margin=margin(t=4,b=4)),
        plot.title=element_text(face='bold',size=10,margin=margin(b=2)),
        plot.subtitle=element_text(size=8,colour='#555555',margin=margin(b=8)),
        plot.caption=element_text(size=7.1,colour='#444444',hjust=0,margin=margin(t=7)),
        legend.position='bottom',legend.text=element_text(size=7.5),
        legend.key.width=grid::unit(.42,'cm'),legend.spacing.x=grid::unit(.16,'cm'),
        panel.spacing.x=grid::unit(.42,'cm'),plot.margin=margin(5,5,5,5))

stem <- file.path(fig_dir,paste0(TAG,'Fig_AUC_repeated_nestedCV'))
ggsave(paste0(stem,'.pdf'),p,width=183/25.4,height=79/25.4,device=cairo_pdf,family='Arial')
svglite::svglite(paste0(stem,'.svg'),width=183/25.4,height=79/25.4)
print(p);dev.off()
ragg::agg_png(paste0(stem,'.png'),width=183/25.4,height=79/25.4,units='in',res=400)
print(p);dev.off()

# Lightweight text/plot-data checks accompany the figure.
qc <- data.table(Check=c('All five repeats present for every panel/model',
                         'Fixed panel sizes reflected in facet labels',
                         'Mean and SD bars use stored repeat summaries',
                         'Plot source data files exist',
                         'PDF, SVG and PNG exports exist'),
                 Passed=c(all(repdat[, .N,by=.(Panel,Model)]$N==5L),
                          identical(levels(repdat$Panel),unname(panel_labels[panel_levels])),
                          all(abs(summary$Lower-(summary$Mean_AUC-summary$SD_between_repeats))<1e-12),
                          file.exists(file.path(data_dir,paste0(TAG,'Fig_AUC_repeat_points.tsv')))&&file.exists(file.path(data_dir,paste0(TAG,'Fig_AUC_mean_SD.tsv'))),
                          all(file.exists(paste0(stem,c('.pdf','.svg','.png'))))))
fwrite(qc,file.path(ROOT,'tables',paste0(TAG,'figure_QC.tsv')),sep='\t')
stopifnot(all(qc$Passed))
print(qc)
cat('FIGURE_COMPLETE\n')


