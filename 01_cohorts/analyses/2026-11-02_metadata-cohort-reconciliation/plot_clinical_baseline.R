# Figure contract: quantitative grid; show observed disease-group imbalance
# and separate cohort-adjusted contrasts within Nw/Ow/Ob. No causal claim.
# R-only drawing; existing data and model results are read-only.
library(ggplot2)
library(patchwork)
args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=','',args[grepl('^--file=',args)][1])
out <- dirname(normalizePath(script,winslash='/'))
root <- dirname(dirname(out))
meta_file <- file.path(root,'02.processed_data/02_aligned_data/metadata_aligned_20260208.csv')
m <- read.csv(meta_file,check.names=FALSE,na.strings=c('NA','','NaN'))
pairs <- read.delim(file.path(out,'tables/cohort_adjusted_pairwise.tsv'),check.names=FALSE)
omni <- read.delim(file.path(out,'tables/cohort_adjusted_omnibus.tsv'),check.names=FALSE)
stopifnot(nrow(m)==4385,!anyDuplicated(m$Sample.Source.ID.processed))
write_data <- function(d,name) write.table(d,file.path(out,'data_for_figures',paste0(name,'.tsv')),sep='\t',quote=TRUE,row.names=FALSE,na='NA')
m$BMI_display <- factor(m$BMI_Group1,levels=c('Normal','Ow','Ob'),labels=c('Nw','Ow','Ob'))
m$Disease_display <- factor(m$Disease,levels=c('CTRL','AA','CRC'),labels=c('CTRL','Adenoma','CRC'))
m$BMI_order <- as.integer(m$BMI_display)
m$Disease_order <- as.integer(m$Disease_display)
colors <- c(CTRL='#5785A4',Adenoma='#C39A42',CRC='#B76B72')
theme_set(theme_classic(base_size=12,base_family='Arial')+
 theme(strip.background=element_blank(),strip.text=element_text(face='bold',size=12),
 legend.position='none',axis.title=element_text(size=12),axis.text=element_text(size=10.5,color='black'),
 plot.title=element_text(size=12,face='bold'),plot.subtitle=element_text(size=10),
 plot.caption=element_text(size=9,hjust=0),plot.margin=margin(7,10,6,7)))

# Raw observations are not altered by cohort adjustment. Missing values are
# excluded per outcome, with counts retained in the figure data and labels.
age <- m[!is.na(m$Age),c('Sample.Source.ID.processed','Study.name','BMI_display','Disease_display','BMI_order','Disease_order','Age')]
age <- age[order(age$BMI_order,age$Disease_order,age$Age),]
age_counts <- aggregate(Age~BMI_display+Disease_display,data=age,FUN=length)
names(age_counts)[3] <- 'N_available'
age_counts$label <- paste0('n=',age_counts$N_available)
age_counts$label_y <- max(age$Age)+4
age$N_available <- age_counts$N_available[match(paste(age$BMI_display,age$Disease_display),paste(age_counts$BMI_display,age_counts$Disease_display))]
age$label <- paste0('n=',age$N_available)
age$label_y <- max(age$Age)+4
age$BMI_levels <- 'Nw|Ow|Ob'; age$Disease_levels <- 'CTRL|Adenoma|CRC'
write_data(age,'Fig_clinical_raw_age')
age_plot <- ggplot(age,aes(Disease_display,Age,fill=Disease_display))+
 geom_boxplot(width=.58,outlier.shape=16,outlier.size=.7,linewidth=.45)+
 geom_text(data=age_counts,aes(x=Disease_display,y=label_y,label=label),inherit.aes=FALSE,size=3.1)+
 facet_wrap(~BMI_display,nrow=1)+scale_fill_manual(values=colors)+
 scale_y_continuous(expand=expansion(mult=c(.02,.08)))+
 labs(title='a  Observed age distributions',subtitle='All available observations; boxes show median and IQR',x=NULL,y='Age (years)')+
 theme(axis.text.x=element_text(angle=35,hjust=1))

sex <- do.call(rbind,lapply(levels(m$BMI_display),function(b) do.call(rbind,lapply(levels(m$Disease_display),function(d) {
 z <- m[m$BMI_display==b & m$Disease_display==d,,drop=FALSE]
 valid <- z$Sex[!is.na(z$Sex)]
 n <- length(valid); count <- sum(valid=='Male'); ph <- count/n
 # Wilson 95% interval for descriptive observed proportions; not cohort adjusted.
 z95 <- qnorm(.975); den <- 1+z95^2/n
 ctr <- (ph+z95^2/(2*n))/den
 half <- z95*sqrt(ph*(1-ph)/n+z95^2/(4*n^2))/den
 data.frame(BMI_display=b,Disease_display=d,BMI_order=match(b,levels(m$BMI_display)),
 Disease_order=match(d,levels(m$Disease_display)),Male_N=count,N_available=n,N_total=nrow(z),N_missing=nrow(z)-n,
 Male_percent=100*ph,CI95_low=100*(ctr-half),CI95_high=100*(ctr+half),
 label=sprintf('%.1f%%\n%d/%d',100*ph,count,n),label_y=100*(ctr+half)+4)
}))))
sex$BMI_display <- factor(sex$BMI_display,levels=levels(m$BMI_display))
sex$Disease_display <- factor(sex$Disease_display,levels=levels(m$Disease_display))
sex$BMI_levels <- 'Nw|Ow|Ob';sex$Disease_levels <- 'CTRL|Adenoma|CRC'
write_data(sex,'Fig_clinical_raw_sex')
age$Panel <- 'observed_age';sex$Panel <- 'observed_male_proportion'
all_columns <- union(names(age),names(sex))
for(nm in setdiff(all_columns,names(age))) age[[nm]] <- NA
for(nm in setdiff(all_columns,names(sex))) sex[[nm]] <- NA
write_data(rbind(age[,all_columns],sex[,all_columns]),'Fig_clinical_raw_distributions')
sex_plot <- ggplot(sex,aes(Disease_display,Male_percent,color=Disease_display))+
 geom_errorbar(aes(ymin=CI95_low,ymax=CI95_high),width=.13,linewidth=.7)+geom_point(size=2.8)+
 geom_text(aes(y=label_y,label=label),color='black',size=3.05,lineheight=.95,vjust=0)+
 facet_wrap(~BMI_display,nrow=1)+scale_color_manual(values=colors)+
 scale_y_continuous(limits=c(0,88),breaks=c(0,20,40,60,80),expand=expansion(mult=c(.01,.02)))+
 labs(title='b  Observed male proportions',subtitle='Points and Wilson 95% CI; labels show % and male/available n',x=NULL,y='Male (%)')+
 theme(axis.text.x=element_text(angle=35,hjust=1))
box_stats <- ggplot_build(age_plot)$data[[1]]
box_stats$outliers <- vapply(box_stats$outliers,function(x) paste(x,collapse='|'),character(1))
box_stats$BMI_display <- levels(m$BMI_display)[as.integer(box_stats$PANEL)]
box_stats$Disease_display <- levels(m$Disease_display)[as.integer(box_stats$x)]
write_data(box_stats,'Fig_clinical_raw_age_boxplot_statistics')
raw_plot <- age_plot/sex_plot + plot_annotation(caption='Raw distributions use all available data. Cohort-adjusted estimates are shown separately.')

# Exact model results plus factor/annotation ordering are saved as plotting data.
adj <- pairs[pairs$Model=='cohort_adjusted',]
adj$BMI_display <- factor(adj$BMI,levels=c('Normal','Ow','Ob'),labels=c('Nw','Ow','Ob'))
adj$Comparison <- paste(ifelse(adj$Group_1=='AA','Adenoma',adj$Group_1),'vs',ifelse(adj$Group_2=='AA','Adenoma',adj$Group_2))
comparison_levels <- c('CRC vs Adenoma','CRC vs CTRL','Adenoma vs CTRL')
adj$Comparison <- factor(adj$Comparison,levels=comparison_levels)
adj$Comparison_order <- as.integer(adj$Comparison)
adj$Comparison_levels <- paste(comparison_levels,collapse='|')
adj$BMI_order <- as.integer(adj$BMI_display)
adj$BMI_levels <- 'Nw|Ow|Ob'
key <- paste(omni$BMI,omni$Variable,omni$Model)
i <- match(paste(adj$BMI,adj$Variable,adj$Model),key)
adj$Overall_P <- omni$P[i];adj$Overall_BH_q <- omni$BH_q[i]
adj$N_CTRL <- omni$N_CTRL[i];adj$N_AA <- omni$N_AA[i];adj$N_CRC <- omni$N_CRC[i]
adj$Cohorts <- omni$Cohorts[i]
adj$Facet_label <- paste0(adj$BMI_display,'\noverall q=',formatC(adj$Overall_BH_q,format='g',digits=2),'\nn=',adj$N)
adj$Pairwise_q_label <- paste0('q=',formatC(adj$BH_q,format='g',digits=2))
adj <- adj[order(adj$Variable,adj$BMI_order,-adj$Comparison_order),]
write_data(adj,'Fig_clinical_cohort_adjusted_forest')
forest <- function(variable,title,xlab,logscale=FALSE) {
 d <- adj[adj$Variable==variable,]
 strips <- setNames(unique(d$Facet_label),unique(as.character(d$BMI_display)))
 p <- ggplot(d,aes(x=Estimate,y=Comparison))+
  geom_vline(xintercept=if(logscale) 1 else 0,linetype=2,color='grey55',linewidth=.45)+
  geom_segment(aes(x=CI95_low,xend=CI95_high,yend=Comparison),linewidth=.7,color='#5785A4')+
  geom_point(size=2.6,color='#315D79')+
  facet_wrap(~BMI_display,nrow=1,labeller=labeller(BMI_display=strips))+
  labs(title=title,x=xlab,y=NULL)+theme(axis.text.y=element_text(size=10),strip.text=element_text(size=10.5),panel.spacing=grid::unit(10,'pt'))
 if(logscale) p <- p+scale_x_log10(breaks=c(.5,1,2,4),limits=c(.3,4)) else p <- p+scale_x_continuous(breaks=c(-2,0,4,8),limits=c(-3,9))
 p
}
age_forest <- forest('age','a  Age differences after cohort adjustment','Mean age difference (years)')
sex_forest <- forest('male','b  Sex differences after cohort adjustment','Male odds ratio (log scale)',TRUE)
forest_plot <- age_forest/sex_forest+plot_annotation(caption=paste('Points and nominal 95% CI; effects are first group relative to second group.',
 'Overall q: BH across 6 tests. Models use cohorts with within-cohort information.',sep='\n'))

save_figure <- function(p,stem,width,height) {
  grDevices::cairo_pdf(file.path(out,'figures',paste0(stem,'.pdf')),width=width,height=height,family='Arial')
 print(p);dev.off()
 svglite::svglite(file.path(out,'figures',paste0(stem,'.svg')),width=width,height=height)
 print(p);dev.off()
 ggsave(file.path(out,'figures',paste0(stem,'.png')),plot=p,device=ragg::agg_png,width=width,height=height,units='in',dpi=300,bg='white')
}
save_figure(raw_plot,'Fig_clinical_raw_distributions',7.2,6.0)
save_figure(forest_plot,'Fig_clinical_cohort_adjusted_forest',8.1,5.9)
stopifnot(nrow(age)==4180,nrow(sex)==9,nrow(adj)==18)
cat('Created two figures in PDF and PNG, with exact plotting data.\n')




