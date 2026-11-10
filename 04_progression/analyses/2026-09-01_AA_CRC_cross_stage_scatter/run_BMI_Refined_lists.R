# Targeted rerun: authoritative Refined BMI lists; linear fit and Spearman statistics.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(ggrepel)})
base_dir <- "H:/My Drive/proj_ob_AA_CRC"
analysis_dir <- file.path(base_dir, "analyses/2026-09-01_AA_CRC_cross_stage_scatter")
list_dir <- file.path(base_dir, "analyses/2026-07-13_bmi_union_overlap_heatmaps/tables")
bmi_groups <- c("Normal", "Ow", "Ob")
p_cutoff <- 0.01
fdr_cutoff <- 0.1
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1, args[1] %in% c("full_list", "strict_in_list"))
selection_mode <- args[1]
output_stem <- paste0("D_BMI_Refined_", selection_mode, "_AA_CRC_spearman")
table_stem <- paste0("BMI_Refined_", selection_mode)
width_in <- 8.9
height_in <- 4.3
for (subdir in c("figures", "tables", "data_for_figures")) dir.create(file.path(analysis_dir, subdir), showWarnings=FALSE)
write_tsv <- function(x, path) write.table(x, path, sep="\t", quote=TRUE, row.names=FALSE, na="NA", fileEncoding="UTF-8")
all_data <- bind_rows(lapply(bmi_groups, function(g) {
  path <- file.path(list_dir, paste0("Feature_List_",g,"_Refined.csv"))
  d <- read.csv(path, stringsAsFactors=FALSE, check.names=FALSE)
  stopifnot(all(c("feature","clean_name","overlap","coef_crc","pval_crc","qval_crc","coef_aa","pval_aa") %in% names(d)), !anyDuplicated(d$feature))
  # Use the paired effect sizes supplied in each authoritative Refined list.
  d %>% transmute(SGB=feature, clean_name=clean_name, original_overlap=overlap,
    coef_CRC=coef_crc, pval_CRC=pval_crc, FDR_CRC=qval_crc, coef_AA=coef_aa, pval_AA=pval_aa,
    BMI_Group=g, input_file=basename(path), input_row=seq_len(nrow(d)),
    passes_raw_P=!is.na(pval_CRC) & pval_CRC < p_cutoff,
    passes_FDR=!is.na(FDR_CRC) & FDR_CRC < fdr_cutoff,
    selected=if(selection_mode=="full_list") TRUE else passes_raw_P | passes_FDR,
    finite_pair=is.finite(coef_AA) & is.finite(coef_CRC))
}))
source_check <- all_data %>% group_by(BMI_Group) %>% summarise(n_listed=n(),n_unique=n_distinct(SGB),
  n_passes_strict=sum(passes_raw_P | passes_FDR),n_selected=sum(selected),n_missing_effect=sum(!finite_pair),.groups="drop")
stopifnot(all(source_check$n_missing_effect==0))
membership <- all_data %>% group_by(SGB) %>% summarise(
  n_associated_BMI=sum(selected), associated_BMI=paste(BMI_Group[selected], collapse="+"), .groups="drop")
plot_data <- all_data %>% filter(selected, finite_pair) %>% left_join(membership, by="SGB") %>%
  mutate(BMI_Group=factor(BMI_Group, levels=bmi_groups), BMI_order=as.integer(BMI_Group),
    association_class=factor(paste(n_associated_BMI, ifelse(n_associated_BMI==1, "BMI group", "BMI groups")),
      levels=c("1 BMI group", "2 BMI groups", "3 BMI groups")),
    direction_concordant=sign(coef_AA)==sign(coef_CRC), label="") %>% arrange(BMI_Group, SGB)
old_labels <- read.delim(file.path(analysis_dir, "tables/scatter_labels_used.tsv"), stringsAsFactors=FALSE)
for (i in seq_len(nrow(old_labels))) {
  ix <- plot_data$SGB==old_labels$SGB[i] & as.character(plot_data$BMI_Group)==old_labels$BMI_Group[i]
  plot_data$label[ix] <- old_labels$label[i]
}
summary_table <- bind_rows(lapply(bmi_groups, function(g) {
  d <- filter(plot_data, BMI_Group==g)
  full <- filter(all_data, BMI_Group==g)
  test <- if(nrow(d)>=3) cor.test(d$coef_AA,d$coef_CRC,method="spearman",exact=FALSE) else NULL
  data.frame(BMI_Group=g, n_tested=nrow(full), n_raw_P=sum(full$passes_raw_P), n_FDR=sum(full$passes_FDR),
    n_selected=sum(full$selected), n_excluded_missing_effect=sum(full$selected & !full$finite_pair),
    n_taxa=nrow(d), spearman_rho=if(is.null(test)) NA_real_ else unname(test$estimate),
    spearman_p=if(is.null(test)) NA_real_ else test$p.value, direction_concordance=mean(d$direction_concordant))
}))
summary_table$spearman_FDR_across_3_BMI <- p.adjust(summary_table$spearman_p,method="BH")
format_p <- function(p) ifelse(p<0.0001, "P < 0.0001", sprintf("P = %.3g",p))
summary_table$annotation <- sprintf("\u03c1 = %.2f; %s\nn = %d; concordant = %.1f%%", summary_table$spearman_rho,
  format_p(summary_table$spearman_p),summary_table$n_taxa,100*summary_table$direction_concordance)
plot_data <- left_join(plot_data, summary_table %>% select(BMI_Group,spearman_rho,spearman_p,spearman_FDR_across_3_BMI,annotation), by="BMI_Group")
plot_data$BMI_Group <- factor(plot_data$BMI_Group,levels=bmi_groups)
annotation_data <- summary_table %>% mutate(BMI_Group=factor(BMI_Group, levels=bmi_groups), placement="facet strip")
strip_labels <- setNames(paste(summary_table$BMI_Group, summary_table$annotation, sep="\n"), summary_table$BMI_Group)
base_plot <- ggplot(plot_data,aes(coef_AA,coef_CRC)) +
  geom_hline(yintercept=0,linetype="dashed",colour="grey30",linewidth=0.35) +
  geom_vline(xintercept=0,linetype="dashed",colour="grey30",linewidth=0.35) +
  facet_wrap(~BMI_Group,nrow=1,scales="free",labeller=labeller(BMI_Group=strip_labels)) +
  scale_x_continuous(breaks=scales::breaks_pretty(n=3)) + scale_y_continuous(breaks=scales::breaks_pretty(n=4)) +
  labs(x="AA vs CTRL effect size",y="CRC vs CTRL effect size") +
  theme_classic(base_size=12,base_family="Arial") +
  theme(strip.background=element_blank(),strip.text=element_text(size=10.5,lineheight=1.15),
    legend.position="bottom",legend.title=element_text(size=11),legend.text=element_text(size=11),
    panel.spacing=grid::unit(0.7,"lines"),plot.margin=margin(5.5,8,5.5,5.5))
colours <- c("1 BMI group"="#E67E22","2 BMI groups"="#2D8BBB","3 BMI groups"="#27AE91")
shapes <- c("1 BMI group"=17,"2 BMI groups"=15,"3 BMI groups"=18)
label_layer <- geom_text_repel(data=filter(plot_data,label!=""),aes(label=label),seed=42,size=3.3,
  min.segment.length=0,box.padding=0.6,point.padding=0.2,max.overlaps=Inf,force=3,max.time=3,max.iter=100000)
plain <- base_plot + geom_smooth(method="lm",formula=y~x,se=TRUE,colour="black",fill="grey80",linewidth=0.65) + geom_point(colour="#2D8BBB",size=2.2,alpha=0.85) + label_layer
smooth <- base_plot + geom_smooth(method="lm",formula=y~x,se=TRUE,colour="black",fill="grey80",linewidth=0.65) +
  geom_point(aes(colour=association_class,shape=association_class),size=2.3,alpha=0.9) +
  scale_colour_manual(values=colours,drop=FALSE) + scale_shape_manual(values=shapes,drop=FALSE) +
  labs(colour="Refined-list membership in",shape="Refined-list membership in") + label_layer
save_plot <- function(p,stem) {
  grDevices::cairo_pdf(file.path(analysis_dir,"figures",paste0(stem,".pdf")),width=width_in,height=height_in,family="Arial")
  print(p); grDevices::dev.off()
  ggsave(file.path(analysis_dir,"figures",paste0(stem,".png")),p,width=width_in,height=height_in,dpi=300,bg="white")
  svglite::svglite(file.path(analysis_dir,"figures",paste0(stem,".svg")),width=width_in,height=height_in)
  print(p); grDevices::dev.off()
  write_tsv(plot_data,file.path(analysis_dir,"data_for_figures",paste0(stem,".tsv")))
  write_tsv(annotation_data,file.path(analysis_dir,"data_for_figures",paste0(stem,"_annotations.tsv")))
  # Each exported figure has its own exact lm curve and confidence-band data.
  fit_data <- ggplot_build(p)$data[[3]]
  fit_data$BMI_Group <- bmi_groups[as.integer(fit_data$PANEL)]
  write_tsv(fit_data,file.path(analysis_dir,"data_for_figures",paste0(stem,"_fit.tsv")))
}
save_plot(plain,paste0(output_stem,"_by_BMI_linear"))
save_plot(smooth,paste0(output_stem,"_by_BMI_color_linear"))
write_tsv(summary_table,file.path(analysis_dir,paste0("tables/",table_stem,"_spearman_summary.tsv")))
write_tsv(all_data %>% left_join(membership,by="SGB") %>% arrange(factor(BMI_Group,levels=bmi_groups),SGB),
  file.path(analysis_dir,paste0("tables/",table_stem,"_selection_audit.tsv")))
write_tsv(filter(plot_data,label!=""),file.path(analysis_dir,paste0("tables/",table_stem,"_labels_used.tsv")))
write_tsv(source_check,file.path(analysis_dir,paste0("tables/",table_stem,"_source_validation.tsv")))
capture.output(sessionInfo(),file=file.path(analysis_dir,paste0("tables/",table_stem,"_sessionInfo.txt")))
print(summary_table)






