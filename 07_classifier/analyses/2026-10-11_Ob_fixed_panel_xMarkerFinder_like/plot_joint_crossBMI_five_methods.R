#!/usr/bin/env Rscript
# Plot-only subset of saved Ob-trained cross-BMI ROC results; no model refitting.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
ROOT <- "H:/My Drive/proj_ob_AA_CRC/analyses/2026-10-11_Ob_fixed_panel_xMarkerFinder_like"
IN_STEM <- "restricted_ObAssoc_L1L2_crossBMI_otherAlgorithms_"
OUT_STEM <- "restricted_ObAssoc_L1L2_crossBMI_joint_five_methods_ROC"
METHODS <- c("Random forest", "LASSO", "Elastic Net", "Radial SVM", "XGBoost")
GROUPS <- c("Ob_test", "Nw", "Ow")
LABEL <- c(Ob_test="Obese", Nw="Nw", Ow="Ow")
COLORS <- c(Obese="#C45144", Nw="#497CAE", Ow="#E4A351")
curve_path <- file.path(ROOT, "data_for_figures", paste0(IN_STEM, "Fig_transfer_ROC_representative_curves.tsv"))
auc_path <- file.path(ROOT, "tables", paste0(IN_STEM, "AUC_summary.tsv"))
nested_curve_path <- file.path(ROOT, "data_for_figures",
                               "restricted_ObAssoc_L1L2_nestedCV_Fig_ROC_representative_repeat.tsv")
nested_auc_path <- file.path(ROOT, "tables", "restricted_ObAssoc_L1L2_nestedCV_summary.tsv")
stopifnot(all(file.exists(c(curve_path, auc_path, nested_curve_path, nested_auc_path))))
stopifnot(requireNamespace("svglite", quietly=TRUE), requireNamespace("ragg", quietly=TRUE))
# Obese uses the requested 477-sample nested-CV figure; Nw/Ow retain transfer curves.
transfer_curves <- fread(curve_path)[Panel=="joint" & Model %in% METHODS &
                                      BMI_evaluation %in% c("Nw","Ow"),
  .(Panel,Model,BMI_evaluation,Representative_repeat,FPR,TPR,Mean_AUC,N)]
transfer_curves[, Source := "Cross-BMI transfer"]
nested_curves <- fread(nested_curve_path)[Record_type=="ROC point" & Panel=="joint" &
                                          Model %in% METHODS,
  .(Panel,Model,Representative_repeat,FPR,TPR,Mean_AUC,N)]
nested_curves[, BMI_evaluation := "Ob_test"]
nested_curves[, Source := "Ob nested CV"]
curves <- rbindlist(list(nested_curves,transfer_curves),use.names=TRUE,fill=TRUE)
transfer_auc <- fread(auc_path)[Panel=="joint" & Model %in% METHODS &
                                 BMI_evaluation %in% c("Nw","Ow"),
  .(Panel,Model,BMI_evaluation,Mean_AUC,N_repeats,N)]
transfer_auc[, Source := "Cross-BMI transfer"]
nested_auc <- fread(nested_auc_path)[Panel=="joint" & Model %in% METHODS,
  .(Panel,Model,BMI_evaluation="Ob_test",Mean_AUC,N_repeats,N=N_samples)]
nested_auc[, Source := "Ob nested CV"]
auc <- rbindlist(list(nested_auc,transfer_auc),use.names=TRUE,fill=TRUE)
stopifnot(nrow(auc)==15L, uniqueN(auc[, .(Model,BMI_evaluation)])==15L)
stopifnot(all(auc$N_repeats==5L), all(is.finite(curves$FPR)), all(is.finite(curves$TPR)))
stopifnot(all(curves$FPR>=0 & curves$FPR<=1), all(curves$TPR>=0 & curves$TPR<=1))
check <- merge(curves[, .(Model,BMI_evaluation,Mean_AUC,N)], auc[, .(Model,BMI_evaluation,Mean_AUC,N)],
               by=c("Model","BMI_evaluation"), suffixes=c("_curve","_summary"))
stopifnot(nrow(check)==nrow(curves))
stopifnot(all(abs(check$Mean_AUC_curve-check$Mean_AUC_summary)<1e-10), all(check$N_curve==check$N_summary))
curves[, Method := factor(Model, levels=METHODS)]
curves[, BMI := factor(LABEL[BMI_evaluation], levels=unname(LABEL[GROUPS]))]
setorder(curves, Method, BMI, FPR, TPR)
auc[, Method := factor(Model, levels=METHODS)]
auc[, BMI := factor(LABEL[BMI_evaluation], levels=unname(LABEL[GROUPS]))]
auc[, Label := sprintf("%s  %.3f", LABEL[BMI_evaluation], Mean_AUC)]
auc[, Label_x := 0.98]
auc[, Label_y := c(Ob_test=0.25, Nw=0.16, Ow=0.07)[BMI_evaluation]]
setorder(auc, Method, BMI)
points <- curves[, .(Record_type="ROC point",Panel,Model,Method,BMI_evaluation,BMI,
                     Representative_repeat,FPR,TPR,Mean_AUC,N,Source,Label=NA_character_)]
labels <- auc[, .(Record_type="AUC annotation",Panel,Model,Method,BMI_evaluation,BMI,
                   Representative_repeat=NA_integer_,FPR=Label_x,TPR=Label_y,Mean_AUC,N,Source,Label)]
plot_data <- rbindlist(list(points,labels), use.names=TRUE, fill=TRUE)
data_path <- file.path(ROOT,"data_for_figures",paste0(OUT_STEM,".tsv"))
fwrite(plot_data,data_path,sep="\t",na="NA",quote=TRUE)
make_panel <- function(method_name, first=FALSE) {
  d <- curves[Model==method_name]
  a <- auc[Model==method_name]
  ggplot(d,aes(FPR,TPR,colour=BMI,group=BMI)) +
    geom_abline(slope=1,intercept=0,linetype="dashed",colour="#B8B8B8",linewidth=0.35) +
    geom_step(data=d[BMI_evaluation=="Ob_test"],direction="hv",linewidth=0.75,alpha=0.96) +
    geom_step(data=d[BMI_evaluation %in% c("Nw","Ow")],direction="vh",linewidth=0.75,alpha=0.96) +
    geom_text(data=a,aes(x=Label_x,y=Label_y,label=Label,colour=BMI),
              inherit.aes=FALSE,hjust=1,size=2.45,fontface="bold",show.legend=FALSE) +
    scale_colour_manual(values=COLORS,drop=FALSE) +
    scale_x_continuous(breaks=c(0,0.5,1),limits=c(0,1),expand=c(0,0)) +
    scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1),expand=c(0,0)) +
    coord_equal() +
    labs(title=method_name,x="False-positive rate",y=if(first) "True-positive rate" else NULL) +
    theme_classic(base_size=9,base_family="Arial") +
    theme(axis.line=element_line(linewidth=0.35,colour="#333333"),
          axis.ticks=element_line(linewidth=0.3,colour="#333333"),
          axis.text=element_text(size=7.2,colour="#333333"),
          axis.title=element_text(size=8.3,colour="#333333"),
          plot.title=element_text(size=9.3,face="bold",hjust=0.5),
          legend.position="none",plot.margin=margin(5,5,5,5))
}
panels <- Map(function(m,i) make_panel(m,i==1L),METHODS,seq_along(METHODS))
fig <- wrap_plots(panels,nrow=1) +
  plot_annotation(title="Ob-trained species + SV model across BMI groups",
                  subtitle="Obese: nested CV; Nw/Ow: transfer evaluation | labels: five-repeat mean AUC",
                  theme=theme(plot.title=element_text(face="bold",size=13),
                              plot.subtitle=element_text(size=9,colour="#555555"),
                              plot.margin=margin(6,6,6,6)))
stem <- file.path(ROOT,"figures",OUT_STEM)
pdf_path <- paste0(stem,".pdf")
svg_path <- paste0(stem,".svg")
png_path <- paste0(stem,".png")
grDevices::cairo_pdf(pdf_path,width=11.1,height=3.4,family="Arial"); print(fig); dev.off()
svglite::svglite(svg_path,width=11.1,height=3.4); print(fig); dev.off()
ragg::agg_png(png_path,width=11.1,height=3.4,units="in",res=400); print(fig); dev.off()
endpoints <- curves[, .(Start=any(FPR==0 & TPR==0),End=any(FPR==1 & TPR==1)),
                    by=.(Model,BMI_evaluation)]
qc <- data.table(Check=c("Five methods and three BMI groups",
                         "Obese nested CV and Nw/Ow transfer sources",
                         "Five repeats per AUC label",
                         "ROC endpoints present","Plotting data and exports present"),
                 Passed=c(uniqueN(curves[,.(Model,BMI_evaluation)])==15L,
                          all(curves[BMI_evaluation=="Ob_test"]$N==477L) &&
                          all(auc[BMI_evaluation=="Ob_test"]$Source=="Ob nested CV") &&
                          all(curves[BMI_evaluation %in% c("Nw","Ow")]$Source=="Cross-BMI transfer"),
                          all(auc$N_repeats==5L),all(endpoints$Start & endpoints$End),
                          all(file.exists(c(data_path,pdf_path,svg_path,png_path)))))
fwrite(qc,file.path(ROOT,"tables",paste0(OUT_STEM,"_QC.tsv")),sep="\t")
stopifnot(all(qc$Passed))
cat("JOINT_CROSS_BMI_FIVE_METHODS_COMPLETE\n")

