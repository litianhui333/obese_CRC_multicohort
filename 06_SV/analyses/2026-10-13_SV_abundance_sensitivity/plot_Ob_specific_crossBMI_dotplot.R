.libPaths(c("D:/R/library",.libPaths()))
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(svglite);library(ragg)})
setDTthreads(2L)
ROOT<-"H:/My Drive/proj_ob_AA_CRC/analyses/2026-10-13_SV_abundance_sensitivity"
STEM<-"Fig_Ob_specific_abundanceAdjusted_crossBMI"
TARGET<-"Ob"
DIRECTION_ORDER<-"direction_order" %in% commandArgs(TRUE)
HORIZONTAL<-"horizontal" %in% commandArgs(TRUE)
BMI_LEVELS<-c("Nw","Ow","Ob")
COLORS<-c(Nw="#4C78A8",Ow="#8C8C8C",Ob="#D55E00")
OFFSETS<-c(Nw=.23,Ow=0,Ob=-.23)
inputs<-file.path(ROOT,"tables",c("formal_effect_specific_candidate_reporting.tsv","run_model_results_long.tsv"))
hashes<-tools::md5sum(inputs)
existing<-list.files(file.path(ROOT,"figures"),full.names=TRUE)
protected<-existing[!startsWith(basename(existing),STEM)]
old_hash<-tools::md5sum(protected)
z<-fread(inputs[1])[BMI_enriched_group==TARGET]
stopifnot(nrow(z)==31L,z[,sum(SV_type=="dSGV")]==21L,z[,sum(SV_type=="vSGV")]==10L,
 all(z$P_difference_1<.05),all(z$P_difference_2<.05),!anyDuplicated(z$BMI_SV_ID))
models<-fread(inputs[2])[Model=="abundance_adjusted" & BMI_SV_ID %in% z$BMI_SV_ID & BMI %in% BMI_LEVELS]
stopifnot(nrow(models)==93L,!anyDuplicated(models[,.(BMI_SV_ID,BMI)]),all(is.finite(models$Disease_Beta)),
 all(is.finite(models$Disease_SE)&models$Disease_SE>0),all(is.na(models$Error)|!nzchar(as.character(models$Error))))
# Cross-check every displayed effect against the formal candidate reporting table.
for(g in BMI_LEVELS){a<-models[BMI==g][match(z$BMI_SV_ID,BMI_SV_ID)];
 stopifnot(max(abs(a$Disease_Beta-z[[paste0("Disease_Beta_",g)]]))<1e-12,
 max(abs(a$Disease_SE-z[[paste0("Disease_SE_",g)]]))<1e-12,
 max(abs(a$Disease_P-z[[paste0("Disease_P_",g)]]))<1e-12)}
z[,Display_species:=Species]
z[grepl("^[A-Z][a-z]+ [a-z][a-z-]+$",Species),Display_species:=sub("^([A-Z])[a-z]+ ","\\1. ",Species)]
# Coordinate labels follow the compact heatmap convention; internal IDs stay in source data only.
dash<-"-" # ASCII range separator prevents locale escape artifacts in plotmath labels
z[,Region_label_full:=gsub(";","; ",gsub("_",dash,sub("^.*:","",SV_feature),fixed=TRUE),fixed=TRUE)]
z[,N_segments:=lengths(strsplit(sub("^.*:","",SV_feature),";",fixed=TRUE))]
compact_region<-function(s){p<-trimws(strsplit(s,";",fixed=TRUE)[[1]]);if(length(p)==1L)return(paste0(p," kb"));if(length(p)==2L)return(paste0(paste(p,collapse="; ")," kb"));nums<-as.numeric(unlist(regmatches(p,gregexpr("[0-9]+",p))));paste0(min(nums),dash,max(nums)," kb (",length(p)," segments)")}
z[,Region_label:=vapply(Region_label_full,compact_region,character(1))]
z[,Plot_label:=paste0(Display_species," ",sub("^t__","",SGB)," : ",Region_label)]
quoted<-function(s)encodeString(s,quote="\"")
label_expr<-function(name,sgb,region){
 suffix<-paste0(sub("^t__","",sgb)," : ",region)
 if(grepl(" SGB[0-9]+(_group)?$",name))name<-sub(" SGB[0-9]+(_group)?$","",name)
 if(grepl(" sp ",name,fixed=TRUE)){
  bits<-strsplit(name," sp ",fixed=TRUE)[[1]]
  return(paste0("italic(",quoted(bits[1]),")~",quoted(paste0("sp ",bits[2]," ",suffix))))
 }
 paste0("italic(",quoted(name),")~",quoted(suffix))
}
z[,Plot_label_expression:=mapply(label_expr,Display_species,SGB,Region_label,USE.NAMES=FALSE)]
stopifnot(!any(grepl("ObSV_",z$Plot_label)),all(nzchar(z$Region_label)))
z[,Reference_ID:=sub(":.*$","",SV_feature)]
z[,SV_intervals_kb:=sub("^.*:","",SV_feature)]
z[,Ob_direction:=fifelse(Disease_Beta_Ob>0,"Higher in Ob CRC","Lower in Ob CRC")]
z[,Direction_order:=fifelse(Disease_Beta_Ob>0,1L,2L)]
if(DIRECTION_ORDER)setorder(z,SV_type,Direction_order,Species,SGB,BMI_SV_ID)else setorder(z,SV_type,Species,SGB,BMI_SV_ID)
z[,Row_order:=seq_len(.N),by=SV_type]
z[,N_rows:=.N,by=SV_type]
z[,Y_row:=N_rows-Row_order+1]
keys<-z[,.(BMI_SV_ID,Species,Display_species,Plot_label,Plot_label_expression,Region_label,Region_label_full,N_segments,Reference_ID,SV_intervals_kb,Ob_direction,Direction_order,Row_order,N_rows,Y_row,
 P_difference_1,P_difference_2,Other_BMI_1,Other_BMI_2,P_joint,Q_BH)]
d<-merge(models,keys,by="BMI_SV_ID",all.x=TRUE,sort=FALSE)
d[,BMI_order:=match(BMI,BMI_LEVELS)]
d[,BMI_color:=unname(COLORS[BMI])]
d[,Y_plot:=Y_row+unname(OFFSETS[BMI])]
d[,X_plot:=Row_order-unname(OFFSETS[BMI])]
d[,Plot_orientation:=if(HORIZONTAL)"SV on x; effect on y; labels 45 degrees"else"effect on x; SV on y"]
d[,CI_low:=Disease_Beta-qnorm(.975)*Disease_SE]
d[,CI_high:=Disease_Beta+qnorm(.975)*Disease_SE]
d[,Association:=fifelse(Disease_P<.01,"P < 0.01","P >= 0.01")]
d[,SV_panel:=fifelse(SV_type=="dSGV","dSV","vSV")]
d[,X_axis:=fifelse(SV_type=="dSGV","Log odds of deletion (CRC vs CTRL)","vSV coefficient (CRC vs CTRL)")]
setorder(d,SV_type,Row_order,BMI_order)
write_tsv<-function(a,area,name)fwrite(a,file.path(ROOT,area,name),sep="\t",na="NA")
make_panel<-function(type){a<-copy(d[SV_type==type]);rows<-unique(a[,.(Y_row,Plot_label,Plot_label_expression)]);setorder(rows,Y_row)
 boundary<-if(DIRECTION_ORDER && uniqueN(a$Direction_order)==2L)min(a[Direction_order==1L,Y_row])-.5 else numeric()
 ggplot(a,aes(Disease_Beta,Y_plot,color=BMI))+
  geom_hline(yintercept=sort(unique(a$Y_row))+.5,color="#EEEEEE",linewidth=.25)+
  geom_hline(yintercept=boundary,color="#999999",linewidth=.6)+
  geom_vline(xintercept=0,linetype="dashed",color="#666666",linewidth=.35)+
  geom_segment(aes(x=CI_low,xend=CI_high,yend=Y_plot),linewidth=.5,alpha=.85)+
  geom_point(aes(shape=Association),size=2.2,stroke=.75)+
  scale_color_manual(values=COLORS,breaks=BMI_LEVELS,labels=c("Nw","Ow","Ob"),name="BMI group")+
  scale_shape_manual(values=c("P < 0.01"=16,"P >= 0.01"=1),breaks=c("P < 0.01","P >= 0.01"),name="CRC association")+
  scale_y_continuous(breaks=rows$Y_row,labels=parse(text=rows$Plot_label_expression),limits=c(.45,max(a$Y_row)+.55),expand=expansion(mult=0))+
  scale_x_continuous(breaks=scales::breaks_pretty(5),expand=expansion(mult=.06))+
  labs(title=paste0(if(type=="dSGV")"dSV"else"vSV"," (",uniqueN(a$BMI_SV_ID)," candidates)"),x=unique(a$X_axis),y=NULL)+
  guides(color=guide_legend(order=1,override.aes=list(shape=16,linewidth=0)),shape=guide_legend(order=2,override.aes=list(color="#333333")))+
  theme_classic(base_size=12,base_family="Arial")+
  theme(axis.text.y=element_text(size=10,color="#222222"),axis.text.x=element_text(size=10,color="#222222"),
   axis.title.x=element_text(size=11),plot.title=element_text(size=12,face="bold"),
   legend.position="top",legend.text=element_text(size=10),legend.title=element_text(size=10),
   plot.margin=margin(8,10,8,8))
}
export<-function(p,stem,w,h){
 cairo_pdf(file.path(ROOT,"figures",paste0(stem,".pdf")),width=w,height=h,family="Arial");print(p);dev.off()
 svglite(file.path(ROOT,"figures",paste0(stem,".svg")),width=w,height=h);print(p);dev.off()
 agg_png(file.path(ROOT,"figures",paste0(stem,".png")),width=w,height=h,units="in",res=300,background="white");print(p);dev.off()
}
# Changing row order must not change any numerical estimate or confidence limit.
prior_path<-file.path(ROOT,"data_for_figures",paste0(STEM,".tsv"))
if(DIRECTION_ORDER && file.exists(prior_path)){
 prior<-fread(prior_path);keys_compare<-c("BMI_SV_ID","BMI")
 current_check<-copy(d);setorderv(prior,keys_compare);setorderv(current_check,keys_compare)
 for(field in c("Disease_Beta","Disease_SE","Disease_P","CI_low","CI_high"))stopifnot(isTRUE(all.equal(prior[[field]],current_check[[field]],tolerance=1e-12)))
 stopifnot(identical(prior$BMI_SV_ID,current_check$BMI_SV_ID),identical(prior$BMI,current_check$BMI))
}
if(DIRECTION_ORDER)stopifnot(all(z[,all(diff(Direction_order)>=0),by=SV_type]$V1))
make_horizontal<-function(type){
 a<-copy(d[SV_type==type]);rows<-unique(a[,.(Row_order,Plot_label_expression)]);setorder(rows,Row_order)
 boundary<-if(DIRECTION_ORDER && uniqueN(a$Direction_order)==2L)max(a[Direction_order==1L,Row_order])+.5 else numeric()
 ggplot(a,aes(X_plot,Disease_Beta,color=BMI))+
  geom_vline(xintercept=sort(unique(a$Row_order))+.5,color="#EEEEEE",linewidth=.25)+
  geom_vline(xintercept=boundary,color="#999999",linewidth=.6)+
  geom_hline(yintercept=0,linetype="dashed",color="#666666",linewidth=.35)+
  geom_segment(aes(xend=X_plot,y=CI_low,yend=CI_high),linewidth=.5,alpha=.85)+
  geom_point(aes(shape=Association),size=2.2,stroke=.75)+
  scale_color_manual(values=COLORS,breaks=BMI_LEVELS,name="BMI group")+
  scale_shape_manual(values=c("P < 0.01"=16,"P >= 0.01"=1),breaks=c("P < 0.01","P >= 0.01"),name="CRC association")+
  scale_x_continuous(breaks=rows$Row_order,labels=parse(text=rows$Plot_label_expression),limits=c(.45,max(a$Row_order)+.55),expand=expansion(mult=0))+
  scale_y_continuous(breaks=scales::breaks_pretty(5),expand=expansion(mult=.06))+
  labs(title=paste0(if(type=="dSGV")"dSV"else"vSV"," (",uniqueN(a$BMI_SV_ID)," candidates)"),x=NULL,y=unique(a$X_axis))+
  guides(color=guide_legend(order=1,override.aes=list(shape=16,linewidth=0)),shape=guide_legend(order=2,override.aes=list(color="#333333")))+
  theme_classic(base_size=12,base_family="Arial")+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=10,color="#222222"),
   axis.text.y=element_text(size=10,color="#222222"),axis.title.y=element_text(size=11),
   plot.title=element_text(size=12,face="bold"),legend.position="top",legend.text=element_text(size=10),legend.title=element_text(size=10),
   plot.margin=margin(8,10,8,170))
}
if(HORIZONTAL)make_panel<-make_horizontal
pd<-make_panel("dSGV");pv<-make_panel("vSGV")
combined<-(pd|pv)+plot_layout(guides="collect",widths=if(HORIZONTAL)c(21,10)else c(1,1))+plot_annotation(
 title="Putative obesity-specific SVs after abundance adjustment",
 subtitle=if(HORIZONTAL)"CRC effects in Nw, Ow and Ob; higher then lower Ob effects from left to right"else if(DIRECTION_ORDER)"CRC effects in Nw, Ow and Ob; rows grouped by Ob effect: higher above, lower below"else"CRC effects in Nw, Ow and Ob; points and 95% Wald confidence intervals",
 caption="Selected by nominal cross-BMI effect differences; none survived BH correction.",
 tag_levels="a",theme=theme(plot.title=element_text(family="Arial",size=14,face="bold"),
 plot.subtitle=element_text(family="Arial",size=11),plot.caption=element_text(family="Arial",size=10,hjust=0)))
combined<-combined & theme(legend.position="top")
export(combined,STEM,if(HORIZONTAL)18 else 14,if(HORIZONTAL)7.5 else 8.4)
export(pd,paste0(STEM,"_dSV"),if(HORIZONTAL)11.8 else 7.5,if(HORIZONTAL)7.5 else 8.4)
pv_export<-if(HORIZONTAL)pv+theme(legend.box="vertical",legend.text=element_text(size=9),legend.title=element_text(size=9))else pv
export(pv_export,paste0(STEM,"_vSV"),if(HORIZONTAL)6.2 else 7.5,if(HORIZONTAL)7.5 else 5.5)
write_tsv(d,"data_for_figures",paste0(STEM,".tsv"))
write_tsv(d[SV_type=="dSGV"],"data_for_figures",paste0(STEM,"_dSV.tsv"))
write_tsv(d[SV_type=="vSGV"],"data_for_figures",paste0(STEM,"_vSV.tsv"))
write_tsv(z,"tables","Ob_specific_crossBMI_candidate_locator.tsv")
stopifnot(nrow(d)==93L,uniqueN(d$BMI_SV_ID)==31L,all(d[,uniqueN(BMI),by=BMI_SV_ID]$V1==3L),
 identical(unname(hashes),unname(tools::md5sum(inputs))),identical(unname(old_hash),unname(tools::md5sum(protected))))
qa<-data.table(Check=c("31 candidates","21 dSV and 10 vSV","93 points; 3 BMI groups each","Matches adjusted model/reporting coefficients",
 "No input modified","Unrelated figures preserved"),Pass=TRUE)
write_tsv(qa,"tables","Ob_specific_crossBMI_dotplot_QA.tsv")
legend<-c("# Abundance-adjusted cross-BMI effect dot plot","", "Status: exploratory", "",
 "Conclusion/evidence contract: directly compare Nw, Ow and Ob CRC effects for the 31 retained nominal Ob-targeted effect-difference candidates. Quantitative two-panel grid; all candidates and all three effects are shown. No models are refitted.",
 "Each position represents one SV, labelled by italicized species name, upright SGB and genomic interval in kb. Internal candidate IDs are retained only in source data. Genus names are italicized but sp./strain strings and SGB identifiers remain upright. Two-segment SVs list both component intervals; SVs with more than two segments use the heatmap compact convention: bounding interval plus explicit segment count; this bounding interval is not a single contiguous deletion. Region_label_full and SV_feature preserve every original component interval. The locator table retains full representative-genome ID and complete SV intervals. Row order is recorded in the exact plotting table. With direction_order: positive Ob CRC effects precede negative effects, with a horizontal separator; within each block sort species, SGB and candidate ID. Without that argument, preserve the earlier alphabetical order. No overlapping SVs are merged.",
 "Colors identify BMI strata. Each point is the CRC-versus-CTRL coefficient from the abundance-adjusted mixed model including age, sex, study random intercept and standardized log1p parent-SGB relative abundance. Model sample sizes and QC flags remain in exact plotting data. Nw/Ow/Ob fits use their own eligible complete-case participants.",
 "Intervals are normal-approximation 95% Wald intervals: beta +/- 1.959964 SE. dSV effects are log odds of deletion; vSV effects are linear coefficients in raw vSV units. The two panels have separate x-axis units and ranges. No clipping or effect standardization is applied.",
 "Filled/open points denote within-stratum CRC association P <0.01 / >=0.01 from the stored model outputs, not significance of BMI specificity. Exact coefficients, SEs, association P values and two Ob-versus-other effect-difference P values are provided in the plotting table. Hollow points are retained regardless of significance.",
 "Candidate inclusion is taken unchanged from formal_effect_specific_candidate_reporting.tsv: 31 nominal Ob candidates (21 dSV,10 vSV), both target-versus-other effect-difference P values <0.05. None survives the reported cross-BMI BH correction. These are putative/exploratory obesity-specific candidates, not confirmed specificity. Selection is restricted to the previous 260-candidate sensitivity library, not genome-wide rediscovery.",
 "Exports: combined 14 x8.4 in, individual dSV 7.5 x8.4 in, vSV 7.5 x5.5 in; PDF/SVG editable text and 300 dpi PNG. For manuscript assembly use individual vector panels. Every figure has its exact plotted rows, grouping, offsets, confidence limits, colors and ordering in data_for_figures.")
writeLines(legend,file.path(ROOT,"Ob_specific_crossBMI_dotplot_legend.md"))
append_unique<-function(path,entries){old<-if(file.exists(path))readLines(path,warn=FALSE)else character();new<-setdiff(entries,old);if(length(new))cat(paste0(new,collapse="\n"),"\n",file=path,append=TRUE,sep="")}
append_unique(file.path(ROOT,"input_files.txt"),inputs)
outputs<-c(file.path(ROOT,"figures",as.vector(outer(c(STEM,paste0(STEM,"_dSV"),paste0(STEM,"_vSV")),c(".pdf",".svg",".png"),paste0))),
 file.path(ROOT,"data_for_figures",paste0(c(STEM,paste0(STEM,"_dSV"),paste0(STEM,"_vSV")),".tsv")),
 file.path(ROOT,"tables",c("Ob_specific_crossBMI_candidate_locator.tsv","Ob_specific_crossBMI_dotplot_QA.tsv")),file.path(ROOT,"Ob_specific_crossBMI_dotplot_legend.md"))
append_unique(file.path(ROOT,"output_files.txt"),outputs)
cat("Completed: 31 SVs, 93 adjusted effects and confidence intervals.\n")
print(d[,.(SVs=uniqueN(BMI_SV_ID),Points=.N,N_P_lt_01=sum(Disease_P<.01)),by=.(SV_type,BMI)])
if(DIRECTION_ORDER){
 write_tsv(z[,.(N_SVs=.N),by=.(SV_type,Ob_direction)],"tables","Ob_specific_crossBMI_direction_order_counts.tsv")
 cat("Direction grouping checked; all beta/SE/P/CI values unchanged.\n")
 print(z[,.(N_SVs=.N),by=.(SV_type,Ob_direction)])
}
if(HORIZONTAL){
 cat("\nHorizontal layout: SVs along x, effects/95% CI along y; labels rotated 45 degrees. Combined panel plotting widths have 21:10 ratio. Export dimensions: combined 18x7.5 in; dSV 11.8x7.5; vSV 6.2x7.5. These supersede the earlier vertical-layout dimensions above. Within each panel higher Ob effects precede lower effects left-to-right. Reproduce with Rscript plot_Ob_specific_crossBMI_dotplot.R direction_order horizontal.\n",file=file.path(ROOT,"Ob_specific_crossBMI_dotplot_legend.md"),append=TRUE)
}