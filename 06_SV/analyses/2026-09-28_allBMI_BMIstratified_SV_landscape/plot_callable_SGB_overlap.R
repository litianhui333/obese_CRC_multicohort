.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages(library(VennDiagram))
library(grid)
arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
root <- dirname(normalizePath(sub("^--file=", "", arg), winslash="/"))
input <- file.path(root,"tables/03_BMIstratified_retained_SGBs.tsv")
dat <- read.delim(input, stringsAsFactors=FALSE)
groups <- c("Normal","Ow","Ob")
sets <- setNames(lapply(groups,function(g) unique(dat$SGB[dat$BMI_Group==g & dat$Retained])),groups)
ids <- sort(unique(unlist(sets)))
membership <- data.frame(SGB=ids,Normal=ids %in% sets$Normal,Ow=ids %in% sets$Ow,Ob=ids %in% sets$Ob)
membership$Region <- apply(membership[,groups],1,function(x) paste(groups[x],collapse=" + "))
regions <- as.data.frame(table(membership$Region), stringsAsFactors=FALSE)
names(regions) <- c("Region","N_SGB")
n <- lengths(sets)
n12 <- length(intersect(sets[[1]],sets[[2]])); n13 <- length(intersect(sets[[1]],sets[[3]]))
n23 <- length(intersect(sets[[2]],sets[[3]])); n123 <- length(Reduce(intersect,sets))
stopifnot(identical(as.integer(n),c(45L,51L,55L)),n123==42L,length(ids)==60L,sum(regions$N_SGB)==60L)
dir.create(file.path(root,"data_for_figures"),showWarnings=FALSE)
stem <- "Fig_callable_SGB_overlap_allBMI"
write.table(membership,file.path(root,"data_for_figures",paste0(stem,"_membership.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
write.table(regions,file.path(root,"data_for_figures",paste0(stem,"_regions.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
summary <- data.frame(Group=groups,Total_SGB=as.integer(n),Shared_all_three=n123,Shared_fraction=n123/as.integer(n),Color=c("#5AAE61","#FF8C00","#9970AB"))
write.table(summary,file.path(root,"data_for_figures",paste0(stem,"_groups.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
draw <- function(){
 grid.newpage()
 grid.text("Callable SGB overlap across BMI groups",y=.96,gp=gpar(fontsize=15,fontface="bold"))
 grid.text("Common all-BMI SV collection",y=.91,gp=gpar(fontsize=11,col="#555555"))
 pushViewport(viewport(y=.49,height=.75,width=.94))
 v <- draw.triple.venn(area1=n[1],area2=n[2],area3=n[3],n12=n12,n23=n23,n13=n13,n123=n123,
  category=c("Normal (45)","Overweight (51)","Obese (55)"),fill=summary$Color,alpha=.34,
  col=summary$Color,lwd=1.5,cex=1.5,fontfamily="sans",cat.cex=1.1,cat.fontfamily="sans",
  cat.col="#333333",scaled=FALSE,euler.d=FALSE,ind=FALSE,margin=.13)
 grid.draw(v);popViewport()
 grid.text("42 of 60 SGBs shared by all three groups (70%)",y=.075,gp=gpar(fontsize=12,fontface="bold"))
 grid.text("Numbers denote SGBs passing group-specific callability thresholds; areas are schematic.",y=.025,gp=gpar(fontsize=8,col="#555555"))
}
svglite::svglite(
 file.path(root,"plots",paste0(stem,".svg")),
 width=6.5,height=6.5
)
draw();dev.off()
grDevices::cairo_pdf(
 file.path(root,"plots",paste0(stem,".pdf")),
 width=6.5,height=6.5,family="Arial"
)
draw();dev.off()
ragg::agg_png(
 file.path(root,"plots",paste0(stem,".png")),
 width=6.5,height=6.5,units="in",res=300
)
draw();dev.off()
print(regions)
