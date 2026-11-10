#!/usr/bin/env Rscript
.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages(library(data.table))
root <- normalizePath(dirname(sub("^--file=", "",
  commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])),
  winslash="/", mustWork=TRUE)
project <- normalizePath(file.path(root, "../.."), winslash="/", mustWork=TRUE)
panel_file <- file.path(root, "tables/restricted_ObAssoc_L1L2_locked_marker_panels.tsv")
tax_file <- file.path(project,
  "analyses/2026-09-29_BMI_enriched_SV_function_modelQC/02_circos_inputs/01_SGB_sector_labels.tsv")
function_file <- file.path(project,
  "analyses/2026-09-29_BMI_enriched_SV_function_modelQC/01_function_classification/02_CDS_function_categories_long.tsv")
stopifnot(all(file.exists(c(panel_file, tax_file, function_file))))
p <- fread(panel_file)
tax <- fread(tax_file)
fun <- fread(function_file)[BMI_enriched_group=="Ob"]
stopifnot(setequal(unique(p$Panel), c("species", "sv", "joint")),
          !anyDuplicated(p[, .(Panel, Feature)]),
          all(p[Type=="SV", SV_type] %in% c("dSGV", "vSGV")))
pick <- rbind(p[Panel=="joint"], p[Panel=="species"], p[Panel=="sv"])
out <- pick[!duplicated(Feature),
  .(Feature, Original_feature, Type, SV_type, SGB, Ob_species_class,
    Ob_L1, L2_nominal, L2_selected_BH)]
out[Type=="species", SGB := sub("^.*\\|(t__SGB[0-9]+)$", "\\1", Original_feature)]
stopifnot(all(grepl("^t__SGB[0-9]+$", out$SGB)))
out[, Taxon := tax$Species[match(SGB, tax$SGB)]]
out[Type=="species" & is.na(Taxon), Taxon := sub("^s__([^|]+).*", "\\1", Original_feature)]
out[Type=="species", Taxon := gsub("_", " ", Taxon, fixed=TRUE)]
out[, Region := fifelse(Type=="SV", sub("^.*:", "", Original_feature), NA_character_)]
for (q in c("species", "sv", "joint")) {
  col <- paste0("In_", q, "_panel")
  rank_col <- paste0(q, "_rank")
  out[, (col) := Feature %in% p[Panel==q, Feature]]
  out[, (rank_col) := p[Panel==q, Rank][match(Feature, p[Panel==q, Feature])]]
}
nonempty <- function(x) unique(x[!is.na(x) & nzchar(trimws(x)) & x!="NA"])
for (col in c("Function_category", "Preferred_name", "Description")) {
  name <- switch(col, Function_category="Annotated_function_categories",
                 Preferred_name="Annotated_gene_symbols",
                 Description="Annotated_description_examples")
  out[, (name) := vapply(Original_feature, function(z) {
    if (!z %in% fun$SV_feature) return(NA_character_)
    v <- nonempty(fun[SV_feature==z, get(col)])
    if (col=="Function_category") v <- setdiff(v, "Others")
    if (col=="Description") v <- head(v, 3L)
    if (!length(v)) NA_character_ else paste(v, collapse="; ")
  }, character(1))]
}
out[, Annotated_CDS_count := vapply(Original_feature, function(z)
  if (z %in% fun$SV_feature) uniqueN(fun[SV_feature==z, Protein_key]) else 0L,
  integer(1))]
setorder(out, -In_joint_panel, joint_rank, -In_species_panel,
         species_rank, -In_sv_panel, sv_rank, Feature)
joint <- out[In_joint_panel==TRUE]
joint_sgb <- joint[Type=="species", SGB]
sv_sgb <- unique(joint[Type=="SV", SGB])
qc <- data.table(Metric=c("all_unique_features", "species_only_members",
  "SV_only_members", "joint_members", "joint_species", "joint_SV",
  "joint_dSV", "joint_vSV", "joint_SV_host_SGBs",
  "joint_SV_host_overlap_with_abundance_SGBs",
  "joint_species_overlap_with_species_only", "joint_SV_overlap_with_SV_only"),
  Value=c(nrow(out), sum(out$In_species_panel), sum(out$In_sv_panel),
    nrow(joint), sum(joint$Type=="species"), sum(joint$Type=="SV"),
    sum(joint$SV_type=="dSGV", na.rm=TRUE),
    sum(joint$SV_type=="vSGV", na.rm=TRUE), length(sv_sgb),
    length(intersect(sv_sgb, joint_sgb)),
    sum(joint$Type=="species" & joint$In_species_panel),
    sum(joint$Type=="SV" & joint$In_sv_panel)))
stopifnot(identical(qc$Value, c(36L,12L,20L,22L,10L,12L,6L,6L,8L,0L,9L,9L)))
fwrite(out, file.path(root,
  "tables/restricted_ObAssoc_L1L2_panel_membership_supplement.tsv"),
  sep="\t", na="")
fwrite(qc, file.path(root,
  "tables/restricted_ObAssoc_L1L2_panel_membership_supplement_QC.tsv"),
  sep="\t")
print(qc)
