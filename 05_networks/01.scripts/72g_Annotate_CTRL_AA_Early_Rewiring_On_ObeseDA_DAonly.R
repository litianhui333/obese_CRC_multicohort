rm(list = ls())

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

base_dir <- "G:/My Drive/proj_ob_AA_CRC"
results_dir <- file.path(
  base_dir,
  "03.results/06_network_analysis/DGCA_Results/Union_DAonly_ObeseDA_Participating_CTRL_CRC_Rewiring_p01"
)
cyto_dir <- file.path(results_dir, "Cytoscape_Excluding_NwSameTrajectory_FullGainLoss")

edge_file <- file.path(
  results_dir,
  "Union_DAonly_ObeseDA_Participating_Rewiring_Edges_Excluding_NwSameTrajectory_p01.csv"
)
ob_ctrl_aa_file <- file.path(
  base_dir,
  "03.results/06_network_analysis/DGCA_Results/CTRL_Ob_vs_AA_Ob/table/DGCA_all_CTRL_Ob_vs_AA_Ob.csv"
)
nw_ctrl_aa_file <- file.path(
  base_dir,
  "03.results/06_network_analysis/DGCA_Results/CTRL_Normal_vs_AA_Normal/table/DGCA_all_CTRL_Normal_vs_AA_Normal.csv"
)

out_dir <- file.path(results_dir, "CTRL_AA_Early_Rewiring_Annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(edge_file), file.exists(ob_ctrl_aa_file), file.exists(nw_ctrl_aa_file))

make_edge_id <- function(g1, g2) {
  paste(pmin(g1, g2), pmax(g1, g2), sep = "___")
}

pretty_name <- function(x) {
  clean_name <- gsub(".*s__", "s__", x)
  parts <- strsplit(clean_name, "\\|")[[1]]
  sp_part <- gsub("^s__", "", parts[1])
  sgb_id <- if (length(parts) > 1) gsub("^t__", "", parts[2]) else ""

  if (grepl("_sp_", sp_part)) {
    res <- gsub("_sp_", " sp. ", sp_part)
    res <- gsub("_", " ", res)
  } else if (grepl("SGB", sp_part)) {
    res <- gsub("_", " ", sp_part)
  } else {
    sub_parts <- strsplit(sp_part, "_")[[1]]
    if (length(sub_parts) >= 2) {
      genus <- sub_parts[1]
      species <- paste(sub_parts[-1], collapse = " ")
      genus_abbr <- paste0(substr(genus, 1, 1), ".")
      res <- paste(genus_abbr, species)
    } else {
      res <- gsub("_", " ", sp_part)
    }
  }

  if (sgb_id != "") paste(res, sgb_id) else res
}

ctrl_aa_relation <- function(ob_sig, nw_sig) {
  case_when(
    ob_sig & !nw_sig ~ "Ob early rewiring only",
    ob_sig & nw_sig ~ "Shared early rewiring",
    !ob_sig & nw_sig ~ "Nw early rewiring only",
    TRUE ~ "No CTRL-AA rewiring"
  )
}

edges <- read_csv(edge_file, show_col_types = FALSE) %>%
  mutate(
    display_edge = paste(pretty_name(Gene1), pretty_name(Gene2), sep = " -- "),
    ob_stage_pattern = paste0(state__CTRL__Ob, state__AA__Ob, state__CRC__Ob),
    nw_stage_pattern = paste0(state__CTRL__Normal, state__AA__Normal, state__CRC__Normal)
  )

ob_ctrl_aa <- read_csv(ob_ctrl_aa_file, show_col_types = FALSE) %>%
  transmute(
    edge_id = make_edge_id(Gene1, Gene2),
    ob_ctrl_aa_CTRL_cor = CTRL_Ob_cor,
    ob_ctrl_aa_AA_cor = AA_Ob_cor,
    ob_ctrl_aa_z = zScoreDiff,
    ob_ctrl_aa_p = pValDiff,
    ob_ctrl_aa_p_adj = pValDiff_adj,
    ob_ctrl_aa_class = Classes
  ) %>%
  distinct(edge_id, .keep_all = TRUE)

nw_ctrl_aa <- read_csv(nw_ctrl_aa_file, show_col_types = FALSE) %>%
  transmute(
    edge_id = make_edge_id(Gene1, Gene2),
    nw_ctrl_aa_CTRL_cor = CTRL_Normal_cor,
    nw_ctrl_aa_AA_cor = AA_Normal_cor,
    nw_ctrl_aa_z = zScoreDiff,
    nw_ctrl_aa_p = pValDiff,
    nw_ctrl_aa_p_adj = pValDiff_adj,
    nw_ctrl_aa_class = Classes
  ) %>%
  distinct(edge_id, .keep_all = TRUE)

annotated_edges <- edges %>%
  left_join(ob_ctrl_aa, by = "edge_id") %>%
  left_join(nw_ctrl_aa, by = "edge_id") %>%
  mutate(
    normal_ctrl_crc_diff_sig_p01 = !is.na(pValDiff_Normal) & pValDiff_Normal < 0.01,
    nw_progression_relation = case_when(
      normal_ctrl_crc_diff_sig_p01 & trajectory__Normal != trajectory__Ob ~ "Shared progression, distinct trajectory",
      normal_ctrl_crc_diff_sig_p01 & trajectory__Normal == trajectory__Ob ~ "Shared progression, same trajectory",
      TRUE ~ "Ob-specific progression"
    ),
    ob_ctrl_aa_sig_p01 = !is.na(ob_ctrl_aa_p) & ob_ctrl_aa_p < 0.01,
    nw_ctrl_aa_sig_p01 = !is.na(nw_ctrl_aa_p) & nw_ctrl_aa_p < 0.01,
    ctrl_aa_rewiring_relation = ctrl_aa_relation(ob_ctrl_aa_sig_p01, nw_ctrl_aa_sig_p01),
    ob_aa_emergent_pattern = state__CTRL__Ob == "0" & state__AA__Ob != "0",
    nw_aa_emergent_pattern = state__CTRL__Normal == "0" & state__AA__Normal != "0",
    early_obesity_signal = ob_ctrl_aa_sig_p01 & !nw_ctrl_aa_sig_p01,
    shared_distinct_and_ob_early = nw_progression_relation == "Shared progression, distinct trajectory" &
      ob_ctrl_aa_sig_p01 & !nw_ctrl_aa_sig_p01
  ) %>%
  arrange(ctrl_aa_rewiring_relation, nw_progression_relation, pValDiff, ob_ctrl_aa_p)

write_csv(
  annotated_edges,
  file.path(out_dir, "Union_DAonly_ObeseDA_Rewiring_Edges_ExcludingNwSameTrajectory_with_CTRLAA_Annotation.csv")
)

summary_all <- annotated_edges %>%
  count(ctrl_aa_rewiring_relation, nw_progression_relation, edge_class, name = "n_edges") %>%
  arrange(ctrl_aa_rewiring_relation, nw_progression_relation, desc(n_edges))

write_csv(summary_all, file.path(out_dir, "CTRL_AA_Rewiring_Summary_AllEdges.csv"))

summary_main <- annotated_edges %>%
  count(ctrl_aa_rewiring_relation, nw_progression_relation, name = "n_edges") %>%
  arrange(desc(n_edges))

write_csv(summary_main, file.path(out_dir, "CTRL_AA_Rewiring_Summary_Main.csv"))

shared_distinct <- annotated_edges %>%
  filter(nw_progression_relation == "Shared progression, distinct trajectory")

write_csv(shared_distinct, file.path(out_dir, "SharedDistinctTrajectory_Edges_with_CTRLAA_Annotation.csv"))

shared_distinct_summary <- shared_distinct %>%
  count(ctrl_aa_rewiring_relation, edge_class, name = "n_edges") %>%
  arrange(ctrl_aa_rewiring_relation, desc(n_edges))

write_csv(
  shared_distinct_summary,
  file.path(out_dir, "SharedDistinctTrajectory_CTRL_AA_Summary.csv")
)

shared_distinct_ob_early <- shared_distinct %>%
  filter(ctrl_aa_rewiring_relation == "Ob early rewiring only") %>%
  select(
    edge_id, display_edge, Gene1, Gene2,
    edge_class, trajectory__Normal, trajectory__Ob,
    nw_progression_relation, ctrl_aa_rewiring_relation,
    ob_stage_pattern, nw_stage_pattern,
    ob_ctrl_aa_class, ob_ctrl_aa_p, nw_ctrl_aa_class, nw_ctrl_aa_p,
    source_Union = da_side, ob_da_taxon
  )

write_csv(
  shared_distinct_ob_early,
  file.path(out_dir, "SharedDistinctTrajectory_ObEarlyOnly_CTRLAA_Edges.csv")
)

top_taxa_shared_ob_early <- shared_distinct_ob_early %>%
  transmute(Gene1, Gene2) %>%
  tidyr::pivot_longer(cols = everything(), values_to = "feature") %>%
  count(feature, sort = TRUE, name = "n_edges") %>%
  mutate(display_name = vapply(feature, pretty_name, character(1)))

write_csv(
  top_taxa_shared_ob_early,
  file.path(out_dir, "SharedDistinctTrajectory_ObEarlyOnly_CTRLAA_TopTaxa.csv")
)

for (set_name in c("Full", "Gain", "Loss")) {
  file_in <- file.path(cyto_dir, paste0("Cytoscape_", set_name, "_Edges_Long.csv"))
  if (!file.exists(file_in)) next

  edges_long <- read_csv(file_in, show_col_types = FALSE) %>%
    left_join(
      annotated_edges %>%
        select(
          edge_id,
          ob_ctrl_aa_class, ob_ctrl_aa_p, ob_ctrl_aa_p_adj, ob_ctrl_aa_sig_p01,
          nw_ctrl_aa_class, nw_ctrl_aa_p, nw_ctrl_aa_p_adj, nw_ctrl_aa_sig_p01,
          ctrl_aa_rewiring_relation, ob_aa_emergent_pattern, nw_aa_emergent_pattern,
          early_obesity_signal, shared_distinct_and_ob_early
        ),
      by = "edge_id"
    )

  write_csv(
    edges_long,
    file.path(cyto_dir, paste0("Cytoscape_", set_name, "_Edges_Long_with_CTRLAA_Annotation.csv"))
  )
}

cat("Saved CTRL->AA annotation outputs to:\n", out_dir, "\n")
cat("Main summary:\n")
print(summary_main, n = Inf)
cat("\nShared progression distinct trajectory summary:\n")
print(shared_distinct_summary, n = Inf)
