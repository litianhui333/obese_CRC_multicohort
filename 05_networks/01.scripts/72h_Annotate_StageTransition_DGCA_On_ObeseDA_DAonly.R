rm(list = ls())

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

base_dir <- "G:/My Drive/proj_ob_AA_CRC"
results_dir <- file.path(
  base_dir,
  "03.results/06_network_analysis/DGCA_Results/Union_DAonly_ObeseDA_Participating_CTRL_CRC_Rewiring_p01"
)

edge_file <- file.path(
  results_dir,
  "Union_DAonly_ObeseDA_Participating_Rewiring_Edges_Excluding_NwSameTrajectory_p01.csv"
)

cyto_dir <- file.path(results_dir, "Cytoscape_Excluding_NwSameTrajectory_FullGainLoss")
out_dir <- file.path(results_dir, "StageTransition_DGCA_Annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ctrl_aa_ob_file <- file.path(base_dir, "03.results/06_network_analysis/DGCA_Results/CTRL_Ob_vs_AA_Ob/table/DGCA_all_CTRL_Ob_vs_AA_Ob.csv")
ctrl_aa_nw_file <- file.path(base_dir, "03.results/06_network_analysis/DGCA_Results/CTRL_Normal_vs_AA_Normal/table/DGCA_all_CTRL_Normal_vs_AA_Normal.csv")
aa_crc_ob_file <- file.path(base_dir, "03.results/06_network_analysis/DGCA_Results/AA_Ob_vs_CRC_Ob/table/DGCA_all_AA_Ob_vs_CRC_Ob.csv")
aa_crc_nw_file <- file.path(base_dir, "03.results/06_network_analysis/DGCA_Results/AA_Normal_vs_CRC_Normal/table/DGCA_all_AA_Normal_vs_CRC_Normal.csv")

stopifnot(
  file.exists(edge_file),
  file.exists(ctrl_aa_ob_file),
  file.exists(ctrl_aa_nw_file),
  file.exists(aa_crc_ob_file),
  file.exists(aa_crc_nw_file)
)

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

transition_relation <- function(ob_sig, nw_sig) {
  case_when(
    ob_sig & !nw_sig ~ "Ob-only significant",
    ob_sig & nw_sig ~ "Both significant",
    !ob_sig & nw_sig ~ "Nw-only significant",
    TRUE ~ "Neither significant"
  )
}

edges <- read_csv(edge_file, show_col_types = FALSE) %>%
  mutate(
    display_edge = paste(pretty_name(Gene1), pretty_name(Gene2), sep = " -- "),
    normal_ctrl_crc_diff_sig_p01 = !is.na(pValDiff_Normal) & pValDiff_Normal < 0.01,
    nw_progression_relation = case_when(
      normal_ctrl_crc_diff_sig_p01 & trajectory__Normal != trajectory__Ob ~ "Shared progression, distinct trajectory",
      normal_ctrl_crc_diff_sig_p01 & trajectory__Normal == trajectory__Ob ~ "Shared progression, same trajectory",
      TRUE ~ "Ob-specific progression"
    )
  )

read_transition <- function(path, prefix, cor1, cor2) {
  read_csv(path, show_col_types = FALSE) %>%
    transmute(
      edge_id = make_edge_id(Gene1, Gene2),
      !!paste0(prefix, "_cor_from") := .data[[cor1]],
      !!paste0(prefix, "_cor_to") := .data[[cor2]],
      !!paste0(prefix, "_z") := zScoreDiff,
      !!paste0(prefix, "_p") := pValDiff,
      !!paste0(prefix, "_p_adj") := pValDiff_adj,
      !!paste0(prefix, "_class") := Classes
    ) %>%
    distinct(edge_id, .keep_all = TRUE)
}

ctrl_aa_ob <- read_transition(ctrl_aa_ob_file, "ob_ctrl_aa", "CTRL_Ob_cor", "AA_Ob_cor")
ctrl_aa_nw <- read_transition(ctrl_aa_nw_file, "nw_ctrl_aa", "CTRL_Normal_cor", "AA_Normal_cor")
aa_crc_ob <- read_transition(aa_crc_ob_file, "ob_aa_crc", "AA_Ob_cor", "CRC_Ob_cor")
aa_crc_nw <- read_transition(aa_crc_nw_file, "nw_aa_crc", "AA_Normal_cor", "CRC_Normal_cor")

annotated <- edges %>%
  left_join(ctrl_aa_ob, by = "edge_id") %>%
  left_join(ctrl_aa_nw, by = "edge_id") %>%
  left_join(aa_crc_ob, by = "edge_id") %>%
  left_join(aa_crc_nw, by = "edge_id") %>%
  mutate(
    ob_ctrl_aa_sig_p01 = !is.na(ob_ctrl_aa_p) & ob_ctrl_aa_p < 0.01,
    nw_ctrl_aa_sig_p01 = !is.na(nw_ctrl_aa_p) & nw_ctrl_aa_p < 0.01,
    ob_aa_crc_sig_p01 = !is.na(ob_aa_crc_p) & ob_aa_crc_p < 0.01,
    nw_aa_crc_sig_p01 = !is.na(nw_aa_crc_p) & nw_aa_crc_p < 0.01,
    ctrl_aa_relation = transition_relation(ob_ctrl_aa_sig_p01, nw_ctrl_aa_sig_p01),
    aa_crc_relation = transition_relation(ob_aa_crc_sig_p01, nw_aa_crc_sig_p01),
    stage_shift_pattern = case_when(
      ob_ctrl_aa_sig_p01 & !nw_ctrl_aa_sig_p01 & !ob_aa_crc_sig_p01 & nw_aa_crc_sig_p01 ~ "Ob early, Nw late",
      ob_ctrl_aa_sig_p01 & !nw_ctrl_aa_sig_p01 & nw_aa_crc_sig_p01 ~ "Ob early, Nw later/also late",
      ob_ctrl_aa_sig_p01 & !nw_ctrl_aa_sig_p01 ~ "Ob early only",
      !ob_ctrl_aa_sig_p01 & nw_aa_crc_sig_p01 ~ "Nw late only",
      TRUE ~ "Other/none"
    )
  ) %>%
  arrange(nw_progression_relation, edge_class, ob_ctrl_aa_p, ob_aa_crc_p)

write_csv(
  annotated,
  file.path(out_dir, "Union_DAonly_ObeseDA_Rewiring_Edges_with_CTRLAA_and_AACRC_DGCA.csv")
)

summary_shared <- annotated %>%
  filter(nw_progression_relation == "Shared progression, distinct trajectory") %>%
  count(ctrl_aa_relation, aa_crc_relation, stage_shift_pattern, edge_class, name = "n_edges") %>%
  arrange(desc(n_edges))

write_csv(summary_shared, file.path(out_dir, "SharedDistinctTrajectory_StageTransition_Summary.csv"))

shared_ob_early_nw_late <- annotated %>%
  filter(
    nw_progression_relation == "Shared progression, distinct trajectory",
    ob_ctrl_aa_sig_p01,
    !nw_ctrl_aa_sig_p01,
    nw_aa_crc_sig_p01
  ) %>%
  select(
    edge_id, display_edge, Gene1, Gene2, edge_class,
    trajectory__Normal, trajectory__Ob,
    ob_ctrl_aa_class, ob_ctrl_aa_p,
    nw_ctrl_aa_class, nw_ctrl_aa_p,
    ob_aa_crc_class, ob_aa_crc_p,
    nw_aa_crc_class, nw_aa_crc_p,
    ctrl_aa_relation, aa_crc_relation, stage_shift_pattern,
    da_side, ob_da_taxon
  )

write_csv(
  shared_ob_early_nw_late,
  file.path(out_dir, "SharedDistinctTrajectory_ObEarly_NwLate_Edges.csv")
)

shared_ob_early_nw_late_taxa <- shared_ob_early_nw_late %>%
  transmute(Gene1, Gene2) %>%
  pivot_longer(cols = everything(), values_to = "feature") %>%
  count(feature, sort = TRUE, name = "n_edges") %>%
  mutate(display_name = vapply(feature, pretty_name, character(1)))

write_csv(
  shared_ob_early_nw_late_taxa,
  file.path(out_dir, "SharedDistinctTrajectory_ObEarly_NwLate_TopTaxa.csv")
)

main_counts <- annotated %>%
  summarise(
    total_edges = n(),
    shared_distinct_edges = sum(nw_progression_relation == "Shared progression, distinct trajectory"),
    shared_distinct_ob_early = sum(nw_progression_relation == "Shared progression, distinct trajectory" & ob_ctrl_aa_sig_p01 & !nw_ctrl_aa_sig_p01),
    shared_distinct_ob_early_nw_late = sum(nw_progression_relation == "Shared progression, distinct trajectory" & ob_ctrl_aa_sig_p01 & !nw_ctrl_aa_sig_p01 & nw_aa_crc_sig_p01)
  ) %>%
  pivot_longer(cols = everything(), names_to = "metric", values_to = "value")

write_csv(main_counts, file.path(out_dir, "StageTransition_MainCounts.csv"))

for (set_name in c("Full", "Gain", "Loss")) {
  file_in <- file.path(cyto_dir, paste0("Cytoscape_", set_name, "_Edges_Long_with_CTRLAA_Annotation.csv"))
  if (!file.exists(file_in)) next

  edges_long <- read_csv(file_in, show_col_types = FALSE) %>%
    left_join(
      annotated %>%
        select(
          edge_id,
          ob_aa_crc_class, ob_aa_crc_p, ob_aa_crc_p_adj, ob_aa_crc_sig_p01,
          nw_aa_crc_class, nw_aa_crc_p, nw_aa_crc_p_adj, nw_aa_crc_sig_p01,
          aa_crc_relation, stage_shift_pattern
        ),
      by = "edge_id"
    )

  write_csv(
    edges_long,
    file.path(cyto_dir, paste0("Cytoscape_", set_name, "_Edges_Long_with_StageTransition_DGCA.csv"))
  )
}

cat("Saved stage-transition DGCA annotation outputs to:\n", out_dir, "\n")
cat("\nMain counts:\n")
print(main_counts, n = Inf)
cat("\nShared distinct transition summary:\n")
print(summary_shared, n = 50)
