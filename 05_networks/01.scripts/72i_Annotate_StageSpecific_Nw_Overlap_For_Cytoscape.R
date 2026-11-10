rm(list = ls())

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
})

base_dir <- "G:/My Drive/proj_ob_AA_CRC"
results_dir <- file.path(
  base_dir,
  "03.results/06_network_analysis/DGCA_Results/Union_DAonly_ObeseDA_Participating_CTRL_CRC_Rewiring_p01"
)
cyto_dir <- file.path(results_dir, "Cytoscape_Excluding_NwSameTrajectory_FullGainLoss")
figure_table_dir <- file.path(base_dir, "Figure_final/network/tables")

dir.create(figure_table_dir, recursive = TRUE, showWarnings = FALSE)

wide_edge_file <- file.path(
  results_dir,
  "Union_DAonly_ObeseDA_Participating_Rewiring_Edges_Excluding_NwSameTrajectory_p01.csv"
)

stopifnot(file.exists(wide_edge_file), dir.exists(cyto_dir))

stage_levels <- c("CTRL", "AA", "CRC")
set_names <- c("Full", "Gain", "Loss")

collapse_true_stages <- function(values, stages = stage_levels) {
  out <- stages[values]
  if (length(out) == 0) "None" else paste(out, collapse = ";")
}

build_overlap_lookup <- function(edge_df) {
  edge_df %>%
    mutate(
      nw_state_CTRL = state__CTRL__Normal,
      nw_state_AA = state__AA__Normal,
      nw_state_CRC = state__CRC__Normal,
      ob_state_CTRL = state__CTRL__Ob,
      ob_state_AA = state__AA__Ob,
      ob_state_CRC = state__CRC__Ob,
      nw_has_signal_CTRL = nw_state_CTRL != "0",
      nw_has_signal_AA = nw_state_AA != "0",
      nw_has_signal_CRC = nw_state_CRC != "0",
      nw_stage_overlap_CTRL = ob_state_CTRL != "0" & ob_state_CTRL == nw_state_CTRL,
      nw_stage_overlap_AA = ob_state_AA != "0" & ob_state_AA == nw_state_AA,
      nw_stage_overlap_CRC = ob_state_CRC != "0" & ob_state_CRC == nw_state_CRC,
      nw_any_significant_stage = nw_has_signal_CTRL | nw_has_signal_AA | nw_has_signal_CRC,
      nw_any_stage_overlap = nw_stage_overlap_CTRL | nw_stage_overlap_AA | nw_stage_overlap_CRC,
      n_nw_significant_stages = nw_has_signal_CTRL + nw_has_signal_AA + nw_has_signal_CRC,
      n_nw_overlap_stages = nw_stage_overlap_CTRL + nw_stage_overlap_AA + nw_stage_overlap_CRC,
      nw_significant_stages = pmap_chr(
        list(nw_has_signal_CTRL, nw_has_signal_AA, nw_has_signal_CRC),
        ~ collapse_true_stages(c(..1, ..2, ..3))
      ),
      nw_overlap_stages = pmap_chr(
        list(nw_stage_overlap_CTRL, nw_stage_overlap_AA, nw_stage_overlap_CRC),
        ~ collapse_true_stages(c(..1, ..2, ..3))
      )
    ) %>%
    select(
      edge_id,
      trajectory__Normal,
      trajectory__Ob,
      nw_state_CTRL,
      nw_state_AA,
      nw_state_CRC,
      ob_state_CTRL,
      ob_state_AA,
      ob_state_CRC,
      nw_has_signal_CTRL,
      nw_has_signal_AA,
      nw_has_signal_CRC,
      nw_stage_overlap_CTRL,
      nw_stage_overlap_AA,
      nw_stage_overlap_CRC,
      nw_any_significant_stage,
      nw_any_stage_overlap,
      n_nw_significant_stages,
      n_nw_overlap_stages,
      nw_significant_stages,
      nw_overlap_stages
    )
}

annotate_stage_specific_overlap <- function(long_df, lookup_df) {
  long_df %>%
    select(-matches(
      "^(trajectory__Normal|trajectory__Ob|nw_state_|ob_state_|nw_has_signal_|nw_stage_overlap_|nw_any_|n_nw_|nw_significant_stages$|nw_overlap_stages$|ob_state_same_stage$|nw_state_same_stage$|nw_same_stage_overlap$|nw_same_stage_overlap_type$|potential_mislabel_)"
    )) %>%
    left_join(lookup_df, by = "edge_id") %>%
    mutate(
      across(
        c(state_CTRL, state_AA, state_CRC,
          nw_state_CTRL, nw_state_AA, nw_state_CRC),
        as.character
      ),
      ob_state_same_stage = case_when(
        Stage == "CTRL" ~ state_CTRL,
        Stage == "AA" ~ state_AA,
        Stage == "CRC" ~ state_CRC,
        TRUE ~ NA_character_
      ),
      nw_state_same_stage = case_when(
        Stage == "CTRL" ~ nw_state_CTRL,
        Stage == "AA" ~ nw_state_AA,
        Stage == "CRC" ~ nw_state_CRC,
        TRUE ~ NA_character_
      ),
      nw_same_stage_overlap = ob_state_same_stage != "0" &
        !is.na(ob_state_same_stage) &
        ob_state_same_stage == nw_state_same_stage,
      nw_same_stage_overlap_type = case_when(
        nw_same_stage_overlap & ob_state_same_stage == "+" ~ "Positive overlap",
        nw_same_stage_overlap & ob_state_same_stage == "-" ~ "Negative overlap",
        TRUE ~ "No stage-specific overlap"
      ),
      potential_mislabel_if_any_nw_significant_used =
        nw_any_significant_stage & !nw_same_stage_overlap,
      potential_mislabel_if_any_overlap_stage_used =
        nw_any_stage_overlap & !nw_same_stage_overlap,
      potential_mislabel_if_normal_mode_used =
        normal_edge_class != "Not significant in Nw" & !nw_same_stage_overlap
    )
}

summarize_set <- function(df, set_name) {
  tibble(
    set_name = set_name,
    total_rows = nrow(df),
    total_edges = n_distinct(df$edge_id),
    stage_specific_overlap_rows = sum(df$nw_same_stage_overlap, na.rm = TRUE),
    stage_specific_overlap_edges = n_distinct(df$edge_id[df$nw_same_stage_overlap]),
    potential_mislabel_rows_if_any_nw_significant_used =
      sum(df$potential_mislabel_if_any_nw_significant_used, na.rm = TRUE),
    potential_mislabel_edges_if_any_nw_significant_used =
      n_distinct(df$edge_id[df$potential_mislabel_if_any_nw_significant_used]),
    potential_mislabel_rows_if_any_overlap_stage_used =
      sum(df$potential_mislabel_if_any_overlap_stage_used, na.rm = TRUE),
    potential_mislabel_edges_if_any_overlap_stage_used =
      n_distinct(df$edge_id[df$potential_mislabel_if_any_overlap_stage_used]),
    potential_mislabel_rows_if_normal_mode_used =
      sum(df$potential_mislabel_if_normal_mode_used, na.rm = TRUE),
    potential_mislabel_edges_if_normal_mode_used =
      n_distinct(df$edge_id[df$potential_mislabel_if_normal_mode_used])
  )
}

wide_edges <- read_csv(wide_edge_file, show_col_types = FALSE)
overlap_lookup <- build_overlap_lookup(wide_edges)

wide_edges_annotated <- wide_edges %>%
  left_join(overlap_lookup, by = c("edge_id", "trajectory__Normal", "trajectory__Ob"))

wide_out_results <- file.path(
  results_dir,
  "Union_DAonly_ObeseDA_Participating_Rewiring_Edges_Excluding_NwSameTrajectory_p01_with_StageSpecificNwOverlap.csv"
)
wide_out_figure <- file.path(
  figure_table_dir,
  "Union_DAonly_ObeseDA_Participating_Rewiring_Edges_Excluding_NwSameTrajectory_p01_with_StageSpecificNwOverlap.csv"
)

write_csv(wide_edges_annotated, wide_out_results)
write_csv(wide_edges_annotated, wide_out_figure)

annotated_sets <- list()
set_summaries <- list()
stage_summaries <- list()
mislabel_rows <- list()

for (set_name in set_names) {
  file_name <- paste0("Cytoscape_", set_name, "_Edges_Long_with_StageTransition_DGCA.csv")
  file_in <- file.path(cyto_dir, file_name)
  if (!file.exists(file_in)) next

  long_df <- read_csv(file_in, show_col_types = FALSE)
  annotated_df <- annotate_stage_specific_overlap(long_df, overlap_lookup)

  write_csv(annotated_df, file_in)
  write_csv(annotated_df, file.path(figure_table_dir, file_name))

  annotated_sets[[set_name]] <- annotated_df
  set_summaries[[set_name]] <- summarize_set(annotated_df, set_name)
  stage_summaries[[set_name]] <- annotated_df %>%
    group_by(Stage) %>%
    summarise(
      set_name = set_name,
      total_rows = n(),
      stage_specific_overlap_rows = sum(nw_same_stage_overlap, na.rm = TRUE),
      stage_specific_overlap_edges = n_distinct(edge_id[nw_same_stage_overlap]),
      potential_mislabel_rows_if_any_nw_significant_used =
        sum(potential_mislabel_if_any_nw_significant_used, na.rm = TRUE),
      potential_mislabel_edges_if_any_nw_significant_used =
        n_distinct(edge_id[potential_mislabel_if_any_nw_significant_used]),
      .groups = "drop"
    ) %>%
    select(set_name, everything())

  mislabel_rows[[set_name]] <- annotated_df %>%
    filter(potential_mislabel_if_any_nw_significant_used) %>%
    arrange(source, target, factor(Stage, levels = stage_levels)) %>%
    select(
      edge_id,
      source,
      target,
      Stage,
      edge_class,
      trajectory,
      trajectory_Normal,
      ob_state_same_stage,
      nw_state_same_stage,
      nw_significant_stages,
      nw_overlap_stages,
      nw_same_stage_overlap,
      potential_mislabel_if_any_nw_significant_used,
      normal_edge_class,
      nw_progression_relation
    )
}

set_summary_df <- bind_rows(set_summaries)
stage_summary_df <- bind_rows(stage_summaries)

write_csv(
  set_summary_df,
  file.path(figure_table_dir, "StageSpecific_NwOverlap_Summary_BySet.csv")
)
write_csv(
  stage_summary_df,
  file.path(figure_table_dir, "StageSpecific_NwOverlap_Summary_ByStage.csv")
)

if (!is.null(annotated_sets[["Full"]])) {
  full_df <- annotated_sets[["Full"]]

  full_mislabel_rows <- mislabel_rows[["Full"]]
  full_mislabel_edges <- full_df %>%
    filter(potential_mislabel_if_any_nw_significant_used) %>%
    group_by(edge_id, source, target, edge_class, trajectory, trajectory_Normal,
             nw_significant_stages, nw_overlap_stages) %>%
    summarise(
      mislabel_stages = paste(Stage, collapse = ";"),
      n_mislabel_stages = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(n_mislabel_stages), source, target)

  write_csv(
    full_mislabel_rows,
    file.path(figure_table_dir, "Cytoscape_Full_PotentialMislabel_Rows_if_AnyNwSignificantUsed.csv")
  )
  write_csv(
    full_mislabel_edges,
    file.path(figure_table_dir, "Cytoscape_Full_PotentialMislabel_Edges_if_AnyNwSignificantUsed.csv")
  )

  full_overlap_rows <- full_df %>%
    filter(nw_same_stage_overlap) %>%
    arrange(source, target, factor(Stage, levels = stage_levels)) %>%
    select(
      edge_id,
      source,
      target,
      Stage,
      edge_class,
      trajectory,
      trajectory_Normal,
      ob_state_same_stage,
      nw_state_same_stage,
      nw_same_stage_overlap_type,
      nw_significant_stages,
      nw_overlap_stages
    )

  write_csv(
    full_overlap_rows,
    file.path(figure_table_dir, "Cytoscape_Full_StageSpecific_NwOverlap_Rows.csv")
  )
}

cat("Saved annotated Cytoscape tables to:\n")
cat("  ", cyto_dir, "\n")
cat("Saved figure-ready tables to:\n")
cat("  ", figure_table_dir, "\n\n")

cat("Summary by set:\n")
print(set_summary_df, n = Inf)

if (!is.null(stage_summary_df) && nrow(stage_summary_df) > 0) {
  cat("\nSummary by stage:\n")
  print(stage_summary_df, n = Inf)
}
