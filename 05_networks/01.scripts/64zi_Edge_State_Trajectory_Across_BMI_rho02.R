# ==============================================================================
# Script Name: 64zi_Edge_State_Trajectory_Across_BMI_rho02.R
# Function:
#   For each bacterial pair, summarize the stage-wise correlation state
#   trajectory (CTRL -> AA -> CRC) within each BMI stratum using |rho| > 0.2
#   as the state definition threshold.
#
# State coding:
#   "+" = rho > 0.2
#   "-" = rho < -0.2
#   "0" = -0.2 <= rho <= 0.2
# ==============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(purrr)
})

base_dir <- "G:/My Drive/proj_ob_AA_CRC"
dgca_root <- file.path(base_dir, "03.results/06_network_analysis/DGCA_Results")
out_dir <- file.path(dgca_root, "Edge_State_Trajectory_Across_BMI_rho02")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rho_cutoff <- 0.2
bmis <- c("Normal", "Ow", "Ob")
stages <- c("CTRL", "AA", "CRC")

safe_write_csv <- function(x, filename) {
  tryCatch({
    write.csv(x, filename, row.names = FALSE)
    filename
  }, error = function(e) {
    fallback_dir <- file.path(tempdir(), "Edge_State_Trajectory_Across_BMI_rho02")
    if (!dir.exists(fallback_dir)) dir.create(fallback_dir, recursive = TRUE, showWarnings = FALSE)
    fallback_file <- file.path(fallback_dir, basename(filename))
    write.csv(x, fallback_file, row.names = FALSE)
    message("Write failed for project directory, saved to temp instead: ", fallback_file)
    fallback_file
  })
}

make_edge_id <- function(g1, g2) {
  paste(pmin(g1, g2), pmax(g1, g2), sep = "___")
}

state_from_rho <- function(cor_val, rho_cutoff = 0.2) {
  case_when(
    is.na(cor_val) ~ "0",
    cor_val > rho_cutoff ~ "+",
    cor_val < -rho_cutoff ~ "-",
    TRUE ~ "0"
  )
}

get_stage_from_pair <- function(stage, bmi) {
  pair_file <- switch(
    stage,
    "CTRL" = file.path(dgca_root, paste0("CTRL_", bmi, "_vs_AA_", bmi), "table",
                       paste0("DGCA_all_CTRL_", bmi, "_vs_AA_", bmi, ".csv")),
    "AA"   = file.path(dgca_root, paste0("AA_", bmi, "_vs_CRC_", bmi), "table",
                       paste0("DGCA_all_AA_", bmi, "_vs_CRC_", bmi, ".csv")),
    "CRC"  = file.path(dgca_root, paste0("AA_", bmi, "_vs_CRC_", bmi), "table",
                       paste0("DGCA_all_AA_", bmi, "_vs_CRC_", bmi, ".csv"))
  )

  if (!file.exists(pair_file)) {
    stop("Missing DGCA file for ", stage, "_", bmi, ": ", pair_file)
  }

  df <- read_csv(pair_file, show_col_types = FALSE)

  cor_col <- paste0(stage, "_", bmi, "_cor")
  p_col <- paste0(stage, "_", bmi, "_pVal")

  df %>%
    transmute(
      BMI = bmi,
      Stage = stage,
      Gene1 = Gene1,
      Gene2 = Gene2,
      edge_id = make_edge_id(Gene1, Gene2),
      cor = .data[[cor_col]],
      p_value = .data[[p_col]],
      state = state_from_rho(.data[[cor_col]], rho_cutoff)
    )
}

stage_df <- map_dfr(bmis, function(bmi) {
  map_dfr(stages, function(stage) get_stage_from_pair(stage, bmi))
})

trajectory_df <- stage_df %>%
  select(BMI, Stage, Gene1, Gene2, edge_id, cor, p_value, state) %>%
  pivot_wider(
    id_cols = c(BMI, edge_id, Gene1, Gene2),
    names_from = Stage,
    values_from = c(cor, p_value, state),
    names_sep = "__"
  ) %>%
  mutate(
    trajectory = paste0(state__CTRL, state__AA, state__CRC),
    n_nonzero_stages = (state__CTRL != "0") + (state__AA != "0") + (state__CRC != "0")
  )

trajectory_across_bmi <- trajectory_df %>%
  select(BMI, edge_id, Gene1, Gene2, trajectory,
         state__CTRL, state__AA, state__CRC,
         cor__CTRL, cor__AA, cor__CRC,
         p_value__CTRL, p_value__AA, p_value__CRC,
         n_nonzero_stages) %>%
  pivot_wider(
    id_cols = c(edge_id, Gene1, Gene2),
    names_from = BMI,
    values_from = c(trajectory,
                    state__CTRL, state__AA, state__CRC,
                    cor__CTRL, cor__AA, cor__CRC,
                    p_value__CTRL, p_value__AA, p_value__CRC,
                    n_nonzero_stages),
    names_sep = "__"
  ) %>%
  mutate(
    trajectory_signature = paste0(
      "Nw:", trajectory__Normal,
      " | Ow:", trajectory__Ow,
      " | Ob:", trajectory__Ob
    ),
    is_bmi_trajectory_different =
      !(trajectory__Normal == trajectory__Ow &
          trajectory__Ow == trajectory__Ob),
    n_distinct_trajectories = pmap_int(
      list(trajectory__Normal, trajectory__Ow, trajectory__Ob),
      ~ dplyr::n_distinct(c(...))
    )
  ) %>%
  arrange(desc(n_distinct_trajectories), desc(is_bmi_trajectory_different))

trajectory_pattern_counts <- trajectory_across_bmi %>%
  count(trajectory__Normal, trajectory__Ow, trajectory__Ob, sort = TRUE, name = "n_edges")

interesting_edges <- trajectory_across_bmi %>%
  filter(is_bmi_trajectory_different) %>%
  arrange(desc(n_distinct_trajectories))

examples_by_pattern <- interesting_edges %>%
  group_by(trajectory__Normal, trajectory__Ow, trajectory__Ob) %>%
  slice_head(n = 10) %>%
  ungroup()

pattern_file <- safe_write_csv(trajectory_pattern_counts,
                               file.path(out_dir, "Edge_Trajectory_Pattern_Counts_Across_BMI_rho02.csv"))
wide_file <- safe_write_csv(trajectory_across_bmi,
                            file.path(out_dir, "Edge_Trajectory_Wide_Across_BMI_rho02.csv"))
example_file <- safe_write_csv(examples_by_pattern,
                               file.path(out_dir, "Edge_Trajectory_Examples_Across_BMI_rho02.csv"))

cat("\n=== Top trajectory patterns across BMI (|rho| > 0.2) ===\n")
print(trajectory_pattern_counts %>% slice_head(n = 20), n = Inf, width = Inf)

cat("\n=== Example BMI-different edges (|rho| > 0.2) ===\n")
print(
  interesting_edges %>%
    select(Gene1, Gene2, trajectory__Normal, trajectory__Ow, trajectory__Ob, trajectory_signature) %>%
    slice_head(n = 20),
  n = Inf, width = Inf
)

message("\nPattern count file: ", pattern_file)
message("Wide trajectory file: ", wide_file)
message("Example file: ", example_file)
