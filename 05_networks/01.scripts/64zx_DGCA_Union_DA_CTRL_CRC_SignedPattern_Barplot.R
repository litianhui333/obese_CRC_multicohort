rm(list = ls())

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

# ==============================================================================
# Script Name: 64zx_DGCA_Union_DA_CTRL_CRC_SignedPattern_Barplot.R
# Function:
#   Summarize DGCA CTRL-to-CRC signed edge-state transitions among union DA taxa
#   and draw stacked barplots for pattern proportions.
#
# Pattern examples:
#   +/0: positive edge in CTRL, absent in CRC
#   -/0: negative edge in CTRL, absent in CRC
#   0/+: absent in CTRL, positive edge in CRC
#   0/-: absent in CTRL, negative edge in CRC
# ==============================================================================

base_dir <- "G:/My Drive/proj_ob_AA_CRC"
dgca_root <- file.path(base_dir, "03.results/06_network_analysis/DGCA_Results")
heatmap_dir <- file.path(
  base_dir,
  "03.results/04_differential_abundance/heatmap_final/BMI_Stratified_Disease_Heatmaps"
)
trajectory_file <- file.path(
  dgca_root,
  "Edge_State_Trajectory_Across_BMI_rho02_p05/Edge_Trajectory_Wide_Across_BMI_rho02_p05.csv"
)
out_dir <- file.path(dgca_root, "Union_DA_CTRL_CRC_SignedPattern_Barplot")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

da_files <- c(
  Normal = "Feature_List_Normal_Refined.csv",
  Ow = "Feature_List_Ow_Refined.csv",
  Ob = "Feature_List_Ob_Refined.csv"
)

stopifnot(file.exists(trajectory_file))
stopifnot(all(file.exists(file.path(heatmap_dir, da_files))))

bmi_levels <- c("Normal", "Ow", "Ob")
bmi_labels <- c("Normal" = "Nw", "Ow" = "Ow", "Ob" = "Ob")

pattern_levels <- c(
  "+/+", "-/-", "+/-", "-/+",
  "+/0", "-/0",
  "0/+", "0/-",
  "0/0"
)

pattern_labels <- c(
  "+/+" = "Positive retained",
  "-/-" = "Negative retained",
  "+/-" = "Positive to negative",
  "-/+" = "Negative to positive",
  "+/0" = "Positive lost",
  "-/0" = "Negative lost",
  "0/+" = "Positive gained",
  "0/-" = "Negative gained",
  "0/0" = "Absent in both"
)

pattern_colors <- c(
  "+/+" = "#E76F51",
  "-/-" = "#3A86B8",
  "+/-" = "#8A5FBF",
  "-/+" = "#2A9D8F",
  "+/0" = "#F4A261",
  "-/0" = "#8ECAE6",
  "0/+" = "#F6C85F",
  "0/-" = "#90BE6D",
  "0/0" = "#EAEAEA"
)

safe_ggsave <- function(filename, plot, ...) {
  tryCatch(
    {
      ggsave(filename, plot, ...)
      filename
    },
    error = function(e) {
      message("Could not write figure: ", filename)
      message("  ", conditionMessage(e))
      NA_character_
    }
  )
}

da_long <- bind_rows(lapply(names(da_files), function(bmi_i) {
  read_csv(file.path(heatmap_dir, da_files[[bmi_i]]), show_col_types = FALSE) %>%
    mutate(source_bmi = bmi_i)
}))
union_da <- unique(da_long$feature)

traj_df <- read_csv(trajectory_file, show_col_types = FALSE)

union_edges <- traj_df %>%
  filter(Gene1 %in% union_da, Gene2 %in% union_da)

message("Union DA taxa: ", length(union_da))
message("Union DA-DA DGCA edge universe: ", nrow(union_edges))

state_long <- bind_rows(lapply(bmi_levels, function(bmi_i) {
  ctrl_col <- paste0("state__CTRL__", bmi_i)
  crc_col <- paste0("state__CRC__", bmi_i)
  cor_ctrl_col <- paste0("cor_round__CTRL__", bmi_i)
  aa_col <- paste0("state__AA__", bmi_i)
  cor_aa_col <- paste0("cor_round__AA__", bmi_i)
  cor_crc_col <- paste0("cor_round__CRC__", bmi_i)

  union_edges %>%
    transmute(
      BMI = bmi_i,
      edge_id,
      Gene1,
      Gene2,
      state_CTRL = .data[[ctrl_col]],
      state_AA = .data[[aa_col]],
      state_CRC = .data[[crc_col]],
      cor_CTRL = .data[[cor_ctrl_col]],
      cor_AA = .data[[cor_aa_col]],
      cor_CRC = .data[[cor_crc_col]],
      pattern = paste0(state_CTRL, "/", state_CRC),
      trajectory_3stage = paste0(state_CTRL, state_AA, state_CRC),
      abs_delta_rho = abs(cor_CRC) - abs(cor_CTRL),
      signed_delta_rho = cor_CRC - cor_CTRL
    )
})) %>%
  mutate(
    BMI = factor(BMI, levels = bmi_levels),
    pattern = factor(pattern, levels = pattern_levels),
    pattern_label = factor(pattern_labels[as.character(pattern)], levels = pattern_labels[pattern_levels]),
    transition_class = case_when(
      pattern %in% c("+/0", "-/0") ~ "Lost by CRC",
      pattern %in% c("0/+", "0/-") ~ "Gained by CRC",
      pattern %in% c("+/+", "-/-") ~ "Retained same sign",
      pattern %in% c("+/-", "-/+") ~ "Sign switched",
      pattern == "0/0" ~ "Absent in both",
      TRUE ~ "Other"
    )
  )

aa_role_summary_detail <- state_long %>%
  filter(pattern != "0/0") %>%
  mutate(
    aa_opposite_to_ctrl_or_crc = case_when(
      state_AA == "0" ~ FALSE,
      state_CTRL != "0" & state_AA != state_CTRL ~ TRUE,
      state_CRC != "0" & state_AA != state_CRC ~ TRUE,
      TRUE ~ FALSE
    ),
    aa_role = case_when(
      aa_opposite_to_ctrl_or_crc ~ "AA opposite sign",
      state_CTRL != "0" & state_CRC != "0" & state_CTRL == state_CRC & state_AA == state_CTRL ~ "AA same as CTRL/CRC",
      state_CTRL != "0" & state_CRC != "0" & state_CTRL == state_CRC & state_AA == "0" ~ "AA transiently absent",
      state_CTRL != "0" & state_CRC == "0" & state_AA == state_CTRL ~ "AA retains CTRL edge",
      state_CTRL != "0" & state_CRC == "0" & state_AA == "0" ~ "AA already absent",
      state_CTRL == "0" & state_CRC != "0" & state_AA == state_CRC ~ "AA already gained",
      state_CTRL == "0" & state_CRC != "0" & state_AA == "0" ~ "AA not yet gained",
      TRUE ~ "Other"
    ),
    include_main_aa_plot = !aa_opposite_to_ctrl_or_crc
  )

aa_role_summary <- aa_role_summary_detail %>%
  filter(include_main_aa_plot) %>%
  count(BMI, pattern, trajectory_3stage, aa_role, transition_class, name = "n_edges") %>%
  group_by(BMI) %>%
  mutate(percent_nonzero_edges = 100 * n_edges / sum(n_edges)) %>%
  ungroup() %>%
  arrange(BMI, pattern, desc(n_edges))

aa_role_by_transition <- aa_role_summary_detail %>%
  filter(include_main_aa_plot) %>%
  count(BMI, transition_class, aa_role, name = "n_edges") %>%
  group_by(BMI, transition_class) %>%
  mutate(percent_within_transition = 100 * n_edges / sum(n_edges)) %>%
  ungroup() %>%
  arrange(BMI, transition_class, desc(n_edges))

aa_opposite_excluded_summary <- aa_role_summary_detail %>%
  filter(aa_opposite_to_ctrl_or_crc) %>%
  count(BMI, pattern, trajectory_3stage, name = "n_edges") %>%
  arrange(BMI, pattern, desc(n_edges))

loss_gain_timing_detail <- aa_role_summary_detail %>%
  filter(
    !aa_opposite_to_ctrl_or_crc,
    transition_class %in% c("Lost by CRC", "Gained by CRC")
  ) %>%
  mutate(
    timing_class = case_when(
      transition_class == "Lost by CRC" & state_AA == "0" ~ "Early loss at AA",
      transition_class == "Lost by CRC" & state_AA == state_CTRL ~ "Late loss at CRC",
      transition_class == "Gained by CRC" & state_AA == state_CRC ~ "Early gain at AA",
      transition_class == "Gained by CRC" & state_AA == "0" ~ "Late gain at CRC",
      TRUE ~ "Other"
    ),
    timing_label = case_when(
      timing_class == "Early loss at AA" ~ "Early loss at AA (+00 / -00)",
      timing_class == "Late loss at CRC" ~ "Late loss at CRC (++0 / --0)",
      timing_class == "Early gain at AA" ~ "Early gain at AA (0++ / 0--)",
      timing_class == "Late gain at CRC" ~ "Late gain at CRC (00+ / 00-)",
      TRUE ~ "Other"
    ),
    timing_simple = case_when(
      timing_class %in% c("Early loss at AA", "Early gain at AA") ~ "AA already changed",
      timing_class %in% c("Late loss at CRC", "Late gain at CRC") ~ "Changed at CRC",
      TRUE ~ "Other"
    )
  )

loss_gain_timing_summary <- loss_gain_timing_detail %>%
  count(BMI, transition_class, timing_class, timing_label, timing_simple, name = "n_edges") %>%
  group_by(BMI, transition_class) %>%
  mutate(
    total_transition_edges = sum(n_edges),
    percent_within_transition = 100 * n_edges / total_transition_edges
  ) %>%
  ungroup() %>%
  arrange(BMI, transition_class, desc(n_edges))

compare_nw_ob_timing <- function(df, transition_i) {
  test_df <- df %>%
    filter(BMI %in% c("Normal", "Ob"), transition_class == transition_i) %>%
    select(BMI, timing_simple, n_edges) %>%
    mutate(
      BMI = as.character(BMI),
      timing_simple = factor(timing_simple, levels = c("AA already changed", "Changed at CRC"))
    )

  tab <- test_df %>%
    tidyr::complete(BMI = c("Normal", "Ob"), timing_simple, fill = list(n_edges = 0)) %>%
    tidyr::pivot_wider(names_from = timing_simple, values_from = n_edges) %>%
    arrange(match(BMI, c("Normal", "Ob")))

  mat <- as.matrix(tab[, c("AA already changed", "Changed at CRC")])
  rownames(mat) <- tab$BMI

  fisher_res <- fisher.test(mat)
  chisq_res <- suppressWarnings(chisq.test(mat, correct = FALSE))

  tibble(
    transition_class = transition_i,
    Nw_AA_already_changed = mat["Normal", "AA already changed"],
    Nw_changed_at_CRC = mat["Normal", "Changed at CRC"],
    Ob_AA_already_changed = mat["Ob", "AA already changed"],
    Ob_changed_at_CRC = mat["Ob", "Changed at CRC"],
    Nw_AA_already_changed_percent = 100 * mat["Normal", "AA already changed"] / sum(mat["Normal", ]),
    Ob_AA_already_changed_percent = 100 * mat["Ob", "AA already changed"] / sum(mat["Ob", ]),
    fisher_p = fisher_res$p.value,
    chisq_p = chisq_res$p.value,
    odds_ratio = unname(fisher_res$estimate)
  )
}

loss_gain_timing_nw_ob_tests <- bind_rows(lapply(
  c("Lost by CRC", "Gained by CRC"),
  function(x) compare_nw_ob_timing(loss_gain_timing_summary, x)
)) %>%
  mutate(
    fisher_fdr = p.adjust(fisher_p, method = "fdr"),
    chisq_fdr = p.adjust(chisq_p, method = "fdr")
  )

pattern_summary_all <- state_long %>%
  count(BMI, pattern, pattern_label, transition_class, name = "n_edges") %>%
  group_by(BMI) %>%
  mutate(percent_all_edges = 100 * n_edges / sum(n_edges)) %>%
  ungroup() %>%
  arrange(BMI, pattern)

pattern_summary_active <- state_long %>%
  filter(pattern != "0/0") %>%
  count(BMI, pattern, pattern_label, transition_class, name = "n_edges") %>%
  group_by(BMI) %>%
  mutate(percent_nonzero_transition_edges = 100 * n_edges / sum(n_edges)) %>%
  ungroup() %>%
  arrange(BMI, pattern)

class_summary_all <- state_long %>%
  count(BMI, transition_class, name = "n_edges") %>%
  group_by(BMI) %>%
  mutate(percent_all_edges = 100 * n_edges / sum(n_edges)) %>%
  ungroup() %>%
  arrange(BMI, desc(n_edges))

class_summary_active <- state_long %>%
  filter(pattern != "0/0") %>%
  count(BMI, transition_class, name = "n_edges") %>%
  group_by(BMI) %>%
  mutate(percent_nonzero_transition_edges = 100 * n_edges / sum(n_edges)) %>%
  ungroup() %>%
  arrange(BMI, desc(n_edges))

write_csv(state_long, file.path(out_dir, "DGCA_Union_DA_CTRL_CRC_SignedPattern_EdgeDetail.csv"))
write_csv(pattern_summary_all, file.path(out_dir, "DGCA_Union_DA_CTRL_CRC_SignedPattern_Summary_AllEdges.csv"))
write_csv(pattern_summary_active, file.path(out_dir, "DGCA_Union_DA_CTRL_CRC_SignedPattern_Summary_Exclude_00.csv"))
write_csv(class_summary_all, file.path(out_dir, "DGCA_Union_DA_CTRL_CRC_TransitionClass_Summary_AllEdges.csv"))
write_csv(class_summary_active, file.path(out_dir, "DGCA_Union_DA_CTRL_CRC_TransitionClass_Summary_Exclude_00.csv"))
write_csv(aa_role_summary_detail, file.path(out_dir, "DGCA_Union_DA_CTRL_CRC_AA_Trajectory_EdgeDetail.csv"))
write_csv(aa_role_summary, file.path(out_dir, "DGCA_Union_DA_CTRL_CRC_AA_Trajectory_Summary_NoOpposite.csv"))
write_csv(aa_role_by_transition, file.path(out_dir, "DGCA_Union_DA_CTRL_CRC_AA_Role_ByTransition_NoOpposite.csv"))
write_csv(aa_opposite_excluded_summary, file.path(out_dir, "DGCA_Union_DA_CTRL_CRC_AA_Opposite_Excluded_Summary.csv"))
write_csv(loss_gain_timing_detail, file.path(out_dir, "DGCA_Union_DA_CTRL_CRC_LossGain_Timing_EdgeDetail_NoOpposite.csv"))
write_csv(loss_gain_timing_summary, file.path(out_dir, "DGCA_Union_DA_CTRL_CRC_LossGain_Timing_Summary_NoOpposite.csv"))
write_csv(loss_gain_timing_nw_ob_tests, file.path(out_dir, "DGCA_Union_DA_CTRL_CRC_LossGain_Timing_Nw_vs_Ob_Tests.csv"))

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
      legend.title = element_blank(),
      legend.position = "right"
    )
}

plot_stacked <- function(plot_df, percent_col, title, subtitle) {
  ggplot(plot_df, aes(x = BMI, y = .data[[percent_col]], fill = pattern)) +
    geom_col(width = 0.68, color = "white", linewidth = 0.25) +
    scale_x_discrete(labels = bmi_labels) +
    scale_fill_manual(
      values = pattern_colors,
      breaks = pattern_levels,
      labels = pattern_levels,
      drop = FALSE
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = "Edges (%)"
    ) +
    theme_pub()
}

p_all <- plot_stacked(
  pattern_summary_all,
  "percent_all_edges",
  "DGCA Union-DA CTRL-to-CRC Signed Edge Patterns",
  "Denominator includes all Union DA-DA edges, including 0/0"
)

p_active <- plot_stacked(
  pattern_summary_active,
  "percent_nonzero_transition_edges",
  "DGCA Union-DA CTRL-to-CRC Signed Edge Patterns",
  "0/0 excluded to emphasize retained, gained, lost, and switched edges"
)

trajectory_colors <- c(
  "+++" = "#E76F51",
  "---" = "#3A86B8",
  "++0" = "#F4A261",
  "--0" = "#8ECAE6",
  "+00" = "#F6C85F",
  "-00" = "#90BE6D",
  "0++" = "#E9C46A",
  "0--" = "#4D908E",
  "00+" = "#F9844A",
  "00-" = "#577590",
  "+0+" = "#B56576",
  "-0-" = "#43AA8B"
)

aa_role_colors <- c(
  "AA same as CTRL/CRC" = "#E76F51",
  "AA transiently absent" = "#B56576",
  "AA retains CTRL edge" = "#F4A261",
  "AA already absent" = "#F6C85F",
  "AA already gained" = "#2A9D8F",
  "AA not yet gained" = "#90BE6D",
  "Other" = "#BDBDBD"
)

loss_gain_timing_colors <- c(
  "AA already changed" = "#4E79A7",
  "Changed at CRC" = "#F28E2B",
  "Other" = "#BDBDBD"
)

p_aa_trajectory <- ggplot(
  aa_role_summary,
  aes(x = BMI, y = percent_nonzero_edges, fill = trajectory_3stage)
) +
  geom_col(width = 0.68, color = "white", linewidth = 0.25) +
  scale_x_discrete(labels = bmi_labels) +
  scale_fill_manual(values = trajectory_colors, drop = FALSE) +
  labs(
    title = "DGCA Union-DA CTRL-AA-CRC Trajectories",
    subtitle = "0/0 and AA-opposite trajectories excluded",
    x = NULL,
    y = "Edges (%)"
  ) +
  theme_pub()

p_aa_role <- ggplot(
  aa_role_by_transition,
  aes(x = BMI, y = percent_within_transition, fill = aa_role)
) +
  geom_col(width = 0.68, color = "white", linewidth = 0.25) +
  facet_wrap(~transition_class, nrow = 1) +
  scale_x_discrete(labels = bmi_labels) +
  scale_fill_manual(values = aa_role_colors, drop = FALSE) +
  labs(
    title = "AA Participation in DGCA Union-DA CTRL-to-CRC Networks",
    subtitle = "Includes retained, lost, and gained CTRL-to-CRC edges; AA-opposite trajectories excluded",
    x = NULL,
    y = "Edges within transition (%)"
  ) +
  theme_pub(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

p_loss_gain_timing <- ggplot(
  loss_gain_timing_summary,
  aes(x = BMI, y = percent_within_transition, fill = timing_simple)
) +
  geom_col(width = 0.68, color = "white", linewidth = 0.25) +
  geom_text(
    aes(label = if_else(percent_within_transition >= 8, sprintf("%.0f%%", percent_within_transition), "")),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "black"
  ) +
  facet_wrap(~transition_class, nrow = 1) +
  scale_x_discrete(labels = bmi_labels) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = loss_gain_timing_colors, drop = FALSE) +
  labs(
    title = "Timing of Lost and Gained DGCA Edges",
    subtitle = "Blue: already changed at AA; orange: changed only at CRC",
    x = NULL,
    y = "Edges within transition (%)"
  ) +
  theme_pub(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

safe_ggsave(
  file.path(out_dir, "DGCA_Union_DA_CTRL_CRC_SignedPattern_StackedBar_AllEdges.pdf"),
  p_all,
  width = 3,
  height = 3
)
safe_ggsave(
  file.path(out_dir, "DGCA_Union_DA_CTRL_CRC_SignedPattern_StackedBar_AllEdges.png"),
  p_all,
  width = 3,
  height = 3,
  dpi = 300
)
safe_ggsave(
  file.path(out_dir, "DGCA_Union_DA_CTRL_CRC_SignedPattern_StackedBar_Exclude_00.pdf"),
  p_active,
  width = 3,
  height = 3
)
safe_ggsave(
  file.path(out_dir, "DGCA_Union_DA_CTRL_CRC_SignedPattern_StackedBar_Exclude_00.png"),
  p_active,
  width = 3,
  height = 3,
  dpi = 300
)
safe_ggsave(
  file.path(out_dir, "DGCA_Union_DA_CTRL_AA_CRC_Trajectory_StackedBar_NoOpposite.pdf"),
  p_aa_trajectory,
  width = 3.4,
  height = 3.2
)
safe_ggsave(
  file.path(out_dir, "DGCA_Union_DA_CTRL_AA_CRC_Trajectory_StackedBar_NoOpposite.png"),
  p_aa_trajectory,
  width = 3.4,
  height = 3.2,
  dpi = 300
)
safe_ggsave(
  file.path(out_dir, "DGCA_Union_DA_AA_Role_By_CTRL_CRC_Transition_NoOpposite.pdf"),
  p_aa_role,
  width = 7.2,
  height = 3.2
)
safe_ggsave(
  file.path(out_dir, "DGCA_Union_DA_AA_Role_By_CTRL_CRC_Transition_NoOpposite.png"),
  p_aa_role,
  width = 7.2,
  height = 3.2,
  dpi = 300
)
safe_ggsave(
  file.path(out_dir, "DGCA_Union_DA_LossGain_Timing_StackedBar_NoOpposite.pdf"),
  p_loss_gain_timing,
  width = 6.2,
  height = 3.2
)
safe_ggsave(
  file.path(out_dir, "DGCA_Union_DA_LossGain_Timing_StackedBar_NoOpposite.png"),
  p_loss_gain_timing,
  width = 6.2,
  height = 3.2,
  dpi = 300
)

message("\nSaved outputs to: ", out_dir)
message("\nPattern summary excluding 0/0:")
print(pattern_summary_active)
message("\nTransition class summary excluding 0/0:")
print(class_summary_active)
message("\nAA trajectory summary excluding AA-opposite trajectories:")
print(aa_role_summary)
message("\nAA-opposite trajectories excluded:")
print(aa_opposite_excluded_summary)
message("\nLoss/gain timing summary:")
print(loss_gain_timing_summary)
message("\nNw vs Ob timing tests:")
print(loss_gain_timing_nw_ob_tests)
