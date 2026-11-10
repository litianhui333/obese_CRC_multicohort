suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(ggpubr)
})

base_dir <- "G:/My Drive/proj_ob_AA_CRC"
analysis_dir <- file.path(base_dir, "analyses/2026-05-29_dysbiosis_score_fdr_replot")
fig_dir <- file.path(analysis_dir, "figures")
tab_dir <- file.path(analysis_dir, "tables")
data_dir <- file.path(analysis_dir, "data_for_figures")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

input_score_file <- file.path(
  base_dir,
  "03.results/04_distance_analysis/statistics_subgroup/data_used_for_statistics_multi_metrics.csv"
)

group_levels <- c(
  "CTRL_Nw", "CTRL_Ow", "CTRL_Ob",
  "AA_Nw", "AA_Ow", "AA_Ob",
  "CRC_Nw", "CRC_Ow", "CRC_Ob"
)

group_label_map <- c(
  CTRL_Nw = "Nw_CTRL", CTRL_Ow = "Ow_CTRL", CTRL_Ob = "Ob_CTRL",
  AA_Nw = "Nw_AA", AA_Ow = "Ow_AA", AA_Ob = "Ob_AA",
  CRC_Nw = "Nw_CRC", CRC_Ow = "Ow_CRC", CRC_Ob = "Ob_CRC"
)

comparison_pairs <- tribble(
  ~ref,       ~target,   ~comparison_type,   ~plot_label,
  "CTRL_Nw", "CTRL_Ow", "BMI within CTRL",  "CTRL_Ow vs CTRL_Nw",
  "CTRL_Nw", "CTRL_Ob", "BMI within CTRL",  "CTRL_Ob vs CTRL_Nw",
  "AA_Nw",   "AA_Ow",   "BMI within AA",    "AA_Ow vs AA_Nw",
  "AA_Nw",   "AA_Ob",   "BMI within AA",    "AA_Ob vs AA_Nw",
  "CRC_Nw",  "CRC_Ow",  "BMI within CRC",   "CRC_Ow vs CRC_Nw",
  "CRC_Nw",  "CRC_Ob",  "BMI within CRC",   "CRC_Ob vs CRC_Nw",
  "CTRL_Nw", "AA_Nw",   "Disease within Nw", "AA_Nw vs CTRL_Nw",
  "CTRL_Nw", "CRC_Nw",  "Disease within Nw", "CRC_Nw vs CTRL_Nw",
  "CTRL_Ow", "AA_Ow",   "Disease within Ow", "AA_Ow vs CTRL_Ow",
  "CTRL_Ow", "CRC_Ow",  "Disease within Ow", "CRC_Ow vs CTRL_Ow",
  "CTRL_Ob", "AA_Ob",   "Disease within Ob", "AA_Ob vs CTRL_Ob",
  "CTRL_Ob", "CRC_Ob",  "Disease within Ob", "CRC_Ob vs CTRL_Ob"
)

plot_comparisons <- c(
  "CTRL_Ob vs CTRL_Nw",
  "CRC_Nw vs CTRL_Nw",
  "CRC_Ow vs CTRL_Ow",
  "CRC_Ob vs CTRL_Ob"
)

format_q <- function(q) {
  ifelse(is.na(q), "NA", formatC(q, format = "g", digits = 2))
}

sig_from_q <- function(q) {
  case_when(
    is.na(q) ~ "ns",
    q < 0.001 ~ "***",
    q < 0.01 ~ "**",
    q < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

shorten_study <- function(x) {
  cleaned <- x %>%
    str_replace("^Public_study__", "") %>%
    str_replace("^This_study_cohort[0-9]+__", "") %>%
    str_replace("^ONCOBIOME_", "") %>%
    str_replace_all("__", "_") %>%
    str_replace_all("_", "-") %>%
    str_squish()

  recode(
    cleaned,
    "2017-YuJ" = "YuJ_2017",
    "2021-YuJ" = "inhouse-HK_1",
    "CHN-WF-CRC" = "SunYY_2025",
    "2023-Yunnan" = "inhouse-YN",
    "2025-HK" = "inhouse-HK_2",
    .default = cleaned
  ) %>%
    str_trunc(width = 18)
}

run_lmm <- function(df, region_name, ref_group, target_group, comparison_type, adjusted = TRUE) {
  sub_df <- df %>%
    filter(Group %in% c(ref_group, target_group)) %>%
    mutate(Group = factor(Group, levels = c(ref_group, target_group)))

  if (nrow(sub_df) < 10 ||
      n_distinct(sub_df$Group) < 2 ||
      n_distinct(sub_df$Study) < 2) {
    return(tibble())
  }

  formula_str <- if (adjusted) {
    "Distance_Score ~ Group + Age + Sex + (1 | Study)"
  } else {
    "Distance_Score ~ Group + (1 | Study)"
  }

  fit <- tryCatch(
    lmer(as.formula(formula_str), data = sub_df, REML = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(tibble())

  coef_tbl <- summary(fit)$coefficients
  coef_name <- paste0("Group", target_group)
  p_col <- match("Pr(>|t|)", colnames(coef_tbl))

  if (!coef_name %in% rownames(coef_tbl) || is.na(p_col)) return(tibble())

  tibble(
    Region = region_name,
    Method = "LMM",
    Adjustment = ifelse(adjusted, "Adjusted(Age+Sex)", "Unadjusted"),
    Comparison = paste(target_group, "vs", ref_group),
    Comparison_Type = comparison_type,
    Ref = ref_group,
    Target = target_group,
    Estimate = unname(coef_tbl[coef_name, "Estimate"]),
    Std_Error = unname(coef_tbl[coef_name, "Std. Error"]),
    P_Value = unname(coef_tbl[coef_name, p_col]),
    N_Samples = nrow(sub_df),
    N_Studies = n_distinct(sub_df$Study)
  )
}

raw_df <- read.csv(input_score_file, check.names = FALSE)

score_df <- raw_df %>%
  transmute(
    SampleID = .data[["Sample.Source.ID.processed"]],
    Study = .data[["Study.name"]],
    Continent = .data[["Continent"]],
    Disease = .data[["Disease"]],
    BMI_Group1 = .data[["BMI_Group1"]],
    Group = recode(.data[["Group_Label"]],
                   CTRL_Normal = "CTRL_Nw",
                   AA_Normal = "AA_Nw",
                   CRC_Normal = "CRC_Nw",
                   .default = .data[["Group_Label"]]),
    Age = suppressWarnings(as.numeric(.data[["Age"]])),
    Sex = as.factor(.data[["Sex"]]),
    Distance_Score = suppressWarnings(as.numeric(.data[["Score_Median"]]))
  ) %>%
  filter(
    !is.na(Distance_Score),
    !is.na(Study),
    !is.na(Continent),
    Group %in% group_levels
  ) %>%
  mutate(Group = factor(Group, levels = group_levels))

study_level <- score_df %>%
  group_by(Region = "Overall", Study, Continent, Group) %>%
  summarise(
    Distance_Score = median(Distance_Score, na.rm = TRUE),
    N_Samples = n(),
    .groups = "drop"
  ) %>%
  bind_rows(
    score_df %>%
      group_by(Region = Continent, Study, Continent, Group) %>%
      summarise(
        Distance_Score = median(Distance_Score, na.rm = TRUE),
        N_Samples = n(),
        .groups = "drop"
      )
  ) %>%
  mutate(Group = factor(Group, levels = group_levels))

write.csv(score_df, file.path(data_dir, "dysbiosis_distance_sample_level_score_median.csv"), row.names = FALSE)
write.csv(study_level, file.path(data_dir, "dysbiosis_distance_study_level_median_plot_data.csv"), row.names = FALSE)

stats_input <- score_df %>%
  select(Region = Continent, everything()) %>%
  bind_rows(score_df %>% mutate(Region = "Overall"))

stats_all <- map_dfr(seq_len(nrow(comparison_pairs)), function(i) {
  pair <- comparison_pairs[i, ]

  map_dfr(c("Overall", "Asia", "Europe"), function(region_name) {
    df_i <- if (region_name == "Overall") score_df else filter(score_df, Continent == region_name)
    bind_rows(
      run_lmm(df_i, region_name, pair$ref, pair$target, pair$comparison_type, adjusted = FALSE),
      run_lmm(df_i, region_name, pair$ref, pair$target, pair$comparison_type, adjusted = TRUE)
    )
  })
}) %>%
  group_by(Region, Adjustment) %>%
  mutate(
    FDR = p.adjust(P_Value, method = "BH"),
    Significance_FDR = sig_from_q(FDR),
    Label_FDR = format_q(FDR)
  ) %>%
  ungroup()

write.csv(stats_all, file.path(tab_dir, "dysbiosis_distance_lmm_all_comparisons_rawP_and_FDR.csv"), row.names = FALSE)

plot_stats <- stats_all %>%
  filter(
    Adjustment == "Adjusted(Age+Sex)",
    Comparison %in% plot_comparisons
  ) %>%
  mutate(
    group1 = Ref,
    group2 = Target,
    label = Label_FDR
  )

write.csv(plot_stats, file.path(tab_dir, "dysbiosis_distance_plot_annotations_FDR.csv"), row.names = FALSE)

make_annotations <- function(plot_df, stat_df) {
  max_df <- plot_df %>%
    group_by(Region) %>%
    summarise(y_base = max(Distance_Score, na.rm = TRUE), .groups = "drop")

  stat_df %>%
    left_join(max_df, by = "Region") %>%
    group_by(Region) %>%
    arrange(match(Comparison, plot_comparisons), .by_group = TRUE) %>%
    mutate(y.position = y_base + row_number() * 0.018) %>%
    ungroup()
}

theme_main <- theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.title = element_text(color = "black"),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold", color = "black"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.32, "cm"),
    legend.spacing.x = unit(0.12, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    plot.subtitle = element_blank(),
    plot.margin = margin(8, 14, 8, 8)
  )

plot_panel <- function(plot_df, stat_df, title, facet = FALSE) {
  ann <- make_annotations(plot_df, stat_df)

  p <- ggplot(plot_df, aes(x = Group, y = Distance_Score)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.55) +
    geom_jitter(aes(size = N_Samples, color = Study), width = 0.16, alpha = 0.78) +
    stat_summary(fun = median, geom = "point", shape = 18, size = 3.8, color = "#D7301F") +
    stat_pvalue_manual(
      ann,
      label = "label",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = 0.35,
      size = 3.0,
      hide.ns = FALSE
    ) +
    scale_x_discrete(labels = group_label_map) +
    scale_color_discrete(
      name = "Cohort",
      labels = shorten_study,
      guide = guide_legend(ncol = 4, override.aes = list(size = 2.0, alpha = 1))
    ) +
    scale_size_continuous(range = c(1.8, 4.2)) +
    guides(size = "none") +
    labs(
      title = title,
      x = NULL,
      y = "Dysbiosis score"
    ) +
    coord_cartesian(clip = "off") +
    theme_main

  if (facet) {
    p <- p + facet_wrap(~Region, nrow = 1)
  }

  p
}

overall_plot_data <- study_level %>% filter(Region == "Overall")
overall_plot_stats <- plot_stats %>% filter(Region == "Overall")

continent_plot_data <- study_level %>% filter(Region %in% c("Asia", "Europe"))
continent_plot_stats <- plot_stats %>% filter(Region %in% c("Asia", "Europe"))

p_overall <- plot_panel(
  overall_plot_data,
  overall_plot_stats,
  "Overall dysbiosis distance score",
  facet = FALSE
)

p_continent <- plot_panel(
  continent_plot_data,
  continent_plot_stats,
  "Regional dysbiosis distance score",
  facet = TRUE
)

ggsave(file.path(fig_dir, "MainFig_Beta_C_Overall_dysbiosis_distance_FDR.pdf"), p_overall, width = 5.3, height = 6.0)
ggsave(file.path(fig_dir, "MainFig_Beta_C_Overall_dysbiosis_distance_FDR.png"), p_overall, width = 5.3, height = 6.0, dpi = 300)

ggsave(file.path(fig_dir, "SuppFig_Beta_C_Continent_dysbiosis_distance_FDR.pdf"), p_continent, width = 9.5, height = 6.0)
ggsave(file.path(fig_dir, "SuppFig_Beta_C_Continent_dysbiosis_distance_FDR.png"), p_continent, width = 9.5, height = 6.0, dpi = 300)

manifest <- tibble(
  file = c(
    "figures/MainFig_Beta_C_Overall_dysbiosis_distance_FDR.pdf",
    "figures/MainFig_Beta_C_Overall_dysbiosis_distance_FDR.png",
    "figures/SuppFig_Beta_C_Continent_dysbiosis_distance_FDR.pdf",
    "figures/SuppFig_Beta_C_Continent_dysbiosis_distance_FDR.png",
    "tables/dysbiosis_distance_lmm_all_comparisons_rawP_and_FDR.csv",
    "tables/dysbiosis_distance_plot_annotations_FDR.csv",
    "data_for_figures/dysbiosis_distance_sample_level_score_median.csv",
    "data_for_figures/dysbiosis_distance_study_level_median_plot_data.csv"
  ),
  description = c(
    "Overall study-level dysbiosis distance plot with FDR q-value annotations.",
    "PNG copy of overall plot.",
    "Asia/Europe study-level dysbiosis distance plot with FDR q-value annotations.",
    "PNG copy of regional plot.",
    "All LMM comparisons with raw P values and BH/FDR q values.",
    "Subset of FDR annotations shown on the plots.",
    "Sample-level median distance to CTRL_Nw used for statistics.",
    "Study-level median data frame used directly for plotting."
  )
)
write.csv(manifest, file.path(tab_dir, "output_manifest.csv"), row.names = FALSE)

message("Done: dysbiosis distance FDR replot outputs written to ", analysis_dir)
