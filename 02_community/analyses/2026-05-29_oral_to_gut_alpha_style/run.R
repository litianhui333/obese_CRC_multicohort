analysis_dir <- "analyses/2026-05-29_oral_to_gut_alpha_style"
source_dir <- "analyses/2026-05-28_alpha_diversity_confirmed"

fig_dir <- file.path(analysis_dir, "figures")
tab_dir <- file.path(analysis_dir, "tables")
data_fig_dir <- file.path(analysis_dir, "data_for_figures")
log_dir <- file.path(analysis_dir, "logs")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

copy_one <- function(from, to) {
  if (!file.exists(from)) stop("Missing source file: ", from)
  if (file.exists(to)) return(invisible(TRUE))
  invisible(file.copy(from, to, overwrite = TRUE))
}

table_files <- c(
  "oral_to_gut_feature_match.tsv",
  "oral_to_gut_pairwise_pooled_LMM_all.tsv",
  "oral_to_gut_pairwise_pooled_LMM_all.csv",
  "oral_to_gut_pairwise_cohort_specific_LM_all.tsv",
  "oral_to_gut_direction_concordance_adjusted.tsv",
  "Supplementary_Table_oral_to_gut_pairwise_statistics.tsv",
  "Supplementary_Table_oral_to_gut_pairwise_statistics.csv"
)

for (x in table_files) {
  copy_one(file.path(source_dir, "tables", x), file.path(tab_dir, x))
}
copy_one(
  file.path(source_dir, "tables", "metadata_with_alpha_confirmed.tsv"),
  file.path(tab_dir, "metadata_with_oral_to_gut_alpha_style.tsv")
)

fig_files <- list.files(file.path(source_dir, "figures"), pattern = "^Oral_to_Gut.*\\.(pdf|png)$", full.names = TRUE)
for (x in fig_files) {
  copy_one(x, file.path(fig_dir, basename(x)))
}

data_files <- list.files(file.path(source_dir, "data_for_figures"), pattern = "^Oral_to_Gut.*\\.tsv$", full.names = TRUE)
for (x in data_files) {
  copy_one(x, file.path(data_fig_dir, basename(x)))
}

summary_lines <- c(
  "Oral-to-gut alpha-style result archive",
  paste("Source directory:", source_dir),
  paste("Tables copied:", length(table_files) + 1),
  paste("Figures copied:", length(fig_files)),
  paste("Figure data files copied:", length(data_files)),
  "",
  "The source workflow calculates log_Oral_Score and Oral_Richness, then applies the same pairwise LMM/LM and plotting format as the alpha-diversity analysis."
)
summary_path <- file.path(tab_dir, "archive_summary.txt")
if (!file.exists(summary_path)) {
  writeLines(summary_lines, summary_path)
}

suppressWarnings(suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
}))

theme_oral <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black", size = 10),
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold", size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30"),
      legend.position = "right",
      legend.key.size = unit(0.4, "cm"),
      legend.text = element_text(size = 9)
    )
}

colors_disease <- c("CTRL" = "#2E9FDF", "AA" = "#F39C12", "CRC" = "#E74C3C")
colors_bmi <- c("Normal" = "#5AAE61", "Ow" = "#FF8C00", "Ob" = "#B07AA1")

format_sig_p <- function(p) {
  ifelse(p < 0.01, "FDR<0.01", paste0("FDR=", sprintf("%.2f", p)))
}

make_zero_rug_y <- function(df) {
  ref <- df %>%
    filter(!(Alpha_Index == "log_Oral_Score" & Alpha_Value <= -5.999)) %>%
    group_by(Alpha_Label) %>%
    summarise(
      y_min = min(Alpha_Value, na.rm = TRUE),
      y_max = max(Alpha_Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(y_range = ifelse(y_max == y_min, abs(y_max) * 0.1 + 1, y_max - y_min),
           zero_y = y_min - y_range * 0.08)

  df %>% left_join(ref %>% select(Alpha_Label, zero_y), by = "Alpha_Label")
}

make_brackets_xbmi <- function(stats, plot_df) {
  if (nrow(stats) == 0) {
    return(stats %>% mutate(x = numeric(), xend = numeric(), y = numeric(), ytick = numeric(), xmid = numeric()))
  }
  y_ref <- plot_df %>%
    group_by(Alpha_Label) %>%
    summarise(
      y_max = max(Alpha_Value, na.rm = TRUE),
      y_range = diff(range(Alpha_Value, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(y_range = ifelse(y_range == 0, abs(y_max) * 0.1 + 1, y_range))

  disease_x <- c("CTRL" = 1, "AA" = 2, "CRC" = 3)
  bmi_offset <- c("Normal" = -0.27, "Ow" = 0, "Ob" = 0.27)

  stats %>%
    left_join(y_ref, by = "Alpha_Label") %>%
    group_by(Alpha_Label, Context) %>%
    arrange(match(Comparison, c("Normal vs Ow", "Ow vs Ob", "Normal vs Ob")), .by_group = TRUE) %>%
    mutate(
      bracket_i = row_number(),
      x = unname(disease_x[Context]) + unname(bmi_offset[Ref_Group]),
      xend = unname(disease_x[Context]) + unname(bmi_offset[Target_Group]),
      y = y_max + y_range * (0.08 + 0.12 * bracket_i),
      ytick = y - y_range * 0.025,
      xmid = (x + xend) / 2
    ) %>%
    ungroup()
}

make_brackets_xdisease <- function(stats, plot_df) {
  if (nrow(stats) == 0) {
    return(stats %>% mutate(x = numeric(), xend = numeric(), y = numeric(), ytick = numeric(), xmid = numeric()))
  }
  y_ref <- plot_df %>%
    group_by(Alpha_Label) %>%
    summarise(
      y_max = max(Alpha_Value, na.rm = TRUE),
      y_range = diff(range(Alpha_Value, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(y_range = ifelse(y_range == 0, abs(y_max) * 0.1 + 1, y_range))

  bmi_x <- c("Normal" = 1, "Ow" = 2, "Ob" = 3)
  disease_offset <- c("CTRL" = -0.27, "AA" = 0, "CRC" = 0.27)

  stats %>%
    left_join(y_ref, by = "Alpha_Label") %>%
    group_by(Alpha_Label, Context) %>%
    arrange(match(Comparison, c("CTRL vs AA", "AA vs CRC", "CTRL vs CRC")), .by_group = TRUE) %>%
    mutate(
      bracket_i = row_number(),
      x = unname(bmi_x[Context]) + unname(disease_offset[Ref_Group]),
      xend = unname(bmi_x[Context]) + unname(disease_offset[Target_Group]),
      y = y_max + y_range * (0.08 + 0.12 * bracket_i),
      ytick = y - y_range * 0.025,
      xmid = (x + xend) / 2
    ) %>%
    ungroup()
}

make_zero_rug_boxplot <- function(subgroup_val, mode, filename, width = 6.4, height = 3.2) {
  plot_df <- read_tsv(file.path(data_fig_dir, paste0(filename, ".tsv")), show_col_types = FALSE) %>%
    mutate(
      Alpha_Label = factor(Alpha_Label, levels = c("Oral Score", "Oral Richness")),
      Disease = factor(Disease, levels = c("CTRL", "AA", "CRC")),
      BMI_Display = factor(BMI_Display, levels = c("Normal", "Ow", "Ob")),
      Is_Zero_Oral_Score = Alpha_Index == "log_Oral_Score" & Alpha_Value <= -5.999
    ) %>%
    make_zero_rug_y()

  dist_df <- plot_df %>% filter(!Is_Zero_Oral_Score)
  zero_df <- plot_df %>% filter(Is_Zero_Oral_Score)
  write_tsv(plot_df, file.path(data_fig_dir, paste0(filename, "_zeroRug.tsv")))

  stats <- read_tsv(file.path(tab_dir, "oral_to_gut_pairwise_pooled_LMM_all.tsv"), show_col_types = FALSE) %>%
    filter(
      Adjustment == "Adjusted_Age_Sex",
      Subgroup_Val == subgroup_val,
      BMI_Group_Column == "BMI_Group1",
      Alpha_Index %in% c("log_Oral_Score", "Oral_Richness")
    ) %>%
    mutate(
      Alpha_Label = factor(ifelse(Alpha_Index == "log_Oral_Score", "Oral Score", "Oral Richness"),
                           levels = c("Oral Score", "Oral Richness")),
      label = ifelse(Significance == "ns", "", format_sig_p(P_FDR))
    ) %>%
    filter(label != "")

  if (mode == "xBMI") {
    stat_df <- stats %>% filter(Type == "BMI_Comparison") %>% make_brackets_xbmi(dist_df)
    p <- ggplot(dist_df, aes(x = Disease, y = Alpha_Value, fill = BMI_Display)) +
      geom_violin(position = position_dodge(0.8), alpha = 0.3, color = NA, trim = FALSE) +
      geom_boxplot(position = position_dodge(0.8), width = 0.2, outlier.size = 0.3,
                   outlier.alpha = 0.6, alpha = 0.8, color = "black", linewidth = 0.2) +
      geom_point(data = zero_df, aes(y = zero_y, color = BMI_Display),
                 position = position_dodge(0.8), shape = 124, size = 2.6,
                 alpha = 0.5, inherit.aes = TRUE, show.legend = FALSE) +
      geom_segment(data = stat_df, aes(x = x, xend = xend, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = 0.3) +
      geom_segment(data = stat_df, aes(x = x, xend = x, y = y, yend = ytick),
                   inherit.aes = FALSE, linewidth = 0.3) +
      geom_segment(data = stat_df, aes(x = xend, xend = xend, y = y, yend = ytick),
                   inherit.aes = FALSE, linewidth = 0.3) +
      geom_text(data = stat_df, aes(x = xmid, y = y, label = label),
                inherit.aes = FALSE, size = 2.25, vjust = -0.15) +
      facet_wrap(~Alpha_Label, scales = "free_y", nrow = 1) +
      scale_fill_manual(values = colors_bmi, drop = TRUE, name = "BMI") +
      scale_color_manual(values = colors_bmi, drop = TRUE, guide = "none") +
      labs(title = paste0(subgroup_val, ": oral-to-gut indices"),
           subtitle = "x-axis: Disease; bottom ticks mark zero oral-score samples",
           x = NULL, y = NULL, fill = "BMI Group") +
      theme_oral(12)
  } else {
    stat_df <- stats %>% filter(Type == "Disease_Comparison") %>% make_brackets_xdisease(dist_df)
    p <- ggplot(dist_df, aes(x = BMI_Display, y = Alpha_Value, fill = Disease)) +
      geom_violin(position = position_dodge(0.8), alpha = 0.3, color = NA, trim = FALSE) +
      geom_boxplot(position = position_dodge(0.8), width = 0.2, outlier.size = 0.3,
                   outlier.alpha = 0.6, alpha = 0.8, color = "black", linewidth = 0.2) +
      geom_point(data = zero_df, aes(y = zero_y, color = Disease),
                 position = position_dodge(0.8), shape = 124, size = 2.6,
                 alpha = 0.5, inherit.aes = TRUE, show.legend = FALSE) +
      geom_segment(data = stat_df, aes(x = x, xend = xend, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = 0.3) +
      geom_segment(data = stat_df, aes(x = x, xend = x, y = y, yend = ytick),
                   inherit.aes = FALSE, linewidth = 0.3) +
      geom_segment(data = stat_df, aes(x = xend, xend = xend, y = y, yend = ytick),
                   inherit.aes = FALSE, linewidth = 0.3) +
      geom_text(data = stat_df, aes(x = xmid, y = y, label = label),
                inherit.aes = FALSE, size = 2.25, vjust = -0.15) +
      facet_wrap(~Alpha_Label, scales = "free_y", nrow = 1) +
      scale_fill_manual(values = colors_disease, drop = TRUE, name = "Disease") +
      scale_color_manual(values = colors_disease, drop = TRUE, guide = "none") +
      labs(title = paste0(subgroup_val, ": oral-to-gut indices"),
           subtitle = "x-axis: BMI; bottom ticks mark zero oral-score samples",
           x = NULL, y = NULL, fill = "Disease") +
      theme_oral(12)
  }

  out_base <- file.path(fig_dir, paste0(filename, "_zeroRug"))
  ggsave(paste0(out_base, ".pdf"), p, width = width, height = height, device = cairo_pdf)
  ggsave(paste0(out_base, ".png"), p, width = width, height = height, dpi = 300)
}

for (region in c("Overall", "Asia", "Europe")) {
  make_zero_rug_boxplot(region, "xBMI", paste0("Oral_to_Gut_", region, "_confirmed_LMM_adjP_xBMI_BMI_Group1"))
  make_zero_rug_boxplot(region, "xDisease", paste0("Oral_to_Gut_", region, "_confirmed_LMM_adjP_xDisease_BMI_Group1"))
}

message("Done")
