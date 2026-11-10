suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

base_dir <- "G:/My Drive/proj_ob_AA_CRC"
analysis_dir <- file.path(base_dir, "analyses/2026-06-02_bmi_specific_overall_concordance")
figures_dir <- file.path(analysis_dir, "figures")
tables_dir <- file.path(analysis_dir, "tables")
data_for_figures_dir <- file.path(analysis_dir, "data_for_figures")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_for_figures_dir, recursive = TRUE, showWarnings = FALSE)

res_file <- file.path(base_dir, "03.results/04_differential_abundance/microbiome_analysis_all_results_BMI_Group1_0.7_20260123.rds")
bmi_source_dir <- file.path(base_dir, "03.results/04_differential_abundance/heatmap_final/BMI_Stratified_Disease_Heatmaps")

all_results <- readRDS(res_file)
overall <- all_results$overall$disease_analysis$analysis_results$disease_analysis_adjust_BMI$all_results$CTRL_vs_CRC$meta_fits
if (!"feature" %in% colnames(overall)) {
  overall <- tibble::rownames_to_column(overall, "feature")
}

overall <- overall %>%
  transmute(
    feature,
    overall_adjusted_coef_crc = coef,
    overall_adjusted_pval_crc = pval,
    overall_adjusted_qval_crc = qval.fdr,
    overall_adjusted_strict_sig = (!is.na(pval) & pval < 0.01) |
      (!is.na(qval.fdr) & qval.fdr < 0.1)
  )

read_specific <- function(file_name, group_label, specific_label) {
  read.csv(file.path(bmi_source_dir, file_name), stringsAsFactors = FALSE) %>%
    filter(overlap == specific_label) %>%
    transmute(
      group = group_label,
      feature,
      clean_name,
      specific_coef_crc = coef_crc,
      specific_pval_crc = pval_crc,
      specific_qval_crc = qval_crc
    )
}

detail <- bind_rows(
  read_specific("Feature_List_Normal_Refined.csv", "Nw-specific", "Nw-specific"),
  read_specific("Feature_List_Ow_Refined.csv", "Ow-specific", "Ow-specific"),
  read_specific("Feature_List_Ob_Refined.csv", "Ob-specific", "Ob-specific")
) %>%
  left_join(overall, by = "feature") %>%
  mutate(
    direction_concordant_with_overall = sign(specific_coef_crc) == sign(overall_adjusted_coef_crc),
    concordant_and_overall_significant = direction_concordant_with_overall & overall_adjusted_strict_sig,
    status = if_else(
      concordant_and_overall_significant,
      "Concordant and overall significant",
      "Other"
    )
  )

summary_tbl <- detail %>%
  count(group, status, name = "n_features") %>%
  group_by(group) %>%
  mutate(
    total_features = sum(n_features),
    proportion = n_features / total_features
  ) %>%
  ungroup()

summary_wide <- summary_tbl %>%
  filter(status == "Concordant and overall significant") %>%
  transmute(
    group,
    concordant_overall_significant_n = n_features,
    total_features,
    proportion
  )

pairwise_groups <- list(
  c("Nw-specific", "Ow-specific"),
  c("Nw-specific", "Ob-specific"),
  c("Ow-specific", "Ob-specific")
)

pairwise_chisq <- bind_rows(lapply(pairwise_groups, function(groups_i) {
  dat_i <- summary_wide %>%
    filter(group %in% groups_i) %>%
    mutate(other_n = total_features - concordant_overall_significant_n) %>%
    arrange(match(group, groups_i))

  test_mat <- as.matrix(dat_i[, c("concordant_overall_significant_n", "other_n")])
  rownames(test_mat) <- dat_i$group

  chisq_res <- suppressWarnings(chisq.test(test_mat, correct = FALSE))
  fisher_res <- fisher.test(test_mat)

  tibble(
    comparison = paste(groups_i, collapse = " vs "),
    group_1 = groups_i[1],
    group_2 = groups_i[2],
    group_1_concordant_overall_significant_n = dat_i$concordant_overall_significant_n[1],
    group_1_total_features = dat_i$total_features[1],
    group_1_proportion = dat_i$proportion[1],
    group_2_concordant_overall_significant_n = dat_i$concordant_overall_significant_n[2],
    group_2_total_features = dat_i$total_features[2],
    group_2_proportion = dat_i$proportion[2],
    chisq_statistic = unname(chisq_res$statistic),
    chisq_df = unname(chisq_res$parameter),
    chisq_p_value = chisq_res$p.value,
    min_expected_count = min(chisq_res$expected),
    fisher_p_value = fisher_res$p.value
  )
})) %>%
  mutate(chisq_p_value_bh = p.adjust(chisq_p_value, method = "BH"))

plot_df <- summary_tbl %>%
  mutate(
    group = factor(group, levels = c("Nw-specific", "Ow-specific", "Ob-specific")),
    status = factor(status, levels = c("Other", "Concordant and overall significant"))
  )

label_df <- summary_wide %>%
  mutate(
    group = factor(group, levels = c("Nw-specific", "Ow-specific", "Ob-specific")),
    label = sprintf("%d/%d\n%.1f%%", concordant_overall_significant_n, total_features, proportion * 100)
  )

write_tsv(detail, file.path(tables_dir, "bmi_specific_overall_concordance_detail.tsv"))
write_tsv(summary_wide, file.path(tables_dir, "bmi_specific_overall_concordance_summary.tsv"))
write_tsv(pairwise_chisq, file.path(tables_dir, "bmi_specific_overall_concordance_pairwise_chisq.tsv"))
write_tsv(plot_df, file.path(data_for_figures_dir, "bmi_specific_overall_concordance_stacked_bar.tsv"))

status_colors <- c(
  "Other" = "#D9D9D9",
  "Concordant and overall significant" = "#4C78A8"
)

p <- ggplot(plot_df, aes(x = group, y = proportion, fill = status)) +
  geom_col(width = 0.68, color = "black", linewidth = 0.35) +
  geom_text(
    data = label_df,
    aes(x = group, y = 0.92, label = label),
    inherit.aes = FALSE,
    size = 3.4,
    lineheight = 0.95
  ) +
  scale_fill_manual(
    values = status_colors,
    breaks = c("Concordant and overall significant", "Other"),
    labels = c("Yes", "No"),
    name = "Captured by BMI-na\u00efve analysis"
  ) +  scale_y_continuous(
    labels = function(x) paste0(round(x * 100), "%"),
    breaks = seq(0, 1, by = 0.25),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = "proportion"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(color = "black", face = "bold", size = 10),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title.y = element_text(color = "black", face = "bold", size = 10),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9)
  ) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1))

out_base <- file.path(figures_dir, "bmi_specific_overall_concordance_stacked_bar")
ggsave(paste0(out_base, ".pdf"), p, width = 4.2, height = 3.4, device = grDevices::cairo_pdf)
ggsave(paste0(out_base, ".png"), p, width = 4.2, height = 3.4, dpi = 300)

print(summary_wide)
print(pairwise_chisq)

