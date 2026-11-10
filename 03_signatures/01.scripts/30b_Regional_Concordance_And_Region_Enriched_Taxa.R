rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
})

# ==============================================================================
# Purpose:
#   Regional sensitivity analysis for BMI-stratified disease-associated taxa.
#
#   1. Test whether BMI-stratified disease-associated taxa show concordant
#      disease-effect directions across regions.
#   2. Identify region-enriched taxa using an explicit heterogeneity criterion,
#      rather than calling a taxon region-specific simply because it is
#      significant in one region and not significant in another.
#
# Main claim supported:
#   Most BMI-stratified disease-associated taxa showed concordant directions
#   across regional strata, while a smaller subset showed region-enriched effects.
# ==============================================================================

base_dir <- "G:/My Drive/proj_ob_AA_CRC"
input_dir <- file.path(base_dir, "03.results/04_differential_abundance")
heatmap_dir <- file.path(
  base_dir,
  "03.results/04_differential_abundance/heatmap_final/BMI_Stratified_Disease_Heatmaps"
)
out_dir <- file.path(
  base_dir,
  "03.results/04_differential_abundance/regional_concordance_refined"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

analysis_rds <- sort(
  list.files(input_dir, pattern = "microbiome_analysis_all_results_BMI_Group1_0.7.*\\.rds", full.names = TRUE),
  decreasing = TRUE
)[1]
if (is.na(analysis_rds) || !file.exists(analysis_rds)) {
  analysis_rds <- sort(
    list.files(heatmap_dir, pattern = "microbiome_analysis_all_results_BMI_Group1.*\\.rds", full.names = TRUE),
    decreasing = TRUE
  )[1]
}
if (is.na(analysis_rds) || !file.exists(analysis_rds)) {
  stop("Cannot find microbiome_analysis_all_results_BMI_Group1 RDS.")
}

all_analysis_results <- readRDS(analysis_rds)

target_comps <- c("CTRL_vs_CRC", "CTRL_vs_AA", "AA_vs_CRC")
bmi_groups <- c("Normal", "Ow", "Ob")
regions <- c("Asia", "Europe")

robust_p_cutoff <- 0.01
robust_q_cutoff <- 0.10
nominal_p_cutoff <- 0.05
het_p_cutoff <- 0.05
het_fdr_cutoff <- 0.10
min_abs_coef_diff <- 0.002

bmi_label <- c("Normal" = "Nw", "Ow" = "Ow", "Ob" = "Ob")

format_sgb_names <- function(names_vec) {
  sapply(names_vec, function(x) {
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
        res <- paste(paste0(substr(genus, 1, 1), "."), species)
      } else {
        res <- gsub("_", " ", sp_part)
      }
    }

    final_res <- trimws(res)
    if (sgb_id != "" && !grepl(sgb_id, final_res)) {
      final_res <- paste(final_res, sgb_id)
    }
    final_res
  })
}

read_feature_list <- function(bmi_group) {
  file_name <- switch(
    bmi_group,
    Normal = "Feature_List_Normal_Refined.csv",
    Ow = "Feature_List_Ow_Refined.csv",
    Ob = "Feature_List_Ob_Refined.csv"
  )
  read_csv(file.path(heatmap_dir, file_name), show_col_types = FALSE) %>%
    mutate(
      BMI_Group = bmi_group,
      BMI_Label = bmi_label[[bmi_group]],
      clean_name = if_else(is.na(clean_name) | clean_name == "", format_sgb_names(feature), clean_name)
    ) %>%
    select(feature, clean_name, heatmap_overlap = overlap, BMI_Group, BMI_Label)
}

extract_meta_fits <- function(region, bmi_group, comp) {
  list_name <- paste0("disease_analysis_", bmi_group)
  obj <- tryCatch(
    all_analysis_results[[region]][["disease_analysis"]][["analysis_results"]][[list_name]][["all_results"]][[comp]][["meta_fits"]],
    error = function(e) NULL
  )
  if (is.null(obj)) {
    return(tibble())
  }

  tibble(
    feature = rownames(obj),
    region = region,
    BMI_Group = bmi_group,
    Comparison = comp,
    coef = obj$coef,
    stderr = obj$stderr,
    pval = obj$pval,
    qval = if ("qval.fdr" %in% names(obj)) obj$qval.fdr else NA_real_
  )
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  cor(x[ok], y[ok])
}

safe_ggsave <- function(filename, plot, width, height, dpi = NULL) {
  tryCatch({
    if (is.null(dpi)) {
      ggsave(filename, plot, width = width, height = height)
    } else {
      ggsave(filename, plot, width = width, height = height, dpi = dpi)
    }
  }, error = function(e) {
    fallback <- sub("(\\.[^.]+)$", "_new\\1", filename)
    message("Could not overwrite ", filename, "; writing ", fallback, " instead.")
    if (is.null(dpi)) {
      ggsave(fallback, plot, width = width, height = height)
    } else {
      ggsave(fallback, plot, width = width, height = height, dpi = dpi)
    }
  })
}

feature_universe <- bind_rows(
  read_feature_list("Normal"),
  read_feature_list("Ow"),
  read_feature_list("Ob")
) %>%
  distinct(feature, BMI_Group, .keep_all = TRUE)

regional_long <- bind_rows(lapply(target_comps, function(comp) {
  bind_rows(lapply(bmi_groups, function(bmi) {
    bind_rows(lapply(regions, function(region) {
      extract_meta_fits(region, bmi, comp)
    }))
  }))
})) %>%
  semi_join(feature_universe, by = c("feature", "BMI_Group"))

regional_wide <- regional_long %>%
  select(feature, BMI_Group, Comparison, region, coef, stderr, pval, qval) %>%
  pivot_wider(
    names_from = region,
    values_from = c(coef, stderr, pval, qval)
  ) %>%
  left_join(feature_universe, by = c("feature", "BMI_Group")) %>%
  mutate(
    BMI_Label = bmi_label[BMI_Group],
    coef_diff_Asia_minus_Europe = coef_Asia - coef_Europe,
    abs_coef_diff = abs(coef_diff_Asia_minus_Europe),
    Het_Z = coef_diff_Asia_minus_Europe / sqrt(stderr_Asia^2 + stderr_Europe^2),
    Het_P = 2 * pnorm(-abs(Het_Z)),
    same_direction = sign(coef_Asia) == sign(coef_Europe),
    both_nominal = pval_Asia < nominal_p_cutoff & pval_Europe < nominal_p_cutoff,
    any_nominal = pval_Asia < nominal_p_cutoff | pval_Europe < nominal_p_cutoff,
    asia_robust = pval_Asia < robust_p_cutoff & qval_Asia < robust_q_cutoff,
    europe_robust = pval_Europe < robust_p_cutoff & qval_Europe < robust_q_cutoff,
    effect_direction = case_when(
      coef_Asia > 0 & coef_Europe > 0 ~ "Both enriched in CRC",
      coef_Asia < 0 & coef_Europe < 0 ~ "Both depleted in CRC",
      TRUE ~ "Direction discordant"
    )
  ) %>%
  group_by(Comparison, BMI_Group) %>%
  mutate(Het_FDR = p.adjust(Het_P, method = "fdr")) %>%
  ungroup() %>%
  mutate(
    regional_class = case_when(
      is.na(coef_Asia) | is.na(coef_Europe) ~ "Missing region estimate",
      !same_direction & Het_P < het_p_cutoff ~ "Direction-discordant",
      Het_FDR < het_fdr_cutoff & abs_coef_diff >= min_abs_coef_diff &
        abs(coef_Asia) > abs(coef_Europe) &
        (asia_robust | pval_Asia < nominal_p_cutoff) ~ "Asia-enriched",
      Het_FDR < het_fdr_cutoff & abs_coef_diff >= min_abs_coef_diff &
        abs(coef_Europe) > abs(coef_Asia) &
        (europe_robust | pval_Europe < nominal_p_cutoff) ~ "Europe-enriched",
      same_direction & any_nominal ~ "Region-concordant",
      TRUE ~ "Weak/uncertain"
    ),
    regional_class = factor(
      regional_class,
      levels = c(
        "Region-concordant",
        "Asia-enriched",
        "Europe-enriched",
        "Direction-discordant",
        "Weak/uncertain",
        "Missing region estimate"
      )
    ),
    BMI_Label = factor(BMI_Label, levels = c("Nw", "Ow", "Ob")),
    BMI_Group = factor(BMI_Group, levels = c("Normal", "Ow", "Ob"))
  )

concordance_summary <- regional_wide %>%
  filter(!is.na(coef_Asia), !is.na(coef_Europe)) %>%
  group_by(Comparison, BMI_Group, BMI_Label) %>%
  summarise(
    n_taxa = n(),
    n_same_direction = sum(same_direction, na.rm = TRUE),
    prop_same_direction = n_same_direction / n_taxa,
    n_any_nominal = sum(any_nominal, na.rm = TRUE),
    n_same_direction_any_nominal = sum(same_direction & any_nominal, na.rm = TRUE),
    prop_same_direction_any_nominal = n_same_direction_any_nominal / n_any_nominal,
    asia_europe_coef_cor = safe_cor(coef_Asia, coef_Europe),
    n_region_enriched = sum(regional_class %in% c("Asia-enriched", "Europe-enriched"), na.rm = TRUE),
    n_direction_discordant = sum(regional_class == "Direction-discordant", na.rm = TRUE),
    .groups = "drop"
  )

class_summary <- regional_wide %>%
  count(Comparison, BMI_Group, BMI_Label, regional_class, name = "n_taxa") %>%
  group_by(Comparison, BMI_Group) %>%
  mutate(prop = n_taxa / sum(n_taxa)) %>%
  ungroup()

region_enriched_taxa <- regional_wide %>%
  filter(regional_class %in% c("Asia-enriched", "Europe-enriched", "Direction-discordant")) %>%
  arrange(Comparison, BMI_Group, regional_class, Het_FDR, desc(abs_coef_diff)) %>%
  select(
    feature, clean_name, BMI_Group, BMI_Label, Comparison, heatmap_overlap,
    regional_class, effect_direction,
    coef_Asia, stderr_Asia, pval_Asia, qval_Asia,
    coef_Europe, stderr_Europe, pval_Europe, qval_Europe,
    coef_diff_Asia_minus_Europe, abs_coef_diff, Het_Z, Het_P, Het_FDR
  )

write_csv(regional_wide, file.path(out_dir, "Regional_Effect_Concordance_All_Taxa.csv"))
write_csv(concordance_summary, file.path(out_dir, "Regional_Concordance_Summary.csv"))
write_csv(class_summary, file.path(out_dir, "Regional_Class_Composition.csv"))
write_csv(region_enriched_taxa, file.path(out_dir, "Region_Enriched_And_Discordant_Taxa.csv"))

plot_ctrl_crc <- regional_wide %>%
  filter(Comparison == "CTRL_vs_CRC", !is.na(coef_Asia), !is.na(coef_Europe))

cor_labels <- concordance_summary %>%
  filter(Comparison == "CTRL_vs_CRC") %>%
  mutate(
    label = sprintf(
      "R = %.2f\nsame dir = %.0f%%",
      asia_europe_coef_cor,
      100 * prop_same_direction
    ),
    x = -Inf,
    y = Inf
  )

p_scatter <- ggplot(plot_ctrl_crc, aes(x = coef_Asia, y = coef_Europe)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey55") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey55") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey65") +
  geom_point(
    data = filter(plot_ctrl_crc, regional_class %in% c("Weak/uncertain", "Missing region estimate")),
    color = "grey70",
    alpha = 0.40,
    size = 1.2
  ) +
  geom_point(
    data = filter(plot_ctrl_crc, !(regional_class %in% c("Weak/uncertain", "Missing region estimate"))),
    aes(color = regional_class),
    alpha = 0.88,
    size = 2.0
  ) +
  geom_text(
    data = cor_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = -0.08,
    vjust = 1.15,
    fontface = "bold",
    size = 3.5
  ) +
  facet_wrap(~BMI_Label, nrow = 1) +
  scale_color_manual(
    values = c(
      "Region-concordant" = "#2CA02C",
      "Asia-enriched" = "#D95F02",
      "Europe-enriched" = "#1B9E77",
      "Direction-discordant" = "#7570B3"
    ),
    drop = FALSE
  ) +
  labs(
    title = "Regional concordance of BMI-stratified CRC-associated taxa",
    x = "Asia effect size",
    y = "Europe effect size",
    color = "Regional class"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

safe_ggsave(file.path(out_dir, "Regional_Concordance_CTRL_vs_CRC_Scatter.pdf"), p_scatter, width = 10.2, height = 4.1)
safe_ggsave(file.path(out_dir, "Regional_Concordance_CTRL_vs_CRC_Scatter.png"), p_scatter, width = 10.2, height = 4.1, dpi = 300)

p_stack <- class_summary %>%
  filter(Comparison == "CTRL_vs_CRC") %>%
  mutate(
    BMI_Label = factor(BMI_Label, levels = c("Nw", "Ow", "Ob")),
    BMI_Group = factor(BMI_Group, levels = c("Normal", "Ow", "Ob"))
  ) %>%
  ggplot(aes(x = BMI_Label, y = prop, fill = regional_class)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(
    values = c(
      "Region-concordant" = "#2CA02C",
      "Asia-enriched" = "#D95F02",
      "Europe-enriched" = "#1B9E77",
      "Direction-discordant" = "#7570B3",
      "Weak/uncertain" = "grey78",
      "Missing region estimate" = "grey45"
    ),
    drop = FALSE
  ) +
  labs(
    title = "Regional evidence classes among CRC-associated taxa",
    x = NULL,
    y = "Proportion of taxa",
    fill = "Regional class"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

safe_ggsave(file.path(out_dir, "Regional_Class_Composition_CTRL_vs_CRC_StackedBar.pdf"), p_stack, width = 7.2, height = 4.5)
safe_ggsave(file.path(out_dir, "Regional_Class_Composition_CTRL_vs_CRC_StackedBar.png"), p_stack, width = 7.2, height = 4.5, dpi = 300)

top_region_enriched <- region_enriched_taxa %>%
  filter(Comparison == "CTRL_vs_CRC", BMI_Group == "Ob") %>%
  group_by(regional_class) %>%
  arrange(Het_FDR, desc(abs_coef_diff)) %>%
  slice_head(n = 12) %>%
  ungroup() %>%
  mutate(
    clean_name = factor(clean_name, levels = rev(unique(clean_name))),
    regional_class = droplevels(regional_class)
  )

if (nrow(top_region_enriched) > 0) {
  p_dumbbell <- ggplot(top_region_enriched) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey55") +
    geom_segment(
      aes(y = clean_name, yend = clean_name, x = coef_Asia, xend = coef_Europe),
      color = "grey78",
      linewidth = 1.0
    ) +
    geom_errorbarh(
      aes(y = clean_name, xmin = coef_Asia - 1.96 * stderr_Asia, xmax = coef_Asia + 1.96 * stderr_Asia),
      height = 0.16,
      color = "#D95F02",
      alpha = 0.45
    ) +
    geom_errorbarh(
      aes(y = clean_name, xmin = coef_Europe - 1.96 * stderr_Europe, xmax = coef_Europe + 1.96 * stderr_Europe),
      height = 0.16,
      color = "#1B9E77",
      alpha = 0.45
    ) +
    geom_point(aes(y = clean_name, x = coef_Asia, fill = "Asia"), shape = 21, size = 3.0) +
    geom_point(aes(y = clean_name, x = coef_Europe, fill = "Europe"), shape = 21, size = 3.0) +
    facet_grid(regional_class ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = c("Asia" = "#D95F02", "Europe" = "#1B9E77")) +
    labs(
      title = "Region-enriched Obese CRC-associated taxa",
      x = "Effect size with 95% CI",
      y = NULL,
      fill = "Region"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top",
      axis.text.y = element_text(face = "italic")
    )

  safe_ggsave(file.path(out_dir, "Obese_CTRL_vs_CRC_Region_Enriched_Taxa_Dumbbell.pdf"), p_dumbbell, width = 7.2, height = 5.4)
  safe_ggsave(file.path(out_dir, "Obese_CTRL_vs_CRC_Region_Enriched_Taxa_Dumbbell.png"), p_dumbbell, width = 7.2, height = 5.4, dpi = 300)
}

cat("Saved regional concordance analysis to:\n")
cat(out_dir, "\n\n")
cat("CTRL_vs_CRC concordance summary:\n")
print(concordance_summary %>% filter(Comparison == "CTRL_vs_CRC"), n = 20)
cat("\nRegion-enriched/discordant CTRL_vs_CRC taxa counts:\n")
print(region_enriched_taxa %>% filter(Comparison == "CTRL_vs_CRC") %>% count(BMI_Label, regional_class), n = 20)
