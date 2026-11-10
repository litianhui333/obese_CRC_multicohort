#!/usr/bin/env Rscript
# MetaPhlAn strain-level (t__) paired longitudinal analysis: LN versus OB.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(vegan); library(pheatmap); library(patchwork)
})

root <- normalizePath(file.path(getwd(), "analyses", "2026-08-06_GF_mouse_metaphlan_t_longitudinal"), mustWork = TRUE)
input_file <- normalizePath(file.path(getwd(), "00.rawdata", "08_GF_mouse", "metaphlan_merged.txt"), mustWork = TRUE)
results_dir <- file.path(root, "results"); figures_dir <- file.path(root, "figures")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

ln_col <- "#3C78A8"; ob_col <- "#CC4125"; wk3_col <- "#6B7280"; wk24_col <- "#7A3E9D"
theme_set(theme_classic(base_size = 11, base_family = "Arial") + theme(panel.grid = element_blank(), plot.title = element_text(face = "bold")))
save_plot <- function(p, name, w = 8, h = 5) {
  ggsave(file.path(figures_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
  ggsave(file.path(figures_dir, paste0(name, ".pdf")), p, width = w, height = h, device = grDevices::cairo_pdf)
  ggsave(file.path(figures_dir, paste0(name, ".svg")), p, width = w, height = h, device = svglite::svglite)
  ggsave(file.path(figures_dir, paste0(name, ".tiff")), p, width = w, height = h, dpi = 600, device = ragg::agg_tiff, bg = "white")
}
perm_p <- function(x, g) {
  # Exact two-sided label permutation test for the difference in mean paired changes.
  stopifnot(length(x) == 10L, sum(g == "OB") == 5L, sum(g == "LN") == 5L)
  obs <- mean(x[g == "OB"]) - mean(x[g == "LN"])
  ix <- combn(seq_along(x), 5)
  null <- apply(ix, 2, function(j) mean(x[j]) - mean(x[-j]))
  mean(abs(null) >= abs(obs) - 1e-12)
}
paired_p <- function(wk3, wk24) {
  if (all(abs(wk24 - wk3) < 1e-15)) return(1)
  suppressWarnings(wilcox.test(wk24, wk3, paired = TRUE, exact = TRUE)$p.value)
}
unpaired_p <- function(x, group) {
  suppressWarnings(wilcox.test(x[group == "OB"], x[group == "LN"], exact = TRUE)$p.value)
}
unpaired_time_p <- function(wk3, wk24) {
  if (all(abs(wk24 - wk3) < 1e-15)) return(1)
  p <- suppressWarnings(wilcox.test(wk24, wk3, exact = TRUE)$p.value)
  ifelse(is.na(p), 1, p)
}

# 1. Import and validate samples.
raw <- read_tsv(input_file, show_col_types = FALSE, name_repair = "minimal")
stopifnot("clade_name" %in% names(raw))
sample_cols <- setdiff(names(raw), "clade_name")
meta <- tibble(SampleID = sample_cols) %>%
  extract(SampleID, into = c("Diet", "MouseID", "Week"),
          regex = "^(LN|OB)_([0-9]+)_WK(3|24)d_metagenome$", remove = FALSE) %>%
  mutate(Week = as.integer(Week),
         Time = factor(if_else(Week == 3, "WK3_CTRL", "WK24_progression"), levels = c("WK3_CTRL", "WK24_progression")),
         Diet = factor(Diet, levels = c("LN", "OB"))) %>%
  arrange(Diet, MouseID, Week)
stopifnot(nrow(meta) == 20L, all(!is.na(meta$Diet)), sum(meta$Diet == "LN") == 10L, sum(meta$Diet == "OB") == 10L)
pair_check <- meta %>% count(Diet, MouseID) %>% filter(n != 2)
stopifnot(nrow(pair_check) == 0L, n_distinct(meta$MouseID) == 10L)
write_csv(meta, file.path(results_dir, "sample_metadata.csv"))

# 2. Extract t__ rows and tidy abundance matrix (values are MetaPhlAn percentages).
strain <- raw %>% filter(str_detect(clade_name, "(^|\\\\|)t__")) %>%
  mutate(strain_id = str_extract(clade_name, "t__[^|]+$"),
         species = str_extract(clade_name, "s__[^|]+"),
         display_name = if_else(!is.na(species), paste(species, strain_id, sep = " | "), strain_id))
stopifnot(nrow(strain) == 484L, !anyDuplicated(strain$strain_id))
abund_long <- strain %>% select(clade_name, strain_id, display_name, all_of(sample_cols)) %>%
  pivot_longer(all_of(sample_cols), names_to = "SampleID", values_to = "RelativeAbundance_pct") %>%
  left_join(meta, by = "SampleID") %>% mutate(RelativeAbundance_pct = as.numeric(RelativeAbundance_pct))
stopifnot(all(is.finite(abund_long$RelativeAbundance_pct)), all(abund_long$RelativeAbundance_pct >= 0))
write_csv(abund_long, file.path(results_dir, "t_strain_abundance_long_unfiltered.csv"))

# 3. Feature filter; preserve its audit table.
feature_filter <- abund_long %>% group_by(clade_name, strain_id, display_name) %>%
  summarise(detection_n = sum(RelativeAbundance_pct > 0), max_abundance_pct = max(RelativeAbundance_pct), .groups = "drop") %>%
  mutate(retain = detection_n >= 4 & max_abundance_pct >= 0.01,
         filter_reason = if_else(retain, "retained", "detection <4 or maximum abundance <0.01%"))
write_csv(feature_filter, file.path(results_dir, "t_strain_filter_audit.csv"))
keep_ids <- filter(feature_filter, retain) %>% pull(strain_id)
filtered_long <- filter(abund_long, strain_id %in% keep_ids)
write_csv(filtered_long, file.path(results_dir, "t_strain_abundance_long_filtered.csv"))

# 4. Alpha diversity from all t__ features; beta diversity uses closed, zero-replaced t__ composition.
mat_all <- abund_long %>% select(strain_id, SampleID, RelativeAbundance_pct) %>% pivot_wider(names_from = strain_id, values_from = RelativeAbundance_pct) %>% arrange(match(SampleID, meta$SampleID))
sample_ids <- mat_all$SampleID; x_all <- as.matrix(select(mat_all, -SampleID)); rownames(x_all) <- sample_ids
alpha <- tibble(SampleID = sample_ids, Richness = rowSums(x_all > 0), Shannon = diversity(x_all, index = "shannon")) %>% left_join(meta, by = "SampleID")
write_csv(alpha, file.path(results_dir, "alpha_diversity_t_strains.csv"))

zrep <- function(x) {
  x <- x / sum(x)
  pos <- x[x > 0]
  replacement <- min(pos) / 2
  x[x == 0] <- replacement
  x / sum(x)
}
comp <- t(apply(x_all, 1, zrep)); clr <- log(comp) - rowMeans(log(comp))
pca <- prcomp(clr, center = FALSE, scale. = FALSE)
pca_df <- as_tibble(pca$x[, 1:2], rownames = "SampleID") %>% left_join(meta, by = "SampleID")
var_exp <- 100 * pca$sdev^2 / sum(pca$sdev^2)
write_csv(pca_df, file.path(results_dir, "t_strain_aitchison_pca_scores.csv"))

# 5. Two transformations are analysed identically: log10 and AST (arcsine square-root).
pseudocount_pct <- 0.001
make_transformed <- function(dat, transform) {
  dat %>% mutate(value = if (transform == "log10") log10(RelativeAbundance_pct + pseudocount_pct) else asin(sqrt(RelativeAbundance_pct / 100)))
}
run_feature_tests <- function(transform) {
  d <- make_transformed(filtered_long, transform)
  # SampleID identifies the individual time-point; omit it so WK3/WK24 combine into one mouse row.
  wide <- d %>% select(strain_id, display_name, clade_name, Diet, MouseID, Week, value, RelativeAbundance_pct) %>%
    pivot_wider(names_from = Week, values_from = c(value, RelativeAbundance_pct), names_sep = "_WK")
  wide <- wide %>% mutate(delta = value_WK24 - value_WK3)
  stats <- wide %>% group_by(strain_id, display_name, clade_name) %>% group_modify(~{
    z <- .x
    wk3 <- z %>% transmute(value = value_WK3, Diet)
    wk24 <- z %>% transmute(value = value_WK24, Diet)
    tibble(
      transform = transform,
      interaction_effect_OB_minus_LN = mean(z$delta[z$Diet == "OB"]) - mean(z$delta[z$Diet == "LN"]),
      interaction_p = perm_p(z$delta, as.character(z$Diet)),
      paired_LN_effect = mean(z$delta[z$Diet == "LN"]), paired_LN_p = paired_p(z$value_WK3[z$Diet == "LN"], z$value_WK24[z$Diet == "LN"]),
      paired_OB_effect = mean(z$delta[z$Diet == "OB"]), paired_OB_p = paired_p(z$value_WK3[z$Diet == "OB"], z$value_WK24[z$Diet == "OB"]),
      # Sensitivity analysis requested by the user: treat WK3 and WK24 in OB as independent groups.
      unpaired_OB_time_effect = mean(z$value_WK24[z$Diet == "OB"]) - mean(z$value_WK3[z$Diet == "OB"]),
      unpaired_OB_time_p = unpaired_time_p(z$value_WK3[z$Diet == "OB"], z$value_WK24[z$Diet == "OB"]),
      baseline_OB_minus_LN = mean(wk3$value[wk3$Diet == "OB"]) - mean(wk3$value[wk3$Diet == "LN"]), baseline_p = unpaired_p(wk3$value, as.character(wk3$Diet)),
      endpoint_OB_minus_LN = mean(wk24$value[wk24$Diet == "OB"]) - mean(wk24$value[wk24$Diet == "LN"]), endpoint_p = unpaired_p(wk24$value, as.character(wk24$Diet)),
      LN_WK3_mean_pct = mean(z$RelativeAbundance_pct_WK3[z$Diet == "LN"]), LN_WK3_median_pct = median(z$RelativeAbundance_pct_WK3[z$Diet == "LN"]),
      LN_WK24_mean_pct = mean(z$RelativeAbundance_pct_WK24[z$Diet == "LN"]), LN_WK24_median_pct = median(z$RelativeAbundance_pct_WK24[z$Diet == "LN"]),
      OB_WK3_mean_pct = mean(z$RelativeAbundance_pct_WK3[z$Diet == "OB"]), OB_WK3_median_pct = median(z$RelativeAbundance_pct_WK3[z$Diet == "OB"]),
      OB_WK24_mean_pct = mean(z$RelativeAbundance_pct_WK24[z$Diet == "OB"]), OB_WK24_median_pct = median(z$RelativeAbundance_pct_WK24[z$Diet == "OB"]),
      detection_n = sum(c(z$RelativeAbundance_pct_WK3, z$RelativeAbundance_pct_WK24) > 0)
    )
  }) %>% ungroup() %>% mutate(interaction_fdr = p.adjust(interaction_p, "BH"), paired_LN_fdr = p.adjust(paired_LN_p, "BH"), paired_OB_fdr = p.adjust(paired_OB_p, "BH"), unpaired_OB_time_fdr = p.adjust(unpaired_OB_time_p, "BH"),
                                interaction_label = case_when(interaction_fdr < 0.05 ~ "FDR < 0.05", interaction_fdr < 0.10 ~ "FDR < 0.10", TRUE ~ "not significant")) %>%
    arrange(interaction_fdr, interaction_p, desc(abs(interaction_effect_OB_minus_LN)))
  list(stats = stats, wide = wide)
}
log_res <- run_feature_tests("log10"); ast_res <- run_feature_tests("AST")
write_csv(log_res$stats, file.path(results_dir, "t_strain_differential_log10.csv"))
write_csv(ast_res$stats, file.path(results_dir, "t_strain_differential_AST.csv"))
combined <- bind_rows(log_res$stats, ast_res$stats)
write_csv(combined, file.path(results_dir, "t_strain_differential_both_transformations.csv"))

# 6. Alpha diversity statistics and source table.
alpha_stats <- alpha %>% pivot_longer(c(Richness, Shannon), names_to = "metric", values_to = "value") %>% group_by(metric) %>% group_modify(~{
  # Pair two time points by mouse, not by the time-point-specific SampleID.
  z <- .x %>% select(Diet, MouseID, Week, value) %>% pivot_wider(names_from = Week, values_from = value, names_prefix = "WK") %>% mutate(delta = WK24 - WK3)
  tibble(interaction_effect_OB_minus_LN = mean(z$delta[z$Diet == "OB"]) - mean(z$delta[z$Diet == "LN"]), interaction_p = perm_p(z$delta, as.character(z$Diet)),
         paired_LN_p = paired_p(z$WK3[z$Diet == "LN"], z$WK24[z$Diet == "LN"]), paired_OB_p = paired_p(z$WK3[z$Diet == "OB"], z$WK24[z$Diet == "OB"]))
}) %>% ungroup() %>% mutate(interaction_fdr = p.adjust(interaction_p, "BH"))
write_csv(alpha_stats, file.path(results_dir, "alpha_diversity_statistics.csv"))

# 7. Figures: alpha, PCA, abundance overview, interaction concordance and candidate detail.
p_alpha <- alpha %>% pivot_longer(c(Richness, Shannon), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(Time, Value, group = MouseID, colour = Diet)) + geom_line(alpha = .6) + geom_point(size = 2.4) +
  facet_wrap(~Metric, scales = "free_y") + scale_colour_manual(values = c(LN = ln_col, OB = ob_col)) + labs(title = "t__ strain alpha diversity", x = NULL, y = NULL, colour = "Diet")
save_plot(p_alpha, "01_alpha_diversity_paired", 8, 4.5)

p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = Diet, shape = Time, group = MouseID)) + geom_line(colour = "grey65", linewidth = .45) + geom_point(size = 3) +
  scale_colour_manual(values = c(LN = ln_col, OB = ob_col)) + scale_shape_manual(values = c(WK3_CTRL = 16, WK24_progression = 17)) +
  labs(title = "Aitchison PCA of t__ strain composition", x = sprintf("PC1 (%.1f%%)", var_exp[1]), y = sprintf("PC2 (%.1f%%)", var_exp[2]), colour = "Diet", shape = "Time")
save_plot(p_pca, "02_aitchison_pca_paired", 7, 5.5)

top_ids <- abund_long %>% group_by(strain_id, display_name) %>% summarise(mean_pct = mean(RelativeAbundance_pct), .groups = "drop") %>% slice_max(mean_pct, n = 15) %>% pull(strain_id)
bar_df <- abund_long %>% filter(strain_id %in% top_ids) %>% mutate(label = display_name) %>% group_by(SampleID, Diet, Time, label) %>% summarise(value = sum(RelativeAbundance_pct), .groups = "drop") %>%
  bind_rows(abund_long %>% filter(!strain_id %in% top_ids) %>% group_by(SampleID, Diet, Time) %>% summarise(value = sum(RelativeAbundance_pct), .groups = "drop") %>% mutate(label = "Other t__ strains")) %>%
  left_join(meta %>% select(SampleID, MouseID), by = "SampleID") %>% arrange(Diet, MouseID, Time)
bar_df$SampleID <- factor(bar_df$SampleID, levels = unique(bar_df$SampleID))
p_bar <- ggplot(bar_df, aes(SampleID, value, fill = label)) + geom_col(width = .85) + facet_grid(~Diet, scales = "free_x", space = "free_x") +
  labs(title = "Top t__ strain composition", x = NULL, y = "Relative abundance (%)", fill = "Strain") + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6), legend.text = element_text(size = 6))
save_plot(p_bar, "03_top_strain_composition", 11, 5.5)

rank_log <- log_res$stats %>% select(strain_id, display_name, log_effect = interaction_effect_OB_minus_LN, log_p = interaction_p, log_fdr = interaction_fdr)
rank_ast <- ast_res$stats %>% select(strain_id, ast_effect = interaction_effect_OB_minus_LN, ast_p = interaction_p, ast_fdr = interaction_fdr)
conc <- left_join(rank_log, rank_ast, by = "strain_id") %>% mutate(label = if_else(rank(pmin(log_fdr, ast_fdr), ties.method = "first") <= 6, display_name, NA_character_), label_hjust = if_else(log_effect < 0, 0, 1))
p_conc <- ggplot(conc, aes(log_effect, ast_effect, colour = pmin(log_fdr, ast_fdr) < .10)) + geom_hline(yintercept = 0, colour = "grey75") + geom_vline(xintercept = 0, colour = "grey75") + geom_point(size = 2, alpha = .8) +
  geom_text(data = filter(conc, !is.na(label)), aes(label = label, hjust = label_hjust), size = 2.5, show.legend = FALSE, check_overlap = TRUE) +
  scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "#7A3E9D"), name = "Either FDR < 0.10") + labs(title = "Interaction-effect concordance: log10 versus AST", x = "log10: OB minus LN change", y = "AST: OB minus LN change") + theme(plot.margin = margin(5.5, 35, 5.5, 5.5))
save_plot(p_conc, "04_transform_concordance", 8, 5.5)

candidate_ids <- conc %>% arrange(pmin(log_fdr, ast_fdr), pmin(log_p, ast_p)) %>% slice_head(n = 12) %>% pull(strain_id)
candidate_names <- feature_filter %>% filter(strain_id %in% candidate_ids) %>% select(strain_id, display_name)
candidate_long <- filtered_long %>% filter(strain_id %in% candidate_ids) %>% left_join(candidate_names, by = "strain_id", suffix = c("", "_candidate")) %>% mutate(display_plot = coalesce(display_name_candidate, display_name))
p_candidates <- ggplot(candidate_long, aes(Time, RelativeAbundance_pct, group = MouseID, colour = Diet)) + geom_line(alpha = .55) + geom_point(size = 1.8) +
  facet_wrap(~display_plot, scales = "free_y", ncol = 3) + scale_colour_manual(values = c(LN = ln_col, OB = ob_col)) +
  labs(title = "Top interaction candidates: paired relative abundances", x = NULL, y = "Relative abundance (%)", colour = "Diet") + theme(strip.text = element_text(size = 7))
save_plot(p_candidates, "05_top_candidates_paired", 10, 8)

heat_mat <- candidate_long %>% select(strain_id, SampleID, RelativeAbundance_pct) %>% pivot_wider(names_from = SampleID, values_from = RelativeAbundance_pct) %>% tibble::column_to_rownames("strain_id") %>% as.matrix()
ordered_samples <- meta %>% arrange(Diet, MouseID, Week) %>% pull(SampleID)
heat_mat <- log10(heat_mat[, ordered_samples, drop = FALSE] + pseudocount_pct)
ann <- meta %>% arrange(match(SampleID, ordered_samples)) %>% select(Diet, Time, MouseID) %>% as.data.frame(); rownames(ann) <- ordered_samples
grDevices::cairo_pdf(file.path(figures_dir, "06_top_candidates_heatmap.pdf"), width = 10, height = 6)
pheatmap(heat_mat, annotation_col = ann, cluster_cols = FALSE, cluster_rows = TRUE, main = "Top interaction candidates (log10 abundance)", fontsize_row = 7)
dev.off()
svglite::svglite(file.path(figures_dir, "06_top_candidates_heatmap.svg"), width = 10, height = 6)
pheatmap(heat_mat, annotation_col = ann, cluster_cols = FALSE, cluster_rows = TRUE, main = "Top interaction candidates (log10 abundance)", fontsize_row = 7)
dev.off()
ragg::agg_tiff(file.path(figures_dir, "06_top_candidates_heatmap.tiff"), width = 10, height = 6, units = "in", res = 600)
pheatmap(heat_mat, annotation_col = ann, cluster_cols = FALSE, cluster_rows = TRUE, main = "Top interaction candidates (log10 abundance)", fontsize_row = 7)
dev.off()
png(file.path(figures_dir, "06_top_candidates_heatmap.png"), width = 6000, height = 3600, res = 600)
pheatmap(heat_mat, annotation_col = ann, cluster_cols = FALSE, cluster_rows = TRUE, main = "Top interaction candidates (log10 abundance)", fontsize_row = 7)
dev.off()

# 8. Manual audit for the top log10 candidate; the same values used in the interaction test are exported.
audit_id <- log_res$stats$strain_id[1]
audit <- log_res$wide %>% filter(strain_id == audit_id) %>% select(strain_id, display_name, Diet, MouseID, value_WK3, value_WK24, delta, RelativeAbundance_pct_WK3, RelativeAbundance_pct_WK24)
write_csv(audit, file.path(results_dir, "manual_audit_top_log10_candidate.csv"))

summary_txt <- c(
  "MetaPhlAn t__ paired longitudinal analysis",
  sprintf("Input t__ features: %d; retained features: %d.", nrow(feature_filter), sum(feature_filter$retain)),
  sprintf("Samples: %d (LN=%d, OB=%d); paired mice: %d.", nrow(meta), sum(meta$Diet == "LN"), sum(meta$Diet == "OB"), n_distinct(meta$MouseID)),
  sprintf("Pseudo-count for log10 transformation: %.3f percentage points.", pseudocount_pct),
  "Primary test: exact two-sided permutation test comparing WK24-WK3 changes in OB versus LN (252 label allocations).",
  "Sensitivity transformations: log10(relative abundance + pseudo-count) and AST=asin(sqrt(relative abundance / 100)).",
  sprintf("Log10: %d features at FDR<0.05; %d at FDR<0.10.", sum(log_res$stats$interaction_fdr < .05), sum(log_res$stats$interaction_fdr < .10)),
  sprintf("AST: %d features at FDR<0.05; %d at FDR<0.10.", sum(ast_res$stats$interaction_fdr < .05), sum(ast_res$stats$interaction_fdr < .10))
)
writeLines(summary_txt, file.path(results_dir, "analysis_summary.txt"))
writeLines(c("Core conclusion: quantify whether longitudinal t__ strain changes differ between obese and lean mice.", "Figure archetype: quantitative grid.", "Hero evidence: exact-permutation interaction effects, shown for both log10 and AST transformations.", "Supporting evidence: paired alpha diversity, Aitchison PCA, composition overview and candidate trajectories.", "Statistics: two-sided exact permutation for interaction; paired Wilcoxon tests within diet; BH correction per transformation.", "Data integrity: all 20 samples and all 484 t__ rows retained in source outputs; formal filter is explicitly audited."), file.path(root, "figure_contract.txt"))
message("Analysis complete: ", root)
