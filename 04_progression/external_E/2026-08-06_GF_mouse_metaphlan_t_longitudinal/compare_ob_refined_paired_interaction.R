#!/usr/bin/env Rscript
# Primary analysis: exact test of the Diet × Time interaction, implemented as
# comparison of each donor's paired WK24-WK3 change (OB versus LN).  Counts of
# OB donors whose within-donor direction matches the historical CRC direction
# are descriptive only (3/5, 4/5, 5/5), not additional significance tests.
suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(stringr) })

root <- normalizePath(file.path(getwd(), "analyses", "2026-08-06_GF_mouse_metaphlan_t_longitudinal"), mustWork = TRUE)
results_dir <- file.path(root, "results")
ref <- read_csv(file.path(getwd(), "analyses", "2026-07-13_bmi_union_overlap_heatmaps", "tables", "Feature_List_Ob_Refined.csv"), show_col_types = FALSE) %>%
  mutate(strain_id = str_extract(feature, "t__[^|]+$"), old_crc_direction = if_else(coef_crc > 0, "CRC enriched", "CRC depleted"))
stats <- read_csv(file.path(results_dir, "t_strain_differential_both_transformations.csv"), show_col_types = FALSE)
abundance <- read_csv(file.path(results_dir, "t_strain_abundance_long_unfiltered.csv"), show_col_types = FALSE)

ob_deltas <- abundance %>% filter(Diet == "OB") %>%
  select(strain_id, MouseID, Week, RelativeAbundance_pct) %>%
  crossing(transform = c("AST", "log10")) %>%
  mutate(value = if_else(transform == "AST", asin(sqrt(RelativeAbundance_pct / 100)), log10(RelativeAbundance_pct + 0.001))) %>%
  select(-RelativeAbundance_pct) %>% pivot_wider(names_from = Week, values_from = value, names_prefix = "WK") %>%
  mutate(ob_donor_delta = WK24 - WK3)

comparison <- inner_join(ref, stats, by = "strain_id") %>%
  transmute(strain_id, clean_name, overlap, coef_crc, qval_crc, old_crc_direction, transform,
            interaction_effect_OB_minus_LN, interaction_p, interaction_fdr,
            interaction_significant_user_rule = interaction_p < 0.01 | interaction_fdr < 0.10,
            interaction_direction_consistent = sign(coef_crc) == sign(interaction_effect_OB_minus_LN)) %>%
  left_join(ob_deltas, by = c("strain_id", "transform")) %>%
  mutate(ob_donor_direction_consistent = sign(ob_donor_delta) == sign(coef_crc)) %>%
  group_by(strain_id, clean_name, overlap, coef_crc, qval_crc, old_crc_direction, transform,
           interaction_effect_OB_minus_LN, interaction_p, interaction_fdr,
           interaction_significant_user_rule, interaction_direction_consistent) %>%
  summarise(n_OB_donors_concordant = sum(ob_donor_direction_consistent, na.rm = TRUE), n_OB_donors = sum(!is.na(ob_donor_delta)),
            donor_direction_support = paste0(n_OB_donors_concordant, "/", n_OB_donors), .groups = "drop") %>%
  arrange(transform, interaction_fdr, interaction_p, desc(n_OB_donors_concordant))

summary <- comparison %>% group_by(transform) %>%
  summarise(overlap_n = n(), interaction_significant_n = sum(interaction_significant_user_rule),
            interaction_significant_and_direction_consistent_n = sum(interaction_significant_user_rule & interaction_direction_consistent),
            interaction_direction_consistent_all_n = sum(interaction_direction_consistent),
            descriptive_3of5_n = sum(n_OB_donors_concordant >= 3), descriptive_4of5_n = sum(n_OB_donors_concordant >= 4), descriptive_5of5_n = sum(n_OB_donors_concordant == 5),
            significant_consistent_with_3of5_n = sum(interaction_significant_user_rule & interaction_direction_consistent & n_OB_donors_concordant >= 3),
            significant_consistent_with_4of5_n = sum(interaction_significant_user_rule & interaction_direction_consistent & n_OB_donors_concordant >= 4),
            significant_consistent_with_5of5_n = sum(interaction_significant_user_rule & interaction_direction_consistent & n_OB_donors_concordant == 5), .groups = "drop")

write_csv(comparison, file.path(results_dir, "overlap_Ob_Refined_vs_paired_interaction_with_OB_donor_support.csv"))
write_csv(summary, file.path(results_dir, "overlap_Ob_Refined_vs_paired_interaction_summary.csv"))

# Strict biological-direction table: CRC coefficient and the *OB within-donor*
# WK24-WK3 mean change must agree in both AST and log10. This excludes taxa
# whose interaction is driven only by a change in LN (for example, SGB1867).
ob_direction <- abundance %>% filter(Diet == "OB") %>%
  select(strain_id, MouseID, Week, RelativeAbundance_pct) %>%
  crossing(transform = c("AST", "log10")) %>%
  mutate(value = if_else(transform == "AST", asin(sqrt(RelativeAbundance_pct / 100)), log10(RelativeAbundance_pct + 0.001))) %>%
  select(-RelativeAbundance_pct) %>% pivot_wider(names_from = Week, values_from = value, names_prefix = "WK") %>%
  mutate(ob_delta = WK24 - WK3) %>% group_by(strain_id, transform) %>%
  summarise(ob_mean_delta = mean(ob_delta), ob_donors_up = sum(ob_delta > 0), ob_donors_down = sum(ob_delta < 0), ob_donors_zero = sum(ob_delta == 0), .groups = "drop")

strict <- comparison %>% left_join(ob_direction, by = c("strain_id", "transform")) %>%
  mutate(ob_crc_direction_consistent = sign(coef_crc) == sign(ob_mean_delta)) %>%
  group_by(strain_id, clean_name, overlap, coef_crc, qval_crc, old_crc_direction) %>%
  summarise(AST_interaction_effect = first(interaction_effect_OB_minus_LN[transform == "AST"]), AST_interaction_p = first(interaction_p[transform == "AST"]), AST_interaction_fdr = first(interaction_fdr[transform == "AST"]), AST_OB_mean_delta = first(ob_mean_delta[transform == "AST"]), AST_OB_up = first(ob_donors_up[transform == "AST"]), AST_OB_down = first(ob_donors_down[transform == "AST"]), AST_direction_consistent = first(ob_crc_direction_consistent[transform == "AST"]),
            log10_interaction_effect = first(interaction_effect_OB_minus_LN[transform == "log10"]), log10_interaction_p = first(interaction_p[transform == "log10"]), log10_interaction_fdr = first(interaction_fdr[transform == "log10"]), log10_OB_mean_delta = first(ob_mean_delta[transform == "log10"]), log10_OB_up = first(ob_donors_up[transform == "log10"]), log10_OB_down = first(ob_donors_down[transform == "log10"]), log10_direction_consistent = first(ob_crc_direction_consistent[transform == "log10"]),
            interaction_significant = any(interaction_significant_user_rule), .groups = "drop") %>%
  mutate(selected_strict_OB_CRC_direction = interaction_significant & AST_direction_consistent & log10_direction_consistent,
         exclusion_reason = case_when(strain_id == "t__SGB1867" ~ "Excluded: OB is unchanged (5/5 zero); interaction is driven by LN decrease", strain_id == "t__SGB15132" ~ "Excluded: 4/5 OB donors decrease despite positive CRC coefficient", TRUE ~ NA_character_))
write_csv(strict, file.path(results_dir, "Ob_Refined_paired_interaction_strict_OB_CRC_direction.csv"))
write_csv(filter(strict, selected_strict_OB_CRC_direction), file.path(results_dir, "selected_12_strict_OB_CRC_direction_taxa.csv"))
print(summary)
