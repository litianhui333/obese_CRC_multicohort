rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

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

comparison <- "CTRL_vs_CRC"
regions <- c("Asia", "Europe")
bmi_groups <- c("Normal", "Ow", "Ob")
bmi_short <- c("Normal" = "Nw", "Ow" = "Ow", "Ob" = "Ob")

robust_p_cutoff <- 0.01
robust_q_cutoff <- 0.10
nominal_p_cutoff <- 0.05

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

safe_write_csv <- function(x, filename) {
  tryCatch({
    write_csv(x, filename)
  }, error = function(e) {
    fallback <- sub("(\\.[^.]+)$", "_new\\1", filename)
    message("Could not overwrite ", filename, "; writing ", fallback, " instead.")
    write_csv(x, fallback)
  })
}

extract_meta <- function(region, bmi_group, features) {
  list_name <- paste0("disease_analysis_", bmi_group)
  res <- tryCatch(
    all_analysis_results[[region]][["disease_analysis"]][["analysis_results"]][[list_name]][["all_results"]][[comparison]][["meta_fits"]],
    error = function(e) NULL
  )
  if (is.null(res)) return(tibble())

  df <- as.data.frame(res)
  if (!"feature" %in% names(df)) {
    df <- df %>% tibble::rownames_to_column("feature")
  }
  q_col <- if ("qval.fdr" %in% names(df)) "qval.fdr" else if ("qval" %in% names(df)) "qval" else NA_character_

  df %>%
    filter(feature %in% features) %>%
    transmute(
      feature,
      region,
      BMI_Group = bmi_group,
      coef = coef,
      pval = pval,
      qval = if (!is.na(q_col)) .data[[q_col]] else NA_real_
    )
}

classify_pattern <- function(robust_vec, nominal_vec) {
  robust_vec[is.na(robust_vec)] <- FALSE
  nominal_vec[is.na(nominal_vec)] <- FALSE

  r_nw <- robust_vec["Normal"]
  r_ow <- robust_vec["Ow"]
  r_ob <- robust_vec["Ob"]
  n_nw <- nominal_vec["Normal"]
  n_ow <- nominal_vec["Ow"]
  n_ob <- nominal_vec["Ob"]

  case_when(
    r_nw & !n_ow & !n_ob ~ "Nw-specific",
    r_ow & !n_nw & !n_ob ~ "Ow-specific",
    r_ob & !n_nw & !n_ow ~ "Ob-specific",
    r_nw & r_ow & !n_ob ~ "Nw+Ow shared",
    r_nw & r_ob & !n_ow ~ "Nw+Ob shared",
    r_ow & r_ob & !n_nw ~ "Ow+Ob shared",
    r_nw & r_ow & r_ob ~ "BMI-independent",
    (n_nw | n_ow | n_ob) ~ "Weak/non-target signal",
    TRUE ~ "No BMI signal"
  )
}

feature_files <- c(
  Normal = "Feature_List_Normal_Refined.csv",
  Ow = "Feature_List_Ow_Refined.csv",
  Ob = "Feature_List_Ob_Refined.csv"
)

overall_long <- bind_rows(lapply(names(feature_files), function(bmi) {
  read_csv(file.path(heatmap_dir, feature_files[[bmi]]), show_col_types = FALSE) %>%
    transmute(feature, clean_name, overall_source_BMI = bmi, overlap, coef_crc = coef_crc)
}))

overall_sets <- overall_long %>%
  distinct(feature, clean_name, overall_source_BMI) %>%
  group_by(feature, clean_name) %>%
  summarise(
    overall_BMI_set = paste(sort(bmi_short[overall_source_BMI]), collapse = "+"),
    .groups = "drop"
  ) %>%
  mutate(
    overall_pattern = case_when(
      overall_BMI_set == "Nw" ~ "Nw-specific",
      overall_BMI_set == "Ow" ~ "Ow-specific",
      overall_BMI_set == "Ob" ~ "Ob-specific",
      overall_BMI_set == "Nw+Ow" ~ "Nw+Ow shared",
      overall_BMI_set == "Nw+Ob" ~ "Nw+Ob shared",
      overall_BMI_set == "Ob+Ow" ~ "Ow+Ob shared",
      overall_BMI_set == "Nw+Ob+Ow" ~ "BMI-independent",
      TRUE ~ overall_BMI_set
    )
  )

target_features <- overall_sets$feature

regional_long <- bind_rows(lapply(regions, function(region) {
  bind_rows(lapply(bmi_groups, function(bmi) {
    extract_meta(region, bmi, target_features)
  }))
})) %>%
  mutate(
    robust = (!is.na(pval) & pval < robust_p_cutoff) | (!is.na(qval) & qval < robust_q_cutoff),
    nominal = !is.na(pval) & pval < nominal_p_cutoff,
    nominal_only_support = nominal,
    coef_direction = case_when(coef > 0 ~ "positive", coef < 0 ~ "negative", TRUE ~ "zero_or_missing")
  )

regional_patterns <- regional_long %>%
  group_by(feature, region) %>%
  summarise(
    regional_pattern = classify_pattern(setNames(robust, BMI_Group), setNames(nominal, BMI_Group)),
    robust_BMI_set = paste(bmi_short[BMI_Group[robust]], collapse = "+"),
    nominal_BMI_set = paste(bmi_short[BMI_Group[nominal]], collapse = "+"),
    .groups = "drop"
  )

regional_support_metrics <- regional_long %>%
  left_join(overall_sets %>% select(feature, overall_pattern), by = "feature") %>%
  left_join(
    overall_long %>%
      select(feature, overall_source_BMI, coef_crc) %>%
      mutate(overall_source_BMI = as.character(overall_source_BMI)) %>%
      distinct(),
    by = c("feature", "BMI_Group" = "overall_source_BMI")
  ) %>%
  group_by(feature, region) %>%
  summarise(
    overall_pattern = first(overall_pattern),
    coef_Normal = coef[BMI_Group == "Normal"][1],
    coef_Ow = coef[BMI_Group == "Ow"][1],
    coef_Ob = coef[BMI_Group == "Ob"][1],
    pval_Normal = pval[BMI_Group == "Normal"][1],
    pval_Ow = pval[BMI_Group == "Ow"][1],
    pval_Ob = pval[BMI_Group == "Ob"][1],
    robust_Normal = robust[BMI_Group == "Normal"][1],
    robust_Ow = robust[BMI_Group == "Ow"][1],
    robust_Ob = robust[BMI_Group == "Ob"][1],
    nominal_Normal = nominal[BMI_Group == "Normal"][1],
    nominal_Ow = nominal[BMI_Group == "Ow"][1],
    nominal_Ob = nominal[BMI_Group == "Ob"][1],
    overall_coef_Normal = coef_crc[BMI_Group == "Normal"][1],
    overall_coef_Ow = coef_crc[BMI_Group == "Ow"][1],
    overall_coef_Ob = coef_crc[BMI_Group == "Ob"][1],
    .groups = "drop"
  ) %>%
  mutate(
    target_bmis = case_when(
      overall_pattern == "Nw-specific" ~ "Normal",
      overall_pattern == "Ow-specific" ~ "Ow",
      overall_pattern == "Ob-specific" ~ "Ob",
      overall_pattern == "Nw+Ow shared" ~ "Normal;Ow",
      overall_pattern == "Nw+Ob shared" ~ "Normal;Ob",
      overall_pattern == "Ow+Ob shared" ~ "Ow;Ob",
      overall_pattern == "BMI-independent" ~ "Normal;Ow;Ob",
      TRUE ~ NA_character_
    )
  ) %>%
  rowwise() %>%
  mutate(
    coef_vec = list(c(Normal = coef_Normal, Ow = coef_Ow, Ob = coef_Ob)),
    overall_coef_vec = list(c(Normal = overall_coef_Normal, Ow = overall_coef_Ow, Ob = overall_coef_Ob)),
    robust_vec = list(c(Normal = robust_Normal, Ow = robust_Ow, Ob = robust_Ob)),
    nominal_vec = list(c(Normal = nominal_Normal, Ow = nominal_Ow, Ob = nominal_Ob)),
    target_vec = list(strsplit(target_bmis, ";")[[1]]),
    non_target_vec = list(setdiff(bmi_groups, target_vec)),
    max_abs_coef = suppressWarnings(max(abs(unlist(coef_vec)), na.rm = TRUE)),
    target_coef_values = list(unlist(coef_vec)[target_vec]),
    target_signs = list(sign(target_coef_values)),
    target_has_nominal = all(unlist(nominal_vec)[target_vec], na.rm = TRUE),
    target_has_signal = any(unlist(nominal_vec)[target_vec], na.rm = TRUE),
    target_dominant = ifelse(
      is.finite(max_abs_coef) && max_abs_coef > 0,
      max(abs(target_coef_values), na.rm = TRUE) >= 0.7 * max_abs_coef,
      FALSE
    ),
    target_same_direction = {
      vals <- target_coef_values
      vals <- vals[!is.na(vals) & vals != 0]
      length(vals) > 0 && length(unique(sign(vals))) == 1
    },
    non_target_same_direction_robust = {
      if (length(non_target_vec) == 0) {
        FALSE
      } else {
        non_vals <- unlist(coef_vec)[non_target_vec]
        non_rob <- unlist(robust_vec)[non_target_vec]
        tar_sign <- unique(sign(target_coef_values[!is.na(target_coef_values) & target_coef_values != 0]))
        if (length(tar_sign) != 1) {
          FALSE
        } else {
          any(non_rob & sign(non_vals) == tar_sign, na.rm = TRUE)
        }
      }
    },
    non_target_same_direction_nominal = {
      if (length(non_target_vec) == 0) {
        FALSE
      } else {
        non_vals <- unlist(coef_vec)[non_target_vec]
        non_nom <- unlist(nominal_vec)[non_target_vec]
        tar_sign <- unique(sign(target_coef_values[!is.na(target_coef_values) & target_coef_values != 0]))
        if (length(tar_sign) != 1) {
          FALSE
        } else {
          any(non_nom & sign(non_vals) == tar_sign, na.rm = TRUE)
        }
      }
    },
    relaxed_pattern_support = case_when(
      overall_pattern %in% c("Nw-specific", "Ow-specific", "Ob-specific") ~
        target_has_signal & target_dominant & target_same_direction & !non_target_same_direction_robust,
      overall_pattern %in% c("Nw+Ow shared", "Nw+Ob shared", "Ow+Ob shared") ~
        target_has_signal & target_dominant & target_same_direction & !non_target_same_direction_robust,
      overall_pattern == "BMI-independent" ~
        target_has_signal & sum(unlist(nominal_vec), na.rm = TRUE) >= 2,
      TRUE ~ FALSE
    ),
    relaxed_pattern_support_p05 = case_when(
      overall_pattern %in% c("Nw-specific", "Ow-specific", "Ob-specific") ~ {
        target_nom <- all(unlist(nominal_vec)[target_vec], na.rm = TRUE)
        non_target_non_nom <- all(!unlist(nominal_vec)[non_target_vec], na.rm = TRUE)
        target_same_as_overall <- all(
          sign(unlist(coef_vec)[target_vec]) == sign(unlist(overall_coef_vec)[target_vec]),
          na.rm = TRUE
        )
        target_nom & non_target_non_nom & target_same_as_overall
      },
      overall_pattern %in% c("Nw+Ow shared", "Nw+Ob shared", "Ow+Ob shared") ~ {
        target_nom <- all(unlist(nominal_vec)[target_vec], na.rm = TRUE)
        non_target_non_nom <- all(!unlist(nominal_vec)[non_target_vec], na.rm = TRUE)
        target_same_as_overall <- all(
          sign(unlist(coef_vec)[target_vec]) == sign(unlist(overall_coef_vec)[target_vec]),
          na.rm = TRUE
        )
        target_nom & non_target_non_nom & target_same_as_overall
      },
      overall_pattern == "BMI-independent" ~ {
        target_nom <- all(unlist(nominal_vec)[target_vec], na.rm = TRUE)
        target_same_as_overall <- all(
          sign(unlist(coef_vec)[target_vec]) == sign(unlist(overall_coef_vec)[target_vec]),
          na.rm = TRUE
        )
        target_nom & target_same_as_overall
      },
      TRUE ~ FALSE
    )
  ) %>%
  ungroup() %>%
  select(
    feature, region, relaxed_pattern_support, relaxed_pattern_support_p05, target_bmis, max_abs_coef,
    coef_Normal, coef_Ow, coef_Ob,
    pval_Normal, pval_Ow, pval_Ob,
    starts_with("nominal_"), starts_with("robust_")
  )

pattern_consistency <- overall_sets %>%
  left_join(regional_patterns, by = "feature") %>%
  left_join(regional_support_metrics, by = c("feature", "region")) %>%
  mutate(
    exact_pattern_match = overall_pattern == regional_pattern,
    broad_pattern_support = case_when(
      overall_pattern == "Nw-specific" ~ regional_pattern %in% c("Nw-specific", "Nw+Ow shared", "Nw+Ob shared", "BMI-independent"),
      overall_pattern == "Ow-specific" ~ regional_pattern %in% c("Ow-specific", "Nw+Ow shared", "Ow+Ob shared", "BMI-independent"),
      overall_pattern == "Ob-specific" ~ regional_pattern %in% c("Ob-specific", "Nw+Ob shared", "Ow+Ob shared", "BMI-independent"),
      overall_pattern == "Nw+Ow shared" ~ regional_pattern %in% c("Nw+Ow shared", "BMI-independent"),
      overall_pattern == "Nw+Ob shared" ~ regional_pattern %in% c("Nw+Ob shared", "BMI-independent"),
      overall_pattern == "Ow+Ob shared" ~ regional_pattern %in% c("Ow+Ob shared", "BMI-independent"),
      overall_pattern == "BMI-independent" ~ regional_pattern == "BMI-independent",
      TRUE ~ FALSE
    )
  )

summary_by_pattern <- pattern_consistency %>%
  group_by(overall_pattern, region) %>%
  summarise(
    n_taxa = n(),
    exact_match = sum(exact_pattern_match, na.rm = TRUE),
    exact_match_pct = exact_match / n_taxa,
    broad_support = sum(broad_pattern_support, na.rm = TRUE),
    broad_support_pct = broad_support / n_taxa,
    relaxed_support = sum(relaxed_pattern_support, na.rm = TRUE),
    relaxed_support_pct = relaxed_support / n_taxa,
    relaxed_p05_support = sum(relaxed_pattern_support_p05, na.rm = TRUE),
    relaxed_p05_support_pct = relaxed_p05_support / n_taxa,
    no_bmi_signal = sum(regional_pattern == "No BMI signal", na.rm = TRUE),
    weak_or_other = sum(regional_pattern == "Weak/non-target signal", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(
    factor(overall_pattern, levels = c("Nw-specific", "Ow-specific", "Ob-specific", "Nw+Ow shared", "Nw+Ob shared", "Ow+Ob shared", "BMI-independent")),
    region
  )

summary_any_region <- pattern_consistency %>%
  group_by(feature, clean_name, overall_pattern) %>%
  summarise(
    exact_in_any_region = any(exact_pattern_match, na.rm = TRUE),
    broad_in_any_region = any(broad_pattern_support, na.rm = TRUE),
    relaxed_in_any_region = any(relaxed_pattern_support, na.rm = TRUE),
    relaxed_p05_in_any_region = any(relaxed_pattern_support_p05, na.rm = TRUE),
    exact_regions = paste(region[exact_pattern_match], collapse = ";"),
    broad_regions = paste(region[broad_pattern_support], collapse = ";"),
    relaxed_regions = paste(region[relaxed_pattern_support], collapse = ";"),
    relaxed_p05_regions = paste(region[relaxed_pattern_support_p05], collapse = ";"),
    .groups = "drop"
  ) %>%
  group_by(overall_pattern) %>%
  summarise(
    n_taxa = n(),
    exact_any = sum(exact_in_any_region, na.rm = TRUE),
    exact_any_pct = exact_any / n_taxa,
    broad_any = sum(broad_in_any_region, na.rm = TRUE),
    broad_any_pct = broad_any / n_taxa,
    relaxed_any = sum(relaxed_in_any_region, na.rm = TRUE),
    relaxed_any_pct = relaxed_any / n_taxa,
    relaxed_p05_any = sum(relaxed_p05_in_any_region, na.rm = TRUE),
    relaxed_p05_any_pct = relaxed_p05_any / n_taxa,
    .groups = "drop"
  ) %>%
  arrange(
    factor(overall_pattern, levels = c("Nw-specific", "Ow-specific", "Ob-specific", "Nw+Ow shared", "Nw+Ob shared", "Ow+Ob shared", "BMI-independent"))
  )

safe_write_csv(pattern_consistency, file.path(out_dir, "Regional_Pattern_Consistency_With_Overall_Long.csv"))
safe_write_csv(summary_by_pattern, file.path(out_dir, "Regional_Pattern_Consistency_By_Overall_Pattern.csv"))
safe_write_csv(summary_any_region, file.path(out_dir, "Regional_Pattern_Consistency_Any_Region_Summary.csv"))

cat("Pattern consistency by overall pattern and region:\n")
print(summary_by_pattern, n = 50)
cat("\nPattern consistency in at least one region:\n")
print(summary_any_region, n = 50)
