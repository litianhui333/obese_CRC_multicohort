suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(lmerTest)
})

analysis_dir <- "analyses/2026-05-28_alpha_diversity_confirmed"
fig_dir <- file.path(analysis_dir, "figures")
tab_dir <- file.path(analysis_dir, "tables")
data_fig_dir <- file.path(analysis_dir, "data_for_figures")
log_dir <- file.path(analysis_dir, "logs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

metadata_path <- "02.processed_data/02_aligned_data/metadata_aligned_20260115.rds"
species_path <- "02.processed_data/03_metaphlan_split_data/metaphlan_species_renormalized_20260115.rds"
oral_list_path <- "05.docs/oral_microbiome.csv"

alpha_indices <- c("Shannon", "Simpson", "Observed")
oral_indices <- c("log_Oral_Score", "Oral_Richness")
bmi_vars <- c(
  "BMI_Group1",
  "BMI_Group3_Ow",
  "BMI_Group1_China_as_OtherAsia",
  "BMI_Group3_Ow_China_as_OtherAsia"
)

metric_label <- function(x) {
  recode(
    x,
    "Observed" = "Richness",
    "log_Oral_Score" = "Oral Score",
    "Oral_Richness" = "Oral Richness",
    .default = x
  )
}

metric_levels <- c("Shannon", "Richness", "Simpson", "Oral Score", "Oral Richness")

display_cohort <- function(x) {
  out <- x
  out <- gsub("^Public_study__", "", out)
  out <- gsub("^This_study_cohort[0-9]+__", "", out)
  recode(
    out,
    "2017_YuJ" = "YuJ_2017",
    "2021_YuJ" = "inhouse-HK_1",
    "CHN_WF-CRC" = "SunYY_2025",
    "2023_Yunnan" = "inhouse-YN",
    "2025_HK" = "inhouse-HK_2",
    .default = out
  )
}

comparison_direction <- function(g1, g2) {
  priority <- c(
    "CTRL" = 1, "Control" = 1, "Healthy" = 1,
    "Normal" = 1, "Normal weight" = 1, "NW" = 1,
    "AA" = 2, "Adenoma" = 2,
    "Ow" = 2, "Overweight" = 2,
    "CRC" = 3, "Cancer" = 3,
    "Ob" = 3, "Obese" = 3,
    "Overweight/Obese" = 3
  )
  p1 <- ifelse(g1 %in% names(priority), unname(priority[g1]), 99)
  p2 <- ifelse(g2 %in% names(priority), unname(priority[g2]), 99)
  if (p1 <= p2) {
    list(ref = g1, target = g2)
  } else {
    list(ref = g2, target = g1)
  }
}

extract_group_coef <- function(fit, target) {
  coefs <- summary(fit)$coefficients
  hit <- grep(target, rownames(coefs), fixed = TRUE)
  if (length(hit) == 0) return(NULL)
  row <- coefs[hit[1], , drop = FALSE]
  p_col <- grep("^Pr\\(", colnames(row), value = TRUE)
  p_val <- if (length(p_col) > 0) row[1, p_col[1]] else NA_real_
  tibble(
    Beta = unname(row[1, "Estimate"]),
    SE = unname(row[1, "Std. Error"]),
    P_Value = unname(p_val)
  )
}

fit_pair <- function(data, alpha_index, group_col, ref, target, adjusted, random_study = TRUE) {
  needed <- c(alpha_index, group_col, "Study.name")
  if (adjusted) needed <- c(needed, "Age", "Sex")
  clean <- data %>%
    filter(if_all(all_of(needed), ~ !is.na(.x))) %>%
    filter(.data[[group_col]] %in% c(ref, target)) %>%
    mutate(
      Group_Factor = factor(.data[[group_col]], levels = c(ref, target)),
      Sex = as.factor(Sex)
    )
  if (nrow(clean) < 10 || length(unique(clean$Group_Factor)) < 2) return(NULL)
  n_ref <- sum(clean$Group_Factor == ref)
  n_target <- sum(clean$Group_Factor == target)
  if (n_ref < 3 || n_target < 3) return(NULL)
  n_studies <- length(unique(clean$Study.name))
  rhs <- if (adjusted) "Group_Factor + Age + Sex" else "Group_Factor"
  formula <- if (random_study && n_studies >= 2) {
    as.formula(paste(alpha_index, "~", rhs, "+ (1|Study.name)"))
  } else {
    as.formula(paste(alpha_index, "~", rhs))
  }
  fit <- tryCatch(
    if (random_study && n_studies >= 2) lmer(formula, data = clean) else lm(formula, data = clean),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  res <- extract_group_coef(fit, target)
  if (is.null(res)) return(NULL)
  res %>%
    mutate(
      Method = ifelse(random_study && n_studies >= 2, "LMM", "LM"),
      Adjustment = ifelse(adjusted, "Adjusted_Age_Sex", "Unadjusted"),
      Ref_Group = ref,
      Target_Group = target,
      n_ref = n_ref,
      n_target = n_target,
      n_samples = nrow(clean),
      n_studies = n_studies
    )
}

run_pairwise_set <- function(data, alpha_index, bmi_col, subgroup_var, subgroup_val) {
  out <- list()
  bmi_groups <- unique(na.omit(as.character(data[[bmi_col]])))
  diseases <- intersect(c("CTRL", "AA", "CRC"), unique(as.character(data$Disease)))

  if (length(bmi_groups) >= 2) {
    for (dz in diseases) {
      dz_data <- data %>% filter(Disease == dz)
      for (comb in combn(bmi_groups, 2, simplify = FALSE)) {
        dir <- comparison_direction(comb[1], comb[2])
        for (adj in c(FALSE, TRUE)) {
          res <- fit_pair(dz_data, alpha_index, bmi_col, dir$ref, dir$target, adj, TRUE)
          if (!is.null(res)) {
            out[[length(out) + 1]] <- res %>%
              mutate(
                Type = "BMI_Comparison",
                Context = dz,
                Comparison = paste(dir$ref, "vs", dir$target),
                BMI_Group_Column = bmi_col,
                Alpha_Index = alpha_index,
                Subgroup_Var = subgroup_var,
                Subgroup_Val = subgroup_val
              )
          }
        }
      }
    }
  }

  disease_groups <- unique(na.omit(as.character(data$Disease)))
  if (length(disease_groups) >= 2) {
    for (bmi in bmi_groups) {
      bmi_data <- data %>% filter(.data[[bmi_col]] == bmi)
      for (comb in combn(disease_groups, 2, simplify = FALSE)) {
        dir <- comparison_direction(comb[1], comb[2])
        for (adj in c(FALSE, TRUE)) {
          res <- fit_pair(bmi_data, alpha_index, "Disease", dir$ref, dir$target, adj, TRUE)
          if (!is.null(res)) {
            out[[length(out) + 1]] <- res %>%
              mutate(
                Type = "Disease_Comparison",
                Context = bmi,
                Comparison = paste(dir$ref, "vs", dir$target),
                BMI_Group_Column = bmi_col,
                Alpha_Index = alpha_index,
                Subgroup_Var = subgroup_var,
                Subgroup_Val = subgroup_val
              )
          }
        }
      }
    }
  }
  bind_rows(out)
}

run_cohort_pairs <- function(data, alpha_index, bmi_col, subgroup_var, subgroup_val) {
  out <- list()
  bmi_groups <- unique(na.omit(as.character(data[[bmi_col]])))
  diseases <- intersect(c("CTRL", "AA", "CRC"), unique(as.character(data$Disease)))

  for (study in sort(unique(data$Study.name))) {
    study_data <- data %>% filter(Study.name == study)
    if (length(bmi_groups) >= 2) {
      for (dz in diseases) {
        dz_data <- study_data %>% filter(Disease == dz)
        for (comb in combn(bmi_groups, 2, simplify = FALSE)) {
          dir <- comparison_direction(comb[1], comb[2])
          for (adj in c(FALSE, TRUE)) {
            res <- fit_pair(dz_data, alpha_index, bmi_col, dir$ref, dir$target, adj, FALSE)
            if (!is.null(res)) {
              out[[length(out) + 1]] <- res %>%
                mutate(
                  Type = "BMI_Comparison",
                  Context = dz,
                  Comparison = paste(dir$ref, "vs", dir$target),
                  BMI_Group_Column = bmi_col,
                  Alpha_Index = alpha_index,
                  Subgroup_Var = subgroup_var,
                  Subgroup_Val = subgroup_val,
                  Study.name = study,
                  Cohort_Display = display_cohort(study)
                )
            }
          }
        }
      }
    }
    disease_groups <- unique(na.omit(as.character(study_data$Disease)))
    if (length(disease_groups) >= 2) {
      for (bmi in bmi_groups) {
        bmi_data <- study_data %>% filter(.data[[bmi_col]] == bmi)
        for (comb in combn(disease_groups, 2, simplify = FALSE)) {
          dir <- comparison_direction(comb[1], comb[2])
          for (adj in c(FALSE, TRUE)) {
            res <- fit_pair(bmi_data, alpha_index, "Disease", dir$ref, dir$target, adj, FALSE)
            if (!is.null(res)) {
              out[[length(out) + 1]] <- res %>%
                mutate(
                  Type = "Disease_Comparison",
                  Context = bmi,
                  Comparison = paste(dir$ref, "vs", dir$target),
                  BMI_Group_Column = bmi_col,
                  Alpha_Index = alpha_index,
                  Subgroup_Var = subgroup_var,
                  Subgroup_Val = subgroup_val,
                  Study.name = study,
                  Cohort_Display = display_cohort(study)
                )
            }
          }
        }
      }
    }
  }
  bind_rows(out)
}

add_fdr <- function(df) {
  df %>%
    group_by(Method, Adjustment, Alpha_Index, Type, Subgroup_Val, BMI_Group_Column, Context) %>%
    mutate(
      P_FDR = p.adjust(P_Value, method = "fdr"),
      Significance = case_when(
        P_FDR < 0.001 ~ "***",
        P_FDR < 0.01 ~ "**",
        P_FDR < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) %>%
    ungroup() %>%
    mutate(
      CI_Lower = Beta - 1.96 * SE,
      CI_Upper = Beta + 1.96 * SE,
      Alpha_Label = metric_label(Alpha_Index),
      Contrast_Direction = paste0(Target_Group, " - ", Ref_Group)
    )
}

message("Reading input data")
metadata <- readRDS(metadata_path)
species_mat <- readRDS(species_path)
common_ids <- intersect(metadata$Sample.Source.ID.processed, colnames(species_mat))
metadata <- metadata %>% filter(Sample.Source.ID.processed %in% common_ids)
species_mat <- species_mat[, metadata$Sample.Source.ID.processed, drop = FALSE]

message("Calculating alpha-diversity metrics")
metadata <- metadata %>%
  mutate(
    Age = as.numeric(Age),
    Shannon = vegan::diversity(t(species_mat), index = "shannon"),
    Simpson = vegan::diversity(t(species_mat), index = "simpson"),
    Observed = colSums(species_mat > 0),
    Cohort_Display = display_cohort(Study.name)
  )

message("Calculating oral-to-gut metrics")
oral_microbiome <- read_csv(oral_list_path, show_col_types = FALSE)
oral_features <- oral_microbiome %>%
  transmute(
    Oral_Feature = as.character(feature),
    Oral_Species_ID = sub("\\|.*$", "", Oral_Feature),
    Oral_Species_ID = ifelse(grepl("^s__", Oral_Species_ID), Oral_Species_ID, paste0("s__", Oral_Species_ID))
  ) %>%
  filter(!is.na(Oral_Feature), Oral_Feature != "", !is.na(Oral_Species_ID), Oral_Species_ID != "") %>%
  distinct()

species_feature_map <- tibble(
  Species_Feature = rownames(species_mat),
  Species_ID = sub("^.*\\|s__", "s__", Species_Feature)
)

oral_feature_match <- oral_features %>%
  left_join(species_feature_map, by = c("Oral_Species_ID" = "Species_ID")) %>%
  mutate(Matched = !is.na(Species_Feature))

valid_oral_bugs <- unique(na.omit(oral_feature_match$Species_Feature))
if (length(valid_oral_bugs) == 0) {
  stop("No oral microbiome features matched the species table.")
}

oral_sub_mat <- species_mat[valid_oral_bugs, , drop = FALSE]
metadata <- metadata %>%
  mutate(
    Oral_Score = colSums(oral_sub_mat, na.rm = TRUE),
    Oral_Richness = colSums(oral_sub_mat > 0, na.rm = TRUE),
    log_Oral_Score = log10(Oral_Score + 1e-6)
  )

write_tsv(oral_feature_match, file.path(tab_dir, "oral_to_gut_feature_match.tsv"))

saveRDS(metadata, file.path(tab_dir, "metadata_with_alpha_confirmed.rds"))
write_tsv(
  metadata %>%
    select(Sample.Source.ID.processed, Study.name, Cohort_Display, Disease, all_of(intersect(bmi_vars, names(metadata))), Continent, Age, Sex, Shannon, Simpson, Observed, Oral_Score, log_Oral_Score, Oral_Richness),
  file.path(tab_dir, "metadata_with_alpha_confirmed.tsv")
)

write_tsv(
  tibble(
    Study.name = sort(unique(metadata$Study.name)),
    Cohort_Display = display_cohort(sort(unique(metadata$Study.name)))
  ),
  file.path(tab_dir, "cohort_display_label_map.tsv")
)

message("Running pooled LMM/LM pairwise analyses")
pooled_results <- list()
for (idx in alpha_indices) {
  for (bmi_col in intersect(bmi_vars, names(metadata))) {
    pooled_results[[length(pooled_results) + 1]] <- run_pairwise_set(metadata, idx, bmi_col, "Dataset", "Overall")
    for (cont in c("Asia", "Europe")) {
      sub_meta <- metadata %>% filter(Continent == cont)
      pooled_results[[length(pooled_results) + 1]] <- run_pairwise_set(sub_meta, idx, bmi_col, "Continent", cont)
    }
  }
}
pooled_results <- bind_rows(pooled_results) %>% add_fdr()
write_tsv(pooled_results, file.path(tab_dir, "alpha_pairwise_pooled_LMM_all.tsv"))
write_csv(pooled_results, file.path(tab_dir, "alpha_pairwise_pooled_LMM_all.csv"))

message("Running cohort-specific LM pairwise analyses")
cohort_results <- list()
for (idx in alpha_indices) {
  for (bmi_col in intersect(bmi_vars, names(metadata))) {
    cohort_results[[length(cohort_results) + 1]] <- run_cohort_pairs(metadata, idx, bmi_col, "Dataset", "Overall")
    for (cont in c("Asia", "Europe")) {
      sub_meta <- metadata %>% filter(Continent == cont)
      cohort_results[[length(cohort_results) + 1]] <- run_cohort_pairs(sub_meta, idx, bmi_col, "Continent", cont)
    }
  }
}
cohort_results <- bind_rows(cohort_results) %>%
  mutate(
    CI_Lower = Beta - 1.96 * SE,
    CI_Upper = Beta + 1.96 * SE,
    Alpha_Label = metric_label(Alpha_Index),
    Contrast_Direction = paste0(Target_Group, " - ", Ref_Group)
  )
write_tsv(cohort_results, file.path(tab_dir, "alpha_pairwise_cohort_specific_LM_all.tsv"))

message("Running oral-to-gut pooled LMM/LM pairwise analyses")
oral_pooled_results <- list()
for (idx in oral_indices) {
  for (bmi_col in intersect(bmi_vars, names(metadata))) {
    oral_pooled_results[[length(oral_pooled_results) + 1]] <- run_pairwise_set(metadata, idx, bmi_col, "Dataset", "Overall")
    for (cont in c("Asia", "Europe")) {
      sub_meta <- metadata %>% filter(Continent == cont)
      oral_pooled_results[[length(oral_pooled_results) + 1]] <- run_pairwise_set(sub_meta, idx, bmi_col, "Continent", cont)
    }
  }
}
oral_pooled_results <- bind_rows(oral_pooled_results) %>% add_fdr()
write_tsv(oral_pooled_results, file.path(tab_dir, "oral_to_gut_pairwise_pooled_LMM_all.tsv"))
write_csv(oral_pooled_results, file.path(tab_dir, "oral_to_gut_pairwise_pooled_LMM_all.csv"))

message("Running oral-to-gut cohort-specific LM pairwise analyses")
oral_cohort_results <- list()
for (idx in oral_indices) {
  for (bmi_col in intersect(bmi_vars, names(metadata))) {
    oral_cohort_results[[length(oral_cohort_results) + 1]] <- run_cohort_pairs(metadata, idx, bmi_col, "Dataset", "Overall")
    for (cont in c("Asia", "Europe")) {
      sub_meta <- metadata %>% filter(Continent == cont)
      oral_cohort_results[[length(oral_cohort_results) + 1]] <- run_cohort_pairs(sub_meta, idx, bmi_col, "Continent", cont)
    }
  }
}
oral_cohort_results <- bind_rows(oral_cohort_results) %>%
  mutate(
    CI_Lower = Beta - 1.96 * SE,
    CI_Upper = Beta + 1.96 * SE,
    Alpha_Label = metric_label(Alpha_Index),
    Contrast_Direction = paste0(Target_Group, " - ", Ref_Group)
  )
write_tsv(oral_cohort_results, file.path(tab_dir, "oral_to_gut_pairwise_cohort_specific_LM_all.tsv"))

oral_concordance <- oral_pooled_results %>%
  filter(Adjustment == "Adjusted_Age_Sex") %>%
  select(
    Type, Context, Comparison, BMI_Group_Column, Alpha_Index, Alpha_Label,
    Subgroup_Var, Subgroup_Val, Pooled_Beta = Beta, Pooled_P_FDR = P_FDR,
    Pooled_Significance = Significance
  ) %>%
  left_join(
    oral_cohort_results %>%
      filter(Adjustment == "Adjusted_Age_Sex") %>%
      select(Type, Context, Comparison, BMI_Group_Column, Alpha_Index, Subgroup_Var, Subgroup_Val, Study.name, Cohort_Display, Beta),
    by = c("Type", "Context", "Comparison", "BMI_Group_Column", "Alpha_Index", "Subgroup_Var", "Subgroup_Val")
  ) %>%
  filter(!is.na(Beta), Pooled_Beta != 0, Beta != 0) %>%
  group_by(Type, Context, Comparison, BMI_Group_Column, Alpha_Index, Alpha_Label, Subgroup_Var, Subgroup_Val, Pooled_Beta, Pooled_P_FDR, Pooled_Significance) %>%
  summarise(
    n_same_direction = sum(sign(Beta) == sign(Pooled_Beta)),
    n_total_direction = n(),
    direction_concordance = n_same_direction / n_total_direction,
    concordance_label = paste0(n_same_direction, "/", n_total_direction),
    concordance_class = ifelse(direction_concordance >= 0.7, ">=70% same direction", "<70% same direction"),
    .groups = "drop"
  )
write_tsv(oral_concordance, file.path(tab_dir, "oral_to_gut_direction_concordance_adjusted.tsv"))

oral_supp_table <- oral_pooled_results %>%
  left_join(
    oral_concordance,
    by = c("Type", "Context", "Comparison", "BMI_Group_Column", "Alpha_Index", "Alpha_Label", "Subgroup_Var", "Subgroup_Val"),
    suffix = c("", "_Concordance")
  ) %>%
  arrange(Alpha_Index, BMI_Group_Column, Type, Subgroup_Val, Context, Comparison, Adjustment) %>%
  select(
    Oral_Metric = Alpha_Label, Alpha_Index, Type, Subgroup_Var, Subgroup_Val,
    BMI_Group_Column, Context, Comparison, Contrast_Direction,
    Method, Adjustment, Beta, SE, CI_Lower, CI_Upper, P_Value, P_FDR, Significance,
    n_ref, n_target, n_samples, n_studies,
    n_same_direction, n_total_direction, direction_concordance, concordance_label, concordance_class
  )
write_tsv(oral_supp_table, file.path(tab_dir, "Supplementary_Table_oral_to_gut_pairwise_statistics.tsv"))
write_csv(oral_supp_table, file.path(tab_dir, "Supplementary_Table_oral_to_gut_pairwise_statistics.csv"))

message("Calculating direction concordance")
concordance <- pooled_results %>%
  filter(Adjustment == "Adjusted_Age_Sex") %>%
  select(
    Type, Context, Comparison, BMI_Group_Column, Alpha_Index, Alpha_Label,
    Subgroup_Var, Subgroup_Val, Pooled_Beta = Beta, Pooled_P_FDR = P_FDR,
    Pooled_Significance = Significance
  ) %>%
  left_join(
    cohort_results %>%
      filter(Adjustment == "Adjusted_Age_Sex") %>%
      select(Type, Context, Comparison, BMI_Group_Column, Alpha_Index, Subgroup_Var, Subgroup_Val, Study.name, Cohort_Display, Beta),
    by = c("Type", "Context", "Comparison", "BMI_Group_Column", "Alpha_Index", "Subgroup_Var", "Subgroup_Val")
  ) %>%
  filter(!is.na(Beta), Pooled_Beta != 0, Beta != 0) %>%
  group_by(Type, Context, Comparison, BMI_Group_Column, Alpha_Index, Alpha_Label, Subgroup_Var, Subgroup_Val, Pooled_Beta, Pooled_P_FDR, Pooled_Significance) %>%
  summarise(
    n_same_direction = sum(sign(Beta) == sign(Pooled_Beta)),
    n_total_direction = n(),
    direction_concordance = n_same_direction / n_total_direction,
    concordance_label = paste0(n_same_direction, "/", n_total_direction),
    concordance_class = ifelse(direction_concordance >= 0.7, ">=70% same direction", "<70% same direction"),
    .groups = "drop"
  )
write_tsv(concordance, file.path(tab_dir, "alpha_direction_concordance_adjusted.tsv"))

supp_table <- pooled_results %>%
  left_join(
    concordance,
    by = c("Type", "Context", "Comparison", "BMI_Group_Column", "Alpha_Index", "Alpha_Label", "Subgroup_Var", "Subgroup_Val"),
    suffix = c("", "_Concordance")
  ) %>%
  arrange(Alpha_Index, BMI_Group_Column, Type, Subgroup_Val, Context, Comparison, Adjustment) %>%
  select(
    Alpha_Metric = Alpha_Label, Alpha_Index, Type, Subgroup_Var, Subgroup_Val,
    BMI_Group_Column, Context, Comparison, Contrast_Direction,
    Method, Adjustment, Beta, SE, CI_Lower, CI_Upper, P_Value, P_FDR, Significance,
    n_ref, n_target, n_samples, n_studies,
    n_same_direction, n_total_direction, direction_concordance, concordance_label, concordance_class
  )
write_tsv(supp_table, file.path(tab_dir, "Supplementary_Table_alpha_diversity_pairwise_statistics.tsv"))
write_csv(supp_table, file.path(tab_dir, "Supplementary_Table_alpha_diversity_pairwise_statistics.csv"))

theme_alpha <- function(base_size = 12) {
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
colors_bmi <- c(
  "Normal" = "#5AAE61",
  "Ow" = "#FF8C00",
  "Ob" = "#B07AA1",
  "Obese" = "#B07AA1",
  "Overweight/Obese" = "#B07AA1"
)

star_y <- function(df, alpha_col) {
  df %>%
    group_by(Alpha_Label) %>%
    summarise(y = max(.data[[alpha_col]], na.rm = TRUE) * 1.07, .groups = "drop")
}

format_sig_p <- function(p) {
  ifelse(p < 0.01, "FDR<0.01", paste0("FDR=", sprintf("%.2f", p)))
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
    mutate(y_range = ifelse(y_range == 0, y_max * 0.1, y_range))

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
    mutate(y_range = ifelse(y_range == 0, y_max * 0.1, y_range))

  bmi_x <- c("Normal" = 1, "Ow" = 2, "Ob" = 3, "Overweight/Obese" = 2)
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

make_boxplot <- function(subgroup_val, bmi_col, mode, metrics, filename, width, height,
                         stats_results = pooled_results, title_suffix = "alpha diversity") {
  sub <- if (subgroup_val == "Overall") metadata else metadata %>% filter(Continent == subgroup_val)
  plot_df <- sub %>%
    select(Sample.Source.ID.processed, Study.name, Cohort_Display, Continent, Disease, all_of(bmi_col), all_of(metrics)) %>%
    pivot_longer(all_of(metrics), names_to = "Alpha_Index", values_to = "Alpha_Value") %>%
    mutate(
      Alpha_Label = factor(metric_label(Alpha_Index), levels = metric_levels),
      Disease = factor(Disease, levels = c("CTRL", "AA", "CRC")),
      BMI_Display = factor(.data[[bmi_col]], levels = c("Normal", "Ow", "Ob", "Overweight/Obese"))
    ) %>%
    filter(!is.na(Disease), !is.na(BMI_Display), !is.na(Alpha_Value))

  write_tsv(plot_df, file.path(data_fig_dir, paste0(filename, ".tsv")))

  if (mode == "xBMI") {
    stat_df <- stats_results %>%
      filter(
        Adjustment == "Adjusted_Age_Sex",
        Subgroup_Val == subgroup_val,
        BMI_Group_Column == bmi_col,
        Type == "BMI_Comparison",
        Alpha_Index %in% metrics
      ) %>%
      mutate(
        Alpha_Label = factor(metric_label(Alpha_Index), levels = metric_levels),
        label = ifelse(Significance == "ns", "", format_sig_p(P_FDR))
      ) %>%
      filter(label != "") %>%
      make_brackets_xbmi(plot_df)

    p <- ggplot(plot_df, aes(x = Disease, y = Alpha_Value, fill = BMI_Display)) +
      geom_violin(position = position_dodge(0.8), alpha = 0.3, color = NA, trim = FALSE) +
      geom_boxplot(position = position_dodge(0.8), width = 0.2,
                   outlier.size = 0.3, outlier.alpha = 0.6,
                   alpha = 0.8, color = "black", linewidth = 0.2) +
      geom_segment(data = stat_df, aes(x = x, xend = xend, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = 0.3) +
      geom_segment(data = stat_df, aes(x = x, xend = x, y = y, yend = ytick),
                   inherit.aes = FALSE, linewidth = 0.3) +
      geom_segment(data = stat_df, aes(x = xend, xend = xend, y = y, yend = ytick),
                   inherit.aes = FALSE, linewidth = 0.3) +
      geom_text(
        data = stat_df,
        aes(x = xmid, y = y, label = label),
        inherit.aes = FALSE,
        size = 2.25,
        vjust = -0.15
      ) +
      facet_wrap(~Alpha_Label, scales = "free_y", nrow = 1) +
      scale_fill_manual(values = colors_bmi, drop = TRUE, name = "BMI") +
      labs(
        title = paste0(subgroup_val, ": ", title_suffix),
        subtitle = "x-axis: Disease",
        x = NULL, y = NULL, fill = "BMI Group"
      ) +
      theme_alpha(12)
  } else {
    stat_df <- stats_results %>%
      filter(
        Adjustment == "Adjusted_Age_Sex",
        Subgroup_Val == subgroup_val,
        BMI_Group_Column == bmi_col,
        Type == "Disease_Comparison",
        Alpha_Index %in% metrics
      ) %>%
      mutate(
        Alpha_Label = factor(metric_label(Alpha_Index), levels = metric_levels),
        label = ifelse(Significance == "ns", "", format_sig_p(P_FDR))
      ) %>%
      filter(label != "") %>%
      make_brackets_xdisease(plot_df)

    p <- ggplot(plot_df, aes(x = BMI_Display, y = Alpha_Value, fill = Disease)) +
      geom_violin(position = position_dodge(0.8), alpha = 0.3, color = NA, trim = FALSE) +
      geom_boxplot(position = position_dodge(0.8), width = 0.2,
                   outlier.size = 0.3, outlier.alpha = 0.6,
                   alpha = 0.8, color = "black", linewidth = 0.2) +
      geom_segment(data = stat_df, aes(x = x, xend = xend, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = 0.3) +
      geom_segment(data = stat_df, aes(x = x, xend = x, y = y, yend = ytick),
                   inherit.aes = FALSE, linewidth = 0.3) +
      geom_segment(data = stat_df, aes(x = xend, xend = xend, y = y, yend = ytick),
                   inherit.aes = FALSE, linewidth = 0.3) +
      geom_text(
        data = stat_df,
        aes(x = xmid, y = y, label = label),
        inherit.aes = FALSE,
        size = 2.25,
        vjust = -0.15
      ) +
      facet_wrap(~Alpha_Label, scales = "free_y", nrow = 1) +
      scale_fill_manual(values = colors_disease, drop = TRUE, name = "Disease") +
      labs(
        title = paste0(subgroup_val, ": ", title_suffix),
        subtitle = "x-axis: BMI",
        x = NULL, y = NULL, fill = "Disease"
      ) +
      theme_alpha(12)
  }

  ggsave(file.path(fig_dir, paste0(filename, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(filename, ".png")), p, width = width, height = height, dpi = 300)
}

make_forest <- function(subgroup_vals, bmi_col, filename, width, height) {
  target_metrics <- c("Shannon", "Observed")
  target_comparisons <- if (bmi_col == "BMI_Group1") {
    c("Normal vs Ow", "Normal vs Ob", "Ow vs Ob")
  } else {
    "Normal vs Overweight/Obese"
  }

  pooled_df <- pooled_results %>%
    filter(
      Adjustment == "Adjusted_Age_Sex",
      BMI_Group_Column == bmi_col,
      Type == "BMI_Comparison",
      Comparison %in% target_comparisons,
      Alpha_Index %in% target_metrics,
      Subgroup_Val %in% subgroup_vals
    ) %>%
    mutate(
      Layer = "Pooled LMM",
      Cohort_Display = "Pooled LMM",
      Alpha_Label = metric_label(Alpha_Index),
      Facet_Label = paste(Alpha_Label, Context, Comparison, Subgroup_Val, sep = " | ")
    )

  cohort_df <- cohort_results %>%
    filter(
      Adjustment == "Adjusted_Age_Sex",
      BMI_Group_Column == bmi_col,
      Type == "BMI_Comparison",
      Comparison %in% target_comparisons,
      Alpha_Index %in% target_metrics,
      Subgroup_Val %in% subgroup_vals
    ) %>%
    mutate(
      Layer = "Cohort LM",
      Alpha_Label = metric_label(Alpha_Index),
      Facet_Label = paste(Alpha_Label, Context, Comparison, Subgroup_Val, sep = " | ")
    )

  plot_df <- bind_rows(
    cohort_df %>% select(Facet_Label, Subgroup_Val, Context, Alpha_Label, Cohort_Display, Layer, Beta, CI_Lower, CI_Upper, P_Value),
    pooled_df %>% select(Facet_Label, Subgroup_Val, Context, Alpha_Label, Cohort_Display, Layer, Beta, CI_Lower, CI_Upper, P_Value, P_FDR, Significance)
  ) %>%
    mutate(
      Context = factor(Context, levels = c("CTRL", "AA", "CRC")),
      Alpha_Label = factor(Alpha_Label, levels = metric_levels),
      Layer = factor(Layer, levels = c("Cohort LM", "Pooled LMM"))
    )

  write_tsv(plot_df, file.path(data_fig_dir, paste0(filename, ".tsv")))

  p <- ggplot(plot_df, aes(x = Beta, y = Cohort_Display)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    geom_errorbar(
      data = plot_df %>% filter(Layer == "Pooled LMM"),
      aes(xmin = CI_Lower, xmax = CI_Upper),
      orientation = "y",
      width = 0.18,
      linewidth = 0.5,
      color = "black"
    ) +
    geom_point(aes(shape = Layer, fill = Layer, color = Layer), size = 2.1, stroke = 0.35) +
    facet_wrap(~Facet_Label, scales = "free", ncol = 2) +
    scale_shape_manual(values = c("Cohort LM" = 21, "Pooled LMM" = 23), name = NULL) +
    scale_fill_manual(values = c("Cohort LM" = "grey78", "Pooled LMM" = "#E45756"), name = NULL) +
    scale_color_manual(values = c("Cohort LM" = "grey55", "Pooled LMM" = "#B22222"), name = NULL) +
    labs(
      x = "Adjusted beta (Target - Ref)",
      y = NULL,
      title = "Cohort-specific alpha-diversity estimates"
    ) +
    theme_alpha(9) +
    theme(
      legend.position = "top",
      axis.text.y = element_text(size = 6.5),
      strip.text = element_text(size = 7.5)
    )

  ggsave(file.path(fig_dir, paste0(filename, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(filename, ".png")), p, width = width, height = height, dpi = 300)
}

make_significant_forest <- function(bmi_col, filename, width, height) {
  sig_keys <- pooled_results %>%
    filter(
      Adjustment == "Adjusted_Age_Sex",
      Subgroup_Val == "Overall",
      BMI_Group_Column == bmi_col,
      Type == "BMI_Comparison",
      Alpha_Index %in% c("Shannon", "Observed"),
      P_FDR < 0.05
    ) %>%
    mutate(
      Alpha_Label = metric_label(Alpha_Index),
      Facet_Label = paste(Alpha_Label, Context, Comparison, sep = " | ")
    )

  if (nrow(sig_keys) == 0) return(invisible(NULL))

  cohort_df <- cohort_results %>%
    filter(
      Adjustment == "Adjusted_Age_Sex",
      Subgroup_Val == "Overall",
      BMI_Group_Column == bmi_col,
      Type == "BMI_Comparison",
      Alpha_Index %in% unique(sig_keys$Alpha_Index)
    ) %>%
    inner_join(
      sig_keys %>% select(Alpha_Index, Context, Comparison, Facet_Label),
      by = c("Alpha_Index", "Context", "Comparison")
    ) %>%
    mutate(Layer = "Cohort LM")

  pooled_df <- sig_keys %>%
    mutate(
      Layer = "Pooled LMM",
      Cohort_Display = "Pooled LMM"
    )

  plot_df <- bind_rows(
    cohort_df %>% select(Facet_Label, Alpha_Label, Context, Comparison, Cohort_Display, Layer, Beta, CI_Lower, CI_Upper, P_Value),
    pooled_df %>% select(Facet_Label, Alpha_Label, Context, Comparison, Cohort_Display, Layer, Beta, CI_Lower, CI_Upper, P_Value, P_FDR, Significance)
  ) %>%
    mutate(
      Layer = factor(Layer, levels = c("Cohort LM", "Pooled LMM")),
      Facet_Label = factor(Facet_Label, levels = unique(sig_keys$Facet_Label))
    )

  write_tsv(plot_df, file.path(data_fig_dir, paste0(filename, ".tsv")))

  p <- ggplot(plot_df, aes(x = Beta, y = Cohort_Display)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    geom_errorbar(
      data = plot_df %>% filter(Layer == "Pooled LMM"),
      aes(xmin = CI_Lower, xmax = CI_Upper),
      orientation = "y",
      width = 0.18,
      linewidth = 0.5,
      color = "black"
    ) +
    geom_point(aes(shape = Layer, fill = Layer, color = Layer), size = 2.1, stroke = 0.35) +
    facet_wrap(~Facet_Label, scales = "free", ncol = 2) +
    scale_shape_manual(values = c("Cohort LM" = 21, "Pooled LMM" = 23), name = NULL) +
    scale_fill_manual(values = c("Cohort LM" = "grey78", "Pooled LMM" = "#E74C3C"), name = NULL) +
    scale_color_manual(values = c("Cohort LM" = "grey55", "Pooled LMM" = "#B22222"), name = NULL) +
    labs(
      x = "Adjusted beta (Target - Ref)",
      y = NULL,
      title = "Cohort-specific estimates for significant pooled BMI comparisons"
    ) +
    theme_alpha(9) +
    theme(
      legend.position = "top",
      axis.text.y = element_text(size = 6.5),
      strip.text = element_text(size = 7.5)
    )

  ggsave(file.path(fig_dir, paste0(filename, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(filename, ".png")), p, width = width, height = height, dpi = 300)
}

make_legacy_style_bmi_forest <- function(subgroup_vals, bmi_col, filename, title_text, width, height,
                                         stats_results = pooled_results,
                                         cohort_stats_results = cohort_results,
                                         target_metrics = c("Shannon", "Observed"),
                                         target_metric_levels = c("Shannon", "Richness")) {
  dataset_map <- c("Overall" = "Overall", "Asia" = "Asia", "Europe" = "nonAsia")
  # ggplot places the first y-axis level at the bottom, so this level order
  # displays top-to-bottom as Nw vs Ob, Nw vs Ow, Ow vs Ob.
  target_comparisons <- c("Ow vs Ob", "Normal vs Ow", "Normal vs Ob")
  comparison_labels <- c(
    "Normal vs Ob" = "Nw vs Ob",
    "Normal vs Ow" = "Nw vs Ow",
    "Ow vs Ob" = "Ow vs Ob"
  )
  disease_levels <- c("CTRL", "AA", "CRC")

  cohort_df <- cohort_stats_results %>%
    filter(
      Adjustment == "Adjusted_Age_Sex",
      Subgroup_Val %in% subgroup_vals,
      BMI_Group_Column == bmi_col,
      Type == "BMI_Comparison",
      Alpha_Index %in% target_metrics,
      Comparison %in% target_comparisons,
      Context %in% disease_levels
    ) %>%
    mutate(
      Dataset = recode(Subgroup_Val, !!!dataset_map),
      Alpha_Label = factor(metric_label(Alpha_Index), levels = target_metric_levels),
      Context = factor(Context, levels = rev(disease_levels)),
      Comparison = factor(Comparison, levels = target_comparisons),
      Study_Label = Cohort_Display
    )

  pooled_df <- stats_results %>%
    filter(
      Adjustment == "Adjusted_Age_Sex",
      Subgroup_Val %in% subgroup_vals,
      BMI_Group_Column == bmi_col,
      Type == "BMI_Comparison",
      Alpha_Index %in% target_metrics,
      Comparison %in% target_comparisons,
      Context %in% disease_levels
    ) %>%
    mutate(
      Dataset = recode(Subgroup_Val, !!!dataset_map),
      Alpha_Label = factor(metric_label(Alpha_Index), levels = target_metric_levels),
      Context = factor(Context, levels = rev(disease_levels)),
      Comparison = factor(Comparison, levels = target_comparisons),
      Sig_Class = if_else(P_FDR < 0.05, "Significant", "Not significant")
    )

  if (nrow(cohort_df) == 0 || nrow(pooled_df) == 0) return(invisible(NULL))

  if (identical(subgroup_vals, "Overall")) {
    cohort_df <- cohort_df %>%
      mutate(Row_Label = factor(as.character(Context), levels = rev(disease_levels)))
    pooled_df <- pooled_df %>%
      mutate(Row_Label = factor(as.character(Context), levels = rev(disease_levels)))
  } else {
    row_levels <- c(
      "Asia | CRC", "Asia | AA", "Asia | CTRL",
      "nonAsia | CRC", "nonAsia | AA", "nonAsia | CTRL"
    )
    cohort_df <- cohort_df %>%
      mutate(Row_Label = factor(paste(Dataset, Context, sep = " | "), levels = row_levels))
    pooled_df <- pooled_df %>%
      mutate(Row_Label = factor(paste(Dataset, Context, sep = " | "), levels = row_levels))
  }

  study_levels <- sort(unique(cohort_df$Study.name))
  shape_pool <- c(0:18, 21:25)
  shape_values <- setNames(rep(shape_pool, length.out = length(study_levels)), study_levels)
  shape_labels <- setNames(display_cohort(study_levels), study_levels)

  plot_df <- bind_rows(
    cohort_df %>%
      mutate(Layer = "Cohort", CI_Lower = NA_real_, CI_Upper = NA_real_, P_FDR = NA_real_, Sig_Class = NA_character_),
    pooled_df %>%
      mutate(Layer = "Pooled", Study.name = NA_character_, Study_Label = NA_character_)
  )
  write_tsv(plot_df, file.path(data_fig_dir, paste0(filename, ".tsv")))

  p <- ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.4) +
    geom_point(
      data = cohort_df,
      aes(x = Beta, y = Comparison, shape = Study.name),
      color = "grey65",
      alpha = 0.55,
      size = 1.9,
      stroke = 0.5,
      position = position_jitter(width = 0, height = 0.10)
    ) +
    geom_errorbar(
      data = pooled_df,
      aes(y = Comparison, xmin = CI_Lower, xmax = CI_Upper, color = Sig_Class),
      orientation = "y",
      width = 0.16,
      linewidth = 0.75
    ) +
    geom_point(
      data = pooled_df,
      aes(x = Beta, y = Comparison, color = Sig_Class),
      size = 3.2
    ) +
    facet_grid(rows = vars(Row_Label), cols = vars(Alpha_Label), scales = "free_x", switch = "y") +
    scale_y_discrete(labels = comparison_labels) +
    scale_shape_manual(values = shape_values, labels = shape_labels, name = "Cohort") +
    scale_color_manual(
      values = c("Significant" = "#C62828", "Not significant" = "black"),
      name = "Pooled LMM"
    ) +
    labs(
      title = title_text,
      subtitle = "Grey points: cohort-level adjusted LM; colored points and CI: pooled adjusted LMM (red = FDR < 0.05)",
      x = "Effect size (beta)",
      y = NULL
    ) +
    theme_bw(base_size = 8.5) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.35, "lines"),
      panel.spacing.y = unit(0.15, "lines"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.title.y = element_blank(),
      axis.text = element_text(size = 8),
      axis.title.x = element_text(size = 9),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(face = "bold", size = 7.5),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.35, "cm"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 9.5),
      plot.subtitle = element_text(hjust = 0.5, size = 7.5),
      strip.placement = "outside"
    ) +
    guides(
      shape = guide_legend(ncol = 4, override.aes = list(alpha = 0.9, color = "grey35", size = 2.2)),
      color = guide_legend(nrow = 1)
    )

  ggsave(file.path(fig_dir, paste0(filename, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(filename, ".png")), p, width = width, height = height, dpi = 300)
}

message("Making figures and source data")
make_boxplot("Overall", "BMI_Group1", "xBMI", c("Shannon", "Observed"), "Alpha_Overall_confirmed_LMM_adjP_xBMI_BMI_Group1", 6.4, 3.2)
make_boxplot("Asia", "BMI_Group1", "xBMI", c("Shannon", "Observed"), "Alpha_Asia_confirmed_LMM_adjP_xBMI_BMI_Group1", 6.4, 3.2)
make_boxplot("Europe", "BMI_Group1", "xBMI", c("Shannon", "Observed"), "Alpha_Europe_confirmed_LMM_adjP_xBMI_BMI_Group1", 6.4, 3.2)
make_boxplot("Overall", "BMI_Group1", "xDisease", c("Shannon", "Observed"), "Alpha_Overall_confirmed_LMM_adjP_xDisease_BMI_Group1", 6.4, 3.2)
make_boxplot("Asia", "BMI_Group1", "xDisease", c("Shannon", "Observed"), "Alpha_Asia_confirmed_LMM_adjP_xDisease_BMI_Group1", 6.4, 3.2)
make_boxplot("Europe", "BMI_Group1", "xDisease", c("Shannon", "Observed"), "Alpha_Europe_confirmed_LMM_adjP_xDisease_BMI_Group1", 6.4, 3.2)
make_legacy_style_bmi_forest("Overall", "BMI_Group1", "Alpha_LMM_CohortOverlay_Overall_BMI_within_Disease_confirmed", "Overall: BMI comparisons within each disease", 3.8, 4.8)
make_legacy_style_bmi_forest(c("Asia", "Europe"), "BMI_Group1", "Alpha_LMM_CohortOverlay_Asia_nonAsia_BMI_within_Disease_confirmed", "Asia and nonAsia: BMI comparisons within each disease", 3.8, 6)

make_boxplot(
  "Overall", "BMI_Group1", "xBMI", oral_indices,
  "Oral_to_Gut_Overall_confirmed_LMM_adjP_xBMI_BMI_Group1", 6.4, 3.2,
  oral_pooled_results, "oral-to-gut indices"
)
make_boxplot(
  "Asia", "BMI_Group1", "xBMI", oral_indices,
  "Oral_to_Gut_Asia_confirmed_LMM_adjP_xBMI_BMI_Group1", 6.4, 3.2,
  oral_pooled_results, "oral-to-gut indices"
)
make_boxplot(
  "Europe", "BMI_Group1", "xBMI", oral_indices,
  "Oral_to_Gut_Europe_confirmed_LMM_adjP_xBMI_BMI_Group1", 6.4, 3.2,
  oral_pooled_results, "oral-to-gut indices"
)
make_boxplot(
  "Overall", "BMI_Group1", "xDisease", oral_indices,
  "Oral_to_Gut_Overall_confirmed_LMM_adjP_xDisease_BMI_Group1", 6.4, 3.2,
  oral_pooled_results, "oral-to-gut indices"
)
make_boxplot(
  "Asia", "BMI_Group1", "xDisease", oral_indices,
  "Oral_to_Gut_Asia_confirmed_LMM_adjP_xDisease_BMI_Group1", 6.4, 3.2,
  oral_pooled_results, "oral-to-gut indices"
)
make_boxplot(
  "Europe", "BMI_Group1", "xDisease", oral_indices,
  "Oral_to_Gut_Europe_confirmed_LMM_adjP_xDisease_BMI_Group1", 6.4, 3.2,
  oral_pooled_results, "oral-to-gut indices"
)
make_legacy_style_bmi_forest(
  "Overall", "BMI_Group1",
  "Oral_to_Gut_LMM_CohortOverlay_Overall_BMI_within_Disease_confirmed",
  "Overall: oral-to-gut BMI comparisons within each disease", 3.8, 4.8,
  oral_pooled_results, oral_cohort_results, oral_indices, c("Oral Score", "Oral Richness")
)
make_legacy_style_bmi_forest(
  c("Asia", "Europe"), "BMI_Group1",
  "Oral_to_Gut_LMM_CohortOverlay_Asia_nonAsia_BMI_within_Disease_confirmed",
  "Asia and nonAsia: oral-to-gut BMI comparisons within each disease", 3.8, 6,
  oral_pooled_results, oral_cohort_results, oral_indices, c("Oral Score", "Oral Richness")
)

summary_lines <- c(
  "Alpha-diversity confirmation analysis",
  paste("Input metadata:", metadata_path),
  paste("Input species table:", species_path),
  paste("Input oral feature list:", oral_list_path),
  paste("Samples:", nrow(metadata)),
  paste("Species features:", nrow(species_mat)),
  paste("Matched oral features:", length(valid_oral_bugs), "of", nrow(oral_features)),
  "",
  "Primary statistics use pairwise LMM/LM contrasts with beta defined as Target - Ref.",
  "Significance is based on FDR within method/adjustment/metric/type/subgroup/BMI/context families.",
  "Oral-to-gut metrics follow the same pairwise model and figure format as the alpha-diversity metrics.",
  "Cohort display labels were harmonized only for figures and tables; raw sample IDs and feature IDs were not altered."
)
writeLines(summary_lines, file.path(tab_dir, "analysis_summary.txt"))

writeLines(c(metadata_path, species_path, oral_list_path), file.path(analysis_dir, "input_files.txt"))
writeLines(
  c(
    file.path(tab_dir, "metadata_with_alpha_confirmed.tsv"),
    file.path(tab_dir, "cohort_display_label_map.tsv"),
    file.path(tab_dir, "oral_to_gut_feature_match.tsv"),
    file.path(tab_dir, "alpha_pairwise_pooled_LMM_all.tsv"),
    file.path(tab_dir, "alpha_pairwise_cohort_specific_LM_all.tsv"),
    file.path(tab_dir, "alpha_direction_concordance_adjusted.tsv"),
    file.path(tab_dir, "Supplementary_Table_alpha_diversity_pairwise_statistics.tsv"),
    file.path(tab_dir, "oral_to_gut_pairwise_pooled_LMM_all.tsv"),
    file.path(tab_dir, "oral_to_gut_pairwise_cohort_specific_LM_all.tsv"),
    file.path(tab_dir, "oral_to_gut_direction_concordance_adjusted.tsv"),
    file.path(tab_dir, "Supplementary_Table_oral_to_gut_pairwise_statistics.tsv"),
    fig_dir,
    data_fig_dir
  ),
  file.path(analysis_dir, "output_files.txt")
)

message("Done")
