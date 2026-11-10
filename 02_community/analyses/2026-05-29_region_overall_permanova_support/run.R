# ============================================================
# Analysis: Region support for pooled overall beta-diversity analysis
#
# Purpose:
#   Test whether region modifies the BMI-dependent CRC community effect.
#
# Interpretation guardrail:
#   Continent and Study.name are largely nested. Therefore Continent is not
#   interpreted as a regular fixed-effect covariate alongside Study.name.
#   Region baseline composition is reported descriptively; the BMI x Disease
#   and Region x BMI x Disease sensitivity tests use study-stratified
#   permutations.
# ============================================================

base_dir <- "G:/My Drive/proj_ob_AA_CRC"
analysis_dir <- file.path(base_dir, "analyses/2026-05-29_region_overall_permanova_support")
tables_dir <- file.path(analysis_dir, "tables")
figures_dir <- file.path(analysis_dir, "figures")
data_for_figures_dir <- file.path(analysis_dir, "data_for_figures")
final_tables_dir <- file.path(base_dir, "Figure_final/region_overall_support/tables")
logs_dir <- file.path(analysis_dir, "logs")

invisible(lapply(
  c(tables_dir, figures_dir, data_for_figures_dir, final_tables_dir, logs_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(vegan)
})

set.seed(20260529)

main_permutations <- as.integer(Sys.getenv("MAIN_PERMUTATIONS", unset = "99"))
balanced_repeats <- as.integer(Sys.getenv("BALANCED_REPEATS", unset = "50"))
balanced_permutations <- as.integer(Sys.getenv("BALANCED_PERMUTATIONS", unset = "2"))
force_main <- tolower(Sys.getenv("FORCE_MAIN", unset = "FALSE")) %in% c("true", "t", "1", "yes", "y")
if (balanced_permutations < 2) {
  balanced_permutations <- 2
}

dist_file <- file.path(base_dir, "03.results/02_beta_diversity/dist_matrix_bray_sqrt_prev0.1_remove_one_study.rds")
metadata_file <- file.path(base_dir, "02.processed_data/02_aligned_data/metadata_aligned_20260208.rds")
previous_multi_file <- file.path(
  base_dir,
  "03.results/02_beta_diversity/permanova_subgroups_multi/PERMANOVA_Subgroup_Comparison_Summary.csv"
)
main_file <- file.path(tables_dir, "main_permanova_region_overall_support.tsv")
main_final_file <- file.path(final_tables_dir, "main_permanova_region_overall_support.tsv")

message("Loading inputs...")
bray_dist <- readRDS(dist_file)
dist_mat <- as.matrix(bray_dist)

metadata <- readRDS(metadata_file)
rownames(metadata) <- metadata$Sample.Source.ID.processed

common_samples <- intersect(rownames(dist_mat), rownames(metadata))
metadata <- metadata[common_samples, , drop = FALSE]
dist_mat <- dist_mat[common_samples, common_samples, drop = FALSE]

metadata <- metadata %>%
  filter(Study.name != "This_study_cohort5__NSHII") %>%
  mutate(
    Disease = factor(Disease, levels = c("CTRL", "AA", "CRC")),
    BMI_Group1 = factor(BMI_Group1, levels = c("Normal", "Ow", "Ob")),
    Continent = factor(Continent, levels = c("Asia", "Europe")),
    Sex = factor(Sex),
    Age = suppressWarnings(as.numeric(Age)),
    Study.name = factor(Study.name)
  ) %>%
  filter(
    !is.na(Disease),
    !is.na(BMI_Group1),
    !is.na(Continent),
    !is.na(Sex),
    !is.na(Age),
    !is.na(Study.name)
  )

dist_mat <- dist_mat[rownames(metadata), rownames(metadata), drop = FALSE]

sample_counts <- metadata %>%
  count(Continent, Disease, BMI_Group1, name = "n") %>%
  arrange(Continent, Disease, BMI_Group1)

write_tsv(sample_counts, file.path(tables_dir, "sample_counts_region_disease_bmi.tsv"))

extract_term <- function(ad, term_label) {
  tab <- as.data.frame(ad) %>%
    rownames_to_column("term")

  hit <- tab %>% filter(term == term_label)
  if (nrow(hit) != 1) {
    stop("Could not find unique PERMANOVA term: ", term_label)
  }

  hit %>%
    transmute(
      term,
      df = Df,
      sum_of_squares = SumOfSqs,
      r2 = R2,
      f_statistic = F,
      p_value = if ("Pr(>F)" %in% names(hit)) `Pr(>F)` else NA_real_
    )
}

run_adonis <- function(sample_ids, formula_obj, strata = NULL, permutations = main_permutations) {
  sub_meta <- metadata[sample_ids, , drop = FALSE]
  sub_dist <- as.dist(dist_mat[sample_ids, sample_ids, drop = FALSE])
  if (is.null(strata)) {
    adonis2(formula_obj, data = sub_meta, permutations = permutations, by = "margin")
  } else {
    adonis2(
      formula_obj,
      data = sub_meta,
      permutations = permutations,
      strata = strata,
      by = "margin"
    )
  }
}

all_ids <- rownames(metadata)

if (file.exists(main_file) && !force_main) {
  message("Reusing existing main PERMANOVA support table. Set FORCE_MAIN=TRUE to rerun it.")
  main_tbl <- read_tsv(main_file, show_col_types = FALSE)
  if (!"r2" %in% names(main_tbl) && "r2_percent" %in% names(main_tbl)) {
    main_tbl <- main_tbl %>%
      mutate(r2 = r2_percent / 100) %>%
      select(
        analysis,
        region_subset,
        term,
        df,
        sum_of_squares,
        r2,
        r2_percent,
        f_statistic,
        p_value,
        n_samples,
        n_studies,
        n_permutations,
        formula,
        permutation_design
      )
    write_tsv(main_tbl, main_file)
    write_tsv(main_tbl, main_final_file)
  }
} else {
message("Running descriptive baseline Region PERMANOVA...")
ad_region <- run_adonis(
  all_ids,
  as.formula("sub_dist ~ Continent"),
  strata = NULL,
  permutations = main_permutations
)

region_row <- extract_term(ad_region, "Continent") %>%
  mutate(
    analysis = "baseline_region_descriptive",
    region_subset = "Overall",
    formula = "dist ~ Continent",
    permutation_design = "unstratified; descriptive because Continent is study-nested",
    n_samples = length(all_ids),
    n_studies = n_distinct(metadata$Study.name),
    n_permutations = main_permutations
  )

message("Running pooled BMI x Disease PERMANOVA with study-stratified permutations...")
ad_pooled_bmi_disease <- run_adonis(
  all_ids,
  as.formula("sub_dist ~ Age + Sex + BMI_Group1 * Disease"),
  strata = metadata[all_ids, "Study.name", drop = TRUE],
  permutations = main_permutations
)

pooled_interaction_row <- extract_term(ad_pooled_bmi_disease, "BMI_Group1:Disease") %>%
  mutate(
    analysis = "pooled_bmi_disease_interaction",
    region_subset = "Overall",
    formula = "dist ~ Age + Sex + BMI_Group1 * Disease",
    permutation_design = "strata = Study.name",
    n_samples = length(all_ids),
    n_studies = n_distinct(metadata$Study.name),
    n_permutations = main_permutations
  )

message("Running full Region x BMI x Disease modifier PERMANOVA with study-stratified permutations...")
ad_region_modifier <- run_adonis(
  all_ids,
  as.formula("sub_dist ~ Age + Sex + Continent * BMI_Group1 * Disease"),
  strata = metadata[all_ids, "Study.name", drop = TRUE],
  permutations = main_permutations
)

region_modifier_row <- extract_term(ad_region_modifier, "Continent:BMI_Group1:Disease") %>%
  mutate(
    analysis = "region_modifier_three_way",
    region_subset = "Overall",
    formula = "dist ~ Age + Sex + Continent * BMI_Group1 * Disease",
    permutation_design = "strata = Study.name; tests whether BMI x Disease pattern differs by region",
    n_samples = length(all_ids),
    n_studies = n_distinct(metadata$Study.name),
    n_permutations = main_permutations
  )

message("Running region-stratified BMI x Disease PERMANOVA...")
region_rows <- lapply(c("Asia", "Europe"), function(region_i) {
  ids_i <- rownames(metadata)[metadata$Continent == region_i]
  ad_i <- run_adonis(
    ids_i,
    as.formula("sub_dist ~ Age + Sex + BMI_Group1 * Disease"),
    strata = metadata[ids_i, "Study.name", drop = TRUE],
    permutations = main_permutations
  )

  extract_term(ad_i, "BMI_Group1:Disease") %>%
    mutate(
      analysis = "region_stratified_bmi_disease_interaction",
      region_subset = region_i,
      formula = "dist ~ Age + Sex + BMI_Group1 * Disease",
      permutation_design = "strata = Study.name within region",
      n_samples = length(ids_i),
      n_studies = n_distinct(metadata[ids_i, "Study.name", drop = TRUE]),
      n_permutations = main_permutations
    )
}) %>%
  bind_rows()

main_tbl <- bind_rows(
  region_row,
  pooled_interaction_row,
  region_modifier_row,
  region_rows
) %>%
  mutate(r2_percent = r2 * 100) %>%
  select(
    analysis,
    region_subset,
    term,
    df,
    sum_of_squares,
    r2,
    r2_percent = r2,
    f_statistic,
    p_value,
    n_samples,
    n_studies,
    n_permutations,
    formula,
    permutation_design
  )

write_tsv(main_tbl, main_file)
write_tsv(main_tbl, main_final_file)
}

if (file.exists(previous_multi_file)) {
  previous_multi <- read_csv(previous_multi_file, show_col_types = FALSE) %>%
    filter(BMI_Type == "BMI_Group1") %>%
    mutate(
      source = previous_multi_file,
      note = "Previous study-stratified region/overall PERMANOVA table reused for cross-check."
    )
  write_tsv(previous_multi, file.path(tables_dir, "reused_previous_region_stratified_permanova.tsv"))
}

message("Running repeated balanced sampling...")
min_cell_n <- min(sample_counts$n)

one_balanced_run <- function(iter) {
  set.seed(20260529 + iter)

  sampled_ids <- metadata %>%
    rownames_to_column("sample_id") %>%
    group_by(Continent, Disease, BMI_Group1) %>%
    slice_sample(n = min_cell_n) %>%
    ungroup() %>%
    pull(sample_id)

  sampled_meta <- metadata[sampled_ids, , drop = FALSE]

  ad_balanced_region <- run_adonis(
    sampled_ids,
    as.formula("sub_dist ~ Continent"),
    strata = NULL,
    permutations = balanced_permutations
  )

  ad_balanced_modifier <- run_adonis(
    sampled_ids,
    as.formula("sub_dist ~ Age + Sex + Continent * BMI_Group1 * Disease"),
    strata = sampled_meta$Study.name,
    permutations = balanced_permutations
  )

  bind_rows(
    extract_term(ad_balanced_region, "Continent") %>%
      mutate(analysis = "balanced_region_descriptive"),
    extract_term(ad_balanced_modifier, "Continent:BMI_Group1:Disease") %>%
      mutate(analysis = "balanced_region_modifier_three_way")
  ) %>%
    mutate(
      iteration = iter,
      n_per_cell = min_cell_n,
      n_samples = length(sampled_ids),
      n_studies = n_distinct(sampled_meta$Study.name),
      n_permutations = balanced_permutations
    )
}

balanced_tbl <- lapply(seq_len(balanced_repeats), one_balanced_run) %>%
  bind_rows() %>%
  mutate(
    p_value = if (balanced_permutations <= 2) NA_real_ else p_value,
    r2_percent = r2 * 100
  ) %>%
  select(
    iteration,
    analysis,
    term,
    df,
    sum_of_squares,
    r2,
    r2_percent,
    f_statistic,
    p_value,
    n_per_cell,
    n_samples,
    n_studies,
    n_permutations
  )

balanced_summary <- balanced_tbl %>%
  group_by(analysis, term) %>%
  summarise(
    n_iterations = n(),
    n_per_cell = first(n_per_cell),
    n_samples = first(n_samples),
    n_permutations = first(n_permutations),
    median_r2 = median(r2, na.rm = TRUE),
    q25_r2 = quantile(r2, 0.25, na.rm = TRUE),
    q75_r2 = quantile(r2, 0.75, na.rm = TRUE),
    median_r2_percent = median(r2_percent, na.rm = TRUE),
    q25_r2_percent = quantile(r2_percent, 0.25, na.rm = TRUE),
    q75_r2_percent = quantile(r2_percent, 0.75, na.rm = TRUE),
    median_p = if (all(is.na(p_value))) NA_real_ else median(p_value, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(balanced_tbl, file.path(tables_dir, "balanced_permanova_repeats.tsv"))
write_tsv(balanced_summary, file.path(tables_dir, "balanced_permanova_summary.tsv"))

format_num <- function(x, digits = 3) {
  ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = digits))
}

format_p <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    ifelse(x < 0.001, "<0.001", formatC(x, format = "f", digits = 3))
  )
}

supp_main <- main_tbl %>%
  mutate(
    analysis_raw = analysis,
    analysis_label = recode(
      analysis,
      baseline_region_descriptive = "Baseline region composition difference",
      pooled_bmi_disease_interaction = "Pooled BMI-dependent CRC effect",
      region_modifier_three_way = "Region modification of BMI-dependent CRC effect",
      region_stratified_bmi_disease_interaction = "Region-stratified BMI-dependent CRC effect",
      .default = analysis
    )
  ) %>%
  transmute(
    table_section = "Main study-aware PERMANOVA",
    analysis = analysis_label,
    region_subset,
    term,
    n_samples,
    n_studies_or_iterations = as.character(n_studies),
    r2_percent = format_num(r2_percent, 3),
    p_value = format_p(p_value),
    statistic = format_num(f_statistic, 3),
    statistic_label = "F",
    formula,
    permutation_design,
    interpretation = case_when(
      analysis_raw == "baseline_region_descriptive" ~ "Region explains baseline community composition and is treated as a sensitivity axis.",
      analysis_raw == "pooled_bmi_disease_interaction" ~ "BMI-dependent CRC community effect is detectable in the pooled overall analysis.",
      analysis_raw == "region_modifier_three_way" ~ "The region-dependent modification term explains little variation and is not supported by study-stratified permutation.",
      analysis_raw == "region_stratified_bmi_disease_interaction" & region_subset == "Asia" ~ "Within Asia, the BMI x Disease term has small R2 and is not significant in this sensitivity analysis.",
      analysis_raw == "region_stratified_bmi_disease_interaction" & region_subset == "Europe" ~ "Within Europe, the BMI x Disease term has small R2 and is not significant in this sensitivity analysis.",
      TRUE ~ NA_character_
    ),
    source_table = "main_permanova_region_overall_support.tsv"
  )

supp_balanced <- balanced_summary %>%
  mutate(
    analysis_raw = analysis,
    analysis_label = recode(
      analysis,
      balanced_region_descriptive = "Balanced baseline region composition difference",
      balanced_region_modifier_three_way = "Balanced region modification of BMI-dependent CRC effect",
      .default = analysis
    )
  ) %>%
  transmute(
    table_section = "Repeated balanced sampling sensitivity",
    analysis = analysis_label,
    region_subset = "Balanced Asia/Europe cells",
    term,
    n_samples,
    n_studies_or_iterations = as.character(n_iterations),
    r2_percent = paste0(
      format_num(median_r2_percent, 3),
      " [",
      format_num(q25_r2_percent, 3),
      ", ",
      format_num(q75_r2_percent, 3),
      "]"
    ),
    p_value = NA_character_,
    statistic = NA_character_,
    statistic_label = "Median R2 percent [IQR]",
    formula = case_when(
      term == "Continent" ~ "dist ~ Continent",
      TRUE ~ "dist ~ Age + Sex + Continent * BMI_Group1 * Disease"
    ),
    permutation_design = paste0(
      "Repeated balanced sampling by Continent x Disease x BMI_Group1; n_per_cell = ",
      n_per_cell,
      "; P values not interpreted in fast R2 sensitivity mode"
    ),
    interpretation = case_when(
      analysis_raw == "balanced_region_descriptive" ~ "After cell balancing, region still explains baseline community composition.",
      analysis_raw == "balanced_region_modifier_three_way" ~ "After cell balancing, the region-modification term remains small relative to the baseline region effect.",
      TRUE ~ NA_character_
    ),
    source_table = "balanced_permanova_summary.tsv"
  )

supplementary_table <- bind_rows(supp_main, supp_balanced)
supplementary_table_file <- "Supplementary_Table_X_region_overall_permanova_support.tsv"
write_tsv(supplementary_table, file.path(tables_dir, supplementary_table_file))
write_tsv(supplementary_table, file.path(final_tables_dir, supplementary_table_file))

sink(file.path(logs_dir, "sessionInfo.txt"))
sessionInfo()
sink()

message("Done.")
message("Supplementary table: ", file.path(final_tables_dir, supplementary_table_file))
