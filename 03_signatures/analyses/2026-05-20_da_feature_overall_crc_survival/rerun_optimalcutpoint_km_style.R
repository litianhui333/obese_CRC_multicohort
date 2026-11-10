suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(survival)
    library(survminer)
})

# Targeted rerun only: the full run.R would regenerate unrelated outputs.
base_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
analysis_dir <- file.path(base_dir, "analyses/2026-05-20_da_feature_overall_crc_survival")
figure_dir <- file.path(analysis_dir, "figures")
figure_data_dir <- file.path(analysis_dir, "data_for_figures")
surv_path <- file.path(base_dir, "00.rawdata/00_metadata/survival_time_info.csv")
abd_path <- file.path(base_dir, "02.processed_data/03_metaphlan_split_data/metaphlan_strain_renormalized_20260115.rds")

target_stems <- c(
    "ObSpecific_OptimalCutpoint_KM_Single_B._glucerasea_SGB4816",
    "ObSpecific_OptimalCutpoint_KM_Single_E._lavalensis_SGB4725",
    "ObSpecific_OptimalCutpoint_KM_Single_GGB3034_SGB4030",
    "ObSpecific_OptimalCutpoint_KM_Single_L._gasseri_SGB7038"
)

transform_ast_percent <- function(x) {
    x_prop <- pmin(pmax(x / 100, 0), 1)
    asin(sqrt(x_prop))
}

clean_survival <- function(surv_df) {
    surv_df %>%
        filter(!is.na(alive), !is.na(surv.time..Y.)) %>%
        mutate(
            status = ifelse(alive == FALSE, 1, 0),
            BMI = as.numeric(as.character(BMI)),
            BMI_Group = cut(
                BMI, c(0, 18.5, 24, 28, 100),
                c("Underweight", "Normal", "Overweight", "Obese"), right = FALSE
            ),
            Stage = case_when(
                grepl("^I$", AJCC) ~ "I",
                grepl("^II[A-C]?$", AJCC) ~ "II",
                grepl("^III[A-C]?$", AJCC) ~ "III",
                grepl("^IV[A-C]?$", AJCC) ~ "IV",
                TRUE ~ NA_character_
            ),
            Gender = factor(Gender, levels = c("F", "M")),
            BMI_Group = factor(BMI_Group, levels = c("Normal", "Underweight", "Overweight", "Obese")),
            Stage = factor(Stage, levels = c("I", "II", "III", "IV"))
        ) %>%
        drop_na(status, surv.time..Y., Age, Gender, BMI_Group)
}

fmt_p <- function(x) {
    ifelse(x < 0.001, format(x, scientific = TRUE, digits = 2), sprintf("%.3f", x))
}

required_inputs <- c(surv_path, abd_path, file.path(figure_data_dir, paste0(target_stems, ".tsv")))
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) stop("Missing required input(s):\n", paste(missing_inputs, collapse = "\n"))

surv_clean <- clean_survival(read_csv(surv_path, show_col_types = FALSE, name_repair = "universal"))
abd_mat <- readRDS(abd_path)
if (grepl("^s__", rownames(abd_mat)[1])) abd_mat <- t(abd_mat)
common_samples <- intersect(surv_clean$Sample, rownames(abd_mat))
surv_base <- surv_clean %>% filter(Sample %in% common_samples)
abd_final <- abd_mat[surv_base$Sample, , drop = FALSE]
abd_final <- abd_final[match(surv_base$Sample, rownames(abd_final)), , drop = FALSE]

for (stem in target_stems) {
    source_data_file <- file.path(figure_data_dir, paste0(stem, ".tsv"))
    prior_plot_data <- read_tsv(source_data_file, show_col_types = FALSE)
    required_columns <- c(
        "Feature", "Feature_Name", "Cutpoint_AST", "Exploratory_Cutpoint_Cox_HR",
        "Exploratory_Cutpoint_Cox_HR_Lower", "Exploratory_Cutpoint_Cox_HR_Upper",
        "Exploratory_Cutpoint_Cox_P", "Exploratory_Logrank_P",
        "Cutpoint_Raw_Abundance_Percent", "Cutpoint_Group_Counts", "Cutpoint_Selection_Method"
    )
    missing_columns <- setdiff(required_columns, names(prior_plot_data))
    if (length(missing_columns) > 0) stop(stem, " is missing: ", paste(missing_columns, collapse = ", "))

    meta <- prior_plot_data %>% slice(1)
    feature <- meta$Feature[[1]]
    if (!feature %in% colnames(abd_final)) stop("Feature not present in abundance matrix: ", feature)
    bacteria_val <- as.numeric(abd_final[, feature])
    n_before_stage_filter <- nrow(surv_base)

    plot_dat <- surv_base %>%
        mutate(Bacteria_Val = bacteria_val, Bacteria_AST = transform_ast_percent(bacteria_val)) %>%
        drop_na(Stage) %>%
        mutate(
            Bacteria_Group = factor(
                ifelse(Bacteria_AST > meta$Cutpoint_AST[[1]] + 1e-12, "High", "Low"),
                levels = c("Low", "High")
            ),
            Feature = feature,
            Feature_Name = meta$Feature_Name[[1]],
            Cutpoint_AST = meta$Cutpoint_AST[[1]],
            Cutpoint_Raw_Abundance_Percent = meta$Cutpoint_Raw_Abundance_Percent[[1]],
            Cutpoint_Group_Counts = meta$Cutpoint_Group_Counts[[1]],
            Exploratory_Cutpoint_Cox_HR = meta$Exploratory_Cutpoint_Cox_HR[[1]],
            Exploratory_Cutpoint_Cox_HR_Lower = meta$Exploratory_Cutpoint_Cox_HR_Lower[[1]],
            Exploratory_Cutpoint_Cox_HR_Upper = meta$Exploratory_Cutpoint_Cox_HR_Upper[[1]],
            Exploratory_Cutpoint_Cox_P = meta$Exploratory_Cutpoint_Cox_P[[1]],
            Exploratory_Logrank_P = meta$Exploratory_Logrank_P[[1]],
            Cutpoint_Selection_Method = meta$Cutpoint_Selection_Method[[1]],
            Stage_Filter_Before_N = n_before_stage_filter,
            Stage_Filter_After_N = n(),
            Plot_Panel = "KM curve and number-at-risk table"
        )

    observed_counts <- plot_dat %>%
        group_by(Bacteria_Group) %>%
        summarise(N = n(), Events = sum(status == 1), .groups = "drop") %>%
        mutate(label = paste0(Bacteria_Group, ":N", N, "_E", Events)) %>%
        pull(label) %>%
        paste(collapse = "; ")
    if (!identical(observed_counts, meta$Cutpoint_Group_Counts[[1]])) {
        stop("Group-count validation failed for ", stem, ". Expected '", meta$Cutpoint_Group_Counts[[1]], "' but reconstructed '", observed_counts, "'.")
    }

    km_fit <- survfit(Surv(surv.time..Y., status) ~ Bacteria_Group, data = plot_dat)
    title_str <- paste0(
        meta$Feature_Name[[1]], "\nOb-specific | Optimal cutpoint | WithStage",
        "\nHR=", sprintf("%.2f", meta$Exploratory_Cutpoint_Cox_HR[[1]]),
        " (", sprintf("%.2f", meta$Exploratory_Cutpoint_Cox_HR_Lower[[1]]),
        "-", sprintf("%.2f", meta$Exploratory_Cutpoint_Cox_HR_Upper[[1]]), ")",
        ", P=", fmt_p(meta$Exploratory_Cutpoint_Cox_P[[1]])
    )

    p_km <- ggsurvplot(
        km_fit, data = plot_dat,
        pval = paste0("Cox P = ", fmt_p(meta$Exploratory_Cutpoint_Cox_P[[1]])),
        conf.int = FALSE, risk.table = TRUE,
        palette = c("#6A8CAF", "#D94F45"),
        title = title_str,
        xlab = "Time (years)", ylab = "Overall survival probability",
        legend.title = "Optimal cutpoint", legend.labs = c("Low", "High"),
        ggtheme = theme_classic(base_size = 12, base_family = "Arial"),
        risk.table.height = 0.24,
        font.title = c(12, "plain", "black"),
        font.x = c(12, "plain", "black"), font.y = c(12, "plain", "black"),
        font.tickslab = c(10, "plain", "black"), font.legend = c(10, "plain", "black"),
        risk.table.fontsize = 5
    )

    write_tsv(plot_dat, source_data_file)
    grDevices::cairo_pdf(file.path(figure_dir, paste0(stem, ".pdf")), width = 5.8, height = 6.8, family = "Arial", onefile = FALSE)
    print(p_km)
    dev.off()
    png(file.path(figure_dir, paste0(stem, ".png")), width = 5.8, height = 6.8, units = "in", res = 300)
    print(p_km)
    dev.off()
    message("Updated: ", stem, " (stage-complete N=", nrow(plot_dat), ")")
}

message("Targeted KM style rerun complete. No unrelated outputs were regenerated.")
