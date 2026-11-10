suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(survival)
    library(survminer)
    library(ggplot2)
})

base_dir <- "G:/My Drive/proj_ob_AA_CRC"
analysis_dir <- file.path(base_dir, "analyses/2026-05-20_da_feature_overall_crc_survival")
table_dir <- file.path(analysis_dir, "tables")
figure_dir <- file.path(analysis_dir, "figures")
figure_data_dir <- file.path(analysis_dir, "data_for_figures")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_data_dir, recursive = TRUE, showWarnings = FALSE)

surv_path <- file.path(base_dir, "00.rawdata/00_metadata/survival_time_info.csv")
abd_path <- file.path(base_dir, "02.processed_data/03_metaphlan_split_data/metaphlan_strain_renormalized_20260115.rds")
feature_path <- file.path(
    base_dir,
    "03.results/04_differential_abundance/heatmap_final/BMI_Stratified_Disease_Heatmaps/Feature_List_Ob_Refined.csv"
)

safe_name <- function(x) {
    gsub("[|/\\\\: *?\"<>]", "_", x)
}

format_sgb_names <- function(x) {
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
            res <- paste0(substr(genus, 1, 1), ". ", species)
        } else {
            res <- gsub("_", " ", sp_part)
        }
    }

    final_res <- trimws(res)
    if (sgb_id != "" && !grepl(sgb_id, final_res, fixed = TRUE)) {
        final_res <- paste(final_res, sgb_id)
    }
    final_res
}

transform_ast_percent <- function(x) {
    x_prop <- x / 100
    x_prop <- pmin(pmax(x_prop, 0), 1)
    asin(sqrt(x_prop))
}

clean_survival <- function(surv_df) {
    surv_df %>%
        filter(!is.na(alive), !is.na(surv.time..Y.)) %>%
        mutate(
            status = ifelse(alive == FALSE, 1, 0),
            BMI = as.numeric(as.character(BMI)),
            BMI_Group = cut(
                BMI,
                breaks = c(0, 18.5, 24, 28, 100),
                labels = c("Underweight", "Normal", "Overweight", "Obese"),
                right = FALSE
            ),
            Stage_Raw = AJCC,
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

extract_group_counts <- function(dat, group_col = "Bacteria_Group") {
    dat %>%
        group_by(.data[[group_col]]) %>%
        summarise(
            N_Group = n(),
            N_Event_Group = sum(status == 1, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        mutate(label = paste0(.data[[group_col]], ":N", N_Group, "_E", N_Event_Group)) %>%
        pull(label) %>%
        paste(collapse = "; ")
}

group_is_valid <- function(dat, group_col = "Bacteria_Group", min_group_n = 5, min_group_events = 2) {
    group_summary <- dat %>%
        group_by(.data[[group_col]]) %>%
        summarise(
            N_Group = n(),
            N_Event_Group = sum(status == 1, na.rm = TRUE),
            .groups = "drop"
        )

    nrow(group_summary) == 2 &&
        min(group_summary$N_Group) >= min_group_n &&
        min(group_summary$N_Event_Group) >= min_group_events
}

run_cox <- function(dat, formula_obj, feature, model_type, grouping_type, feature_meta) {
    warn_msgs <- character()
    fit <- tryCatch(
        withCallingHandlers(
            coxph(formula_obj, data = dat, control = coxph.control(iter.max = 50)),
            warning = function(w) {
                warn_msgs <<- c(warn_msgs, conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        ),
        error = function(e) NULL
    )

    if (is.null(fit)) {
        return(NULL)
    }

    sum_fit <- summary(fit)
    coef_names <- rownames(sum_fit$coefficients)
    target_row <- grep("^Bacteria_Group", coef_names)
    if (length(target_row) == 0) {
        return(NULL)
    }

    coef_row <- sum_fit$coefficients[target_row[1], , drop = FALSE]
    conf_row <- sum_fit$conf.int[target_row[1], , drop = FALSE]

    zph_p <- tryCatch({
        z <- cox.zph(fit)
        z$table["GLOBAL", "p"]
    }, error = function(e) NA_real_)
    meta <- feature_meta %>%
        filter(Feature == feature) %>%
        slice(1)

    data.frame(
        Feature = feature,
        Feature_Name = meta$Feature_Name,
        Ob_DA_Class = meta$Ob_DA_Class,
        Candidate_Set = meta$Candidate_Set,
        Grouping = grouping_type,
        Model_Type = model_type,
        Term = rownames(coef_row)[1],
        HR = unname(conf_row[1, "exp(coef)"]),
        HR_Lower = unname(conf_row[1, "lower .95"]),
        HR_Upper = unname(conf_row[1, "upper .95"]),
        P_Value = unname(coef_row[1, "Pr(>|z|)"]),
        N = sum_fit$n,
        N_Event = sum_fit$nevent,
        Group_Counts = extract_group_counts(dat),
        PH_Global_P = zph_p,
        Cox_Warnings = paste(unique(warn_msgs), collapse = " | "),
        Median_Abundance = median(dat$Bacteria_Val, na.rm = TRUE),
        Prevalence = mean(dat$Bacteria_Val > 0, na.rm = TRUE),
        Ob_CRC_Coef = meta$Ob_CRC_Coef,
        Ob_CRC_P = meta$Ob_CRC_P,
        Ob_CRC_Q = meta$Ob_CRC_Q,
        stringsAsFactors = FALSE
    )
}

run_cox_continuous_ast <- function(dat, feature, feature_meta) {
    warn_msgs <- character()
    fit <- tryCatch(
        withCallingHandlers(
            coxph(
                Surv(surv.time..Y., status) ~ Bacteria_AST + Age + Gender + BMI_Group + Stage,
                data = dat,
                control = coxph.control(iter.max = 50)
            ),
            warning = function(w) {
                warn_msgs <<- c(warn_msgs, conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        ),
        error = function(e) NULL
    )

    if (is.null(fit)) {
        return(NULL)
    }

    sum_fit <- summary(fit)
    if (!"Bacteria_AST" %in% rownames(sum_fit$coefficients)) {
        return(NULL)
    }

    coef_row <- sum_fit$coefficients["Bacteria_AST", , drop = FALSE]
    beta <- unname(coef_row[1, "coef"])
    se <- unname(coef_row[1, "se(coef)"])
    p_val <- unname(coef_row[1, "Pr(>|z|)"])
    scale_step <- 0.1
    ast_sd <- sd(dat$Bacteria_AST, na.rm = TRUE)
    ast_iqr <- IQR(dat$Bacteria_AST, na.rm = TRUE)

    zph_p <- tryCatch({
        z <- cox.zph(fit)
        z$table["GLOBAL", "p"]
    }, error = function(e) NA_real_)
    zph_taxon_p <- tryCatch({
        z <- cox.zph(fit)
        z$table["Bacteria_AST", "p"]
    }, error = function(e) NA_real_)

    meta <- feature_meta %>%
        filter(Feature == feature) %>%
        slice(1)

    data.frame(
        Feature = feature,
        Feature_Name = meta$Feature_Name,
        Ob_DA_Class = meta$Ob_DA_Class,
        Candidate_Set = meta$Candidate_Set,
        Model_Type = "Continuous_AST_WithStage",
        Term = "Bacteria_AST",
        Coef_AST = beta,
        SE_AST = se,
        HR_Per_1_AST = exp(beta),
        HR_Per_0.1_AST = exp(beta * scale_step),
        HR_Per_0.1_AST_Lower = exp((beta - 1.96 * se) * scale_step),
        HR_Per_0.1_AST_Upper = exp((beta + 1.96 * se) * scale_step),
        HR_Per_SD_AST = exp(beta * ast_sd),
        HR_Per_SD_AST_Lower = exp((beta - 1.96 * se) * ast_sd),
        HR_Per_SD_AST_Upper = exp((beta + 1.96 * se) * ast_sd),
        HR_Per_IQR_AST = exp(beta * ast_iqr),
        HR_Per_IQR_AST_Lower = exp((beta - 1.96 * se) * ast_iqr),
        HR_Per_IQR_AST_Upper = exp((beta + 1.96 * se) * ast_iqr),
        P_Value = p_val,
        N = sum_fit$n,
        N_Event = sum_fit$nevent,
        PH_Taxon_P = zph_taxon_p,
        PH_Global_P = zph_p,
        Cox_Warnings = paste(unique(warn_msgs), collapse = " | "),
        Raw_Median_Abundance = median(dat$Bacteria_Val, na.rm = TRUE),
        Raw_Prevalence = mean(dat$Bacteria_Val > 0, na.rm = TRUE),
        AST_Mean = mean(dat$Bacteria_AST, na.rm = TRUE),
        AST_SD = ast_sd,
        AST_IQR = ast_iqr,
        AST_Min = min(dat$Bacteria_AST, na.rm = TRUE),
        AST_Max = max(dat$Bacteria_AST, na.rm = TRUE),
        Ob_CRC_Coef = meta$Ob_CRC_Coef,
        Ob_CRC_P = meta$Ob_CRC_P,
        Ob_CRC_Q = meta$Ob_CRC_Q,
        stringsAsFactors = FALSE
    )
}

km_curve_data <- function(km_fit, feature, feature_name, grouping, model_type) {
    s <- summary(km_fit)
    strata <- if (is.null(s$strata)) rep(NA_character_, length(s$time)) else as.character(s$strata)

    tibble(
        Feature = feature,
        Feature_Name = feature_name,
        Grouping = grouping,
        Model_Type = model_type,
        strata = strata,
        time = s$time,
        n_risk = s$n.risk,
        n_event = s$n.event,
        n_censor = s$n.censor,
        survival = s$surv,
        std_error = s$std.err,
        lower = s$lower,
        upper = s$upper
    )
}

plot_km <- function(dat, feature, grouping, model_type, model_res) {
    km_fit <- survfit(Surv(surv.time..Y., status) ~ Bacteria_Group, data = dat)
    prefix <- paste0("KM_", grouping, "_", model_type, "_", safe_name(feature))

    write_tsv(
        km_curve_data(km_fit, feature, model_res$Feature_Name, grouping, model_type),
        file.path(figure_data_dir, paste0(prefix, ".tsv"))
    )

    title_str <- paste0(
        model_res$Feature_Name,
        "\n", model_res$Ob_DA_Class, " | ", grouping, " | ", model_type,
        "\nHR=", sprintf("%.2f", model_res$HR),
        " (", sprintf("%.2f", model_res$HR_Lower), "-", sprintf("%.2f", model_res$HR_Upper), ")",
        ", P=", sprintf("%.3g", model_res$P_Value)
    )

    pal <- if (grouping == "PA") c("grey55", "#D39200") else c("#6A8CAF", "#D94F45")

    p_km <- ggsurvplot(
        km_fit,
        data = dat,
        pval = paste0("Cox P = ", sprintf("%.3g", model_res$P_Value)),
        conf.int = FALSE,
        risk.table = TRUE,
        palette = pal,
        title = title_str,
        xlab = "Time (years)",
        ylab = "Overall survival probability",
        legend.title = grouping,
        legend.labs = levels(dat$Bacteria_Group),
        ggtheme = theme_classic(base_size = 12),
        risk.table.height = 0.24,
        font.title = c(12, "plain", "black"),
        font.x = c(12, "plain", "black"),
        font.y = c(12, "plain", "black"),
        font.tickslab = c(10, "plain", "black"),
        font.legend = c(10, "plain", "black"),
        risk.table.fontsize = 3.5
    )

    pdf_file <- file.path(figure_dir, paste0(prefix, ".pdf"))
    png_file <- file.path(figure_dir, paste0(prefix, ".png"))

    pdf(pdf_file, width = 5.8, height = 6.8, onefile = FALSE)
    print(p_km)
    dev.off()

    png(png_file, width = 5.8, height = 6.8, units = "in", res = 300)
    print(p_km)
    dev.off()
}

message("=== 1. Load inputs ===")
surv_df <- read.csv(surv_path, stringsAsFactors = FALSE)
abd_mat_raw <- readRDS(abd_path)
feature_raw <- read.csv(feature_path, stringsAsFactors = FALSE, check.names = FALSE)

if (grepl("^s__", rownames(abd_mat_raw)[1])) {
    abd_mat_raw <- t(abd_mat_raw)
}

message("=== 2. Prepare survival cohort ===")
surv_clean <- clean_survival(surv_df)
common_samples <- intersect(surv_clean$Sample, rownames(abd_mat_raw))
surv_base <- surv_clean %>%
    filter(Sample %in% common_samples)
abd_final <- abd_mat_raw[surv_base$Sample, , drop = FALSE]
abd_final <- abd_final[match(surv_base$Sample, rownames(abd_final)), , drop = FALSE]

cohort_summary <- surv_base %>%
    summarise(
        N_NoStage_Model_Base = n(),
        N_Event_NoStage_Model_Base = sum(status == 1, na.rm = TRUE),
        N_Censored_NoStage_Model_Base = sum(status == 0, na.rm = TRUE),
        N_With_Stage = sum(!is.na(Stage)),
        N_Event_With_Stage = sum(status == 1 & !is.na(Stage), na.rm = TRUE),
        Median_Followup_Years = median(surv.time..Y., na.rm = TRUE),
        N_Normal_BMI = sum(BMI_Group == "Normal", na.rm = TRUE),
        N_Underweight_BMI = sum(BMI_Group == "Underweight", na.rm = TRUE),
        N_Overweight_BMI = sum(BMI_Group == "Overweight", na.rm = TRUE),
        N_Obese_BMI = sum(BMI_Group == "Obese", na.rm = TRUE)
    )

write_csv(cohort_summary, file.path(table_dir, "Overall_CRC_Cohort_Summary.csv"))
print(cohort_summary)

message("=== 3. Prepare DA candidate taxa ===")
candidate_df <- feature_raw %>%
    transmute(
        Feature = feature,
        Feature_Name = if ("clean_name" %in% names(feature_raw)) clean_name else vapply(feature, format_sgb_names, character(1)),
        Ob_DA_Class = overlap,
        Candidate_Set = ifelse(overlap == "Ob-specific", "Ob-specific", "Ob-overlap"),
        Ob_CRC_Coef = if ("coef_crc" %in% names(feature_raw)) coef_crc else NA_real_,
        Ob_CRC_P = if ("pval_crc" %in% names(feature_raw)) pval_crc else NA_real_,
        Ob_CRC_Q = if ("qval_crc" %in% names(feature_raw)) qval_crc else NA_real_
    ) %>%
    distinct(Feature, .keep_all = TRUE) %>%
    mutate(
        Available_In_Abundance = Feature %in% colnames(abd_final),
        Feature_Source = "Obese BMI stratum CTRL vs CRC DA feature list"
    )

write_csv(candidate_df, file.path(table_dir, "Overall_CRC_DA_Feature_Candidate_Taxa.csv"))

target_features <- candidate_df %>%
    filter(Available_In_Abundance) %>%
    pull(Feature)

message("Candidate taxa: ", nrow(candidate_df))
message("Available in survival abundance matrix: ", length(target_features))
print(table(candidate_df$Ob_DA_Class, candidate_df$Available_In_Abundance, useNA = "ifany"))

feature_meta <- candidate_df %>%
    select(Feature, Feature_Name, Ob_DA_Class, Candidate_Set, Ob_CRC_Coef, Ob_CRC_P, Ob_CRC_Q)

message("=== 4. Run Cox models ===")
results <- list()
skipped <- list()
continuous_results <- list()
continuous_skipped <- list()

for (feature in target_features) {
    val <- as.numeric(abd_final[, feature])
    dat_base <- surv_base %>%
        mutate(
            Bacteria_Val = val,
            Bacteria_AST = transform_ast_percent(val)
        )

    dat_continuous <- dat_base %>% drop_na(Stage)
    if (sd(dat_continuous$Bacteria_AST, na.rm = TRUE) > 0) {
        res_continuous <- run_cox_continuous_ast(dat_continuous, feature, feature_meta)
        if (is.null(res_continuous)) {
            continuous_skipped[[length(continuous_skipped) + 1]] <- tibble(
                Feature = feature,
                Feature_Name = feature_meta$Feature_Name[match(feature, feature_meta$Feature)],
                Reason = "Continuous AST Cox model failed or term unavailable",
                N = nrow(dat_continuous),
                N_Event = sum(dat_continuous$status == 1, na.rm = TRUE),
                AST_SD = sd(dat_continuous$Bacteria_AST, na.rm = TRUE)
            )
        } else {
            continuous_results[[length(continuous_results) + 1]] <- res_continuous
        }
    } else {
        continuous_skipped[[length(continuous_skipped) + 1]] <- tibble(
            Feature = feature,
            Feature_Name = feature_meta$Feature_Name[match(feature, feature_meta$Feature)],
            Reason = "Continuous AST abundance has no variation",
            N = nrow(dat_continuous),
            N_Event = sum(dat_continuous$status == 1, na.rm = TRUE),
            AST_SD = sd(dat_continuous$Bacteria_AST, na.rm = TRUE)
        )
    }

    med_val <- median(val, na.rm = TRUE)
    dat_median <- dat_base %>%
        mutate(Bacteria_Group = factor(ifelse(Bacteria_Val >= med_val, "High", "Low"), levels = c("Low", "High")))

    dat_pa <- dat_base %>%
        mutate(Bacteria_Group = factor(ifelse(Bacteria_Val > 0, "Present", "Absent"), levels = c("Absent", "Present")))

    model_specs <- list(
        list(
            grouping = "Median",
            dat = dat_median,
            with_stage_formula = Surv(surv.time..Y., status) ~ Bacteria_Group + Age + Gender + BMI_Group + Stage,
            no_stage_formula = Surv(surv.time..Y., status) ~ Bacteria_Group + Age + Gender + BMI_Group
        ),
        list(
            grouping = "PA",
            dat = dat_pa,
            with_stage_formula = Surv(surv.time..Y., status) ~ Bacteria_Group + Age + Gender + BMI_Group + Stage,
            no_stage_formula = Surv(surv.time..Y., status) ~ Bacteria_Group + Age + Gender + BMI_Group
        )
    )

    for (spec in model_specs) {
        dat_no_stage <- spec$dat
        dat_with_stage <- spec$dat %>% drop_na(Stage)

        for (model_type in c("WithStage", "NoStage")) {
            dat_model <- if (model_type == "WithStage") dat_with_stage else dat_no_stage
            valid <- group_is_valid(dat_model)

            if (!valid) {
                skipped[[length(skipped) + 1]] <- tibble(
                    Feature = feature,
                    Feature_Name = feature_meta$Feature_Name[match(feature, feature_meta$Feature)],
                    Grouping = spec$grouping,
                    Model_Type = model_type,
                    Reason = "Group balance or event count below threshold",
                    Group_Counts = extract_group_counts(dat_model)
                )
                next
            }

            formula_obj <- if (model_type == "WithStage") spec$with_stage_formula else spec$no_stage_formula
            res <- run_cox(dat_model, formula_obj, feature, model_type, spec$grouping, feature_meta)

            if (is.null(res)) {
                skipped[[length(skipped) + 1]] <- tibble(
                    Feature = feature,
                    Feature_Name = feature_meta$Feature_Name[match(feature, feature_meta$Feature)],
                    Grouping = spec$grouping,
                    Model_Type = model_type,
                    Reason = "Cox model failed or taxon term unavailable",
                    Group_Counts = extract_group_counts(dat_model)
                )
                next
            }

            results[[length(results) + 1]] <- res

            if (model_type == "WithStage" && res$P_Value < 0.05) {
                plot_km(dat_model, feature, spec$grouping, model_type, res)
            }
        }
    }
}

if (length(results) == 0) {
    stop("No valid Cox models were produced.")
}

final_df <- bind_rows(results) %>%
    group_by(Grouping, Model_Type) %>%
    mutate(P_Adj_Within_Model = p.adjust(P_Value, method = "BH")) %>%
    ungroup() %>%
    mutate(P_Adj_All = p.adjust(P_Value, method = "BH")) %>%
    arrange(P_Value)

write_csv(final_df, file.path(table_dir, "Overall_CRC_DA_Feature_Survival_Cox_Results.csv"))

sig_nominal <- final_df %>%
    filter(P_Value < 0.05) %>%
    arrange(P_Value)
write_csv(sig_nominal, file.path(table_dir, "Overall_CRC_DA_Feature_Survival_Cox_P005.csv"))

sig_fdr_model <- final_df %>%
    filter(P_Adj_Within_Model < 0.10) %>%
    arrange(P_Adj_Within_Model)
write_csv(sig_fdr_model, file.path(table_dir, "Overall_CRC_DA_Feature_Survival_Cox_FDR10_WithinModel.csv"))

summary_table <- final_df %>%
    group_by(Candidate_Set, Ob_DA_Class, Grouping, Model_Type) %>%
    summarise(
        N_Tests = n(),
        N_P005 = sum(P_Value < 0.05, na.rm = TRUE),
        N_FDR10_Within_Model = sum(P_Adj_Within_Model < 0.10, na.rm = TRUE),
        Min_P = min(P_Value, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(Grouping, Model_Type, Candidate_Set, Ob_DA_Class)
write_csv(summary_table, file.path(table_dir, "Overall_CRC_DA_Feature_Survival_Significance_Summary.csv"))

skipped_df <- if (length(skipped) > 0) bind_rows(skipped) else tibble()
write_csv(skipped_df, file.path(table_dir, "Overall_CRC_DA_Feature_Survival_Skipped_Models.csv"))

continuous_df <- if (length(continuous_results) > 0) {
    bind_rows(continuous_results) %>%
        mutate(P_Adj_FDR = p.adjust(P_Value, method = "BH")) %>%
        arrange(P_Value)
} else {
    tibble()
}
write_csv(continuous_df, file.path(table_dir, "Overall_CRC_DA_Feature_Survival_Continuous_AST_WithStage.csv"))

continuous_p005 <- continuous_df %>%
    filter(P_Value < 0.05) %>%
    arrange(P_Value)
write_csv(continuous_p005, file.path(table_dir, "Overall_CRC_DA_Feature_Survival_Continuous_AST_WithStage_P005.csv"))

continuous_fdr10 <- continuous_df %>%
    filter(P_Adj_FDR < 0.10) %>%
    arrange(P_Adj_FDR)
write_csv(continuous_fdr10, file.path(table_dir, "Overall_CRC_DA_Feature_Survival_Continuous_AST_WithStage_FDR10.csv"))

continuous_skipped_df <- if (length(continuous_skipped) > 0) bind_rows(continuous_skipped) else tibble()
write_csv(continuous_skipped_df, file.path(table_dir, "Overall_CRC_DA_Feature_Survival_Continuous_AST_Skipped.csv"))

fmt_hr_ci <- function(hr, lower, upper) {
    ifelse(
        is.na(hr),
        NA_character_,
        sprintf("%.2f (%.2f-%.2f)", hr, lower, upper)
    )
}

fmt_p <- function(p) {
    ifelse(is.na(p), NA_character_, formatC(p, format = "g", digits = 3))
}

binary_with_stage <- final_df %>%
    filter(Model_Type == "WithStage") %>%
    mutate(
        HR_CI = fmt_hr_ci(HR, HR_Lower, HR_Upper),
        P_Text = fmt_p(P_Value),
        FDR_Text = fmt_p(P_Adj_Within_Model)
    )

median_wide <- binary_with_stage %>%
    filter(Grouping == "Median") %>%
    transmute(
        Feature,
        Median_HR = HR,
        Median_HR_Lower = HR_Lower,
        Median_HR_Upper = HR_Upper,
        Median_HR_CI = HR_CI,
        Median_P = P_Value,
        Median_FDR = P_Adj_Within_Model,
        Median_P_Text = P_Text,
        Median_FDR_Text = FDR_Text,
        Median_N = N,
        Median_N_Event = N_Event,
        Median_Group_Counts = Group_Counts
    )

pa_wide <- binary_with_stage %>%
    filter(Grouping == "PA") %>%
    transmute(
        Feature,
        PA_HR = HR,
        PA_HR_Lower = HR_Lower,
        PA_HR_Upper = HR_Upper,
        PA_HR_CI = HR_CI,
        PA_P = P_Value,
        PA_FDR = P_Adj_Within_Model,
        PA_P_Text = P_Text,
        PA_FDR_Text = FDR_Text,
        PA_N = N,
        PA_N_Event = N_Event,
        PA_Group_Counts = Group_Counts
    )

continuous_wide <- continuous_df %>%
    transmute(
        Feature,
        Continuous_HR_Per_SD_AST = HR_Per_SD_AST,
        Continuous_HR_Per_SD_AST_Lower = HR_Per_SD_AST_Lower,
        Continuous_HR_Per_SD_AST_Upper = HR_Per_SD_AST_Upper,
        Continuous_HR_Per_SD_AST_CI = fmt_hr_ci(HR_Per_SD_AST, HR_Per_SD_AST_Lower, HR_Per_SD_AST_Upper),
        Continuous_HR_Per_0.1_AST = HR_Per_0.1_AST,
        Continuous_HR_Per_0.1_AST_Lower = HR_Per_0.1_AST_Lower,
        Continuous_HR_Per_0.1_AST_Upper = HR_Per_0.1_AST_Upper,
        Continuous_HR_Per_0.1_AST_CI = fmt_hr_ci(HR_Per_0.1_AST, HR_Per_0.1_AST_Lower, HR_Per_0.1_AST_Upper),
        Continuous_HR_Per_IQR_AST = HR_Per_IQR_AST,
        Continuous_P = P_Value,
        Continuous_FDR = P_Adj_FDR,
        Continuous_P_Text = fmt_p(P_Value),
        Continuous_FDR_Text = fmt_p(P_Adj_FDR),
        Continuous_N = N,
        Continuous_N_Event = N_Event,
        Raw_Median_Abundance,
        Raw_Prevalence,
        AST_SD,
        AST_IQR,
        Continuous_PH_Taxon_P = PH_Taxon_P,
        Continuous_PH_Global_P = PH_Global_P
    )

integrated_survival <- candidate_df %>%
    select(
        Feature,
        Feature_Name,
        Ob_DA_Class,
        Candidate_Set,
        Ob_CRC_Coef,
        Ob_CRC_P,
        Ob_CRC_Q,
        Available_In_Abundance
    ) %>%
    left_join(continuous_wide, by = "Feature") %>%
    left_join(median_wide, by = "Feature") %>%
    left_join(pa_wide, by = "Feature") %>%
    rowwise() %>%
    mutate(
        DA_Direction = case_when(
            Ob_CRC_Coef > 0 ~ "CRC-enriched in obese DA",
            Ob_CRC_Coef < 0 ~ "CRC-depleted in obese DA",
            TRUE ~ NA_character_
        ),
        Continuous_Direction = case_when(
            is.na(Continuous_HR_Per_SD_AST) ~ NA_character_,
            Continuous_HR_Per_SD_AST > 1 ~ "Higher abundance, higher hazard",
            Continuous_HR_Per_SD_AST < 1 ~ "Higher abundance, lower hazard",
            TRUE ~ "No HR direction"
        ),
        Median_Direction = case_when(
            is.na(Median_HR) ~ NA_character_,
            Median_HR > 1 ~ "High vs low, higher hazard",
            Median_HR < 1 ~ "High vs low, lower hazard",
            TRUE ~ "No HR direction"
        ),
        PA_Direction = case_when(
            is.na(PA_HR) ~ NA_character_,
            PA_HR > 1 ~ "Present vs absent, higher hazard",
            PA_HR < 1 ~ "Present vs absent, lower hazard",
            TRUE ~ "No HR direction"
        ),
        Continuous_Nominal = !is.na(Continuous_P) && Continuous_P < 0.05,
        Median_Nominal = !is.na(Median_P) && Median_P < 0.05,
        PA_Nominal = !is.na(PA_P) && PA_P < 0.05,
        Any_Grouped_Nominal = Median_Nominal || PA_Nominal,
        Any_Nominal = Continuous_Nominal || Any_Grouped_Nominal,
        Any_FDR10 = (!is.na(Continuous_FDR) && Continuous_FDR < 0.10) ||
            (!is.na(Median_FDR) && Median_FDR < 0.10) ||
            (!is.na(PA_FDR) && PA_FDR < 0.10),
        Significant_Models = paste(
            c(
                if (Continuous_Nominal) "Continuous_AST" else NULL,
                if (Median_Nominal) "Median" else NULL,
                if (PA_Nominal) "Presence_absence" else NULL
            ),
            collapse = "; "
        ),
        Significant_Models = ifelse(Significant_Models == "", NA_character_, Significant_Models),
        Best_Grouped_Model = {
            p_vals <- c(Median = Median_P, Presence_absence = PA_P)
            if (all(is.na(p_vals))) {
                NA_character_
            } else {
                names(p_vals)[which.min(replace(p_vals, is.na(p_vals), Inf))]
            }
        },
        Best_Grouped_P = {
            p_vals <- c(Median_P, PA_P)
            if (all(is.na(p_vals))) NA_real_ else min(p_vals, na.rm = TRUE)
        },
        Best_Grouped_HR = case_when(
            Best_Grouped_Model == "Median" ~ Median_HR,
            Best_Grouped_Model == "Presence_absence" ~ PA_HR,
            TRUE ~ NA_real_
        ),
        Best_Grouped_HR_CI = case_when(
            Best_Grouped_Model == "Median" ~ Median_HR_CI,
            Best_Grouped_Model == "Presence_absence" ~ PA_HR_CI,
            TRUE ~ NA_character_
        ),
        Evidence_Class = case_when(
            Any_FDR10 && Continuous_Nominal && Any_Grouped_Nominal ~ "FDR-supported plus concordant nominal models",
            Any_FDR10 ~ "FDR-supported in at least one model",
            Continuous_Nominal && Any_Grouped_Nominal ~ "Continuous and grouped nominal",
            Continuous_Nominal ~ "Continuous-only nominal",
            Any_Grouped_Nominal ~ "Grouped-only nominal",
            TRUE ~ "No nominal survival association"
        ),
        Suggested_Manuscript_Use = case_when(
            Any_FDR10 ~ "Highlight as strongest exploratory survival candidate",
            Continuous_Nominal && Any_Grouped_Nominal ~ "Describe as model-consistent supplementary candidate",
            Continuous_Nominal ~ "Describe as continuous-abundance exploratory candidate; avoid relying on KM dichotomization alone",
            Any_Grouped_Nominal ~ "Describe as dichotomized/KM exploratory candidate; continuous model did not support it",
            TRUE ~ "Do not highlight"
        )
    ) %>%
    ungroup() %>%
    mutate(
        Ob_CRC_P_Text = fmt_p(Ob_CRC_P),
        Ob_CRC_Q_Text = fmt_p(Ob_CRC_Q)
    ) %>%
    arrange(
        desc(Any_FDR10),
        desc(Continuous_Nominal & Any_Grouped_Nominal),
        desc(Continuous_Nominal),
        desc(Any_Grouped_Nominal),
        coalesce(Continuous_P, Best_Grouped_P, 1)
    )

write_csv(
    integrated_survival,
    file.path(table_dir, "Overall_CRC_DA_Feature_Survival_Integrated_All_Taxa.csv")
)

manuscript_candidate_table <- integrated_survival %>%
    filter(Any_Nominal) %>%
    transmute(
        Feature_Name,
        Ob_DA_Class,
        DA_Direction,
        Ob_CRC_Coef,
        Ob_CRC_Q,
        Continuous_HR_Per_SD_AST_CI,
        Continuous_P,
        Continuous_FDR,
        Median_HR_CI,
        Median_P,
        Median_FDR,
        PA_HR_CI,
        PA_P,
        PA_FDR,
        Significant_Models,
        Evidence_Class,
        Suggested_Manuscript_Use,
        Raw_Prevalence,
        AST_SD,
        AST_IQR
    )

write_csv(
    manuscript_candidate_table,
    file.path(table_dir, "Overall_CRC_DA_Feature_Survival_Manuscript_Candidates.csv")
)

manuscript_compact_table <- manuscript_candidate_table %>%
    transmute(
        Taxon = Feature_Name,
        `DA source` = Ob_DA_Class,
        `DA direction` = DA_Direction,
        `DA q` = fmt_p(Ob_CRC_Q),
        `Continuous AST HR per SD (95% CI)` = Continuous_HR_Per_SD_AST_CI,
        `Continuous P` = fmt_p(Continuous_P),
        `Continuous FDR` = fmt_p(Continuous_FDR),
        `Median high vs low HR (95% CI)` = Median_HR_CI,
        `Median P` = fmt_p(Median_P),
        `Presence vs absence HR (95% CI)` = PA_HR_CI,
        `Presence/absence P` = fmt_p(PA_P),
        `Nominal models` = Significant_Models,
        `Interpretation tier` = Evidence_Class
    )

write_csv(
    manuscript_compact_table,
    file.path(table_dir, "Overall_CRC_DA_Feature_Survival_Manuscript_Compact_Table.csv")
)

message("=== 5. Build manuscript-style survival figures ===")

forest_plot_data <- integrated_survival %>%
    filter(Any_Nominal, !is.na(Continuous_HR_Per_SD_AST)) %>%
    mutate(
        DA_Direction_Short = case_when(
            Ob_CRC_Coef > 0 ~ "CRC-enriched",
            Ob_CRC_Coef < 0 ~ "CRC-depleted",
            TRUE ~ "Unknown"
        ),
        Continuous_Status = ifelse(Continuous_Nominal, "Continuous P<0.05", "Continuous ns"),
        Plot_Label = factor(Feature_Name, levels = rev(Feature_Name[order(Continuous_HR_Per_SD_AST)]))
    )

write_tsv(
    forest_plot_data,
    file.path(figure_data_dir, "Continuous_AST_Forest_Manuscript_Candidates.tsv")
)

p_forest <- ggplot(
    forest_plot_data,
    aes(x = Continuous_HR_Per_SD_AST, y = Plot_Label)
) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, color = "grey45") +
    geom_errorbar(
        aes(xmin = Continuous_HR_Per_SD_AST_Lower, xmax = Continuous_HR_Per_SD_AST_Upper),
        orientation = "y",
        width = 0.18,
        linewidth = 0.45,
        color = "grey30"
    ) +
    geom_point(aes(color = DA_Direction_Short, shape = Continuous_Status), size = 2.4, stroke = 0.9) +
    scale_x_log10() +
    scale_color_manual(
        values = c("CRC-depleted" = "#3B6EA8", "CRC-enriched" = "#C44E52", "Unknown" = "grey45")
    ) +
    scale_shape_manual(values = c("Continuous P<0.05" = 16, "Continuous ns" = 1)) +
    labs(
        x = "Hazard ratio per 1 SD AST-transformed abundance",
        y = NULL,
        color = "DA direction",
        shape = "Continuous model",
        title = "Stage-adjusted continuous AST Cox model"
    ) +
    theme_classic(base_size = 12) +
    theme(
        plot.title = element_text(size = 12, face = "plain"),
        axis.text.y = element_text(size = 9),
        legend.position = "right"
    )

ggsave(
    file.path(figure_dir, "Continuous_AST_Forest_Manuscript_Candidates.pdf"),
    p_forest,
    width = 6.9,
    height = 4.7
)
ggsave(
    file.path(figure_dir, "Continuous_AST_Forest_Manuscript_Candidates.png"),
    p_forest,
    width = 6.9,
    height = 4.7,
    dpi = 300
)

evidence_plot_data <- integrated_survival %>%
    filter(Any_Nominal) %>%
    transmute(
        Feature,
        Feature_Name,
        Feature_Order = match(Feature_Name, levels(forest_plot_data$Plot_Label)),
        DA_Direction,
        `Continuous AST` = Continuous_HR_Per_SD_AST,
        `Median split` = Median_HR,
        `Presence/absence` = PA_HR,
        Continuous_AST_P = Continuous_P,
        Median_split_P = Median_P,
        Presence_absence_P = PA_P
    ) %>%
    pivot_longer(
        cols = c(`Continuous AST`, `Median split`, `Presence/absence`),
        names_to = "Model",
        values_to = "HR"
    ) %>%
    mutate(
        P_Value = case_when(
            Model == "Continuous AST" ~ Continuous_AST_P,
            Model == "Median split" ~ Median_split_P,
            Model == "Presence/absence" ~ Presence_absence_P,
            TRUE ~ NA_real_
        ),
        Model = factor(Model, levels = c("Continuous AST", "Median split", "Presence/absence")),
        Feature_Name = factor(
            Feature_Name,
            levels = rev(levels(forest_plot_data$Plot_Label))
        ),
        Signed_NegLog10_P = case_when(
            is.na(P_Value) | is.na(HR) ~ NA_real_,
            HR >= 1 ~ pmin(-log10(P_Value), 3),
            HR < 1 ~ -pmin(-log10(P_Value), 3)
        ),
        Sig_Label = ifelse(!is.na(P_Value) & P_Value < 0.05, "*", "")
    ) %>%
    select(
        Feature,
        Feature_Name,
        DA_Direction,
        Model,
        HR,
        P_Value,
        Signed_NegLog10_P,
        Sig_Label
    )

write_tsv(
    evidence_plot_data,
    file.path(figure_data_dir, "Survival_Model_Evidence_Heatmap_Manuscript_Candidates.tsv")
)

p_evidence <- ggplot(evidence_plot_data, aes(x = Model, y = Feature_Name, fill = Signed_NegLog10_P)) +
    geom_tile(color = "white", linewidth = 0.45) +
    geom_text(aes(label = Sig_Label), size = 4.1, color = "black") +
    scale_fill_gradient2(
        low = "#3B6EA8",
        mid = "white",
        high = "#C44E52",
        midpoint = 0,
        limits = c(-3, 3),
        na.value = "grey92",
        name = "Signed -log10(P)"
    ) +
    labs(
        x = NULL,
        y = NULL,
        title = "Survival evidence across three stage-adjusted models"
    ) +
    theme_classic(base_size = 12) +
    theme(
        plot.title = element_text(size = 12, face = "plain"),
        axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = element_text(size = 9),
        legend.position = "right"
    )

ggsave(
    file.path(figure_dir, "Survival_Model_Evidence_Heatmap_Manuscript_Candidates.pdf"),
    p_evidence,
    width = 6.4,
    height = 4.5
)
ggsave(
    file.path(figure_dir, "Survival_Model_Evidence_Heatmap_Manuscript_Candidates.png"),
    p_evidence,
    width = 6.4,
    height = 4.5,
    dpi = 300
)

mode_level <- function(x) {
    tab <- sort(table(x), decreasing = TRUE)
    names(tab)[1]
}

inverse_ast_percent <- function(x) {
    (sin(x)^2) * 100
}

representative_ast_values <- function(ast_values) {
    q_vals <- quantile(ast_values, probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE)
    low_ast <- q_vals[1]
    high_ast <- q_vals[2]
    method <- "25th vs 75th percentile"

    if (is.na(high_ast) || high_ast <= low_ast) {
        positive_ast <- ast_values[ast_values > 0]
        low_ast <- 0
        if (length(positive_ast) > 0) {
            high_ast <- quantile(positive_ast, probs = 0.75, na.rm = TRUE, names = FALSE)
            method <- "0 vs positive 75th percentile"
        } else {
            high_ast <- max(ast_values, na.rm = TRUE)
            method <- "0 vs maximum"
        }
    }

    c(low_ast = low_ast, high_ast = high_ast, method = method)
}

predict_adjusted_survival_curves <- function(feature, meta_row) {
    val <- as.numeric(abd_final[, feature])
    dat <- surv_base %>%
        mutate(
            Bacteria_Val = val,
            Bacteria_AST = transform_ast_percent(val)
        ) %>%
        drop_na(Stage)

    fit <- coxph(
        Surv(surv.time..Y., status) ~ Bacteria_AST + Age + Gender + BMI_Group + Stage,
        data = dat,
        control = coxph.control(iter.max = 50)
    )

    ast_pair <- representative_ast_values(dat$Bacteria_AST)
    low_ast <- as.numeric(ast_pair["low_ast"])
    high_ast <- as.numeric(ast_pair["high_ast"])
    selection_method <- unname(ast_pair["method"])

    base_newdata <- data.frame(
        Bacteria_AST = c(low_ast, high_ast),
        Age = median(dat$Age, na.rm = TRUE),
        Gender = mode_level(dat$Gender),
        BMI_Group = mode_level(dat$BMI_Group),
        Stage = mode_level(dat$Stage),
        stringsAsFactors = FALSE
    )
    base_newdata$Gender <- factor(base_newdata$Gender, levels = levels(dat$Gender))
    base_newdata$BMI_Group <- factor(base_newdata$BMI_Group, levels = levels(dat$BMI_Group))
    base_newdata$Stage <- factor(base_newdata$Stage, levels = levels(dat$Stage))

    make_curve <- function(row_id, curve_role) {
        sf <- survfit(fit, newdata = base_newdata[row_id, , drop = FALSE])
        s <- summary(sf)
        lower <- if (is.null(s$lower)) rep(NA_real_, length(s$time)) else s$lower
        upper <- if (is.null(s$upper)) rep(NA_real_, length(s$time)) else s$upper

        tibble(
            Feature = feature,
            Feature_Name = meta_row$Feature_Name,
            DA_Direction = meta_row$DA_Direction,
            Curve_Role = curve_role,
            time = c(0, s$time),
            survival = c(1, s$surv),
            lower = c(1, lower),
            upper = c(1, upper),
            AST_Value = base_newdata$Bacteria_AST[row_id],
            Raw_Abundance_Percent = inverse_ast_percent(base_newdata$Bacteria_AST[row_id]),
            Representative_Value_Method = selection_method,
            HR_Per_SD_AST = meta_row$Continuous_HR_Per_SD_AST,
            HR_Per_SD_AST_Lower = meta_row$Continuous_HR_Per_SD_AST_Lower,
            HR_Per_SD_AST_Upper = meta_row$Continuous_HR_Per_SD_AST_Upper,
            P_Value = meta_row$Continuous_P,
            P_Adj_FDR = meta_row$Continuous_FDR,
            Reference_Age = base_newdata$Age[row_id],
            Reference_Gender = as.character(base_newdata$Gender[row_id]),
            Reference_BMI_Group = as.character(base_newdata$BMI_Group[row_id]),
            Reference_Stage = as.character(base_newdata$Stage[row_id])
        )
    }

    bind_rows(
        make_curve(1, "Low"),
        make_curve(2, "High")
    )
}

adjusted_curve_features <- integrated_survival %>%
    filter(Continuous_Nominal) %>%
    arrange(Ob_CRC_Coef > 0, Continuous_P)

adjusted_curve_data <- bind_rows(lapply(seq_len(nrow(adjusted_curve_features)), function(i) {
    predict_adjusted_survival_curves(adjusted_curve_features$Feature[i], adjusted_curve_features[i, ])
})) %>%
    mutate(
        Curve_Display = case_when(
            Curve_Role == "Low" ~ "Low reference",
            HR_Per_SD_AST >= 1 ~ "High abundance (HR>1)",
            HR_Per_SD_AST < 1 ~ "High abundance (HR<1)",
            TRUE ~ "High abundance"
        ),
        Plot_Label = paste0(
            Feature_Name,
            " | HR=", sprintf("%.2f", HR_Per_SD_AST),
            ", P=", fmt_p(P_Value)
        ),
        Plot_Label = factor(Plot_Label, levels = unique(Plot_Label)),
        Curve_ID = paste(Feature, Curve_Role, sep = "__")
    )

write_tsv(
    adjusted_curve_data,
    file.path(figure_data_dir, "Continuous_AST_Adjusted_Survival_Selected_Candidates.tsv")
)

adjusted_label_map <- adjusted_curve_data %>%
    distinct(Plot_Label, Feature_Name, HR_Per_SD_AST, P_Value) %>%
    mutate(
        Feature_Label = vapply(Feature_Name, function(x) paste(strwrap(x, width = 26), collapse = "\n"), character(1)),
        Facet_Label = paste0(
            Feature_Label,
            "\nHR=", sprintf("%.2f", HR_Per_SD_AST),
            ", P=", fmt_p(P_Value)
        )
    )
adjusted_label_map <- setNames(adjusted_label_map$Facet_Label, as.character(adjusted_label_map$Plot_Label))

p_adjusted <- ggplot(
    adjusted_curve_data,
    aes(x = time, y = survival, group = Curve_ID, color = Curve_Display, linetype = Curve_Display)
) +
    geom_step(linewidth = 0.75) +
    facet_wrap(~Plot_Label, ncol = 3, labeller = as_labeller(adjusted_label_map)) +
    scale_color_manual(
        values = c(
            "Low reference" = "grey45",
            "High abundance (HR<1)" = "#3B6EA8",
            "High abundance (HR>1)" = "#C44E52",
            "High abundance" = "#C44E52"
        )
    ) +
    scale_linetype_manual(
        values = c(
            "Low reference" = "dashed",
            "High abundance (HR<1)" = "solid",
            "High abundance (HR>1)" = "solid",
            "High abundance" = "solid"
        )
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
        x = "Time (years)",
        y = "Predicted overall survival probability",
        color = NULL,
        linetype = NULL,
        title = "Adjusted survival curves for continuous AST-significant taxa"
    ) +
    theme_classic(base_size = 12) +
    theme(
        plot.title = element_text(size = 12, face = "plain"),
        strip.background = element_blank(),
        strip.text = element_text(size = 9),
        legend.position = "bottom"
    )

ggsave(
    file.path(figure_dir, "Continuous_AST_Adjusted_Survival_Selected_Candidates.pdf"),
    p_adjusted,
    width = 7.2,
    height = 6.0
)
ggsave(
    file.path(figure_dir, "Continuous_AST_Adjusted_Survival_Selected_Candidates.png"),
    p_adjusted,
    width = 7.2,
    height = 6.0,
    dpi = 300
)

best_grouped_candidates <- integrated_survival %>%
    filter(Any_Grouped_Nominal) %>%
    arrange(Best_Grouped_P)

build_best_grouped_km <- function(feature, meta_row) {
    val <- as.numeric(abd_final[, feature])
    dat <- surv_base %>%
        mutate(Bacteria_Val = val)

    if (meta_row$Best_Grouped_Model == "Median") {
        med_val <- median(val, na.rm = TRUE)
        dat <- dat %>%
            mutate(Bacteria_Group = factor(ifelse(Bacteria_Val >= med_val, "High", "Low"), levels = c("Low", "High")))
    } else {
        dat <- dat %>%
            mutate(Bacteria_Group = factor(ifelse(Bacteria_Val > 0, "Present", "Absent"), levels = c("Absent", "Present")))
    }

    dat <- dat %>% drop_na(Stage)
    km_fit <- survfit(Surv(surv.time..Y., status) ~ Bacteria_Group, data = dat)
    s <- summary(km_fit)
    strata <- if (is.null(s$strata)) rep(levels(dat$Bacteria_Group)[1], length(s$time)) else as.character(s$strata)
    strata <- sub("^Bacteria_Group=", "", strata)

    tibble(
        Feature = feature,
        Feature_Name = meta_row$Feature_Name,
        DA_Direction = meta_row$DA_Direction,
        Best_Grouped_Model = meta_row$Best_Grouped_Model,
        Best_Grouped_P = meta_row$Best_Grouped_P,
        Best_Grouped_HR = meta_row$Best_Grouped_HR,
        Best_Grouped_HR_CI = meta_row$Best_Grouped_HR_CI,
        Best_Grouped_Direction_Check = if ("Best_Grouped_Direction_Check" %in% names(meta_row)) {
            meta_row$Best_Grouped_Direction_Check
        } else {
            NA_character_
        },
        strata = strata,
        time = s$time,
        n_risk = s$n.risk,
        n_event = s$n.event,
        n_censor = s$n.censor,
        survival = s$surv,
        lower = s$lower,
        upper = s$upper,
        Group_Counts = if (meta_row$Best_Grouped_Model == "Median") meta_row$Median_Group_Counts else meta_row$PA_Group_Counts
    )
}

best_grouped_km_data <- bind_rows(lapply(seq_len(nrow(best_grouped_candidates)), function(i) {
    build_best_grouped_km(best_grouped_candidates$Feature[i], best_grouped_candidates[i, ])
})) %>%
    group_by(Feature, strata) %>%
    arrange(time, .by_group = TRUE) %>%
    group_modify(~bind_rows(
        .x[1, ] %>% mutate(time = 0, n_event = 0, n_censor = 0, survival = 1, lower = 1, upper = 1),
        .x
    )) %>%
    ungroup() %>%
    mutate(
        Best_Grouped_Model_Display = ifelse(Best_Grouped_Model == "Presence_absence", "PA", Best_Grouped_Model),
        Plot_Label = paste0(
            Feature_Name,
            " | ", Best_Grouped_Model_Display,
            " HR=", sprintf("%.2f", Best_Grouped_HR),
            ", P=", fmt_p(Best_Grouped_P)
        ),
        Plot_Label = factor(Plot_Label, levels = unique(Plot_Label)),
        Curve_ID = paste(Feature, strata, sep = "__")
    )

write_tsv(
    best_grouped_km_data,
    file.path(figure_data_dir, "KM_Best_Grouped_Manuscript_Candidates.tsv")
)

best_km_label_map <- best_grouped_km_data %>%
    distinct(Plot_Label, Feature_Name, Best_Grouped_Model_Display, Best_Grouped_HR, Best_Grouped_P) %>%
    mutate(
        Feature_Label = vapply(Feature_Name, function(x) paste(strwrap(x, width = 28), collapse = "\n"), character(1)),
        Facet_Label = paste0(
            Feature_Label,
            "\n", Best_Grouped_Model_Display,
            " HR=", sprintf("%.2f", Best_Grouped_HR),
            ", P=", fmt_p(Best_Grouped_P)
        )
    )
best_km_label_map <- setNames(best_km_label_map$Facet_Label, as.character(best_km_label_map$Plot_Label))

p_best_km <- ggplot(best_grouped_km_data, aes(x = time, y = survival, group = Curve_ID, color = strata)) +
    geom_step(linewidth = 0.75) +
    facet_wrap(~Plot_Label, ncol = 3, labeller = as_labeller(best_km_label_map)) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_color_manual(values = c("Absent" = "grey55", "Present" = "#D39200", "Low" = "#6A8CAF", "High" = "#D94F45")) +
    labs(
        x = "Time (years)",
        y = "Overall survival probability",
        color = NULL,
        title = "Kaplan-Meier curves for best grouped survival candidates"
    ) +
    theme_classic(base_size = 12) +
    theme(
        plot.title = element_text(size = 12, face = "plain"),
        strip.background = element_blank(),
        strip.text = element_text(size = 8.5),
        legend.position = "bottom"
    )

ggsave(
    file.path(figure_dir, "KM_Best_Grouped_Manuscript_Candidates.pdf"),
    p_best_km,
    width = 8.0,
    height = 6.4
)
ggsave(
    file.path(figure_dir, "KM_Best_Grouped_Manuscript_Candidates.png"),
    p_best_km,
    width = 8.0,
    height = 6.4,
    dpi = 300
)

message("=== 6. Build Ob-specific-only survival summaries ===")

adjust_fdr_no_na <- function(p) {
    out <- rep(NA_real_, length(p))
    idx <- !is.na(p)
    out[idx] <- p.adjust(p[idx], method = "BH")
    out
}

classify_direction <- function(da_coef, hr) {
    if (is.na(da_coef) || is.na(hr)) {
        return(NA_character_)
    }
    if (da_coef > 0 && hr > 1) {
        "Concordant adverse: CRC-enriched and higher abundance/hit has higher hazard"
    } else if (da_coef > 0 && hr < 1) {
        "Discordant protective: CRC-enriched but higher abundance/hit has lower hazard"
    } else if (da_coef < 0 && hr < 1) {
        "Concordant protective: CRC-depleted and higher abundance/hit has lower hazard"
    } else if (da_coef < 0 && hr > 1) {
        "Discordant adverse: CRC-depleted but higher abundance/hit has higher hazard"
    } else {
        "No clear direction"
    }
}

short_direction <- function(x) {
    case_when(
        grepl("^Concordant adverse", x) ~ "Concordant adverse",
        grepl("^Concordant protective", x) ~ "Concordant protective",
        grepl("^Discordant protective", x) ~ "Discordant protective",
        grepl("^Discordant adverse", x) ~ "Discordant adverse",
        TRUE ~ x
    )
}

ob_specific_survival <- integrated_survival %>%
    filter(Ob_DA_Class == "Ob-specific") %>%
    mutate(
        Continuous_FDR_ObSpecific = adjust_fdr_no_na(Continuous_P),
        Median_FDR_ObSpecific = adjust_fdr_no_na(Median_P),
        PA_FDR_ObSpecific = adjust_fdr_no_na(PA_P),
        Continuous_Direction_Check = vapply(
            seq_along(Ob_CRC_Coef),
            function(i) classify_direction(Ob_CRC_Coef[i], Continuous_HR_Per_SD_AST[i]),
            character(1)
        ),
        Median_Direction_Check = vapply(
            seq_along(Ob_CRC_Coef),
            function(i) classify_direction(Ob_CRC_Coef[i], Median_HR[i]),
            character(1)
        ),
        PA_Direction_Check = vapply(
            seq_along(Ob_CRC_Coef),
            function(i) classify_direction(Ob_CRC_Coef[i], PA_HR[i]),
            character(1)
        ),
        Continuous_Nominal_ObSpecific = !is.na(Continuous_P) & Continuous_P < 0.05,
        Median_Nominal_ObSpecific = !is.na(Median_P) & Median_P < 0.05,
        PA_Nominal_ObSpecific = !is.na(PA_P) & PA_P < 0.05,
        Continuous_FDR10_ObSpecific = !is.na(Continuous_FDR_ObSpecific) & Continuous_FDR_ObSpecific < 0.10,
        Median_FDR10_ObSpecific = !is.na(Median_FDR_ObSpecific) & Median_FDR_ObSpecific < 0.10,
        PA_FDR10_ObSpecific = !is.na(PA_FDR_ObSpecific) & PA_FDR_ObSpecific < 0.10,
        Any_Nominal_ObSpecific = Continuous_Nominal_ObSpecific | Median_Nominal_ObSpecific | PA_Nominal_ObSpecific,
        Any_FDR10_ObSpecific = Continuous_FDR10_ObSpecific | Median_FDR10_ObSpecific | PA_FDR10_ObSpecific,
        Best_Grouped_FDR_ObSpecific = case_when(
            Best_Grouped_Model == "Median" ~ Median_FDR_ObSpecific,
            Best_Grouped_Model == "Presence_absence" ~ PA_FDR_ObSpecific,
            TRUE ~ NA_real_
        ),
        Best_Grouped_Direction_Check = case_when(
            Best_Grouped_Model == "Median" ~ Median_Direction_Check,
            Best_Grouped_Model == "Presence_absence" ~ PA_Direction_Check,
            TRUE ~ NA_character_
        ),
        ObSpecific_Evidence_Class = case_when(
            Any_FDR10_ObSpecific & Continuous_Nominal_ObSpecific & (Median_Nominal_ObSpecific | PA_Nominal_ObSpecific) ~
                "Ob-specific FDR-supported plus grouped support",
            Any_FDR10_ObSpecific ~ "Ob-specific FDR-supported",
            Continuous_Nominal_ObSpecific & (Median_Nominal_ObSpecific | PA_Nominal_ObSpecific) ~
                "Ob-specific continuous and grouped nominal",
            Continuous_Nominal_ObSpecific ~ "Ob-specific continuous-only nominal",
            Median_Nominal_ObSpecific | PA_Nominal_ObSpecific ~ "Ob-specific grouped-only nominal",
            TRUE ~ "No Ob-specific nominal survival association"
        )
    ) %>%
    rowwise() %>%
    mutate(
        Significant_Models_ObSpecific = paste(
            c(
                if (Continuous_Nominal_ObSpecific) "Continuous_AST" else NULL,
                if (Median_Nominal_ObSpecific) "Median" else NULL,
                if (PA_Nominal_ObSpecific) "Presence_absence" else NULL
            ),
            collapse = "; "
        ),
        Significant_Models_ObSpecific = ifelse(Significant_Models_ObSpecific == "", NA_character_, Significant_Models_ObSpecific),
        FDR10_Models_ObSpecific = paste(
            c(
                if (Continuous_FDR10_ObSpecific) "Continuous_AST" else NULL,
                if (Median_FDR10_ObSpecific) "Median" else NULL,
                if (PA_FDR10_ObSpecific) "Presence_absence" else NULL
            ),
            collapse = "; "
        ),
        FDR10_Models_ObSpecific = ifelse(FDR10_Models_ObSpecific == "", NA_character_, FDR10_Models_ObSpecific)
    ) %>%
    ungroup() %>%
    arrange(
        desc(Any_FDR10_ObSpecific),
        desc(Continuous_FDR10_ObSpecific),
        coalesce(Continuous_P, Best_Grouped_P, 1)
    )

write_csv(
    ob_specific_survival,
    file.path(table_dir, "Overall_CRC_DA_Feature_Survival_ObSpecific_Integrated.csv")
)

ob_specific_candidates <- ob_specific_survival %>%
    filter(Any_Nominal_ObSpecific) %>%
    transmute(
        Feature_Name,
        DA_Direction,
        Ob_CRC_Coef,
        Ob_CRC_Q,
        Continuous_HR_Per_SD_AST_CI,
        Continuous_P,
        Continuous_FDR_ObSpecific,
        Continuous_Direction_Check,
        Median_HR_CI,
        Median_P,
        Median_FDR_ObSpecific,
        Median_Direction_Check,
        PA_HR_CI,
        PA_P,
        PA_FDR_ObSpecific,
        PA_Direction_Check,
        Significant_Models_ObSpecific,
        FDR10_Models_ObSpecific,
        ObSpecific_Evidence_Class,
        Raw_Prevalence,
        AST_SD,
        AST_IQR
    )

write_csv(
    ob_specific_candidates,
    file.path(table_dir, "Overall_CRC_DA_Feature_Survival_ObSpecific_Candidates.csv")
)

ob_specific_direction_audit <- ob_specific_survival %>%
    select(
        Feature,
        Feature_Name,
        DA_Direction,
        Ob_CRC_Coef,
        Ob_CRC_Q,
        Continuous_HR_Per_SD_AST,
        Continuous_P,
        Continuous_FDR_ObSpecific,
        Continuous_Direction_Check,
        Median_HR,
        Median_P,
        Median_FDR_ObSpecific,
        Median_Direction_Check,
        PA_HR,
        PA_P,
        PA_FDR_ObSpecific,
        PA_Direction_Check,
        Significant_Models_ObSpecific,
        FDR10_Models_ObSpecific,
        ObSpecific_Evidence_Class
    )

write_csv(
    ob_specific_direction_audit,
    file.path(table_dir, "Overall_CRC_DA_Feature_Survival_ObSpecific_Direction_Audit.csv")
)

ob_cont_fdr10 <- ob_specific_survival %>%
    filter(Continuous_FDR10_ObSpecific) %>%
    mutate(
        Continuous_FDR = Continuous_FDR_ObSpecific,
        DA_Direction_Short = case_when(
            Ob_CRC_Coef > 0 ~ "CRC-enriched",
            Ob_CRC_Coef < 0 ~ "CRC-depleted",
            TRUE ~ "Unknown"
        ),
        Direction_Check_Short = short_direction(Continuous_Direction_Check),
        Plot_Label = factor(Feature_Name, levels = rev(Feature_Name[order(Continuous_HR_Per_SD_AST)]))
    )

if (nrow(ob_cont_fdr10) > 0) {
    write_tsv(
        ob_cont_fdr10,
        file.path(figure_data_dir, "ObSpecific_Continuous_AST_FDR10_Forest.tsv")
    )

    p_ob_forest <- ggplot(ob_cont_fdr10, aes(x = Continuous_HR_Per_SD_AST, y = Plot_Label)) +
        geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, color = "grey45") +
        geom_errorbar(
            aes(xmin = Continuous_HR_Per_SD_AST_Lower, xmax = Continuous_HR_Per_SD_AST_Upper),
            orientation = "y",
            width = 0.18,
            linewidth = 0.45,
            color = "grey30"
        ) +
        geom_point(aes(color = DA_Direction_Short), size = 2.8) +
        scale_x_log10() +
        scale_color_manual(
            values = c("CRC-depleted" = "#3B6EA8", "CRC-enriched" = "#C44E52", "Unknown" = "grey45")
        ) +
        labs(
            x = "Hazard ratio per 1 SD AST-transformed abundance",
            y = NULL,
            color = "Obese DA direction",
            title = "Ob-specific taxa with continuous AST FDR < 0.10"
        ) +
        theme_classic(base_size = 12) +
        theme(
            plot.title = element_text(size = 12, face = "plain"),
            axis.text.y = element_text(size = 10),
            legend.position = "right"
        )

    ggsave(
        file.path(figure_dir, "ObSpecific_Continuous_AST_FDR10_Forest.pdf"),
        p_ob_forest,
        width = 6.2,
        height = 2.8
    )
    ggsave(
        file.path(figure_dir, "ObSpecific_Continuous_AST_FDR10_Forest.png"),
        p_ob_forest,
        width = 6.2,
        height = 2.8,
        dpi = 300
    )

    ob_adjusted_curve_data <- bind_rows(lapply(seq_len(nrow(ob_cont_fdr10)), function(i) {
        predict_adjusted_survival_curves(ob_cont_fdr10$Feature[i], ob_cont_fdr10[i, ])
    })) %>%
        mutate(
            Curve_Display = case_when(
                Curve_Role == "Low" ~ "Low reference",
                HR_Per_SD_AST >= 1 ~ "High abundance (HR>1)",
                HR_Per_SD_AST < 1 ~ "High abundance (HR<1)",
                TRUE ~ "High abundance"
            ),
            Plot_Label = paste0(
                Feature_Name,
                " | HR=", sprintf("%.2f", HR_Per_SD_AST),
                ", FDR=", fmt_p(P_Adj_FDR)
            ),
            Plot_Label = factor(Plot_Label, levels = unique(Plot_Label)),
            Curve_ID = paste(Feature, Curve_Role, sep = "__")
        )

    write_tsv(
        ob_adjusted_curve_data,
        file.path(figure_data_dir, "ObSpecific_Continuous_AST_FDR10_Adjusted_Survival.tsv")
    )

    ob_adjusted_label_map <- ob_adjusted_curve_data %>%
        distinct(Plot_Label, Feature_Name, HR_Per_SD_AST, P_Value, P_Adj_FDR) %>%
        mutate(
            Feature_Label = vapply(Feature_Name, function(x) paste(strwrap(x, width = 26), collapse = "\n"), character(1)),
            Facet_Label = paste0(
                Feature_Label,
                "\nHR=", sprintf("%.2f", HR_Per_SD_AST),
                ", P=", fmt_p(P_Value),
                ", FDR=", fmt_p(P_Adj_FDR)
            )
        )
    ob_adjusted_label_map <- setNames(ob_adjusted_label_map$Facet_Label, as.character(ob_adjusted_label_map$Plot_Label))

    p_ob_adjusted <- ggplot(
        ob_adjusted_curve_data,
        aes(x = time, y = survival, group = Curve_ID, color = Curve_Display, linetype = Curve_Display)
    ) +
        geom_step(linewidth = 0.75) +
        facet_wrap(~Plot_Label, ncol = 2, labeller = as_labeller(ob_adjusted_label_map)) +
        scale_color_manual(
            values = c(
                "Low reference" = "grey45",
                "High abundance (HR<1)" = "#3B6EA8",
                "High abundance (HR>1)" = "#C44E52",
                "High abundance" = "#C44E52"
            )
        ) +
        scale_linetype_manual(
            values = c(
                "Low reference" = "dashed",
                "High abundance (HR<1)" = "solid",
                "High abundance (HR>1)" = "solid",
                "High abundance" = "solid"
            )
        ) +
        coord_cartesian(ylim = c(0, 1)) +
        labs(
            x = "Time (years)",
            y = "Predicted overall survival probability",
            color = NULL,
            linetype = NULL,
            title = "Adjusted survival curves for Ob-specific continuous AST FDR taxa"
        ) +
        theme_classic(base_size = 12) +
        theme(
            plot.title = element_text(size = 12, face = "plain"),
            strip.background = element_blank(),
            strip.text = element_text(size = 9),
            legend.position = "bottom"
        )

    ggsave(
        file.path(figure_dir, "ObSpecific_Continuous_AST_FDR10_Adjusted_Survival.pdf"),
        p_ob_adjusted,
        width = 6.8,
        height = 5.6
    )
    ggsave(
        file.path(figure_dir, "ObSpecific_Continuous_AST_FDR10_Adjusted_Survival.png"),
        p_ob_adjusted,
        width = 6.8,
        height = 5.6,
        dpi = 300
    )
}

ob_specific_evidence_data <- ob_specific_survival %>%
    filter(Any_Nominal_ObSpecific) %>%
    transmute(
        Feature,
        Feature_Name,
        DA_Direction,
        `Continuous AST` = Continuous_HR_Per_SD_AST,
        `Median split` = Median_HR,
        `Presence/absence` = PA_HR,
        Continuous_AST_P = Continuous_P,
        Median_split_P = Median_P,
        Presence_absence_P = PA_P,
        Continuous_AST_FDR = Continuous_FDR_ObSpecific,
        Median_split_FDR = Median_FDR_ObSpecific,
        Presence_absence_FDR = PA_FDR_ObSpecific
    ) %>%
    pivot_longer(
        cols = c(`Continuous AST`, `Median split`, `Presence/absence`),
        names_to = "Model",
        values_to = "HR"
    ) %>%
    mutate(
        P_Value = case_when(
            Model == "Continuous AST" ~ Continuous_AST_P,
            Model == "Median split" ~ Median_split_P,
            Model == "Presence/absence" ~ Presence_absence_P,
            TRUE ~ NA_real_
        ),
        FDR_ObSpecific = case_when(
            Model == "Continuous AST" ~ Continuous_AST_FDR,
            Model == "Median split" ~ Median_split_FDR,
            Model == "Presence/absence" ~ Presence_absence_FDR,
            TRUE ~ NA_real_
        ),
        Model = factor(Model, levels = c("Continuous AST", "Median split", "Presence/absence")),
        Feature_Name = factor(Feature_Name, levels = rev(unique(Feature_Name))),
        Nominal = !is.na(P_Value) & P_Value < 0.05,
        FDR10 = !is.na(FDR_ObSpecific) & FDR_ObSpecific < 0.10,
        Signed_NegLog10_P = case_when(
            !Nominal ~ NA_real_,
            HR >= 1 ~ pmin(-log10(P_Value), 3),
            HR < 1 ~ -pmin(-log10(P_Value), 3)
        ),
        Sig_Label = case_when(
            FDR10 ~ "\u2020",
            Nominal ~ "*",
            TRUE ~ ""
        )
    ) %>%
    select(Feature, Feature_Name, DA_Direction, Model, HR, P_Value, FDR_ObSpecific, Nominal, FDR10, Signed_NegLog10_P, Sig_Label)

write_tsv(
    ob_specific_evidence_data,
    file.path(figure_data_dir, "ObSpecific_Significant_Model_Evidence_Heatmap.tsv")
)

p_ob_evidence <- ggplot(ob_specific_evidence_data, aes(x = Model, y = Feature_Name, fill = Signed_NegLog10_P)) +
    geom_tile(color = "white", linewidth = 0.45) +
    geom_text(aes(label = Sig_Label), size = 4.1, color = "black") +
    scale_fill_gradient2(
        low = "#3B6EA8",
        mid = "white",
        high = "#C44E52",
        midpoint = 0,
        limits = c(-3, 3),
        na.value = "grey94",
        name = "Signed -log10(P)"
    ) +
    labs(
        x = NULL,
        y = NULL,
        title = "Ob-specific taxa with nominal survival associations"
    ) +
    theme_classic(base_size = 12) +
    theme(
        plot.title = element_text(size = 12, face = "plain"),
        axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = element_text(size = 9),
        legend.position = "right"
    )

ggsave(
    file.path(figure_dir, "ObSpecific_Significant_Model_Evidence_Heatmap.pdf"),
    p_ob_evidence,
    width = 6.2,
    height = 3.7
)
ggsave(
    file.path(figure_dir, "ObSpecific_Significant_Model_Evidence_Heatmap.png"),
    p_ob_evidence,
    width = 6.2,
    height = 3.7,
    dpi = 300
)

ob_grouped_nominal <- ob_specific_survival %>%
    filter(Median_Nominal_ObSpecific | PA_Nominal_ObSpecific) %>%
    arrange(Best_Grouped_P)

if (nrow(ob_grouped_nominal) > 0) {
    ob_grouped_km_data <- bind_rows(lapply(seq_len(nrow(ob_grouped_nominal)), function(i) {
        build_best_grouped_km(ob_grouped_nominal$Feature[i], ob_grouped_nominal[i, ])
    })) %>%
        group_by(Feature, strata) %>%
        arrange(time, .by_group = TRUE) %>%
        group_modify(~bind_rows(
            .x[1, ] %>% mutate(time = 0, n_event = 0, n_censor = 0, survival = 1, lower = 1, upper = 1),
            .x
        )) %>%
        ungroup() %>%
        mutate(
            Best_Grouped_Model_Display = ifelse(Best_Grouped_Model == "Presence_absence", "PA", Best_Grouped_Model),
            Direction_Check_Short = short_direction(Best_Grouped_Direction_Check),
            Plot_Label = paste0(
                Feature_Name,
                " | ", Best_Grouped_Model_Display,
                " HR=", sprintf("%.2f", Best_Grouped_HR),
                ", P=", fmt_p(Best_Grouped_P)
            ),
            Plot_Label = factor(Plot_Label, levels = unique(Plot_Label)),
            Curve_ID = paste(Feature, strata, sep = "__")
        )

    write_tsv(
        ob_grouped_km_data,
        file.path(figure_data_dir, "ObSpecific_KM_Grouped_Nominal_Candidates.tsv")
    )

    ob_grouped_label_map <- ob_grouped_km_data %>%
        distinct(Plot_Label, Feature_Name, Best_Grouped_Model_Display, Best_Grouped_HR, Best_Grouped_P, Direction_Check_Short) %>%
        mutate(
            Feature_Label = vapply(Feature_Name, function(x) paste(strwrap(x, width = 28), collapse = "\n"), character(1)),
            Direction_Label = ifelse(grepl("^Discordant", Direction_Check_Short), "discordant", "consistent"),
            Facet_Label = paste0(
                Feature_Label,
                "\n", Best_Grouped_Model_Display,
                " HR=", sprintf("%.2f", Best_Grouped_HR),
                ", P=", fmt_p(Best_Grouped_P),
                "\n", Direction_Label
            )
        )
    ob_grouped_label_map <- setNames(ob_grouped_label_map$Facet_Label, as.character(ob_grouped_label_map$Plot_Label))

    p_ob_km <- ggplot(ob_grouped_km_data, aes(x = time, y = survival, group = Curve_ID, color = strata)) +
        geom_step(linewidth = 0.75) +
        facet_wrap(~Plot_Label, ncol = 3, labeller = as_labeller(ob_grouped_label_map)) +
        coord_cartesian(ylim = c(0, 1)) +
        scale_color_manual(values = c("Absent" = "grey55", "Present" = "#D39200", "Low" = "#6A8CAF", "High" = "#D94F45")) +
        labs(
            x = "Time (years)",
            y = "Overall survival probability",
            color = NULL,
            title = "Ob-specific grouped nominal KM candidates"
        ) +
        theme_classic(base_size = 12) +
        theme(
            plot.title = element_text(size = 12, face = "plain"),
            strip.background = element_blank(),
            strip.text = element_text(size = 8.2),
            legend.position = "bottom"
        )

    ggsave(
        file.path(figure_dir, "ObSpecific_KM_Grouped_Nominal_Candidates.pdf"),
        p_ob_km,
        width = 7.5,
        height = 4.8
    )
    ggsave(
        file.path(figure_dir, "ObSpecific_KM_Grouped_Nominal_Candidates.png"),
        p_ob_km,
        width = 7.5,
        height = 4.8,
        dpi = 300
    )
}

message("=== 7. Build exploratory optimal-cutpoint KM for selected Ob-specific continuous taxa ===")

optimal_cutpoint_features <- c(
    "s__Lactobacillus_gasseri|t__SGB7038",
    "s__GGB3034_SGB4030|t__SGB4030",
    "s__Blautia_glucerasea|t__SGB4816",
    "s__Enterocloster_lavalensis|t__SGB4725"
)

find_exploratory_cutpoint <- function(dat, min_group_n = 10, min_group_events = 2) {
    candidate_cuts <- sort(unique(dat$Bacteria_AST))
    candidate_cuts <- candidate_cuts[candidate_cuts < max(dat$Bacteria_AST, na.rm = TRUE)]

    scan_res <- lapply(candidate_cuts, function(cut_ast) {
        group <- factor(ifelse(dat$Bacteria_AST > cut_ast, "High", "Low"), levels = c("Low", "High"))
        group_tab <- dat %>%
            mutate(Cutpoint_Group = group) %>%
            group_by(Cutpoint_Group) %>%
            summarise(
                N_Group = n(),
                N_Event_Group = sum(status == 1, na.rm = TRUE),
                .groups = "drop"
            )

        if (
            nrow(group_tab) != 2 ||
                min(group_tab$N_Group) < min_group_n ||
                min(group_tab$N_Event_Group) < min_group_events
        ) {
            return(NULL)
        }

        logrank_p <- tryCatch({
            sdif <- survdiff(Surv(surv.time..Y., status) ~ group, data = dat)
            pchisq(sdif$chisq, length(sdif$n) - 1, lower.tail = FALSE)
        }, error = function(e) NA_real_)

        if (is.na(logrank_p)) {
            return(NULL)
        }

        tibble(
            Cutpoint_AST = cut_ast,
            Exploratory_Logrank_P = logrank_p,
            Low_N = group_tab$N_Group[group_tab$Cutpoint_Group == "Low"],
            High_N = group_tab$N_Group[group_tab$Cutpoint_Group == "High"],
            Low_Event = group_tab$N_Event_Group[group_tab$Cutpoint_Group == "Low"],
            High_Event = group_tab$N_Event_Group[group_tab$Cutpoint_Group == "High"]
        )
    })

    scan_df <- bind_rows(scan_res)
    if (nrow(scan_df) == 0) {
        stop("No candidate cutpoint met minimum group/event constraints.")
    }

    scan_df %>%
        arrange(Exploratory_Logrank_P, desc(High_N)) %>%
        slice(1)
}

make_optimal_cutpoint_km <- function(feature, meta_row, min_group_n = 10, min_group_events = 2) {
    val <- as.numeric(abd_final[, feature])
    dat <- surv_base %>%
        mutate(
            Bacteria_Val = val,
            Bacteria_AST = transform_ast_percent(val)
        ) %>%
        drop_na(Stage)

    cut_best <- find_exploratory_cutpoint(
        dat,
        min_group_n = min_group_n,
        min_group_events = min_group_events
    )
    cut_ast <- cut_best$Cutpoint_AST
    cut_raw <- inverse_ast_percent(cut_ast)

    dat <- dat %>%
        mutate(
            Cutpoint_Group = factor(
                ifelse(Bacteria_AST > cut_ast, "High", "Low"),
                levels = c("Low", "High")
            )
        )

    group_summary <- dat %>%
        group_by(Cutpoint_Group) %>%
        summarise(
            N_Group = n(),
            N_Event_Group = sum(status == 1, na.rm = TRUE),
            .groups = "drop"
        )

    cutpoint_group_counts <- group_summary %>%
        mutate(label = paste0(Cutpoint_Group, ":N", N_Group, "_E", N_Event_Group)) %>%
        pull(label) %>%
        paste(collapse = "; ")

    cox_cut <- coxph(
        Surv(surv.time..Y., status) ~ Cutpoint_Group + Age + Gender + BMI_Group + Stage,
        data = dat,
        control = coxph.control(iter.max = 50)
    )
    cox_sum <- summary(cox_cut)
    cut_term <- grep("^Cutpoint_Group", rownames(cox_sum$coefficients), value = TRUE)[1]
    cut_hr <- unname(cox_sum$conf.int[cut_term, "exp(coef)"])
    cut_lower <- unname(cox_sum$conf.int[cut_term, "lower .95"])
    cut_upper <- unname(cox_sum$conf.int[cut_term, "upper .95"])
    cut_p <- unname(cox_sum$coefficients[cut_term, "Pr(>|z|)"])

    logrank_p <- cut_best$Exploratory_Logrank_P

    km_fit <- survfit(Surv(surv.time..Y., status) ~ Cutpoint_Group, data = dat)
    s <- summary(km_fit)
    strata <- if (is.null(s$strata)) rep(levels(dat$Cutpoint_Group)[1], length(s$time)) else as.character(s$strata)
    strata <- sub("^Cutpoint_Group=", "", strata)

    km_data <- tibble(
        Feature = feature,
        Feature_Name = meta_row$Feature_Name,
        DA_Direction = meta_row$DA_Direction,
        strata = strata,
        time = s$time,
        n_risk = s$n.risk,
        n_event = s$n.event,
        n_censor = s$n.censor,
        survival = s$surv,
        lower = s$lower,
        upper = s$upper,
        Cutpoint_AST = cut_ast,
        Cutpoint_Raw_Abundance_Percent = cut_raw,
        Cutpoint_Group_Counts = cutpoint_group_counts,
        Exploratory_Cutpoint_Cox_HR = cut_hr,
        Exploratory_Cutpoint_Cox_HR_Lower = cut_lower,
        Exploratory_Cutpoint_Cox_HR_Upper = cut_upper,
        Exploratory_Cutpoint_Cox_P = cut_p,
        Exploratory_Logrank_P = logrank_p,
        Cutpoint_Selection_Method = "Exploratory minimum log-rank P over AST cutpoints with min group N >= 10 and min events >= 2",
        Continuous_HR_Per_SD_AST = meta_row$Continuous_HR_Per_SD_AST,
        Continuous_HR_Per_SD_AST_Lower = meta_row$Continuous_HR_Per_SD_AST_Lower,
        Continuous_HR_Per_SD_AST_Upper = meta_row$Continuous_HR_Per_SD_AST_Upper,
        Continuous_P = meta_row$Continuous_P,
        Continuous_FDR_ObSpecific = meta_row$Continuous_FDR_ObSpecific,
        PA_HR = meta_row$PA_HR,
        PA_P = meta_row$PA_P,
        PA_FDR_ObSpecific = meta_row$PA_FDR_ObSpecific,
        Raw_Prevalence = meta_row$Raw_Prevalence,
        AST_SD = meta_row$AST_SD,
        AST_IQR = meta_row$AST_IQR
    )

    km_data <- km_data %>%
        group_by(Feature, strata) %>%
        arrange(time, .by_group = TRUE) %>%
        group_modify(~bind_rows(
            .x[1, ] %>% mutate(time = 0, n_event = 0, n_censor = 0, survival = 1, lower = 1, upper = 1),
            .x
        )) %>%
        ungroup()

    summary_row <- tibble(
        Feature = feature,
        Feature_Name = meta_row$Feature_Name,
        DA_Direction = meta_row$DA_Direction,
        Cutpoint_AST = cut_ast,
        Cutpoint_Raw_Abundance_Percent = cut_raw,
        Cutpoint_Group_Counts = cutpoint_group_counts,
        Exploratory_Cutpoint_Cox_HR_CI = fmt_hr_ci(cut_hr, cut_lower, cut_upper),
        Exploratory_Cutpoint_Cox_HR = cut_hr,
        Exploratory_Cutpoint_Cox_HR_Lower = cut_lower,
        Exploratory_Cutpoint_Cox_HR_Upper = cut_upper,
        Exploratory_Cutpoint_Cox_P = cut_p,
        Exploratory_Logrank_P = logrank_p,
        Cutpoint_Selection_Method = "Exploratory minimum log-rank P over AST cutpoints with min group N >= 10 and min events >= 2",
        Continuous_HR_Per_SD_AST_CI = meta_row$Continuous_HR_Per_SD_AST_CI,
        Continuous_HR_Per_SD_AST = meta_row$Continuous_HR_Per_SD_AST,
        Continuous_P = meta_row$Continuous_P,
        Continuous_FDR_ObSpecific = meta_row$Continuous_FDR_ObSpecific,
        PA_HR_CI = meta_row$PA_HR_CI,
        PA_P = meta_row$PA_P,
        PA_FDR_ObSpecific = meta_row$PA_FDR_ObSpecific,
        Raw_Prevalence = meta_row$Raw_Prevalence,
        Note = "Exploratory optimal cutpoint is used for visualization; primary inference is continuous AST Cox."
    )

    list(km_data = km_data, summary_row = summary_row)
}

optimal_cutpoint_meta <- ob_specific_survival %>%
    filter(Feature %in% optimal_cutpoint_features) %>%
    arrange(match(Feature, optimal_cutpoint_features))

optimal_cutpoint_results <- lapply(seq_len(nrow(optimal_cutpoint_meta)), function(i) {
    make_optimal_cutpoint_km(optimal_cutpoint_meta$Feature[i], optimal_cutpoint_meta[i, ])
})

optimal_cutpoint_summary <- bind_rows(lapply(optimal_cutpoint_results, `[[`, "summary_row"))
optimal_cutpoint_km_data <- bind_rows(lapply(optimal_cutpoint_results, `[[`, "km_data")) %>%
    mutate(
        Curve_Display = case_when(
            strata == "Low" ~ "Low reference",
            strata == "High" & Continuous_HR_Per_SD_AST < 1 ~ "High abundance (HR<1)",
            strata == "High" & Continuous_HR_Per_SD_AST >= 1 ~ "High abundance (HR>1)",
            TRUE ~ strata
        ),
        Plot_Label = paste0(
            Feature_Name,
            " | cont HR=", sprintf("%.2f", Continuous_HR_Per_SD_AST),
            ", FDR=", fmt_p(Continuous_FDR_ObSpecific),
            " | cut HR=", sprintf("%.2f", Exploratory_Cutpoint_Cox_HR)
        ),
        Plot_Label = factor(Plot_Label, levels = unique(Plot_Label)),
        Curve_ID = paste(Feature, strata, sep = "__")
    )

write_csv(
    optimal_cutpoint_summary,
    file.path(table_dir, "Overall_CRC_DA_Feature_Survival_ObSpecific_OptimalCutpoint_Selected4.csv")
)

write_tsv(
    optimal_cutpoint_km_data,
    file.path(figure_data_dir, "ObSpecific_OptimalCutpoint_KM_Selected4.tsv")
)

optimal_cutpoint_label_map <- optimal_cutpoint_km_data %>%
    distinct(
        Plot_Label,
        Feature_Name,
        Continuous_HR_Per_SD_AST,
        Continuous_P,
        Continuous_FDR_ObSpecific,
        Cutpoint_Raw_Abundance_Percent,
        Exploratory_Cutpoint_Cox_HR,
        Exploratory_Cutpoint_Cox_P
    ) %>%
    mutate(
        Feature_Label = vapply(Feature_Name, function(x) paste(strwrap(x, width = 26), collapse = "\n"), character(1)),
        Facet_Label = paste0(
            Feature_Label,
            "\nContinuous HR=", sprintf("%.2f", Continuous_HR_Per_SD_AST),
            ", P=", fmt_p(Continuous_P),
            ", FDR=", fmt_p(Continuous_FDR_ObSpecific),
            "\nExploratory cut=", signif(Cutpoint_Raw_Abundance_Percent, 3),
            "%, cut-HR=", sprintf("%.2f", Exploratory_Cutpoint_Cox_HR)
        )
    )
optimal_cutpoint_label_map <- setNames(optimal_cutpoint_label_map$Facet_Label, as.character(optimal_cutpoint_label_map$Plot_Label))

p_optimal_cutpoint <- ggplot(
    optimal_cutpoint_km_data,
    aes(x = time, y = survival, group = Curve_ID, color = Curve_Display)
) +
    geom_step(linewidth = 0.75) +
    facet_wrap(~Plot_Label, ncol = 2, labeller = as_labeller(optimal_cutpoint_label_map)) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_color_manual(
        values = c(
            "Low reference" = "grey55",
            "High abundance (HR<1)" = "#3B6EA8",
            "High abundance (HR>1)" = "#C44E52"
        )
    ) +
    labs(
        x = "Time (years)",
        y = "Overall survival probability",
        color = NULL,
        title = "Exploratory optimal-cutpoint KM for selected Ob-specific taxa"
    ) +
    theme_classic(base_size = 12) +
    theme(
        plot.title = element_text(size = 12, face = "plain"),
        strip.background = element_blank(),
        strip.text = element_text(size = 8.2),
        legend.position = "bottom"
    )

ggsave(
    file.path(figure_dir, "ObSpecific_OptimalCutpoint_KM_Selected4.pdf"),
    p_optimal_cutpoint,
    width = 8.2,
    height = 5.8
)
ggsave(
    file.path(figure_dir, "ObSpecific_OptimalCutpoint_KM_Selected4.png"),
    p_optimal_cutpoint,
    width = 8.2,
    height = 5.8,
    dpi = 300
)

single_optimal_plot_data <- optimal_cutpoint_km_data %>%
    mutate(
        Plot_Group = factor(strata, levels = c("Low", "High")),
        Single_Title = paste0(
            Feature_Name,
            " | HR=", sprintf("%.2f", Continuous_HR_Per_SD_AST),
            ", P=", fmt_p(Continuous_P)
        )
    )

for (feature_i in unique(single_optimal_plot_data$Feature)) {
    one_plot_data <- single_optimal_plot_data %>%
        filter(Feature == feature_i)
    one_feature_name <- unique(one_plot_data$Feature_Name)[1]
    one_title <- paste0(
        one_feature_name,
        "\nHR=", sprintf("%.2f", unique(one_plot_data$Continuous_HR_Per_SD_AST)[1]),
        ", P=", fmt_p(unique(one_plot_data$Continuous_P)[1])
    )
    one_prefix <- paste0("ObSpecific_OptimalCutpoint_KM_Single_", safe_name(one_feature_name))

    write_tsv(
        one_plot_data,
        file.path(figure_data_dir, paste0(one_prefix, ".tsv"))
    )

    p_single <- ggplot(
        one_plot_data,
        aes(x = time, y = survival, color = Plot_Group, group = Curve_ID)
    ) +
        geom_step(linewidth = 0.9) +
        coord_cartesian(ylim = c(0, 1)) +
        scale_color_manual(
            values = c("Low" = "grey55", "High" = "#C44E52"),
            breaks = c("Low", "High")
        ) +
        labs(
            x = "Time (years)",
            y = "Overall survival",
            color = NULL,
            title = one_title
        ) +
        theme_classic(base_size = 13) +
        theme(
            plot.title = element_text(size = 12, face = "plain", hjust = 0.5, lineheight = 0.95),
            axis.title = element_text(size = 13),
            axis.text = element_text(size = 11),
            legend.position = "top",
            legend.text = element_text(size = 11),
            legend.key.width = unit(0.75, "cm"),
            plot.margin = margin(6, 7, 5, 6)
        )

    ggsave(
        file.path(figure_dir, paste0(one_prefix, ".pdf")),
        p_single,
        width = 3.25,
        height = 3.15
    )
    ggsave(
        file.path(figure_dir, paste0(one_prefix, ".png")),
        p_single,
        width = 3.25,
        height = 3.15,
        dpi = 300
    )
}

capture.output(sessionInfo(), file = file.path(table_dir, "session_info.txt"))

all_outputs <- c(
    list.files(table_dir, full.names = TRUE),
    list.files(figure_dir, full.names = TRUE),
    list.files(figure_data_dir, full.names = TRUE)
)
writeLines(
    sub(paste0("^", gsub("\\\\", "/", analysis_dir), "/?"), "", gsub("\\\\", "/", all_outputs)),
    con = file.path(analysis_dir, "output_files.txt")
)

message("=== 5. Done ===")
message("Results written to: ", analysis_dir)
print(summary_table)
if (nrow(continuous_df) > 0) {
    message("Continuous AST WithStage models: ", nrow(continuous_df))
    print(continuous_df %>% select(Feature_Name, Ob_DA_Class, HR_Per_SD_AST, HR_Per_IQR_AST, P_Value, P_Adj_FDR) %>% head(10))
}
