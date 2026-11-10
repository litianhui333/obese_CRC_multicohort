# ==============================================================================
# 脚本名称: 29l_Nw_and_TripleOverlap_Trajectory_Barplots.R
# Script Name: 29l_Nw_and_TripleOverlap_Trajectory_Barplots.R
#
# 功能描述 | Description:
#   基于三个 BMI 组（Nw/Ow/Ob）的批次校正丰度数据，绘制：
#   Based on batch-corrected abundance data across three BMI groups, generate:
#
#   PDF 1: Nw-specific Enriched (↑, 来自 Normal 热图中 Nw-specific 且 coef_crc>0)
#   PDF 2: Nw-specific Depleted (↓, 来自 Normal 热图中 Nw-specific 且 coef_crc<0)
#   PDF 3: Triple-Overlap Enriched (↑, Nw+Ow+Ob 三组共有，coef_crc>0)
#   PDF 4: Triple-Overlap Depleted (↓, Nw+Ow+Ob 三组共有，coef_crc<0)
#
# 依赖 | Dependencies:
#   - 已运行 29_BMI_Stratified_Disease_Heatmap.R（生成 Feature_List_*_Refined.csv）
#   - 已运行或将重跑 MMUPHin 批次校正（abd_corrected）
#   - 函数 plot_bmi_trajectory_barplots() 从 29k 脚本中复制/提炼
#
# 输出 | Output:
#   Nw_Specific_Trajectory_3BMI_Enriched.pdf
#   Nw_Specific_Trajectory_3BMI_Depleted.pdf
#   Triple_Overlap_Trajectory_3BMI_Enriched.pdf
#   Triple_Overlap_Trajectory_3BMI_Depleted.pdf
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(MMUPHin)
  library(tibble)
  library(patchwork)
})

# ==============================================================================
# 0. 路径 & 工具函数 | Paths & utility functions
# ==============================================================================
base_dir   <- "G:/My Drive/proj_ob_AA_CRC"
source(file.path(base_dir, "01.scripts/DA_analysis_figure_utils_functions.R"))
source(file.path(base_dir, "01.scripts/DA_function.R"))

input_dir  <- file.path(base_dir, "03.results/04_differential_abundance")
data_dir   <- file.path(base_dir, "02.processed_data/02_aligned_data")
sp_dir     <- file.path(base_dir, "02.processed_data/03_metaphlan_split_data")
heat_dir   <- file.path(base_dir, "03.results/04_differential_abundance/heatmap_final/BMI_Stratified_Disease_Heatmaps")
output_dir <- heat_dir
cache_dir  <- file.path(base_dir, "03.results/04_differential_abundance/heatmap_final/cache")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(cache_dir))  dir.create(cache_dir,  recursive = TRUE)

# ==============================================================================
# 1. 读取菌名单 + 方向匹配过滤 | Read feature lists + direction match filter
# ==============================================================================
message("=== 1. 读取菌名单 + 方向匹配 (Feature_List_*_Refined.csv) | Reading feature lists with direction match ===")

res_file  <- sort(list.files(input_dir, pattern = "microbiome_analysis_all_results_BMI_Group1_0.7.*\\.rds", full.names = TRUE), decreasing = TRUE)[1]
all_analysis_results <- readRDS(res_file)
data_info <- all_analysis_results$overall$disease_analysis$data_info
fa <- data_info$filtered_data
meta <- data_info$metadata_filtered

calculate_sample_direction <- function(features, bmi_group, fa, meta) {
  bmi_mask <- meta$BMI_Group1 == bmi_group
  ctrl_idx <- which(bmi_mask & meta$Disease == "CTRL")
  crc_idx <- which(bmi_mask & meta$Disease == "CRC")
  
  result <- c()
  for (feat in features) {
    if (feat %in% rownames(fa)) {
      ctrl_vals <- as.numeric(fa[feat, ctrl_idx])
      crc_vals <- as.numeric(fa[feat, crc_idx])
      ctrl_pos <- ctrl_vals[!is.na(ctrl_vals) & ctrl_vals > 0]
      crc_pos <- crc_vals[!is.na(crc_vals) & crc_vals > 0]
      
      if (length(ctrl_pos) > 0 && length(crc_pos) > 0) {
        direction <- ifelse(mean(crc_pos) > mean(ctrl_pos), "UP", "DOWN")
        result[feat] <- direction
      }
    }
  }
  result
}

# ---- 1a. Nw-specific ----
nw_csv <- read.csv(file.path(heat_dir, "Feature_List_Normal_Refined.csv"), stringsAsFactors = FALSE)
nw_sample_dir <- calculate_sample_direction(nw_csv$feature, "Normal", fa, meta)
nw_csv$sample_dir <- nw_sample_dir[nw_csv$feature]
nw_csv$direction_match <- ifelse(
  (nw_csv$coef_crc > 0 & nw_csv$sample_dir == "UP") | 
  (nw_csv$coef_crc < 0 & nw_csv$sample_dir == "DOWN"),
  "Match", "Mismatch"
)

nw_en_df <- nw_csv[grepl("Nw-specific", nw_csv$overlap) & nw_csv$direction_match == "Match" & nw_csv$sample_dir == "UP", ]
nw_en_full <- nw_en_df$feature

nw_de_df <- nw_csv[grepl("Nw-specific", nw_csv$overlap) & nw_csv$direction_match == "Match" & nw_csv$sample_dir == "DOWN", ]
nw_de_full <- nw_de_df$feature

message(sprintf("   Nw-specific (matched): Enriched: %d | Depleted: %d", length(nw_en_full), length(nw_de_full)))

# ---- 1b. Ow-specific ----
ow_csv <- read.csv(file.path(heat_dir, "Feature_List_Ow_Refined.csv"), stringsAsFactors = FALSE)
ow_sample_dir <- calculate_sample_direction(ow_csv$feature, "Ow", fa, meta)
ow_csv$sample_dir <- ow_sample_dir[ow_csv$feature]
ow_csv$direction_match <- ifelse(
  (ow_csv$coef_crc > 0 & ow_csv$sample_dir == "UP") | 
  (ow_csv$coef_crc < 0 & ow_csv$sample_dir == "DOWN"),
  "Match", "Mismatch"
)

ow_en_df <- ow_csv[grepl("Ow-specific", ow_csv$overlap) & ow_csv$direction_match == "Match" & ow_csv$sample_dir == "UP", ]
ow_en_full <- ow_en_df$feature

ow_de_df <- ow_csv[grepl("Ow-specific", ow_csv$overlap) & ow_csv$direction_match == "Match" & ow_csv$sample_dir == "DOWN", ]
ow_de_full <- ow_de_df$feature

message(sprintf("   Ow-specific (matched): Enriched: %d | Depleted: %d", length(ow_en_full), length(ow_de_full)))

# ---- 1c. Triple-overlap (from Normal CSV) ----
triple_all <- nw_csv[grepl("Overlap with both Ow and Ob", nw_csv$overlap), ]
triple_en_df <- triple_all[triple_all$direction_match == "Match" & triple_all$sample_dir == "UP", ]
triple_en_full <- triple_en_df$feature

triple_de_df <- triple_all[triple_all$direction_match == "Match" & triple_all$sample_dir == "DOWN", ]
triple_de_full <- triple_de_df$feature

message(sprintf("   Triple-overlap (matched): Enriched: %d | Depleted: %d", length(triple_en_full), length(triple_de_full)))

# ---- 1d. 合并所有需要丰度数据的菌 ----
all_csv_features <- unique(c(nw_en_full, nw_de_full, ow_en_full, ow_de_full, triple_en_full, triple_de_full))

# 状态标签映射 | Status label mapping
status_map <- c(
  setNames(rep("Enriched",        length(nw_en_full)),     nw_en_full),
  setNames(rep("Depleted_NwSpec", length(nw_de_full)),     nw_de_full),
  setNames(rep("Enriched",        length(ow_en_full)),     ow_en_full),
  setNames(rep("Depleted_OwSpec", length(ow_de_full)),     ow_de_full),
  setNames(rep("Enriched",        length(triple_en_full)), triple_en_full),  # shared Enriched
  setNames(rep("Depleted_Triple", length(triple_de_full)), triple_de_full)   # shared Depleted
)
# 类型标签（用于区分 Nw-specific vs Ow-specific vs Triple-overlap）| Category label
category_map <- c(
  setNames(rep("Nw-specific",    length(nw_en_full)),     nw_en_full),
  setNames(rep("Nw-specific",    length(nw_de_full)),     nw_de_full),
  setNames(rep("Ow-specific",    length(ow_en_full)),     ow_en_full),
  setNames(rep("Ow-specific",    length(ow_de_full)),     ow_de_full),
  setNames(rep("Triple-Overlap", length(triple_en_full)), triple_en_full),
  setNames(rep("Triple-Overlap", length(triple_de_full)), triple_de_full)
)

# ==============================================================================
# 2. 加载原始数据 | Load raw data
# ==============================================================================
message("=== 2. 加载数据 | Loading data ===")
meta_file <- sort(list.files(data_dir, pattern = "metadata_aligned.*\\.rds",           full.names = TRUE), decreasing = TRUE)[1]
sp_file   <- sort(list.files(sp_dir,   pattern = "metaphlan_strain_renormalized.*\\.rds", full.names = TRUE), decreasing = TRUE)[1]

metadata             <- readRDS(meta_file)
metaphlan            <- readRDS(sp_file)

rownames(metadata)   <- metadata$Sample.Source.ID.processed
metadata$Disease     <- factor(metadata$Disease,    levels = c("CTRL", "AA", "CRC"))
metadata$BMI_Group1  <- factor(metadata$BMI_Group1, levels = c("Normal", "Ow", "Ob"))
sex_col <- if ("Sex" %in% colnames(metadata)) "Sex" else "Gender"

# ==============================================================================
# 3. MMUPHin 批次校正 (与 Script 29 / 29k 一致) | Batch correction
# ==============================================================================

# 只从 data_info 获取 filtered_data，不依赖缓存
metaphlan_filtered <- data_info$filtered_data
meta_clean <- data_info$metadata_filtered

# 运行批次校正 | Run batch correction
cache_file <- file.path(cache_dir, "abd_corrected_meta_clean.rds")

if (file.exists(cache_file)) {
  message("=== 3. 加载批次校正缓存 | Loading cached batch-corrected data ===")
  cache <- readRDS(cache_file)
  abd_corrected <- cache$abd_corrected
  meta_clean    <- cache$meta_clean
} else {
  message("=== 3. 运行 MMUPHin 批次校正 | Running MMUPHin batch correction ===")

  # 与 Script 29 / 29k 完全一致的函数
  get_robust_features_detailed <- function(data_source, region, bmi_group, comp,
                                           p_cutoff = 0.01, q_cutoff = 0.1, weight_limit = 0.7) {
    list_name <- paste0("disease_analysis_", bmi_group)
    res_obj   <- tryCatch(
      data_source[[region]][["disease_analysis"]][["analysis_results"]][[list_name]][["all_results"]][[comp]],
      error = function(e) NULL)
    if (is.null(res_obj)) return(NULL)
    df_meta    <- res_obj[["meta_fits"]]; df_maaslin <- res_obj[["maaslin_fits"]]
    if (is.null(df_meta) | is.null(df_maaslin)) return(NULL)
    weight_cols <- grep("^weight_", colnames(df_meta), value = TRUE)
    w_mat <- as.matrix(df_meta[, weight_cols]); w_mat[is.na(w_mat)] <- 0
    perform_loo_rescue <- function(f_name, weight_row, original_coef_sign) {
      tryCatch({
        max_idx <- which.max(weight_row); if (length(max_idx) == 0) return(list(saved = FALSE))
        exclude_study    <- gsub("^weight_", "", names(weight_row)[max_idx])
        remaining_studies <- names(df_maaslin)[names(df_maaslin) != exclude_study]
        coefs <- c(); ses <- c()
        for (s in remaining_studies) {
          df_s <- df_maaslin[[s]]
          if ("feature" %in% colnames(df_s) && f_name %in% df_s$feature) {
            idx <- which(df_s$feature == f_name); val_c <- df_s$coef[idx]; val_se <- df_s$stderr[idx]
            if (!is.na(val_c) && !is.na(val_se) && val_se > 1e-10) { coefs <- c(coefs, val_c); ses <- c(ses, val_se) }
          }
        }
        if (length(coefs) < 1) return(list(saved = FALSE))
        weights <- 1/(ses^2); valid_w <- is.finite(weights) & !is.na(weights)
        if (sum(valid_w) == 0) return(list(saved = FALSE))
        weights <- weights[valid_w]; coefs <- coefs[valid_w]; sum_w <- sum(weights)
        if (sum_w == 0) return(list(saved = FALSE))
        combined_coef <- sum(coefs * weights)/sum_w; combined_se <- sqrt(1/sum_w)
        combined_p    <- 2*(1 - pnorm(abs(combined_coef/combined_se)))
        if (is.na(combined_p) || sign(combined_coef) != original_coef_sign) return(list(saved = FALSE))
        if (combined_p < 0.05) return(list(saved = TRUE, coef = combined_coef, pval = combined_p))
        return(list(saved = FALSE))
      }, error = function(e) list(saved = FALSE))
    }
    audit_df <- data.frame(
      feature    = rownames(df_meta),
      Max_Weight = apply(w_mat, 1, function(x) { s <- sum(x, na.rm = TRUE); if (s <= 1e-10) 1 else max(x, na.rm = TRUE)/s }),
      pval       = df_meta$pval, qval = df_meta$qval.fdr, coef = df_meta$coef, stringsAsFactors = FALSE)
    result_list <- list()
    for (i in seq_len(nrow(audit_df))) {
      row <- audit_df[i, ]; f <- row$feature
      is_candidate <- !is.na(row$pval) && (row$pval < p_cutoff | (!is.na(row$qval) && row$qval < q_cutoff))
      if (!is_candidate) next
      status <- "Discard"; final_coef <- row$coef; mw <- row$Max_Weight; if (is.na(mw)) mw <- 1
      if (mw < weight_limit) { status <- "Robust" } else {
        loo_res <- perform_loo_rescue(f, w_mat[f, ], sign(row$coef))
        if (loo_res$saved) { status <- "Rescued"; final_coef <- loo_res$coef }
      }
      if (status != "Discard")
        result_list[[f]] <- data.frame(feature = f, pval = row$pval, qval = row$qval, coef = final_coef, robust_status = status)
    }
    if (length(result_list) == 0) return(NULL)
    return(do.call(rbind, result_list))
  }

  df_robust_Nw <- get_robust_features_detailed(all_analysis_results, "overall", "Normal", "CTRL_vs_CRC")
  df_robust_Ow <- get_robust_features_detailed(all_analysis_results, "overall", "Ow",     "CTRL_vs_CRC")
  df_robust_Ob <- get_robust_features_detailed(all_analysis_results, "overall", "Ob",     "CTRL_vs_CRC")
  pass_nw      <- if (!is.null(df_robust_Nw)) df_robust_Nw$feature else character(0)
  pass_ow      <- if (!is.null(df_robust_Ow)) df_robust_Ow$feature else character(0)
  pass_ob      <- if (!is.null(df_robust_Ob)) df_robust_Ob$feature else character(0)
  all_targets  <- unique(c(pass_nw, pass_ow, pass_ob))

  config <- create_analysis_config(
    transform_method = "AST", analysis_method = "HS", BMI_Group_var = "BMI_Group1",
    covariates = c("Age", "Sex"), study_prevalence_cut = 0.1, abundance_threshold = 0.0001,
    min_studies_ratio = 0.5, pval_type = "both", sig_criterion = "or",
    pval_threshold = 0.02, qval_threshold = 0.2, is_inhouse = FALSE,
    force_keep = all_targets)

  filtering_result   <- smart_data_filtering(t(metaphlan), metadata, config)
  metaphlan_filtered <- filtering_result$filtered_data

  common_samples <- intersect(rownames(metadata), colnames(metaphlan_filtered))
  meta_clean     <- metadata[common_samples, ] %>%
    drop_na(all_of(c("Study.name", "Disease", "BMI_Group1", "Age", sex_col)))
  abd_clean      <- metaphlan_filtered[, rownames(meta_clean)]
  meta_clean$Study.name <- as.factor(meta_clean$Study.name)

  try(unlink("adjust_batch_diagnostic.pdf", force = TRUE))
  fit_adjust    <- adjust_batch(
    feature_abd = abd_clean / 100, batch = "Study.name",
    covariates  = c("Disease", "BMI_Group1", "Age", sex_col),
    data        = meta_clean, control = list(verbose = FALSE))
  abd_corrected <- fit_adjust$feature_abd_adj * 100

# 保存缓存 | Save cache for reuse
saveRDS(list(abd_corrected = abd_corrected, meta_clean = meta_clean, metaphlan_filtered = metaphlan_filtered), cache_file)
  message("   >>> 批次校正完成，已保存缓存 | Batch correction done, cached.")
}

# ==============================================================================
# 4. 名字对齐: CSV 全名 -> 批次矩阵行名 | Name alignment
# ==============================================================================
message("=== 4. 名字对齐 | Name alignment ===")

map_to_matrix <- function(targets, mat_rownames) {
  # 三级匹配策略 | Three-tier matching: exact → SGB ID → species name
  result <- list()
  for (t in targets) {
    if (t %in% mat_rownames) { result[[t]] <- t; next }
    sgb <- regmatches(t, regexpr("SGB[0-9]+", t))
    if (length(sgb) > 0) {
      m <- grep(sgb, mat_rownames, value = TRUE, fixed = TRUE)
      if (length(m) >= 1) { result[[t]] <- m[1]; next }
    }
    sp_key <- gsub("\\|t__.*", "", gsub(".*\\|s__", "s__", t))
    m <- grep(sp_key, mat_rownames, value = TRUE, fixed = TRUE)
    if (length(m) >= 1) { result[[t]] <- m[1] }
  }
  return(result)
}

mapping     <- map_to_matrix(all_csv_features, rownames(abd_corrected))
found       <- names(mapping)
missing     <- setdiff(all_csv_features, found)
matched_rows <- unlist(mapping); names(matched_rows) <- found

if (length(missing) > 0)
  message("   >>> 未匹配（跳过）: ", paste(format_sgb_names(missing), collapse = "; "))
message(sprintf("   >>> 匹配: %d / %d 个目标菌", length(found), length(all_csv_features)))

# ==============================================================================
# 5. 提取三组数据，构建长格式 data frame | Extract per-BMI data
# ==============================================================================
message("=== 5. 提取三组数据 | Extracting per-BMI data ===")
long_list <- list()

for (bmi in c("Normal", "Ow", "Ob")) {
  meta_bmi  <- meta_clean[meta_clean$BMI_Group1 == bmi, ]
  samps_bmi <- intersect(rownames(meta_bmi), colnames(metaphlan_filtered))
  meta_bmi  <- meta_bmi[samps_bmi, ]

  for (csv_name in found) {
    mat_row   <- matched_rows[csv_name]
    raw_abund <- as.numeric(metaphlan_filtered[mat_row, samps_bmi])
    abund_vec <- asin(sqrt(pmin(pmax(raw_abund, 0) / 100, 1)))
    disp_name <- format_sgb_names(csv_name)

    tmp <- data.frame(
      SampleID = samps_bmi,
      Disease  = meta_bmi$Disease,
      BMI_Group = bmi,
      Species   = disp_name,
      Feature   = csv_name,
      Abundance = abund_vec,
      Status    = status_map[csv_name],
      Category  = category_map[csv_name],
      stringsAsFactors = FALSE
    )
    long_list[[paste(bmi, csv_name)]] <- tmp
  }
  message(sprintf("   %s: %d 样本 | %s: %d samples", bmi, length(samps_bmi), bmi, length(samps_bmi)))
}

plot_df <- as.data.frame(do.call(rbind, long_list))
plot_df$Disease   <- factor(plot_df$Disease,   levels = c("CTRL", "AA", "CRC"))
plot_df$BMI_Group <- factor(plot_df$BMI_Group, levels = c("Normal", "Ow", "Ob"))

# ==============================================================================
# 6. 核心绘图函数 (与 29k 完全相同) | Core plotting function (same as 29k)
# ==============================================================================
# ==============================================================================
# 6. 核心绘图函数 — CNS 极简风格 (可复用) | Core Plotting Function — CNS Style
# ==============================================================================

#' plot_bmi_trajectory_barplots (CNS style)
#'
#' 【功能】绘制菌群在三个 BMI 组（Nw/Ow/Ob）内的疾病轨迹柱状图 (CNS 极简分幅风格)。
#' 【Function】Plot microbial disease-trajectory barplots across three BMI groups
#'   using CNS-style per-species subplots assembled via patchwork.
#'
#' @param df            长格式数据框 | Long-format data frame
#' @param status_filter 用于过滤 Status 列的字符串 | Status filter ("Enriched" or others)
#' @param category_filter 用于过滤 Category 列的字符串 | Category filter ("Nw-specific" or "Triple-Overlap")
#' @param filename      输出文件名 | Output filename
#' @param output_dir    输出目录 | Output directory
#' @param sort_by_bmi   排序依据 BMI 组 | BMI group for ranking species
#' @param disease_colors 疾病阶段配色 | Disease stage colors
#' @param page_title    总标题 | Page title (NULL for auto)
#' @param ncol_plots    每行小图数量 | Subplots per row (default 3)
#' @param subplot_w     每列小图宽度 (英寸) | Subplot width
#' @param subplot_h     每行小图高度 (英寸) | Subplot height

plot_bmi_trajectory_barplots <- function(
    df,
    status_filter,
    category_filter,
    filename,
    output_dir,
    sort_by_bmi    = "Ob",
    disease_colors = c("CTRL" = "#2E9FDF", "AA" = "#F39C12", "CRC" = "#E74C3C"),
    page_title     = NULL,
    ncol_plots     = 3,
    subplot_w      = 3.2,
    subplot_h      = 2.2
) {
  # ---- 1. 子集化 | Subset by category and direction ----
  sub_df <- df[df$Category == category_filter, ]
  if (status_filter == "Enriched") {
    sub_df <- sub_df[sub_df$Status == "Enriched", , drop = FALSE]
  } else {
    sub_df <- sub_df[sub_df$Status != "Enriched", , drop = FALSE]
  }
  sub_df <- as.data.frame(sub_df)

  if (nrow(sub_df) == 0) {
    message(sprintf("   >>> [SKIP] No %s taxa in %s.", status_filter, category_filter))
    return(invisible(NULL))
  }

  n_species <- length(unique(sub_df$Species))
  message(sprintf("   >>> Plotting %s [%s]: %d taxa (CNS style)", category_filter, status_filter, n_species))

  # ---- 2. 菌名排序 | Species ordering ----
  sort_fn       <- if (status_filter == "Enriched") dplyr::desc else identity
  species_order <- sub_df[sub_df$BMI_Group == sort_by_bmi & sub_df$Disease == "CRC", ] %>%
    group_by(Species) %>%
    summarise(m = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(sort_fn(m)) %>% pull(Species)
  remaining     <- setdiff(unique(sub_df$Species), species_order)
  sub_df$Species <- factor(sub_df$Species, levels = c(species_order, remaining))

  # ---- 3. 标签 | Labels ----
  bmi_labs <- c("Normal" = "Normal", "Ow" = "Overweight", "Ob" = "Obese")

  if (is.null(page_title)) {
    dir_word   <- if (status_filter == "Enriched") "Enriched (+)" else "Depleted (-)"
    cate_label <- switch(category_filter,
      "Nw-specific"    = "Nw-Specific Signatures",
      "Triple-Overlap" = "Triple-Overlap (Nw + Ow + Ob) Signatures",
      category_filter)
    page_title <- paste0(cate_label, " - ", dir_word, " (n = ", n_species, ")")
  }

  # ---- 4. 逐菌绘图 | Per-taxon subplot loop ----
  plot_list <- list()
  for (sp in levels(sub_df$Species)) {
    sp_df <- sub_df[sub_df$Species == sp, ]

    p <- ggplot(sp_df, aes(x = Disease, y = Abundance, fill = Disease)) +
      # 半透明均值柱 | Semi-transparent bar
      stat_summary(fun = mean, geom = "bar",
                   color = "black", linewidth = 0.35, alpha = 0.65, width = 0.75) +
      # 无帽误差棒 | Capless error bar
      stat_summary(fun.data = mean_se, geom = "errorbar",
                   width = 0, linewidth = 0.5, color = "black") +
      # firebrick 实线轨迹 | firebrick solid trajectory line
      stat_summary(fun = mean, geom = "line", aes(group = 1),
                   color = "firebrick", linewidth = 0.7, linetype = "solid") +
      # 红色边框 + 填充白点 | Dots
      stat_summary(fun = mean, geom = "point",
                   size = 2.2, color = "firebrick", shape = 21, fill = "white", stroke = 1.1) +
      # 固定 Y 轴 | Fixed scale
      facet_wrap(~ BMI_Group, scales = "fixed", labeller = labeller(BMI_Group = bmi_labs), ncol = 3) +
      scale_fill_manual(values = disease_colors) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(title = sp, x = NULL, y = "Rel. Abundance (%)") +
      theme_classic(base_size = 10) +
      theme(
        plot.title       = element_text(face = "bold.italic", size = 10.5, hjust = 0.5, margin = margin(b = 8)),
        axis.line        = element_line(color = "black", linewidth = 0.4),
        axis.ticks       = element_line(color = "black", linewidth = 0.4),
        axis.text.x      = element_text(color = "black", size = 7.5, angle = 45, hjust = 1, vjust = 1),
        axis.text.y      = element_text(color = "black", size = 7.5),
        axis.title.y     = element_text(color = "black", face = "bold", size = 8.5),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 8.5, margin = margin(b = 4)),
        legend.position  = "none",
        panel.spacing.x  = unit(0.2, "lines")
      )
    plot_list[[sp]] <- p
  }

  # ---- 5. Patchwork 拼接 | Assemble with patchwork ----
  final_plot <- wrap_plots(plot_list, ncol = ncol_plots) +
    plot_annotation(
      title    = page_title,
      subtitle = paste0("Mean \u00b1 SEM  |  Batch-Corrected Relative Abundance (%)  |  Y-axis fixed across BMI groups"),
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(size = 9,  color = "grey40", hjust = 0.5)
      )
    )

  # ---- 6. 动态尺寸 & 保存 | Dynamic sizing and save ----
  n_rows   <- ceiling(n_species / ncol_plots)
  fig_w    <- ncol_plots * subplot_w + 0.5
  fig_h    <- max(3, n_rows * subplot_h + 1.5)
  out_path <- file.path(output_dir, filename)

  ggsave(out_path, final_plot, width = fig_w, height = fig_h,
         limitsize = FALSE, dpi = 300, device = grDevices::cairo_pdf)
  message(sprintf("   >>> Saved: %s  (%.1f x %.1f in, %d rows x %d cols)",
                  basename(out_path), fig_w, fig_h, n_rows, ncol_plots))
  return(invisible(NULL))
}

# ==============================================================================
# 7. 生成 PDF | Generate PDFs
# ==============================================================================
message("=== 6. 生成 PDF | Generating PDFs ===")

# Nw-specific
plot_bmi_trajectory_barplots(
  df              = plot_df,
  status_filter   = "Enriched",
  category_filter = "Nw-specific",
  filename        = "Nw_Specific_Trajectory_3BMI_Enriched_Matched.pdf",
  output_dir      = output_dir,
  sort_by_bmi     = "Normal"
)

plot_bmi_trajectory_barplots(
  df              = plot_df,
  status_filter   = "Depleted",
  category_filter = "Nw-specific",
  filename        = "Nw_Specific_Trajectory_3BMI_Depleted_Matched.pdf",
  output_dir      = output_dir,
  sort_by_bmi     = "Normal"
)

# Ow-specific
plot_bmi_trajectory_barplots(
  df              = plot_df,
  status_filter   = "Enriched",
  category_filter = "Ow-specific",
  filename        = "Ow_Specific_Trajectory_3BMI_Enriched_Matched.pdf",
  output_dir      = output_dir,
  sort_by_bmi     = "Ow"
)

plot_bmi_trajectory_barplots(
  df              = plot_df,
  status_filter   = "Depleted",
  category_filter = "Ow-specific",
  filename        = "Ow_Specific_Trajectory_3BMI_Depleted_Matched.pdf",
  output_dir      = output_dir,
  sort_by_bmi     = "Ow"
)

# Triple-Overlap
plot_bmi_trajectory_barplots(
  df              = plot_df,
  status_filter   = "Enriched",
  category_filter = "Triple-Overlap",
  filename        = "Triple_Overlap_Trajectory_3BMI_Enriched_Matched.pdf",
  output_dir      = output_dir,
  sort_by_bmi     = "Ob"
)

plot_bmi_trajectory_barplots(
  df              = plot_df,
  status_filter   = "Depleted",
  category_filter = "Triple-Overlap",
  filename        = "Triple_Overlap_Trajectory_3BMI_Depleted_Matched.pdf",
  output_dir      = output_dir,
  sort_by_bmi     = "Ob"
)

message("✅ 完成! | Done!")
