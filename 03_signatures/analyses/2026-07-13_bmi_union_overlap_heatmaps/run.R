# ==============================================================================
# 脚本名称: 29_BMI_Stratified_Disease_Heatmap.R
# 功能描述:
#   参考 28_Universal_Core_Trajectory_Clustering_core.R 的数据处理流程
#   为 Nw, Ow, Ob 三个 BMI 组分别生成一张热图, 每张热图展示:
#     (1) 个体丰度热图 (该 BMI 组的样本, 按 Disease CTRL/AA/CRC 分割列)
#     (2) 该 BMI 组内按 Disease 的均值热图 (3列)
#     (3) 效应值热图 (CTRL vs AA, CTRL vs CRC), 标注星号
#   行首先按 CTRL_vs_CRC 系数方向拆分 (升高/降低),
#   再在各方向组内按 CTRL/AA/CRC 均值的 Pearson 距离独立聚类
#   最左侧行注释标注: 该菌在另外两组里的 overlap 情况
# ==============================================================================

library(dplyr)
library(tidyr)
library(MMUPHin)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(tibble)

# ==============================================================================
# 0. 加载工具函数
# ==============================================================================
source("G:/My Drive/proj_ob_AA_CRC/01.scripts/DA_analysis_figure_utils_functions.R")
source("G:/My Drive/proj_ob_AA_CRC/01.scripts/DA_function.R")

# ==============================================================================
# 0.5 核心函数: 稳健性审计 (直接移植自 Script 28)
# ==============================================================================
get_robust_features_detailed <- function(data_source, region, bmi_group, comp, 
                                         p_cutoff = 0.01, q_cutoff = 0.1, weight_limit = 0.7) {
  list_name <- paste0("disease_analysis_", bmi_group)
  res_obj <- tryCatch(
    data_source[[region]][["disease_analysis"]][["analysis_results"]][[list_name]][["all_results"]][[comp]],
    error = function(e) NULL
  )
  if (is.null(res_obj)) return(NULL)

  df_meta <- res_obj[["meta_fits"]]
  df_maaslin <- res_obj[["maaslin_fits"]]
  if (is.null(df_meta) | is.null(df_maaslin)) return(NULL)

  weight_cols <- grep("^weight_", colnames(df_meta), value = TRUE)
  w_mat <- as.matrix(df_meta[, weight_cols])
  w_mat[is.na(w_mat)] <- 0

  perform_loo_rescue <- function(f_name, weight_row, original_coef_sign) {
    tryCatch({
      max_idx <- which.max(weight_row)
      if (length(max_idx) == 0) return(list(saved = FALSE))
      exclude_study <- gsub("^weight_", "", names(weight_row)[max_idx])
      remaining_studies <- names(df_maaslin)[names(df_maaslin) != exclude_study]
      
      coefs <- c(); ses <- c()
      for (s in remaining_studies) {
        df_s <- df_maaslin[[s]]
        if ("feature" %in% colnames(df_s) && f_name %in% df_s$feature) {
          idx <- which(df_s$feature == f_name)
          val_c <- df_s$coef[idx]; val_se <- df_s$stderr[idx]
          if (!is.na(val_c) && !is.na(val_se) && val_se > 1e-10) {
            coefs <- c(coefs, val_c); ses <- c(ses, val_se)
          }
        }
      }
      if (length(coefs) < 1) return(list(saved = FALSE))
      weights <- 1 / (ses^2)
      valid_w <- is.finite(weights) & !is.na(weights)
      if (sum(valid_w) == 0) return(list(saved = FALSE))
      weights <- weights[valid_w]; coefs <- coefs[valid_w]
      sum_w <- sum(weights)
      if (sum_w == 0) return(list(saved = FALSE))
      combined_coef <- sum(coefs * weights) / sum_w
      combined_se <- sqrt(1 / sum_w)
      combined_p <- 2 * (1 - pnorm(abs(combined_coef / combined_se)))
      
      if (is.na(combined_p) || sign(combined_coef) != original_coef_sign) return(list(saved = FALSE))
      if (combined_p < 0.05) { # Rescue follows original P < 0.05 for consistency
        return(list(saved = TRUE, coef = combined_coef, pval = combined_p))
      }
      return(list(saved = FALSE))
    }, error = function(e) list(saved = FALSE))
  }

  audit_df <- data.frame(
    feature = rownames(df_meta),
    Max_Weight = apply(w_mat, 1, function(x) {
      s <- sum(x, na.rm = TRUE)
      if (s <= 1e-10) 1 else max(x, na.rm = TRUE) / s
    }),
    pval = df_meta$pval,
    qval = df_meta$qval.fdr,
    coef = df_meta$coef,
    stringsAsFactors = FALSE
  )

  result_list <- list()
  for (i in seq_len(nrow(audit_df))) {
    row <- audit_df[i, ]
    f <- row$feature
    # Stringent Candidate logic: P < p_cutoff OR Q < q_cutoff
    is_candidate <- !is.na(row$pval) && (row$pval < p_cutoff | (!is.na(row$qval) && row$qval < q_cutoff))
    if (!is_candidate) next
    
    status <- "Discard"
    final_coef <- row$coef
    mw <- row$Max_Weight; if (is.na(mw)) mw <- 1
    
    if (mw < weight_limit) { status <- "Robust" } 
    else {
      loo_res <- perform_loo_rescue(f, w_mat[f, ], sign(row$coef))
      if (loo_res$saved) {
        status <- "Rescued"; final_coef <- loo_res$coef
      }
    }
    
    if (status != "Discard") {
      result_list[[f]] <- data.frame(
        feature = f, pval = row$pval, qval = row$qval, 
        coef = final_coef, robust_status = status
      )
    }
  }
  if (length(result_list) == 0) return(NULL)
  return(do.call(rbind, result_list))
}

# ==============================================================================
# 1. 路径设置与数据加载
# ==============================================================================
base_dir <- "G:/My Drive/proj_ob_AA_CRC"
analysis_dir <- file.path(base_dir, "analyses/2026-07-13_bmi_union_overlap_heatmaps")
input_dir <- file.path(base_dir, "03.results/04_differential_abundance")
data_dir <- file.path(base_dir, "02.processed_data/02_aligned_data")
sp_dir <- file.path(base_dir, "02.processed_data/03_metaphlan_split_data")

output_dir <- file.path(analysis_dir, "figures")
table_dir <- file.path(analysis_dir, "tables")
data_fig_dir <- file.path(analysis_dir, "data_for_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_fig_dir, recursive = TRUE, showWarnings = FALSE)

message("=== 1. 加载数据 ===")
meta_file <- sort(list.files(data_dir, pattern = "metadata_aligned.*\\.rds", full.names = TRUE), decreasing = TRUE)[1]
sp_file <- sort(list.files(sp_dir, pattern = "metaphlan_strain_renormalized.*\\.rds", full.names = TRUE), decreasing = TRUE)[1]
res_file <- sort(list.files(input_dir, pattern = "microbiome_analysis_all_results_BMI_Group1_0.7.*\\.rds", full.names = TRUE), decreasing = TRUE)[1]

metadata <- readRDS(meta_file)
metaphlan <- readRDS(sp_file)
all_analysis_results <- readRDS(res_file)

rownames(metadata) <- metadata$Sample.Source.ID.processed
metadata$BMI <- as.numeric(as.character(metadata$BMI))
metadata$Disease <- factor(metadata$Disease, levels = c("CTRL", "AA", "CRC"))
metadata$BMI_Group1 <- factor(metadata$BMI_Group1, levels = c("Normal", "Ow", "Ob"))

# ==============================================================================
# 2. 全局批次校正 (只做一次)
# ==============================================================================
message("=== 2. 全局过滤与 MMUPHin 批次校正 ===")

sex_col <- if ("Sex" %in% colnames(metadata)) "Sex" else "Gender"

# 先确定所有需要保留的特征 (三个 BMI 组的并集, 严格标准: P < 0.01 | FDR < 0.1)
message("   正在识别三组里各自显著的菌群 (严格筛选: P < 0.01 | FDR < 0.1)...")
df_robust_Nw <- get_robust_features_detailed(all_analysis_results, "overall", "Normal", "CTRL_vs_CRC", p_cutoff = 0.01, q_cutoff = 0.1)
df_robust_Ow <- get_robust_features_detailed(all_analysis_results, "overall", "Ow",     "CTRL_vs_CRC", p_cutoff = 0.01, q_cutoff = 0.1)
df_robust_Ob <- get_robust_features_detailed(all_analysis_results, "overall", "Ob",     "CTRL_vs_CRC", p_cutoff = 0.01, q_cutoff = 0.1)

pass_nw <- if (!is.null(df_robust_Nw)) df_robust_Nw$feature else character(0)
pass_ow <- if (!is.null(df_robust_Ow)) df_robust_Ow$feature else character(0)
pass_ob <- if (!is.null(df_robust_Ob)) df_robust_Ob$feature else character(0)
all_targets <- unique(c(pass_nw, pass_ow, pass_ob))

cat(sprintf(
  "   Nw=%d  Ow=%d  Ob=%d  Union=%d\n",
  length(pass_nw), length(pass_ow), length(pass_ob), length(all_targets)
))

# 提前准备全局 P 值与 FDR 矩阵，用于后续的具体 Specificity 判定 (P > 0.05 逻辑)
# 这些值来自 robust_stat_overall 的缓存结果 (v5_Pure)
rds_dir_pure <- file.path(base_dir, "03.results/04_differential_abundance/heatmap_final/Robust_Analysis_v5_Pure")
pure_file <- list.files(rds_dir_pure, pattern = "Robust_Stat_overall_CTRL_vs_CRC.*\\.rds", full.names = TRUE)[1]
dat_pure <- readRDS(pure_file)
m_p_all_pure <- as.matrix(dat_pure$matrix_data$pval)

# 如果 RDS 中没有 qval，则手动计算 FDR
if ("qval" %in% names(dat_pure$matrix_data)) {
  m_q_all_pure <- as.matrix(dat_pure$matrix_data$qval)
} else {
  message("   RDS 中未发现 qval 矩阵，正在根据 pval 矩阵计算 FDR...")
  m_q_all_pure <- apply(m_p_all_pure, 2, p.adjust, method = "fdr")
  rownames(m_q_all_pure) <- rownames(m_p_all_pure)
}

# 适配列名 (RDS 里是 "Normal", "Overweight", "Obese")
if (ncol(m_p_all_pure) >= 2) {
  colnames(m_p_all_pure)[colnames(m_p_all_pure) == "Overweight"] <- "Ow"
  colnames(m_p_all_pure)[colnames(m_p_all_pure) == "Obese"]      <- "Ob"
  colnames(m_q_all_pure)[colnames(m_q_all_pure) == "Overweight"] <- "Ow"
  colnames(m_q_all_pure)[colnames(m_q_all_pure) == "Obese"]      <- "Ob"
}

# Build the display set from the union of differential taxa.
# A taxon is shown in every BMI heatmap where it is either robust in that BMI
# group or has the overlap-level CTRL-vs-CRC signal (P < 0.05) in that BMI group.
bmi_short <- c("Normal" = "Nw", "Ow" = "Ow", "Ob" = "Ob")
get_display_bmis <- function(f) {
  robust_bmis <- c(
    if (f %in% pass_nw) "Normal",
    if (f %in% pass_ow) "Ow",
    if (f %in% pass_ob) "Ob"
  )
  p_bmis <- c()
  for (bg in c("Normal", "Ow", "Ob")) {
    pv <- m_p_all_pure[f, bg]
    if (!is.na(pv) && pv < 0.05) p_bmis <- c(p_bmis, bg)
  }
  intersect(c("Normal", "Ow", "Ob"), unique(c(robust_bmis, p_bmis)))
}
display_bmis_by_feature <- lapply(all_targets, get_display_bmis)
names(display_bmis_by_feature) <- all_targets

display_plan <- data.frame(
  feature = all_targets,
  display_bmis = vapply(display_bmis_by_feature, function(x) paste(bmi_short[x], collapse = ";"), character(1)),
  display_class = vapply(display_bmis_by_feature, function(x) {
    if (length(x) == 3) return("BMI-independent")
    if (length(x) == 2) return(paste0("Overlap_", paste(bmi_short[x], collapse = "_")))
    paste0(bmi_short[x], "-specific")
  }, character(1)),
  robust_in_Nw = all_targets %in% pass_nw,
  robust_in_Ow = all_targets %in% pass_ow,
  robust_in_Ob = all_targets %in% pass_ob,
  p_Nw = m_p_all_pure[all_targets, "Normal"],
  p_Ow = m_p_all_pure[all_targets, "Ow"],
  p_Ob = m_p_all_pure[all_targets, "Ob"],
  stringsAsFactors = FALSE
)
write.table(display_plan, file.path(table_dir, "display_feature_plan_union_overlap.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)

config <- create_analysis_config(
  transform_method     = "AST",
  analysis_method      = "HS",
  BMI_Group_var        = "BMI_Group1",
  covariates           = c("Age", "Sex"),
  study_prevalence_cut = 0.1,
  abundance_threshold  = 0.0001,
  min_studies_ratio    = 0.5,
  pval_type            = "both",
  sig_criterion        = "or",
  pval_threshold       = 0.02,
  qval_threshold       = 0.2,
  is_inhouse           = FALSE,
  force_keep           = all_targets
)

filtering_result <- smart_data_filtering(t(metaphlan), metadata, config)
metaphlan_filtered <- filtering_result$filtered_data

common_samples <- intersect(rownames(metadata), colnames(metaphlan_filtered))
meta_clean <- metadata[common_samples, ] %>%
  drop_na(all_of(c("Study.name", "Disease", "BMI_Group1", "Age", sex_col)))
abd_clean <- metaphlan_filtered[, rownames(meta_clean)]
meta_clean$Study.name <- as.factor(meta_clean$Study.name)

message("   正在执行 MMUPHin 批次校正...")
try(unlink("adjust_batch_diagnostic.pdf", force = TRUE))
fit_adjust <- adjust_batch(
  feature_abd = abd_clean / 100,
  batch       = "Study.name",
  covariates  = c("Disease", "BMI_Group1", "Age", sex_col),
  data        = meta_clean,
  control     = list(verbose = FALSE)
)
abd_corrected <- fit_adjust$feature_abd_adj
abd_ast <- asin(sqrt(pmax(fit_adjust$feature_abd_adj, 0)))

# ==============================================================================
# 3. 绘图循环 (一个 BMI 组绘制一张热图)
# ==============================================================================
bmi_label_map <- c("Normal" = "Nw (Normal Weight)", "Ow" = "Ow (Overweight)", "Ob" = "Ob (Obese)")
bmi_color_map <- c("Normal" = "#5AAE61", "Ow" = "#FF8C00", "Ob" = "#9970AB")

for (BMI_GRP in c("Normal", "Ow", "Ob")) {
  message(sprintf("\n\n====== 正在处理 %s 组 ======", BMI_GRP))

  # 该组的显著菌名单
  pass_this <- switch(BMI_GRP,
    Normal = pass_nw,
    Ow = pass_ow,
    Ob = pass_ob
  )
  if (length(pass_this) == 0) {
    message("   [跳过] 该 BMI 组无显著特征。")
    next
  }
  display_this <- all_targets[vapply(all_targets, function(f) {
    BMI_GRP %in% display_bmis_by_feature[[f]]
  }, logical(1))]
  available_targets <- intersect(display_this, rownames(abd_ast))
  if (length(available_targets) == 0) {
    message("   [跳过] 所有特征批次校正后丢失。")
    next
  }
  message(sprintf("   展示特征: own robust=%d, union/overlap display=%d, available=%d",
    length(pass_this), length(display_this), length(available_targets)))

  # 仅限该 BMI 组的样本
  meta_bmi <- meta_clean %>% filter(BMI_Group1 == BMI_GRP)
  samples_bmi <- meta_bmi %>%
    arrange(Disease) %>%
    pull(Sample.Source.ID.processed)
  samples_bmi <- intersect(samples_bmi, colnames(abd_ast))
  meta_bmi_ord <- meta_bmi[samples_bmi, ]

  # ============================================================
  # A. 均值矩阵 (仅限该 BMI 组, 3列: CTRL, AA, CRC)
  # ============================================================
  df_mean_bmi <- t(abd_ast[available_targets, samples_bmi, drop = FALSE]) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("SampleID") %>%
    inner_join(meta_bmi_ord %>% select(Sample.Source.ID.processed, Disease),
      by = c("SampleID" = "Sample.Source.ID.processed")
    ) %>%
    pivot_longer(cols = all_of(available_targets), names_to = "Feature", values_to = "val") %>%
    group_by(Feature, Disease) %>%
    summarise(mean_val = mean(val, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Disease, values_from = mean_val)

  mat_mean_bmi <- as.data.frame(df_mean_bmi)
  rownames(mat_mean_bmi) <- mat_mean_bmi$Feature
  mat_mean_bmi$Feature <- NULL
  # 确保列顺序
  available_diseases <- intersect(c("CTRL", "AA", "CRC"), colnames(mat_mean_bmi))
  mat_mean_bmi <- mat_mean_bmi[, available_diseases, drop = FALSE]

  # Z-score
  mat_mean_bmi_z <- t(apply(mat_mean_bmi, 1, scale))
  colnames(mat_mean_bmi_z) <- available_diseases
  rownames(mat_mean_bmi_z) <- rownames(mat_mean_bmi)

  # ============================================================
  # B. 行方向判断 + 方向分组内 Pearson 聚类
  # ============================================================
  # ============================================================
  # D. 效应值矩阵 (先计算，用于聚类和展示)
  # ============================================================
  get_stats <- function(comp, feats) {
    extract_model_stats_final(all_analysis_results, "overall", "disease", BMI_GRP, comp, feats)
  }
  df_ca_all <- get_stats("CTRL_vs_AA", available_targets)
  df_cc_all <- get_stats("CTRL_vs_CRC", available_targets)

  # 得到已经审计过的菌的 coef 覆盖
  patch_vec_all <- function(df_raw, df_audit, target_features) {
    coef_vec <- df_raw[match(target_features, df_raw$feature), "coef"]
    if (!is.null(df_audit)) {
      common <- intersect(target_features, df_audit$feature)
      coef_vec[match(common, target_features)] <- df_audit[match(common, df_audit$feature), "coef"]
    }
    coef_vec
  }
  df_audit_this <- switch(BMI_GRP,
    Normal = df_robust_Nw,
    Ow = df_robust_Ow,
    Ob = df_robust_Ob
  )

  m_c_eff_all <- cbind(
    "CTRL vs AA"  = df_ca_all[match(available_targets, df_ca_all$feature), "coef"],
    "CTRL vs CRC" = patch_vec_all(df_cc_all, df_audit_this, available_targets)
  )
  rownames(m_c_eff_all) <- available_targets

  # ============================================================
  # B. 行排序与拆分逻辑 (根据用户最新要求: 4路拆分 + Overlap标签归类 + AA排序)
  #   1. 按 CRC 和 AA 效应值正负拆分 (4 个趋势组)
  #   2. 在每个小组内部，按 Overlap 标签归类 (不让相同标签的菌分散)
  #   3. 在最终小组内，按 AA 效应值从小到大排序
  # ============================================================
  # 提前计算 Overlap 标签用于排序
  other_bmis <- setdiff(c("Normal", "Ow", "Ob"), BMI_GRP)
  pass_other1 <- switch(other_bmis[1], Normal = pass_nw, Ow = pass_ow, Ob = pass_ob)
  pass_other2 <- switch(other_bmis[2], Normal = pass_nw, Ow = pass_ow, Ob = pass_ob)
  bmi_short <- c("Normal" = "Nw", "Ow" = "Ow", "Ob" = "Ob")

  other1_label <- bmi_short[other_bmis[1]]
  other2_label <- bmi_short[other_bmis[2]]

  other1_full <- paste0("Overlap with ", other1_label)
  other2_full <- paste0("Overlap with ", other2_label)
  both_full   <- paste0("Overlap with both ", other1_label, " and ", other2_label)
  specific_full <- paste0(bmi_short[BMI_GRP], "-specific")

  get_overlap_label <- function(f) {
    # 新逻辑：
    # 1. 目标组需要通过robust标准 (P<0.01 | FDR<0.1)，已在pass_*里
    # 2. 看其他组的P值：
    #    - 其他两组都 > 0.05 → specific
    #    - 一组 < 0.05，另一组 > 0.05 → overlap with 那一组
    #    - 两组都 < 0.05 → BMI-independent (三组都显著)
    
    in_nw <- f %in% pass_nw
    in_ow <- f %in% pass_ow
    in_ob <- f %in% pass_ob
    
    # 原始P值
    p_nw <- m_p_all_pure[f, "Normal"]; if(is.na(p_nw)) p_nw <- 1
    p_ow <- m_p_all_pure[f, "Ow"];     if(is.na(p_ow)) p_ow <- 1
    p_ob <- m_p_all_pure[f, "Ob"];     if(is.na(p_ob)) p_ob <- 1
    
    # 判断当前处理组
    if (BMI_GRP == "Normal") {
      # Nw组显著，检查Ow和Ob
      if (p_ow > 0.05 & p_ob > 0.05) {
        return(specific_full)  # Nw-specific
      } else if (p_ow < 0.05 & p_ob > 0.05) {
        return(other1_full)   # Overlap with Ow
      } else if (p_ow > 0.05 & p_ob < 0.05) {
        return(other2_full)   # Overlap with Ob
      } else {
        return(both_full)     # BMI-independent (三组都显著)
      }
    } else if (BMI_GRP == "Ow") {
      # Ow组显著，检查Nw和Ob
      if (p_nw > 0.05 & p_ob > 0.05) {
        return(specific_full)  # Ow-specific
      } else if (p_nw < 0.05 & p_ob > 0.05) {
        return(other1_full)   # Overlap with Nw
      } else if (p_nw > 0.05 & p_ob < 0.05) {
        return(other2_full)   # Overlap with Ob
      } else {
        return(both_full)     # BMI-independent
      }
    } else if (BMI_GRP == "Ob") {
      # Ob组显著，检查Nw和Ow
      if (p_nw > 0.05 & p_ow > 0.05) {
        return(specific_full)  # Ob-specific
      } else if (p_nw < 0.05 & p_ow > 0.05) {
        return(other1_full)   # Overlap with Nw
      } else if (p_nw > 0.05 & p_ow < 0.05) {
        return(other2_full)   # Overlap with Ow
      } else {
        return(both_full)     # BMI-independent
      }
    }
    
    return(specific_full)  # 兜底
  }

  df_sorting <- data.frame(
    feature = available_targets,
    coef_crc = m_c_eff_all[, "CTRL vs CRC"],
    coef_aa = m_c_eff_all[, "CTRL vs AA"],
    pval_aa = df_ca_all$pval[match(available_targets, df_ca_all$feature)],
    stringsAsFactors = FALSE
  )

  df_sorting$overlap_label <- sapply(df_sorting$feature, get_overlap_label)
  df_sorting$overlap_fixed <- factor(df_sorting$overlap_label, levels = c(specific_full, other1_full, other2_full, both_full))

  # 分组颜色 (Enriched vs Depleted)
  category_pal <- c(
    "Enriched in CRC" = "#B2182B", # Red
    "Depleted in CRC" = "#2166AC"  # Blue
  )

  # 定义趋势小组 (仅分 Enriched vs Depleted)
  df_sorting <- df_sorting %>%
    mutate(
      joint_group = factor(
        ifelse(coef_crc > 0, "Enriched in CRC", "Depleted in CRC"),
        levels = c("Enriched in CRC", "Depleted in CRC")
      )
    ) %>%
    arrange(joint_group, overlap_fixed, coef_aa)

  feature_order <- df_sorting$feature
  direction_split <- df_sorting$joint_group
  names(direction_split) <- feature_order

  mat_mean_bmi_z <- mat_mean_bmi_z[feature_order, ]

  # ============================================================
  # C. 个体矩阵 (Study-specific Z-score)
  # ============================================================
  df_indi_long <- t(abd_ast[feature_order, samples_bmi, drop = FALSE]) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("SampleID") %>%
    inner_join(meta_bmi_ord %>% select(Sample.Source.ID.processed, Study.name),
      by = c("SampleID" = "Sample.Source.ID.processed")
    ) %>%
    pivot_longer(cols = all_of(feature_order), names_to = "Feature", values_to = "Abundance") %>%
    group_by(Feature, Study.name) %>%
    mutate(Abundance_Scaled = as.vector(scale(Abundance))) %>%
    ungroup() %>%
    mutate(Abundance_Scaled = replace_na(Abundance_Scaled, 0))

  mat_indi_z <- df_indi_long %>%
    select(SampleID, Feature, Abundance_Scaled) %>%
    pivot_wider(names_from = SampleID, values_from = Abundance_Scaled) %>%
    column_to_rownames("Feature") %>%
    as.matrix()
  mat_indi_z <- mat_indi_z[feature_order, samples_bmi, drop = FALSE]
  #mat_indi_z[mat_indi_z > 3.0] <- 3.0
  #mat_indi_z[mat_indi_z < -3.0] <- -3.0

  # ============================================================
  # D. 效应值矩阵 (CTRL vs AA 和 CTRL vs CRC, 该 BMI 组)
  # ============================================================
  m_c_eff <- m_c_eff_all[feature_order, ]

  # 注意：df_ca 和 df_cc 这里只用于 pval 和 qval，coef 我们用上面算好的 m_c_eff_all
  df_ca <- get_stats("CTRL_vs_AA", feature_order)
  df_cc <- get_stats("CTRL_vs_CRC", feature_order)

  m_p_eff <- cbind(
    "CTRL vs AA"  = df_ca[match(feature_order, df_ca$feature), "pval"],
    "CTRL vs CRC" = df_cc[match(feature_order, df_cc$feature), "pval"]
  )
  m_q_eff <- cbind(
    "CTRL vs AA"  = df_ca[match(feature_order, df_ca$feature), "qval.fdr"],
    "CTRL vs CRC" = df_cc[match(feature_order, df_cc$feature), "qval.fdr"]
  )
  rownames(m_c_eff) <- rownames(m_p_eff) <- rownames(m_q_eff) <- feature_order

  # mask: CTRL vs CRC 仅显示通过 LOO 的菌
  m_p_eff[!feature_order %in% pass_this, "CTRL vs CRC"] <- 1.0
  m_q_eff[!feature_order %in% pass_this, "CTRL vs CRC"] <- 1.0

  # ============================================================
  # E. 左侧行注释 (直接复用排序时算好的 Overlap 标签)
  # ============================================================
  overlap_vec <- df_sorting$overlap_label
  names(overlap_vec) <- feature_order

  # 1. 定义 Overlap 颜色 (动态标签映射)
  overlap_levels <- c(specific_full, other1_full, other2_full, both_full)
  overlap_pal <- setNames(
    c("#832a3dff", "#1c5579ff", "#8b8909ff", "black"),
    overlap_levels
  )

  # 2. 定义趋势颜色 (方块)
  category_pal <- c(
    "Enriched in CRC" = "#B2182B",
    "Depleted in CRC" = "#2166AC"
  )

  # 3. 菌名颜色按照 Overlap 状态染色
  row_name_colors <- overlap_pal[as.character(overlap_vec)]

  left_ha <- rowAnnotation(
    Trend   = direction_split[feature_order],
    Overlap = overlap_vec[feature_order],
    col = list(
      Trend   = category_pal,
      Overlap = overlap_pal
    ),
    show_annotation_name = TRUE,
    annotation_name_side = "top",
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    simple_anno_size = unit(4, "mm"),
    show_legend = c(TRUE, TRUE)
  )

  # ============================================================
  # F. 热图定义
  # ============================================================
  clean_names <- format_sgb_names(feature_order)

col_fun_indi <- colorRamp2(
    c(-1.5,-1,-0.8, -0.5,0, 0.5, 1.5, 3),
    c("#08306B", "#2166AC",  "#92C5DE", "#BFDFF0","white", "#F4A582", "#D6604D", "#e4051f")
  )
  col_fun_mean <- colorRamp2(c(-1.5, 0, 1.5), c("#2166AC", "white", "#B2182B"))
  col_fun_eff <- colorRamp2(
    c(-0.05, -0.005, 0, 0.005, 0.05),
    c("#313695", "#Aeb6e8", "white", "#ffadad", "#a50026")
  )

  top_ha_indi <- HeatmapAnnotation(
    Disease = meta_bmi_ord$Disease,
    col = list(Disease = c("CTRL" = "#2E9FDF", "AA" = "#F39C12", "CRC" = "#E74C3C")),
    show_annotation_name = FALSE
  )

  # 个体热图 — 按 direction 拆分，组内按效应值趋势聚类
  clean_names_reorder <- format_sgb_names(feature_order)

  ht_indi <- Heatmap(
    mat_indi_z[feature_order, samples_bmi, drop = FALSE],
    name = "Z-score\n(Individual)",
    column_split = factor(meta_bmi_ord$Disease, levels = c("CTRL", "AA", "CRC")),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    row_split = direction_split,
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 12, fontface = "bold"),
    show_row_names = FALSE,
    show_column_names = FALSE,
    top_annotation = top_ha_indi,
    left_annotation = left_ha,
    col = col_fun_indi,
    row_gap = unit(1.0, "mm"),
    column_gap = unit(0.5, "mm"),
    border = TRUE,
    width = unit(11, "cm"),
    use_raster = FALSE
  )


  # 效应值热图
  m_c_eff_reorder <- m_c_eff[feature_order, , drop = FALSE]
  m_p_eff_reorder <- m_p_eff[feature_order, , drop = FALSE]
  m_q_eff_reorder <- m_q_eff[feature_order, , drop = FALSE]

  # 效应值热图 — 使用 layer_fun (兼容 row_split)
  ht_eff <- Heatmap(
    m_c_eff_reorder,
    name = "Effect Size\n(LogFC)",
    col = col_fun_eff,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_names = TRUE,
    row_labels = clean_names_reorder,
    row_names_gp = gpar(fontface = "italic", fontsize = 9, col = row_name_colors),
    row_names_side = "right",
    row_names_max_width = unit(10, "cm"), # 关键：增加行名最大宽度
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 11),
    column_names_rot = 45,
    layer_fun = function(j, i, x, y, w, h, fill) {
      pv <- pindex(m_p_eff_reorder, i, j)
      qv <- pindex(m_q_eff_reorder, i, j)
      is_strong <- (!is.na(qv) & qv < 0.1) | (!is.na(pv) & pv < 0.01)
      is_weak <- (!is.na(pv) & pv < 0.05) & !is_strong
      if (any(is_strong)) grid.text("**", x[is_strong], y[is_strong], gp = gpar(fontsize = 10))
      if (any(is_weak)) grid.text("*", x[is_weak], y[is_weak], gp = gpar(fontsize = 10, col = "grey30"))
    },
    border = TRUE,
    width = unit(2, "cm")
  )

  ht_list <- ht_indi + ht_eff # 移除均值热图

  # ============================================================
  # G. 输出 PDF
  # ============================================================
  n_feat <- length(feature_order)
  pdf_height <- max(3, n_feat * 0.12 + 2.5) # 稍微增加高度以适应更简单的分组
  pdf_width <- 25 # 宽度大幅减小
  pdf_file <- file.path(output_dir, paste0("Disease_Heatmap_", BMI_GRP, "_Refined_v5_Strict.pdf"))

  pdf(pdf_file, width = pdf_width, height = pdf_height)
  draw(ht_list,
    main_heatmap = "Z-score\n(Individual)",
    column_title = paste0(bmi_label_map[BMI_GRP], " — CTRL vs CRC Significant Taxa (n=", n_feat, ")"),
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    padding = unit(c(5, 5, 5, 40), "mm") # 增加右侧 padding
  )
  dev.off()

  message(sprintf("✅ %s 组热图已输出: %s", BMI_GRP, pdf_file))

  # 同时保存特征信息 CSV
  csv_out <- data.frame(
    feature    = feature_order,
    clean_name = clean_names_reorder,
    overlap    = overlap_vec[feature_order],
    coef_crc   = m_c_eff_reorder[, "CTRL vs CRC"],
    pval_crc   = m_p_eff_reorder[, "CTRL vs CRC"],
    qval_crc   = m_q_eff_reorder[, "CTRL vs CRC"],
    coef_aa    = m_c_eff_reorder[, "CTRL vs AA"],
    pval_aa    = m_p_eff_reorder[, "CTRL vs AA"]
  )
  write.csv(csv_out, file.path(table_dir, paste0("Feature_List_", BMI_GRP, "_Refined.csv")), row.names = FALSE)

  individual_plot_data <- df_indi_long %>%
    inner_join(meta_bmi_ord %>% select(Sample.Source.ID.processed, Disease, BMI_Group1),
      by = c("SampleID" = "Sample.Source.ID.processed")
    ) %>%
    mutate(
      Feature = factor(Feature, levels = feature_order),
      SampleID = factor(SampleID, levels = samples_bmi),
      trend = as.character(direction_split[as.character(Feature)]),
      overlap = as.character(overlap_vec[as.character(Feature)])
    ) %>%
    arrange(Feature, SampleID)
  write.table(individual_plot_data,
    file.path(data_fig_dir, paste0("Disease_Heatmap_", BMI_GRP, "_individual_z.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE)

  effect_plot_data <- data.frame(
    feature = rep(feature_order, times = ncol(m_c_eff_reorder)),
    clean_name = rep(clean_names_reorder, times = ncol(m_c_eff_reorder)),
    comparison = rep(colnames(m_c_eff_reorder), each = length(feature_order)),
    coef = as.vector(m_c_eff_reorder),
    pval = as.vector(m_p_eff_reorder),
    qval = as.vector(m_q_eff_reorder),
    trend = rep(as.character(direction_split[feature_order]), times = ncol(m_c_eff_reorder)),
    overlap = rep(as.character(overlap_vec[feature_order]), times = ncol(m_c_eff_reorder)),
    stringsAsFactors = FALSE
  )
  write.table(effect_plot_data,
    file.path(data_fig_dir, paste0("Disease_Heatmap_", BMI_GRP, "_effect.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE)
}

# ==============================================================================
# 4. 新增: Triple Overlap (核心共有菌) 总结热图 (Individual Samples 版)
# ==============================================================================
message("\n====== 正在生成 Triple Overlap (核心共有菌) 总结热图 (Individual) ======")

triple_targets <- intersect(pass_nw, intersect(pass_ow, pass_ob))

if (FALSE && length(triple_targets) > 0) {
  # A. 提取所有三组样本并排序
  meta_triple <- meta_clean %>%
    filter(BMI_Group1 %in% c("Normal", "Ow", "Ob")) %>%
    mutate(
      BMI_Group1 = factor(BMI_Group1, levels = c("Normal", "Ow", "Ob")),
      Disease = factor(Disease, levels = c("CTRL", "AA", "CRC")),
      # 创建用于分列的复合因子
      Split_Col = factor(paste0(BMI_Group1, "_", Disease), levels = c(
        "Normal_CTRL", "Normal_AA", "Normal_CRC",
        "Ow_CTRL", "Ow_AA", "Ow_CRC",
        "Ob_CTRL", "Ob_AA", "Ob_CRC"
      ))
    ) %>%
    arrange(Split_Col)
  
  samples_triple <- meta_triple$Sample.Source.ID.processed
  
  # B. 提取并标准化丰度 (按 Cohort 内部进行 Z-score)
  mat_triple_raw <- abd_ast[triple_targets, samples_triple, drop = FALSE]
  
  # 分 Cohort 进行标准化
  mat_triple_z_list <- list()
  for (bg in c("Normal", "Ow", "Ob")) {
    s_bg <- meta_triple$Sample.Source.ID.processed[meta_triple$BMI_Group1 == bg]
    if (length(s_bg) > 0) {
      m_bg <- mat_triple_raw[, s_bg, drop = FALSE]
      z_bg <- t(apply(m_bg, 1, function(x) { s <- sd(x, na.rm=TRUE); if(s==0) x else (x - mean(x, na.rm=TRUE))/s }))
      mat_triple_z_list[[bg]] <- z_bg
    }
  }
  mat_triple_z <- do.call(cbind, mat_triple_z_list)
  mat_triple_z <- mat_triple_z[triple_targets, samples_triple] # 恢复原本列顺序
  
  # 截断离群值
  mat_triple_z[mat_triple_z > 3] <- 3
  mat_triple_z[mat_triple_z < -3] <- -3

  # C. 提取 3 组效应值 (逻辑不变)
  get_eff_triple <- function(bmi_grp) {
    ca <- extract_model_stats_final(all_analysis_results, "overall", "disease", bmi_grp, "CTRL_vs_AA", triple_targets)
    cc <- extract_model_stats_final(all_analysis_results, "overall", "disease", bmi_grp, "CTRL_vs_CRC", triple_targets)
    df_robust <- switch(bmi_grp, Normal = df_robust_Nw, Ow = df_robust_Ow, Ob = df_robust_Ob)
    cc$coef <- df_robust$coef[match(cc$feature, df_robust$feature)]
    list(ca = ca, cc = cc)
  }
  eff_nw <- get_eff_triple("Normal")
  eff_ow <- get_eff_triple("Ow")
  eff_ob <- get_eff_triple("Ob")

  m_coef <- cbind(
    Nw_AA = eff_nw$ca$coef, Nw_CRC = eff_nw$cc$coef,
    Ow_AA = eff_ow$ca$coef, Ow_CRC = eff_ow$cc$coef,
    Ob_AA = eff_ob$ca$coef, Ob_CRC = eff_ob$cc$coef
  )
  rownames(m_coef) <- triple_targets

  m_pval <- cbind(
    Nw_AA = eff_nw$ca$pval, Nw_CRC = eff_nw$cc$pval,
    Ow_AA = eff_ow$ca$pval, Ow_CRC = eff_ow$cc$pval,
    Ob_AA = eff_ob$ca$pval, Ob_CRC = eff_ob$cc$pval
  )
  m_qval <- cbind(
    Nw_AA = eff_nw$ca$qval.fdr, Nw_CRC = eff_nw$cc$qval.fdr,
    Ow_AA = eff_ow$ca$qval.fdr, Ow_CRC = eff_ow$cc$qval.fdr,
    Ob_AA = eff_ob$ca$qval.fdr, Ob_CRC = eff_ob$cc$qval.fdr
  )
  rownames(m_pval) <- rownames(m_qval) <- triple_targets

  # D. 排序逻辑 (按平均 CRC 趋势)
  avg_crc <- rowMeans(m_coef[, c("Nw_CRC", "Ow_CRC", "Ob_CRC")])
  direction <- ifelse(avg_crc > 0, "Consistent Increase", "Consistent Decrease")
  avg_aa <- rowMeans(m_coef[, c("Nw_AA", "Ow_AA", "Ob_AA")])
  feature_order <- triple_targets[order(direction, avg_aa)]

  # E. 绘图
  clean_names <- format_sgb_names(feature_order)
  col_abd <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
  col_eff_fun <- colorRamp2(c(-0.05, 0, 0.05), c("#313695", "white", "#a50026"))

  ht_main <- Heatmap(mat_triple_z[feature_order, ], 
    name = "Individual\nZ-score", 
    col = col_abd,
    cluster_columns = FALSE, cluster_rows = FALSE,
    column_split = meta_triple$Split_Col,
    row_split = direction[match(feature_order, triple_targets)],
    row_title_rot = 0, row_title_gp = gpar(fontsize = 8, fontface = "bold"),
    show_row_names = FALSE, show_column_names = FALSE,
    column_gap = unit(0.5, "mm"),
    border = TRUE, width = unit(12, "cm"),
    use_raster = FALSE)

  draw_eff <- function(m_c, m_p, m_q, title, show_names = FALSE) {
    Heatmap(m_c[feature_order, ], col = col_eff_fun, cluster_columns = FALSE, cluster_rows = FALSE,
      column_labels = c("AA", "CRC"), column_title = title, column_title_gp = gpar(fontsize = 8),
      show_row_names = show_names, row_labels = clean_names, 
      row_names_gp = gpar(fontsize = 7, fontface = "italic"),
      layer_fun = function(j, i, x, y, w, h, fill) {
        pv <- pindex(m_p[feature_order, ], i, j); qv <- pindex(m_q[feature_order, ], i, j)
        is_strong <- (!is.na(qv) & qv < 0.1) | (!is.na(pv) & pv < 0.01)
        is_weak <- (!is.na(pv) & pv < 0.05) & !is_strong
        if (any(is_strong)) grid.text("**", x[is_strong], y[is_strong], gp = gpar(fontsize = 10))
        if (any(is_weak)) grid.text("*", x[is_weak], y[is_weak], gp = gpar(fontsize = 10, col = "grey30"))
      }, border = TRUE, width = unit(1.8, "cm"))
  }

  ht_nw <- draw_eff(m_coef[, 1:2], m_pval[, 1:2], m_qval[, 1:2], "Nw")
  ht_ow <- draw_eff(m_coef[, 3:4], m_pval[, 3:4], m_qval[, 3:4], "Ow")
  ht_ob <- draw_eff(m_coef[, 5:6], m_pval[, 5:6], m_qval[, 5:6], "Ob", show_names = TRUE)

  final_triple <- ht_main + ht_nw + ht_ow + ht_ob
  pdf_file <- file.path(output_dir, "Triple_Overlap_Summary_Heatmap_Individual.pdf")
  pdf(pdf_file, width = 25, height = max(4, length(feature_order) * 0.08 + 2.5))
  draw(final_triple, main_heatmap = "Individual\nZ-score", 
       column_title = paste0("Triple Overlap Core Microbiota (n=", length(feature_order), ") - Individual Samples"),
       padding = unit(c(5, 5, 5, 40), "mm"))
  dev.off()
  message("✅ Triple Overlap 总结热图 (Individual) 已输出: ", pdf_file)
} else {
  message("⚠️ 未识别到 Triple Overlap 特征。")
}

message("\n🎉 全部处理完成！")
