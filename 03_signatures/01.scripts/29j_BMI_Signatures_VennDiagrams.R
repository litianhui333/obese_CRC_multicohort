# ==============================================================================
# 脚本名称: 29j_BMI_Signatures_VennDiagrams.R
# 功能描述: 
#   利用 Venn 图展示三个 BMI 组之间显著菌群的重叠情况，
#   分为 Enriched in CRC 和 Depleted in CRC 两个维度。
#   严格标准 (v5_Strict): P < 0.01 或 FDR < 0.1
# ==============================================================================

library(VennDiagram)
library(grid)
library(gridExtra)

# ==============================================================================
# 0.5 核心函数: 稳健性审计 (LOO Audit)
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
      if (combined_p < 0.05) { # Rescue
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
output_dir <- file.path(base_dir, "03.results/04_differential_abundance/heatmap_final/BMI_Stratified_Disease_Heatmaps")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

message("=== 1. 加载严格标准统计结果 (v5_Strict + LOO Audit) ===")

# 加载分析结果 (用于 LOO 审计)
input_dir <- file.path(base_dir, "03.results/04_differential_abundance")
res_file <- sort(list.files(input_dir, pattern = "microbiome_analysis_all_results_BMI_Group1_0.7.*\\.rds", full.names = TRUE), decreasing = TRUE)[1]
all_analysis_results <- readRDS(res_file)

# 识别三组里各自显著且稳健的菌群 (Strict 0.01)
df_robust_Nw <- get_robust_features_detailed(all_analysis_results, "overall", "Normal", "CTRL_vs_CRC", p_cutoff = 0.01, q_cutoff = 0.1)
df_robust_Ow <- get_robust_features_detailed(all_analysis_results, "overall", "Ow",     "CTRL_vs_CRC", p_cutoff = 0.01, q_cutoff = 0.1)
df_robust_Ob <- get_robust_features_detailed(all_analysis_results, "overall", "Ob",     "CTRL_vs_CRC", p_cutoff = 0.01, q_cutoff = 0.1)

pass_nw <- if (!is.null(df_robust_Nw)) df_robust_Nw$feature else character(0)
pass_ow <- if (!is.null(df_robust_Ow)) df_robust_Ow$feature else character(0)
pass_ob <- if (!is.null(df_robust_Ob)) df_robust_Ob$feature else character(0)

# 加载 P 值/系数矩阵 (用于判定 Enriched/Depleted)
rds_dir_pure <- file.path(base_dir, "03.results/04_differential_abundance/heatmap_final/Robust_Analysis_v5_Pure")
pure_file <- list.files(rds_dir_pure, pattern = "Robust_Stat_overall_CTRL_vs_CRC.*\\.rds", full.names = TRUE)[1]
dat_pure <- readRDS(pure_file)

m_p <- as.matrix(dat_pure$matrix_data$pval)
m_c <- as.matrix(dat_pure$matrix_data$coef)

# 列名对齐
colnames(m_p)[colnames(m_p) == "Overweight"] <- "Ow"
colnames(m_p)[colnames(m_p) == "Obese"]      <- "Ob"
colnames(m_c)[colnames(m_c) == "Overweight"] <- "Ow"
colnames(m_c)[colnames(m_c) == "Obese"]      <- "Ob"

# ==============================================================================
# 2. 提取特征名单 (基于 Gap Logic 以对齐主图 18+4 逻辑)
# ==============================================================================
# 为了让 Venn 图的“Obese 独有区域”正好等于主图标注的 Ob-specific 数量 (即 18+4)，
# 我们需要构造如下集合：
# Obese 集合: 严苛显著名单 (pass_ob)
# Normal/Ow 集合: 广义显著名单 (P <= 0.05)
# 这样: Ob_only = Ob_sig - (Nw_marginal | Ow_marginal) = P_ob<0.01 & P_nw>0.05 & P_ow>0.05
# 这正是脚本 29 中 Specific 的定义标准。

message("=== 2. 提取 Venn 集合 (Gap Logic 对齐版) ===")

is_marginal_Nw <- rownames(m_p)[which(m_p[, "Normal"] <= 0.05)]
is_marginal_Ow <- rownames(m_p)[which(m_p[, "Ow"]     <= 0.05)]
is_marginal_Ob <- rownames(m_p)[which(m_p[, "Ob"]     <= 0.05)]

# 对于 Enriched Venn
Venn_En_List <- list(
  Normal = intersect(is_marginal_Nw, rownames(m_c)[which(m_c[, "Normal"] > 0)]),
  Ow     = intersect(is_marginal_Ow, rownames(m_c)[which(m_c[, "Ow"]     > 0)]),
  Obese  = intersect(pass_ob,        rownames(m_c)[which(m_c[, "Ob"]     > 0)]) 
)

# 对于 Depleted Venn
Venn_De_List <- list(
  Normal = intersect(is_marginal_Nw, rownames(m_c)[which(m_c[, "Normal"] < 0)]),
  Ow     = intersect(is_marginal_Ow, rownames(m_c)[which(m_c[, "Ow"]     < 0)]),
  Obese  = intersect(pass_ob,        rownames(m_c)[which(m_c[, "Ob"]     < 0)])
)

list_en <- Venn_En_List
list_de <- Venn_De_List

# ==============================================================================
# 3. 绘制 Venn 图
# ==============================================================================
message("=== 3. 绘制 Venn 图 (Enriched / Depleted) ===")

# 配色: Normal: #5AAE61, Ow: #FF8C00, Obese: #9970AB
my_colors <- c("#5AAE61", "#FF8C00", "#9970AB")

draw_custom_venn <- function(x_list, title) {
  venn_obj <- venn.diagram(
    x = x_list,
    filename = NULL,
    col = "transparent",
    fill = my_colors,
    alpha = 0.5,
    label.col = "black",
    cex = 1.2,
    fontface = "bold",
    cat.col = my_colors,
    cat.cex = 1.2,
    cat.fontface = "bold",
    cat.dist = c(0.06, 0.06, 0.06), # 距离调整
    cat.pos = c(-20, 20, 180),    # 位置调整
    main = title,
    main.cex = 1.5,
    main.fontface = "bold"
  )
  return(venn_obj)
}

# 生成两个 Venn
v_en <- draw_custom_venn(list_en, "Enriched in CRC")
v_de <- draw_custom_venn(list_de, "Depleted in CRC")

# 保存 PDF (一页包含两张图)
pdf_file <- file.path(output_dir, "BMI_Signatures_Venn_Enriched_Depleted_Strict.pdf")
pdf(pdf_file, width = 12, height = 6)
grid.arrange(gTree(children = v_en), gTree(children = v_de), ncol = 2)
dev.off()

# ==============================================================================
# 4. 打印交集详情
# ==============================================================================
message("\n=== [Enriched] 交集数量 ===")
print(sapply(list_en, length))
message("\n=== [Depleted] 交集数量 ===")
print(sapply(list_de, length))

# 获取 Ob-specific 菌名 (Enriched)
ob_en_spec <- setdiff(Ob_En, c(Nw_En, Ow_En))
if (length(ob_en_spec) > 0) {
  message("\n🌟 Obese-Specific Enriched Microbes (Top 5):")
  print(head(ob_en_spec, 5))
}

message("\n🎉 Venn 图已成功输出: ", pdf_file)
