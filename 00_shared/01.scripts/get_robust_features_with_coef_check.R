# ==============================================================================
# 稳健性审计函数 (增强版：检查系数变化与符号反转)
# ==============================================================================
get_robust_features_with_coef_check <- function(data_source, region, bmi_group, comp, 
                                                p_thresh, fdr_thresh, 
                                                weight_limit = 0.7) {
  
  list_name <- paste0("disease_analysis_", bmi_group)
  # 增加安全性检查
  res_obj <- tryCatch(data_source[[region]][["disease_analysis"]][["analysis_results"]][[list_name]][["all_results"]][[comp]], error=function(e) NULL)
  if(is.null(res_obj)) return(NULL)
  
  df_meta <- res_obj[["meta_fits"]]
  df_maaslin <- res_obj[["maaslin_fits"]]
  
  if (is.null(df_meta) | is.null(df_maaslin)) return(NULL)
  
  # 权重提取
  weight_cols <- grep("^weight_", colnames(df_meta), value = TRUE)
  w_mat <- as.matrix(df_meta[, weight_cols])
  w_mat[is.na(w_mat)] <- 0
  
  # --- 内部函数：LOO 计算并返回详细数据 ---
  calc_loo_stats <- function(f_name, weight_row) {
    max_idx <- which.max(weight_row)
    if(length(max_idx)==0) return(list(coef=NA, pval=NA, note="No Weight"))
    
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
    
    if (length(coefs) < 1) return(list(coef=NA, pval=NA, note="Too few studies"))
    
    # 逆方差加权计算
    weights <- 1 / (ses^2)
    sum_w <- sum(weights)
    if(sum_w == 0) return(list(coef=NA, pval=NA, note="Zero Weight"))
    
    combined_coef <- sum(coefs * weights) / sum_w
    combined_se <- sqrt(1 / sum_w)
    combined_p <- 2 * (1 - pnorm(abs(combined_coef / combined_se)))
    
    return(list(coef = combined_coef, pval = combined_p, note = "Calculated"))
  }
  
  # --- 主循环：构建详细报表 ---
  # 先筛选出原始显著的候选菌
  candidates <- rownames(df_meta)
  
  results_list <- list()
  
  for (f in candidates) {
    row <- df_meta[f, ]
    
    # 1. 基础显著性检查
    is_sig <- !is.na(row$pval) && (row$pval < 0.05) && (row$pval < p_thresh | row$qval.fdr < fdr_thresh)
    if (!is_sig) next # 不显著的直接跳过，不浪费计算资源
    
    # 计算最大权重
    w_row <- w_mat[f, ]
    max_w <- if(sum(w_row)<=0) 1 else max(w_row)/sum(w_row)
    
    # 2. 判断是否需要拯救
    orig_coef <- row$coef
    loo_coef <- NA
    loo_p <- NA
    status <- "Discard"
    note <- ""
    
    if (max_w < weight_limit) {
      # 权重分布健康，不需要 LOO，直接通过
      status <- "Robust"
      loo_coef <- orig_coef # 此时 LOO Coef 视为与原始一致 (或者 NA)
      loo_p <- row$pval
      note <- "Low Weight Ratio"
    } else {
      # 权重过高，执行 LOO
      loo_res <- calc_loo_stats(f, w_row)
      loo_coef <- loo_res$coef
      loo_p <- loo_res$pval
      
      if (!is.na(loo_p) && loo_p < 0.05) {
        status <- "Rescued"
        note <- "Significant after LOO"
      } else {
        status <- "Discard"
        note <- "Not significant after LOO"
      }
    }
    
    # 3. 符号反转检查 (Sign Flip Check)
    # 只有当两个 Coef 都有值时才检查
    sign_flipped <- FALSE
    if (!is.na(orig_coef) && !is.na(loo_coef)) {
      if (sign(orig_coef) != sign(loo_coef)) {
        sign_flipped <- TRUE
      }
    }
    
    # 如果通过了筛选，但符号反转了，这非常危险，建议标记或剔除
    if (status %in% c("Robust", "Rescued") && sign_flipped && max_w >= weight_limit) {
      # 注意：如果是 Robust (权重低)，通常 coef 不会变。
      # 只有 Rescued (去掉了大权重) 的时候，符号反转才是有意义的风险。
      note <- paste0(note, " [WARNING: Sign Flipped]")
    }
    
    results_list[[f]] <- data.frame(
      feature = f,
      Original_P = row$pval,
      Original_Coef = orig_coef,
      Max_Weight = max_w,
      LOO_P = loo_p,
      LOO_Coef = loo_coef,
      Sign_Flipped = sign_flipped, # TRUE/FALSE
      Status = status,
      Note = note,
      stringsAsFactors = FALSE
    )
  }
  
  if(length(results_list) == 0) return(NULL)
  return(do.call(rbind, results_list))
}

# ==============================================================================
# 4. 运行并检查
# ==============================================================================
# 运行新函数
df_check_norm <- get_robust_features_with_coef_check(all_analysis_results, reg, "Normal", comp, CUSTOM_PVAL, CUSTOM_FDR)
df_check_ob   <- get_robust_features_with_coef_check(all_analysis_results, reg, "Ob",     comp, CUSTOM_PVAL, CUSTOM_FDR)

# 查看那些 "被拯救" (Rescued) 的菌，重点看 Coef 变化
message("--- Normal 组被拯救的菌 (检查 Coef 变化) ---")
if(!is.null(df_check_norm)) {
  rescued_norm <- df_check_norm %>% filter(Status == "Rescued")
  print(rescued_norm %>% select(feature, Original_Coef, LOO_Coef, Sign_Flipped, Note))
}

message("\n--- Obese 组被拯救的菌 (检查 Coef 变化) ---")
if(!is.null(df_check_ob)) {
  rescued_ob <- df_check_ob %>% filter(Status == "Rescued")
  print(rescued_ob %>% select(feature, Original_Coef, LOO_Coef, Sign_Flipped, Note))
}

# 自动获取过滤后的最终名单 (排除掉符号反转的)
# 逻辑：必须是 Robust 或 Rescued，且如果是 Rescued，不能发生符号反转
get_clean_features <- function(df) {
  if(is.null(df)) return(c())
  df %>% 
    filter(Status %in% c("Robust", "Rescued")) %>%
    filter(!(Status == "Rescued" & Sign_Flipped == TRUE)) %>% # 剔除拯救后反转的
    pull(feature)
}

s_norm_final <- get_clean_features(df_check_norm)
s_ob_final   <- get_clean_features(df_check_ob)

message(paste("\n最终保留数量: Normal =", length(s_norm_final), ", Ob =", length(s_ob_final)))