# ==============================================================================
# 菌群 Meta 分析与可视化综合工具箱
# ==============================================================================

# --- 1. 加载必要的 R 包 ---
library(dplyr)          # 数据清洗与操作
library(tidyr)          # 数据转换 (replace_na, pivot)
library(tibble)         # 行名处理 (rownames_to_column)
library(purrr)          # 迭代与合并 (reduce, map)
library(ComplexHeatmap) # 核心热图绘制
library(circlize)       # 热图颜色映射 (colorRamp2)
library(grid)           # 绘图底层图形元素
library(ggplot2)        # 基础绘图 (轨迹图、森林图)
library(ggrepel)        # 标签防重叠
library(scales)         # 百分比与比例显示
library(VennDiagram)    # 韦恩图绘制
library(UpSetR)         # UpSet 图绘制
library(futile.logger)  # 日志管理 (用于静默 VennDiagram)

# 临时关闭 VennDiagram 的日志输出以保持控制台整洁
flog.threshold(ERROR)

# ==============================================================================
# --- 2. 基础数据清洗函数 ---
# ==============================================================================

# 处理 NA 的 P 值，将其视为不显著 (1)
clean_p <- function(p) ifelse(is.na(p), 1, p)

# 处理 NA 的系数 (Coefficient)，将其视为 0
clean_c <- function(c) ifelse(is.na(c), 0, c)

# ==============================================================================
# --- 3. 核心统计数据提取函数 ---
# ==============================================================================

# 从复杂的 Meta 分析结果列表中提取系数、P值和 FDR
extract_model_stats_final <- function(root_list, region_name, analysis_type, sub_group, comparison = NULL, features) {
  if (length(features) == 0) return(data.frame(feature = character(), coef = numeric(), pval = numeric(), qval.fdr = numeric()))
  
  res_df <- NULL
  if (analysis_type == "disease") {
    list_name <- if(!grepl("disease_analysis_", sub_group)) paste0("disease_analysis_", sub_group) else sub_group
    try({ res_df <- root_list[[region_name]][["disease_analysis"]][["analysis_results"]][[list_name]][["all_results"]][[comparison]][["meta_fits"]] }, silent = TRUE)
  } else if (analysis_type == "obesity") {
    list_name <- if(!grepl("obesity_analysis_", sub_group)) paste0("obesity_analysis_", sub_group) else sub_group
    try({ res_df <- root_list[[region_name]][["obesity_analysis"]][["analysis_results"]][[list_name]][["all_results"]][[comparison]][["meta_fits"]] }, silent = TRUE)
  } else if (analysis_type == "bmi") {
    try({ res_df <- root_list[[region_name]][["bmi_cor_analysis"]][["results_by_group"]][[sub_group]][["meta_result"]][["meta_fits"]] }, silent = TRUE)
  }
  
  if (is.null(res_df) || (is.data.frame(res_df) && nrow(res_df) == 0)) {
    return(data.frame(feature = features, coef = 0, pval = 1, qval.fdr = 1))
  }
  if (!is.data.frame(res_df)) res_df <- as.data.frame(res_df)
  if (!"feature" %in% colnames(res_df)) res_df$feature <- rownames(res_df)
  
  out <- data.frame(feature = features) %>%
    left_join(res_df %>% distinct(feature, .keep_all = TRUE), by = "feature") %>%
    select(feature, coef, pval, qval.fdr) %>%
    mutate(coef = replace_na(coef, 0), pval = replace_na(pval, 1), qval.fdr = replace_na(qval.fdr, 1))
  return(out)
}

# 针对 BMI 分层提取不同病程阶段的轨迹数据 (C-A, A-C, C-CRC)
extract_stratified_data <- function(bmi, comps, features) {
  results <- list()
  for(cp in comps) {
    bmi_key <- if(!grepl("disease_analysis_", bmi)) paste0("disease_analysis_", bmi) else bmi
    res <- tryCatch({
      all_analysis_results[["overall"]][["disease_analysis"]][["analysis_results"]][[bmi_key]][["all_results"]][[cp]][["meta_fits"]]
    }, error = function(e) NULL)
    
    if(!is.null(res)) {
      results[[cp]] <- res %>% filter(feature %in% features) %>%
        select(feature, coef, pval, qval.fdr) %>%
        rename(!!paste0(bmi, "_", cp, "_C") := coef, 
               !!paste0(bmi, "_", cp, "_P") := pval, 
               !!paste0(bmi, "_", cp, "_Q") := qval.fdr)
    } else {
      results[[cp]] <- data.frame(feature = features) %>%
        mutate(!!paste0(bmi, "_", cp, "_C") := NA, 
               !!paste0(bmi, "_", cp, "_P") := NA, !!paste0(bmi, "_", cp, "_Q") := NA)
    }
  }
  reduce(results, full_join, by="feature")
}

# ==============================================================================
# --- 4. 菌群分类与特征筛选函数 ---
# ==============================================================================

# 获取符合显著性阈值的特征列表
get_sig_feats <- function(data_source, region, bmi_group, comparison, p_cutoff, q_cutoff) {
  list_name <- paste0("disease_analysis_", bmi_group)
  df <- tryCatch({
    data_source[[region]][["disease_analysis"]][["analysis_results"]][[list_name]][["all_results"]][[comparison]][["meta_fits"]]
  }, error = function(e) return(NULL))
  
  if (is.null(df)) return(character(0))
  if (!"feature" %in% colnames(df)) df <- df %>% as.data.frame() %>% rownames_to_column("feature")
  
  sig_df <- df %>% filter(pval < 0.05) %>% filter(pval < p_cutoff | qval.fdr < q_cutoff)
  return(sig_df$feature)
}

# 将特征根据 BMI 组的变化方向进行分类
get_reg_cat_lists <- function(all_analysis_results, reg, comp_name, p_val, fdr_val) {
  s_norm <- get_sig_feats(all_analysis_results, reg, "Normal", comp_name, p_val, fdr_val)
  s_ob   <- get_sig_feats(all_analysis_results, reg, "Ob",     comp_name, p_val, fdr_val)
  all_sig <- unique(c(s_norm, s_ob))
  
  if(length(all_sig) == 0) return(list(uni=c(), norm=c(), ob=c()))
  
  d_n <- extract_model_stats_final(all_analysis_results, reg, "disease", "Normal", comp_name, all_sig)
  d_o <- extract_model_stats_final(all_analysis_results, reg, "disease", "Ob",     comp_name, all_sig)
  
  res <- data.frame(feature = all_sig) %>%
    mutate(c_n = d_n$coef[match(feature, d_n$feature)], c_o = d_o$coef[match(feature, d_o$feature)]) %>%
    mutate(cat = case_when(
      feature %in% s_norm & feature %in% s_ob ~ ifelse(sign(c_n) == sign(c_o), "Universal Core", "Normal-driven"),
      feature %in% s_norm ~ "Normal-driven",
      feature %in% s_ob   ~ "Obesity-specific",
      TRUE ~ "Other"
    ))
  
  return(list(uni=res$feature[res$cat=="Universal Core"], norm=res$feature[res$cat=="Normal-driven"], ob=res$feature[res$cat=="Obesity-specific"]))
}

# ==============================================================================
# --- 5. 逻辑修正与清洗函数 (用于 Venn 和 UpSet) ---
# ==============================================================================

# 根据系数方向清洗 Venn 集合，确保交集中的菌变化方向一致
clean_sets_by_direction <- function(s_norm, s_ow, s_ob, df_coef) {
  s_norm_f <- s_norm; s_ow_f <- s_ow; s_ob_f <- s_ob
  all_f <- unique(c(s_norm, s_ow, s_ob))
  for (f in all_f) {
    row <- df_coef[df_coef$feature == f, ]
    if (nrow(row) == 0) next
    # 若方向不一，从特定组剔除，使其在 Venn 图中显示为 Unique 而非 Overlap
    if (f %in% s_norm_f && f %in% s_ob_f && sign(row$coef_norm) != sign(row$coef_ob)) {
      if (row$coef_norm > 0) s_ob_f <- setdiff(s_ob_f, f) else s_norm_f <- setdiff(s_norm_f, f)
    }
  }
  return(list(norm = s_norm_f, ow = s_ow_f, ob = s_ob_f))
}

# ==============================================================================
# --- 6. 可视化辅助函数 ---
# ==============================================================================

# 将复杂的物种 ID (s__xxx|t__SGBxxx) 格式化为简称
format_sgb_names <- function(names_vec) {
  sapply(names_vec, function(x) {
    clean_name <- gsub(".*s__", "s__", x)
    parts <- strsplit(clean_name, "\\|")[[1]]
    sp_part <- gsub("^s__", "", parts[1])
    sgb_id <- if(length(parts) > 1) gsub("^t__", "", parts[2]) else ""
    res <- gsub("_", " ", sp_part)
    if (sgb_id != "" && !grepl(sgb_id, res)) res <- paste(res, sgb_id)
    return(res)
  })
}

# 准备热图矩阵
prep_mat <- function(dat_list, field, species_vector) {
  mat <- do.call(rbind, lapply(dat_list, function(x) x[[field]]))
  rownames(mat) <- names(dat_list)
  colnames(mat) <- species_vector 
  return(mat)
}

# 生成热图显著性星号标注层
make_layer_fun <- function(p_mat, q_mat) {
  function(j, i, x, y, width, height, fill) {
    v_p <- pindex(p_mat, i, j); v_q <- pindex(q_mat, i, j)
    is_sig <- (v_p < 0.05) & (v_p < CUSTOM_PVAL | v_q < CUSTOM_FDR)
    l <- ifelse(is_sig, ifelse(v_p < 0.01, "**", "*"), "")
    if(any(l != "")) grid.text(l[l != ""], x[l != ""], y[l != ""], gp = gpar(fontsize = 8))
  }
}

# 处理 NA 的 P 值，将其视为不显著 (1)
clean_p <- function(p) {
  ifelse(is.na(p), 1, p)
}

# 处理 NA 的系数 (Coefficient)，将其视为 0
clean_c <- function(c) {
  ifelse(is.na(c), 0, c)
}


get_reg_cat_lists <- function(all_analysis_results, reg, comp_name, p_val, fdr_val) {
  # 提取显著特征
  s_norm <- get_sig_feats(all_analysis_results, reg, "Normal", comp_name, p_val, fdr_val)
  s_ob   <- get_sig_feats(all_analysis_results, reg, "Ob",     comp_name, p_val, fdr_val)
  all_sig <- unique(c(s_norm, s_ob))
  
  if(length(all_sig) == 0) return(list(uni=c(), norm=c(), ob=c()))
  
  # 提取系数进行方向判定
  d_n <- extract_model_stats_final(all_analysis_results, reg, "disease", "Normal", comp_name, all_sig)
  d_o <- extract_model_stats_final(all_analysis_results, reg, "disease", "Ob",     comp_name, all_sig)
  
  df_check <- data.frame(
    feature = all_sig,
    c_n = d_n$coef[match(all_sig, d_n$feature)],
    c_o = d_o$coef[match(all_sig, d_o$feature)]
  )
  
  # 分类归纳
  res <- df_check %>%
    mutate(cat = case_when(
      feature %in% s_norm & feature %in% s_ob ~ ifelse(sign(c_n) == sign(c_o), "Universal Core", 
                                                       ifelse(c_n > 0, "Normal-driven", "Obesity-specific")),
      feature %in% s_norm ~ "Normal-driven",
      feature %in% s_ob   ~ "Obesity-specific",
      TRUE ~ "Other"
    ))
  
  return(list(
    uni  = res$feature[res$cat == "Universal Core"],
    norm = res$feature[res$cat == "Normal-driven"],
    ob   = res$feature[res$cat == "Obesity-specific"]
  ))
}

extract_full_metrics_multi <- function(data_source, region, bmi_group, comparison, target_features) {
  list_name <- paste0("disease_analysis_", bmi_group)
  df <- tryCatch({
    data_source[[region]][["disease_analysis"]][["analysis_results"]][[list_name]][["all_results"]][[comparison]][["meta_fits"]]
  }, error = function(e) return(NULL))
  
  if (is.null(df)) return(NULL)
  if (!"feature" %in% colnames(df)) df <- df %>% as.data.frame() %>% rownames_to_column("feature")
  
  # 提取目标特征，并重命名列 (添加 region 和 bmi_group 后缀)
  df %>%
    filter(feature %in% target_features) %>%
    select(feature, coef, pval, qval.fdr) %>%
    rename_with(~paste0(., "_", region, "_", bmi_group), -feature)
}

# [cite_start]从 Meta 分析结果中提取系数、P值和 FDR [cite: 168, 250, 251]
extract_model_stats_final <- function(root_list, region_name, analysis_type, sub_group, comparison = NULL, features) {
  res_df <- NULL
  if (analysis_type == "disease") {
    list_name <- if(!grepl("disease_analysis_", sub_group)) paste0("disease_analysis_", sub_group) else sub_group
    try({ res_df <- root_list[[region_name]][["disease_analysis"]][["analysis_results"]][[list_name]][["all_results"]][[comparison]][["meta_fits"]] }, silent = TRUE)
  } else if (analysis_type == "obesity") {
    list_name <- if(!grepl("obesity_analysis_", sub_group)) paste0("obesity_analysis_", sub_group) else sub_group
    try({ res_df <- root_list[[region_name]][["obesity_analysis"]][["analysis_results"]][[list_name]][["all_results"]][[comparison]][["meta_fits"]] }, silent = TRUE)
  } else if (analysis_type == "bmi") {
    try({ res_df <- root_list[[region_name]][["bmi_cor_analysis"]][["results_by_group"]][[sub_group]][["meta_result"]][["meta_fits"]] }, silent = TRUE)
  }
  
  if (is.null(res_df) || (is.data.frame(res_df) && nrow(res_df) == 0)) {
    return(data.frame(feature = features, coef = 0, pval = 1, qval.fdr = 1))
  }
  if (!is.data.frame(res_df)) res_df <- as.data.frame(res_df)
  if (!"feature" %in% colnames(res_df)) res_df$feature <- rownames(res_df)
  
  out <- data.frame(feature = features) %>%
    left_join(res_df %>% distinct(feature, .keep_all = TRUE), by = "feature") %>%
    select(feature, coef, pval, qval.fdr) %>%
    mutate(coef = replace_na(coef, 0), pval = replace_na(pval, 1), qval.fdr = replace_na(qval.fdr, 1))
  return(out)
}

# [cite_start]针对 BMI 分层提取轨迹数据 [cite: 251]
extract_stratified_data <- function(bmi, comps, features) {
  results <- list()
  for(cp in comps) {
    bmi_key <- if(!grepl("disease_analysis_", bmi)) paste0("disease_analysis_", bmi) else bmi
    res <- tryCatch({
      all_analysis_results[["overall"]][["disease_analysis"]][["analysis_results"]][[bmi_key]][["all_results"]][[cp]][["meta_fits"]]
    }, error = function(e) NULL)
    
    if(!is.null(res)) {
      results[[cp]] <- res %>% 
        filter(feature %in% features) %>%
        select(feature, coef, pval, qval.fdr) %>%
        rename(!!paste0(bmi, "_", cp, "_C") := coef, 
               !!paste0(bmi, "_", cp, "_P") := pval, 
               !!paste0(bmi, "_", cp, "_Q") := qval.fdr)
    } else {
      results[[cp]] <- data.frame(feature = features) %>%
        mutate(!!paste0(bmi, "_", cp, "_C") := NA, 
               !!paste0(bmi, "_", cp, "_P") := NA, 
               !!paste0(bmi, "_", cp, "_Q") := NA)
    }
  }
  reduce(results, full_join, by="feature")
}

# [cite_start]提取森林图所需的单个研究 (Study) 统计量 [cite: 10]
extract_study_stats <- function(study_name, fits_list, target_bugs_vector) {
  obj <- fits_list[[study_name]]
  res_df <- if ("results" %in% names(obj)) obj$results else obj
  
  target_rows <- res_df %>% filter(feature %in% target_bugs_vector) 
  if (nrow(target_rows) == 0) return(NULL)
  
  data.frame(
    Study   = study_name,
    Feature = target_rows$feature,
    Beta    = target_rows$coef,
    SE      = target_rows$stderr,
    Pval    = target_rows$pval,
    N       = target_rows$N,
    Type    = "Study"
  )
}

# [cite_start]获取显著特征列表 [cite: 168]
get_sig_feats <- function(data_source, region, bmi_group, comparison, p_cutoff, q_cutoff) {
  list_name <- paste0("disease_analysis_", bmi_group)
  df <- tryCatch({
    data_source[[region]][["disease_analysis"]][["analysis_results"]][[list_name]][["all_results"]][[comparison]][["meta_fits"]]
  }, error = function(e) return(NULL))
  
  if (is.null(df)) return(character(0))
  if (!"feature" %in% colnames(df)) df <- df %>% as.data.frame() %>% rownames_to_column("feature")
  
  sig_df <- df %>%
    filter(pval < 0.05) %>%
    filter(pval < p_cutoff | qval.fdr < q_cutoff)
  return(sig_df$feature)
}

# [cite_start]将特征分类 (Universal Core, Normal-driven, Obesity-specific) [cite: 251]
get_reg_cat_lists <- function(all_analysis_results, reg, comp_name, p_val, fdr_val) {
  s_norm <- get_sig_feats(all_analysis_results, reg, "Normal", comp_name, p_val, fdr_val)
  s_ob   <- get_sig_feats(all_analysis_results, reg, "Ob",     comp_name, p_val, fdr_val)
  all_sig <- unique(c(s_norm, s_ob))
  
  if(length(all_sig) == 0) return(list(uni=c(), norm=c(), ob=c()))
  
  d_n <- extract_model_stats_final(all_analysis_results, reg, "disease", "Normal", comp_name, all_sig)
  d_o <- extract_model_stats_final(all_analysis_results, reg, "disease", "Ob",     comp_name, all_sig)
  
  df_check <- data.frame(
    feature = all_sig,
    c_n = d_n$coef[match(all_sig, d_n$feature)],
    c_o = d_o$coef[match(all_sig, d_o$feature)]
  )
  
  res <- df_check %>%
    mutate(cat = case_when(
      feature %in% s_norm & feature %in% s_ob ~ ifelse(sign(c_n) == sign(c_o), "Universal Core", 
                                                       ifelse(c_n > 0, "Normal-driven", "Obesity-specific")),
      feature %in% s_norm ~ "Normal-driven",
      feature %in% s_ob   ~ "Obesity-specific",
      TRUE ~ "Other"
    ))
  
  return(list(
    uni  = res$feature[res$cat == "Universal Core"],
    norm = res$feature[res$cat == "Normal-driven"],
    ob   = res$feature[res$cat == "Obesity-specific"]
  ))
}

# [cite_start]将 SGB ID 格式化为更具可读性的种名 [cite: 168, 251]
format_sgb_names <- function(names_vec) {
  sapply(names_vec, function(x) {
    # 1. 预处理：去掉路径前缀，仅保留 s__ 开始的部分
    clean_name <- gsub(".*s__", "s__", x)
    
    # 2. 分解 species 部分和 strain/SGB 部分
    parts <- strsplit(clean_name, "\\|")[[1]]
    sp_part <- gsub("^s__", "", parts[1]) # Genus_species...
    sgb_id <- if(length(parts) > 1) gsub("^t__", "", parts[2]) else ""
    
    # 逻辑 A：处理包含 _sp_ 的情况 (如 Clostridium_sp_AF34_13)
    if (grepl("_sp_", sp_part)) {
      res <- gsub("_sp_", " sp. ", sp_part)
      res <- gsub("_", " ", res)
    } 
    # 逻辑 B：处理种名本身就是 SGB 的情况 (如 Clostridium_SGB4751)
    else if (grepl("SGB", sp_part)) {
      res <- gsub("_", " ", sp_part)
    } 
    # 逻辑 C：标准属名_种名 (如 Blautia_wexlerae) -> B. wexlerae
    else {
      sub_parts <- strsplit(sp_part, "_")[[1]]
      if (length(sub_parts) >= 2) {
        genus <- sub_parts[1]
        species <- paste(sub_parts[-1], collapse = " ")
        genus_abbr <- paste0(substr(genus, 1, 1), ".")
        res <- paste(genus_abbr, species)
      } else {
        res <- gsub("_", " ", sp_part)
      }
    }
    
    # 合并结果，并去除可能的重复 SGB ID (针对 s__...SGBxxx|t__SGBxxx 情况)
    final_res <- trimws(res)
    if (sgb_id != "" && !grepl(sgb_id, final_res)) {
      final_res <- paste(final_res, sgb_id)
    }
    return(final_res)
  })
}

get_robust_sig_feats_rescue <- function(data_source, region, bmi_group, comp, p_thresh, fdr_thresh, weight_limit = 0.7, min_studies = 2) {
  list_name <- paste0("disease_analysis_", bmi_group)
  df_meta <- data_source[[region]][["disease_analysis"]][["analysis_results"]][[list_name]][["all_results"]][[comp]][["meta_fits"]]
  df_maaslin <- data_source[[region]][["disease_analysis"]][["analysis_results"]][[list_name]][["all_results"]][[comp]][["maaslin_fits"]]
  
  if (is.null(df_meta) | is.null(df_maaslin)) return(character(0))
  
  # 权重提取
  weight_cols <- grep("^weight_", colnames(df_meta), value = TRUE)
  w_mat <- as.matrix(df_meta[, weight_cols])
  w_mat[is.na(w_mat)] <- 0
  
  # 内部拯救函数：自动剔除最大权重研究并重算 P 值
  perform_loo_rescue <- function(f_name, weight_row) {
    # 找到该菌权重最大的研究名
    max_idx <- which.max(weight_row)
    exclude_study <- gsub("^weight_", "", names(weight_row)[max_idx])
    
    remaining_studies <- names(df_maaslin)[names(df_maaslin) != exclude_study]
    coefs <- c(); ses <- c()
    
    for (s in remaining_studies) {
      df_s <- df_maaslin[[s]]
      # 【核心修复】：通过 feature 列匹配，不依赖 rownames
      if ("feature" %in% colnames(df_s) && f_name %in% df_s$feature) {
        idx <- which(df_s$feature == f_name)
        val_c <- df_s$coef[idx]; val_se <- df_s$stderr[idx]
        if (!is.na(val_c) && !is.na(val_se)) { coefs <- c(coefs, val_c); ses <- c(ses, val_se) }
      }
    }
    
    if (length(coefs) < 1) return(FALSE) # 孤证不立
    
    # 简单的逆方差汇总 (LOO Meta)
    weights <- 1 / (ses^2)
    combined_coef <- sum(coefs * weights) / sum(weights)
    combined_se <- sqrt(1 / sum(weights))
    combined_p <- 2 * (1 - pnorm(abs(combined_coef / combined_se)))
    
    return(combined_p < 0.05) # 拯救标准：剔除后依然显著
  }
  
  audit_df <- data.frame(
    feature = rownames(df_meta),
    Max_Weight = apply(w_mat, 1, function(x) { s=sum(x); if(s<=0) return(1) else max(x)/s }),
    pval = df_meta$pval,
    qval = df_meta$qval.fdr
  )
  
  pass_features <- c()
  for (i in 1:nrow(audit_df)) {
    f <- as.character(audit_df$feature[i])
    # 1. 必须显著
    is_sig <- !is.na(audit_df$pval[i]) && (audit_df$pval[i] < 0.05) && (audit_df$pval[i] < p_thresh | audit_df$qval[i] < fdr_thresh)
    if (!is_sig) next
    
    # 2. 稳健性判断
    if (audit_df$Max_Weight[i] < weight_limit) {
      pass_features <- c(pass_features, f) # 权重分布均匀，直接通过
    } else if (perform_loo_rescue(f, w_mat[f, ])) {
      pass_features <- c(pass_features, f) # 权重虽大但通过拯救测试，找回！
      message(">>> Rescued: ", f)
    }
  }
  return(pass_features)
}

# [cite_start]准备热图所需的矩阵 [cite: 168]
prep_mat <- function(layer_list, value_col, features) {
  mat_list <- lapply(layer_list, function(df) {
    res <- df[match(features, df$feature), value_col]
    return(res)
  })
  mat <- do.call(cbind, mat_list) # 结果: 行=Features, 列=Groups
  rownames(mat) <- features
  return(mat)
}

format_sgb_names <- function(names_vec) {
  sapply(names_vec, function(x) {
    # 1. 预处理：去掉路径前缀，仅保留 s__ 开始的部分
    clean_name <- gsub(".*s__", "s__", x)
    
    # 2. 分解 species 部分和 strain/SGB 部分
    parts <- strsplit(clean_name, "\\|")[[1]]
    sp_part <- gsub("^s__", "", parts[1]) # Genus_species...
    sgb_id <- if(length(parts) > 1) gsub("^t__", "", parts[2]) else ""
    
    # 逻辑 A：处理包含 _sp_ 的情况 (如 Clostridium_sp_AF34_13)
    if (grepl("_sp_", sp_part)) {
      res <- gsub("_sp_", " sp. ", sp_part)
      res <- gsub("_", " ", res)
    } 
    # 逻辑 B：处理种名本身就是 SGB 的情况 (如 Clostridium_SGB4751)
    else if (grepl("SGB", sp_part)) {
      res <- gsub("_", " ", sp_part)
    } 
    # 逻辑 C：标准属名_种名 (如 Blautia_wexlerae) -> B. wexlerae
    else {
      sub_parts <- strsplit(sp_part, "_")[[1]]
      if (length(sub_parts) >= 2) {
        genus <- sub_parts[1]
        species <- paste(sub_parts[-1], collapse = " ")
        genus_abbr <- paste0(substr(genus, 1, 1), ".")
        res <- paste(genus_abbr, species)
      } else {
        res <- gsub("_", " ", sp_part)
      }
    }
    
    # 合并结果，并去除可能的重复 SGB ID (针对 s__...SGBxxx|t__SGBxxx 情况)
    final_res <- trimws(res)
    if (sgb_id != "" && !grepl(sgb_id, final_res)) {
      final_res <- paste(final_res, sgb_id)
    }
    return(final_res)
  })
}

# [cite_start]生成热图单元格标注函数 (* 或 **) [cite: 168]
make_layer_fun <- function(p_mat, q_mat) {
  force(p_mat); force(q_mat)
  function(j, i, x, y, width, height, fill) {
    v_p <- pindex(p_mat, i, j)
    v_q <- pindex(q_mat, i, j)
    is_sig <- (v_p < 0.05) & (v_p < CUSTOM_PVAL | v_q < CUSTOM_FDR)
    l <- ifelse(is_sig, ifelse(v_p < 0.01, "**", "*"), "")
    if(any(l != "")) grid.text(l[l != ""], x[l != ""], y[l != ""], gp = gpar(fontsize = 8))
  }
}

# [cite_start]提取轨迹绘图统计量 [cite: 251]
get_trajectory_stats <- function(feature_list) {
  bmi_groups <- c("Normal", "Ob")
  comparisons <- c("CTRL_vs_AA", "CTRL_vs_CRC") 
  results <- list()
  for(bmi in bmi_groups) {
    results[[paste0(bmi, "_base")]] <- data.frame(BMI=bmi, Stage="CTRL", Coef=0, SE=0, feature=feature_list)
    for(comp in comparisons) {
      bmi_key <- paste0("disease_analysis_", bmi)
      res <- all_analysis_results[["overall"]][["disease_analysis"]][["analysis_results"]][[bmi_key]][["all_results"]][[comp]][["meta_fits"]]
      if(!is.null(res)) {
        tmp <- res %>% filter(feature %in% feature_list) %>% 
          select(feature, coef, stderr) %>%
          mutate(BMI=bmi, Stage=ifelse(comp=="CTRL_vs_AA", "AA", "CRC")) %>%
          rename(Coef=coef, SE=stderr)
        results[[paste0(bmi, "_", comp)]] <- tmp
      }
    }
  }
  bind_rows(results)
}

# 提取全量指标并自动生成带后缀的列名
extract_full_metrics_multi <- function(data_source, region, bmi_group, comparison, target_features) {
  list_name <- paste0("disease_analysis_", bmi_group)
  df <- tryCatch({
    data_source[[region]][["disease_analysis"]][["analysis_results"]][[list_name]][["all_results"]][[comparison]][["meta_fits"]]
  }, error = function(e) return(NULL))
  
  if (is.null(df)) return(NULL)
  if (!"feature" %in% colnames(df)) df <- df %>% as.data.frame() %>% tibble::rownames_to_column("feature")
  
  # 提取目标特征，并重命名列 (添加 region 和 bmi_group 后缀)
  df %>%
    filter(feature %in% target_features) %>%
    select(feature, coef, pval, qval.fdr) %>%
    rename_with(~paste0(., "_", region, "_", bmi_group), -feature)
}

# 为 UpSet 图准备分类好的特征数据
get_categorized_features <- function(region, comp_name) {
  # 调用之前提炼的 get_sig_feats
  s_norm <- get_sig_feats(all_analysis_results, region, "Normal", comp_name, CUSTOM_PVAL, CUSTOM_FDR)
  s_ob   <- get_sig_feats(all_analysis_results, region, "Ob",     comp_name, CUSTOM_PVAL, CUSTOM_FDR)
  all_feats <- unique(c(s_norm, s_ob))
  
  if(length(all_feats) == 0) return(NULL)
  
  # 提取系数用于方向判断
  df_norm <- extract_model_stats_final(all_analysis_results, region, "disease", "Normal", comp_name, all_feats)
  df_ob   <- extract_model_stats_final(all_analysis_results, region, "disease", "Ob",     comp_name, all_feats)
  
  dir_check <- data.frame(feature = all_feats) %>%
    left_join(df_norm %>% select(feature, coef_norm = coef), by = "feature") %>%
    left_join(df_ob %>% select(feature, coef_ob = coef), by = "feature") %>%
    mutate(coef_norm = tidyr::replace_na(coef_norm, 0),
           coef_ob   = tidyr::replace_na(coef_ob, 0))
  
  # 应用热图分类逻辑
  categorized <- dir_check %>%
    mutate(Category = case_when(
      feature %in% s_norm & feature %in% s_ob ~ case_when(
        sign(coef_norm) == sign(coef_ob) ~ "Universal Core",
        coef_norm > 0 ~ "Normal-driven",
        coef_ob > 0   ~ "Obesity-specific",
        TRUE ~ "Normal-driven"
      ),
      feature %in% s_norm ~ "Normal-driven",
      feature %in% s_ob   ~ "Obesity-specific",
      TRUE ~ "Other"
    ))
  return(categorized)
}
# 根据系数方向清洗集合，剔除方向不一致的虚假重叠
clean_sets_by_direction <- function(s_norm, s_ow, s_ob, df_coef) {
  s_norm_final <- s_norm; s_ow_final <- s_ow; s_ob_final <- s_ob
  all_feats <- unique(c(s_norm, s_ow, s_ob))
  
  for (f in all_feats) {
    row_idx <- which(df_coef$feature == f)
    if (length(row_idx) == 0) next
    c_norm <- df_coef$coef_norm[row_idx]; c_ow <- df_coef$coef_ow[row_idx]; c_ob <- df_coef$coef_ob[row_idx]
    
    # Norm vs Ob 冲突处理
    if (f %in% s_norm_final && f %in% s_ob_final) {
      if (sign(c_norm) != sign(c_ob) && c_norm != 0 && c_ob != 0) {
        if (c_norm > 0) s_ob_final <- setdiff(s_ob_final, f) else s_norm_final <- setdiff(s_norm_final, f)
      }
    }
    # Norm vs Ow 冲突处理
    if (f %in% s_norm_final && f %in% s_ow_final) {
      if (sign(c_norm) != sign(c_ow) && c_norm != 0 && c_ow != 0) {
        if (c_norm > 0) s_ow_final <- setdiff(s_ow_final, f) else s_norm_final <- setdiff(s_norm_final, f)
      }
    }
    # Ow vs Ob 冲突处理
    if (f %in% s_ow_final && f %in% s_ob_final) {
      if (sign(c_ow) != sign(c_ob) && c_ow != 0 && c_ob != 0) {
        if (c_ow > 0) s_ob_final <- setdiff(s_ob_final, f) else s_ow_final <- setdiff(s_ow_final, f)
      }
    }
  }
  return(list(norm = s_norm_final, ow = s_ow_final, ob = s_ob_final))
}

# 获取跨地区交集且变化方向完全一致的特征
get_consistent_overlap <- function(comp_name, cat_name, region_category_data) {
  # 提取三地的该类别数据列表
  list_3reg <- lapply(region_category_data[[comp_name]], `[[`, cat_name)
  
  # 1. 取 feature 的交集
  common_feats <- Reduce(intersect, lapply(list_3reg, function(x) x$feature))
  if(length(common_feats) == 0) return(data.frame())
  
  # 2. 检查跨地区方向一致性 (所有地区方向必须跟第一个地区 overall 一致)
  consistent_mask <- sapply(common_feats, function(f) {
    dirs <- sapply(list_3reg, function(df) sign(df$coef[df$feature == f]))
    all(dirs == dirs[1]) 
  })
  
  final_feats <- common_feats[consistent_mask]
  return(list_3reg$overall %>% filter(feature %in% final_feats))
}

