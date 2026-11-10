# ==============================================================================
# 脚本名称: alpha_analysis_functions.R
# 功能: 定义 Alpha 多样性分析的核心函数 (Meta-analysis & LMM)
# 更新: 
# 1. 强制统一返回单一数据框 (Data Frame)，解决 bind_rows 报错问题
# 2. 支持 Unadjusted 和 Adjusted (Age+Sex) 两种模式
# 3. Adjusted Meta 采用 Two-stage IPD 方法 (Study内回归 -> Pool)
# ==============================================================================

library(tidyverse)
library(vegan)
library(lmerTest)
library(metafor)

#' 辅助函数: 判定参考组
get_comparison_direction <- function(g1, g2) {
  priority <- c(
    "CTRL" = 1, "Control" = 1, "Healthy" = 1,
    "Normal" = 1, "Normal weight" = 1, "NW" = 1,
    "AA" = 2, "Adenoma" = 2,
    "Ow" = 2, "Overweight" = 2,
    "CRC" = 3, "Cancer" = 3,
    "Ob" = 3, "Obese" = 2, "Overweight/Obese" = 2
  )
  
  p1 <- ifelse(g1 %in% names(priority), priority[g1], 99)
  p2 <- ifelse(g2 %in% names(priority), priority[g2], 99)
  
  if (p1 <= p2) {
    return(list(ref = g1, target = g2))
  } else {
    return(list(ref = g2, target = g1))
  }
}

#' 执行 Alpha 多样性分析 (包含 Unadjusted 和 Adjusted)
perform_analysis <- function(data, alpha_index, bmi_group_var, 
                             subgroup_var = "All", subgroup_value = "All", 
                             metaphlan_data_species = NULL) {
  
  # 1. 动态计算 Alpha 指数
  if (!alpha_index %in% names(data)) {
    if (is.null(metaphlan_data_species)) stop("需要提供物种丰度表来计算Alpha指数")
    common_samples <- intersect(rownames(data), colnames(metaphlan_data_species))
    if (length(common_samples) == 0) return(NULL)
    
    sub_abu <- metaphlan_data_species[, common_samples]
    if (grepl("Shannon", alpha_index)) {
      data[[alpha_index]] <- vegan::diversity(t(sub_abu), index = "shannon")
    } else if (grepl("Simpson", alpha_index)) {
      data[[alpha_index]] <- vegan::diversity(t(sub_abu), index = "simpson")
    } else if (grepl("Observed", alpha_index)) {
      data[[alpha_index]] <- colSums(sub_abu > 0)
    }
  }
  
  results_list <- list()
  
  # 移除NA值
  if(!bmi_group_var %in% names(data)) return(NULL)
  bmi_groups <- unique(na.omit(as.character(data[[bmi_group_var]])))
  
  # 内部辅助函数：处理单次比较并记录结果
  run_comparison <- function(sub_data, group_col, ref, tgt, context_val, type_label) {
    
    res_rows <- list()
    
    # 清洗数据：取交集确保 Unadj 和 Adj 使用相同的样本 (公平比较)
    clean_data <- sub_data %>% 
      filter(!is.na(.data[[group_col]]), !is.na(Age), !is.na(Sex)) %>%
      filter(.data[[group_col]] %in% c(ref, tgt)) %>%
      mutate(
        Group_Factor = factor(.data[[group_col]], levels = c(ref, tgt)), # Ref is baseline
        Sex = as.factor(Sex)
      )
    
    if(nrow(clean_data) < 10) return(NULL) 
    
    n_studies <- length(unique(clean_data$Study.name))
    
    # --------------------------------------------------------------------------
    # 1. LMM (Linear Mixed Models)
    # --------------------------------------------------------------------------
    
    # A. LMM Unadjusted
    formula_unadj <- if(n_studies >= 2) {
      as.formula(paste(alpha_index, "~ Group_Factor + (1|Study.name)"))
    } else {
      as.formula(paste(alpha_index, "~ Group_Factor"))
    }
    
    fit_unadj <- tryCatch({
      if(n_studies >= 2) lmer(formula_unadj, data = clean_data) else lm(formula_unadj, data = clean_data)
    }, error = function(e) NULL)
    
    if(!is.null(fit_unadj)) {
      coefs <- summary(fit_unadj)$coefficients
      target_idx <- grep(tgt, rownames(coefs), fixed = TRUE)
      if(length(target_idx) > 0) {
        # p-value 提取兼容性处理
        p_val <- coefs[target_idx[1], match("Pr(>|t|)", colnames(coefs))]
        if(is.na(p_val) || is.null(p_val)) {
          t_val <- coefs[target_idx[1], 3] 
          df_resid <- if(class(fit_unadj) == "lmerMod") df.residual(fit_unadj) else fit_unadj$df.residual
          p_val <- 2 * (1 - pt(abs(t_val), df_resid))
        }
        res_rows[[length(res_rows)+1]] <- data.frame(
          Method = if(n_studies >=2) "LMM" else "LM",
          Adjustment = "Unadjusted",
          Beta = coefs[target_idx[1], 1], SE = coefs[target_idx[1], 2], P_Value = p_val,
          stringsAsFactors = FALSE
        )
      }
    }
    
    # B. LMM Adjusted (Age + Sex)
    formula_adj <- if(n_studies >= 2) {
      as.formula(paste(alpha_index, "~ Group_Factor + Age + Sex + (1|Study.name)"))
    } else {
      as.formula(paste(alpha_index, "~ Group_Factor + Age + Sex"))
    }
    
    fit_adj <- tryCatch({
      if(n_studies >= 2) lmer(formula_adj, data = clean_data) else lm(formula_adj, data = clean_data)
    }, error = function(e) NULL)
    
    if(!is.null(fit_adj)) {
      coefs <- summary(fit_adj)$coefficients
      target_idx <- grep(tgt, rownames(coefs), fixed = TRUE)
      if(length(target_idx) > 0) {
        p_val <- coefs[target_idx[1], match("Pr(>|t|)", colnames(coefs))]
        if(is.na(p_val) || is.null(p_val)) {
          t_val <- coefs[target_idx[1], 3]
          df_resid <- if(class(fit_adj) == "lmerMod") df.residual(fit_adj) else fit_adj$df.residual
          p_val <- 2 * (1 - pt(abs(t_val), df_resid))
        }
        res_rows[[length(res_rows)+1]] <- data.frame(
          Method = if(n_studies >=2) "LMM" else "LM",
          Adjustment = "Adjusted_Age_Sex",
          Beta = coefs[target_idx[1], 1], SE = coefs[target_idx[1], 2], P_Value = p_val,
          stringsAsFactors = FALSE
        )
      }
    }
    
    # --------------------------------------------------------------------------
    # 2. Meta-Analysis
    # --------------------------------------------------------------------------
    if (n_studies >= 2) {
      
      # A. Meta Unadjusted (SMD based on Summary Stats)
      # -----------------------------------------------
      stats_summary <- clean_data %>%
        group_by(Study.name) %>%
        summarise(
          n1 = sum(Group_Factor == tgt, na.rm=T),
          m1 = mean(.data[[alpha_index]][Group_Factor == tgt], na.rm=T),
          sd1 = sd(.data[[alpha_index]][Group_Factor == tgt], na.rm=T),
          n2 = sum(Group_Factor == ref, na.rm=T),
          m2 = mean(.data[[alpha_index]][Group_Factor == ref], na.rm=T),
          sd2 = sd(.data[[alpha_index]][Group_Factor == ref], na.rm=T),
          .groups = "drop"
        ) %>% filter(n1 >= 3, n2 >= 3)
      
      if(nrow(stats_summary) >= 2) {
        es_unadj <- escalc(measure = "SMD", 
                           m1i=m1, sd1i=sd1, n1i=n1, 
                           m2i=m2, sd2i=sd2, n2i=n2, 
                           data = stats_summary)
        
        rma_unadj <- tryCatch({rma(yi, vi, data = es_unadj, method = "REML", test = "knha")}, error=function(e) NULL)
        
        if(!is.null(rma_unadj)) {
          res_rows[[length(res_rows)+1]] <- data.frame(
            Method = "Meta_Analysis",
            Adjustment = "Unadjusted",
            Beta = rma_unadj$beta[1], SE = rma_unadj$se, P_Value = rma_unadj$pval,
            stringsAsFactors = FALSE
          )
        }
      }
      
      # B. Meta Adjusted (Two-Stage: LM per study -> Pool)
      # --------------------------------------------------
      study_list <- unique(clean_data$Study.name)
      meta_adj_df <- data.frame()
      
      for(s in study_list) {
        s_data <- clean_data %>% filter(Study.name == s)
        if(nrow(s_data) < 6 || length(unique(s_data$Group_Factor)) < 2) next
        
        fit_s <- tryCatch({
          lm(as.formula(paste(alpha_index, "~ Group_Factor + Age + Sex")), data = s_data)
        }, error = function(e) NULL)
        
        if(!is.null(fit_s)) {
          coefs_s <- summary(fit_s)$coefficients
          target_idx_s <- grep(tgt, rownames(coefs_s), fixed=TRUE)
          if(length(target_idx_s) > 0) {
            meta_adj_df <- rbind(meta_adj_df, data.frame(
              yi = coefs_s[target_idx_s[1], 1], # Beta
              vi = (coefs_s[target_idx_s[1], 2])^2 # SE^2 = Variance
            ))
          }
        }
      }
      
      if(nrow(meta_adj_df) >= 2) {
        rma_adj <- tryCatch({rma(yi, vi, data = meta_adj_df, method = "REML", test = "knha")}, error=function(e) NULL)
        
        if(!is.null(rma_adj)) {
          res_rows[[length(res_rows)+1]] <- data.frame(
            Method = "Meta_Analysis",
            Adjustment = "Adjusted_Age_Sex",
            Beta = rma_adj$beta[1], SE = rma_adj$se, P_Value = rma_adj$pval,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    
    # 整理并返回数据框
    if(length(res_rows) > 0) {
      final_df <- bind_rows(res_rows)
      final_df$Type <- type_label
      final_df$Context <- context_val
      final_df$Comparison <- paste(ref, "vs", tgt)
      final_df$Ref_Group <- ref
      final_df$Target_Group <- tgt
      final_df$BMI_Group_Column <- bmi_group_var
      final_df$Alpha_Index <- alpha_index
      final_df$Subgroup_Var <- subgroup_var
      final_df$Subgroup_Val <- subgroup_value
      return(final_df)
    } else {
      return(NULL)
    }
  }
  
  # ==========================================================================
  # 执行比较循环
  # ==========================================================================
  
  # A. BMI Comparison
  if (length(bmi_groups) >= 2) {
    disease_groups <- intersect(c("CTRL", "AA", "CRC"), unique(data$Disease))
    if(length(disease_groups) == 0) disease_groups <- unique(data$Disease)
    
    bmi_combs <- combn(bmi_groups, 2, simplify = FALSE)
    for (dz in disease_groups) {
      dz_data <- data %>% filter(Disease == dz)
      for (comb in bmi_combs) {
        dir <- get_comparison_direction(comb[1], comb[2])
        res <- run_comparison(dz_data, bmi_group_var, dir$ref, dir$target, 
                              context_val = dz, type_label = "BMI_Comparison")
        if(!is.null(res)) results_list[[length(results_list)+1]] <- res
      }
    }
  }
  
  # B. Disease Comparison
  valid_diseases <- unique(as.character(data$Disease))
  if(length(valid_diseases) >= 2) {
    dz_combs <- combn(valid_diseases, 2, simplify = FALSE)
    for (bmi in bmi_groups) {
      bmi_data <- data %>% filter(.data[[bmi_group_var]] == bmi)
      for (comb in dz_combs) {
        dir <- get_comparison_direction(comb[1], comb[2])
        res <- run_comparison(bmi_data, "Disease", dir$ref, dir$target, 
                              context_val = bmi, type_label = "Disease_Comparison")
        if(!is.null(res)) results_list[[length(results_list)+1]] <- res
      }
    }
  }
  
  return(bind_rows(results_list)) # 返回单一数据框，解决 bind_rows 报错
}