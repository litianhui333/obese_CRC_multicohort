#### 微生物组分析参数配置系统 ####

#' 创建微生物组分析配置
#'
#' @param transform_method 数据转换方法，"AST"或"LOG"
#' @param analysis_method 分析方法，"LM"或"HS"
#' @param covariates 需要校正的协变量，默认c("Age", "Sex")
#' @param min_studies_ratio 最小研究数比例，默认1/3
#' @param prevalence_threshold 出现率阈值，默认0.1(10%)
#' @param abundance_threshold 丰度阈值，默认0.01(0.01%)
#' @param is_inhouse 是否为inhouse数据，默认FALSE
#' @return 包含所有分析参数的列表

create_analysis_config <- function(
  transform_method = "AST",
  analysis_method = "LM",
  covariates = c("Age", "Sex"),
  min_studies_ratio = 1 / 3,
  BMI_Group_var = "BMI_Group1_China_as_OtherAsia",
  study_prevalence_cut = 0.1,
  abundance_threshold = 1e-4, # 建议设为 1e-4 (0.01%)，0.01 (1%) 有点太高了
  pval_type = "both",
  sig_criterion = "or",
  pval_threshold = 0.05,
  qval_threshold = 0.2, # 根据刚才讨论，建议默认放宽到 0.2
  is_inhouse = FALSE,
  force_keep = NULL # <--- 新增：白名单参数
) {
  # === 参数验证 ===
  if (!transform_method %in% c("AST", "LOG")) {
    stop("错误：transform_method必须是'AST'或'LOG'")
  }

  if (!analysis_method %in% c("LM", "HS")) {
    stop("错误：analysis_method必须是'LM'或'HS'")
  }

  if (min_studies_ratio <= 0 || min_studies_ratio > 1) {
    stop("错误：min_studies_ratio必须在0到1之间")
  }

  if (study_prevalence_cut < 0 || study_prevalence_cut > 1) {
    stop("错误：study_prevalence_cut必须在0到1之间")
  }

  if (abundance_threshold <= 0) {
    stop("错误：abundance_threshold必须大于0")
  }

  if (!pval_type %in% c("pval", "qval.fdr", "both")) {
    stop("错误：pval_type必须是 'pval','qval.fdr' 或 'both'")
  }

  if (!sig_criterion %in% c("or", "and")) {
    stop("错误：sig_criterion 'or' 或 'and'")
  }

  if (pval_threshold > 1) { # 放宽一点检查，允许用户设置宽松阈值
    stop("错误：pval_threshold 必须 <= 1")
  }

  if (qval_threshold > 1) {
    stop("错误：qval_threshold 必须 <= 1 ")
  }

  # 验证 force_keep
  if (!is.null(force_keep) && !is.character(force_keep)) {
    stop("错误：force_keep 必须是字符向量 (character vector)")
  }

  # === 创建配置列表 ===
  config <- list(
    transform_method = transform_method,
    analysis_method = analysis_method,
    covariates = covariates,
    BMI_Group_var = BMI_Group_var,
    min_studies_ratio = min_studies_ratio,
    study_prevalence_cut = study_prevalence_cut,
    abundance_threshold = abundance_threshold,
    pval_type = pval_type,
    sig_criterion = sig_criterion,
    pval_threshold = pval_threshold,
    qval_threshold = qval_threshold,
    is_inhouse = is_inhouse,
    force_keep = force_keep # <--- 存入配置
  )

  # === 添加配置摘要信息 ===
  # 处理白名单的显示文本
  force_keep_info <- if (is.null(force_keep)) {
    "无"
  } else {
    paste0(length(force_keep), "个特征 (", paste(head(force_keep, 2), collapse = ","), "...)")
  }

  config$summary <- sprintf(
    "分析配置摘要:\n  - 转换/分析: %s / %s\n  - 协变量: %s (BMI分组: %s)\n  - 过滤标准: 研究比例=%.2f, 出现率=%.1f%%, 丰度阈值=%.4f%%\n  - 统计阈值: P<%.2f, Q<%.2f (%s / %s)\n  - 白名单(强制保留): %s",
    transform_method,
    analysis_method,
    paste(covariates, collapse = ","),
    BMI_Group_var,
    min_studies_ratio,
    study_prevalence_cut * 100,
    abundance_threshold * 100,
    pval_threshold,
    qval_threshold,
    pval_type,
    sig_criterion,
    force_keep_info # <--- 显示白名单信息
  )

  cat("=== 微生物组分析配置已更新 ===\n")
  cat(config$summary, "\n")

  return(config)
}

#' 获取研究信息（不进行数量验证）
#'
#' 只统计研究信息，不验证研究数量是否足够
#'
#' @param metadata 样本元数据，必须包含Study.name列
#' @param config 分析配置列表
#' @return 包含研究统计信息的列表

get_study_info <- function(metadata, config) {
  cat("=== 获取研究信息 ===\n")

  # 检查metadata结构
  if (!"Study.name" %in% colnames(metadata)) {
    stop("错误：metadata必须包含Study.name列")
  }

  # 获取研究信息
  studies <- unique(metadata$Study.name)
  total_studies <- length(studies)

  cat("检测到的研究数量:", total_studies, "\n")
  cat("研究名称:", paste(studies, collapse = ", "), "\n")

  # 计算用于条件3的最小研究数（但不验证是否足够）
  if (config$is_inhouse) {
    min_studies_required <- 1
    cat("inhouse数据模式: 最小研究数要求 = 1\n")
  } else {
    min_studies_required <- max(1, floor(total_studies * config$min_studies_ratio))
    cat(sprintf(
      "多研究模式: 最小研究数要求 = floor(%d × %.2f) = %d\n",
      total_studies, config$min_studies_ratio, min_studies_required
    ))
  }

  # 计算每个研究的样本数（仅用于信息显示，不用于验证）
  study_sample_counts <- table(metadata$Study.name)
  cat("\n各研究样本数量:\n")
  print(study_sample_counts)

  # 返回研究信息（不包含验证状态）
  study_info <- list(
    total_studies = total_studies,
    study_names = studies,
    min_studies_required = min_studies_required, # 仅用于条件3计算
    study_sample_counts = study_sample_counts
  )

  cat("✅ 研究信息获取完成\n")

  return(study_info)
}

#' 创建MMUPHin分析控制参数
#'
#' @param config 分析配置列表
#' @return 用于MMUPHin::lm_meta的control参数列表

create_mmupihin_control <- function(config) {
  control_params <- list(
    transform = config$transform_method,
    analysis_method = "LM", # MMUPHin固定参数
    rma_method = config$analysis_method,
    verbose = FALSE
  )

  cat("=== MMUPHin控制参数已创建 ===\n")
  cat(sprintf(
    "转换方法: %s, 元分析方法: %s\n",
    control_params$transform, control_params$rma_method
  ))

  return(control_params)
}

#### 智能数据过滤模块 ####

#' 智能数据过滤函数
#'
#' 根据配置参数过滤低丰度和低出现率的微生物特征
#'
#' @param metaphlan_data 微生物丰度矩阵（行是物种，列是样本）
#' @param metadata 样本元数据（必须包含Sample.Source.ID和Study.name）
#' @param config 分析配置列表
#' @return 包含过滤后数据和过滤统计的列表

smart_data_filtering <- function(metaphlan_data, metadata, config) {
  cat("=== 开始智能数据过滤 ===\n")

  # 自动检测方向并对齐
  common_cols <- intersect(colnames(metaphlan_data), metadata$Sample.Source.ID.processed)
  common_rows <- intersect(rownames(metaphlan_data), metadata$Sample.Source.ID.processed)

  if (length(common_cols) >= length(common_rows) && length(common_cols) > 0) {
    metaphlan_data <- metaphlan_data[, common_cols, drop = FALSE]
    metadata <- metadata[metadata$Sample.Source.ID.processed %in% common_cols, ]
  } else if (length(common_rows) > 0) {
    metaphlan_data <- t(metaphlan_data[common_rows, , drop = FALSE])
    metadata <- metadata[metadata$Sample.Source.ID.processed %in% common_rows, ]
  } else {
    stop("错误: 无法对齐样本ID，请检查数据。")
  }

  total_features <- nrow(metaphlan_data)
  studies <- unique(metadata$Study.name)

  study_presence_count <- rep(0, total_features)
  for (s in studies) {
    s_samples <- metadata$Sample.Source.ID.processed[metadata$Study.name == s]
    if (length(s_samples) > 0) {
      s_data <- metaphlan_data[, s_samples, drop = FALSE]
      prev <- rowSums(s_data > config$abundance_threshold) / length(s_samples)
      study_presence_count <- study_presence_count + (prev > config$study_prevalence_cut)
    }
  }

  min_studies <- round(config$min_studies_ratio * length(studies))
  keep_breadth <- study_presence_count >= min_studies

  keep_force <- rep(FALSE, total_features)
  if (!is.null(config$force_keep)) {
    pattern <- paste(config$force_keep, collapse = "|")
    keep_force <- grepl(pattern, rownames(metaphlan_data), ignore.case = TRUE)
  }

  final_keep <- keep_breadth | keep_force
  filtered_data <- metaphlan_data[final_keep, , drop = FALSE]

  cat(sprintf("   特征过滤: 原始 %d -> 保留 %d\n", total_features, nrow(filtered_data)))
  return(list(filtered_data = filtered_data, keep_indices = final_keep))
}


#' 生成比较组（修正版）
#'
#' 根据分组自动生成所有可能的两两比较
#' 如果没有指定参考组，则生成所有两两比较
#'
#' @param groups 所有分组名称
#' @param reference_group 参考组（可选，如不指定则生成所有两两比较）
#' @return 比较组列表

generate_comparisons_fixed <- function(groups, reference_group = NULL, fixed_order = NULL) {
  # 如果指定了固定顺序，使用固定顺序
  if (!is.null(fixed_order)) {
    # 确保所有组都在固定顺序中
    if (!all(groups %in% names(fixed_order))) {
      stop("固定顺序必须包含所有分组")
    }
    groups <- groups[order(match(groups, names(fixed_order)))]
  }

  # 如果指定了参考组，确保参考组在前面
  if (!is.null(reference_group)) {
    if (!reference_group %in% groups) {
      stop("参考组不在分组中")
    }
    # 将参考组移到最前面
    groups <- c(reference_group, setdiff(groups, reference_group))
  }

  # 生成两两比较，确保顺序
  comparisons <- list()
  for (i in 1:(length(groups) - 1)) {
    for (j in (i + 1):length(groups)) {
      comp_name <- paste0(groups[i], "_vs_", groups[j])
      comparisons[[comp_name]] <- c(groups[i], groups[j])
    }
  }

  return(comparisons)
}

#' 修改后的通用分组差异分析函数
#'
#' 支持固定比较方向

perform_group_analysis_fixed <- function(filtered_data, metadata, config,
                                         group_var, reference_group = NULL,
                                         comparisons = NULL,
                                         fixed_order = NULL) {
  cat("=== 通用分组差异分析（固定方向） ===\n")
  cat("分组变量:", group_var, "\n")

  # 1. 数据准备和验证
  cat("1. 准备分析数据...\n")

  if (!group_var %in% colnames(metadata)) {
    stop("错误：分组变量'", group_var, "'不存在于元数据中")
  }

  groups <- unique(metadata[[group_var]])
  groups <- groups[!is.na(groups)]
  cat("   检测到的分组: ", paste(groups, collapse = ", "), "\n")
  cat("   分组样本数:\n")
  print(table(metadata[[group_var]]))

  # 2. 自动生成比较组（使用固定顺序）
  cat("2. 生成比较组（固定顺序）...\n")

  if (is.null(comparisons)) {
    comparisons <- generate_comparisons_fixed(groups, reference_group, fixed_order)
  } else {
    # 如果用户提供了比较组，验证并标准化
    comparisons <- standardize_comparisons(comparisons, groups, reference_group, fixed_order)
  }

  cat("   比较组数量: ", length(comparisons), "\n")
  cat("   比较组顺序:\n")
  for (comp_name in names(comparisons)) {
    cat("     ", comp_name, ": ", comparisons[[comp_name]][1], " vs ", comparisons[[comp_name]][2], "\n")
  }

  # 获取汇总统计函数
  summary_func <- if (is.null(config$summary_stat) || config$summary_stat == "mean") {
    function(x) rowMeans(x, na.rm = TRUE)
  } else {
    function(x) apply(x, 1, median, na.rm = TRUE)
  }

  # 3. 执行差异分析
  cat("3. 执行差异分析...\n")

  all_results <- list()

  for (comp_name in names(comparisons)) {
    cat("   分析比较: ", comp_name, "\n")

    comp_groups <- comparisons[[comp_name]]
    reference <- comp_groups[1]
    comparison <- comp_groups[2]
    cat("     参考组:", reference, " vs 比较组:", comparison, "\n")

    # 筛选当前比较的样本
    comp_samples <- metadata$Sample.Source.ID.processed[metadata[[group_var]] %in% comp_groups]
    comp_metadata <- metadata[metadata$Sample.Source.ID.processed %in% comp_samples, ]
    comp_data <- filtered_data[, comp_samples, drop = FALSE]

    # 执行单次比较分析
    comp_result <- analyze_single_comparison_fixed(
      comp_data, comp_metadata, config, group_var, comp_groups, comp_name
    )

    if (is.null(comp_result)) {
      cat("     分析失败，跳过\n")
      next
    }

    # 计算每组均值/中位数
    cat("     计算各组汇总统计量...\n")
    for (group in comp_groups) {
      group_samples <- comp_metadata$Sample.Source.ID.processed[comp_metadata[[group_var]] == group]

      if (length(group_samples) > 0) {
        group_summary <- summary_func(comp_data[, group_samples, drop = FALSE])
        col_name <- paste0("mean_", group)
        comp_result$meta_fits[[col_name]] <- group_summary[rownames(comp_result$meta_fits)]
      }
    }

    # 添加解释性列（始终基于固定顺序）
    comp_result$meta_fits$direction <- ifelse(
      comp_result$meta_fits$coef > 0,
      paste0("up_in_", comparison), # coef>0表示在比较组中上调
      paste0("down_in_", comparison) # coef<0表示在比较组中下调
    )

    # 计算组间差值（比较组 - 参考组）
    ref_mean_col <- paste0("mean_", reference)
    comp_mean_col <- paste0("mean_", comparison)

    if (ref_mean_col %in% colnames(comp_result$meta_fits) &&
      comp_mean_col %in% colnames(comp_result$meta_fits)) {
      comp_result$meta_fits$mean_diff <- comp_result$meta_fits[[comp_mean_col]] - comp_result$meta_fits[[ref_mean_col]]
    }

    # 添加比较信息
    comp_result$comparison_info <- list(
      comparison_name = comp_name,
      reference_group = reference,
      comparison_group = comparison,
      interpretation = paste0("coef>0: 在", comparison, "中更高; coef<0: 在", reference, "中更高"),
      fixed_direction = TRUE
    )

    all_results[[comp_name]] <- comp_result
  }

  # 4. 结果整合
  cat("4. 整合分析结果...\n")
  significant_features <- extract_significant_features_fixed(all_results, config)

  # 5. 生成分析报告
  analysis_report <- generate_analysis_report(all_results, significant_features, group_var)

  # 6. 返回完整结果
  final_result <- list(
    group_variable = group_var,
    groups = groups,
    comparisons = comparisons,
    all_results = all_results,
    significant_features = significant_features,
    report = analysis_report,
    config = config
  )

  cat("✅ 分组分析完成: ", group_var, "\n")
  cat("   发现显著特征数: ", nrow(significant_features), "\n")

  return(final_result)
}


#' 标准化比较组
standardize_comparisons <- function(comparisons, groups, reference_group = NULL, fixed_order = NULL) {
  standardized <- list()

  for (comp_name in names(comparisons)) {
    comp_groups <- comparisons[[comp_name]]

    if (length(comp_groups) != 2) {
      stop("比较组必须包含两个分组")
    }

    # 确定顺序
    if (!is.null(fixed_order)) {
      # 使用固定顺序
      order1 <- match(comp_groups[1], names(fixed_order))
      order2 <- match(comp_groups[2], names(fixed_order))

      if (order1 > order2) {
        # 交换顺序
        comp_groups <- rev(comp_groups)
        comp_name <- paste0(comp_groups[1], "_vs_", comp_groups[2])
      }
    } else if (!is.null(reference_group)) {
      # 确保参考组在前
      if (comp_groups[2] == reference_group) {
        comp_groups <- rev(comp_groups)
        comp_name <- paste0(comp_groups[1], "_vs_", comp_groups[2])
      }
    }

    standardized[[comp_name]] <- comp_groups
  }

  return(standardized)
}

#' 修改后的单次比较分析
analyze_single_comparison_fixed <- function(comp_data, comp_metadata, config,
                                            group_var, comp_groups, comp_name) {
  reference <- comp_groups[1]
  comparison <- comp_groups[2]

  # 创建二分类变量（参考组=0，比较组=1）
  comp_metadata$group_binary <- ifelse(comp_metadata[[group_var]] == reference, 0, 1)

  # 样本数检查
  n_reference <- sum(comp_metadata$group_binary == 0)
  n_comparison <- sum(comp_metadata$group_binary == 1)

  cat("     样本数: ", reference, "=", n_reference, ", ", comparison, "=", n_comparison, "\n")

  if (n_reference < 3 | n_comparison < 3) {
    warning("比较组", comp_name, "样本数不足，结果可能不可靠")
    return(NULL)
  }

  library(MMUPHin)

  tryCatch(
    {
      studies <- unique(comp_metadata$Study.name)

      if (length(studies) > 1) {
        result <- lm_meta(
          feature_abd = comp_data,
          batch = "Study.name",
          exposure = "group_binary",
          covariates = config$covariates,
          data = comp_metadata,
          control = list(
            transform = config$transform_method,
            analysis_method = "LM",
            rma_method = config$analysis_method,
            verbose = FALSE
          )
        )
      } else {
        result <- perform_single_study_analysis(comp_data, comp_metadata, config, "group_binary")
      }

      result$comparison_info <- list(
        comparison_name = comp_name,
        groups = comp_groups,
        sample_sizes = c(n_reference, n_comparison),
        reference_group = reference,
        comparison_group = comparison
      )

      return(result)
    },
    error = function(e) {
      warning("比较组", comp_name, "分析失败: ", e$message)
      return(NULL)
    }
  )
}

#' 提取显著特征
#'
#' 从所有比较结果中提取显著差异特征
#'
#' @param all_results 所有比较结果
#' @param config 分析配置
#' @return 显著特征数据框

#' 提取显著特征（固定方向版本）
#' 提取显著特征（固定方向版本）- 可配置阈值
extract_significant_features_fixed <- function(all_results, config) {
  library(dplyr)

  # 从config中读取阈值参数，设置默认值
  pval_threshold <- if (!is.null(config$pval_threshold)) config$pval_threshold else 0.05
  qval_threshold <- if (!is.null(config$qval_threshold)) config$qval_threshold else 0.05

  # 读取要使用的p值类型（默认使用qval.fdr）
  pval_type <- if (!is.null(config$pval_type)) config$pval_type else "qval.fdr"

  # 读取显著性标准（"or" 或 "and"，默认"or"）
  sig_criterion <- if (!is.null(config$sig_criterion)) config$sig_criterion else "or"

  significant_list <- list()

  for (comp_name in names(all_results)) {
    result <- all_results[[comp_name]]

    if (!is.null(result) && !is.null(result$meta_fits)) {
      # 提取当前比较的显著特征
      comp_features <- result$meta_fits
      comp_features$comparison <- comp_name

      # 从比较信息中获取参考组和比较组
      ref_group <- result$comparison_info$reference_group
      comp_group <- result$comparison_info$comparison_group

      # 应用可配置的显著性阈值
      if (pval_type == "pval" || pval_type == "pval.raw") {
        # 使用原始p值
        pval_col <- if ("pval" %in% colnames(comp_features)) "pval" else if ("pval.raw" %in% colnames(comp_features)) "pval.raw" else "pval"
        if (pval_col %in% colnames(comp_features)) {
          sig_features <- comp_features[comp_features[[pval_col]] < pval_threshold, ]
        } else {
          warning("在比较 ", comp_name, " 中未找到原始p值列")
          next
        }
      } else if (pval_type == "qval.fdr") {
        # 使用FDR校正的q值
        qval_col <- if ("qval.fdr" %in% colnames(comp_features)) "qval.fdr" else if ("pval.adj" %in% colnames(comp_features)) "pval.adj" else if ("qval" %in% colnames(comp_features)) "qval" else "qval.fdr"

        if (qval_col %in% colnames(comp_features)) {
          sig_features <- comp_features[comp_features[[qval_col]] < qval_threshold, ]
        } else {
          warning("在比较 ", comp_name, " 中未找到FDR校正的q值列")
          next
        }
      } else if (pval_type == "both") {
        # 同时使用原始p值和校正q值（"或"的关系）
        pval_col <- if ("pval" %in% colnames(comp_features)) "pval" else if ("pval.raw" %in% colnames(comp_features)) "pval.raw" else NULL
        qval_col <- if ("qval.fdr" %in% colnames(comp_features)) "qval.fdr" else if ("pval.adj" %in% colnames(comp_features)) "pval.adj" else if ("qval" %in% colnames(comp_features)) "qval" else NULL

        if (is.null(pval_col) || is.null(qval_col)) {
          warning("在比较 ", comp_name, " 中未找到所需的p值和q值列")
          next
        }

        if (sig_criterion == "or") {
          # 或关系：满足任意一个阈值
          sig_features <- comp_features[
            comp_features[[pval_col]] < pval_threshold |
              comp_features[[qval_col]] < qval_threshold,
          ]
        } else if (sig_criterion == "and") {
          # 与关系：同时满足两个阈值
          sig_features <- comp_features[
            comp_features[[pval_col]] < pval_threshold &
              comp_features[[qval_col]] < qval_threshold,
          ]
        } else {
          warning("不支持的sig_criterion: ", sig_criterion, "，使用默认的'or'")
          sig_features <- comp_features[
            comp_features[[pval_col]] < pval_threshold |
              comp_features[[qval_col]] < qval_threshold,
          ]
        }
      } else {
        warning("不支持的pval_type: ", pval_type, "，使用默认的'qval.fdr'")
        # 回退到默认
        qval_col <- if ("qval.fdr" %in% colnames(comp_features)) "qval.fdr" else "qval"
        if (qval_col %in% colnames(comp_features)) {
          sig_features <- comp_features[comp_features[[qval_col]] < qval_threshold, ]
        } else {
          warning("在比较 ", comp_name, " 中未找到q值列")
          next
        }
      }

      # 添加比较信息 (skip if empty)
      if (nrow(sig_features) > 0) {
        sig_features$reference_group <- ref_group
        sig_features$comparison_group <- comp_group
      }

      significant_list[[comp_name]] <- sig_features
    }
  }

  # 合并所有显著特征
  if (length(significant_list) > 0) {
    all_significant <- dplyr::bind_rows(significant_list, .id = "source_comparison")

    # 移除可能的重复行
    all_significant <- all_significant[!duplicated(all_significant), ]

    return(all_significant)
  } else {
    return(data.frame()) # 空数据框
  }
}

#' 执行BMI连续变量与菌群的MMUPHin分析
#'
#' @param filtered_data 过滤后的菌群丰度数据（行是特征，列是样本）
#' @param metadata 元数据，必须包含BMI列
#' @param config 配置列表
#' @param bmi_var BMI变量名，默认为"BMI"
#' @param covariates 协变量，从config获取
#' @param min_samples_per_study 每个研究的最小样本数
#'
#' @return 分析结果列表
perform_bmi_analysis <- function(filtered_data, metadata, config,
                                 bmi_var = "BMI",
                                 by_disease_group = FALSE,
                                 disease_var = "Group",
                                 include_overall = TRUE, # 新增：是否包含全样本
                                 min_samples_per_study = 1) {
  # --- 增加防御性修复 ---
  # 1. 确保 bmi_var 只是一个字符串
  if (length(bmi_var) > 1) {
    warning("检测到 bmi_var 包含多个值，仅使用第一个：", bmi_var[1])
    bmi_var <- bmi_var[1]
  }

  # 2. 确保 metadata 是标准的 data.frame（解决 tibble 兼容性）
  metadata <- as.data.frame(metadata)


  cat("=== BMI连续变量与菌群关联分析 ===\n")

  # 1. 数据验证
  cat("1. 验证数据...\n")

  if (!bmi_var %in% colnames(metadata)) {
    stop("错误：BMI变量'", bmi_var, "'不存在于元数据中")
  }

  if (!is.numeric(metadata[[bmi_var]])) {
    warning("BMI变量不是数值型，尝试转换为数值型...")
    metadata[[bmi_var]] <- as.numeric(as.character(metadata[[bmi_var]]))
  }

  covariates <- config$covariates

  # 移除缺失的BMI值
  valid_samples <- metadata$Sample.Source.ID.processed[!is.na(metadata[[bmi_var]])]
  cat("   有效样本数（有BMI数据）: ", length(valid_samples), "\n")

  if (length(valid_samples) < 10) {
    stop("错误：BMI有效样本数不足（<10）")
  }

  # 筛选数据
  bmi_data <- filtered_data[, colnames(filtered_data) %in% valid_samples, drop = FALSE]
  bmi_metadata <- metadata[metadata$Sample.Source.ID.processed %in% valid_samples, ]
  bmi_data <- bmi_data[, bmi_metadata$Sample.Source.ID.processed] # 确保顺序

  # 2. 判断是否按疾病分组分析
  if (by_disease_group) {
    cat("2. 将进行分层分析...\n")

    if (!disease_var %in% colnames(bmi_metadata)) {
      stop("错误：疾病分组变量'", disease_var, "'不存在于元数据中")
    }

    # 获取所有疾病组
    disease_groups <- unique(bmi_metadata[[disease_var]])
    disease_groups <- disease_groups[!is.na(disease_groups)]

    cat("   找到疾病组: ", paste(disease_groups, collapse = ", "), "\n")

    # 存储所有结果
    all_results <- list()

    # 首先分析全样本（如果要求包含）
    if (include_overall) {
      cat("\n--- 分析全样本 ---\n")
      cat("   样本数: ", nrow(bmi_metadata), "\n")

      overall_result <- perform_bmi_analysis_single_dataset(
        data = bmi_data,
        metadata = bmi_metadata,
        bmi_var = bmi_var,
        covariates = covariates,
        config = config,
        min_samples_per_study = min_samples_per_study
      )

      all_results[["Overall"]] <- overall_result
    }

    # 然后分析每个疾病亚组
    for (group in disease_groups) {
      cat("\n--- 分析疾病组: ", group, " ---\n")

      # 筛选该疾病组的样本
      group_samples <- bmi_metadata$Sample.Source.ID.processed[bmi_metadata[[disease_var]] == group]
      group_data <- bmi_data[, colnames(bmi_data) %in% group_samples, drop = FALSE]
      group_metadata <- bmi_metadata[bmi_metadata$Sample.Source.ID.processed %in% group_samples, ]

      cat("   样本数: ", nrow(group_metadata), "\n")

      if (nrow(group_metadata) >= 10) {
        # 直接调用原有的分析流程
        result <- perform_bmi_analysis_single_dataset(
          data = group_data,
          metadata = group_metadata,
          bmi_var = bmi_var,
          covariates = covariates,
          config = config,
          min_samples_per_study = min_samples_per_study
        )

        all_results[[group]] <- result
      } else {
        cat("   样本数不足，跳过该组分析\n")
        all_results[[group]] <- list(
          group = group,
          n_samples = nrow(group_metadata),
          message = "样本数不足，跳过分析"
        )
      }
    }

    # 3. 比较不同组的结果
    cat("\n3. 比较不同组的结果...\n")
    comparison <- compare_group_results(all_results)

    # 4. 返回所有结果
    final_result <- list(
      analysis_type = "continuous_bmi_stratified",
      bmi_variable = bmi_var,
      stratify_by = disease_var,
      results_by_group = all_results,
      comparison_results = comparison,
      config = config
    )
  } else {
    # 不在疾病亚组中分析，只做全样本
    cat("2. 在全样本中进行分析...\n")

    # 直接调用原有的分析流程
    result <- perform_bmi_analysis_single_dataset(
      data = bmi_data,
      metadata = bmi_metadata,
      bmi_var = bmi_var,
      covariates = covariates,
      config = config,
      min_samples_per_study = min_samples_per_study
    )

    final_result <- result
  }

  cat("\n✅ BMI分析完成!\n")
  return(final_result)
}

# 比较不同组的结果
compare_group_results <- function(all_results) {
  cat("   比较不同组的结果...\n")

  # 收集所有组的显著特征
  all_sig_features <- list()
  all_effects <- data.frame()

  for (group_name in names(all_results)) {
    group_result <- all_results[[group_name]]

    if (!is.null(group_result$significant_features) && nrow(group_result$significant_features) > 0) {
      sig_features <- group_result$significant_features
      sig_features$group <- group_name
      sig_features$feature <- rownames(sig_features)

      all_sig_features[[group_name]] <- sig_features
    }
  }

  # 找出各组共同的显著特征
  common_features <- NULL
  if (length(all_sig_features) > 1) {
    feature_sets <- lapply(all_sig_features, function(x) unique(x$feature))
    common_features <- Reduce(intersect, feature_sets)
    cat("   各组共同的显著特征数量: ", length(common_features), "\n")
  }

  # 找出各组特有的显著特征
  unique_features <- list()
  for (group_name in names(all_sig_features)) {
    other_groups <- setdiff(names(all_sig_features), group_name)
    other_features <- unique(unlist(lapply(all_sig_features[other_groups], function(x) x$feature)))
    group_features <- all_sig_features[[group_name]]$feature
    unique_features[[group_name]] <- setdiff(group_features, other_features)
    cat("   ", group_name, "组特有的显著特征数量: ", length(unique_features[[group_name]]), "\n")
  }

  return(list(
    common_features = common_features,
    unique_features = unique_features
  ))
}


# 辅助函数：对单个数据集执行完整的BMI分析
perform_bmi_analysis_single_dataset <- function(data, metadata, bmi_var, covariates,
                                                config, min_samples_per_study) {
  # 检查研究数量和样本数
  studies <- unique(metadata$Study.name)
  cat("   研究数量: ", length(studies), "\n")

  study_sample_counts <- table(metadata$Study.name)
  cat("   各研究样本数:\n")
  print(study_sample_counts)

  # 如果提供了协变量，则使用
  if (!is.null(covariates)) {
    cat("   使用协变量: ", paste(covariates, collapse = ", "), "\n")
  }

  result <- list(meta_fits = NULL)

  # 检查是否有足够的研究进行meta分析
  if (length(studies) > 1) {
    cat("   检测到多个研究，执行meta分析...\n")

    # 确保有足够的研究
    valid_studies <- names(study_sample_counts[study_sample_counts >= min_samples_per_study])

    if (length(valid_studies) < 2) {
      cat("   警告：只有", length(valid_studies), "个研究有足够样本，改用单研究分析\n")
      # 如果没有定义单研究函数，这里需要定义或者直接处理
      result <- analyze_single_study_bmi(data, metadata, bmi_var, config, covariates)
    } else {
      # 只保留有足够样本的研究
      metadata_valid <- metadata[metadata$Study.name %in% valid_studies, ]
      data_valid <- data[, metadata_valid$Sample.Source.ID.processed]

      cat("   用于meta分析的研究: ", paste(valid_studies, collapse = ", "), "\n")

      # --- [核心修复] 直接调用 MMUPHin::lm_meta，而不是 perform_bmi_analysis ---
      library(MMUPHin)

      # 准备 Control 参数
      control_params <- list(
        transform = config$transform_method,
        analysis_method = "LM",
        rma_method = config$analysis_method,
        verbose = FALSE
      )

      tryCatch(
        {
          fit_res <- lm_meta(
            feature_abd = data_valid,
            batch = "Study.name",
            exposure = bmi_var,
            covariates = covariates,
            data = metadata_valid,
            control = control_params
          )
          result <- list(
            meta_fits = fit_res$meta_fits,
            maaslin_fits = fit_res$maaslin_fits
          )
        },
        error = function(e) {
          warning("   MMUPHin meta分析失败: ", e$message)
          result <- list(meta_fits = NULL)
        }
      )
      # -------------------------------------------------------------
    }
  } else {
    cat("   只有一个研究，执行单研究分析...\n")
    result <- analyze_single_study_bmi(data, metadata, bmi_var, config, covariates)
  }

  # 结果整理
  cat("   整理分析结果...\n")

  if (!is.null(result$meta_fits)) {
    # 添加特征信息
    result$meta_fits$feature_id <- rownames(result$meta_fits)
    result$meta_fits$feature_name <- rownames(result$meta_fits)

    # 计算FDR (如果lm_meta没算或需要重算)
    if (!"qval.fdr" %in% colnames(result$meta_fits)) {
      result$meta_fits$qval.fdr <- p.adjust(result$meta_fits$pval, method = "fdr")
    }

    # 排序结果
    result$meta_fits <- result$meta_fits[order(result$meta_fits$pval), ]

    # 识别显著特征
    sig_cutoff <- ifelse(!is.null(config$qval_threshold), config$qval_threshold, 0.2)
    sig_features <- result$meta_fits[result$meta_fits$qval.fdr < sig_cutoff, ]

    cat("   总共分析特征数: ", nrow(result$meta_fits), "\n")
    cat("   显著特征数 (FDR < ", sig_cutoff, "): ", nrow(sig_features), "\n")

    # 添加方向解释
    result$meta_fits$direction <- ifelse(result$meta_fits$coef > 0,
      "positive_association",
      "negative_association"
    )
    result$meta_fits$interpretation <- ifelse(
      result$meta_fits$coef > 0,
      paste("BMI每增加1单位，菌群丰度增加", round(abs(result$meta_fits$coef), 4)),
      paste("BMI每增加1单位，菌群丰度减少", round(abs(result$meta_fits$coef), 4))
    )
  }

  # 生成分析报告 (确保 generate_bmi_report 存在，如果不存在可暂时注释掉)
  analysis_report <- NULL
  if (exists("generate_bmi_report")) {
    analysis_report <- generate_bmi_report(result, metadata, bmi_var, config)
  }

  # 返回结果
  final_result <- list(
    analysis_type = "continuous_bmi",
    n_samples = nrow(metadata),
    n_studies = length(studies),
    meta_result = result,
    significant_features = if (exists("sig_features")) sig_features else NULL,
    report = analysis_report,
    config = config
  )

  return(final_result)
}

#' 执行单研究 BMI 分析 (作为 Meta 分析的 fallback 或单中心分析)
analyze_single_study_bmi <- function(data, metadata, bmi_var, config, covariates) {
  library(MMUPHin) # 虽然是单研究，复用 Maaslin2 逻辑或简单 LM

  cat("   [单研究模式] 使用线性模型分析...\n")

  # 构造公式
  # Feature ~ BMI + Age + Sex
  # 注意：单研究不需要 random effect (batch)

  # 这里我们依然可以使用 MMUPHin::lm_meta，但 batch 设为 NULL 或不传?
  # MMUPHin::lm_meta 要求 batch。如果是单研究，我们可以伪造一个 batch 列，或者直接用 base R lm / Maaslin2
  # 为了保持一致性，我们这里手动循环做 LM

  features <- rownames(data)
  n_features <- length(features)

  # 预分配结果
  pvals <- numeric(n_features)
  coefs <- numeric(n_features)
  se <- numeric(n_features)

  # 转换数据 (AST 或 Log)
  if (config$transform_method == "AST") {
    data_trans <- asin(sqrt(data))
  } else {
    data_trans <- log(data + min(data[data > 0]) / 2)
  }

  # 循环每个特征
  for (i in 1:n_features) {
    y <- as.numeric(data_trans[i, ])

    # 构建数据框
    df_lm <- data.frame(Y = y, BMI = metadata[[bmi_var]])

    # 添加协变量
    for (cov in covariates) {
      if (cov %in% colnames(metadata)) {
        df_lm[[cov]] <- metadata[[cov]]
      }
    }

    # 构建公式
    formula_str <- paste("Y ~ BMI", paste(covariates[covariates %in% colnames(df_lm)], collapse = " + "), sep = " + ")

    # 拟合
    tryCatch(
      {
        fit <- lm(as.formula(formula_str), data = df_lm)
        res <- summary(fit)

        # 提取 BMI 的系数 (假设 BMI 是第二个系数，Intercept 是第一个)
        # 更严谨的提取方法：
        coef_idx <- grep("BMI", rownames(res$coefficients))
        if (length(coef_idx) > 0) {
          coefs[i] <- res$coefficients[coef_idx, "Estimate"]
          pvals[i] <- res$coefficients[coef_idx, "Pr(>|t|)"]
          se[i] <- res$coefficients[coef_idx, "Std. Error"]
        } else {
          pvals[i] <- NA
          coefs[i] <- 0
        }
      },
      error = function(e) {
        pvals[i] <- NA
        coefs[i] <- 0
      }
    )
  }

  # 整理结果表 (模拟 lm_meta 的输出格式)
  res_df <- data.frame(
    coef = coefs,
    pval = pvals,
    stderr = se,
    row.names = features
  )

  # 移除 NA
  res_df <- na.omit(res_df)

  return(list(meta_fits = res_df))
}

# 比较不同疾病组的结果
compare_disease_group_results <- function(all_results) {
  cat("   比较不同疾病组的结果...\n")

  # 收集所有组的显著特征
  all_sig_features <- list()
  all_effects <- data.frame()

  for (group_name in names(all_results)) {
    group_result <- all_results[[group_name]]

    if (!is.null(group_result$significant_features) && nrow(group_result$significant_features) > 0) {
      sig_features <- group_result$significant_features
      sig_features$group <- group_name
      sig_features$feature <- rownames(sig_features)

      all_sig_features[[group_name]] <- sig_features
    }
  }

  # 找出各组共同的显著特征
  common_features <- NULL
  if (length(all_sig_features) > 1) {
    feature_sets <- lapply(all_sig_features, function(x) unique(x$feature))
    common_features <- Reduce(intersect, feature_sets)
    cat("   各组共同的显著特征数量: ", length(common_features), "\n")
  }

  # 找出各组特有的显著特征
  unique_features <- list()
  for (group_name in names(all_sig_features)) {
    other_groups <- setdiff(names(all_sig_features), group_name)
    other_features <- unique(unlist(lapply(all_sig_features[other_groups], function(x) x$feature)))
    group_features <- all_sig_features[[group_name]]$feature
    unique_features[[group_name]] <- setdiff(group_features, other_features)
    cat("   ", group_name, "组特有的显著特征数量: ", length(unique_features[[group_name]]), "\n")
  }

  return(list(
    common_features = common_features,
    unique_features = unique_features
  ))
}

#' 提取显著特征并进行分类
#'
#' @param result_obj 结果对象，如disease_result
#' @param result_name 结果名称，用于变量命名，如"adjBMI"或"unadjBMI"
#' @param data_name 数据集名称，用于变量命名，如"cohort1"或"cohort2"
#' @return 返回一个包含所有特征列表的对象
extract_disease_features <- function(result_obj, result_name = "", data_name = "") {
  # 构建变量名前缀
  prefix <- ""
  if (data_name != "") {
    prefix <- paste0(data_name, "_")
  }
  if (result_name != "") {
    if (prefix != "") {
      prefix <- paste0(prefix, result_name, "_")
    } else {
      prefix <- paste0(result_name, "_")
    }
  }

  # 初始化所有变量
  all_features <- list()

  # 检查结果对象
  if (!is.null(result_obj) && !is.null(result_obj[["significant_features"]])) {
    # 提取比较组特征
    all_features[[paste0(prefix, "CRC_vs_CTRL_features")]] <-
      result_obj[["significant_features"]]$feature[
        result_obj[["significant_features"]]$source_comparison %in% c("CTRL_vs_CRC")
      ]

    all_features[[paste0(prefix, "CTRL_vs_AA_features")]] <-
      result_obj[["significant_features"]]$feature[
        result_obj[["significant_features"]]$source_comparison == "CTRL_vs_AA"
      ]

    all_features[[paste0(prefix, "CRC_vs_AA_features")]] <-
      result_obj[["significant_features"]]$feature[
        result_obj[["significant_features"]]$source_comparison == "AA_vs_CRC"
      ]

    # 上下调特征
    all_features[[paste0(prefix, "CRC_vs_CTRL_up")]] <-
      result_obj[["significant_features"]]$feature[
        result_obj[["significant_features"]]$source_comparison %in% c("CTRL_vs_CRC") &
          result_obj[["significant_features"]]$direction == "up_in_CRC"
      ]

    all_features[[paste0(prefix, "CTRL_vs_AA_up")]] <-
      result_obj[["significant_features"]]$feature[
        result_obj[["significant_features"]]$source_comparison == "CTRL_vs_AA" &
          result_obj[["significant_features"]]$direction == "up_in_AA"
      ]

    all_features[[paste0(prefix, "CRC_vs_AA_up")]] <-
      result_obj[["significant_features"]]$feature[
        result_obj[["significant_features"]]$source_comparison == "AA_vs_CRC" &
          result_obj[["significant_features"]]$direction == "up_in_CRC"
      ]

    all_features[[paste0(prefix, "CRC_vs_CTRL_down")]] <-
      result_obj[["significant_features"]]$feature[
        result_obj[["significant_features"]]$source_comparison %in% c("CTRL_vs_CRC") &
          result_obj[["significant_features"]]$direction == "down_in_CRC"
      ]

    all_features[[paste0(prefix, "CTRL_vs_AA_down")]] <-
      result_obj[["significant_features"]]$feature[
        result_obj[["significant_features"]]$source_comparison == "CTRL_vs_AA" &
          result_obj[["significant_features"]]$direction == "down_in_AA"
      ]

    all_features[[paste0(prefix, "CRC_vs_AA_down")]] <-
      result_obj[["significant_features"]]$feature[
        result_obj[["significant_features"]]$source_comparison == "AA_vs_CRC" &
          result_obj[["significant_features"]]$direction == "down_in_CRC"
      ]

    # 输出统计信息
    cat("=== ", ifelse(data_name != "", paste0(data_name, " "), ""),
      ifelse(result_name != "", paste0("(", result_name, ") "), ""),
      "===\n",
      sep = ""
    )
    cat("CRC vs CTRL特征数:", length(all_features[[paste0(prefix, "CRC_vs_CTRL_features")]]), "\n")
    cat("CTRL vs AA特征数:", length(all_features[[paste0(prefix, "CTRL_vs_AA_features")]]), "\n")
    cat("CRC vs AA特征数:", length(all_features[[paste0(prefix, "CRC_vs_AA_features")]]), "\n")
    cat("\n")
  } else {
    # 初始化空变量
    variable_names <- c(
      "CRC_vs_CTRL_features", "CTRL_vs_AA_features", "CRC_vs_AA_features",
      "CRC_vs_CTRL_up", "CTRL_vs_AA_up", "CRC_vs_AA_up",
      "CRC_vs_CTRL_down", "CTRL_vs_AA_down", "CRC_vs_AA_down"
    )

    for (var_name in variable_names) {
      all_features[[paste0(prefix, var_name)]] <- character(0)
    }

    cat("=== ", ifelse(data_name != "", paste0(data_name, " "), ""),
      ifelse(result_name != "", paste0("(", result_name, ") "), ""),
      "===\n",
      sep = ""
    )
    cat("无显著特征\n\n")
  }

  return(all_features)
}


# 另一种简化版本：如果只需要处理两种固定情况
extract_disease_comparison_features <- function(disease_result) {
  # 这是原始代码的封装版本
  if (!is.null(disease_result) && !is.null(disease_result[["significant_features"]])) {
    # 提取比较组特征
    CRC_vs_CTRL_features <- disease_result[["significant_features"]]$feature[
      disease_result[["significant_features"]]$source_comparison %in% c("CTRL_vs_CRC")
    ]

    CTRL_vs_AA_features <- disease_result[["significant_features"]]$feature[
      disease_result[["significant_features"]]$source_comparison == "CTRL_vs_AA"
    ]

    CRC_vs_AA_features <- disease_result[["significant_features"]]$feature[
      disease_result[["significant_features"]]$source_comparison == "AA_vs_CRC"
    ]

    # 上下调特征
    CRC_vs_CTRL_up <- disease_result[["significant_features"]]$feature[
      disease_result[["significant_features"]]$source_comparison %in% c("CTRL_vs_CRC") &
        disease_result[["significant_features"]]$direction == "up_in_CRC"
    ]

    CTRL_vs_AA_up <- disease_result[["significant_features"]]$feature[
      disease_result[["significant_features"]]$source_comparison == "CTRL_vs_AA" &
        disease_result[["significant_features"]]$direction == "up_in_AA"
    ]

    CRC_vs_AA_up <- disease_result[["significant_features"]]$feature[
      disease_result[["significant_features"]]$source_comparison == "AA_vs_CRC" &
        disease_result[["significant_features"]]$direction == "up_in_CRC"
    ]

    CRC_vs_CTRL_down <- disease_result[["significant_features"]]$feature[
      disease_result[["significant_features"]]$source_comparison %in% c("CTRL_vs_CRC") &
        disease_result[["significant_features"]]$direction == "down_in_CRC"
    ]

    CTRL_vs_AA_down <- disease_result[["significant_features"]]$feature[
      disease_result[["significant_features"]]$source_comparison == "CTRL_vs_AA" &
        disease_result[["significant_features"]]$direction == "down_in_AA"
    ]

    CRC_vs_AA_down <- disease_result[["significant_features"]]$feature[
      disease_result[["significant_features"]]$source_comparison == "AA_vs_CRC" &
        disease_result[["significant_features"]]$direction == "down_in_CRC"
    ]

    cat("CRC vs CTRL特征数:", length(CRC_vs_CTRL_features), "\n")
    cat("CTRL vs AA特征数:", length(CTRL_vs_AA_features), "\n")
    cat("CRC vs AA特征数:", length(CRC_vs_AA_features), "\n")

    return(list(
      CRC_vs_CTRL_features = CRC_vs_CTRL_features,
      CTRL_vs_AA_features = CTRL_vs_AA_features,
      CRC_vs_AA_features = CRC_vs_AA_features,
      CRC_vs_CTRL_up = CRC_vs_CTRL_up,
      CTRL_vs_AA_up = CTRL_vs_AA_up,
      CRC_vs_AA_up = CRC_vs_AA_up,
      CRC_vs_CTRL_down = CRC_vs_CTRL_down,
      CTRL_vs_AA_down = CTRL_vs_AA_down,
      CRC_vs_AA_down = CRC_vs_AA_down
    ))
  } else {
    # 初始化空变量
    return(list(
      CRC_vs_CTRL_features = character(0),
      CTRL_vs_AA_features = character(0),
      CRC_vs_AA_features = character(0),
      CRC_vs_CTRL_up = character(0),
      CTRL_vs_AA_up = character(0),
      CRC_vs_AA_up = character(0),
      CRC_vs_CTRL_down = character(0),
      CTRL_vs_AA_down = character(0),
      CRC_vs_AA_down = character(0)
    ))
  }
}


#' 提取肥胖相关显著特征并进行分类
#'
#' @param result_obj 结果对象，如obesity_result
#' @param result_name 结果名称，用于变量命名，如"adjBMI"或"unadjBMI"
#' @param data_name 数据集名称，用于变量命名，如"cohort1"或"cohort2"
#' @return 返回一个包含所有特征列表的对象
extract_obesity_features <- function(result_obj, result_name = "", data_name = "") {
  # 构建变量名前缀
  prefix <- ""
  if (data_name != "") {
    prefix <- paste0(data_name, "_")
  }
  if (result_name != "") {
    if (prefix != "") {
      prefix <- paste0(prefix, result_name, "_")
    } else {
      prefix <- paste0(result_name, "_")
    }
  }

  # 初始化所有变量
  all_features <- list()

  # 检查结果对象
  if (!is.null(result_obj) && !is.null(result_obj[["significant_features"]])) {
    # 定义比较组
    comparisons <- c("Normal_vs_Ob", "Normal_vs_Ow", "Ow_vs_Ob")

    # 提取每个比较组的特征
    for (comparison in comparisons) {
      # 提取该比较组的全部特征
      comparison_key <- gsub("_vs_", "_vs_", comparison) # 保持原样
      var_name <- gsub("_vs_", "_vs_", comparison) # 作为变量名

      all_features[[paste0(prefix, var_name, "_features")]] <-
        result_obj[["significant_features"]]$feature[
          result_obj[["significant_features"]]$source_comparison == comparison
        ]

      # 确定方向标签
      if (comparison == "Normal_vs_Ow") {
        up_direction <- "up_in_Ow"
        down_direction <- "down_in_Ow"
      } else {
        # Normal_vs_Ob 和 Ow_vs_Ob
        up_direction <- "up_in_Ob"
        down_direction <- "down_in_Ob"
      }

      # 提取上调特征
      all_features[[paste0(prefix, var_name, "_up")]] <-
        result_obj[["significant_features"]]$feature[
          result_obj[["significant_features"]]$source_comparison == comparison &
            result_obj[["significant_features"]]$direction == up_direction
        ]

      # 提取下调特征
      all_features[[paste0(prefix, var_name, "_down")]] <-
        result_obj[["significant_features"]]$feature[
          result_obj[["significant_features"]]$source_comparison == comparison &
            result_obj[["significant_features"]]$direction == down_direction
        ]
    }

    # 输出统计信息
    cat("=== ", ifelse(data_name != "", paste0(data_name, " "), ""),
      ifelse(result_name != "", paste0("(", result_name, ") "), ""),
      "===\n",
      sep = ""
    )

    for (comparison in comparisons) {
      var_name <- gsub("_vs_", "_vs_", comparison)
      total_var <- paste0(prefix, var_name, "_features")
      up_var <- paste0(prefix, var_name, "_up")
      down_var <- paste0(prefix, var_name, "_down")

      cat(paste0(comparison, ":\n"))
      cat("  全部特征数: ", length(all_features[[total_var]]), "\n")
      cat("  上调特征数: ", length(all_features[[up_var]]), "\n")
      cat("  下调特征数: ", length(all_features[[down_var]]), "\n")
    }
    cat("\n")
  } else {
    # 初始化空变量
    comparisons <- c("Normal_vs_Ob", "Normal_vs_Ow", "Ow_vs_Ob")
    suffixes <- c("_features", "_up", "_down")

    for (comparison in comparisons) {
      var_base <- gsub("_vs_", "_vs_", comparison)
      for (suffix in suffixes) {
        all_features[[paste0(prefix, var_base, suffix)]] <- character(0)
      }
    }

    cat("=== ", ifelse(data_name != "", paste0(data_name, " "), ""),
      ifelse(result_name != "", paste0("(", result_name, ") "), ""),
      "===\n",
      sep = ""
    )
    cat("肥胖分析无显著特征\n\n")
  }

  return(all_features)
}

#' 提取肥胖相关显著特征（简化版本，保持原始命名）
#'
#' @param obesity_result 肥胖结果对象
#' @return 返回一个包含所有特征列表的对象
extract_obesity_comparison_features <- function(obesity_result) {
  if (!is.null(obesity_result) && !is.null(obesity_result[["significant_features"]])) {
    # 定义所有比较组
    comparisons <- c("Normal_vs_Ob", "Normal_vs_Ow", "Ow_vs_Ob")

    # 初始化结果列表
    result_list <- list()

    for (comparison in comparisons) {
      # 提取该比较组的全部特征
      comparison_key <- comparison
      var_name <- gsub("_vs_", "_vs_", comparison)

      result_list[[paste0(var_name, "_features")]] <-
        obesity_result[["significant_features"]]$feature[
          obesity_result[["significant_features"]]$source_comparison == comparison
        ]

      # 确定方向标签
      if (comparison == "Normal_vs_Ow") {
        up_direction <- "up_in_Ow"
        down_direction <- "down_in_Ow"
      } else {
        # Normal_vs_Ob 和 Ow_vs_Ob
        up_direction <- "up_in_Ob"
        down_direction <- "down_in_Ob"
      }

      # 提取上调特征
      result_list[[paste0(var_name, "_up")]] <-
        obesity_result[["significant_features"]]$feature[
          obesity_result[["significant_features"]]$source_comparison == comparison &
            obesity_result[["significant_features"]]$direction == up_direction
        ]

      # 提取下调特征
      result_list[[paste0(var_name, "_down")]] <-
        obesity_result[["significant_features"]]$feature[
          obesity_result[["significant_features"]]$source_comparison == comparison &
            obesity_result[["significant_features"]]$direction == down_direction
        ]
    }

    # 输出统计信息
    cat("肥胖相关特征统计:\n")
    for (comparison in comparisons) {
      var_name <- gsub("_vs_", "_vs_", comparison)
      cat(paste0(comparison, ":\n"))
      cat("  全部特征数: ", length(result_list[[paste0(var_name, "_features")]]), "\n")
      cat("  上调特征数: ", length(result_list[[paste0(var_name, "_up")]]), "\n")
      cat("  下调特征数: ", length(result_list[[paste0(var_name, "_down")]]), "\n")
    }
    cat("\n")

    return(result_list)
  } else {
    # 初始化空变量
    comparisons <- c("Normal_vs_Ob", "Normal_vs_Ow", "Ow_vs_Ob")
    result_list <- list()

    for (comparison in comparisons) {
      var_name <- gsub("_vs_", "_vs_", comparison)
      result_list[[paste0(var_name, "_features")]] <- character(0)
      result_list[[paste0(var_name, "_up")]] <- character(0)
      result_list[[paste0(var_name, "_down")]] <- character(0)
    }

    cat("肥胖分析无显著特征\n")
    return(result_list)
  }
}


# 创建完整的分析函数，同时处理疾病和肥胖结果
analyze_all_features <- function(disease_result = NULL, obesity_result = NULL,
                                 result_name = "", data_name = "") {
  # 初始化结果列表
  all_results <- list()

  # 处理疾病结果
  if (!is.null(disease_result)) {
    disease_features <- extract_disease_features(
      result_obj = disease_result,
      result_name = result_name,
      data_name = data_name
    )
    all_results <- c(all_results, disease_features)
  }

  # 处理肥胖结果
  if (!is.null(obesity_result)) {
    obesity_features <- extract_obesity_features(
      result_obj = obesity_result,
      result_name = result_name,
      data_name = data_name
    )
    all_results <- c(all_results, obesity_features)
  }

  return(all_results)
}

run_comprehensive_analysis <- function(
  metaphlan_data, # metaphlan数据矩阵
  metadata, # 元数据
  analysis_type = c("obesity", "disease"), # 分析类型
  config, # 配置列表
  data_name = "custom", # 数据集名称
  disease_order = c("CTRL" = 1, "AA" = 2, "CRC" = 3), # 疾病顺序
  obesity_order = c("Normal" = 1, "Ow" = 2, "Ob" = 3) # 肥胖等级顺序
) {
  # 初始化结果列表
  results <- list()
  results$data_info <- list()
  results$analysis_results <- list()
  results$extracted_features <- list()

  # 设置分析类型
  analysis_type <- match.arg(analysis_type)

  cat("\n=== 开始", data_name, "数据集", analysis_type, "分析 ===\n")

  #### 1.1 数据过滤 ####
  cat("1.1 数据过滤...\n")

  filtering_result <- smart_data_filtering(
    metaphlan_data = t(metaphlan_data),
    metadata = metadata,
    config = config
  )

  # 提取过滤后的数据
  filtered_data <- filtering_result$filtered_data
  metadata_filtered <- metadata[metadata$Sample.Source.ID.processed %in% colnames(filtered_data), ]

  # 存储过滤信息
  results$data_info$filtering_result <- filtering_result
  results$data_info$filtered_data <- filtered_data
  results$data_info$metadata_filtered <- metadata_filtered

  cat("过滤后特征数:", nrow(filtered_data), "\n")
  cat("过滤后样本数:", ncol(filtered_data), "\n")

  #### 1.2 肥胖相关菌群特征分析 ####
  if (analysis_type == "obesity") {
    cat("\n1.2 肥胖相关菌群特征分析...\n")

    # 分离各疾病组数据
    disease_groups <- unique(metadata_filtered$Disease)

    for (disease_group in disease_groups) {
      cat("\n  ", disease_group, "人群分析...\n")

      # 提取当前疾病组数据
      disease_samples <- metadata_filtered$Sample.Source.ID.processed[metadata_filtered$Disease == disease_group]
      disease_metadata <- metadata_filtered[metadata_filtered$Sample.Source.ID.processed %in% disease_samples, ]
      disease_data <- filtered_data[, disease_samples]

      # 检查样本数量
      cat("  ", disease_group, "人群样本数:", length(disease_samples), "\n")
      if (!is.null(disease_metadata$BMI_Group1)) {
        cat("    - Normal组:", sum(disease_metadata$BMI_Group1 == "Normal", na.rm = TRUE), "\n")
        cat("    - Ow组:", sum(disease_metadata$BMI_Group1 == "Ow", na.rm = TRUE), "\n")
        cat("    - Ob组:", sum(disease_metadata$BMI_Group1 == "Ob", na.rm = TRUE), "\n")
      }

      # 跳过样本数不足的组
      if (length(disease_samples) < 5) {
        cat("    警告: 样本数不足，跳过", disease_group, "组分析\n")
        next
      }

      # 执行肥胖相关分析
      obesity_result <- perform_group_analysis_fixed(
        filtered_data = disease_data,
        metadata = disease_metadata,
        config = config,
        group_var = config$BMI_Group_var,
        fixed_order = obesity_order
      )

      # 提取特征
      obesity_features <- extract_obesity_features(
        result_obj = obesity_result,
        result_name = paste("obesity_features", disease_group, sep = "_"),
        data_name = paste(data_name, disease_group, "population", sep = "_")
      )

      # 存储结果
      results$analysis_results[[paste("obesity_analysis", disease_group, sep = "_")]] <- obesity_result
      results$extracted_features[[paste("obesity_features", disease_group, sep = "_")]] <- obesity_features

      # 存储样本信息
      results$data_info[[paste("samples", disease_group, sep = "_")]] <- disease_samples
      results$data_info[[paste("metadata", disease_group, sep = "_")]] <- disease_metadata
      results$data_info[[paste("data", disease_group, sep = "_")]] <- disease_data
    }
  }

  #### 1.3 疾病相关菌群特征分析 ####
  if (analysis_type == "disease") {
    cat("\n1.3 疾病相关菌群特征分析...\n")

    # 1.3.1 整体疾病分析
    cat("  1.3.1 整体疾病分析...\n")

    # 校正BMI的疾病分析
    disease_DA_adjust_BMI <- perform_group_analysis_fixed(
      filtered_data = filtered_data,
      metadata = metadata_filtered,
      config = modifyList(config, list(covariates = c("Age", "Sex", "BMI_Group1"))),
      group_var = "Disease",
      fixed_order = disease_order
    )

    # 不校正BMI的疾病分析
    disease_DA_no_BMI <- perform_group_analysis_fixed(
      filtered_data = filtered_data,
      metadata = metadata_filtered,
      config = modifyList(config, list(covariates = c("Age", "Sex"))),
      group_var = "Disease",
      fixed_order = disease_order
    )

    # 提取特征
    disease_features_adjust_BMI <- extract_disease_features(
      result_obj = disease_DA_adjust_BMI,
      result_name = "disease_features_adjust_BMI",
      data_name = paste(data_name, "overall", sep = "_")
    )

    disease_features_no_BMI <- extract_disease_features(
      result_obj = disease_DA_no_BMI,
      result_name = "disease_features_no_BMI",
      data_name = paste(data_name, "overall", sep = "_")
    )

    # 存储结果
    results$analysis_results$disease_analysis_adjust_BMI <- disease_DA_adjust_BMI
    results$analysis_results$disease_analysis_no_BMI <- disease_DA_no_BMI
    results$extracted_features$disease_features_adjust_BMI <- disease_features_adjust_BMI
    results$extracted_features$disease_features_no_BMI <- disease_features_no_BMI

    # 1.3.2 按BMI分组进行疾病分析
    cat("  1.3.2 按BMI分组进行疾病分析...\n")

    # 检查是否有BMI分组信息
    if (!is.null(metadata_filtered$BMI_Group1)) {
      bmi_groups <- unique(metadata_filtered$BMI_Group1)
      bmi_groups <- bmi_groups[!is.na(bmi_groups)]

      for (bmi_group in bmi_groups) {
        cat("    ", bmi_group, "人群分析...\n")

        # 提取当前BMI组数据
        bmi_samples <- metadata_filtered$Sample.Source.ID.processed[metadata_filtered$BMI_Group1 == bmi_group]
        bmi_metadata <- metadata_filtered[metadata_filtered$Sample.Source.ID.processed %in% bmi_samples, ]
        bmi_data <- filtered_data[, bmi_samples]

        # 检查样本数量
        cat("    ", bmi_group, "人群样本数:", length(bmi_samples), "\n")
        disease_counts <- table(bmi_metadata$Disease)
        for (disease in names(disease_counts)) {
          cat("      - ", disease, "组:", disease_counts[disease], "\n")
        }

        # 跳过样本数不足的组
        if (length(bmi_samples) < 5) {
          cat("      警告: 样本数不足，跳过", bmi_group, "组分析\n")
          next
        }

        # 执行疾病分析
        disease_DA_bmi_group <- perform_group_analysis_fixed(
          filtered_data = bmi_data,
          metadata = bmi_metadata,
          config = config,
          group_var = "Disease",
          fixed_order = disease_order
        )

        # 提取特征
        disease_features_bmi_group <- extract_disease_features(
          result_obj = disease_DA_bmi_group,
          result_name = paste("disease_features", bmi_group, sep = "_"),
          data_name = paste(data_name, bmi_group, "population", sep = "_")
        )

        # 存储结果
        results$analysis_results[[paste("disease_analysis", bmi_group, sep = "_")]] <- disease_DA_bmi_group
        results$extracted_features[[paste("disease_features", bmi_group, sep = "_")]] <- disease_features_bmi_group

        # 存储样本信息
        results$data_info[[paste("samples", bmi_group, sep = "_")]] <- bmi_samples
        results$data_info[[paste("metadata", bmi_group, sep = "_")]] <- bmi_metadata
        results$data_info[[paste("data", bmi_group, sep = "_")]] <- bmi_data
      }
    } else {
      cat("    警告: 元数据中缺少BMI分组信息，跳过按BMI分组的疾病分析\n")
    }
  }

  # 添加分析信息
  results$analysis_info <- list(
    data_name = data_name,
    analysis_type = analysis_type,
    analysis_date = Sys.time(),
    original_sample_count = ncol(metaphlan_data),
    filtered_sample_count = ncol(filtered_data),
    original_feature_count = nrow(metaphlan_data),
    filtered_feature_count = nrow(filtered_data)
  )

  cat("\n=== ", data_name, "数据集", analysis_type, "分析完成 ===\n")
  cat("分析结果已存储在列表中，包含以下组件:\n")
  cat("  - data_info: 数据信息和过滤结果\n")
  cat("  - analysis_results: 详细分析结果\n")
  cat("  - extracted_features: 提取的特征结果\n")
  cat("  - analysis_info: 分析元信息\n")

  return(results)
}

#' 生成分析报告
#'
#' @param all_results 所有比较结果
#' @param significant_features 显著特征数据框
#' @param group_var 分组变量名
#' @return 包含分析报告的列表
generate_analysis_report <- function(all_results, significant_features, group_var) {
  cat("生成分析报告...\n")

  # 初始化报告
  report <- list()
  report$analysis_date <- Sys.time()
  report$group_variable <- group_var
  report$comparisons <- names(all_results)
  report$total_significant_features <- nrow(significant_features)

  # 每个比较的显著特征数
  sig_by_comparison <- table(significant_features$source_comparison)
  report$significant_features_by_comparison <- as.list(sig_by_comparison)

  # 每个比较的样本数
  sample_sizes <- list()
  for (comp_name in names(all_results)) {
    if (!is.null(all_results[[comp_name]]$comparison_info)) {
      sample_sizes[[comp_name]] <- all_results[[comp_name]]$comparison_info$sample_sizes
    }
  }
  report$sample_sizes <- sample_sizes

  # 每个比较的上下调特征数
  if (nrow(significant_features) > 0) {
    direction_counts <- table(significant_features$source_comparison, significant_features$direction)
    report$direction_counts <- as.list(direction_counts)
  } else {
    report$direction_counts <- list()
  }

  # 生成文本摘要
  report$summary_text <- paste(
    "分组变量:", group_var, "\n",
    "比较组:", paste(names(all_results), collapse = ", "), "\n",
    "总显著特征数:", nrow(significant_features), "\n",
    "分析日期:", report$analysis_date, "\n"
  )

  # 打印报告摘要
  cat(report$summary_text)

  return(report)
}
