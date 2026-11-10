##### Beta Diversity Analysis Pipeline #####

library(tidyverse)
library(vegan)

# ==============================================================================
# 0. 加载外部函数库
# ==============================================================================
# 确保 06_beta_diversity_function_code.R 在工作目录下
# source("G:/My Drive/proj_ob_AA_CRC/01.scripts/06_beta_diversity_function_code.R")

# ==============================================================================
# 1. 配置路径与参数
# ==============================================================================
base_dir <- "G:/My Drive/proj_ob_AA_CRC"
input_dir <- file.path(base_dir, "02.processed_data/02_aligned_data")
sp_dir <- file.path(base_dir, "02.processed_data/03_metaphlan_split_data")
output_dir <- file.path(base_dir, "03.results/02_beta_diversity")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- 核心参数 (文件名会包含这些信息) ---
prev_threshold <- 0.1 # 过滤阈值 (10% 样本存在)
trans_method <- "sqrt" # 转换方法 (Hellinger)
dist_method <- "bray" # 距离算法
n_cores_perm <- 15 # PERMANOVA 并行核心数

# ======================== 2. 数据加载 ========================
# 动态搜索最新的数据文件
meta_file <- sort(list.files(input_dir, pattern = "metadata_aligned.*\\.rds", full.names = TRUE), decreasing = TRUE)[1]
sp_file <- sort(list.files(sp_dir, pattern = "metaphlan_species_renormalized.*\\.rds", full.names = TRUE), decreasing = TRUE)[1]

if (is.na(meta_file) || is.na(sp_file)) stop("错误: 未找到输入数据文件，请检查路径。")

message(">>> 正在加载数据...")
metadata <- readRDS(meta_file)
metaphlan_raw <- readRDS(sp_file)

# 确保行名一致
rownames(metadata) <- metadata$Sample.Source.ID.processed
metaphlan <- as.data.frame(t(metaphlan_raw)) # 转置为 行=样本, 列=菌

# ======================== 3. 数据预处理与转换 ========================
message(">>> 开始数据预处理...")

# --- 步骤 A: 普及率过滤 (Prevalence Filtering) ---
# 保留在至少 prev_threshold 比例的样本中存在的菌
keep_taxa <- colSums(metaphlan > 0) > (nrow(metaphlan) * prev_threshold)
metaphlan_filt <- metaphlan[, keep_taxa]

cat(sprintf(
  "   特征过滤: %d -> %d (阈值: %.0f%%)\n",
  ncol(metaphlan), ncol(metaphlan_filt), prev_threshold * 100
))

# --- 步骤 B: 强制再次 TSS 标准化 (Re-normalization) ---
# 关键解释：过滤后每行之和可能不再是 1。
# 必须重新除以行和，使剩余的核心菌群总和为 1。
metaphlan_tss <- metaphlan_filt / rowSums(metaphlan_filt)

# --- 步骤 C: 平方根转换 (Square Root / Hellinger) ---
# Hellinger = sqrt(Relative_Abundance)
message("   执行平方根转换 (Hellinger Transformation)...")
metaphlan_final <- sqrt(metaphlan_tss)

# --- 样本对齐 ---
common_samples <- intersect(rownames(metaphlan_final), rownames(metadata))
metaphlan_final <- metaphlan_final[common_samples, ]
metadata_final <- metadata[common_samples, ]

# 构建分组列 (用于画图)
metadata_final <- metadata_final %>%
  mutate(Group_Label = paste(Disease, BMI_Group1, sep = "_")) %>%
  # 设定因子水平顺序 (Control -> AA -> CRC)
  mutate(Group_Label = factor(Group_Label, levels = c(
    "CTRL_Normal", "CTRL_Ow", "CTRL_Ob",
    "AA_Normal",   "AA_Ow",   "AA_Ob",
    "CRC_Normal",  "CRC_Ow",  "CRC_Ob"
  )))

# ======================== 4. 计算/读取距离矩阵 ========================
# 构建带有参数的文件名，确保可追溯
dist_file_name <- sprintf(
  "dist_matrix_%s_%s_prev%s.rds",
  dist_method, trans_method, prev_threshold
)
dist_file_path <- file.path(output_dir, dist_file_name)

if (file.exists(dist_file_path)) {
  message(sprintf(">>> 检测到已存在的距离矩阵: %s", dist_file_name))
  message("    直接读取，跳过计算...")
  bray_dist <- readRDS(dist_file_path)

  # 简单检查样本是否匹配
  if (!all(rownames(as.matrix(bray_dist)) %in% rownames(metadata_final))) {
    warning("注意: 读取的矩阵样本与当前Metadata不完全匹配，建议删除旧文件重新计算。")
  }
} else {
  message(sprintf(">>> 计算新的距离矩阵 (%s)...", dist_method))
  # 计算 Bray-Curtis
  bray_dist <- vegdist(metaphlan_final, method = dist_method)

  # 保存矩阵
  message(sprintf("    保存距离矩阵到: %s", dist_file_name))
  saveRDS(bray_dist, dist_file_path)
}


# 添加 Overall 分组方便后续分析
metadata_final$Overall <- "Overall"

# ==============================================================================
# 终极调试版函数: 储存完整的 adonis2 原始结果对象
# ==============================================================================

analyze_disease_by_bmi_strata_raw <- function(bray_matrix,
                                              metadata,
                                              disease_factor,
                                              bmi_factors,
                                              subgroup_var,
                                              adjust_vars = c("Age", "Sex"),
                                              strata_var = NULL) {
  require(vegan)
  require(tidyverse)
  require(utils)

  # 初始化空列表
  results_list <- list()

  for (bmi_col in bmi_factors) {
    subgroups <- unique(as.character(metadata[[subgroup_var]]))
    subgroups <- subgroups[!is.na(subgroups)]

    for (sub in subgroups) {
      sub_samples_all <- rownames(metadata)[metadata[[subgroup_var]] == sub & !is.na(metadata[[subgroup_var]])]
      current_bmi_levels <- unique(as.character(metadata[sub_samples_all, bmi_col]))
      current_bmi_levels <- current_bmi_levels[!is.na(current_bmi_levels)]

      for (bmi_lvl in current_bmi_levels) {
        target_samples <- rownames(metadata)[
          metadata[[subgroup_var]] == sub &
            metadata[[bmi_col]] == bmi_lvl &
            !is.na(metadata[[bmi_col]]) &
            rownames(metadata) %in% sub_samples_all
        ]

        if (length(target_samples) < 5) next
        curr_meta <- metadata[target_samples, ]

        diseases <- unique(as.character(curr_meta[[disease_factor]]))
        diseases <- diseases[!is.na(diseases)]
        if (length(diseases) < 2) next

        pairs <- combn(diseases, 2, simplify = FALSE)

        for (pair in pairs) {
          g1 <- pair[1]
          g2 <- pair[2]

          # === Step 1: 数据清洗 (去NA) ===
          pair_samples <- rownames(curr_meta)[curr_meta[[disease_factor]] %in% c(g1, g2)]
          if (length(pair_samples) < 5) next

          pair_meta_temp <- curr_meta[pair_samples, ]

          # 检查 NA
          cols_to_check <- c(disease_factor, adjust_vars)
          if (!is.null(strata_var)) cols_to_check <- c(cols_to_check, strata_var)

          valid_rows <- complete.cases(pair_meta_temp[, cols_to_check])
          pair_meta_clean <- pair_meta_temp[valid_rows, ]
          pair_samples_clean <- rownames(pair_meta_clean)

          if (nrow(pair_meta_clean) < 5) next

          # 同步距离矩阵
          pair_dist_mat <- as.matrix(bray_matrix)[pair_samples_clean, pair_samples_clean]
          pair_dist <- as.dist(pair_dist_mat)

          # === Step 2: 构建公式与 Strata ===
          covariates_str <- paste(adjust_vars, collapse = " + ")
          formula_str <- paste("pair_dist ~", covariates_str, "+", disease_factor)

          perm_strata <- NULL
          if (!is.null(strata_var)) {
            perm_strata <- factor(pair_meta_clean[[strata_var]])
            if (nlevels(perm_strata) < 2) perm_strata <- NULL
          }

          # === Step 3: 运行并储存完整对象 ===
          # 定义一个唯一的名称，方便后续检索
          list_key <- paste(sub, bmi_lvl, g1, "vs", g2, sep = "_")

          tryCatch(
            {
              # 运行模型
              if (!is.null(perm_strata)) {
                permanova_res <- adonis2(as.formula(formula_str),
                  data = pair_meta_clean,
                  permutations = 999,
                  strata = perm_strata,
                  by = "margin"
                )
              } else {
                permanova_res <- adonis2(as.formula(formula_str),
                  data = pair_meta_clean,
                  permutations = 999,
                  by = "margin"
                )
              }

              # === Step 4: 运行 Betadisper (Dispersion) ===
              # 直接算两组之间的 Dispersion 差异
              disp_mod <- betadisper(pair_dist, pair_meta_clean[[disease_factor]])
              disp_perm <- permutest(disp_mod, permutations = 999)

              # 储存结果：包含“信息”, “PERMANOVA原始对象” 和 “Dispersion结果”
              results_list[[list_key]] <- list(
                Info = data.frame(
                  BMI_Criteria = bmi_col,
                  Subgroup_Var = subgroup_var,
                  Subgroup_Val = sub,
                  BMI_Level = bmi_lvl,
                  Group1 = g1,
                  Group2 = g2,
                  N_Samples = nrow(pair_meta_clean),
                  Formula = formula_str,
                  stringsAsFactors = FALSE
                ),
                Raw_Model = permanova_res,
                Dispersion_Perm = disp_perm
              )
            },
            error = function(e) {
              # 如果出错，只存错误信息
              results_list[[list_key]] <- list(
                Info = NULL,
                Raw_Model = NULL,
                Dispersion_Perm = NULL,
                Error = e$message
              )
            }
          )
        }
      }
    }
  }

  return(results_list)
}

# ==============================================================================
# 4. 执行分析 Part 1: 在 BMI 分层下进行疾病两两比较
# ==============================================================================
message("\n=== Part 1: 复杂分层分析 (Disease within BMI strata) ===")

# 请确保上面定义的 analyze_disease_by_bmi_strata 函数已经被运行加载

# 1.1 按 Continent 分层 (例如：China, Europe 下的 Normal/Obese 分别做 CRC vs Control)
results_continent_bmi_disease <- analyze_disease_by_bmi_strata_raw(
  bray_matrix = bray_dist,
  metadata = metadata_final,
  disease_factor = "Disease",
  bmi_factors = c("BMI_Group1"),
  subgroup_var = "Continent",
  adjust_vars = c("Age", "Sex"),
  strata_var = "Study.name"
)

print("按洲分层结果预览：")
print(results_continent_bmi_disease)

# 1.2 整体分析 (Overall) - 修正了你原始代码中的变量名错误
results_overall_bmi_disease <- analyze_disease_by_bmi_strata_raw(
  bray_matrix = bray_dist,
  metadata = metadata_final,
  disease_factor = "Disease",
  bmi_factors = c("BMI_Group1"),
  subgroup_var = "Overall",
  adjust_vars = c("Age", "Sex"),
  strata_var = "Study.name"
)

print("整体分析结果预览：")
print(results_overall_bmi_disease)
library(tidyverse)

# ==============================================================================
# 任务 1: 提取第三行 (Disease) 的结果
# ==============================================================================
# 定义提取函数
extract_disease_only <- function(res_list) {
  # 遍历列表中的每一个元素
  results <- lapply(names(res_list), function(key) {
    # 获取原始模型对象
    # 注意：根据你的描述，结构是 item[["Raw_Model"]]
    item <- res_list[[key]]

    # 安全检查：确保 item 有内容且有 Raw_Model
    if (is.null(item) || is.null(item$Raw_Model)) {
      return(NULL)
    }

    model <- item$Raw_Model

    # --- 核心提取逻辑：提取第 3 行 ---
    # 通常第3行是主效应 (Disease)，建议同时检查行名是否符合预期
    # 这里我们严格执行"提取第3行"
    target_row <- model[3, , drop = FALSE]

    # 添加辅助列
    target_row$Term <- rownames(target_row) # 将行名变成一列
    target_row$Comparison_ID <- key # 把列表名作为ID

    return(target_row)
  })

  # 合并所有结果
  final_df <- bind_rows(results)

  # 整理列顺序 (ID 放最前)
  final_df <- final_df %>%
    select(Comparison_ID, Term, everything())

  return(final_df)
}

# 执行提取
df_disease_overall <- extract_disease_only(results_overall_bmi_disease)
df_disease_continent <- extract_disease_only(results_continent_bmi_disease)
df_disease_result <- rbind(df_disease_overall, df_disease_continent)
# 可选：将 Comparison_ID 拆分为更细致的列 (Subgroup, BMI, Groups)
# 假设命名格式严格为 "Overall_Ow_CRC_vs_CTRL" (下划线分隔)
# 注意：如果 Disease 名字里带下划线，separate 可能会切错，这里仅做简单切分演示
df_disease_result <- df_disease_result %>%
  separate(Comparison_ID, c("Subgroup", "BMI", "Group1", "vs", "Group2"), sep = "_", remove = FALSE) %>%
  select(-vs) # 去掉 'vs' 这一列

df_disease_result <- df_disease_result %>%
  group_by(Subgroup, BMI) %>% # 按亚组分组
  mutate(P_Adj = p.adjust(`Pr(>F)`, method = "fdr")) %>%
  ungroup()

# 预览
print("=== 任务1: Disease 结果预览 ===")
print(df_disease_result)


# ==============================================================================
# 任务 2: 提取前三行 (包括矫正因素) 的结果
# ==============================================================================
extract_full_model <- function(res_list) {
  results <- lapply(names(res_list), function(key) {
    item <- res_list[[key]]
    if (is.null(item) || is.null(item$Raw_Model)) {
      return(NULL)
    }

    model <- item$Raw_Model

    # --- 核心提取逻辑：提取前 3 行 ---
    # 假设顺序是 Age, Sex, Disease
    target_rows <- model[1:3, , drop = FALSE]

    # 添加辅助列
    target_rows$Term <- rownames(target_rows)
    target_rows$Comparison_ID <- key

    return(target_rows)
  })

  final_df <- bind_rows(results)

  # 整理列顺序
  final_df <- final_df %>%
    select(Comparison_ID, Term, everything())

  return(final_df)
}

# 执行提取
df_full_model_result <- extract_full_model(results_overall_bmi_disease)
df_full_model_result_continent <- extract_full_model(results_continent_bmi_disease)
df_full_model_result <- rbind(df_full_model_result, df_full_model_result_continent)

# 同样可以进行拆分 ID
df_full_model_result <- df_full_model_result %>%
  separate(Comparison_ID, c("Subgroup", "BMI", "Group1", "vs", "Group2"), sep = "_", remove = FALSE) %>%
  select(-vs)


# ==============================================================================
# 任务 3: 提取 Betadisper (Dispersion) 的结果
# ==============================================================================
extract_dispersion_results <- function(res_list) {
  results <- lapply(names(res_list), function(key) {
    item <- res_list[[key]]
    if (is.null(item) || is.null(item$Dispersion_Perm)) {
      return(NULL)
    }

    disp_perm <- item$Dispersion_Perm

    # 从 permutest 结果中提取 p-value
    disp_p <- disp_perm$tab$`Pr(>F)`[1]
    disp_f <- disp_perm$tab$F[1]

    res_df <- data.frame(
      Comparison_ID = key,
      Dispersion_F = disp_f,
      Dispersion_P = disp_p,
      stringsAsFactors = FALSE
    )

    return(res_df)
  })

  final_df <- bind_rows(results)
  return(final_df)
}

# 执行提取
df_disp_overall <- extract_dispersion_results(results_overall_bmi_disease)
df_disp_continent <- extract_dispersion_results(results_continent_bmi_disease)
df_disp_result <- rbind(df_disp_overall, df_disp_continent)

# 拆分 ID 并合并到 Disease 结果中
df_disp_result <- df_disp_result %>%
  separate(Comparison_ID, c("Subgroup", "BMI", "Group1", "vs", "Group2"), sep = "_", remove = FALSE) %>%
  select(-vs)

# 预览
print("=== 任务3: Betadisper 结果预览 ===")
print(head(df_disp_result))

# ==============================================================================
# 保存结果
# ==============================================================================
write.csv(df_disease_result, file = file.path(output_dir, "beta_diversity_disease_only.csv"), row.names = FALSE)
write.csv(df_full_model_result, file.path(output_dir, "beta_diversity_full_model_adjusted.csv"), row.names = FALSE)
write.csv(df_disp_result, file.path(output_dir, "beta_diversity_dispersion_pairwise.csv"), row.names = FALSE)

message("分析完成！结果已保存至: ", output_dir)
