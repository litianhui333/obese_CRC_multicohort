# ==============================================================================
# 脚本名称: Run_PERMANOVA_By_Continent.R
# 功能描述: 
#   1. 加载并预处理数据 (与 06_beta_pairwise...R 保持严格一致)
#   2. 执行分层 PERMANOVA 分析 (按 Continent 分组) + 总体分析 (Overall)
#   3. 可视化结果
# ==============================================================================

library(vegan)
library(tidyverse)
library(dplyr)

# ======================== 1. 核心函数定义 ========================

# --- 函数 A: 核心 PERMANOVA (Robust Version) ---
run_permanova_robust <- function(bray_dist, factors, metadata, n_cores = 15) {
  
  # 确保 metadata 和 distance matrix 的行名/样本名完全一致
  bray_matrix <- as.matrix(bray_dist)
  common_samples <- intersect(rownames(metadata), rownames(bray_matrix))
  
  if(length(common_samples) < nrow(metadata) || length(common_samples) < nrow(bray_matrix)){
    # warning(paste("注意：Metadata 和 距离矩阵样本不完全匹配。将仅使用共有的", length(common_samples), "个样本。"))
  }
  
  metadata <- metadata[common_samples, , drop=FALSE]
  bray_matrix <- bray_matrix[common_samples, common_samples]
  
  # 修正点：将 StringsAsFactors 改为小写的 stringsAsFactors
  results_df <- data.frame(
    Factor = character(),
    R2 = numeric(),
    P_value = numeric(),
    F_Model = numeric(),
    N_Samples = integer(),
    stringsAsFactors = FALSE
  )
  
  factors <- factors[factors %in% colnames(metadata)]
  
  for (factor in factors) {
    # cat("  正在分析:", factor, "...\n")
    valid_samples <- rownames(metadata)[!is.na(metadata[[factor]])]
    
    if (length(valid_samples) < 5) next
    
    sub_metadata <- metadata[valid_samples, , drop = FALSE]
    sub_dist <- as.dist(bray_matrix[valid_samples, valid_samples])
    form <- as.formula(paste("sub_dist ~", factor))
    
    # 智能选择 Strata
    if (factor %in% c("Study.name", "Continent", "Country")) {
      perm_ctrl <- 999 
    } else {
      if ("Study.name" %in% colnames(sub_metadata) && length(unique(sub_metadata$Study.name)) > 1) {
        perm_ctrl <- how(nperm = 999, blocks = as.factor(sub_metadata$Study.name))
      } else {
        perm_ctrl <- 999
      }
    }
    
    tryCatch({
      res <- adonis2(form, data = sub_metadata, permutations = perm_ctrl, parallel = n_cores)
      stats <- data.frame(
        Factor = factor,
        R2 = res$R2[1],
        P_value = res$`Pr(>F)`[1],
        F_Model = res$F[1],
        N_Samples = length(valid_samples)
      )
      results_df <- rbind(results_df, stats)
    }, error = function(e) { })
  }
  
  if(nrow(results_df) > 0) {
    results_df$P_adj <- p.adjust(results_df$P_value, method = "BH")
  }
  
  return(results_df)
}

# --- 函数 B: 分层运行 Wrapper (Wrapper Function) ---
run_permanova_by_region <- function(bray_dist, factors, metadata, 
                                    split_col = "Continent", 
                                    target_regions = c("Asia", "Europe"),
                                    include_overall = TRUE,
                                    n_cores = 15) {
  
  all_region_results <- list()
  bray_matrix_full <- as.matrix(bray_dist)
  
  # --- 1. Overall 分析 (全数据集) ---
  if (include_overall) {
    cat("\n=======================================\n")
    cat("正在启动总体分析 (Overall)\n")
    cat("=======================================\n")
    
    # 检查样本量
    if (nrow(metadata) >= 10) {
      # 对于 Overall 分析，使用所有传入的 factors
      # 如果 split_col (如 Continent) 在 factors 里，也会被计算，这是合理的
      overall_res <- run_permanova_robust(
        bray_dist = as.dist(bray_matrix_full),
        factors = factors,
        metadata = metadata,
        n_cores = n_cores
      )
      
      if (!is.null(overall_res) && nrow(overall_res) > 0) {
        overall_res$Region <- "Overall"
        all_region_results[["Overall"]] <- overall_res
      }
    } else {
      cat("  -> 警告：总体样本量过少 (<10)，跳过 Overall 分析。\n")
    }
  }
  
  # --- 2. 分层分析 (Regional) ---
  for (region in target_regions) {
    cat("\n=======================================\n")
    cat("正在启动分层分析: ", region, "\n")
    cat("=======================================\n")
    
    target_idx <- which(metadata[[split_col]] == region)
    target_samples <- rownames(metadata)[target_idx]
    
    if (length(target_samples) < 10) {
      cat("  -> 警告：样本量过少 (<10)，跳过该地区。\n")
      next
    }
    
    sub_metadata <- metadata[target_samples, , drop = FALSE]
    sub_dist_mat <- bray_matrix_full[target_samples, target_samples]
    sub_dist <- as.dist(sub_dist_mat)
    
    # 动态调整因子列表：移除分组变量本身 + 移除单一水平变量
    current_factors <- setdiff(factors, split_col)
    valid_factors <- c()
    for (f in current_factors) {
      # 移除 NA 后检查唯一值数量
      if (length(unique(na.omit(sub_metadata[[f]]))) > 1) {
        valid_factors <- c(valid_factors, f)
      } else {
        cat("  -> 变量", f, "在该地区只有单一水平（或全NA），已移除。\n")
      }
    }
    
    region_result <- run_permanova_robust(
      bray_dist = sub_dist, 
      factors = valid_factors, 
      metadata = sub_metadata, 
      n_cores = n_cores
    )
    
    if (!is.null(region_result) && nrow(region_result) > 0) {
      region_result$Region <- region
      all_region_results[[region]] <- region_result
    }
  }
  
  final_combined_df <- do.call(rbind, all_region_results)
  
  # 调整 Factor 顺序，让 Overall 排在最前面 (可选)
  if("Overall" %in% final_combined_df$Region) {
    final_combined_df$Region <- factor(final_combined_df$Region, levels = c("Overall", target_regions))
  }
  
  return(final_combined_df)
}

# ======================== 2. 数据加载与预处理 (一致性流程) ========================

# --- 2.1 路径设置 (与上传的脚本一致) ---
base_dir   <- "D:/Tianhui/proj_ob_AA_CRC"
input_dir  <- file.path(base_dir, "02.processed_data/02_aligned_data")
sp_dir     <- file.path(base_dir, "02.processed_data/03_metaphlan_split_data")
output_dir <- file.path(base_dir, "03.results/02_beta_diversity/")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 参数设置
prev_threshold <- 0.1      # 过滤阈值
trans_method   <- "sqrt"   # 转换方法
dist_method    <- "bray"   # 距离算法
n_cores_perm   <- 15

# --- 2.2 加载原始数据 ---
meta_file <- sort(list.files(input_dir, pattern = "metadata_aligned.*\\.rds", full.names = TRUE), decreasing = TRUE)[1]
sp_file   <- sort(list.files(sp_dir, pattern = "metaphlan_species_renormalized.*\\.rds", full.names = TRUE), decreasing = TRUE)[1]

if (is.na(meta_file) || is.na(sp_file)) stop("错误: 未找到输入数据文件，请检查路径。")

message(">>> 正在加载数据...")
metadata <- readRDS(meta_file)
metaphlan_raw <- readRDS(sp_file)

# 设置行名 (保持一致性)
rownames(metadata) <- metadata$Sample.Source.ID.processed
metaphlan <- as.data.frame(t(metaphlan_raw)) # 转置: 行=样本

# --- 2.3 预处理 (Filtering -> TSS -> Hellinger) ---
message(">>> 执行数据预处理...")

# A. Prevalence Filtering
keep_taxa <- colSums(metaphlan > 0) > (nrow(metaphlan) * prev_threshold)
metaphlan_filt <- metaphlan[, keep_taxa]
cat(sprintf("   特征过滤: %d -> %d (阈值: %.0f%%)\n", ncol(metaphlan), ncol(metaphlan_filt), prev_threshold * 100))

# B. Re-normalization (TSS)
metaphlan_tss <- metaphlan_filt / rowSums(metaphlan_filt)

# C. Transformation (Hellinger/Sqrt)
metaphlan_final <- sqrt(metaphlan_tss)

# D. 样本对齐
common_samples  <- intersect(rownames(metaphlan_final), rownames(metadata))
metaphlan_final <- metaphlan_final[common_samples, ]
metadata_final  <- metadata[common_samples, ]

## 去掉这个对列，因为只包含female会导致偏差
metadata_final <- metadata_final %>%
  filter(!Study.name == "This_study_cohort5__NSHII")
# --- 2.4 获取距离矩阵 (读取或计算) ---
dist_file_name <- sprintf("dist_matrix_%s_%s_prev%s_remove_one_study.rds", dist_method, trans_method, prev_threshold)
dist_file_path <- file.path(output_dir, dist_file_name)

if (file.exists(dist_file_path)) {
  message(sprintf(">>> 读取已存在的距离矩阵: %s", dist_file_name))
  bray_dist <- readRDS(dist_file_path)
} else {
  message(sprintf(">>> 计算新的距离矩阵 (%s)...", dist_method))
  bray_dist <- vegdist(metaphlan_final, method = dist_method)
  saveRDS(bray_dist, dist_file_path)
}

# ======================== 3. 运行分层分析 ========================

# 定义分析因子 (建议根据你的实际列名调整)
# 注意：不需要手动移除 "Continent"，函数会自动处理
my_factors <- c("Age", "Sex", "BMI", "BMI_Group1","Disease", "Study.name") 

# 额外处理：如果需要标准化 BMI 可以在这里做
metadata_final$BMI_scaled <- scale(metadata_final$BMI) 

# 运行函数
stratified_results <- run_permanova_by_region(
  bray_dist = bray_dist,
  factors = my_factors,
  metadata = metadata_final,
  split_col = "Continent",          
  target_regions = c("Asia", "Europe"),
  include_overall = TRUE,            # 开启 Overall 计算
  n_cores = n_cores_perm
)
# 创建保存子目录
sub_dir <- file.path(output_dir, "permanova_subgroups_univariate")
if(!dir.exists(sub_dir)) dir.create(sub_dir, recursive = TRUE)
save(stratified_results, file = file.path(sub_dir, "PERMANOVA_By_Continent_Results_with_Overall_remove_one_study_20260122.RData"))
write.csv(stratified_results,file = file.path(sub_dir, "PERMANOVA_By_Continent_Results_with_Overall_remove_one_study_20260122.csv"))
# ======================== 4. 统计验证: 多因素 PERMANOVA ========================
message("\n>>> 正在启动多分层 PERMANOVA 统计验证 (By Continent & BMI type)...")

# ------------------------------------------------------------------------------
# 1. 参数与容器准备
# ------------------------------------------------------------------------------
target_regions <- c("Overall", "Asia", "Europe")
bmi_test_vars  <- c("BMI", "BMI_Group1") # 对比连续型和分类型
permanova_results_full <- list()

# 创建保存子目录
sub_dir <- file.path(output_dir, "permanova_subgroups_multi")
if(!dir.exists(sub_dir)) dir.create(sub_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 2. 嵌套循环分析
# ------------------------------------------------------------------------------
for (region in target_regions) {
  for (bmi_var in bmi_test_vars) {
    
    message(sprintf("\n>>> 正在分析 [区域: %s] | [BMI变量: %s] ...", region, bmi_var))
    
    # 2.1 提取并清洗子集数据
    if (region == "Overall") {
      meta_sub <- metadata_final
    } else {
      meta_sub <- metadata_final %>% filter(Continent == region)
    }
    
    # 定义当前模型需要的列
    current_vars <- c("Study.name", "Age", "Sex", "Disease", bmi_var)
    
    # 剔除缺失值 (Complete cases)
    valid_idx <- complete.cases(meta_sub[, current_vars])
    meta_clean <- meta_sub[valid_idx, ]
    
    if (nrow(meta_clean) < 20) {
      message(sprintf("  跳过 %s: 有效样本量不足 (%d)", region, nrow(meta_clean)))
      next
    }
    
    # 2.2 提取对应的距离矩阵
    # 注意：必须从全局距离矩阵中，根据清洗后的 SampleID 提取子集
    sample_ids <- rownames(meta_clean)
    dist_mat_all <- as.matrix(bray_dist)
    dist_clean <- as.dist(dist_mat_all[sample_ids, sample_ids])
    
    # 2.3 准备 Stratification (批次控制)
    # 如果该子集只包含一个 Study，则不设置 strata
    blocks <- if(length(unique(meta_clean$Study.name)) > 1) as.factor(meta_clean$Study.name) else NULL
    
    # 2.4 构建公式
    # 基础公式：~ Age + Sex + BMI_var + Disease
    # 交互公式：~ Age + Sex + BMI_var * Disease
    formula_main <- as.formula(paste("dist_clean ~ Age + Sex +", bmi_var, "+ Disease"))
    formula_inter <- as.formula(paste("dist_clean ~ Age + Sex +", bmi_var, "* Disease"))
    
    # 2.5 运行 adonis2
    tryCatch({
      # 模型 1: 主效应 (Type III SS)
      perm_main <- adonis2(formula_main, data = meta_clean, permutations = 999,
                           strata = blocks, by = "margin", parallel = n_cores_perm)
      
      # 模型 2: 交互项
      perm_inter <- adonis2(formula_inter, data = meta_clean, permutations = 999,
                            strata = blocks, by = "margin", parallel = n_cores_perm)
      
      # 2.6 提取关键统计量并存储
      res_entry <- data.frame(
        Region = region,
        BMI_Type = bmi_var,
        N_Samples = nrow(meta_clean),
        N_Studies = length(unique(meta_clean$Study.name)),
        P_BMI = perm_main[bmi_var, "Pr(>F)"],
        R2_BMI = perm_main[bmi_var, "R2"],
        P_Disease = perm_main["Disease", "Pr(>F)"],
        R2_Disease = perm_main["Disease", "R2"],
        P_Interaction = perm_inter[paste0(bmi_var, ":Disease"), "Pr(>F)"],
        R2_Interaction = perm_inter[paste0(bmi_var, ":Disease"), "R2"]
      )
      
      permanova_results_full[[paste(region, bmi_var, sep="_")]] <- res_entry
      
      # 保存详细结果到文本文件
      out_txt <- file.path(sub_dir, paste0("PERMANOVA_", region, "_", bmi_var, ".txt"))
      sink(out_txt)
      cat("Region:", region, "\nBMI Variable:", bmi_var, "\n")
      print(perm_main)
      print(perm_inter)
      sink()
      
    }, error = function(e) {
      message(sprintf("  分析 %s_%s 时出错: %s", region, bmi_var, e$message))
    })
  }
}

# ------------------------------------------------------------------------------
# 3. 结果汇总输出
# ------------------------------------------------------------------------------
final_summary_table <- bind_rows(permanova_results_full)

# 添加显著性标记
final_summary_table <- final_summary_table %>%
  mutate(Sig_Inter = case_when(
    P_Interaction < 0.001 ~ "***",
    P_Interaction < 0.01 ~ "**",
    P_Interaction < 0.05 ~ "*",
    TRUE ~ "ns"
  ))

write.csv(final_summary_table, file.path(sub_dir, "PERMANOVA_Subgroup_Comparison_Summary.csv"), row.names = FALSE)

print("--- 所有子组 PERMANOVA 分析完成！汇总表已生成 ---")
print(final_summary_table)

# ======================== 4. 结果展示与可视化 ========================

message("\n>>> 分析完成，结果概览 (按 Region 和 R2 排序):")
## 单因素barplot
library(tidyverse)
library(ggsci)

# 2. 数据清洗与准备
plot_data <- stratified_results %>%
  # 过滤 Study.name (R2 过大)
  filter(!Factor %in% c("Study.name")) %>%
  # 生成显著性标签
  mutate(
    Significance = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01  ~ "**",
      P_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ),
    Label = ifelse(Significance == "ns", "", Significance),
    R2_Percent = R2 * 100
  ) %>%
  # 设置顺序
  mutate(Factor = factor(Factor, levels = c("Disease", "BMI_Group1", "BMI", "Age", "Sex"))) %>%
  mutate(Region = factor(Region, levels = c("Overall", "Asia", "Europe")))

# 3. 绘图
p <- ggplot(plot_data, aes(x = Factor, y = R2_Percent, fill = Factor)) +
  # 柱状图
  geom_bar(stat = "identity", width = 0.7, color = "black", size = 0.3) +
  # 显著性标签
  geom_text(aes(label = Label), vjust = -0.5, size = 5, fontface = "bold") +
  # 分面
  facet_wrap(~ Region, scales = "fixed") + 
  # 配色与主题
  scale_fill_npg() +
  labs(title = "Variance Explained (R2) by Factors",
       subtitle = "PERMANOVA Analysis | Adjusted by Study",
       y = "R2 (%)", x = NULL,
       caption = "*** P < 0.001, ** P < 0.01, * P < 0.05") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    strip.background = element_rect(fill = "#EFEFEF"),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)))

# 保存
 ggsave(paste0(output_dir,"/permanova_subgroups_univariate/PERMANOVA_R2_Plot_uni.pdf"), p, width = 5.4, height = 3.24)
 
##### 多变量画图
 library(tidyverse)
 library(ggsci)
 # 假设你的原始数据框叫 final_summary_table
 plot_data <- final_summary_table %>%
   # 1. 提取所有 R2 和 P 值，转为长格式
   pivot_longer(cols = starts_with(c("R2_", "P_")), 
                names_to = c(".value", "Term"), 
                names_sep = "_") %>%
   # 2. 核心逻辑：根据 BMI_Type 重新定义 Factor 名字
   mutate(Factor = case_when(
     Term == "Disease" ~ "Disease",
     Term == "BMI" & BMI_Type == "BMI" ~ "BMI",
     Term == "BMI" & BMI_Type == "BMI_Group1" ~ "BMI_Group1",
     Term == "Interaction" & BMI_Type == "BMI" ~ "Intersection(Continous)",
     Term == "Interaction" & BMI_Type == "BMI_Group1" ~ "Intersection(Category)",
     TRUE ~ NA_character_
   )) %>%
   # 3. 去掉不需要的重复行 (例如 BMI_Type 为 BMI 时产生的 Disease 行)
   filter(!is.na(Factor)) %>%
   # 如果同一个 Region 下 Disease 出现了两次（因为 BMI 和 BMI_Group1 各贡献一次），取均值即可
   group_by(Region, Factor) %>%
   summarise(R2 = mean(R2), P_value = mean(P), .groups = "drop") %>%
   # 4. 格式化用于绘图
   mutate(
     R2_Percent = R2 * 100,
     Significance = case_when(
       P_value < 0.001 ~ "***",
       P_value < 0.01  ~ "**",
       P_value < 0.05  ~ "*",
       TRUE            ~ "ns"
     ),
     Label = ifelse(Significance == "ns", "", Significance),
     Factor = factor(Factor, levels = c("Disease", "BMI", "BMI_Group1", 
                                        "Intersection(Continous)", "Intersection(Category)")),
     Region = factor(Region, levels = c("Overall", "Asia", "Europe"))
   )
 
 
 p <- ggplot(plot_data, aes(x = Factor, y = R2_Percent, fill = Factor)) +
   # 绘制柱状图
   geom_bar(stat = "identity", width = 0.7, color = "black", size = 0.3) +
   
   # 添加显著性星号 (P < 0.05 才会显示)
   geom_text(aes(label = Label), vjust = -0.5, size = 5, fontface = "bold") +
   
   # 按地域分面
   facet_wrap(~ Region, scales = "fixed") + 
   
   # 使用 NPG 杂志配色
   scale_fill_npg() +
   
   # 坐标轴与标题
   labs(
     title = "Variance Explained (R2) by Clinical Factors",
     subtitle = "Multivariate PERMANOVA | Global vs Regional Comparison",
     y = "R2 (%)", 
     x = NULL,
     caption = "Notes: *** P < 0.001, ** P < 0.01, * P < 0.05. \n'Intersection' represents the interaction effect of BMI and Disease."
   ) +
   
   # 主题美化
   theme_bw() +
   theme(
     # 标题样式
     plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
     plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
     
     # X 轴文字旋转 45 度，防止标签重叠
     axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold", color = "black"),
     axis.text.y = element_text(size = 10, color = "black"),
     axis.title.y = element_text(size = 12, face = "bold"),
     
     # 分面条目样式
     strip.background = element_rect(fill = "#EFEFEF", color = NA),
     strip.text = element_text(face = "bold", size = 12),
     
     # 隐藏图例（因为 X 轴已有标签）
     legend.position = "none",
     
     # 去掉垂直网格线
     panel.grid.major.x = element_blank(),
     panel.grid.minor = element_blank()
   ) +
   
   # 扩展 Y 轴上限，防止显著性标签溢出
   scale_y_continuous(expand = expansion(mult = c(0, 0.2)))
 
 # ======================== 3. 保存与查看 ========================
 
 print(p)
 
 # 根据需要保存为 PDF 或图片
 ggsave(filename = file.path(sub_dir,"PERMANOVA_R2_Comparison_Final.pdf"), plot = p, width = 5.4, height = 4)
 