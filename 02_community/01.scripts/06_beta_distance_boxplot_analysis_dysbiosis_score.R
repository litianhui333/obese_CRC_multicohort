# ==============================================================================
# 脚本名称: 06_beta_distance_boxplot_analysis.R
# 功能描述:
#   1. 预处理: 普及率过滤 -> TSS标准化 (Re-normalization) -> 平方根转化 (Hellinger)
#   2. 距离计算: Bray-Curtis 距离 (支持缓存，避免重复计算)
#   3. 核心分析: 在 Study 内部计算 Dysbiosis Score (距离健康基线的程度)
#   4. 可视化: 绘制跨 Study 的 Boxplot 并做统计检验
#   5. 统计验证: 多因素 PERMANOVA (验证独立效应与交互作用)
# ==============================================================================

library(tidyverse)
library(vegan)
library(ggpubr) # 用于统计检验和绘图
library(ggsci) # 优美的配色

# ======================== 1. 路径与参数设置 ========================
base_dir <- "G:/My Drive/proj_ob_AA_CRC"
input_dir <- file.path(base_dir, "02.processed_data/02_aligned_data")
sp_dir <- file.path(base_dir, "02.processed_data/03_metaphlan_split_data")
output_dir <- file.path(base_dir, "03.results/04_distance_analysis")

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
dim(as.matrix(bray_dist))
# ======================== 5. 核心分析: 计算 Dysbiosis Score ========================
message(">>> 开始计算 Study 内部 Dysbiosis Score...")

# 转换距离为矩阵方便索引
dist_mat <- as.matrix(bray_dist)

# 定义参考组条件: Control & Normal
ref_disease <- "CTRL"
ref_bmi <- "Normal"

# 结果容器
score_list <- list()
all_studies <- unique(metadata_final$Study.name)

for (std in all_studies) {
  # 1. 找出该 Study 下的所有样本
  study_samples <- rownames(metadata_final)[metadata_final$Study.name == std]

  # 2. 找出该 Study 下的参考样本 (Healthy Normal)
  ref_samples <- rownames(metadata_final)[
    metadata_final$Study.name == std &
      metadata_final$Disease == ref_disease &
      metadata_final$BMI_Group1 == ref_bmi
  ]

  # 3. 检查参考组是否足够
  if (length(ref_samples) < 5) {
    message(sprintf("   [跳过] Study: %s - 参考样本不足 (n=%d)", std, length(ref_samples)))
    next
  }

  # 4. 计算每个样本的 Score
  # Score = 该样本 到 所有参考样本距离 的平均值
  study_scores <- sapply(study_samples, function(sid) {
    # 提取距离向量
    dists <- dist_mat[sid, ref_samples]
    # 如果样本自己就是参考样本，去掉自对自距离(0)
    dists <- dists[dists > 1e-6]
    if (length(dists) == 0) {
      return(0)
    } # 极端情况保护
    return(median(dists))
  })

  # 5. 存入结果
  score_list[[std]] <- data.frame(
    SampleID = names(study_scores),
    Dysbiosis_Score = study_scores,
    Study = std
  )
}

# [跳过] Study: This_study_cohort1__ONCOBIOME_AtezoTRIBE - 参考样本不足 (n=0)
# [跳过] Study: This_study_cohort2__ONCOBIOME_COLOBIOME - 参考样本不足 (n=0)

# 合并所有结果
dysbiosis_df <- do.call(rbind, score_list)
message(sprintf(
  ">>> 计算完成，共处理 %d 个研究，生成 %d 个样本分数。",
  length(score_list), nrow(dysbiosis_df)
))

# 合并回元数据用于画图
plot_data <- metadata_final %>%
  inner_join(dysbiosis_df, by = c("Sample.Source.ID.processed" = "SampleID", "Study.name" = "Study"))

# 保存分数表
write.csv(plot_data, file.path(output_dir, "dysbiosis_scores_table.csv"), row.names = FALSE)

# ======================== 6. 可视化: Dysbiosis Score Boxplot ========================
message(">>> 正在绘图...")

# 定义比较组
my_comparisons <- list(
  c("CTRL_Normal", "CTRL_Ob"),
  c("CTRL_Normal", "AA_Normal"),
  c("CTRL_Normal", "CRC_Normal"),
  c("CTRL_Ob", "CRC_Ob")
)

# 计算每个 Study 的平均值 (用于画点)
study_means <- plot_data %>%
  group_by(Study.name, Group_Label, Disease) %>%
  summarise(
    Mean_Score = mean(Dysbiosis_Score, na.rm = TRUE),
    SE_Score = sd(Dysbiosis_Score, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# 主图: Boxplot (样本分布) + Jitter Points (Study Means)
p_box <- ggplot(plot_data, aes(x = Group_Label, y = Dysbiosis_Score)) +

  # 1. Boxplot (展示整体样本分布) - 设为背景，透明度高
  geom_boxplot(aes(fill = Disease), alpha = 0.5, outlier.shape = NA, width = 0.6) +

  # 2. Study Means (每个 Study 一个点) - Color 区分 Study, Shape 区分 Study (可选)
  # 注意: 这里的点代表该 Study 在该 Group 下的平均 Dysbiosis Score
  # 建议: 如果 Study 很多，用颜色区分可能太杂。但用户问 "是否每个study为一个点"，所以必须展示 Study 维度。
  # 这里用 shape 区分 Study (最多支持 ~6-10 种形状，多了会重复)，如果不合适，可以用 point size mapping total N?
  # 用户明确提到 "figure legend 缺失"，说明他希望能看到 legend。
  # 我们尝试用 shape 和 color 同时映射?
  geom_jitter(
    data = study_means,
    aes(y = Mean_Score, shape = Study.name), # 移除 color=Study.name 以免图例太乱，或者保留？用户提到legend缺失
    size = 3, alpha = 0.8, width = 0.2
  ) +

  # 统计检验 (基于样本数据)
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    vjust = 0.5,
    hide.ns = FALSE
  ) +
  scale_fill_npg() +
  scale_shape_manual(values = rep(1:25, times = 5)) + # 确保足够多的形状

  labs(
    title = "Microbial Dysbiosis Score across BMI & Disease",
    subtitle = "Boxplot: Sample Distribution | Points: Study Means",
    y = "Dysbiosis Score (Mean Dist to Healthy Normal)",
    x = NULL,
    shape = "Study"
  ) + # 明确 Legend 标题

  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right", # 强制显示 Legend
    legend.box = "vertical"
  )

print(p_box)

# 保存图片 (PDF & PNG)
ggsave(file.path(output_dir, "boxplot_dysbiosis_score_study_means.pdf"), p_box, width = 10, height = 6) # 稍微加宽以容纳 Legend
ggsave(file.path(output_dir, "boxplot_dysbiosis_score_study_means.png"), p_box, width = 10, height = 6, dpi = 300)

# ======================== 7. 统计验证: 多因素 PERMANOVA ========================
message("\n>>> 正在进行多因素 PERMANOVA 统计验证...")

# 7.1 数据清洗与对齐
# 这一步至关重要：多因素模型要求所有变量都不能有 NA
target_vars <- c("Study.name", "Age", "Sex", "BMI_Group1", "Disease")

# 找出所有目标变量都齐全的样本 (Complete Cases)
metadata_final <- metadata_final %>%
  filter(!Study.name %in% c("This_study_cohort5__NSHII")) # 先剔除未知研究)
valid_idx <- metadata_final[complete.cases(metadata_final[, target_vars]), ]
valid_idx <- valid_idx$Sample.Source.ID.processed
n_original <- nrow(metadata_final)
n_final <- length(valid_idx)

if (n_final < 10) {
  stop("错误：剔除缺失值后样本量太少 (<10)，无法分析！")
}

message(sprintf("   原始样本数: %d", n_original))
message(sprintf(
  "   清洗后样本数: %d (剔除了 %d 个含有缺失值的样本)",
  n_final, n_original - n_final
))

# 提取清洗后的 Metadata
meta_clean <- metadata_final[valid_idx, ]

# 提取清洗后的 距离矩阵
# 逻辑：Dist对象 -> 转矩阵 -> 按索引切片 -> 转回 Dist对象
dist_mat_full <- as.matrix(bray_dist)[valid_idx, valid_idx]
dist_clean <- as.dist(dist_mat_full[valid_idx, valid_idx])

# 准备批次校正的分组因子 (Strata)
blocks_strata <- as.factor(meta_clean$Study.name)

# 7.2 运行模型 1：加法模型 (Additive Model)
# 目的：获取 Age, Sex, BMI, Disease 的独立主效应 P 值
message("\n   正在运行 Model 1: 主效应模型 (Age + Sex + BMI + Disease)...")

perm_main <- adonis2(
  formula = dist_clean ~ Age + Sex + BMI_Group1 + Disease,
  data = meta_clean,
  permutations = 999,
  strata = blocks_strata, # 限制在 Study 内部置换(控制批次)
  by = "margin", # 关键：计算独立效应 (Type III SS)
  parallel = n_cores_perm
)

print("--- Model 1 结果 (Main Effects) ---")
print(perm_main)

# 7.3 运行模型 2：全模型 (Full/Interaction Model)
# 目的：专门检查 BMI * Disease 的交互作用是否显著
message("\n   正在运行 Model 2: 交互模型 (Age + Sex + BMI * Disease)...")

perm_full <- adonis2(
  formula = dist_clean ~ Age + Sex + BMI_Group1 * Disease,
  data = meta_clean,
  permutations = 999,
  strata = blocks_strata,
  by = "margin", # 关键
  parallel = n_cores_perm
)

print("--- Model 2 结果 (Interaction Term) ---")
print(perm_full)

save(perm_main, perm_full, file = file.path(output_dir, "permanova_results_remove_one_study.RData"))

# 7.4 简单的结论提示 & 结果保存
p_bmi <- perm_main["BMI_Group1", "Pr(>F)"]
p_disease <- perm_main["Disease", "Pr(>F)"]
p_inter <- perm_full["BMI_Group1:Disease", "R2"]
r2_inter <- perm_full["BMI_Group1:Disease", "R2"]

# 将结果输出到文件
stats_out_file <- file.path(output_dir, "permanova_stats_report.txt")
sink(stats_out_file)
cat("================ PERMANOVA Analysis Report ================\n")
cat("Date:", as.character(Sys.time()), "\n")
cat("Distance Metric:", dist_method, "(Sqrt Transformed)\n")
cat("Samples included:", n_final, "\n\n")

cat(">>> Model 1: Additive (Main Effects) Results\n")
print(perm_main)
cat("\n")

cat(">>> Model 2: Interaction (Full Model) Results\n")
print(perm_full)
cat("\n")

cat(">>> Summary Interpretation\n")
cat(sprintf("1. BMI Main Effect: P = %.4f\n", p_bmi))
cat(sprintf("2. Disease Main Effect: P = %.4f\n", p_disease))
cat(sprintf("3. Interaction Effect: P = %.4f (R2 = %.5f)\n", p_inter, r2_inter))

if (p_inter > 0.05) {
  cat("\nConclusion: Interaction is NOT significant. Supports additive effect model (Risk Seeding).\n")
} else {
  cat("\nConclusion: Interaction IS significant. Effect of BMI depends on disease status.\n")
}
sink()

message("\n>>> 统计分析完成！详细报告已保存至: ", stats_out_file)
message(">>> 全部流程结束。")
