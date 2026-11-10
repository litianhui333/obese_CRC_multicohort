##### Beta Diversity Dispersion Analysis (9 Groups) #####

library(tidyverse)
library(vegan)

# ==============================================================================
# 1. 配置路径与参数
# ==============================================================================
base_dir <- "G:/My Drive/proj_ob_AA_CRC"
input_dir <- file.path(base_dir, "02.processed_data/02_aligned_data")
sp_dir <- file.path(base_dir, "02.processed_data/03_metaphlan_split_data")
output_dir <- file.path(base_dir, "03.results/02_beta_diversity")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

prev_threshold <- 0.1
trans_method <- "sqrt"
dist_method <- "bray"

# ======================== 2. 数据加载 ========================
meta_file <- sort(list.files(input_dir, pattern = "metadata_aligned.*\\.rds", full.names = TRUE), decreasing = TRUE)[1]
sp_file <- sort(list.files(sp_dir, pattern = "metaphlan_species_renormalized.*\\.rds", full.names = TRUE), decreasing = TRUE)[1]

if (is.na(meta_file) || is.na(sp_file)) stop("错误: 未找到输入数据文件。")

message(">>> 加载数据...")
metadata <- readRDS(meta_file)
metaphlan_raw <- readRDS(sp_file)

rownames(metadata) <- metadata$Sample.Source.ID.processed
metaphlan <- as.data.frame(t(metaphlan_raw))

# ======================== 3. 数据过滤与转换 ========================
message(">>> 数据预处理...")
keep_taxa <- colSums(metaphlan > 0) > (nrow(metaphlan) * prev_threshold)
metaphlan_filt <- metaphlan[, keep_taxa]
metaphlan_tss <- metaphlan_filt / rowSums(metaphlan_filt)
metaphlan_final <- sqrt(metaphlan_tss)

common_samples <- intersect(rownames(metaphlan_final), rownames(metadata))
metaphlan_final <- metaphlan_final[common_samples, ]
metadata_final <- metadata[common_samples, ]

# 核心：生成完整的 9 组标记 (Group_Label)
metadata_final <- metadata_final %>%
  mutate(Group_Label = paste(Disease, BMI_Group1, sep = "_")) %>%
  mutate(Group_Label = factor(Group_Label, levels = c(
    "CTRL_Normal", "CTRL_Ow", "CTRL_Ob",
    "AA_Normal",   "AA_Ow",   "AA_Ob",
    "CRC_Normal",  "CRC_Ow",  "CRC_Ob"
  ))) %>%
  # 滤除带有 NA 分组的非法样本
  filter(!is.na(Group_Label))

# 同步样本
common_samples_clean <- rownames(metadata_final)
metaphlan_final <- metaphlan_final[common_samples_clean, ]

# ======================== 4. 提取/计算距离矩阵 ========================
dist_file_name <- sprintf("dist_matrix_%s_%s_prev%s.rds", dist_method, trans_method, prev_threshold)
dist_file_path <- file.path(output_dir, dist_file_name)

if (file.exists(dist_file_path)) {
  message(">>> 读取已存在的距离矩阵...")
  bray_dist <- readRDS(dist_file_path)
  
  # 对齐到清洗过的样本
  bray_mat <- as.matrix(bray_dist)[common_samples_clean, common_samples_clean]
  bray_dist <- as.dist(bray_mat)
} else {
  message(">>> 计算全新的距离矩阵...")
  bray_dist <- vegdist(metaphlan_final, method = dist_method)
  saveRDS(bray_dist, dist_file_path)
}

# ==============================================================================
# 5. 执行 Betadisper (全 9 组的离散度计算)
# ==============================================================================
message("\n=== 开始计算 9 分组的 Beta Dispersion ===")

# 将所有样本和其分组传入 betadisper
mod_disp <- betadisper(bray_dist, metadata_final$Group_Label)

# 1. 提取每个样本的离散度（到对应分组空间质心的距离）
distances_df <- data.frame(
  SampleID = names(mod_disp$distances),
  Group_Label = metadata_final$Group_Label,
  Disease = metadata_final$Disease,
  BMI_Group = metadata_final$BMI_Group1,
  Distance_to_Centroid = mod_disp$distances
)

# 保存单个样本的离散度用于后续画 Boxplot
out_dist_file <- file.path(output_dir, "beta_dispersion_9groups_distances.csv")
write.csv(distances_df, out_dist_file, row.names = FALSE)
message(">> 样本中心距离已保存: ", out_dist_file)

# 2. 全局显著性检验与成对比较 (Permutation Test)
message(">> 执行 Permutation 统计学检验 (计算中，请稍候)...")
set.seed(123)
perm_res <- permutest(mod_disp, pairwise = TRUE, permutations = 999)

# 提取两两比较的 P 值
pairwise_pvals <- perm_res$pairwise$permuted
pairwise_df <- as.data.frame(pairwise_pvals)
pairwise_df$Comparison <- rownames(pairwise_df)
colnames(pairwise_df)[1] <- "P_Value_Permuted"

pairwise_df <- pairwise_df %>%
  select(Comparison, P_Value_Permuted) %>%
  mutate(FDR_Adjusted = p.adjust(P_Value_Permuted, method = "BH"))

out_pair_file <- file.path(output_dir, "beta_dispersion_9groups_pairwise_pvalues.csv")
write.csv(pairwise_df, out_pair_file, row.names = FALSE)
message(">> 组间成对检验结果已保存: ", out_pair_file)

# ==============================================================================
# 6. 可视化 (可选)
# ==============================================================================
message("\n>>> 生成 Boxplot 以预览 9 组离散度分布...")
library(ggplot2)

p_box <- ggplot(distances_df, aes(x = Group_Label, y = Distance_to_Centroid, fill = Disease)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 0.5) +
  theme_bw() +
  scale_fill_manual(values = c("CTRL" = "#00A08A", "AA" = "#F2A6A2", "CRC" = "#922687")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Beta Diversity Dispersion (Distance to Centroid)",
       y = "Distance to Centroid", x = "Subgroup")

ggsave(file.path(output_dir, "beta_dispersion_9groups_boxplot.pdf"), p_box, width = 8, height = 6)
message(">> 可视化完成。分析结束！")
