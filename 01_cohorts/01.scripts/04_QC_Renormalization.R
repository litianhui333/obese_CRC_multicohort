# ==============================================================================
# 脚本名称: 04_QC_Renormalization.R
# 功能: 
# 1. 检查 Unclassified 比例 (QC)
# 2. 移除 Unclassified 并重新归一化 (Renormalize relative abundance)
# 3. 提取 Species 层级并重新计算 Alpha 多样性
# ==============================================================================

library(tidyverse)
library(vegan)
library(lmerTest)
library(ggpubr)

# 1. 路径设置
# ------------------------------------------------------------------------------
# 输入：对齐后的原始数据 (包含 UNCLASSIFIED)
input_dir  <- "D:/Tianhui/proj_ob_AA_CRC/02.processed_data/02_aligned_data"
# 输出：清洗后的数据 (用于下一步分析)
output_dir <- "D:/Tianhui/proj_ob_AA_CRC/02.processed_data/04_cleaned_data"

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
timestamp <- format(Sys.Date(), "%Y%m%d")

# 2. 读取数据
# ------------------------------------------------------------------------------
message("正在读取对齐后的原始数据...")

# A. 读取 Metadata
meta_file <- list.files(input_dir, pattern = "metadata_aligned.*\\.rds", full.names = TRUE)[1]
metadata <- readRDS(meta_file)
rownames(metadata) <- metadata$Sample.Source.ID.processed

# B. 读取 完整的 Metaphlan 表 (必须包含 UNCLASSIFIED)
tax_file <- list.files(input_dir, pattern = "metaphlan_aligned.*\\.rds", full.names = TRUE)[1]
full_tax_table <- readRDS(tax_file)

# 确保格式正确 (行=Taxa, 列=Sample)
if("clade_name" %in% names(full_tax_table)) {
  full_tax_table <- full_tax_table %>% column_to_rownames("clade_name")
}

# 确保样本对齐
common_ids <- intersect(rownames(metadata), colnames(full_tax_table))
full_tax_table <- full_tax_table[, common_ids]
metadata <- metadata[common_ids, ]

message(paste("分析样本数:", length(common_ids)))

# ==============================================================================
# Part 1: Unclassified QC 分析
# ==============================================================================
message("=== Step 1: Unclassified 比例 QC ===")

# 查找 UNCLASSIFIED 行
unclass_row <- grep("UNCLASSIFIED", rownames(full_tax_table))

if(length(unclass_row) > 0) {
  # 提取 Unclassified 丰度
  # 注意：Metaphlan 输出通常是百分比 (0-100)
  unclass_abun <- as.numeric(full_tax_table[unclass_row, ])
  metadata$Unclassified_Prop <- unclass_abun
  
  # 1.1 绘制 QC 箱线图
  # --------------------
  message("绘制 QC 箱线图...")
  
  plot_qc <- function(x_var, title) {
    ggplot(metadata, aes_string(x = x_var, y = "Unclassified_Prop", fill = x_var)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
      stat_compare_means(label = "p.format", method = "wilcox.test", ref.group = ".all.") +
      labs(title = title, y = "Unclassified %", x = NULL) +
      theme_classic() +
      theme(legend.position = "none")
  }
  
  p1 <- plot_qc("Disease", "Unclassified by Disease")
  p2 <- plot_qc("BMI_Group1", "Unclassified by BMI")
  
  pdf(file.path(output_dir, paste0("QC_Unclassified_Distribution_", timestamp, ".pdf")), width = 8, height = 4)
  print(ggarrange(p1, p2, ncol = 2))
  dev.off()
  
  # 1.2 统计检验 (LMM)
  # --------------------
  message("执行 LMM 检验...")
  # 检查 Unclassified 是否受生物学因素显著影响
  lmm_qc <- tryCatch({
    lmer(Unclassified_Prop ~ Disease + BMI_Group1 + Age + Sex + (1|Study.name), data = metadata)
  }, error = function(e) NULL)
  
  if(!is.null(lmm_qc)) {
    print(summary(lmm_qc))
    # 保存统计结果到文本文件
    capture.output(summary(lmm_qc), file = file.path(output_dir, "QC_Unclassified_LMM_Result.txt"))
  }
  } else {
  message("警告: 未检测到 UNCLASSIFIED 行。假设数据已经清洗过。")
  clean_tax_table <- full_tax_table
  }

# ==============================================================================
# Part 2: 读取拆分数据并重新归一化 (修改部分)
# ==============================================================================
message("\n=== Step 2: Species 和 Strain 分别重新归一化 ===")

# 定义归一化函数 (Sum-scaling to 100)
renormalize_to_100 <- function(df, level_name) {
  message(paste0("正在处理: ", level_name))
  message(paste0("  - 原始维度: ", nrow(df), " x ", ncol(df)))
  
  # 确保全是数值
  # (Metaphlan拆分后的表通常是纯数值矩阵，行名是clade)
  
  # 按列求和并归一化
  df_norm <- apply(df, 2, function(x) {
    s <- sum(x, na.rm = TRUE)
    if(s == 0) return(x) # 防止全0样本报错
    return((x / s) * 100)
  }) %>% as.data.frame()
  
  # 简单检查归一化结果 (检查前几个样本和是否为100)
  check_sums <- colSums(df_norm)
  message(paste0("  - 归一化后样本均值总和: ", round(mean(check_sums), 2), "% (应接近 100)"))
  
  return(df_norm)
}

# --- 处理 Species ---
# 1. 查找最新的 Species 文件
sp_files <- dir_ls(split_data_dir, glob = "*metaphlan_species*.rds")[1]
# 排除掉已经是 renormalized 的文件 (防止重复处理)
sp_files <- sp_files[!grepl("renormalized", sp_files)]

if(length(sp_files) > 0) {
  latest_sp_file <- sort(sp_files, decreasing = TRUE)[1]
  message(paste("读取 Species 文件:", basename(latest_sp_file)))
  
  species_df <- readRDS(latest_sp_file)
  
  # 2. 重新归一化
  species_norm <- renormalize_to_100(species_df, "Species")
  
  # 3. 保存回原目录 (split_data_dir)
  out_name_sp <- paste0("metaphlan_species_renormalized_", timestamp, ".rds")
  out_csv_sp <- paste0("metaphlan_species_renormalized_", timestamp, ".csv")
  
  saveRDS(species_norm, file.path(split_data_dir, out_name_sp))
  write.csv(species_norm, file.path(split_data_dir, out_csv_sp), quote = FALSE)
  message(paste("  -> 已保存:", out_name_sp))
  
} else {
  warning("未找到 Species RDS 文件！")
}

# --- 处理 Strain ---
# 1. 查找最新的 Strain 文件
st_files <- dir_ls(split_data_dir, glob = "*metaphlan_strain*.rds")
st_files <- st_files[!grepl("renormalized", st_files)] # 排除已归一化的

if(length(st_files) > 0) {
  latest_st_file <- sort(st_files, decreasing = TRUE)[1]
  message(paste("\n读取 Strain 文件:", basename(latest_st_file)))
  
  strain_df <- readRDS(latest_st_file)
  
  # 2. 重新归一化
  strain_norm <- renormalize_to_100(strain_df, "Strain")
  
  # 3. 保存回原目录 (split_data_dir)
  out_name_st <- paste0("metaphlan_strain_renormalized_", timestamp, ".rds")
  out_csv_st <- paste0("metaphlan_strain_renormalized_", timestamp, ".csv")
  
  saveRDS(strain_norm, file.path(split_data_dir, out_name_st))
  write.csv(strain_norm, file.path(split_data_dir, out_csv_st), quote = FALSE)
  message(paste("  -> 已保存:", out_name_st))
  
} else {
  warning("未找到 Strain RDS 文件！")
}

message("\n=== 所有归一化工作完成 ===")
message(paste("结果已保存在:", split_data_dir))

