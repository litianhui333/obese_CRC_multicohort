# ==============================================================================
# 脚本名称: id_mapping_and_merge.R
# 功能: 匹配 Metaphlan 丰度表的列名与 Metadata 的样本ID，并进行对齐
# ==============================================================================

library(tidyverse)
library(data.table)
library(fs)

# ==============================================================================
# 1. 定义路径与加载数据
# ==============================================================================

# 定义输入路径
metaphlan_path <- "G:/My Drive/proj_ob_AA_CRC/02.processed_data/01_metaphlan_data_filter/metaphlan4_merged_filtered.tsv"
metadata_path <- "G:/My Drive/proj_ob_AA_CRC/02.processed_data/01_metadata_processed/metadata_processed_20260208_restored.rds"
sra_dir <- "G:/My Drive/proj_ob_AA_CRC/00.rawdata/03_SRA_data/"
output_dir <- "G:/My Drive/proj_ob_AA_CRC/02.processed_data/02_aligned_data"

# 创建输出目录
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

message(">>> 1. 正在加载数据...")

# 1.1 读取 Metaphlan 数据
# check.names = FALSE 极其重要，防止 "SRR123456" 变成 "X.SRR123456"
tax_df <- fread(metaphlan_path, header = TRUE)
# 转换为 tibble 方便操作，且将 clade_name 设为第一列
tax_df <- as_tibble(tax_df)
colnames(tax_df)[-1] <- gsub("d_metagenome", "", colnames(tax_df)[-1])
# 1.2 读取 Metadata
meta_df <- readRDS(metadata_path)

# 检查 Metadata 关键列
target_id_col <- "Sample.Source.ID.processed"
if (!target_id_col %in% names(meta_df)) {
  stop(paste0("错误: Metadata 中找不到目标列 '", target_id_col, "'"))
}

message(paste0("Metaphlan 数据加载完成: ", ncol(tax_df) - 1, " 个样本 (列)"))
message(paste0("Metadata 数据加载完成: ", nrow(meta_df), " 个样本 (行)"))


# ==============================================================================
# 2. 构建 SRA 映射桥接表 (The Bridge)
# ==============================================================================
message("\n>>> 2. 读取 SRA RunInfo 文件并构建映射关系...")

# 获取所有 SRA 信息文件
sra_files <- dir_ls(sra_dir, glob = "*.csv")

# 定义读取函数 (利用 fread 自动识别分隔符)
read_sra_info <- function(file_path) {
  dt <- fread(file_path, select = c("Run", "BioProject", "Sample_Name"), fill = TRUE)
  # 确保列存在，防止某些文件缺失列报错
  if (!all(c("Run", "BioProject", "Sample_Name") %in% names(dt))) {
    return(NULL)
  }
  return(dt)
}

# 批量读取并合并
sra_map_raw <- map_dfr(sra_files, read_sra_info) %>%
  distinct(Run, .keep_all = TRUE) # 去重

# ------------------------------------------------------------------------------
# [关键] 应用特殊命名规则
# ------------------------------------------------------------------------------
# 目标：构造一个 'Run_processed' 列，模拟 Metaphlan 列名可能对应的 SampleID 格式
# 或者构造一个逻辑，将 Metaphlan 的 Run ID 映射到 Metadata 的 Sample Name
# ------------------------------------------------------------------------------

sra_bridge <- sra_map_raw %>%
  mutate(
    # 这里应用你之前定义的规则：将 Run ID 转换为可能在 Metadata 中出现的形式
    # 注意：你的规则似乎是把 Run ID 转换成某种格式，去匹配 Metadata？
    # 或者 Metadata 里存的是 Sample_Name？
    # 下面的逻辑是：Metadata 里的 ID (Sample.Source.ID.processed) 可能等于这些变换后的 Run
    Mapped_ID = case_when(
      BioProject == "PRJEB27928" ~ paste0(Run, "-21-0"),
      BioProject == "PRJEB12449" ~ paste0(Run, "-27-0-0"),
      BioProject == "PRJNA447983" ~ paste0("CRC_MR_", Run, "_17"),
      BioProject == "PRJEB7774" ~ paste0("SID", Run),
      TRUE ~ Run # 默认保持 Run ID
    ),

    # 有时候 Metadata 里直接用的是 Sample_Name，这里做一个备用列
    # 优先用 Sample_Name 匹配，还是用上面的规则匹配？
    # 通常策略：先尝试用 Run 匹配，再尝试用 Sample_Name 匹配
    Raw_Run = Run
  )

message(paste0("SRA 映射表构建完成，涵盖 ", nrow(sra_bridge), " 个 Run ID。"))


# ==============================================================================
# 3. 执行 ID 匹配与替换
# ==============================================================================
message("\n>>> 3. 开始执行 ID 匹配...")

# 获取 Metaphlan 当前的所有列名 (排除 clade_name)
current_tax_cols <- names(tax_df)[-1] # 假设第一列是 clade_name
target_meta_ids <- meta_df[[target_id_col]]

# --- 匹配策略 ---

# 准备一个用来重命名的向量: named vector [Old_Name = New_Name]
# 初始状态：假设不需要改名
rename_map <- setNames(current_tax_cols, current_tax_cols)

# 3.1 情况 A: 直接匹配 (Perfect Match)
# Metaphlan 列名 直接就是 Metadata 里的 ID
direct_matches <- intersect(current_tax_cols, target_meta_ids)

# 3.2 情况 B: 通过 SRA Bridge 匹配
# 剩下的未匹配列
remaining_cols <- setdiff(current_tax_cols, direct_matches)

# 在 Bridge 表中查找这些 remaining_cols (作为 Run)
# 看看它们对应的 'Sample_Name' 或 'Mapped_ID' 是否在 Metadata 中
matched_via_bridge <- sra_bridge %>%
  filter(Run %in% remaining_cols) %>%
  # 优先检查 Sample_Name 是否在 Metadata 中
  mutate(Target_Match = case_when(
    Sample_Name %in% target_meta_ids ~ Sample_Name,
    Mapped_ID %in% target_meta_ids ~ Mapped_ID,
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Target_Match)) %>%
  select(Run, Target_Match)

# 更新重命名映射表
# 格式: rename_map[Old_Column_Name] <- New_Target_Name
if (nrow(matched_via_bridge) > 0) {
  for (i in 1:nrow(matched_via_bridge)) {
    old_col <- matched_via_bridge$Run[i]
    new_col <- matched_via_bridge$Target_Match[i]
    rename_map[old_col] <- new_col
  }
}

# 3.3 应用重命名
# 使用 plyr::rename 或者 dplyr::rename_with (但在向量操作上 data.table 或 base R 更直观)
# 我们直接修改 names(tax_df)
# 只需要修改那些在该映射表中的列

# 找到所有需要修改的列索引
cols_to_rename_indices <- match(names(rename_map), names(tax_df))
# 排除 NA (clade_name 或其他没找到的)
valid_indices <- cols_to_rename_indices[!is.na(cols_to_rename_indices)]

# 执行重命名
# 注意：rename_map 的 names 是旧名字 (Run)，values 是新名字 (Metadata ID)
names(tax_df)[valid_indices] <- rename_map[names(tax_df)[valid_indices]]

message("ID 重命名完成。")

# ==============================================================================
# 4. 最终对齐与筛选 (Intersection)
# ==============================================================================
message("\n>>> 4. 最终数据对齐...")

# 获取重命名后的 Metaphlan 样本名
new_tax_samples <- names(tax_df)[-1]

# 找出交集 (Final Intersection)
common_samples <- intersect(new_tax_samples, target_meta_ids)

# --- 报告 ---
message("--------------------------------------------------")
message(paste0("Metadata 样本数: ", length(target_meta_ids)))
message(paste0("Metaphlan 原始样本数: ", length(current_tax_cols)))
message(paste0("最终成功匹配样本数: ", length(common_samples)))
message("--------------------------------------------------")

# 检查哪些丢失了
lost_meta <- setdiff(target_meta_ids, common_samples)
lost_tax <- setdiff(new_tax_samples, common_samples)

if (length(lost_meta) > 0) {
  message(paste0("警告: 有 ", length(lost_meta), " 个 Metadata 样本在 Metaphlan 中找不到数据。"))
  # print(head(lost_meta))
}

# --- 执行最终取子集 ---

# 1. 过滤 Metadata
final_meta <- meta_df %>%
  filter(.data[[target_id_col]] %in% common_samples) %>%
  arrange(match(.data[[target_id_col]], common_samples)) # 按 common_samples 排序

# 2. 过滤 Metaphlan
# select(clade_name, all_of(common_samples))
final_tax <- tax_df %>%
  select(clade_name, all_of(common_samples))

# 确保两者的样本顺序完全一致 (这对后续分析至关重要)
# Metadata 的行顺序必须等于 Metaphlan 的列顺序(除了clade_name)
if (!all(final_meta[[target_id_col]] == names(final_tax)[-1])) {
  stop("严重错误: 数据对齐后顺序不一致，请检查代码逻辑！")
} else {
  message("检查通过: Metadata 行顺序与 Metaphlan 列顺序完全一致。")
}

# ==============================================================================
# 5. 保存结果
# ==============================================================================
message("\n>>> 5. 保存对齐后的数据...")

# 定义文件名
date_str <- format(Sys.Date(), "%Y%m%d")

# 保存 Metadata
saveRDS(final_meta, file.path(output_dir, paste0("metadata_aligned_", date_str, ".rds")))
write.csv(final_meta, file.path(output_dir, paste0("metadata_aligned_", date_str, ".csv")), row.names = FALSE)

# 保存 Metaphlan 表
# 建议保存为 RDS (R专用，快，保留类型) 和 TSV (通用)
saveRDS(final_tax, file.path(output_dir, paste0("metaphlan_aligned_", date_str, ".rds")))
write.table(final_tax, file.path(output_dir, paste0("metaphlan_aligned_", date_str, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message(paste0("所有文件已保存至: ", output_dir))
message("流程结束。")
