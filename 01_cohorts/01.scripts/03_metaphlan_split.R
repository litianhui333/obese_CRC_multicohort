# ==============================================================================
# 脚本名称: 03_metaphlan_split.R
# 功能: 将 Metaphlan 丰度表拆分为 Species (物种) 和 Strain (菌株) 层级
# 输入: 02_aligned_data 中的 metaphlan_aligned
# 输出: 03_metaphlan_split_data (Species表, Strain表, Taxonomy信息)
# ==============================================================================

library(tidyverse)
library(data.table)
library(fs)

# ======================== 1. 路径设置 ========================
input_dir  <- "D:/Tianhui/proj_ob_AA_CRC/02.processed_data/02_aligned_data"
output_dir <- "D:/Tianhui/proj_ob_AA_CRC/02.processed_data/03_metaphlan_split_data"

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

timestamp <- format(Sys.Date(), "%Y%m%d")

# ======================== 2. 数据读取 ========================
# 读取对齐后的 Metaphlan 数据
tax_files <- dir_ls(input_dir, glob = "*metaphlan_aligned*.rds")
if(length(tax_files) == 0) stop("未找到对齐后的 Metaphlan RDS 文件！")
latest_tax_file <- sort(tax_files, decreasing = TRUE)[1]

message(paste("正在读取 Metaphlan 数据:", latest_tax_file))
metaphlan_data <- readRDS(latest_tax_file)

# 确保 clade_name 存在
if(!"clade_name" %in% names(metaphlan_data)) {
  # 如果行名是 clade_name，转换为列
  metaphlan_data <- metaphlan_data %>% rownames_to_column("clade_name")
}

message(paste0("数据加载完成。维度: ", nrow(metaphlan_data), " 行 x ", ncol(metaphlan_data), " 列"))


# ======================== 3. 层级拆分逻辑 ========================
message("\n=== 开始拆分 Species 和 Strain 层级 ===")

# --- 3.1 提取 Strain (菌株) 层级 ---
# 特征: 包含 "t__" (terminal level for Metaphlan strains)
strain_df <- metaphlan_data %>% 
  filter(grepl("t__", clade_name))

message(paste("提取到 Strain (菌株) 层级行数:", nrow(strain_df)))

# --- 3.2 提取 Species (物种) 层级 ---
# 特征: 包含 "s__" 但不包含 "t__" (即种水平汇总)
species_df <- metaphlan_data %>% 
  filter(grepl("s__", clade_name) & !grepl("t__", clade_name))

message(paste("提取到 Species (物种) 层级行数:", nrow(species_df)))


# ======================== 4. 处理 Strain 的 Taxonomy ========================
message("\n=== 处理 Strain 分类信息 ===")

if (nrow(strain_df) > 0) {
  # 提取 clade_name 进行分割
  # 使用 data.table 的 tstrsplit 或 tidyr::separate 处理字符串会更快
  # 这里为了稳健性，使用 base split
  
  raw_clades <- strain_df$clade_name
  
  # 分割字符串
  split_list <- strsplit(raw_clades, "\\|")
  
  # 转换为数据框 (处理可能长度不一致的情况，虽然 Metaphlan 通常是一致的)
  # 预定义层级名称
  ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species", "strain")
  
  taxonomy_mat <- do.call(rbind, lapply(split_list, function(x) {
    # 补齐或截断到8个层级
    length(x) <- 8 
    return(x)
  }))
  
  taxonomy_df <- as.data.frame(taxonomy_mat)
  colnames(taxonomy_df) <- ranks
  
  # 创建简化的 ID: "Species|Strain"
  # 这种格式在后续分析（如 Panphlan 或 关联分析）中很常用
  taxonomy_df$clade_name_simple <- paste0(taxonomy_df$species, "|", taxonomy_df$strain)
  
  # 将完整 Taxonomy 保存备用
  # 并将 Strain 数据表的 clade_name 替换为简单格式? 
  # 这里保留原始 clade_name 作为行名通常更安全，或者使用简单名。
  # 你的原始代码使用了简单名。我们这里也生成简单名版本。
  
  strain_df_simple <- strain_df
  strain_df_simple$clade_name <- taxonomy_df$clade_name_simple
  
  message("Taxonomy 解析完成。")
} else {
  taxonomy_df <- NULL
  strain_df_simple <- NULL
  warning("警告: 未检测到 Strain 数据。")
}


# ======================== 5. 格式整理与保存 ========================
message("\n=== 保存文件 ===")

save_data <- function(df, prefix, use_simple_name = FALSE) {
  if (is.null(df) || nrow(df) == 0) return()
  
  # 整理行名: 将 clade_name 设为行名，并移除该列
  # 注意：如果有重复的 clade_name (极少见)，这步会报错
  df_save <- df %>% column_to_rownames("clade_name")
  
  # 文件路径
  rds_path <- file.path(output_dir, paste0(prefix, "_", timestamp, ".rds"))
  csv_path <- file.path(output_dir, paste0(prefix, "_", timestamp, ".csv")) # 或 .tsv
  
  # 保存
  saveRDS(df_save, rds_path)
  # write.csv 对大文件较慢且包含行名处理，建议 write.table
  write.table(df_save, csv_path, sep = ",", quote = FALSE, col.names = NA)
  
  message(paste("已保存:", basename(rds_path)))
}

# 1. 保存 Species 数据
save_data(species_df, "metaphlan_species")

# 2. 保存 Strain 数据 (使用简化名版本)
if(!is.null(strain_df_simple)) {
  save_data(strain_df_simple, "metaphlan_strain")
}

# 3. 保存 Taxonomy 信息 (仅针对 Strain)
if (!is.null(taxonomy_df)) {
  tax_rds <- file.path(output_dir, paste0("strain_taxonomy_", timestamp, ".rds"))
  saveRDS(taxonomy_df, tax_rds)
  write.csv(taxonomy_df, file.path(output_dir, paste0("strain_taxonomy_", timestamp, ".csv")), row.names = FALSE)
  message(paste("已保存: Strain Taxonomy 信息"))
}

message(paste0("\n所有拆分数据已保存至: ", output_dir))
