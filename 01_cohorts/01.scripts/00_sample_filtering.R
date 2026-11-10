# ==============================================================================
# 脚本名称: metaphlan_processing.R
# 功能: 读取Metaphlan4数据，合并样本，并进行严格的质量控制和筛选
# ==============================================================================

# 1. 加载必要的库
library(tidyverse)
library(data.table) # 用于高效处理大文件
library(fs)         # 用于文件路径操作

# ==============================================================================
# 2. 定义路径 (根据用户提供的信息)
# 注意：在R中路径通常使用正斜杠 "/" 或双反斜杠 "\\"
# ==============================================================================

metadata_path <- "D:/Tianhui/proj_ob_AA_CRC/00. rawdata/00_mpa4_data/all_studies_metaphlan_summary.csv"
data_dir      <- "D:/Tianhui/proj_ob_AA_CRC/00. rawdata/00_mpa4_data"

# ==============================================================================
# 3. 读取元数据并根据 Total_Reads 筛选
# ==============================================================================

message("正在读取元数据并筛选 Reads > 100,000 的样本...")

# 读取 CSV
# 注意：check.names = FALSE 防止 R 自动修改列名
meta_df <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)

# 假设 CSV 中包含样本名称的列叫 "Sample_ID" 或类似的名称。
# 如果你的列名不同，请在此处修改，例如: valid_samples <- meta_df$Run_accession ...
# 这里我们假设第一列是样本名，或者有一列叫 'sample_id'
# 请根据实际 CSV 文件头修改下面的 'Sample_ID'

# 检查 Reads 数列名，通常是 "Total_Reads" 或 "n_reads"
if (!"Total_Reads" %in% names(meta_df)) {
  stop("错误: 在 CSV 中找不到 'Total_Reads' 列，请检查列名。")
}

# 筛选 Reads > 100000 的样本
filtered_meta <- meta_df %>%
  filter(Total_Reads > 1000000)
filtered_meta <- filtered_meta[!duplicated(filtered_meta$Sample_ID), ]  # 去重
# 获取合格的样本ID列表 (用于后续匹配)
# 假设样本ID在名为 "Sample" 或第一列
# 这里假设样本ID就是文件名(不含扩展名)，或者在CSV中有对应列
valid_sample_ids <- filtered_meta$Sample_ID # 默认取第一列作为样本ID，请根据实际情况调整

message(paste("初步筛选完成，剩余合格样本数:", length(valid_sample_ids)))


# ==============================================================================
# 4. 批量读取 TXT 文件并合并
# ==============================================================================

message("正在读取并合并所有 TXT 文件...")

# 获取目录下所有 txt 文件
files <- dir_ls(data_dir, glob = "*.txt")

# 定义读取单个文件的函数
read_mpa_file <- function(file_path) {
  # data.table::fread 非常智能，skip = "#" 会自动跳过以#开头的行
  # 并且它读取速度比 read.table 快得多
  dt <- fread(file_path, skip = "#", header = TRUE)
  
  # Metaphlan 输出通常是: clade_name  |  abundance_value
  # 列名通常是 clade_name 和 样本名
  # 我们需要确保第二列的列名是纯净的样本名
  
  # 获取当前文件的样本名 (假设文件名即样本名，或取文件的第二列列名)
  # 这里我们直接使用文件内容中的列名
  return(dt)
}

# 读取所有文件列表 (使用 map 读取为 list)
list_of_dfs <- map(files, read_mpa_file)

# 合并所有数据框
# 使用 reduce + full_join 保证所有物种都被保留
merged_df <- list_of_dfs %>%
  reduce(function(x, y) full_join(x, y, by = "clade_name"))

# 将所有 NA 替换为 0 (因为 full_join 会产生 NA)
merged_df[is.na(merged_df)] <- 0

colnames(merged_df)[-1] <- gsub("d_metagenome","", colnames(merged_df)[-1]) # 清理列名中的多余部分

message(paste("合并完成。包含物种数:", nrow(merged_df), "包含列数:", ncol(merged_df)))


# ==============================================================================
# 5. 数据清洗逻辑
# ==============================================================================

# 5.1 对齐样本：只保留在 Step 3 中 Total_Reads 达标的样本
# -----------------------------------------------------------
# 获取当前合并表中的样本列名 (排除 clade_name)
current_samples <- setdiff(names(merged_df), "clade_name")

# 找出既在合并表中，又在 valid_sample_ids 里的样本
# 注意：你需要确保 CSV 里的 ID 和 TXT 文件里的列名是一致的
# 如果文件名和列名不一致，这里可能需要额外的字符串处理
cols_to_keep <- intersect(current_samples, valid_sample_ids)
length(cols_to_keep)

# 如果交集为空，可能是命名格式问题，暂时保留所有列并发出警告
if(length(cols_to_keep) == 0) {
  warning("警告: CSV中的样本ID与TXT文件内的列名不匹配。将跳过基于CSV的列筛选，使用所有读取的样本。")
  cols_to_keep <- current_samples
}

# 更新 merged_df，只保留合格样本
final_df <- merged_df %>%
  select(clade_name, all_of(cols_to_keep))

message(paste("基于Total_Reads筛选后，剩余样本数:", length(cols_to_keep)))


# 5.2 处理 UNCLASSIFIED (筛选 > 50% 的样本)
# -----------------------------------------------------------

# 提取 UNCLASSIFIED 行
unclassified_row <- final_df %>% 
  filter(grepl("UNCLASSIFIED", clade_name, ignore.case = TRUE))

if(nrow(unclassified_row) > 0) {
  # 转置以便更容易处理
  unclass_long <- unclassified_row %>%
    pivot_longer(cols = -clade_name, names_to = "Sample", values_to = "Abundance")
  
  # 找出 Unclassified > 50 的样本 (Metaphlan 默认输出是 0-100)
  bad_samples_unclassified <- unclass_long %>%
    filter(Abundance > 50) %>%
    pull(Sample)
  
  message(paste("发现 Unclassified > 50% 的样本数:", length(bad_samples_unclassified)))
  if(length(bad_samples_unclassified) > 0) {
    print(bad_samples_unclassified)
  }
} else {
  bad_samples_unclassified <- c()
  message("未找到 UNCLASSIFIED 行，跳过此步骤。")
}


# 5.3 去除 UNCLASSIFIED 行，并检查单一优势菌 (> 50%)
# -----------------------------------------------------------

# 先从表中移除 UNCLASSIFIED
df_clean_step1 <- final_df %>%
  filter(!grepl("UNCLASSIFIED", clade_name, ignore.case = TRUE))

# 还要移除之前标记的 Unclassified 过高的样本
samples_step1 <- setdiff(names(df_clean_step1)[-1], bad_samples_unclassified)

df_clean_step2 <- df_clean_step1 %>%
  select(clade_name, all_of(samples_step1))

df_clean_step2 <- df_clean_step2[grep("t__", df_clean_step2$clade_name), ]  # 只保留属级及以下的分类单元
# 检查剩余物种中，是否有某一个菌的丰度超过 50%
# 计算每个样本在去除 Unclassified 后的最大丰度值
# 注意：这里我们直接看原始丰度值。如果 Unclassified 很高但被移除了，
# 剩余菌的相对比例并未重新归一化到100。
# 这里的逻辑是：如果原始数据中（排除Unclassified后），某个菌本身占比就 > 50

# 获取每个样本列的最大值
max_values <- df_clean_step2 %>%
  summarise(across(where(is.numeric), max, na.rm = TRUE)) %>%
  pivot_longer(everything(), names_to = "Sample", values_to = "Max_Abundance")

# 找出单一菌株 > 50 的样本
bad_samples_dominant <- max_values %>%
  filter(Max_Abundance > 50) %>%
  pull(Sample)

message(paste("发现单一菌株丰度 > 50% 的样本数:", length(bad_samples_dominant)))
if(length(bad_samples_dominant) > 0) {
  print(bad_samples_dominant)
}


# ==============================================================================
# 6. 生成最终表格并保存
# ==============================================================================

# 最终需要保留的样本 = (Total_Reads合格) - (Unclass过高) - (单一菌过高)
final_samples_to_keep <- setdiff(samples_step1, bad_samples_dominant)

final_output_table <- final_df %>%
  select(clade_name, all_of(final_samples_to_keep))

message("========================================================")
message(paste("原始样本数 (CSV):", nrow(meta_df)))
message(paste("最终保留样本数:", length(final_samples_to_keep)))
message("========================================================")

# 保存结果
output_dir <- "D:/Tianhui/proj_ob_AA_CRC/02.processed_data/01_metaphlan_data_filter/"
dir.create(output_dir, recursive = TRUE)
output_file <- file.path(output_dir, "metaphlan4_merged_filtered.tsv")
write.table(final_output_table, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

message(paste("结果已保存至:", output_file))
message("处理完成！")
