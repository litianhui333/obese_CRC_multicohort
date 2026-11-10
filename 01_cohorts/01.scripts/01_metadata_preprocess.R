# ==============================================================================
# 脚本名称: 01_metadata_preprocess.R
# 功能: 合并、清洗和处理多中心研究的元数据 (Metadata)
# 原作者: Tianhui Li
# 整理日期: 2025-12-04 (Updated: 2026-01-15)
# ==============================================================================

# ======================== 1. 环境设置与依赖加载 ========================
library(tidyverse)
library(purrr)
library(stringr) # 明确加载 stringr 用于字符串处理

# 设置输入和输出路径
# 设置输入和输出路径
input_dir  <- "G:/My Drive/proj_ob_AA_CRC/00.rawdata/00_metadata_check/"
output_dir <- "G:/My Drive/proj_ob_AA_CRC/02.processed_data/01_metadata_processed/"

# 创建输出目录（如果不存在）
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 设置全局日期标记（所有输出文件使用相同时间戳）
timestamp <- format(Sys.time(), "%Y%m%d")


# ======================== 2. 函数定义 ========================

#' 读取单个元数据文件
#' @param path 文件路径
#' @return 数据框
read_metadata <- function(path) {
  # 根据扩展名选择读取方式
  if (grepl("\\.(csv)$", path, ignore.case = TRUE)) {
    df <- read.csv(path, stringsAsFactors = FALSE)
  } else {
    df <- read.delim(path, stringsAsFactors = FALSE)
  }
  
  # 确保 Age 列为 numeric 类型，防止字符型干扰
  if ("Age" %in% colnames(df)) {
    df$Age <- as.numeric(df$Age)
  }
  
  return(df)
}

#' 计算 BMI 分组
#' @param df 包含 BMI, Country 等列的数据框
#' @return 添加了 BMI 分组列的数据框
calculate_bmi_groups <- function(df) {
  df %>%
    mutate(
      # 定义洲际分组
      Continent = case_when(
        Country == "China" ~ "China",
        Country %in% c("India", "Japan") ~ "Other_Asia",
        TRUE ~ "International"
      ),
      
      # 分组策略 1: 基于地区差异的标准 (Original)
      BMI_Group1 = case_when(
        BMI < 18.5 ~ "Underweight",
        Continent == "China" & BMI < 24 ~ "Normal",
        Continent == "China" & BMI >= 24 & BMI < 28 ~ "Ow",
        Continent == "China" & BMI >= 28 ~ "Ob",
        Continent == "Other_Asia" & BMI < 23 ~ "Normal", 
        Continent == "Other_Asia" & BMI >= 23 & BMI < 25 ~ "Ow",
        Continent == "Other_Asia" & BMI >= 25 ~ "Ob",
        Continent == "International" & BMI < 25 ~ "Normal",
        Continent == "International" & BMI >= 25 & BMI < 30 ~ "Ow",
        Continent == "International" & BMI >= 30 ~ "Ob",
        TRUE ~ NA_character_
      ),
      
      # 分组策略 1 (变体): 中国使用 Other_Asia 标准
      BMI_Group1_China_as_OtherAsia = case_when(
        BMI < 18.5 ~ "Underweight",
        Country == "China" & BMI < 23 ~ "Normal",
        Country == "China" & BMI >= 23 & BMI < 25 ~ "Ow",
        Country == "China" & BMI >= 25 ~ "Ob",
        Continent == "Other_Asia" & BMI < 23 ~ "Normal", 
        Continent == "Other_Asia" & BMI >= 23 & BMI < 25 ~ "Ow",
        Continent == "Other_Asia" & BMI >= 25 ~ "Ob",
        Continent == "International" & BMI < 25 ~ "Normal",
        Continent == "International" & BMI >= 25 & BMI < 30 ~ "Ow",
        Continent == "International" & BMI >= 30 ~ "Ob",
        TRUE ~ NA_character_
      ),
      
      # 分组策略 2: 仅区分 Normal vs Obese (Original)
      BMI_Group2_ob = case_when(
        BMI < 18.5 ~ "Underweight",
        Continent == "China" & BMI < 28 ~ "Normal",
        Continent == "China" & BMI >= 28 ~ "Obese",
        Continent == "Other_Asia" & BMI < 25 ~ "Normal",
        Continent == "Other_Asia" & BMI >= 25 ~ "Obese",
        Continent == "International" & BMI < 30 ~ "Normal",
        Continent == "International" & BMI >= 30 ~ "Obese",
        TRUE ~ NA_character_
      ),
      
      # 分组策略 2 (变体): 中国使用 Other_Asia 标准
      BMI_Group2_ob_China_as_OtherAsia = case_when(
        BMI < 18.5 ~ "Underweight",
        Country == "China" & BMI < 25 ~ "Normal",
        Country == "China" & BMI >= 25 ~ "Obese",
        Continent == "Other_Asia" & BMI < 25 ~ "Normal",
        Continent == "Other_Asia" & BMI >= 25 ~ "Obese",
        Continent == "International" & BMI < 30 ~ "Normal",
        Continent == "International" & BMI >= 30 ~ "Obese",
        TRUE ~ NA_character_
      ),
      
      # 分组策略 3: 合并超重/肥胖 (Original)
      BMI_Group3_Ow = case_when(
        BMI < 18.5 ~ "Underweight",
        BMI_Group1 %in% c("Ow", "Ob") ~ "Overweight/Obese",
        BMI_Group1 == "Normal" ~ "Normal",
        TRUE ~ NA_character_
      ),
      
      # 分组策略 3 (变体): 中国使用 Other_Asia 标准
      BMI_Group3_Ow_China_as_OtherAsia = case_when(
        BMI < 18.5 ~ "Underweight",
        BMI_Group1_China_as_OtherAsia %in% c("Ow", "Ob") ~ "Overweight/Obese",
        BMI_Group1_China_as_OtherAsia == "Normal" ~ "Normal",
        TRUE ~ NA_character_
      )
    ) %>%
    select(-Continent) # 计算完成后移除临时变量 Continent
}


# ======================== 3. 数据读取与初步清洗 ========================
message("正在读取元数据文件...")

# 获取所有文件路径
metadata_paths <- list.files(input_dir, pattern = "\\.(tsv|csv)$", full.names = TRUE)

# 合并所有元数据文件
# .id = "source_file" 会自动添加一列显示来源文件名
combined_metadata <- map_dfr(metadata_paths, read_metadata, .id = "source_file")

message(paste("已合并文件数:", length(metadata_paths), "| 总行数:", nrow(combined_metadata)))

# [NEW] 导出原始 Tumor.Staging.AJCC 分布情况 (Before Cleaning)
if ("Tumor.Staging.AJCC" %in% colnames(combined_metadata)) {
    staging_tab <- table(combined_metadata$Tumor.Staging.AJCC, useNA = "ifany") %>% as.data.frame()
    colnames(staging_tab) <- c("Raw_Stage", "Count")
    write.csv(staging_tab, file.path(output_dir, "Tumor_Staging_Raw_Frequency.csv"), row.names = FALSE)
    print("Raw Tumor Staging frequency table saved.")
}

# 标准化基础字段
combined_metadata$Country[combined_metadata$Country == "USA"] <- "United States"
combined_metadata$Sex[combined_metadata$Sex == ""] <- NA

# 初步计算 BMI 分组 (在未过滤的数据上)
combined_metadata <- calculate_bmi_groups(combined_metadata)


# ======================== 4. 数据过滤 (Filtering) ========================
message("正在进行数据过滤...")

# 1. 重新计算 Continent 用于简单的 Asia/Europe 分类
combined_metadata <- combined_metadata %>%
  mutate(Continent = if_else(Country %in% c("China", "India", "Japan"), "Asia", "Europe")) 

# 2. 生成过滤后的数据集
filtered_metadata <- combined_metadata %>%
  filter(BMI >= 18.5) %>%        # 过滤掉体重不足样本
  calculate_bmi_groups() %>%     # 重新计算分组 (确保分组逻辑基于过滤后的数据)
  mutate(Continent = if_else(Country %in% c("China", "India", "Japan"), "Asia", "Europe")) 


# ======================== 5. 变量类型转换与标准化 ========================
message("正在转换变量类型...")

# 定义研究类型 (Inhouse vs Public)
inhouse_studies <- c("2025_HK", "2023_Yunnan", "2021_YuJ")

combined_metadata$Study_type <- ifelse(combined_metadata$Study.name %in% inhouse_studies, "inhouse", "Public")
filtered_metadata$Study_type <- ifelse(filtered_metadata$Study.name %in% inhouse_studies, "inhouse", "Public")

# 统一 Disease/Condition 命名
filtered_metadata$Study.condition[filtered_metadata$Study.condition == "colorectal carcinoma"] <- "CRC"
combined_metadata$Study.condition[combined_metadata$Study.condition == "colorectal carcinoma"] <- "CRC"

# 因子化和数值化处理 (针对 filtered_metadata)
filtered_metadata <- filtered_metadata %>%
  mutate(
    BMI_Group1 = factor(BMI_Group1, levels = c("Normal", "Ow", "Ob")),
    BMI_Group2_ob = factor(BMI_Group2_ob, levels = c("Normal", "Obese")),
    BMI_Group3_Ow = factor(BMI_Group3_Ow, levels = c("Normal", "Overweight/Obese")),
    
    BMI_Group1_China_as_OtherAsia = factor(BMI_Group1_China_as_OtherAsia, levels = c("Normal", "Ow", "Ob")),
    BMI_Group2_ob_China_as_OtherAsia = factor(BMI_Group2_ob_China_as_OtherAsia, levels = c("Normal", "Obese")),
    BMI_Group3_Ow_China_as_OtherAsia = factor(BMI_Group3_Ow_China_as_OtherAsia, levels = c("Normal", "Overweight/Obese")),
    
    Disease = factor(Disease, levels = c("Control", "Adenoma", "colorectal carcinoma"), labels = c("CTRL", "AA", "CRC")),
    Sex = factor(Sex, levels = c("Male", "Female")),
    Study_type = factor(Study_type, levels = c("inhouse", "Public")),
    Age = as.numeric(Age)
  )

# 临床信息清洗
# 1. 肿瘤位置
filtered_metadata$Primary.Tumor.Location <- ifelse(
  filtered_metadata$Primary.Tumor.Location %in% c("left_sided", "right_sided"),
  filtered_metadata$Primary.Tumor.Location,
  NA
)

# 2. 肿瘤分期 (AJCC) Cleaning
# 目标: 合并子分期 (如 IIA -> II) 并将非标准分期设为 NA
# Step 1: 移除分期末尾的字母 (A, B, C)
filtered_metadata$Tumor.Staging.AJCC <- gsub("[ABC]$", "", filtered_metadata$Tumor.Staging.AJCC)

# Step 2: 仅保留标准的 I, II, III, IV
valid_stages <- c("I", "II", "III", "IV")
filtered_metadata$Tumor.Staging.AJCC <- ifelse(
  filtered_metadata$Tumor.Staging.AJCC %in% valid_stages,
  filtered_metadata$Tumor.Staging.AJCC,
  NA
)

# ======================== 6. 数据保存 ========================
message("正在保存处理后的数据...")

# 构建输出文件名 (包含日期)
final_rds_path <- file.path(output_dir, paste0("metadata_processed_", timestamp, ".rds"))
final_csv_path <- file.path(output_dir, paste0("metadata_processed_", timestamp, ".csv"))

# 保存
saveRDS(filtered_metadata, final_rds_path)
write.csv(filtered_metadata, final_csv_path, row.names = FALSE)

message(paste("处理完成!"))
message(paste("RDS 文件:", final_rds_path))
message(paste("CSV 文件:", final_csv_path))
# 去掉分期后面的字母 (如 IIA -> II)
# Incomplete trailing fragment retained as a comment: filtered_metadata$Tumor.Staging.AJCC <- gsub("^(I{1,