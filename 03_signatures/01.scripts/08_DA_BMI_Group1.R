
# ==============================================================================
# 脚本名称: 08_DA_updated_for_new_data.R
# 功能: 
#   1. 多区域分析 (Overall, Asia, Europe)
#   2. 分层差异分析 (Stratified by BMI & Stratified by Disease)
#   3. 全局疾病分析 (Model A: Adj BMI vs Model B: No Adj BMI)
#   4. BMI 连续变量相关性 (Overall & Stratified)
# ==============================================================================


# 1. 路径设置
# ------------------------------------------------------------------------------
base_dir   <- "G:/My Drive/proj_ob_AA_CRC"
input_dir  <- file.path(base_dir, "02.processed_data/02_aligned_data")
sp_dir     <- file.path(base_dir, "02.processed_data/03_metaphlan_split_data")
output_dir <- file.path(base_dir, "03.results/04_differential_abundance")
script_dir <- file.path(base_dir, "01.scripts")

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
timestamp <- format(Sys.Date(), "%Y%m%d")

library(tidyverse)
library(MMUPHin)

# 加载您的核心函数
# 确保 DA_function.R 位于 script_dir 下
source(file.path(script_dir, "DA_function.R"))

# ======================== 加载必要包 ========================
library(tidyverse)
library(vegan)
library(permute)
library(stringr)

# 2. 数据读取
# ------------------------------------------------------------------------------
message("=== 读取数据 ===")
meta_file <- sort(list.files(input_dir, pattern = "metadata_aligned.*\\.rds", full.names = T), decreasing = T)[1]
sp_file   <- sort(list.files(sp_dir, pattern = "metaphlan_strain_renormalized.*\\.rds", full.names = T), decreasing = T)[1]

if (length(meta_file) == 0 || length(sp_file) == 0) {
  stop("未找到 metadata 或 metaphlan 数据文件，请检查路径。")
}

metadata <- readRDS(meta_file)
metaphlan <- readRDS(sp_file)
rownames(metadata) <- metadata$Sample.Source.ID.processed

# 确保 BMI 是数值型，并确保 BMI_Group1 和 Disease 是字符或因子
metadata$BMI <- as.numeric(as.character(metadata$BMI))

# === [关键修复] 全局指定 Reference Levels ===
# 这将一劳永逸地解决 "Internal Maaslin run error... more than 2 levels" 的问题
# 1. 设置 Disease: 参考组设为 CTRL
metadata$Disease <- factor(metadata$Disease, levels = c("CTRL", "AA", "CRC"))

# 2. 设置 BMI_Group1: 参考组设为 Normal
metadata$BMI_Group1 <- factor(metadata$BMI_Group1, levels = c("Normal", "Ow", "Ob"))

# (如果您跑的是 China_as_OtherAsia 脚本，也要设置那个变量)
if("BMI_Group1_China_as_OtherAsia" %in% names(metadata)){
  metadata$BMI_Group1_China_as_OtherAsia <- factor(metadata$BMI_Group1_China_as_OtherAsia, 
                                                   levels = c("Normal", "Ow", "Ob"))}


# ======================== Study_type分组分析脚本 ========================
# 创建过滤配置
config <- create_analysis_config(
  transform_method = "AST",
  analysis_method = "HS",
  BMI_Group_var = "BMI_Group1", 
  covariates = c("Age", "Sex"),
  study_prevalence_cut = 0.1,
  abundance_threshold = 0.0001,
  min_studies_ratio = 0.7,
  pval_type = "both",
  sig_criterion = "or",
  pval_threshold = 0.02,
  qval_threshold = 0.2,
  is_inhouse = FALSE,
  force_keep = c(
    "s__Fusobacterium_nucleatum", 
    "s__Peptostreptococcus_stomatis", 
    "s__Parvimonas_micra",
    "s__Solobacterium_moorei"
  )
)

# ======================== 按Continent拆分数据 ========================

# 提取Asia和Europe的样本
Asia_samples <- metadata$Sample.Source.ID.processed[metadata$Continent == "Asia" & !is.na(metadata$Continent)]
Europe_samples <- metadata$Sample.Source.ID.processed[metadata$Continent == "Europe" & !is.na(metadata$Continent)]

cat("Asia样本数:", length(Asia_samples), "\n")
cat("Europe样本数:", length(Europe_samples), "\n")


# 创建Asia数据子集
metadata_Asia <- metadata[metadata$Sample.Source.ID.processed %in% Asia_samples, ]
metaphlan_Asia <- as.data.frame(t(metaphlan[,colnames(metaphlan) %in% Asia_samples, ]))

# 创建Europe数据子集
metadata_Europe <- metadata[metadata$Sample.Source.ID.processed %in% Europe_samples, ]
metaphlan_Europe <- as.data.frame(t(metaphlan[,colnames(metaphlan) %in% Europe_samples, ]))

# 使用示例函数调用原代码中的分析
overall_obesity_analysis <- run_comprehensive_analysis(
  metaphlan_data = metaphlan,  # 注意：原代码中metaphlan_inhouse是转置的
  metadata = metadata,
  analysis_type = "obesity",
  config = config,
  data_name = "overall"
)

overall_disease_analysis <- run_comprehensive_analysis(
  metaphlan_data = metaphlan,
  metadata = metadata,
  analysis_type = "disease",
  config = config,
  data_name = "overall"
)

# BMI关联分析
# 假设已有metaphlan_data和metadata
filtering_result <- smart_data_filtering(
  metaphlan_data = t(metaphlan),
  metadata = metadata,
  config = config
)

overall_bmi_cor_analysis <- perform_bmi_analysis(
  filtered_data = filtering_result$filtered_data,
  metadata = metadata,
  config = config,
  bmi_var = "BMI",
  by_disease_group = TRUE,     # 开启疾病分层
  disease_var = "Disease",     # 指定分组变量
  include_overall = TRUE       # 包含全样本（默认）
)



# 1. 定义固定顺序
disease_order <- c("CTRL" = 1, "AA" = 2, "CRC" = 3)
obesity_order <- c("Normal" = 1, "Ow" = 2,"Ob" = 3)


# 使用示例函数调用原代码中的分析
Asia_obesity_analysis <- run_comprehensive_analysis(
  metaphlan_data = metaphlan_Asia,  # 注意：原代码中metaphlan_inhouse是转置的
  metadata = metadata_Asia,
  analysis_type = "obesity",
  config = config,
  data_name = "Asia"
)

Asia_disease_analysis <- run_comprehensive_analysis(
  metaphlan_data = metaphlan_Asia,
  metadata = metadata_Asia,
  analysis_type = "disease",
  config = config,
  data_name = "Asia"
)


# 假设已有metaphlan_data和metadata
filtering_result_Asia <- smart_data_filtering(
  metaphlan_data = t(metaphlan_Asia),
  metadata = metadata_Asia,
  config = config
)


Asia_bmi_cor_analysis <- perform_bmi_analysis(
  filtered_data = filtering_result_Asia$filtered_data,
  metadata = metadata_Asia,
  config = config,
  bmi_var = "BMI",
  by_disease_group = TRUE,     # 开启疾病分层
  disease_var = "Disease",     # 指定分组变量
  include_overall = TRUE       # 包含全样本（默认）
)

Europe_obesity_analysis <- run_comprehensive_analysis(
  metaphlan_data = metaphlan_Europe,  # 注意：原代码中metaphlan_inhouse是转置的
  metadata = metadata_Europe,
  analysis_type = "obesity",
  config = config,
  data_name = "Europe"
)

Europe_disease_analysis <- run_comprehensive_analysis(
  metaphlan_data = metaphlan_Europe,
  metadata = metadata_Europe,
  analysis_type = "disease",
  config = config,
  data_name = "Europe"
)


# 假设已有metaphlan_data和metadata
filtering_result_Europe <- smart_data_filtering(
  metaphlan_data = t(metaphlan_Europe),
  metadata = metadata_Europe,
  config = config
)

Europe_bmi_cor_analysis <- perform_bmi_analysis(
  filtered_data = filtering_result_Europe$filtered_data,
  metadata = metadata_Europe,
  config = config,
  bmi_var = "BMI",
  by_disease_group = TRUE,     # 开启疾病分层
  disease_var = "Disease",     # 指定分组变量
  include_overall = TRUE       # 包含全样本（默认）
)

# ====================== 第二步：更新结构化结果列表（包含BMI相关性分析） ======================
all_analysis_results <- list(
  # 1. 全球（overall）分析结果
  overall = list(
    obesity_analysis = overall_obesity_analysis,
    disease_analysis = overall_disease_analysis,
    disease_DA = list(
      adjust_BMI = overall_disease_analysis$analysis_results$disease_analysis_adjust_BMI,
      no_BMI = overall_disease_analysis$analysis_results$disease_analysis_no_BMI
    ),
    bmi_cor_analysis = overall_bmi_cor_analysis  # 新增BMI相关性分析
  ),
  
  # 2. 亚洲（Asia）分析结果
  Asia = list(
    obesity_analysis = Asia_obesity_analysis,
    disease_analysis = Asia_disease_analysis,
    disease_DA = list(
      adjust_BMI = Asia_disease_analysis$analysis_results$disease_analysis_adjust_BMI,
      no_BMI = Asia_disease_analysis$analysis_results$disease_analysis_no_BMI
    ),
    bmi_cor_analysis = Asia_bmi_cor_analysis  # 新增BMI相关性分析
  ),
  
  # 3. 欧洲（Europe）分析结果
  Europe = list(
    obesity_analysis = Europe_obesity_analysis,
    disease_analysis = Europe_disease_analysis,
    disease_DA = list(
      adjust_BMI = Europe_obesity_analysis$analysis_results$disease_analysis_adjust_BMI,
      no_BMI = Europe_obesity_analysis$analysis_results$disease_analysis_no_BMI
    ),
    bmi_cor_analysis = Europe_bmi_cor_analysis  # 新增BMI相关性分析
  ),
  
  # 可选：保存分析配置信息（便于后续复现）
  config_info = list(
    disease_order = disease_order,
    obesity_order = obesity_order,
    base_config = config,
    bmi_var = "BMI"  # 补充BMI变量名配置
  )
)

# ====================== 第三步：重新保存为RDS文件（覆盖或新命名） ======================
# 保存
message("\n=== 保存结果 ===")
rds_path <- file.path(output_dir, paste0("microbiome_analysis_all_results_BMI_Group1_", timestamp, ".rds"))
saveRDS(all_analysis_results, file = rds_path)
message(paste("结果已保存:", rds_path))

# 打印保存成功提示
cat(paste("所有分析结果（含BMI相关性分析）已保存到：", rds_path, "\n"))

# ====================== 第三步：加载验证（可选） ======================
# 加载保存的RDS文件
all_analysis_results <- readRDS("D:/Tianhui/proj_ob_AA_CRC/03.results/04_differential_abundance/microbiome_analysis_all_results_BMI_Group1_20260116.rds")

# 验证加载结果（检查列表结构）
cat("\n=== 验证加载结果 ===\n")
# 查看主列表的层级
str(all_analysis_results, max.level = 3)

# 示例：提取Asia区域调整BMI的疾病差异分析结果
asia_da_adjust <- loaded_results$Asia$disease_DA$adjust_BMI
cat(paste("\nAsia区域调整BMI的疾病差异分析结果维度：", dim(asia_da_adjust)[1], "行 ×", dim(asia_da_adjust)[2], "列\n"))

