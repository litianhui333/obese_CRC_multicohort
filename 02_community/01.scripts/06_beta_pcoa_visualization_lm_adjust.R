# ==============================================================================
# 脚本名称: 06_beta_pcoa_visualization_pro.R
# 功能描述: 
#   1. 坐标层级批次校正: 使用 PCoA 轴坐标对 Study.name 求残差消除批次影响
#   2. 增强可视化: 包含样本点、95% 置信椭圆，并注明校正后的解释度
# ==============================================================================

# ======================== 1. 路径与环境设置 ========================
base_dir   <- "D:/Tianhui/proj_ob_AA_CRC"
input_dir  <- file.path(base_dir, "02.processed_data/02_aligned_data")
sp_dir     <- file.path(base_dir, "02.processed_data/03_metaphlan_split_data")
output_dir <- file.path(base_dir, "03.results/03_pcoa_plots")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 动态搜索最新的数据文件
meta_file <- sort(list.files(input_dir, pattern = "metadata_aligned.*\\.rds", full.names = TRUE), decreasing = TRUE)[1]
sp_file   <- sort(list.files(sp_dir, pattern = "metaphlan_species_renormalized.*\\.rds", full.names = TRUE), decreasing = TRUE)[1]

if(is.na(meta_file) || is.na(sp_file)) {
  stop("错误: 找不到输入文件，请检查路径。")
}

# ======================== 2. 加载必要包 ========================
library(tidyverse)
library(vegan)
library(ggforce)
library(patchwork)
library(ggrepel)

# ======================== 3. 数据加载与预处理 ========================
message("正在加载并预处理数据...")
metadata <- readRDS(meta_file)
metaphlan_raw <- readRDS(sp_file)
rownames(metadata) <- metadata$Sample.Source.ID.processed
metaphlan <- as.data.frame(t(metaphlan_raw))

# 普及率过滤 (Prevalence Filter: 10%)
prev_threshold <- 0.1
metaphlan_filt <- metaphlan[, colSums(metaphlan > 0) > (nrow(metaphlan) * prev_threshold)]

# 样本对齐
common_samples  <- intersect(rownames(metaphlan_filt), rownames(metadata))
metaphlan_final <- metaphlan_filt[common_samples, ]
metadata_final  <- metadata[common_samples, ]

# ======================== 4. 绘图辅助设置 ========================

get_color_schemes <- function() {
  list(disease = c("CTRL" = "#2E9FDF", "AA" = "#F39C12", "CRC" = "#E74C3C"))
}

cns_theme <- function(base_size = 11) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1.0),
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# ======================== 5. 核心函数：坐标残差校正绘图 ========================

plot_pcoa_residual_adjusted <- function(abundance_matrix, metadata, bmi_col, batch_col = "Study.name",
                                        subset_col = NULL, disease_col = "Disease", colors = get_color_schemes()) {
  
  plot_list <- list()
  
  # 定义分析子集 (Overall 或按 Continent)
  if (is.null(subset_col)) {
    subsets <- list(Overall = rownames(metadata))
  } else {
    sub_vals <- as.character(unique(na.omit(metadata[[subset_col]])))
    subsets <- list()
    for (v in sub_vals) subsets[[v]] <- rownames(metadata[which(metadata[[subset_col]] == v), ])
  }
  
  for (sub_name in names(subsets)) {
    curr_samples <- subsets[[sub_name]]
    meta_sub <- metadata[curr_samples, ]
    abd_sub <- abundance_matrix[curr_samples, ]
    
    # 进一步按 BMI 分层
    bmi_groups <- unique(na.omit(meta_sub[[bmi_col]]))
    
    for (bg in bmi_groups) {
      final_samples <- rownames(meta_sub[which(meta_sub[[bmi_col]] == bg), ])
      if (length(final_samples) < 10) next
      
      # 1. 计算 Bray-Curtis 距离
      local_dist <- vegdist(abd_sub[final_samples, ], method = "bray")
      
      # 2. 计算原始 PCoA 坐标
      pcoa_res <- cmdscale(local_dist, k = 2, eig = TRUE)
      
      # 计算解释度 (基于原始 PCoA 结果)
      eig <- pcoa_res$eig
      var_exp <- eig / sum(eig[eig > 0])
      pc1_var <- round(var_exp[1] * 100, 1)
      pc2_var <- round(var_exp[2] * 100, 1)
      
      df <- as.data.frame(pcoa_res$points)
      colnames(df) <- c("PCo1", "PCo2")
      df$Sample <- rownames(df)
      df <- merge(df, meta_sub[final_samples, ], by.x = "Sample", by.y = "row.names")
      
      # 3. 核心步骤：通过线性模型残差移除 Study 的影响
      if (length(unique(df[[batch_col]])) > 1) {
        df$PCo1_adj <- resid(lm(as.formula(paste("PCo1 ~", batch_col)), data = df))
        df$PCo2_adj <- resid(lm(as.formula(paste("PCo2 ~", batch_col)), data = df))
        
        # 修改标签，注明是校正后的，但保留原始解释度作为参考
        x_lab <- paste0("PCo1 (Study-corrected) [", pc1_var, "%]")
        y_lab <- paste0("PCo2 (Study-corrected) [", pc2_var, "%]")
      } else {
        df$PCo1_adj <- df$PCo1
        df$PCo2_adj <- df$PCo2
        x_lab <- paste0("PCo1 [", pc1_var, "%]")
        y_lab <- paste0("PCo2 [", pc2_var, "%]")
      }
      
      # 4. 绘图 (移除质心和连线)
      p <- ggplot(df, aes(x = PCo1_adj, y = PCo2_adj, color = .data[[disease_col]])) +
        # 绘制 95% 置信椭圆
        stat_ellipse(aes(fill = .data[[disease_col]]), level = 0.68, 
                     geom = "polygon", alpha = 0.1, linewidth = 0.5, linetype = 1.2) +
        # 绘制样本点
        geom_point(size = 0.5, alpha = 0.3) +
        scale_color_manual(values = colors$disease) +
        scale_fill_manual(values = colors$disease) +
        labs(x = x_lab, y = y_lab,
             title = paste0("Residual PCoA: ", sub_name),
             subtitle = paste0("BMI: ", bg, " (n=", nrow(df), ") | 95% CI")) +
        cns_theme()
      
      plot_list[[paste0(sub_name, "_", bg)]] <- p
    }
  }
  return(plot_list)
}

# ======================== 6. 执行分析与保存 ========================

message("开始生成分析图表...")
save_path <- file.path(output_dir, "residual_adjusted_analysis")
if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)

# 1. Overall 分析
res_overall <- plot_pcoa_residual_adjusted(metaphlan_final, metadata_final, "BMI_Group1", subset_col = NULL)

# 2. Continent 分析
res_continent <- plot_pcoa_residual_adjusted(metaphlan_final, metadata_final, "BMI_Group1", subset_col = "Continent")

# 合并保存
all_plots <- c(res_overall, res_continent)

for (n in names(all_plots)) {
  file_name <- paste0(gsub(" ", "_", n), "_residual_adj.pdf")
  ggsave(file.path(save_path, file_name), all_plots[[n]], width = 4, height = 3)
}

message("分析成功完成。所有结果已保存至: ", save_path)


# ======================== 6. 新增函数：微生物组轨迹图 (9组质心 + 置信区间线) ========================

plot_pcoa_trajectory_analysis <- function(abundance_matrix, metadata, bmi_col, batch_col = "Study.name",
                                          subset_col = NULL, disease_col = "Disease", colors = get_color_schemes()) {
  
  plot_list <- list()
  
  # 定义分析子集
  if (is.null(subset_col)) {
    subsets <- list(Overall = rownames(metadata))
  } else {
    sub_vals <- as.character(unique(na.omit(metadata[[subset_col]])))
    subsets <- list()
    for (v in sub_vals) subsets[[v]] <- rownames(metadata[which(metadata[[subset_col]] == v), ])
  }
  
  for (sub_name in names(subsets)) {
    curr_samples <- subsets[[sub_name]]
    if (length(curr_samples) < 15) next # 需要足够的样本量计算9个质心
    
    meta_sub <- metadata[curr_samples, ]
    abd_sub <- abundance_matrix[curr_samples, ]
    
    # 1. 计算局部 PCoA
    local_dist <- vegdist(abd_sub, method = "bray")
    pcoa_res <- cmdscale(local_dist, k = 2, eig = TRUE)
    
    eig <- pcoa_res$eig
    var_exp <- eig / sum(eig[eig > 0])
    pc1_var <- round(var_exp[1] * 100, 1)
    pc2_var <- round(var_exp[2] * 100, 1)
    
    df <- as.data.frame(pcoa_res$points)
    colnames(df) <- c("PCo1", "PCo2")
    df$Sample <- rownames(df)
    df <- merge(df, meta_sub, by.x = "Sample", by.y = "row.names")
    
    # 2. 坐标残差校正 (移除批次影响)
    if (length(unique(df[[batch_col]])) > 1) {
      df$PCo1_adj <- resid(lm(as.formula(paste("PCo1 ~", batch_col)), data = df))
      df$PCo2_adj <- resid(lm(as.formula(paste("PCo2 ~", batch_col)), data = df))
      x_lab <- paste0("PCo1 (Study-corrected) [", pc1_var, "%]")
      y_lab <- paste0("PCo2 (Study-corrected) [", pc2_var, "%]")
    } else {
      df$PCo1_adj <- df$PCo1
      df$PCo2_adj <- df$PCo2
      x_lab <- paste0("PCo1 [", pc1_var, "%]")
      y_lab <- paste0("PCo2 [", pc2_var, "%]")
    }
    
    # 3. 计算 9 个组的质心与标准误
    df$Disease_BMI <- paste0(df[[disease_col]], "_", df[[bmi_col]])
    
    centroids <- df %>%
      group_by(Disease_BMI, .data[[disease_col]], .data[[bmi_col]]) %>%
      summarise(
        Axis.1 = mean(PCo1_adj),
        Axis.2 = mean(PCo2_adj),
        se.1 = sd(PCo1_adj) / sqrt(n()),
        se.2 = sd(PCo2_adj) / sqrt(n()),
        .groups = 'drop'
      )
    
    # 因子排序
    centroids[[bmi_col]] <- factor(centroids[[bmi_col]], levels = c("Normal", "Ow", "Ob"))
    centroids[[disease_col]] <- factor(centroids[[disease_col]], levels = c("CTRL", "AA", "CRC"))
    centroids <- centroids %>% arrange(.data[[disease_col]], .data[[bmi_col]])
    
    # 4. 绘图
    p <- ggplot(centroids, aes(x = Axis.1, y = Axis.2)) +
      # 置信区间线
      geom_errorbar(aes(ymin = Axis.2 - 1.96*se.2, ymax = Axis.2 + 1.96*se.2, color = .data[[disease_col]]), width = 0, alpha = 0.5) +
      geom_errorbarh(aes(xmin = Axis.1 - 1.96*se.1, xmax = Axis.1 + 1.96*se.1, color = .data[[disease_col]]), height = 0, alpha = 0.5) +
      # 轨迹线
      geom_path(aes(group = .data[[disease_col]], color = .data[[disease_col]]), 
                arrow = arrow(length = unit(0.3, "cm"), type = "closed"), 
                linewidth = 1, linetype = "dashed") +
      # 质心点
      geom_point(aes(color = .data[[disease_col]], shape = .data[[bmi_col]]), size = 6) +
      geom_text_repel(aes(label = Disease_BMI, color = .data[[disease_col]]), size = 3, max.overlaps = 20) +
      
      scale_color_manual(values = colors$disease) +
      scale_shape_manual(values = c("Normal"=19, "Ow"=17, "Ob"=15)) + 
      
      labs(x = x_lab, y = y_lab,
           title = paste0("Microbiome Trajectory: ", sub_name),
           subtitle = "Centroids with 95% CI (Mean +/- 1.96*SE) | Batch Corrected") +
      cns_theme()
    
    plot_list[[sub_name]] <- p
  }
  return(plot_list)
}


# ======================== 7. 执行分析与保存 ========================

message("开始生成分析图表...")

# 1. 生成 Overall 轨迹图
traj_overall <- plot_pcoa_trajectory_analysis(metaphlan_final, metadata_final, "BMI_Group1", subset_col = NULL)

# 2. 生成按 Continent 分组的轨迹图
traj_continent <- plot_pcoa_trajectory_analysis(metaphlan_final, metadata_final, "BMI_Group1", subset_col = "Continent")

# 合并并保存所有轨迹图
all_traj_plots <- c(traj_overall, traj_continent)
for (n in names(all_traj_plots)) {
  ggsave(file.path(save_path, paste0("trajectory_", gsub(" ", "_", n), ".pdf")), all_traj_plots[[n]], width = 6, height = 4)
}



message("分析成功完成。轨迹图已包含各 Continent 子集并保存至 residual_adjusted_analysis。")
