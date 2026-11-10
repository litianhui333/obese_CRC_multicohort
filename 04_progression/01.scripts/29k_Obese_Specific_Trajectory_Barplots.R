# ==============================================================================
# 脚本名称: 29k_Obese_Specific_Trajectory_Barplots.R
# 功能描述:
#   对 Ob-specific 菌群绘制三个 BMI 组 (Nw/Ow/Ob) 的 CTRL -> AA -> CRC 轨迹柱状图
#   - 菌名单直接来自主图 Script 29 生成的 Feature_List_Ob_Clustered.csv (保证一致)
#   - Enriched: CSV 中 overlap=="Ob-specific" & coef_crc>0 (24个，实际匹配到多少展示多少)
#   - Depleted: CSV 中 overlap=="Ob-specific" & coef_crc<0 (7个)
#   - 均值柱 + SEM 误差棒 + 均值连接线
#   - 命名使用与 Script 29 完全一致的 format_sgb_names()
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(MMUPHin)
  library(tibble)
  library(patchwork)
})

# ==============================================================================
# 0. 路径 & 工具函数
# ==============================================================================
base_dir  <- "G:/My Drive/proj_ob_AA_CRC"
source(file.path(base_dir, "01.scripts/DA_analysis_figure_utils_functions.R"))
source(file.path(base_dir, "01.scripts/DA_function.R"))

input_dir  <- file.path(base_dir, "03.results/04_differential_abundance")
data_dir   <- file.path(base_dir, "02.processed_data/02_aligned_data")
sp_dir     <- file.path(base_dir, "02.processed_data/03_metaphlan_split_data")
heat_dir   <- file.path(base_dir, "03.results/04_differential_abundance/heatmap_final/BMI_Stratified_Disease_Heatmaps")
output_dir <- heat_dir
cache_dir  <- file.path(base_dir, "03.results/04_differential_abundance/heatmap_final/cache")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(cache_dir))  dir.create(cache_dir,  recursive = TRUE)

# ==============================================================================
# 1. 直接从主图 CSV 读取权威菌名单 (与 Script 29 100% 一致)
# ==============================================================================
message("=== 1. 读取权威菌名单 (Feature_List_Ob_Refined.csv == Script 29 实际写入的文件) ===")
# 注意: Script 29 在第 613 行写入的是 Feature_List_{BMI}_Refined.csv
csv_file <- file.path(heat_dir, "Feature_List_Ob_Refined.csv")
if (!file.exists(csv_file)) stop("❌ 未找到 Feature_List_Ob_Refined.csv，请先运行 Script 29。")

res_file  <- sort(list.files(input_dir, pattern="microbiome_analysis_all_results_BMI_Group1_0.7.*\\.rds", full.names=TRUE), decreasing=TRUE)[1]
all_analysis_results <- readRDS(res_file)
data_info <- all_analysis_results$overall$disease_analysis$data_info
fa <- data_info$filtered_data
meta <- data_info$metadata_filtered

calculate_sample_direction <- function(features, fa, meta) {
  ob_mask <- meta$BMI_Group1 == "Ob"
  ctrl_idx <- which(ob_mask & meta$Disease == "CTRL")
  crc_idx <- which(ob_mask & meta$Disease == "CRC")
  
  result <- c()
  for (feat in features) {
    if (feat %in% rownames(fa)) {
      ctrl_vals <- as.numeric(fa[feat, ctrl_idx])
      crc_vals <- as.numeric(fa[feat, crc_idx])
      ctrl_pos <- ctrl_vals[!is.na(ctrl_vals) & ctrl_vals > 0]
      crc_pos <- crc_vals[!is.na(crc_vals) & crc_vals > 0]
      
      if (length(ctrl_pos) > 0 && length(crc_pos) > 0) {
        ctrl_mean <- mean(ctrl_pos)
        crc_mean <- mean(crc_pos)
        direction <- ifelse(crc_mean > ctrl_mean, "UP", "DOWN")
        result[feat] <- direction
      }
    }
  }
  result
}

feat_df <- read.csv(csv_file, stringsAsFactors=FALSE)
cat("Total rows in CSV:", nrow(feat_df), "\n")
cat("Overlap categories:\n"); print(table(feat_df$overlap))

sample_dir <- calculate_sample_direction(feat_df$feature, fa, meta)

feat_df$sample_dir <- sample_dir[feat_df$feature]
feat_df$direction_match <- ifelse(
  (feat_df$coef_crc > 0 & feat_df$sample_dir == "UP") | 
  (feat_df$coef_crc < 0 & feat_df$sample_dir == "DOWN"),
  "Match", "Mismatch"
)

cat("\n=== Direction check ===\n")
print(table(feat_df$direction_match))

ob_en_df <- feat_df[grepl("Ob-specific", feat_df$overlap) | feat_df$overlap == "Overlap with Ow", ]
ob_en_df <- ob_en_df[ob_en_df$direction_match == "Match" & ob_en_df$sample_dir == "UP", ]
ob_en_full <- ob_en_df$feature

b_de_df <- feat_df[grepl("Ob-specific", feat_df$overlap) | feat_df$overlap == "Overlap with Ow", ]
b_de_df <- b_de_df[b_de_df$direction_match == "Match" & b_de_df$sample_dir == "DOWN", ]
ob_de_full <- b_de_df$feature

message(sprintf("   >>> Enriched (Ob-specific + Overlap Ow + matched + UP): %d | Depleted (Ob-specific + Overlap Ow + matched + DOWN): %d",
                length(ob_en_full), length(ob_de_full)))

# 状态标签映射
status_map <- c(
  setNames(rep("Enriched", length(ob_en_full)), ob_en_full),
  setNames(b_de_df$overlap,                     ob_de_full)  # overlap 标签作为 Status
)

# ==============================================================================
# 2. 加载原始数据 (与 Script 29 完全相同)
# ==============================================================================
message("=== 2. 加载原始数据 ===")
meta_file <- sort(list.files(data_dir, pattern="metadata_aligned.*\\.rds", full.names=TRUE), decreasing=TRUE)[1]
sp_file   <- sort(list.files(sp_dir,   pattern="metaphlan_strain_renormalized.*\\.rds", full.names=TRUE), decreasing=TRUE)[1]

metadata             <- readRDS(meta_file)
metaphlan            <- readRDS(sp_file)

rownames(metadata)   <- metadata$Sample.Source.ID.processed
metadata$Disease     <- factor(metadata$Disease,    levels=c("CTRL","AA","CRC"))
metadata$BMI_Group1  <- factor(metadata$BMI_Group1, levels=c("Normal","Ow","Ob"))
sex_col <- if ("Sex" %in% colnames(metadata)) "Sex" else "Gender"

# ==============================================================================
# 3. 与 Script 29 完全相同的 MMUPHin 批次校正
# ==============================================================================
message("=== 3. MMUPHin 批次校正 (全局，与 Script 29 完全一致) ===")

# 同 Script 29: 用 get_robust_features_detailed 识别 all_targets
get_robust_features_detailed <- function(data_source, region, bmi_group, comp,
                                         p_cutoff=0.01, q_cutoff=0.1, weight_limit=0.7) {
  list_name <- paste0("disease_analysis_", bmi_group)
  res_obj <- tryCatch(
    data_source[[region]][["disease_analysis"]][["analysis_results"]][[list_name]][["all_results"]][[comp]],
    error=function(e) NULL)
  if (is.null(res_obj)) return(NULL)
  df_meta   <- res_obj[["meta_fits"]]; df_maaslin <- res_obj[["maaslin_fits"]]
  if (is.null(df_meta) | is.null(df_maaslin)) return(NULL)
  weight_cols <- grep("^weight_", colnames(df_meta), value=TRUE)
  w_mat <- as.matrix(df_meta[, weight_cols]); w_mat[is.na(w_mat)] <- 0
  perform_loo_rescue <- function(f_name, weight_row, original_coef_sign) {
    tryCatch({
      max_idx <- which.max(weight_row); if(length(max_idx)==0) return(list(saved=FALSE))
      exclude_study <- gsub("^weight_","",names(weight_row)[max_idx])
      remaining_studies <- names(df_maaslin)[names(df_maaslin) != exclude_study]
      coefs <- c(); ses <- c()
      for (s in remaining_studies) {
        df_s <- df_maaslin[[s]]
        if ("feature" %in% colnames(df_s) && f_name %in% df_s$feature) {
          idx <- which(df_s$feature==f_name); val_c <- df_s$coef[idx]; val_se <- df_s$stderr[idx]
          if (!is.na(val_c) && !is.na(val_se) && val_se>1e-10) { coefs <- c(coefs,val_c); ses <- c(ses,val_se) }
        }
      }
      if (length(coefs)<1) return(list(saved=FALSE))
      weights <- 1/(ses^2); valid_w <- is.finite(weights)&!is.na(weights)
      if (sum(valid_w)==0) return(list(saved=FALSE))
      weights <- weights[valid_w]; coefs <- coefs[valid_w]; sum_w <- sum(weights)
      if (sum_w==0) return(list(saved=FALSE))
      combined_coef <- sum(coefs*weights)/sum_w; combined_se <- sqrt(1/sum_w)
      combined_p <- 2*(1-pnorm(abs(combined_coef/combined_se)))
      if (is.na(combined_p)||sign(combined_coef)!=original_coef_sign) return(list(saved=FALSE))
      if (combined_p<0.05) return(list(saved=TRUE,coef=combined_coef,pval=combined_p))
      return(list(saved=FALSE))
    }, error=function(e) list(saved=FALSE))
  }
  audit_df <- data.frame(
    feature=rownames(df_meta),
    Max_Weight=apply(w_mat,1,function(x){s=sum(x,na.rm=TRUE);if(s<=1e-10)1 else max(x,na.rm=TRUE)/s}),
    pval=df_meta$pval, qval=df_meta$qval.fdr, coef=df_meta$coef, stringsAsFactors=FALSE)
  result_list <- list()
  for (i in seq_len(nrow(audit_df))) {
    row <- audit_df[i,]; f <- row$feature
    is_candidate <- !is.na(row$pval) && (row$pval<p_cutoff | (!is.na(row$qval) && row$qval<q_cutoff))
    if (!is_candidate) next
    status <- "Discard"; final_coef <- row$coef; mw <- row$Max_Weight; if(is.na(mw)) mw <- 1
    if (mw<weight_limit) { status <- "Robust" } else {
      loo_res <- perform_loo_rescue(f, w_mat[f,], sign(row$coef))
      if (loo_res$saved) { status <- "Rescued"; final_coef <- loo_res$coef }
    }
    if (status!="Discard") result_list[[f]] <- data.frame(feature=f,pval=row$pval,qval=row$qval,coef=final_coef,robust_status=status)
  }
  if (length(result_list)==0) return(NULL)
  return(do.call(rbind, result_list))
}

df_robust_Nw <- get_robust_features_detailed(all_analysis_results,"overall","Normal","CTRL_vs_CRC")
df_robust_Ow <- get_robust_features_detailed(all_analysis_results,"overall","Ow",    "CTRL_vs_CRC")
df_robust_Ob <- get_robust_features_detailed(all_analysis_results,"overall","Ob",    "CTRL_vs_CRC")

pass_nw <- if(!is.null(df_robust_Nw)) df_robust_Nw$feature else character(0)
pass_ow <- if(!is.null(df_robust_Ow)) df_robust_Ow$feature else character(0)
pass_ob <- if(!is.null(df_robust_Ob)) df_robust_Ob$feature else character(0)
all_targets <- unique(c(pass_nw, pass_ow, pass_ob))

config <- create_analysis_config(
  transform_method="AST", analysis_method="HS", BMI_Group_var="BMI_Group1",
  covariates=c("Age","Sex"), study_prevalence_cut=0.1, abundance_threshold=0.0001,
  min_studies_ratio=0.5, pval_type="both", sig_criterion="or",
  pval_threshold=0.02, qval_threshold=0.2, is_inhouse=FALSE, force_keep=all_targets)

filtering_result   <- smart_data_filtering(t(metaphlan), metadata, config)
metaphlan_filtered <- filtering_result$filtered_data

common_samples <- intersect(rownames(metadata), colnames(metaphlan_filtered))
meta_clean     <- metadata[common_samples,] %>%
  drop_na(all_of(c("Study.name","Disease","BMI_Group1","Age",sex_col)))
abd_clean      <- metaphlan_filtered[, rownames(meta_clean)]
meta_clean$Study.name <- as.factor(meta_clean$Study.name)

try(unlink("adjust_batch_diagnostic.pdf", force=TRUE))
fit_adjust <- adjust_batch(
  feature_abd=abd_clean/100, batch="Study.name",
  covariates=c("Disease","BMI_Group1","Age",sex_col),
  data=meta_clean, control=list(verbose=FALSE))
abd_corrected <- fit_adjust$feature_abd_adj * 100  # 批次校正后的 %

# ==============================================================================
# 4. 名字对齐: 将 CSV 中的全名映射到批次校正矩阵的行名
# ==============================================================================
message("=== 4. 名字对齐 ===")
all_target_csv <- c(ob_en_full, ob_de_full)

# 建立从 CSV 全名 -> 矩阵行名 的映射
map_to_matrix <- function(targets, mat_rownames) {
  result <- list()
  for (t in targets) {
    # 精确匹配
    if (t %in% mat_rownames) { result[[t]] <- t; next }
    # SGB 编号匹配 (最可靠)
    sgb <- regmatches(t, regexpr("SGB[0-9]+", t))
    if (length(sgb) > 0) {
      m <- grep(sgb, mat_rownames, value=TRUE, fixed=TRUE)
      if (length(m) >= 1) { result[[t]] <- m[1]; next }
    }
    # 种名匹配
    sp_name <- gsub("\\|t__.*","", gsub(".*s__","s__",t))
    sp_key  <- gsub("s__","", sp_name)
    m <- grep(sp_key, mat_rownames, value=TRUE, fixed=TRUE)
    if (length(m) >= 1) { result[[t]] <- m[1] }
  }
  return(result)
}

mapping <- map_to_matrix(all_target_csv, rownames(abd_corrected))
found    <- names(mapping)
missing  <- setdiff(all_target_csv, found)

message(sprintf("   >>> CSV %d 个目标菌 -> 批次矩阵匹配: %d 个 (%d 个未匹配)",
                length(all_target_csv), length(found), length(missing)))
if (length(missing) > 0) {
  message("   >>> 未匹配 (将跳过): ", paste(format_sgb_names(missing), collapse="; "))
}

# 最终用于绘图的矩阵行名
matched_rows <- unlist(mapping)  # 矩阵中的实际行名
names(matched_rows) <- found     # CSV 中的完整名 (key)

# ==============================================================================
# 5. 提取三组数据并重塑
# ==============================================================================
message("=== 5. 提取三组数据 ===")
long_list <- list()

for (bmi in c("Normal","Ow","Ob")) {
  meta_bmi  <- meta_clean[meta_clean$BMI_Group1 == bmi, ]
  samps_bmi <- rownames(meta_bmi)
  samps_bmi <- intersect(samps_bmi, colnames(metaphlan_filtered))
  meta_bmi  <- meta_bmi[samps_bmi, ]
  
  for (csv_name in found) {
    mat_row    <- matched_rows[csv_name]
    raw_abund  <- as.numeric(metaphlan_filtered[mat_row, samps_bmi])
    # AST transformation for percentage abundance.
    abund_vec  <- asin(sqrt(pmin(pmax(raw_abund, 0) / 100, 1)))
    # 使用与 Script 29 完全相同的 format_sgb_names() 命名
    disp_name  <- format_sgb_names(csv_name)
    
    tmp <- data.frame(
      SampleID  = samps_bmi,
      Disease   = meta_bmi$Disease,
      BMI_Group = bmi,
      Species   = disp_name,
      Feature   = csv_name,
      Abundance = abund_vec,
      Status    = status_map[csv_name],
      stringsAsFactors = FALSE
    )
    long_list[[paste(bmi, csv_name)]] <- tmp
  }
  message(sprintf("   %s: %d 样本", bmi, length(samps_bmi)))
}

plot_df <- as.data.frame(do.call(rbind, long_list))
plot_df$Disease      <- factor(plot_df$Disease,   levels=c("CTRL","AA","CRC"))
plot_df$BMI_Group    <- factor(plot_df$BMI_Group, levels=c("Normal","Ow","Ob"))
# Status: Enriched / Ob-specific / Overlap with Nw / Overlap with both Nw and Ow
# 用于绘图标注 depleted 菌的 overlap 类型
plot_df$Is_Depleted <- plot_df$Status != "Enriched"

# ==============================================================================
# 6. 核心绘图函数 — CNS 极简风格 (可复用) | Core Plotting Function — CNS Style
# ==============================================================================

#' plot_bmi_trajectory_barplots  (CNS style)
#'
#' 【功能】每个菌独立小图, 三个 BMI 组并排 (固定 Y 轴), Patchwork 拼接输出。
#' 【Function】One subplot per taxon, 3 BMI groups side-by-side (fixed Y-axis),
#'   assembled with patchwork.
#'
#' 图形特征 | Visual style:
#'   - 半透明柱 (α=0.65) + 无帽误差棒 | Semi-transparent bars + capless error bars
#'   - firebrick 实线连接均值点 | Solid firebrick trajectory line + white-fill circles
#'   - theme_classic CNS 极简主题 | CNS minimalist theme_classic
#'   - scales = "fixed": 三组共享 Y 轴范围 | Fixed Y-axis across 3 BMI groups
#'
#' @param df            长格式数据框 | Long-format data frame
#'                      必须列 | Required cols: Disease, BMI_Group, Species, Abundance, Status
#' @param status_type   "Enriched" 或其他任意字符串 (代表降低菌)
#'                      "Enriched" or any other string (for depleted taxa)
#' @param filename      输出文件名 (相对 output_dir) | Output filename
#' @param output_dir    输出目录 | Output directory
#' @param sort_by_bmi   菌排序依据的 BMI 组 | BMI group for species ordering. Default "Ob"
#' @param disease_colors 疾病颜色向量 | Disease stage colors
#' @param page_title    总标题 (NULL 自动) | Overall page title (NULL = auto)
#' @param ncol_plots    每行小图数量 | Subplots per row. Default 3
#' @param subplot_w     每列小图宽度 (英寸) | Width per subplot column. Default 3.2
#' @param subplot_h     每行小图高度 (英寸) | Height per subplot row.  Default 2.2

plot_bmi_trajectory_barplots <- function(
    df,
    status_type,
    filename,
    output_dir,
    sort_by_bmi    = "Ob",
    disease_colors = c("CTRL" = "#2E9FDF", "AA" = "#F39C12", "CRC" = "#E74C3C"),
    page_title     = NULL,
    ncol_plots     = 3,
    subplot_w      = 3.2,
    subplot_h      = 2.2
) {
  # ---- 1. 数据过滤 | Subset by direction ----
  if (status_type == "Enriched") {
    sub_df <- as.data.frame(df[df$Status == "Enriched",  , drop = FALSE])
  } else {
    sub_df <- as.data.frame(df[df$Status != "Enriched",  , drop = FALSE])
  }
  if (nrow(sub_df) == 0) {
    message(sprintf("   >>> [SKIP] No %s taxa.", status_type)); return(invisible(NULL))
  }

  n_species <- length(unique(sub_df$Species))
  message(sprintf("   >>> Plotting %s: %d taxa (CNS style)", status_type, n_species))

  # ---- 2. 菌名排序 | Order species by mean CRC abundance in sort_by_bmi ----
  sort_fn       <- if (status_type == "Enriched") dplyr::desc else identity
  species_order <- sub_df[sub_df$BMI_Group == sort_by_bmi & sub_df$Disease == "CRC", ] %>%
    group_by(Species) %>%
    summarise(m = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(sort_fn(m)) %>% pull(Species)
  remaining     <- setdiff(unique(sub_df$Species), species_order)
  sub_df$Species <- factor(sub_df$Species, levels = c(species_order, remaining))

  # ---- 3. BMI 分面标签 | BMI facet labels ----
  bmi_labs <- c("Normal" = "Normal", "Ow" = "Overweight", "Ob" = "Obese")

  # ---- 4. 总标题 | Overall page title ----
  if (is.null(page_title)) {
    dir_word   <- if (status_type == "Enriched") "Enriched (+)" else "Depleted (-)"
    page_title <- paste0("Ob-Specific Microbial Signatures - ", dir_word, " (n = ", n_species, ")")
  }

  # ---- 5. 逐菌绘图 | Per-taxon subplot loop ----
  plot_list <- list()
  for (sp in levels(sub_df$Species)) {
    sp_df <- sub_df[sub_df$Species == sp, ]

    p <- ggplot(sp_df, aes(x = Disease, y = Abundance, fill = Disease)) +

      # 半透明均值柱 | Semi-transparent mean bar
      stat_summary(fun = mean, geom = "bar",
                   color = "black", linewidth = 0.35, alpha = 0.65, width = 0.75) +

      # 无帽误差棒 (SEM) | Capless SEM error bars
      stat_summary(fun.data = mean_se, geom = "errorbar",
                   width = 0, linewidth = 0.5, color = "black") +

      # firebrick 实线轨迹 | Solid firebrick mean trajectory line
      stat_summary(fun = mean, geom = "line", aes(group = 1),
                   color = "firebrick", linewidth = 0.7, linetype = "solid") +

      # 空心白圆 + 红边 | Open white circles with firebrick stroke
      stat_summary(fun = mean, geom = "point",
                   size = 2.2, color = "firebrick", shape = 21, fill = "white", stroke = 1.1) +

      # 三 BMI 组并排，固定 Y 轴 | 3-BMI side-by-side, FIXED Y-axis
      facet_wrap(~ BMI_Group, scales = "fixed",
                 labeller = labeller(BMI_Group = bmi_labs), ncol = 3) +

      scale_fill_manual(values = disease_colors) +

      # Y 轴向上延展 15% 防裁切 | Expand Y to avoid clipping error bars
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +

      # 菌名为标题 (斜体加粗) | Species name as bold-italic title
      labs(title = sp, x = NULL, y = "Rel. Abundance (%)") +

      # CNS 极简主题 | CNS minimalist theme
      theme_classic(base_size = 10) +
      theme(
        plot.title       = element_text(face = "bold.italic", size = 10.5,
                                        hjust = 0.5, margin = margin(b = 8)),
        axis.line        = element_line(color = "black", linewidth = 0.4),
        axis.ticks       = element_line(color = "black", linewidth = 0.4),
        axis.text.x      = element_text(color = "black", size = 7.5,
                                        angle = 45, hjust = 1, vjust = 1),
        axis.text.y      = element_text(color = "black", size = 7.5),
        axis.title.y     = element_text(color = "black", face = "bold", size = 8.5),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 8.5, margin = margin(b = 4)),
        legend.position  = "none",
        panel.spacing.x  = unit(0.2, "lines")
      )

    plot_list[[sp]] <- p
  }

  # ---- 6. Patchwork 拼接 | Assemble with patchwork ----
  final_plot <- wrap_plots(plot_list, ncol = ncol_plots) +
    plot_annotation(
      title    = page_title,
      subtitle = paste0("Mean \u00b1 SEM  |  Batch-Corrected Relative Abundance (%)  |  Y-axis fixed across BMI groups"),
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(size = 9,  color = "grey40", hjust = 0.5)
      )
    )

  # ---- 7. 动态尺寸 & 保存 | Dynamic sizing and save ----
  n_rows   <- ceiling(n_species / ncol_plots)
  fig_w    <- ncol_plots * subplot_w + 0.5          # +0.5 for margins
  fig_h    <- max(3, n_rows * subplot_h + 1.5)      # +1.5 for title/subtitle
  out_path <- file.path(output_dir, filename)

  ggsave(out_path, final_plot, width = fig_w, height = fig_h,
         limitsize = FALSE, dpi = 300, device = grDevices::cairo_pdf)
  message(sprintf("   >>> Saved: %s  (%.1f x %.1f in, %d rows x %d cols)",
                  basename(out_path), fig_w, fig_h, n_rows, ncol_plots))
  return(invisible(NULL))
}

# ==============================================================================
# 7. 生成输出 | Generate PDF Output
# ==============================================================================
message("=== 6. 生成 PDF (CNS style) | Generating PDFs ===")

# 升高的 Ob-specific 菌 (方向一致) | Ob-specific enriched taxa (direction matched)
plot_bmi_trajectory_barplots(
  df = plot_df, status_type = "Enriched",
  filename   = "Ob_Specific_Trajectory_3BMI_Enriched_Matched.pdf",
  output_dir = output_dir
)

# 降低的菌 (方向一致) | Depleted taxa (direction matched)
plot_bmi_trajectory_barplots(
  df = plot_df, status_type = "Depleted",
  filename   = "Ob_Specific_Trajectory_3BMI_Depleted_Matched.pdf",
  output_dir = output_dir
)

message("✅ 完成! | Done!")

# ==============================================================================
# 8. 保存批次校正缓存，供 29l 等脚本复用 | Cache for reuse by 29l etc.
# ==============================================================================
cache_file <- file.path(cache_dir, "abd_corrected_meta_clean.rds")
saveRDS(list(abd_corrected = abd_corrected, meta_clean = meta_clean), cache_file)



# ==============================================================================
# 8. 保存批次校正缓存，供 29l 等脚本复用 | Cache for reuse by 29l etc.
# ==============================================================================
cache_file <- file.path(cache_dir, "abd_corrected_meta_clean.rds")
saveRDS(list(abd_corrected = abd_corrected, meta_clean = meta_clean), cache_file)
message(">>> 批次校正结果已缓存 | Cached: ", basename(cache_file))
