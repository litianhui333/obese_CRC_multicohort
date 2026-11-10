# ==============================================================================
# Script Name: 64_DGCA_Network_Analysis_9Groups.R
# Function: 
#   1. Load smart filtered bacteria list (from 10_robust configuration).
#   2. Load batch-corrected abundance data from 29k cache.
#   3. Subset corrected data by the robust bacterial features.
#   4. Compute pairwise DGCA between 9 groups (BMI_Disease).
# ==============================================================================

rm(list = ls())

# Suppress warnings
suppressPackageStartupMessages({
  library(tidyverse)
  library(reshape2)
  library(broom)
  library(ggplot2)
  library(RColorBrewer)
  library(DGCA)
})

# Paths
base_dir <- "G:/My Drive/proj_ob_AA_CRC"
script_dir <- file.path(base_dir, "01.scripts")
data_dir <- file.path(base_dir, "02.processed_data/02_aligned_data")
sp_dir <- file.path(base_dir, "02.processed_data/03_metaphlan_split_data")

cache_dir <- file.path(base_dir, "03.results/04_differential_abundance/heatmap_final/cache")
out_dir <- file.path(base_dir, "03.results/06_network_analysis/DGCA_Results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load Utilities
source(file.path(script_dir, "DA_function.R"))

# ==============================================================================
# 1. Feature Filtering (Mimicking 10_robust_meta_network_execution.R)
# ==============================================================================
message("=== 1. Smart Data Filtering (extracting robust features) ===")

meta_file <- sort(list.files(data_dir, pattern = "metadata_aligned.*\\.rds", full.names = T), decreasing = T)[1]
sp_file <- sort(list.files(sp_dir, pattern = "metaphlan_strain_renormalized.*\\.rds", full.names = T), decreasing = T)[1]

metadata <- readRDS(meta_file)
metaphlan <- readRDS(sp_file) 
rownames(metadata) <- metadata$Sample.Source.ID.processed

config <- list(
  transform_method = "AST",
  analysis_method = "HS",
  BMI_Group_var = "BMI_Group1",
  covariates = c("Age", "Sex"),
  study_prevalence_cut = 0.1,
  abundance_threshold = 0.0001,
  min_studies_ratio = 0.7,   # 10_robust uses 0.7 
  pval_type = "both",
  sig_criterion = "or",
  pval_threshold = 0.02,
  qval_threshold = 0.2,
  is_inhouse = FALSE,
  force_keep = NULL
)

filtering_result <- smart_data_filtering(
  metaphlan_data = t(metaphlan), 
  metadata = metadata,
  config = config
)

feature_data_filt <- filtering_result$filtered_data
# When processing transposed metaphlan, smart_data_filtering outputs it with transposed logic
# DA_function.R automatically formats filtered_data with rows as features
filt_features <- rownames(feature_data_filt)

# ==============================================================================
# 2. Corrected Data Sourcing (from 29k_Obese_Specific_Trajectory_Barplots.R)
# ==============================================================================
message("\n=== 2. Load Corrected Abundance Cache ===")

cache_file <- file.path(cache_dir, "abd_corrected_meta_clean.rds")
if (!file.exists(cache_file)) {
    stop("Cache file not found, please run 29k_Obese_Specific_Trajectory_Barplots.R first.")
}

cache_res <- readRDS(cache_file)
abd_corrected <- cache_res$abd_corrected  # Rows = Features, Cols = Samples
meta_clean <- cache_res$meta_clean

intersect_features <- intersect(filt_features, rownames(abd_corrected))
message(sprintf("   >>> %d robust features intersected with corrected target dataset.", length(intersect_features)))

if (length(intersect_features) == 0) {
  stop("   >>> No intersecting features! Check your metaphlan input transpose format or cache consistency.")
}

# Isolate the corrected data of interest
data_mat <- abd_corrected[intersect_features, ] 

# ==============================================================================
# 3. Setup 9 Groups
# ==============================================================================
message("\n=== 3. Setup 9 Groups ===")
meta_clean <- meta_clean[!is.na(meta_clean$Disease) & !is.na(meta_clean$BMI_Group1), ]
meta_clean$Analysis_Group <- paste(meta_clean$Disease, meta_clean$BMI_Group1, sep = "_")

valid_samps <- intersect(rownames(meta_clean), colnames(data_mat))
meta_clean <- meta_clean[valid_samps, ]
data_mat <- data_mat[, valid_samps]

groups <- unique(meta_clean$Analysis_Group)

# Generate biologically meaningful biological comparisons
comparisons <- list()

# 1. Within different BMI groups, compare Disease
bmis <- c("Normal", "Ow", "Ob")
dis_pairs <- list(c("CTRL", "AA"), c("AA", "CRC"), c("CTRL", "CRC"))

for(b in bmis) {
    for(dp in dis_pairs) {
        g1 <- paste0(dp[1], "_", b)
        g2 <- paste0(dp[2], "_", b)
        if(g1 %in% groups && g2 %in% groups) {
            comparisons[[length(comparisons) + 1]] <- c(g1, g2)
        }
    }
}

# 2. Within different Disease groups, compare BMI
diseases <- c("CTRL", "AA", "CRC")
bmi_pairs <- list(c("Normal", "Ow"), c("Ow", "Ob"), c("Normal", "Ob"))

for(d in diseases) {
    for(bp in bmi_pairs) {
        g1 <- paste0(d, "_", bp[1])
        g2 <- paste0(d, "_", bp[2])
        if(g1 %in% groups && g2 %in% groups) {
            comparisons[[length(comparisons) + 1]] <- c(g1, g2)
        }
    }
}

message(sprintf("   >>> Built %d disease-BMI groups, resulting in %d targeted pairwise DGCA combinations.", length(groups), length(comparisons)))

# ==============================================================================
# 4. DGCA Analysis Core Loop
# ==============================================================================
message("\n=== 4. Executing Pairwise Differential Correlation Comparison ===")

mycol <- rev(brewer.pal(9, "Spectral"))
names(mycol) <- c('< -0.2', '-0.2 ~ -0.1', '-0.1 ~ 0',
                  '0 ~ 0.1', '0.1 ~ 0.2', '0.2 ~ 0.3',
                  '0.3 ~ 0.4', '0.4 ~ 0.5', '> 0.5')

save_or_not <- TRUE

for(i in seq_along(comparisons)) {
  comp <- comparisons[[i]]
  comp1 <- comp[1]
  comp2 <- comp[2]
  group_name <- paste0(comp1, "_vs_", comp2)
  
  message(sprintf("\n>>> [%d/%d] Testing: %s", i, length(comparisons), group_name))
  
  # A. Filter subset
  c_samps <- rownames(meta_clean)[meta_clean$Analysis_Group %in% comp]
  c_meta <- meta_clean[c_samps, ]
  c_mat <- data_mat[, c_samps]
  
  # Check if there are enough samples per group
  sz <- table(c_meta$Analysis_Group)
  if(length(sz) < 2 || any(sz < 5)) {
    message("   ... Skipping: insufficient samples in one or both groups (< 5).")
    next
  }
  
  # Remove variance-zero features inside this subset
  keep_f <- apply(c_mat, 1, var) > 0
  if(sum(keep_f) < 5) {
     message("   ... Skipping: Insufficient variable features.")
     next
  }
  c_mat <- c_mat[keep_f, ]
  
  # B. Build comparative design matrix
  # Samples x Conditions
  stage_design_mat <- t(sapply(c_meta$Analysis_Group, function(x) {
    if(x == comp1) return(c(1, 0)) else return(c(0, 1))
  }))
  colnames(stage_design_mat) <- comp
  rownames(stage_design_mat) <- c_samps
  
  # C. Run Differential gene correlation logic
  set.seed(2023)
  tryCatch({
    ddcor_res <- ddcorAll(
      inputMat = c_mat, 
      corrType = 'spearman',
      design = stage_design_mat,
      compare = comp
    )
    
    # Evaluate significance
    all_summary <- ddcor_res
    all_summary$abs_zScoreDiff <- abs(all_summary$zScoreDiff)
    den_fil_df <- all_summary[!is.na(all_summary$empPVals) & all_summary$empPVals < 0.05, ]
    
    if (nrow(den_fil_df) == 0) {
      message("   ... No significant (empPVals < 0.05) robust correlations detected. Outputting base stats, skipping plots.")
    }
    
    if (save_or_not) {
      d_dir <- file.path(out_dir, group_name)
      if(!dir.exists(file.path(d_dir, "fig"))) dir.create(file.path(d_dir, "fig"), recursive = TRUE)
      if(!dir.exists(file.path(d_dir, "table"))) dir.create(file.path(d_dir, "table"), recursive = TRUE)
      
      write.csv(all_summary, file = file.path(d_dir, "table", paste0("DGCA_all_", group_name, ".csv")), row.names = FALSE)
      
      if (nrow(den_fil_df) > 0) {
          write.csv(den_fil_df, file = file.path(d_dir, "table", paste0("DGCA_sig_", group_name, ".csv")), row.names = FALSE)
          
          # D. Plotting
          # Assign unified type
          den_fil_df$type <- "Bacteria"
          
          p <- ggplot(data=den_fil_df, aes(x=zScoreDiff, group=type, fill=type)) +
            geom_density(adjust=1.5, position="stack") +
            theme_classic() + scale_fill_brewer(palette = "Set1")
            
          # Format classes
          bar_df2 <- data.frame(Classes = factor(den_fil_df$Classes), type = factor(den_fil_df$type))
          # Since there's only 1 type, we handle edge cases appropriately
          bar_dcast = reshape2::dcast(bar_df2, type ~ Classes, fun.aggregate = length)
          bar_melt = reshape2::melt(bar_dcast, id.vars = "type", variable.name = 'Classes')
          bar_melt <- bar_melt %>% dplyr::group_by(type) %>% dplyr::mutate(percentage = value/sum(value)*100)
          
          write.csv(bar_melt, file = file.path(d_dir, "table", paste0("DGCA_prop_", group_name, ".csv")), row.names = FALSE)
          
          g_bar <- ggplot(bar_melt, aes(x = type, y = percentage, fill = Classes)) + 
            geom_col(position = "stack") +
            scale_fill_brewer(palette = "Spectral") +
            theme(panel.grid.minor.y = element_line(colour = "grey"),
                  panel.background = element_rect(color = 'black', fill = 'transparent'),
                  axis.text.x = element_text(color="black",size=12,face ="bold"),
                  axis.text.y = element_text(color="black",size=12),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(color="black",size=14),
                  axis.line.y = element_line(color="black", linetype ="solid", linewidth=0.5),
                  axis.line.x = element_line(color="black", linetype ="solid", linewidth=0.5))
                  
          pdf(file = file.path(d_dir, "fig", paste0("DensityPlot-z_score_", group_name, ".pdf")), width = 6, height=7,
              device = grDevices::cairo_pdf)  # cairo_pdf handles fonts better
          print(p)
          dev.off()
          
          pdf(file = file.path(d_dir, "fig", paste0("BarPlot-Classes_", group_name, ".pdf")), width = 5, height=5,
              device = grDevices::cairo_pdf)
          print(g_bar)
          dev.off()
      }
    }
  }, error = function(e){
    message("   ... ERROR in pairwise test: ", e$message)
  })
}

message("\n=== All Comparisons Processed Successfully! ===")
