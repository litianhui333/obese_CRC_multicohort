suppressPackageStartupMessages({
  suppressWarnings(library(dplyr))
  suppressWarnings(library(tidyr))
  suppressWarnings(library(ggplot2))
  suppressWarnings(library(patchwork))
  suppressWarnings(library(readr))
})

base_dir <- "G:/My Drive/proj_ob_AA_CRC"
analysis_dir <- file.path(base_dir, "analyses/2026-05-30_da_trajectory_barplots_with_significance")
figures_dir <- file.path(analysis_dir, "figures")
tables_dir <- file.path(analysis_dir, "tables")
data_for_figures_dir <- file.path(analysis_dir, "data_for_figures")

input_dir <- file.path(base_dir, "03.results/04_differential_abundance")
heat_dir <- file.path(input_dir, "heatmap_final/BMI_Stratified_Disease_Heatmaps")
cache_file <- file.path(input_dir, "heatmap_final/cache/abd_corrected_meta_clean.rds")
res_file <- file.path(input_dir, "microbiome_analysis_all_results_BMI_Group1_0.7_20260123.rds")
robust_file <- file.path(input_dir, "heatmap_final/Robust_Analysis_v5_Pure/Robust_Stat_overall_CTRL_vs_CRC.rds")
robust_transition_files <- c(
  "CTRL_vs_AA" = file.path(input_dir, "heatmap_final/Robust_Analysis_v5_Pure/Robust_Stat_overall_CTRL_vs_AA.rds"),
  "AA_vs_CRC" = file.path(input_dir, "heatmap_final/Robust_Analysis_v5_Pure/Robust_Stat_overall_AA_vs_CRC.rds")
)

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_for_figures_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(cache_file)) stop("Missing cache: ", cache_file)
if (!file.exists(res_file)) stop("Missing MMUPHin result: ", res_file)
if (!file.exists(robust_file)) stop("Missing robust/rescue P-value matrix: ", robust_file)
if (any(!file.exists(robust_transition_files))) {
  stop("Missing robust/rescue transition P-value matrix: ", paste(robust_transition_files[!file.exists(robust_transition_files)], collapse = "; "))
}

format_sgb_names_local <- function(names_vec) {
  vapply(names_vec, function(x) {
    clean_name <- gsub(".*s__", "s__", x)
    parts <- strsplit(clean_name, "\\|")[[1]]
    sp_part <- gsub("^s__", "", parts[1])
    sgb_id <- if (length(parts) > 1) gsub("^t__", "", parts[2]) else ""

    if (grepl("_sp_", sp_part)) {
      label <- gsub("_sp_", " sp. ", sp_part)
      label <- gsub("_", " ", label)
    } else if (grepl("SGB", sp_part)) {
      label <- gsub("_", " ", sp_part)
    } else {
      sub_parts <- strsplit(sp_part, "_")[[1]]
      if (length(sub_parts) >= 2) {
        genus <- sub_parts[1]
        species <- paste(sub_parts[-1], collapse = " ")
        label <- paste0(substr(genus, 1, 1), ". ", species)
      } else {
        label <- gsub("_", " ", sp_part)
      }
    }

    label <- trimws(label)
    if (sgb_id != "" && !grepl(sgb_id, label, fixed = TRUE)) {
      label <- paste(label, sgb_id)
    }
    label
  }, character(1))
}

format_p_star <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

format_p_label <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p >= 0.05) return(NA_character_)
  if (p < 0.0001) return("<0.0001")
  if (p < 0.001) return(sprintf("%.4f", p))
  if (p < 0.01) return(sprintf("%.3f", p))
  sprintf("%.2f", p)
}

bmi_display_levels <- c("Nw", "Ow", "Ob")

format_bmi_display <- function(x) {
  dplyr::recode(as.character(x), "Normal" = "Nw", "Ow" = "Ow", "Ob" = "Ob")
}

calculate_sample_direction <- function(features, bmi_group, fa, meta) {
  bmi_mask <- meta$BMI_Group1 == bmi_group
  ctrl_idx <- which(bmi_mask & meta$Disease == "CTRL")
  crc_idx <- which(bmi_mask & meta$Disease == "CRC")

  result <- c()
  for (feat in features) {
    if (feat %in% rownames(fa)) {
      ctrl_vals <- as.numeric(fa[feat, ctrl_idx])
      crc_vals <- as.numeric(fa[feat, crc_idx])
      ctrl_pos <- ctrl_vals[!is.na(ctrl_vals) & ctrl_vals > 0]
      crc_pos <- crc_vals[!is.na(crc_vals) & crc_vals > 0]

      if (length(ctrl_pos) > 0 && length(crc_pos) > 0) {
        result[feat] <- ifelse(mean(crc_pos) > mean(ctrl_pos), "UP", "DOWN")
      }
    }
  }
  result
}

map_to_matrix <- function(targets, mat_rownames) {
  result <- list()
  for (target in targets) {
    if (target %in% mat_rownames) {
      result[[target]] <- target
      next
    }

    sgb <- regmatches(target, regexpr("SGB[0-9]+", target))
    if (length(sgb) > 0) {
      matched <- grep(sgb, mat_rownames, value = TRUE, fixed = TRUE)
      if (length(matched) >= 1) {
        result[[target]] <- matched[1]
        next
      }
    }

    species_key <- gsub("\\|t__.*", "", gsub(".*\\|s__", "s__", target))
    matched <- grep(species_key, mat_rownames, value = TRUE, fixed = TRUE)
    if (length(matched) >= 1) result[[target]] <- matched[1]
  }
  result
}

extract_p_table <- function(all_results, feature_universe, contrast_name = "CTRL_vs_CRC") {
  bmi_map <- c("Normal" = "Normal", "Ow" = "Ow", "Ob" = "Ob")
  bmi_display_map <- c("Normal" = "Nw", "Ow" = "Ow", "Ob" = "Ob")
  rows <- list()

  for (bmi in names(bmi_map)) {
    list_name <- paste0("disease_analysis_", bmi_map[[bmi]])
    mf <- tryCatch(
      all_results$overall$disease_analysis$analysis_results[[list_name]]$all_results[[contrast_name]]$meta_fits,
      error = function(e) NULL
    )
    if (is.null(mf)) next

    if (!"feature" %in% colnames(mf)) {
      mf <- tibble::rownames_to_column(mf, "feature")
    }

    mf <- mf %>%
      filter(feature %in% feature_universe) %>%
      transmute(
        feature,
        BMI_Group = bmi_display_map[[bmi]],
        contrast = contrast_name,
        coef_raw_mmuphin = coef,
        pval_raw_mmuphin = pval,
        qval_raw_mmuphin = qval.fdr
      )
    rows[[bmi]] <- mf
  }

  bind_rows(rows)
}

extract_robust_p_table <- function(robust_path, feature_universe, contrast_name = "CTRL_vs_CRC") {
  robust <- readRDS(robust_path)
  p_mat <- as.matrix(robust$matrix_data$pval)
  coef_mat <- as.matrix(robust$matrix_data$coef)
  bmi_cols <- c("Nw" = "Normal", "Ow" = "Overweight", "Ob" = "Obese")
  features <- intersect(feature_universe, rownames(p_mat))

  bind_rows(lapply(names(bmi_cols), function(bmi) {
    col_name <- bmi_cols[[bmi]]
    tibble(
      feature = features,
      BMI_Group = bmi,
      contrast = contrast_name,
      coef_mmuphin = coef_mat[features, col_name],
      pval_mmuphin = p_mat[features, col_name],
      qval_mmuphin = NA_real_,
      sig_star = vapply(p_mat[features, col_name], format_p_star, character(1)),
      sig_label = vapply(p_mat[features, col_name], format_p_label, character(1))
    )
  }))
}

all_analysis_results <- readRDS(res_file)
data_info <- all_analysis_results$overall$disease_analysis$data_info
fa <- data_info$filtered_data
meta <- data_info$metadata_filtered

normal_csv <- read.csv(file.path(heat_dir, "Feature_List_Normal_Refined.csv"), stringsAsFactors = FALSE)
ow_csv <- read.csv(file.path(heat_dir, "Feature_List_Ow_Refined.csv"), stringsAsFactors = FALSE)
ob_csv <- read.csv(file.path(heat_dir, "Feature_List_Ob_Refined.csv"), stringsAsFactors = FALSE)

add_direction_match <- function(df, bmi_group) {
  sample_dir <- calculate_sample_direction(df$feature, bmi_group, fa, meta)
  df$sample_dir <- sample_dir[df$feature]
  df$direction_match <- ifelse(
    (df$coef_crc > 0 & df$sample_dir == "UP") |
      (df$coef_crc < 0 & df$sample_dir == "DOWN"),
    "Match",
    "Mismatch"
  )
  df
}

normal_csv <- add_direction_match(normal_csv, "Normal")
ow_csv <- add_direction_match(ow_csv, "Ow")
ob_csv <- add_direction_match(ob_csv, "Ob")

selected_taxa <- bind_rows(
  normal_csv %>%
    filter(grepl("Nw-specific", overlap)) %>%
    mutate(Category = "Nw-specific", DirectionSet = if_else(coef_crc > 0, "Enriched", "Depleted")),
  ow_csv %>%
    filter(grepl("Ow-specific", overlap)) %>%
    mutate(Category = "Ow-specific", DirectionSet = if_else(coef_crc > 0, "Enriched", "Depleted")),
  ob_csv %>%
    filter(grepl("Ob-specific", overlap) | overlap == "Overlap with Ow") %>%
    mutate(Category = "Ob-specific", DirectionSet = if_else(coef_crc > 0, "Enriched", "Depleted")),
  normal_csv %>%
    filter(grepl("Overlap with both Ow and Ob", overlap)) %>%
    mutate(Category = "Triple-Overlap", DirectionSet = if_else(coef_crc > 0, "Enriched", "Depleted"))
) %>%
  distinct(feature, Category, DirectionSet, .keep_all = TRUE) %>%
  mutate(Species = format_sgb_names_local(feature))

write_tsv(selected_taxa, file.path(tables_dir, "selected_da_trajectory_taxa.tsv"))

cache <- readRDS(cache_file)
abd_corrected <- cache$abd_corrected
meta_clean <- cache$meta_clean
meta_clean$Disease <- factor(meta_clean$Disease, levels = c("CTRL", "AA", "CRC"))
meta_clean$BMI_Group1 <- factor(meta_clean$BMI_Group1, levels = c("Normal", "Ow", "Ob"))

mapping <- map_to_matrix(selected_taxa$feature, rownames(abd_corrected))
missing <- setdiff(selected_taxa$feature, names(mapping))
if (length(missing) > 0) {
  warning("Skipped taxa not found in corrected abundance matrix: ", paste(format_sgb_names_local(missing), collapse = "; "))
}

matched_rows <- unlist(mapping)
names(matched_rows) <- names(mapping)

plot_rows <- list()
for (feature in names(matched_rows)) {
  meta_row <- selected_taxa %>% filter(feature == !!feature) %>% slice(1)
  values <- as.numeric(abd_corrected[matched_rows[[feature]], ])

  tmp <- data.frame(
    SampleID = colnames(abd_corrected),
    RelAbundance = values,
    stringsAsFactors = FALSE
  ) %>%
    inner_join(
      meta_clean %>%
        tibble::rownames_to_column("SampleID") %>%
        select(SampleID, Disease, BMI_Group1),
      by = "SampleID"
    ) %>%
    mutate(
      Abundance = asin(sqrt(pmin(pmax(RelAbundance, 0), 100) / 100)),
      Feature = feature,
      MatrixFeature = matched_rows[[feature]],
      Species = meta_row$Species,
      Category = meta_row$Category,
      DirectionSet = meta_row$DirectionSet
    )
  plot_rows[[feature]] <- tmp
}

plot_df <- bind_rows(plot_rows) %>%
  mutate(
    Disease = factor(Disease, levels = c("CTRL", "AA", "CRC")),
    BMI_Group = factor(format_bmi_display(BMI_Group1), levels = bmi_display_levels),
    Species = factor(Species, levels = unique(selected_taxa$Species[selected_taxa$feature %in% names(mapping)]))
  ) %>%
  select(-BMI_Group1)

p_table_raw_crc <- extract_p_table(all_analysis_results, names(mapping), "CTRL_vs_CRC")
p_table_robust_crc <- extract_robust_p_table(robust_file, names(mapping), "CTRL_vs_CRC")
p_table_crc <- p_table_robust_crc %>%
  left_join(p_table_raw_crc, by = c("feature", "BMI_Group", "contrast")) %>%
  mutate(show_annotation = !is.na(pval_mmuphin) & pval_mmuphin < 0.05)

p_table_transition <- bind_rows(lapply(names(robust_transition_files), function(contrast_name) {
  extract_robust_p_table(robust_transition_files[[contrast_name]], names(mapping), contrast_name) %>%
    left_join(
      extract_p_table(all_analysis_results, names(mapping), contrast_name),
      by = c("feature", "BMI_Group", "contrast")
    ) %>%
    mutate(
      show_annotation = !is.na(pval_mmuphin) & pval_mmuphin < 0.05 &
        ((!is.na(pval_raw_mmuphin) & pval_raw_mmuphin < 0.01) |
           (!is.na(qval_raw_mmuphin) & qval_raw_mmuphin < 0.1))
    ) %>%
    filter(show_annotation)
}))

p_table <- bind_rows(p_table_crc, p_table_transition) %>%
  mutate(
    contrast = factor(contrast, levels = c("CTRL_vs_AA", "CTRL_vs_CRC", "AA_vs_CRC")),
    contrast_label = recode(as.character(contrast), "CTRL_vs_AA" = "CTRL-AA", "CTRL_vs_CRC" = "CTRL-CRC", "AA_vs_CRC" = "AA-CRC")
  )

summary_for_sig <- plot_df %>%
  group_by(Feature, Species, Category, DirectionSet, BMI_Group, Disease) %>%
  summarise(
    mean_abundance = mean(Abundance, na.rm = TRUE),
    sem_abundance = sd(Abundance, na.rm = TRUE) / sqrt(sum(!is.na(Abundance))),
    .groups = "drop"
  )

sig_df <- summary_for_sig %>%
  group_by(Feature, Species, Category, DirectionSet, BMI_Group) %>%
  summarise(
    y = max(mean_abundance + sem_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Feature, Species) %>%
  mutate(y_max_species = max(y, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(p_table, by = c("Feature" = "feature", "BMI_Group")) %>%
  arrange(Feature, BMI_Group, contrast) %>%
  group_by(Feature, Species, BMI_Group) %>%
  mutate(
    annotation_rank = row_number(),
    y_position = y + (0.16 + 0.13 * (annotation_rank - 1)) * y_max_species
  ) %>%
  ungroup()

write_tsv(plot_df, file.path(data_for_figures_dir, "da_trajectory_plot_data.tsv"))
write_tsv(sig_df, file.path(data_for_figures_dir, "da_trajectory_significance.tsv"))

disease_colors <- c("CTRL" = "#2E9FDF", "AA" = "#F39C12", "CRC" = "#E74C3C")
bmi_labs <- c("Nw" = "Nw", "Ow" = "Ow", "Ob" = "Ob")

wrap_title <- function(x, width = 34) paste(strwrap(x, width = width), collapse = "\n")

plot_one_set <- function(df, sig, category_filter, direction_filter, filename, sort_by_bmi) {
  sub_df <- df %>%
    filter(Category == category_filter, DirectionSet == direction_filter)
  sub_sig <- sig %>%
    filter(Category == category_filter, DirectionSet == direction_filter)
  sub_df$BMI_Group <- factor(as.character(sub_df$BMI_Group), levels = bmi_display_levels)
  sub_sig$BMI_Group <- factor(as.character(sub_sig$BMI_Group), levels = bmi_display_levels)

  if (nrow(sub_df) == 0) {
    message("Skip empty set: ", category_filter, " ", direction_filter)
    return(invisible(NULL))
  }

  sort_fn <- if (direction_filter == "Enriched") desc else identity
  species_order <- sub_df %>%
    filter(BMI_Group == sort_by_bmi, Disease == "CRC") %>%
    group_by(Species) %>%
    summarise(m = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(sort_fn(m)) %>%
    pull(Species) %>%
    as.character()
  species_order <- c(species_order, setdiff(unique(as.character(sub_df$Species)), species_order))
  species_order <- species_order[!is.na(species_order) & species_order != ""]

  sub_df$Species <- factor(as.character(sub_df$Species), levels = species_order)
  sub_sig$Species <- factor(as.character(sub_sig$Species), levels = species_order)

  plot_list <- list()
  for (species in species_order) {
    species_df <- sub_df %>% filter(Species == species)
    species_sig <- sub_sig %>%
      filter(Species == species, show_annotation, !is.na(sig_label)) %>%
      mutate(
        BMI_Group = factor(as.character(BMI_Group), levels = bmi_display_levels),
        x_start = case_when(
          contrast == "CTRL_vs_AA" ~ 1,
          contrast == "AA_vs_CRC" ~ 2,
          TRUE ~ 1
        ),
        x_end = case_when(
          contrast == "CTRL_vs_AA" ~ 2,
          contrast == "AA_vs_CRC" ~ 3,
          TRUE ~ 3
        ),
        x_mid = (x_start + x_end) / 2
      )
    if (nrow(species_df) == 0) next
    if (nrow(species_sig) == 0) {
      species_sig <- data.frame(
        BMI_Group = factor(character(0), levels = levels(sub_df$BMI_Group)),
        x_start = numeric(0),
        x_end = numeric(0),
        x_mid = numeric(0),
        y_position = numeric(0),
        sig_label = character(0),
        contrast = factor(character(0), levels = c("CTRL_vs_AA", "CTRL_vs_CRC", "AA_vs_CRC"))
      )
    }

    p <- ggplot(species_df, aes(x = Disease, y = Abundance, fill = Disease)) +
      stat_summary(fun = mean, geom = "bar", color = "black", linewidth = 0.35, alpha = 0.65, width = 0.75) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0, linewidth = 0.5, color = "black") +
      stat_summary(fun = mean, geom = "line", aes(group = 1), color = "firebrick", linewidth = 0.7) +
      stat_summary(fun = mean, geom = "point", size = 2.2, color = "firebrick", shape = 21, fill = "white", stroke = 1.1) +
      geom_segment(
        data = species_sig,
        aes(x = x_start, xend = x_end, y = y_position, yend = y_position),
        inherit.aes = FALSE,
        linewidth = 0.35,
        color = "black"
      ) +
      geom_text(
        data = species_sig,
        aes(x = x_mid, y = y_position, label = sig_label),
        inherit.aes = FALSE,
        vjust = -0.25,
        size = 2.5,
        color = "black"
      ) +
      facet_wrap(~ BMI_Group, scales = "fixed", labeller = labeller(BMI_Group = bmi_labs), ncol = 3) +
      scale_fill_manual(values = disease_colors) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.34))) +
      labs(title = wrap_title(species), x = NULL, y = "AST abundance") +
      theme_classic(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold.italic", size = 10.5, hjust = 0.5, margin = margin(b = 8)),
        axis.line = element_line(color = "black", linewidth = 0.4),
        axis.ticks = element_line(color = "black", linewidth = 0.4),
        axis.text.x = element_text(color = "black", size = 7.5, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black", size = 7.5),
        axis.title.y = element_text(color = "black", face = "bold", size = 8.5),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8.5, margin = margin(b = 4)),
        legend.position = "none",
        panel.spacing.x = unit(0.2, "lines")
      )

    plot_list[[species]] <- p
  }

  title_label <- switch(
    category_filter,
    "Nw-specific" = "Nw-specific signatures",
    "Ow-specific" = "Ow-specific signatures",
    "Ob-specific" = "Ob-focused signatures",
    "Triple-Overlap" = "Triple-overlap signatures",
    category_filter
  )
  direction_label <- if (direction_filter == "Enriched") "Enriched (+)" else "Depleted (-)"

  ncol_plots <- 3
  subplot_w <- 3.2
  subplot_h <- 2.25
  n_species <- length(plot_list)
  if (n_species == 0) {
    message("Skip empty plot list: ", category_filter, " ", direction_filter)
    return(invisible(NULL))
  }
  n_rows <- ceiling(n_species / ncol_plots)
  fig_w <- ncol_plots * subplot_w + 0.5
  fig_h <- max(3.2, n_rows * subplot_h + 1.5)

  final_plot <- wrap_plots(plot_list, ncol = ncol_plots) +
    plot_annotation(
      title = paste0(title_label, " - ", direction_label, " (n = ", n_species, ")"),
      subtitle = "Mean +/- SEM | AST-transformed batch-corrected abundance | P-values shown for significant CTRL-AA, CTRL-CRC and AA-CRC contrasts",
      theme = theme(
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5)
      )
    )

  out_base <- file.path(figures_dir, filename)
  ggsave(paste0(out_base, ".pdf"), final_plot, width = fig_w, height = fig_h, limitsize = FALSE, dpi = 300, device = grDevices::cairo_pdf)
  ggsave(paste0(out_base, ".png"), final_plot, width = fig_w, height = fig_h, limitsize = FALSE, dpi = 300)

  message("Saved: ", paste0(out_base, ".pdf"))
  invisible(final_plot)
}

plot_one_set(plot_df, sig_df, "Nw-specific", "Enriched", "Nw_Specific_Trajectory_3BMI_Enriched_Significance", "Nw")
plot_one_set(plot_df, sig_df, "Nw-specific", "Depleted", "Nw_Specific_Trajectory_3BMI_Depleted_Significance", "Nw")
plot_one_set(plot_df, sig_df, "Ow-specific", "Enriched", "Ow_Specific_Trajectory_3BMI_Enriched_Significance", "Ow")
plot_one_set(plot_df, sig_df, "Ow-specific", "Depleted", "Ow_Specific_Trajectory_3BMI_Depleted_Significance", "Ow")
plot_one_set(plot_df, sig_df, "Ob-specific", "Enriched", "Ob_Specific_Trajectory_3BMI_Enriched_Significance", "Ob")
plot_one_set(plot_df, sig_df, "Ob-specific", "Depleted", "Ob_Specific_Trajectory_3BMI_Depleted_Significance", "Ob")
plot_one_set(plot_df, sig_df, "Triple-Overlap", "Enriched", "Triple_Overlap_Trajectory_3BMI_Enriched_Significance", "Ob")
plot_one_set(plot_df, sig_df, "Triple-Overlap", "Depleted", "Triple_Overlap_Trajectory_3BMI_Depleted_Significance", "Ob")

focused_features <- c(
  "s__Bacteroides_caccae|t__SGB1877",
  "s__Ruminococcus_torques|t__SGB4608",
  "s__Hungatella_hathewayi|t__SGB4742",
  "s__Anaerotruncus_colihominis|t__SGB14963"
)

plot_focused_taxa <- function(features, filename) {
  focused_mapping <- map_to_matrix(features, rownames(abd_corrected))
  if (length(focused_mapping) == 0) {
    warning("No focused taxa found in corrected abundance matrix.")
    return(invisible(NULL))
  }

  focused_rows <- list()
  for (feature in names(focused_mapping)) {
    values <- as.numeric(abd_corrected[focused_mapping[[feature]], ])
    focused_rows[[feature]] <- data.frame(
      SampleID = colnames(abd_corrected),
      RelAbundance = values,
      stringsAsFactors = FALSE
    ) %>%
      inner_join(
        meta_clean %>%
          tibble::rownames_to_column("SampleID") %>%
          select(SampleID, Disease, BMI_Group1),
        by = "SampleID"
      ) %>%
      mutate(
        Abundance = asin(sqrt(pmin(pmax(RelAbundance, 0), 100) / 100)),
        Feature = feature,
        MatrixFeature = focused_mapping[[feature]],
        Species = format_sgb_names_local(feature),
        BMI_Group = factor(format_bmi_display(BMI_Group1), levels = bmi_display_levels)
      ) %>%
      select(-BMI_Group1)
  }

  focused_df <- bind_rows(focused_rows) %>%
    mutate(
      Disease = factor(Disease, levels = c("CTRL", "AA", "CRC")),
      Species = factor(Species, levels = format_sgb_names_local(features))
    )

  focused_raw_crc <- extract_p_table(all_analysis_results, names(focused_mapping), "CTRL_vs_CRC")
  focused_robust_crc <- extract_robust_p_table(robust_file, names(focused_mapping), "CTRL_vs_CRC")
  focused_crc <- focused_robust_crc %>%
    left_join(focused_raw_crc, by = c("feature", "BMI_Group", "contrast")) %>%
    mutate(
      show_annotation = !is.na(pval_mmuphin) & pval_mmuphin < 0.05 &
        ((!is.na(pval_raw_mmuphin) & pval_raw_mmuphin < 0.01) |
           (!is.na(qval_raw_mmuphin) & qval_raw_mmuphin < 0.1)),
      sig_label = if_else(show_annotation, vapply(pval_mmuphin, format_p_label, character(1)), NA_character_)
    )

  focused_ctrl_aa <- extract_p_table(all_analysis_results, names(focused_mapping), "CTRL_vs_AA") %>%
    transmute(
      feature,
      BMI_Group,
      contrast,
      coef_mmuphin = coef_raw_mmuphin,
      pval_mmuphin = pval_raw_mmuphin,
      qval_mmuphin = qval_raw_mmuphin,
      coef_raw_mmuphin,
      pval_raw_mmuphin,
      qval_raw_mmuphin,
      sig_star = vapply(pval_raw_mmuphin, format_p_star, character(1)),
      sig_label = vapply(pval_raw_mmuphin, format_p_label, character(1)),
      show_annotation = !is.na(pval_raw_mmuphin) & pval_raw_mmuphin < 0.05
    )

  focused_sig_base <- focused_df %>%
    group_by(Feature, Species, BMI_Group, Disease) %>%
    summarise(
      mean_abundance = mean(Abundance, na.rm = TRUE),
      sem_abundance = sd(Abundance, na.rm = TRUE) / sqrt(sum(!is.na(Abundance))),
      .groups = "drop"
    ) %>%
    group_by(Feature, Species, BMI_Group) %>%
    summarise(y = max(mean_abundance + sem_abundance, na.rm = TRUE), .groups = "drop") %>%
    group_by(Feature, Species) %>%
    mutate(y_max_species = max(y, na.rm = TRUE)) %>%
    ungroup()

  focused_sig <- bind_rows(focused_ctrl_aa, focused_crc) %>%
    mutate(
      Species = format_sgb_names_local(feature),
      contrast = factor(contrast, levels = c("CTRL_vs_AA", "CTRL_vs_CRC")),
      contrast_label = recode(as.character(contrast), "CTRL_vs_AA" = "CTRL-AA", "CTRL_vs_CRC" = "CTRL-CRC")
    ) %>%
    left_join(focused_sig_base, by = c("feature" = "Feature", "Species", "BMI_Group")) %>%
    arrange(feature, BMI_Group, contrast) %>%
    group_by(feature, Species, BMI_Group) %>%
    mutate(
      annotation_rank = row_number(),
      y_position = y + (0.16 + 0.13 * (annotation_rank - 1)) * y_max_species
    ) %>%
    ungroup()

  write_tsv(focused_df, file.path(data_for_figures_dir, paste0(filename, "_plot_data.tsv")))
  write_tsv(focused_sig, file.path(data_for_figures_dir, paste0(filename, "_significance.tsv")))

  plot_list <- list()
  for (species in levels(focused_df$Species)) {
    species_df <- focused_df %>% filter(Species == species)
    species_sig <- focused_sig %>%
      filter(Species == species, show_annotation, !is.na(sig_label), !is.na(y_position)) %>%
      mutate(
        BMI_Group = factor(as.character(BMI_Group), levels = bmi_display_levels),
        x_start = if_else(contrast == "CTRL_vs_AA", 1, 1),
        x_end = if_else(contrast == "CTRL_vs_AA", 2, 3),
        x_mid = (x_start + x_end) / 2
      )
    if (nrow(species_sig) == 0) {
      species_sig <- data.frame(
        BMI_Group = factor(character(0), levels = bmi_display_levels),
        x_start = numeric(0),
        x_end = numeric(0),
        x_mid = numeric(0),
        y_position = numeric(0),
        sig_label = character(0)
      )
    }

    plot_list[[species]] <- ggplot(species_df, aes(x = Disease, y = Abundance, fill = Disease)) +
      stat_summary(fun = mean, geom = "bar", color = "black", linewidth = 0.35, alpha = 0.65, width = 0.75) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = 0, linewidth = 0.5, color = "black") +
      stat_summary(fun = mean, geom = "line", aes(group = 1), color = "firebrick", linewidth = 0.7) +
      stat_summary(fun = mean, geom = "point", size = 2.2, color = "firebrick", shape = 21, fill = "white", stroke = 1.1) +
      geom_segment(
        data = species_sig,
        aes(x = x_start, xend = x_end, y = y_position, yend = y_position),
        inherit.aes = FALSE,
        linewidth = 0.35,
        color = "black"
      ) +
      geom_text(
        data = species_sig,
        aes(x = x_mid, y = y_position, label = sig_label),
        inherit.aes = FALSE,
        vjust = -0.25,
        size = 2.5,
        color = "black"
      ) +
      facet_wrap(~ BMI_Group, scales = "fixed", labeller = labeller(BMI_Group = bmi_labs), ncol = 3) +
      scale_fill_manual(values = disease_colors) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.34))) +
      labs(title = wrap_title(species), x = NULL, y = "AST abundance") +
      theme_classic(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold.italic", size = 10.5, hjust = 0.5, margin = margin(b = 8)),
        axis.line = element_line(color = "black", linewidth = 0.4),
        axis.ticks = element_line(color = "black", linewidth = 0.4),
        axis.text.x = element_text(color = "black", size = 7.5, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black", size = 7.5),
        axis.title.y = element_text(color = "black", face = "bold", size = 8.5),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 8.5, margin = margin(b = 4)),
        legend.position = "none",
        panel.spacing.x = unit(0.2, "lines")
      )
  }

  final_plot <- wrap_plots(plot_list, ncol = 2) +
    plot_annotation(
      title = "Focused BMI-stratified trajectories",
      subtitle = "CTRL-CRC uses strict raw P/FDR plus robust/rescue P; CTRL-AA shown when P < 0.05",
      theme = theme(
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5)
      )
    )
  out_base <- file.path(figures_dir, filename)
  ggsave(paste0(out_base, ".pdf"), final_plot, width = 7.2, height = 7.0, limitsize = FALSE, dpi = 300, device = grDevices::cairo_pdf)
  ggsave(paste0(out_base, ".png"), final_plot, width = 7.2, height = 7.0, limitsize = FALSE, dpi = 300)
  message("Saved: ", paste0(out_base, ".pdf"))
  invisible(final_plot)
}

plot_focused_taxa(focused_features, "Focused_B_caccae_R_torques_3BMI_CTRL_CRC_CTRL_AA_Significance")

message("Done.")
