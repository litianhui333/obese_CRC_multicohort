rm(list = ls())

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
})

# ==============================================================================
# Script Name: 64zs_DGCA_DA_Seed_Subnetwork_Obese_Union_Heatmaps.R
# Function:
#   Extract obese-defined DGCA networks among:
#     1. Ob significant DA taxa
#     2. Union of significant DA taxa from Normal/Ow/Ob
#   and draw 3x3 BMI-by-stage heatmaps without Loss/Gain text labels.
#
# Notes:
#   - Shared BMI trajectories are removed, matching the current 64zp emphasis.
#   - Only DA-DA edges are retained; non-DA nodes are not shown.
#   - Ob significant taxa are ordered first and labeled in both seed sets.
# ==============================================================================

base_dir <- "G:/My Drive/proj_ob_AA_CRC"
dgca_root <- file.path(base_dir, "03.results/06_network_analysis/DGCA_Results")
heatmap_dir <- file.path(
  base_dir,
  "03.results/04_differential_abundance/heatmap_final/BMI_Stratified_Disease_Heatmaps"
)

trajectory_file <- file.path(
  dgca_root,
  "Edge_State_Trajectory_Across_BMI_rho02_p05/Edge_Trajectory_Wide_Across_BMI_rho02_p05.csv"
)
ctrl_crc_ob_file <- file.path(
  dgca_root,
  "CTRL_Ob_vs_CRC_Ob/table/DGCA_all_CTRL_Ob_vs_CRC_Ob.csv"
)

out_dir <- file.path(dgca_root, "DA_only_subnetwork_Obese_Union_Heatmaps")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(trajectory_file), file.exists(ctrl_crc_ob_file))

da_files <- c(
  Normal = "Feature_List_Normal_Refined.csv",
  Ow = "Feature_List_Ow_Refined.csv",
  Ob = "Feature_List_Ob_Refined.csv"
)
missing_da <- da_files[!file.exists(file.path(heatmap_dir, da_files))]
if (length(missing_da) > 0) {
  stop("Missing DA feature files: ", paste(missing_da, collapse = "; "))
}

gain_patterns <- c("00+", "0++", "00-", "0--")
loss_patterns <- c("+00", "++0", "-00", "--0")
bmi_levels <- c("Ob", "Ow", "Normal")
bmi_labels <- c("Ob" = "Ob", "Ow" = "Ow", "Normal" = "Nw")
stage_levels <- c("CTRL", "AA", "CRC")
p_diff_cutoff <- 0.001

stage_col_map <- tribble(
  ~BMI,      ~Stage,  ~cor_col,
  "Normal",  "CTRL",  "cor_round__CTRL__Normal",
  "Normal",  "AA",    "cor_round__AA__Normal",
  "Normal",  "CRC",   "cor_round__CRC__Normal",
  "Ow",      "CTRL",  "cor_round__CTRL__Ow",
  "Ow",      "AA",    "cor_round__AA__Ow",
  "Ow",      "CRC",   "cor_round__CRC__Ow",
  "Ob",      "CTRL",  "cor_round__CTRL__Ob",
  "Ob",      "AA",    "cor_round__AA__Ob",
  "Ob",      "CRC",   "cor_round__CRC__Ob"
)

format_sgb_name <- function(x) {
  x <- str_replace(x, "^s__", "")
  x <- str_replace_all(x, "_", " ")
  x
}

make_edge_id <- function(g1, g2) {
  paste(pmin(g1, g2), pmax(g1, g2), sep = "___")
}

filter_shared_bmi_trajectory <- function(df) {
  df %>%
    filter(!(trajectory__Normal == trajectory__Ow & trajectory__Ow == trajectory__Ob))
}

safe_ggsave <- function(filename, plot, ...) {
  tryCatch(
    {
      ggsave(filename, plot, ...)
      filename
    },
    error = function(e) {
      message("Could not write figure, possibly open/locked: ", filename)
      message("  ", conditionMessage(e))
      NA_character_
    }
  )
}

traj_df <- read_csv(trajectory_file, show_col_types = FALSE)
ctrl_crc_ob_df <- read_csv(ctrl_crc_ob_file, show_col_types = FALSE) %>%
  transmute(
    edge_id = make_edge_id(Gene1, Gene2),
    pValDiff,
    pValDiff_adj,
    cor_ctrl_ob = CTRL_Ob_cor,
    cor_crc_ob = CRC_Ob_cor
  ) %>%
  distinct(edge_id, .keep_all = TRUE)

da_long <- bind_rows(lapply(names(da_files), function(bmi_i) {
  read_csv(file.path(heatmap_dir, da_files[[bmi_i]]), show_col_types = FALSE) %>%
    mutate(source_bmi = bmi_i)
})) %>%
  mutate(clean_name = if_else(is.na(clean_name) | clean_name == "", format_sgb_name(feature), clean_name))

ob_seed_tbl <- da_long %>%
  filter(source_bmi == "Ob") %>%
  distinct(feature, .keep_all = TRUE) %>%
  transmute(
    seed_taxon = feature,
    clean_name,
    seed_source = "Ob significant DA",
    crc_direction = if_else(coef_crc > 0, "CRC-enriched", "CRC-depleted"),
    ob_overlap_class = overlap,
    ob_coef_crc = coef_crc,
    ob_qval_crc = qval_crc
  )

union_seed_tbl <- da_long %>%
  group_by(feature) %>%
  summarise(
    clean_name = first(na.omit(clean_name)),
    seed_source = paste(sort(unique(source_bmi)), collapse = ";"),
    source_overlap = paste(sort(unique(paste(source_bmi, overlap, sep = ":"))), collapse = ";"),
    n_pos = sum(coef_crc > 0, na.rm = TRUE),
    n_neg = sum(coef_crc < 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    crc_direction = case_when(
      n_pos > 0 & n_neg == 0 ~ "CRC-enriched",
      n_neg > 0 & n_pos == 0 ~ "CRC-depleted",
      TRUE ~ "Mixed direction"
    ),
    is_ob_seed = feature %in% ob_seed_tbl$seed_taxon
  )

ob_seed_taxa <- ob_seed_tbl$seed_taxon
union_seed_taxa <- union_seed_tbl$feature

ob_label_lookup <- ob_seed_tbl %>%
  select(node = seed_taxon, clean_name, ob_overlap_class)

edge_sets <- bind_rows(
  traj_df %>%
    filter_shared_bmi_trajectory() %>%
    filter(trajectory__Ob %in% gain_patterns) %>%
    transmute(edge_id, Gene1, Gene2, edge_component = "Gain", heatmap_layer = "Lower"),
  traj_df %>%
    filter_shared_bmi_trajectory() %>%
    filter(trajectory__Ob %in% loss_patterns) %>%
    transmute(edge_id, Gene1, Gene2, edge_component = "Loss", heatmap_layer = "Upper"),
  traj_df %>%
    filter_shared_bmi_trajectory() %>%
    left_join(ctrl_crc_ob_df, by = "edge_id") %>%
    filter(
      trajectory__Ob %in% c("+++", "---"),
      !is.na(pValDiff),
      pValDiff < p_diff_cutoff,
      !is.na(cor_ctrl_ob),
      !is.na(cor_crc_ob),
      abs(cor_crc_ob) != abs(cor_ctrl_ob)
    ) %>%
    mutate(
      edge_component = if_else(abs(cor_crc_ob) > abs(cor_ctrl_ob), "Strengthened", "Weakened"),
      heatmap_layer = if_else(edge_component == "Strengthened", "Lower", "Upper")
    ) %>%
    transmute(edge_id, Gene1, Gene2, edge_component, heatmap_layer)
) %>%
  distinct(edge_id, heatmap_layer, .keep_all = TRUE)

make_da_only_heatmap <- function(seed_taxa, seed_table, set_name, set_title, network_change_label, edge_layer_filter) {
  seed_direction_lookup <- seed_table %>%
    mutate(node = if ("seed_taxon" %in% names(.)) seed_taxon else feature) %>%
    distinct(node, crc_direction)

  subnet_edges <- edge_sets %>%
    filter(
      Gene1 %in% seed_taxa,
      Gene2 %in% seed_taxa,
      heatmap_layer %in% edge_layer_filter
    )

  if (nrow(subnet_edges) == 0) {
    warning("No DA-only subnetwork edges found for ", set_name)
    return(invisible(NULL))
  }

  node_degree <- subnet_edges %>%
    select(Gene1, Gene2) %>%
    pivot_longer(cols = everything(), values_to = "node") %>%
    count(node, name = "degree")

  nodes <- node_degree %>%
    mutate(
      node_group = case_when(
        node %in% ob_seed_taxa ~ "Ob significant DA",
        node %in% seed_taxa ~ "Other DA seed",
        TRUE ~ "Other"
      ),
      group_rank = case_when(
        node_group == "Ob significant DA" ~ 1L,
        node_group == "Other DA seed" ~ 2L,
        TRUE ~ 3L
      )
    ) %>%
    left_join(ob_label_lookup, by = "node") %>%
    left_join(seed_direction_lookup, by = "node") %>%
    mutate(
      display_name = case_when(
        !is.na(clean_name) ~ clean_name,
        TRUE ~ format_sgb_name(node)
      ),
      crc_direction = factor(crc_direction, levels = c("CRC-enriched", "CRC-depleted")),
      direction_rank = case_when(
        crc_direction == "CRC-enriched" ~ 1L,
        crc_direction == "CRC-depleted" ~ 2L,
        TRUE ~ 3L
      )
    ) %>%
    arrange(direction_rank, group_rank, desc(degree), display_name) %>%
    mutate(idx = row_number())

  plot_df <- bind_rows(lapply(seq_len(nrow(stage_col_map)), function(i) {
    bmi_i <- stage_col_map$BMI[i]
    stage_i <- stage_col_map$Stage[i]
    cor_col_i <- stage_col_map$cor_col[i]

    subnet_edges %>%
      left_join(traj_df %>% select(edge_id, all_of(cor_col_i)), by = "edge_id") %>%
      left_join(nodes %>% select(node, i = idx), by = c("Gene1" = "node")) %>%
      left_join(nodes %>% select(node, j = idx), by = c("Gene2" = "node")) %>%
      transmute(
        seed_set = set_name,
        edge_id,
        Gene1,
        Gene2,
        edge_component,
        heatmap_layer,
        BMI = bmi_i,
        Stage = stage_i,
        rho = .data[[cor_col_i]],
        row_idx = pmax(i, j),
        col_idx = pmin(i, j)
      )
  })) %>%
    filter(!is.na(rho)) %>%
    mutate(
      BMI = factor(BMI, levels = bmi_levels),
      Stage = factor(Stage, levels = stage_levels)
    )

  node_ob_crc_shift <- plot_df %>%
    filter(Stage %in% c("CTRL", "CRC"), BMI %in% c("Normal", "Ob")) %>%
    select(edge_id, Gene1, Gene2, BMI, Stage, rho) %>%
    pivot_longer(c(Gene1, Gene2), values_to = "node", names_to = "node_side") %>%
    group_by(node, BMI, Stage) %>%
    summarise(mean_abs_rho = mean(abs(rho), na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = c(BMI, Stage), values_from = mean_abs_rho) %>%
    mutate(ob_crc_shift = `Ob_CRC` - `Normal_CRC`)

  nodes <- nodes %>%
    left_join(node_ob_crc_shift %>% select(node, ob_crc_shift), by = "node") %>%
    arrange(direction_rank, group_rank, desc(ob_crc_shift), desc(degree), display_name) %>%
    mutate(idx = row_number())

  plot_df <- plot_df %>%
    select(-row_idx, -col_idx) %>%
    left_join(nodes %>% select(node, i = idx), by = c("Gene1" = "node")) %>%
    left_join(nodes %>% select(node, j = idx), by = c("Gene2" = "node")) %>%
    mutate(
      row_idx = pmax(i, j),
      col_idx = pmin(i, j)
    )

  label_nodes <- nodes %>%
    filter(node_group == "Ob significant DA") %>%
    transmute(
      row_idx = idx,
      col_idx = -max(nodes$idx) * 0.02,
      label = display_name
    )

  direction_bar <- nodes %>%
    transmute(
      row_idx = idx,
      x = -max(nodes$idx) * 0.008,
      xend = -max(nodes$idx) * 0.008,
      y = idx - 0.45,
      yend = idx + 0.45,
      crc_direction
    )

  bottom_direction_bar <- nodes %>%
    transmute(
      col_idx = idx,
      x = idx - 0.45,
      xend = idx + 0.45,
      y = max(nodes$idx) + 0.85,
      yend = max(nodes$idx) + 0.85,
      crc_direction
    )

  p_base <- ggplot(plot_df, aes(x = col_idx, y = row_idx, fill = rho)) +
    geom_tile(width = 0.98, height = 0.98) +
    geom_segment(
      data = direction_bar,
      aes(x = x, xend = xend, y = y, yend = yend, color = crc_direction),
      inherit.aes = FALSE,
      linewidth = 2.2,
      lineend = "butt"
    ) +
    geom_segment(
      data = bottom_direction_bar,
      aes(x = x, xend = xend, y = y, yend = yend, color = crc_direction),
      inherit.aes = FALSE,
      linewidth = 2.2,
      lineend = "butt"
    ) +
    scale_fill_gradientn(
      colors = c("#2166AC", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#B2182B"),
      values = scales::rescale(c(-0.35, -0.1, 0, 0.1, 0.35)),
      limits = c(-0.35, 0.35),
      oob = scales::squish,
      name = "Correlation"
    ) +
    scale_color_manual(
      values = c("CRC-enriched" = "#B2182B", "CRC-depleted" = "#2166AC"),
      name = "Taxon status in CRC",
      drop = FALSE
    ) +
    scale_y_reverse() +
    coord_fixed(clip = "off") +
    facet_grid(BMI ~ Stage, labeller = labeller(BMI = bmi_labels), switch = "y") +
    labs(
      title = set_title,
      subtitle = paste0(network_change_label, "; all selected DA taxa are shown; red/blue bar marks CRC enriched/depleted taxa"),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      strip.text = element_text(face = "bold"),
      strip.text.y.left = element_text(angle = 0, face = "bold"),
      strip.placement = "outside",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
      legend.position = "right",
      plot.margin = margin(6, 6, 18, 18)
    )

  p_labeled <- p_base +
    geom_text(
      data = label_nodes,
      aes(x = col_idx, y = row_idx, label = label),
      inherit.aes = FALSE,
      hjust = 1,
      size = if_else(nrow(label_nodes) > 35, 1.5, 1.8),
      fontface = "italic",
      color = "grey10"
    ) +
    theme(plot.margin = margin(6, 6, 18, 72))

  width_i <- if_else(nrow(nodes) > 465, 10.8, 9.8)
  height_i <- if_else(nrow(nodes) > 465, 9.2, 8.6)

  safe_ggsave(file.path(out_dir, paste0(set_name, "_DAonly_3x3_RhoHeatmap.pdf")), p_labeled, width = width_i, height = height_i)
  safe_ggsave(file.path(out_dir, paste0(set_name, "_DAonly_3x3_RhoHeatmap.png")), p_labeled, width = width_i, height = height_i, dpi = 300)
  safe_ggsave(file.path(out_dir, paste0(set_name, "_NoObTaxaLabels_DAonly_3x3_RhoHeatmap.pdf")), p_base, width = width_i, height = height_i)
  safe_ggsave(file.path(out_dir, paste0(set_name, "_NoObTaxaLabels_DAonly_3x3_RhoHeatmap.png")), p_base, width = width_i, height = height_i, dpi = 300)

  write_csv(seed_table, file.path(out_dir, paste0(set_name, "_SeedTaxa.csv")))
  write_csv(nodes, file.path(out_dir, paste0(set_name, "_Nodes.csv")))
  write_csv(subnet_edges, file.path(out_dir, paste0(set_name, "_DAonly_Subnetwork_Edges.csv")))
  write_csv(plot_df, file.path(out_dir, paste0(set_name, "_Heatmap_Data.csv")))

  trend_summary <- plot_df %>%
    mutate(abs_rho = abs(rho)) %>%
    group_by(seed_set, BMI, Stage) %>%
    summarise(
      n_edges = n_distinct(edge_id),
      median_abs_rho = median(abs_rho, na.rm = TRUE),
      mean_abs_rho = mean(abs_rho, na.rm = TRUE),
      q25_abs_rho = quantile(abs_rho, 0.25, na.rm = TRUE),
      q75_abs_rho = quantile(abs_rho, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  write_csv(trend_summary, file.path(out_dir, paste0(set_name, "_MedianAbsRho_Trend_By_BMI_Stage.csv")))

  p_trend <- ggplot(trend_summary, aes(x = Stage, y = median_abs_rho, color = BMI, group = BMI)) +
    geom_line(linewidth = 0.75) +
    geom_point(size = 2.2) +
    scale_color_manual(
      values = c("Ob" = "#9970AB", "Ow" = "#FF8C00", "Normal" = "#5AAE61"),
      breaks = bmi_levels,
      labels = bmi_labels
    ) +
    labs(
      title = paste0(set_title, ": median network strength"),
      subtitle = network_change_label,
      x = NULL,
      y = expression("Median |" * rho * "|")
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  safe_ggsave(file.path(out_dir, paste0(set_name, "_MedianAbsRho_Trend_By_BMI_Stage.pdf")), p_trend, width = 4.8, height = 3.4)
  safe_ggsave(file.path(out_dir, paste0(set_name, "_MedianAbsRho_Trend_By_BMI_Stage.png")), p_trend, width = 4.8, height = 3.4, dpi = 300)

  tibble(
    seed_set = set_name,
    network_change = network_change_label,
    n_seed_taxa = length(unique(seed_taxa)),
    n_ob_seed_taxa = length(intersect(unique(seed_taxa), ob_seed_taxa)),
    n_edges = nrow(subnet_edges),
    n_nodes = nrow(nodes),
    n_labeled_ob_taxa = nrow(label_nodes),
    n_lower_edges = sum(subnet_edges$heatmap_layer == "Lower"),
    n_upper_edges = sum(subnet_edges$heatmap_layer == "Upper")
  )
}

summary_tbl <- bind_rows(
  make_da_only_heatmap(
    seed_taxa = ob_seed_tbl$seed_taxon,
    seed_table = ob_seed_tbl,
    set_name = "Ob_Significant_DA_Network_Gain_Strengthened",
    set_title = "Obese DA Taxa Network Gain / Strengthened",
    network_change_label = "Lower network: gain or strengthened edges",
    edge_layer_filter = "Lower"
  ),
  make_da_only_heatmap(
    seed_taxa = ob_seed_tbl$seed_taxon,
    seed_table = ob_seed_tbl,
    set_name = "Ob_Significant_DA_Network_Loss_Weakened",
    set_title = "Obese DA Taxa Network Loss / Weakened",
    network_change_label = "Upper network: loss or weakened edges",
    edge_layer_filter = "Upper"
  ),
  make_da_only_heatmap(
    seed_taxa = union_seed_tbl$feature,
    seed_table = union_seed_tbl,
    set_name = "Union_All_BMI_DA_Network_Gain_Strengthened_Highlight_Ob",
    set_title = "Union DA Taxa Network Gain / Strengthened",
    network_change_label = "Lower network: gain or strengthened edges",
    edge_layer_filter = "Lower"
  ),
  make_da_only_heatmap(
    seed_taxa = union_seed_tbl$feature,
    seed_table = union_seed_tbl,
    set_name = "Union_All_BMI_DA_Network_Loss_Weakened_Highlight_Ob",
    set_title = "Union DA Taxa Network Loss / Weakened",
    network_change_label = "Upper network: loss or weakened edges",
    edge_layer_filter = "Upper"
  )
)

write_csv(summary_tbl, file.path(out_dir, "DAonly_Subnetwork_Heatmap_Summary.csv"))

message("\nSaved DA-only subnetwork heatmaps to: ", out_dir)
print(summary_tbl)
