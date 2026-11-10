rm(list = ls())

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
})

# ==============================================================================
# Purpose:
#   Extract the source-based DGCA rewiring subnetwork defined as:
#     1) both endpoints are within the Union All-BMI DA taxa universe
#     2) at least one endpoint is an obese-DA taxon
#     3) CTRL_Ob vs CRC_Ob DGCA pValDiff < 0.01
#
# Notes:
#   - This is stricter than the earlier "obese-DA-involved" network because all
#     nodes remain inside the Union DA-only network.
#   - Output includes both the full source-based version and a contextualized
#     version excluding edges that show the same rewiring mode in normal-weight.
# ==============================================================================

base_dir <- "G:/My Drive/proj_ob_AA_CRC"
dgca_root <- file.path(base_dir, "03.results/06_network_analysis/DGCA_Results")
heatmap_dir <- file.path(
  base_dir,
  "03.results/04_differential_abundance/heatmap_final/BMI_Stratified_Disease_Heatmaps"
)
tax_file <- file.path(
  base_dir,
  "02.processed_data/03_metaphlan_split_data/strain_taxonomy_20260115.rds"
)

traj_path <- file.path(
  dgca_root,
  "Edge_State_Trajectory_Across_BMI_rho02_p05/Edge_Trajectory_Wide_Across_BMI_rho02_p05.csv"
)
ob_diff_path <- file.path(
  dgca_root,
  "CTRL_Ob_vs_CRC_Ob/table/DGCA_all_CTRL_Ob_vs_CRC_Ob.csv"
)
normal_diff_path <- file.path(
  dgca_root,
  "CTRL_Normal_vs_CRC_Normal/table/DGCA_all_CTRL_Normal_vs_CRC_Normal.csv"
)

feature_files <- c(
  Normal = file.path(heatmap_dir, "Feature_List_Normal_Refined.csv"),
  Ow = file.path(heatmap_dir, "Feature_List_Ow_Refined.csv"),
  Ob = file.path(heatmap_dir, "Feature_List_Ob_Refined.csv")
)

out_dir <- file.path(
  dgca_root,
  "Union_DAonly_ObeseDA_Participating_CTRL_CRC_Rewiring_p01"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  file.exists(traj_path),
  file.exists(ob_diff_path),
  file.exists(normal_diff_path),
  all(file.exists(feature_files))
)

dgca_p_cutoff <- 0.01
top_n_hubs <- 25

edge_class_order <- c(
  "Positive correlation lost",
  "Positive correlation gained",
  "Negative correlation lost",
  "Negative correlation gained",
  "Positive correlation strengthened",
  "Positive correlation weakened",
  "Negative correlation strengthened",
  "Negative correlation weakened"
)

class_short_map <- c(
  "Normal" = "Nw-spec",
  "Ow" = "Ow-spec",
  "Ob" = "Ob-spec",
  "Normal;Ow" = "Nw+Ow",
  "Normal;Ob" = "Nw+Ob",
  "Ob;Ow" = "Ow+Ob",
  "Normal;Ob;Ow" = "BMI-ind"
)

make_edge_id <- function(g1, g2) {
  paste(pmin(g1, g2), pmax(g1, g2), sep = "___")
}

clean_tax_label <- function(x) {
  x %>%
    str_replace("^[a-z]__", "") %>%
    str_replace_all("_", " ") %>%
    na_if("") %>%
    replace_na("Unclassified")
}

format_sgb_names <- function(names_vec) {
  sapply(names_vec, function(x) {
    clean_name <- gsub(".*s__", "s__", x)
    parts <- strsplit(clean_name, "\\|")[[1]]
    sp_part <- gsub("^s__", "", parts[1])
    sgb_id <- if (length(parts) > 1) gsub("^t__", "", parts[2]) else ""

    if (grepl("_sp_", sp_part)) {
      res <- gsub("_sp_", " sp. ", sp_part)
      res <- gsub("_", " ", res)
    } else if (grepl("SGB", sp_part)) {
      res <- gsub("_", " ", sp_part)
    } else {
      sub_parts <- strsplit(sp_part, "_")[[1]]
      if (length(sub_parts) >= 2) {
        genus <- sub_parts[1]
        species <- paste(sub_parts[-1], collapse = " ")
        genus_abbr <- paste0(substr(genus, 1, 1), ".")
        res <- paste(genus_abbr, species)
      } else {
        res <- gsub("_", " ", sp_part)
      }
    }

    if (sgb_id != "") paste(res, sgb_id) else res
  }, USE.NAMES = FALSE)
}

classify_edge <- function(trajectory, ctrl_cor, crc_cor) {
  case_when(
    trajectory %in% c("+00", "++0") ~ "Positive correlation lost",
    trajectory %in% c("00+", "0++") ~ "Positive correlation gained",
    trajectory %in% c("-00", "--0") ~ "Negative correlation lost",
    trajectory %in% c("00-", "0--") ~ "Negative correlation gained",
    trajectory == "+++" & abs(crc_cor) > abs(ctrl_cor) ~ "Positive correlation strengthened",
    trajectory == "+++" & abs(crc_cor) < abs(ctrl_cor) ~ "Positive correlation weakened",
    trajectory == "---" & abs(crc_cor) > abs(ctrl_cor) ~ "Negative correlation strengthened",
    trajectory == "---" & abs(crc_cor) < abs(ctrl_cor) ~ "Negative correlation weakened",
    TRUE ~ "Other"
  )
}

safe_file_label <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
}

feature_long <- bind_rows(lapply(names(feature_files), function(bmi_i) {
  read_csv(feature_files[[bmi_i]], show_col_types = FALSE) %>%
    mutate(source_bmi = bmi_i)
}))

union_meta <- feature_long %>%
  group_by(feature) %>%
  summarise(
    clean_name = dplyr::first(na.omit(clean_name)),
    source_bmi = paste(sort(unique(source_bmi)), collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(
    class_short = recode(source_bmi, !!!class_short_map),
    class_short = if_else(is.na(class_short), source_bmi, class_short)
  )

union_taxa <- union_meta$feature

feature_ob <- read_csv(feature_files[["Ob"]], show_col_types = FALSE) %>%
  rename(seed_taxon = feature, overlap_class = overlap) %>%
  distinct(seed_taxon, .keep_all = TRUE)

ob_da_set <- feature_ob$seed_taxon

display_ob_specific <- feature_ob %>%
  filter(overlap_class == "Ob-specific") %>%
  pull(seed_taxon) %>%
  union("s__Akkermansia_muciniphila|t__SGB9226")

taxonomy <- if (file.exists(tax_file)) {
  readRDS(tax_file) %>%
    transmute(
      node = clade_name_simple,
      phylum = clean_tax_label(phylum),
      family = clean_tax_label(family),
      genus = clean_tax_label(genus)
    ) %>%
    distinct(node, .keep_all = TRUE)
} else {
  tibble(node = character(), phylum = character(), family = character(), genus = character())
}

traj_df <- read_csv(traj_path, show_col_types = FALSE) %>%
  mutate(edge_id = make_edge_id(Gene1, Gene2)) %>%
  filter(Gene1 %in% union_taxa, Gene2 %in% union_taxa)

ob_diff <- read_csv(ob_diff_path, show_col_types = FALSE) %>%
  transmute(
    edge_id = make_edge_id(Gene1, Gene2),
    CTRL_Ob_cor,
    CRC_Ob_cor,
    zScoreDiff,
    pValDiff,
    pValDiff_adj,
    dgca_class = Classes
  ) %>%
  distinct(edge_id, .keep_all = TRUE)

normal_diff <- read_csv(normal_diff_path, show_col_types = FALSE) %>%
  transmute(
    edge_id = make_edge_id(Gene1, Gene2),
    CTRL_Normal_cor,
    CRC_Normal_cor,
    pValDiff_Normal = pValDiff,
    pValDiff_adj_Normal = pValDiff_adj,
    dgca_class_Normal = Classes
  ) %>%
  distinct(edge_id, .keep_all = TRUE)

rewiring_edges <- traj_df %>%
  mutate(
    has_obese_da_taxon = Gene1 %in% ob_da_set | Gene2 %in% ob_da_set,
    both_union_da = Gene1 %in% union_taxa & Gene2 %in% union_taxa,
    da_side = case_when(
      Gene1 %in% ob_da_set & Gene2 %in% ob_da_set ~ "Both taxa are obese-DA",
      Gene1 %in% ob_da_set | Gene2 %in% ob_da_set ~ "One taxon is obese-DA; one is other union-DA",
      TRUE ~ "No obese-DA taxon"
    ),
    ob_da_taxon = case_when(
      Gene1 %in% ob_da_set & Gene2 %in% ob_da_set ~ paste(Gene1, Gene2, sep = "; "),
      Gene1 %in% ob_da_set ~ Gene1,
      Gene2 %in% ob_da_set ~ Gene2,
      TRUE ~ NA_character_
    )
  ) %>%
  filter(both_union_da, has_obese_da_taxon) %>%
  left_join(ob_diff, by = "edge_id") %>%
  left_join(normal_diff, by = "edge_id") %>%
  filter(!is.na(pValDiff), pValDiff < dgca_p_cutoff) %>%
  mutate(
    edge_class = classify_edge(trajectory__Ob, CTRL_Ob_cor, CRC_Ob_cor),
    edge_class_normal = if_else(
      !is.na(pValDiff_Normal) & pValDiff_Normal < dgca_p_cutoff,
      classify_edge(trajectory__Normal, CTRL_Normal_cor, CRC_Normal_cor),
      "Not significant in Nw"
    ),
    same_mode_in_normal = as.character(edge_class) == edge_class_normal,
    delta_cor = CRC_Ob_cor - CTRL_Ob_cor,
    abs_delta_cor = abs(delta_cor),
    ctrl_abs_cor = abs(CTRL_Ob_cor),
    crc_abs_cor = abs(CRC_Ob_cor)
  ) %>%
  filter(edge_class %in% edge_class_order) %>%
  mutate(edge_class = factor(edge_class, levels = edge_class_order)) %>%
  arrange(edge_class, pValDiff, desc(abs_delta_cor))

contextualized_edges <- rewiring_edges %>%
  filter(!same_mode_in_normal)

trajectory_specific_edges <- rewiring_edges %>%
  filter(trajectory__Ob != trajectory__Normal)

node_group_label <- function(nodes) {
  case_when(
    nodes %in% display_ob_specific ~ "Ob-specific obese-DA taxon",
    nodes %in% ob_da_set ~ "Other obese-DA taxon",
    TRUE ~ "Other union-DA taxon"
  )
}

build_node_table <- function(edge_df, network_label) {
  node_degree <- bind_rows(
    edge_df %>% transmute(node = Gene1),
    edge_df %>% transmute(node = Gene2)
  ) %>%
    count(node, name = "degree", sort = TRUE)

  node_strength <- bind_rows(
    edge_df %>% transmute(node = Gene1, abs_cor = abs(CTRL_Ob_cor), abs_cor_crc = abs(CRC_Ob_cor)),
    edge_df %>% transmute(node = Gene2, abs_cor = abs(CTRL_Ob_cor), abs_cor_crc = abs(CRC_Ob_cor))
  ) %>%
    group_by(node) %>%
    summarise(
      strength_ctrl_ob = sum(abs_cor, na.rm = TRUE),
      strength_crc_ob = sum(abs_cor_crc, na.rm = TRUE),
      .groups = "drop"
    )

  node_degree %>%
    left_join(node_strength, by = "node") %>%
    left_join(
      feature_ob %>%
        transmute(
          node = seed_taxon,
          ob_clean_name = clean_name,
          overlap_class,
          coef_crc,
          pval_crc,
          qval_crc,
          coef_aa,
          pval_aa
        ),
      by = "node"
    ) %>%
    left_join(
      union_meta %>%
        transmute(
          node = feature,
          union_clean_name = clean_name,
          source_bmi,
          class_short
        ),
      by = "node"
    ) %>%
    left_join(taxonomy, by = "node") %>%
    mutate(
      network = network_label,
      is_obese_da_taxon = node %in% ob_da_set,
      display_group = node_group_label(node),
      label = format_sgb_names(node),
      clean_name = coalesce(union_clean_name, ob_clean_name, label),
      phylum = replace_na(phylum, "Unclassified"),
      family = replace_na(family, "Unclassified"),
      genus = replace_na(genus, "Unclassified")
    ) %>%
    arrange(desc(degree), desc(strength_crc_ob), node)
}

summarize_network_set <- function(edge_df, set_label, file_tag, plot_title_prefix) {
  nodes_all <- build_node_table(edge_df, paste0(set_label, ": Union DA-only obese-DA-participating rewiring"))

  nodes_by_class <- edge_df %>%
    group_split(edge_class) %>%
    lapply(function(df) build_node_table(df, as.character(unique(df$edge_class)))) %>%
    bind_rows()

  hubs_all <- bind_rows(nodes_all, nodes_by_class)

  class_summary <- edge_df %>%
    count(edge_class, da_side, name = "n_edges") %>%
    group_by(edge_class) %>%
    mutate(
      class_total = sum(n_edges),
      prop_within_class = n_edges / class_total
    ) %>%
    ungroup() %>%
    arrange(desc(class_total), edge_class, da_side)

  class_totals <- edge_df %>%
    count(edge_class, name = "n_edges") %>%
    mutate(
      prop_edges = n_edges / sum(n_edges),
      edge_class = factor(edge_class, levels = edge_class_order)
    ) %>%
    arrange(desc(n_edges))

  hub_summary <- hubs_all %>%
    group_by(network) %>%
    summarise(
      n_nodes = n(),
      n_obese_da_nodes = sum(is_obese_da_taxon),
      pct_obese_da_nodes = n_obese_da_nodes / n_nodes,
      n_ob_specific_nodes = sum(display_group == "Ob-specific obese-DA taxon"),
      .groups = "drop"
    ) %>%
    arrange(desc(n_nodes))

  top_summary <- hubs_all %>%
    group_by(network) %>%
    slice_head(n = top_n_hubs) %>%
    summarise(
      top_n = n(),
      n_obese_da_top_hubs = sum(is_obese_da_taxon),
      pct_obese_da_top_hubs = n_obese_da_top_hubs / top_n,
      n_ob_specific_top_hubs = sum(display_group == "Ob-specific obese-DA taxon"),
      .groups = "drop"
    ) %>%
    arrange(network)

  class_plot <- ggplot(class_totals, aes(x = fct_reorder(edge_class, n_edges), y = n_edges)) +
    geom_col(width = 0.72, fill = "#007A87", color = "black", linewidth = 0.2) +
    geom_text(aes(label = paste0(n_edges, "\n", scales::percent(prop_edges, accuracy = 1))),
              hjust = -0.08, size = 3.2, lineheight = 0.85) +
    coord_flip(clip = "off") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = paste0(plot_title_prefix, ": Rewiring Modes"),
      subtitle = paste0("Union DA-only; at least one endpoint is obese-DA; DGCA pValDiff < ", dgca_p_cutoff),
      x = NULL,
      y = "Edges"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.major.y = element_blank()
    )

  overall_plot_df <- nodes_all %>%
    slice_head(n = top_n_hubs) %>%
    mutate(label = fct_reorder(label, degree))

  hub_plot <- ggplot(overall_plot_df, aes(x = label, y = degree, fill = display_group)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.2) +
    coord_flip() +
    scale_fill_manual(values = c(
      "Ob-specific obese-DA taxon" = "#B2182B",
      "Other obese-DA taxon" = "#4D4D4D",
      "Other union-DA taxon" = "#BDBDBD"
    )) +
    labs(
      title = sprintf("%s: Top %d Hubs", plot_title_prefix, top_n_hubs),
      subtitle = paste0("CTRL vs CRC DGCA pValDiff < ", dgca_p_cutoff),
      x = NULL,
      y = "Degree",
      fill = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top",
      panel.grid.major.y = element_blank()
    )

  write_csv(edge_df, file.path(out_dir, paste0("Union_DAonly_ObeseDA_Participating_Rewiring_Edges_", file_tag, "_p01.csv")))
  write_csv(nodes_all, file.path(out_dir, paste0("Union_DAonly_ObeseDA_Participating_Rewiring_Nodes_", file_tag, "_p01.csv")))
  write_csv(hubs_all, file.path(out_dir, paste0("Union_DAonly_ObeseDA_Participating_Rewiring_Hubs_AllAndByClass_", file_tag, "_p01.csv")))
  write_csv(class_summary, file.path(out_dir, paste0("Union_DAonly_ObeseDA_Participating_Rewiring_Class_By_DA_Side_", file_tag, "_p01.csv")))
  write_csv(class_totals, file.path(out_dir, paste0("Union_DAonly_ObeseDA_Participating_Rewiring_Class_Totals_", file_tag, "_p01.csv")))
  write_csv(hub_summary, file.path(out_dir, paste0("Union_DAonly_ObeseDA_Participating_Rewiring_Hub_DA_Summary_", file_tag, "_p01.csv")))
  write_csv(top_summary, file.path(out_dir, paste0("Union_DAonly_ObeseDA_Participating_Rewiring_TopHub_DA_Summary_", file_tag, "_p01.csv")))

  for (class_i in levels(droplevels(edge_df$edge_class))) {
    class_label <- safe_file_label(class_i)
    edge_i <- edge_df %>% filter(edge_class == class_i)
    node_i <- nodes_by_class %>% filter(network == class_i)

    write_csv(edge_i, file.path(out_dir, paste0("Cytoscape_", file_tag, "_", class_label, "_Edges_p01.csv")))
    write_csv(node_i, file.path(out_dir, paste0("Cytoscape_", file_tag, "_", class_label, "_Nodes_p01.csv")))
  }

  ggsave(file.path(out_dir, paste0("Union_DAonly_ObeseDA_Participating_Rewiring_Class_Composition_", file_tag, "_p01.pdf")),
         class_plot, width = 7.2, height = 4.8)
  ggsave(file.path(out_dir, paste0("Union_DAonly_ObeseDA_Participating_Rewiring_Class_Composition_", file_tag, "_p01.png")),
         class_plot, width = 7.2, height = 4.8, dpi = 300)
  ggsave(file.path(out_dir, paste0("Union_DAonly_ObeseDA_Participating_Rewiring_Top25_Hubs_", file_tag, "_p01.pdf")),
         hub_plot, width = 8.4, height = 6.2)
  ggsave(file.path(out_dir, paste0("Union_DAonly_ObeseDA_Participating_Rewiring_Top25_Hubs_", file_tag, "_p01.png")),
         hub_plot, width = 8.4, height = 6.2, dpi = 300)

  list(
    edges = edge_df,
    nodes = nodes_all,
    hubs = hubs_all,
    class_totals = class_totals,
    top_hub_summary = top_summary
  )
}

same_normal_summary <- rewiring_edges %>%
  count(edge_class, same_mode_in_normal, edge_class_normal, name = "n_edges") %>%
  group_by(edge_class) %>%
  mutate(prop_within_obese_class = n_edges / sum(n_edges)) %>%
  ungroup() %>%
  arrange(edge_class, desc(same_mode_in_normal), desc(n_edges))

write_csv(
  same_normal_summary,
  file.path(out_dir, "Union_DAonly_ObeseDA_Participating_Rewiring_SameModeInNw_Summary_p01.csv")
)

same_trajectory_summary <- rewiring_edges %>%
  count(edge_class, trajectory__Normal, trajectory__Ob, name = "n_edges") %>%
  mutate(
    same_trajectory_as_nw = trajectory__Normal == trajectory__Ob
  ) %>%
  group_by(edge_class) %>%
  mutate(prop_within_obese_class = n_edges / sum(n_edges)) %>%
  ungroup() %>%
  arrange(edge_class, desc(same_trajectory_as_nw), desc(n_edges))

write_csv(
  same_trajectory_summary,
  file.path(out_dir, "Union_DAonly_ObeseDA_Participating_Rewiring_SameTrajectoryAsNw_Summary_p01.csv")
)

full_outputs <- summarize_network_set(
  rewiring_edges,
  set_label = "Full",
  file_tag = "Full",
  plot_title_prefix = "Union DA-only obese-DA participating rewiring"
)

contextualized_outputs <- summarize_network_set(
  contextualized_edges,
  set_label = "Excluding Nw-same",
  file_tag = "Excluding_NwSame",
  plot_title_prefix = "Union DA-only obese-contextualized rewiring"
)

trajectory_specific_outputs <- summarize_network_set(
  trajectory_specific_edges,
  set_label = "Excluding Nw-same trajectory",
  file_tag = "Excluding_NwSameTrajectory",
  plot_title_prefix = "Union DA-only obese-trajectory-specific rewiring"
)

cat("Saved union DA-only obese-DA participating rewiring analysis to:\n", out_dir, "\n")
cat("\nFull network size:\n")
cat("  Edges:", nrow(rewiring_edges), "\n")
cat("  Nodes:", n_distinct(c(rewiring_edges$Gene1, rewiring_edges$Gene2)), "\n")
cat("\nMajor edge-change modes (full):\n")
print(full_outputs$class_totals)
cat("\nEdges remaining after excluding same-mode Nw rewiring:\n")
print(contextualized_outputs$class_totals)
cat("\nEdges remaining after excluding identical Nw trajectories:\n")
print(trajectory_specific_outputs$class_totals)
