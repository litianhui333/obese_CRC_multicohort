rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(igraph)
})

base_dir <- "G:/My Drive/proj_ob_AA_CRC"
analysis_dir <- file.path(base_dir, "analyses", "2026-06-05_overall_and_union_da_topology")
fig_dir <- file.path(analysis_dir, "figures")
tab_dir <- file.path(analysis_dir, "tables")
data_fig_dir <- file.path(analysis_dir, "data_for_figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_fig_dir, recursive = TRUE, showWarnings = FALSE)

dgca_root <- file.path(base_dir, "03.results", "06_network_analysis", "DGCA_Results")
seed_file <- file.path(
  dgca_root,
  "DA_only_subnetwork_Obese_Union_Heatmaps",
  "Union_All_BMI_DA_Highlight_Ob_SeedTaxa.csv"
)

comparisons <- c(
  "CTRL_Normal_vs_CRC_Normal",
  "CTRL_Ow_vs_CRC_Ow",
  "CTRL_Ob_vs_CRC_Ob"
)

bmi_labels <- c(
  "CTRL_Normal_vs_CRC_Normal" = "Normal",
  "CTRL_Ow_vs_CRC_Ow" = "Overweight",
  "CTRL_Ob_vs_CRC_Ob" = "Obese"
)

p_cutoff <- 0.01

read_comparison_pairs <- function(comp_name) {
  file_path <- file.path(dgca_root, comp_name, "table", paste0("DGCA_all_", comp_name, ".csv"))
  df <- read_csv(file_path, show_col_types = FALSE)

  ctrl_group <- sub("_vs_.*$", "", comp_name)
  crc_group <- sub("^.*_vs_", "", comp_name)

  df %>%
    transmute(
      Taxa1 = Gene1,
      Taxa2 = Gene2,
      ctrl_r = .data[[paste0(ctrl_group, "_cor")]],
      ctrl_p = .data[[paste0(ctrl_group, "_pVal")]],
      crc_r = .data[[paste0(crc_group, "_cor")]],
      crc_p = .data[[paste0(crc_group, "_pVal")]],
      BMI = bmi_labels[[comp_name]]
    )
}

all_pairs_raw <- bind_rows(lapply(comparisons, read_comparison_pairs))

overall_taxa <- sort(unique(c(all_pairs_raw$Taxa1, all_pairs_raw$Taxa2)))
union_da_taxa <- read_csv(seed_file, show_col_types = FALSE) %>%
  pull(feature) %>%
  unique() %>%
  sort()

taxa_universes <- list(
  Overall = overall_taxa,
  `Union DA` = union_da_taxa
)

calc_topology <- function(df_edges, all_vertices) {
  if (nrow(df_edges) == 0) {
    g <- make_empty_graph(n = length(all_vertices), directed = FALSE)
    V(g)$name <- all_vertices
  } else {
    g <- graph_from_data_frame(
      df_edges %>% select(Taxa1, Taxa2),
      directed = FALSE,
      vertices = tibble(name = all_vertices)
    ) %>%
      simplify(remove.multiple = TRUE, remove.loops = TRUE)
  }

  clust <- suppressWarnings(transitivity(g, type = "global"))
  if (is.na(clust)) clust <- 0

  mod_score <- if (ecount(g) > 0) {
    tryCatch(modularity(cluster_fast_greedy(g)), error = function(e) NA_real_)
  } else {
    NA_real_
  }

  tibble(
    Nodes = vcount(g),
    Edges = ecount(g),
    Pos_Edges = sum(df_edges$r > 0, na.rm = TRUE),
    Neg_Edges = sum(df_edges$r < 0, na.rm = TRUE),
    Avg_Degree = mean(degree(g)),
    Density = edge_density(g, loops = FALSE),
    Clustering_Coefficient = clust,
    Modularity = mod_score
  )
}

calc_universe_topology <- function(network_name, taxa_set) {
  universe_pairs <- all_pairs_raw %>%
    filter(Taxa1 %in% taxa_set, Taxa2 %in% taxa_set)

  bind_rows(lapply(split(universe_pairs, universe_pairs$BMI), function(df) {
    ctrl_edges <- df %>%
      filter(!is.na(ctrl_p), ctrl_p < p_cutoff) %>%
      transmute(Taxa1, Taxa2, r = ctrl_r)

    crc_edges <- df %>%
      filter(!is.na(crc_p), crc_p < p_cutoff) %>%
      transmute(Taxa1, Taxa2, r = crc_r)

    bind_rows(
      bind_cols(
        tibble(Network = network_name, Rule = "raw p < 0.01", BMI = unique(df$BMI), Disease = "CTRL"),
        calc_topology(ctrl_edges, taxa_set)
      ),
      bind_cols(
        tibble(Network = network_name, Rule = "raw p < 0.01", BMI = unique(df$BMI), Disease = "CRC"),
        calc_topology(crc_edges, taxa_set)
      )
    )
  }))
}

topology_df <- bind_rows(lapply(names(taxa_universes), function(network_name) {
  calc_universe_topology(network_name, taxa_universes[[network_name]])
})) %>%
  mutate(
    Network = factor(Network, levels = c("Overall", "Union DA")),
    BMI = factor(BMI, levels = c("Normal", "Overweight", "Obese")),
    Disease = factor(Disease, levels = c("CTRL", "CRC"))
  ) %>%
  arrange(Network, BMI, Disease)

delta_df <- topology_df %>%
  pivot_wider(
    names_from = Disease,
    values_from = c(
      Nodes, Edges, Pos_Edges, Neg_Edges, Avg_Degree,
      Density, Clustering_Coefficient, Modularity
    )
  ) %>%
  mutate(
    delta_Edges = Edges_CRC - Edges_CTRL,
    delta_Pos_Edges = Pos_Edges_CRC - Pos_Edges_CTRL,
    delta_Neg_Edges = Neg_Edges_CRC - Neg_Edges_CTRL,
    delta_Avg_Degree = Avg_Degree_CRC - Avg_Degree_CTRL,
    delta_Density = Density_CRC - Density_CTRL,
    delta_Clustering_Coefficient = Clustering_Coefficient_CRC - Clustering_Coefficient_CTRL,
    delta_Modularity = Modularity_CRC - Modularity_CTRL
  ) %>%
  arrange(Network, BMI)

plot_df <- topology_df %>%
  select(
    Network, BMI, Disease, Edges, Avg_Degree, Density,
    Clustering_Coefficient, Modularity
  ) %>%
  pivot_longer(
    cols = c(Edges, Avg_Degree, Density, Clustering_Coefficient, Modularity),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("Edges", "Density", "Avg_Degree", "Clustering_Coefficient", "Modularity"),
      labels = c("Edge count", "Density", "Mean degree", "Clustering", "Modularity")
    ),
    Network = factor(Network, levels = c("Overall", "Union DA")),
    BMI = factor(BMI, levels = c("Normal", "Overweight", "Obese")),
    Disease = factor(Disease, levels = c("CTRL", "CRC"))
  )

write_csv(topology_df, file.path(tab_dir, "overall_and_union_DA_CTRL_vs_CRC_topology_p01.csv"))
write_csv(delta_df, file.path(tab_dir, "overall_and_union_DA_CTRL_vs_CRC_topology_delta_p01.csv"))
write_csv(
  topology_df %>% filter(Network == "Overall") %>% select(-Network),
  file.path(tab_dir, "overall_CTRL_vs_CRC_topology_p01.csv")
)
write_csv(
  topology_df %>% filter(Network == "Union DA") %>% select(-Network),
  file.path(tab_dir, "union_DA_CTRL_vs_CRC_topology_p01.csv")
)
write_csv(
  delta_df %>% filter(Network == "Overall") %>% select(-Network),
  file.path(tab_dir, "overall_CTRL_vs_CRC_topology_delta_p01.csv")
)
write_csv(
  delta_df %>% filter(Network == "Union DA") %>% select(-Network),
  file.path(tab_dir, "union_DA_CTRL_vs_CRC_topology_delta_p01.csv")
)
write_tsv(plot_df, file.path(data_fig_dir, "overall_and_union_DA_CTRL_vs_CRC_topology_p01_plotting_data.tsv"))

bmi_colors <- c(
  "Normal" = "#1f6f5f",
  "Overweight" = "#c46a2d",
  "Obese" = "#7a4c8a"
)

base_theme <- theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "#f3f3f3", color = "black", linewidth = 0.3),
    strip.text = element_text(face = "bold", size = 10),
    axis.text = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    panel.spacing = unit(0.9, "lines")
  )

make_panel <- function(df, metric_name, panel_title, y_lab) {
  ggplot(
    df %>% filter(Metric == metric_name),
    aes(x = Disease, y = Value, group = BMI, color = BMI)
  ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.4) +
    scale_color_manual(values = bmi_colors) +
    scale_x_discrete(expand = expansion(mult = c(0.08, 0.08))) +
    labs(title = panel_title, x = NULL, y = y_lab) +
    base_theme
}

build_topology_plot <- function(df, title_prefix = "") {
  p_edges <- make_panel(df, "Edge count", paste0(title_prefix, "A. Edge Count"), "Edges")
  p_density <- make_panel(df, "Density", paste0(title_prefix, "B. Density"), "Density")
  p_degree <- make_panel(df, "Mean degree", paste0(title_prefix, "C. Mean Degree"), "Mean degree")
p_cluster <- make_panel(df, "Clustering", paste0(title_prefix, "D. Clustering Coefficient"), "Clustering coefficient")
  p_mod <- make_panel(df, "Modularity", paste0(title_prefix, "E. Modularity"), "Modularity")

  (
    (p_edges | p_density | p_degree) /
      (p_cluster | p_mod | plot_spacer())
  ) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

overall_plot <- build_topology_plot(plot_df %>% filter(Network == "Overall"))
union_da_plot <- build_topology_plot(plot_df %>% filter(Network == "Union DA"))

combined_plot <- ggplot(
  plot_df,
  aes(x = Disease, y = Value, group = BMI, color = BMI)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.1) +
  facet_wrap(vars(Network, Metric), nrow = 2, scales = "free_y") +
  scale_color_manual(values = bmi_colors) +
  scale_x_discrete(expand = expansion(mult = c(0.08, 0.08))) +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 8)
  )

save_plot <- function(plot, stem, width, height) {
  pdf_file <- file.path(fig_dir, paste0(stem, ".pdf"))
  png_file <- file.path(fig_dir, paste0(stem, ".png"))
  tryCatch(
    ggsave(pdf_file, plot, width = width, height = height, device = cairo_pdf),
    error = function(e) {
      message("cairo_pdf failed for ", basename(pdf_file), "; retrying with grDevices::pdf().")
      ggsave(pdf_file, plot, width = width, height = height, device = grDevices::pdf)
    }
  )
  ggsave(png_file, plot, width = width, height = height, dpi = 300)
  c(pdf = pdf_file, png = png_file)
}

save_plot(overall_plot, "overall_CTRL_vs_CRC_topology_p01", 7.2, 4.8)
save_plot(union_da_plot, "union_DA_CTRL_vs_CRC_topology_p01", 7.2, 4.8)
save_plot(combined_plot, "overall_and_union_DA_CTRL_vs_CRC_topology_p01", 9.2, 4.8)

nonredundant_plot_df <- topology_df %>%
  select(
    Network, BMI, Disease, Density, Clustering_Coefficient,
    Modularity, Pos_Edges, Neg_Edges
  ) %>%
  pivot_longer(
    cols = c(Density, Clustering_Coefficient, Modularity, Pos_Edges, Neg_Edges),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("Density", "Clustering_Coefficient", "Modularity", "Pos_Edges", "Neg_Edges"),
      labels = c("Density", "Clustering", "Modularity", "Positive edges", "Negative edges")
    ),
    Network = factor(Network, levels = c("Overall", "Union DA")),
    BMI = factor(BMI, levels = c("Normal", "Overweight", "Obese")),
    Disease = factor(Disease, levels = c("CTRL", "CRC"))
  )

write_tsv(
  nonredundant_plot_df,
  file.path(data_fig_dir, "overall_and_union_DA_CTRL_vs_CRC_topology_p01_nonredundant_plotting_data.tsv")
)

make_nonredundant_panel <- function(df, metric_name, panel_title, y_lab) {
  ggplot(
    df %>% filter(Metric == metric_name),
    aes(x = Disease, y = Value, group = BMI, color = BMI)
  ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.4) +
    scale_color_manual(values = bmi_colors) +
    scale_x_discrete(expand = expansion(mult = c(0.08, 0.08))) +
    labs(title = panel_title, x = NULL, y = y_lab) +
    base_theme
}

build_nonredundant_plot <- function(df) {
  p_density <- make_nonredundant_panel(df, "Density", "A. Density", "Density")
  p_cluster <- make_nonredundant_panel(df, "Clustering", "B. Clustering", "Clustering coefficient")
  p_mod <- make_nonredundant_panel(df, "Modularity", "C. Modularity", "Modularity")
  p_pos <- make_nonredundant_panel(df, "Positive edges", "D. Positive Edges", "Positive edges")
  p_neg <- make_nonredundant_panel(df, "Negative edges", "E. Negative Edges", "Negative edges")

  (
    (p_density | p_cluster | p_mod) /
      (p_pos | p_neg | plot_spacer())
  ) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

overall_nonredundant_plot <- build_nonredundant_plot(
  nonredundant_plot_df %>% filter(Network == "Overall")
)
union_da_nonredundant_plot <- build_nonredundant_plot(
  nonredundant_plot_df %>% filter(Network == "Union DA")
)

save_plot(overall_nonredundant_plot, "overall_CTRL_vs_CRC_topology_p01_nonredundant", 7.2, 4.8)
save_plot(union_da_nonredundant_plot, "union_DA_CTRL_vs_CRC_topology_p01_nonredundant", 7.2, 4.8)

edge_reduction_df <- delta_df %>%
  transmute(
    Network,
    BMI,
    `Positive edges` = (Pos_Edges_CTRL - Pos_Edges_CRC) / Pos_Edges_CTRL,
    `Negative edges` = (Neg_Edges_CTRL - Neg_Edges_CRC) / Neg_Edges_CTRL
  ) %>%
  pivot_longer(
    cols = c(`Positive edges`, `Negative edges`),
    names_to = "Edge_Type",
    values_to = "Reduction_Fraction"
  ) %>%
  mutate(
    Network = factor(Network, levels = c("Overall", "Union DA")),
    BMI = factor(BMI, levels = c("Normal", "Overweight", "Obese")),
    Edge_Type = factor(Edge_Type, levels = c("Positive edges", "Negative edges"))
  )

compact4_plot_df <- topology_df %>%
  select(Network, BMI, Disease, Density, Clustering_Coefficient, Modularity) %>%
  pivot_longer(
    cols = c(Density, Clustering_Coefficient, Modularity),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("Density", "Clustering_Coefficient", "Modularity"),
      labels = c("Density", "Clustering", "Modularity")
    ),
    Network = factor(Network, levels = c("Overall", "Union DA")),
    BMI = factor(BMI, levels = c("Normal", "Overweight", "Obese")),
    Disease = factor(Disease, levels = c("CTRL", "CRC"))
  )

write_tsv(
  compact4_plot_df,
  file.path(data_fig_dir, "overall_and_union_DA_CTRL_vs_CRC_topology_p01_compact4_line_plotting_data.tsv")
)
write_tsv(
  edge_reduction_df,
  file.path(data_fig_dir, "overall_and_union_DA_CTRL_vs_CRC_edge_reduction_fraction_p01_plotting_data.tsv")
)

build_compact4_plot <- function(network_name) {
  line_df <- compact4_plot_df %>% filter(Network == network_name)
  bar_df <- edge_reduction_df %>% filter(Network == network_name)

  p_density <- make_nonredundant_panel(line_df, "Density", "A. Density", "Density")
  p_cluster <- make_nonredundant_panel(line_df, "Clustering", "B. Clustering", "Clustering coefficient")
  p_mod <- make_nonredundant_panel(line_df, "Modularity", "C. Modularity", "Modularity")

  p_reduction <- ggplot(
    bar_df,
    aes(x = BMI, y = Reduction_Fraction, fill = Edge_Type)
  ) +
    geom_col(position = position_dodge(width = 0.72), width = 0.62, color = "black", linewidth = 0.2) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
    scale_fill_manual(values = c("Positive edges" = "#2b8c67", "Negative edges" = "#b75a5a")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = "D. Edge Reduction",
      x = NULL,
      y = "Reduction from CTRL to CRC",
      fill = NULL
    ) +
    base_theme +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "bottom"
    )

  (
    (p_density | p_cluster) /
      (p_mod | p_reduction)
  ) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

overall_compact4_plot <- build_compact4_plot("Overall")
union_da_compact4_plot <- build_compact4_plot("Union DA")

save_plot(overall_compact4_plot, "overall_CTRL_vs_CRC_topology_p01_compact4", 6.4, 5.2)
save_plot(union_da_compact4_plot, "union_DA_CTRL_vs_CRC_topology_p01_compact4", 6.4, 5.2)

message("Overall and union DA topology outputs saved to: ", analysis_dir)
