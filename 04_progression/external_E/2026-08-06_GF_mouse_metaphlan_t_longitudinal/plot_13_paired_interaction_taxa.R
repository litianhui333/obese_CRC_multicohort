#!/usr/bin/env Rscript
# Figure contract: show taxa whose OB within-donor changes agree with the old
# CRC coefficient and are significant in the paired Diet x Time interaction.
# Each panel shows raw
# relative abundance and all paired trajectories; panel-specific y scales avoid
# suppressing low-abundance taxa.
suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })

root <- normalizePath(file.path(getwd(), "analyses", "2026-08-06_GF_mouse_metaphlan_t_longitudinal"), mustWork = TRUE)
results_dir <- file.path(root, "results"); figures_dir <- file.path(root, "figures")
comparison <- read_csv(file.path(results_dir, "Ob_Refined_paired_interaction_strict_OB_CRC_direction.csv"), show_col_types = FALSE)
abundance <- read_csv(file.path(results_dir, "t_strain_abundance_long_unfiltered.csv"), show_col_types = FALSE)

# Strict selection requires OB directional agreement with CRC in both transformations.
selected <- comparison %>% filter(selected_strict_OB_CRC_direction) %>%
  select(strain_id, clean_name) %>% arrange(clean_name) %>%
  mutate(panel_label = vapply(clean_name, function(x) paste(strwrap(x, width = 23), collapse = "\n"), character(1)))
stopifnot(nrow(selected) == 12L)
write_csv(selected, file.path(results_dir, "selected_12_paired_interaction_strict_OB_CRC_direction_taxa.csv"))

plot_data <- abundance %>% inner_join(selected, by = "strain_id") %>%
  mutate(Diet = factor(Diet, levels = c("LN", "OB"), labels = c("LN", "OB")),
         Time = factor(if_else(Week == 3, "WK3", "WK24"), levels = c("WK3", "WK24")),
         panel_label = factor(panel_label, levels = selected$panel_label))
write_csv(plot_data, file.path(results_dir, "source_data_12_paired_interaction_strict_OB_CRC_direction_taxa.csv"))

theme_set(theme_classic(base_size = 6.5, base_family = "Arial") + theme(
  axis.line = element_line(linewidth = .3), axis.ticks = element_line(linewidth = .3),
  strip.text = element_text(size = 5.8, face = "bold"), plot.title = element_text(size = 8.5, face = "bold"),
  legend.title = element_text(size = 6), legend.text = element_text(size = 5.8), panel.spacing = grid::unit(2, "mm")))

p <- ggplot(plot_data, aes(Time, RelativeAbundance_pct, group = MouseID)) +
  geom_line(aes(colour = Diet), linewidth = .42, alpha = .6) +
  geom_point(aes(colour = Diet, shape = Time), size = 1.55) +
  facet_wrap(~panel_label, ncol = 5, scales = "free_y") +
  scale_colour_manual(values = c(LN = "#3C78A8", OB = "#CC4125")) +
  scale_shape_manual(values = c(WK3 = 16, WK24 = 17)) +
  labs(title = "Paired longitudinal changes of 12 OB–CRC direction-consistent taxa", subtitle = "Paired Diet x Time interaction; raw relative abundance shown", x = NULL, y = "Relative abundance (%)", colour = "BMI", shape = "Time") +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95", colour = "grey35"),
        axis.text.x = element_text(size = 5.5), axis.text.y = element_text(size = 5.5), axis.title.y = element_text(size = 6.5),
        plot.margin = margin(4, 4, 4, 4))

stem <- file.path(figures_dir, "08_12_strict_OB_CRC_direction_taxa_WK3_WK24")
w <- 170 / 25.4; h <- 110 / 25.4
ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = 600, bg = "white")
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h); print(p); dev.off()
svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(p); dev.off()
ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in", res = 600); print(p); dev.off()

writeLines(c(
  "Core conclusion: these 12 taxa are significant under the paired Diet x Time interaction and their OB within-donor changes agree with historical CRC direction in AST and log10.",
  "Panel map: 12 taxon panels in 5 columns x 3 rows; x-axis is WK3/WK24, raw abundance is shown, and lines trace paired donors.",
  "Selection: raw interaction P < 0.01 OR BH-FDR < 0.10; OB mean delta and CRC coefficient must have the same sign in AST and log10.",
  "Source data: results/source_data_12_paired_interaction_strict_OB_CRC_direction_taxa.csv."
), file.path(root, "figure_contract_12_strict_OB_CRC_direction_taxa.txt"))
