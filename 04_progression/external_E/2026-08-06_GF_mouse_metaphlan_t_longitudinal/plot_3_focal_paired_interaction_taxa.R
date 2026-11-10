#!/usr/bin/env Rscript
# Figure contract: three focal paired-interaction taxa in separate horizontal
# facets, with raw abundance and paired donor trajectories shown at WK3/WK24.
suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })

root <- normalizePath(file.path(getwd(), "analyses", "2026-08-06_GF_mouse_metaphlan_t_longitudinal"), mustWork = TRUE)
results_dir <- file.path(root, "results"); figures_dir <- file.path(root, "figures")
focal <- tibble(strain_id = c("t__SGB15078", "t__SGB2303", "t__SGB4608"),
                panel_label = c("D. welbionis\nSGB15078", "A. onderdonkii\nSGB2303", "R. torques\nSGB4608"))
abundance <- read_csv(file.path(results_dir, "t_strain_abundance_long_unfiltered.csv"), show_col_types = FALSE)
plot_data <- abundance %>% inner_join(focal, by = "strain_id") %>%
  mutate(Diet = factor(Diet, levels = c("LN", "OB"), labels = c("LN", "OB")),
         Time = factor(if_else(Week == 3, "WK3", "WK24"), levels = c("WK3", "WK24")),
         panel_label = factor(panel_label, levels = focal$panel_label))
write_csv(plot_data, file.path(results_dir, "source_data_3_focal_paired_interaction_taxa.csv"))

theme_set(theme_classic(base_size = 6.5, base_family = "Arial") + theme(
  axis.line = element_line(linewidth = .3), axis.ticks = element_line(linewidth = .3),
  strip.text = element_text(size = 6.5, face = "bold"), plot.title = element_text(size = 8.5, face = "bold"),
  legend.title = element_text(size = 6), legend.text = element_text(size = 5.8), panel.spacing = grid::unit(3, "mm")))

p <- ggplot(plot_data, aes(Time, RelativeAbundance_pct, group = MouseID)) +
  geom_line(aes(colour = Diet), linewidth = .5, alpha = .65) +
  geom_point(aes(colour = Diet, shape = Time), size = 1.8) +
  facet_wrap(~panel_label, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = c(LN = "#3C78A8", OB = "#CC4125")) +
  scale_shape_manual(values = c(WK3 = 16, WK24 = 17)) +
  labs(title = "Paired longitudinal changes of three focal taxa", subtitle = "Raw relative abundance; colours denote BMI strata", x = NULL, y = "Relative abundance (%)", colour = "BMI", shape = "Time") +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95", colour = "grey35"),
        axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 6.5), plot.margin = margin(4, 4, 4, 4))

stem <- file.path(figures_dir, "10_3_focal_paired_interaction_taxa_WK3_WK24")
w <- 150 / 25.4; h <- 68 / 25.4
ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = 600, bg = "white")
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h); print(p); dev.off()
svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(p); dev.off()
ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in", res = 600); print(p); dev.off()
