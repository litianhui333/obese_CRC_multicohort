#!/usr/bin/env Rscript
# Figure contract: show the remaining 10 taxa after the three individually
# displayed taxa (SGB15078, SGB2303, SGB4608) are removed from the 13-taxon set.
suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })

root <- normalizePath(file.path(getwd(), "analyses", "2026-08-06_GF_mouse_metaphlan_t_longitudinal"), mustWork = TRUE)
results_dir <- file.path(root, "results"); figures_dir <- file.path(root, "figures")
selected <- read_csv(file.path(results_dir, "selected_12_paired_interaction_strict_OB_CRC_direction_taxa.csv"), show_col_types = FALSE) %>%
  filter(!strain_id %in% c("t__SGB15078", "t__SGB2303", "t__SGB4608")) %>% arrange(clean_name) %>%
  mutate(panel_label = vapply(clean_name, function(x) paste(strwrap(x, width = 23), collapse = "\n"), character(1)))
stopifnot(nrow(selected) == 9L)
abundance <- read_csv(file.path(results_dir, "t_strain_abundance_long_unfiltered.csv"), show_col_types = FALSE)
plot_data <- abundance %>% inner_join(selected, by = "strain_id") %>%
  mutate(Diet = factor(Diet, levels = c("LN", "OB"), labels = c("LN", "OB")),
         Time = factor(if_else(Week == 3, "WK3", "WK24"), levels = c("WK3", "WK24")),
         panel_label = factor(panel_label, levels = selected$panel_label))
write_csv(selected, file.path(results_dir, "selected_9_remaining_strict_OB_CRC_direction_taxa.csv"))
write_csv(plot_data, file.path(results_dir, "source_data_9_remaining_strict_OB_CRC_direction_taxa.csv"))

theme_set(theme_classic(base_size = 6.3, base_family = "Arial") + theme(
  axis.line = element_line(linewidth = .3), axis.ticks = element_line(linewidth = .3),
  strip.text = element_text(size = 5.8, face = "bold"), plot.title = element_text(size = 8.5, face = "bold"),
  legend.title = element_text(size = 6), legend.text = element_text(size = 5.8), panel.spacing = grid::unit(2, "mm")))

p <- ggplot(plot_data, aes(Time, RelativeAbundance_pct, group = MouseID)) +
  geom_line(aes(colour = Diet), linewidth = .42, alpha = .6) +
  geom_point(aes(colour = Diet, shape = Time), size = 1.55) +
  facet_wrap(~panel_label, ncol = 5, scales = "free_y") +
  scale_colour_manual(values = c(LN = "#3C78A8", OB = "#CC4125")) +
  scale_shape_manual(values = c(WK3 = 16, WK24 = 17)) +
  labs(title = "Paired longitudinal changes of the remaining 9 taxa", subtitle = "Three focal taxa are shown separately; raw relative abundance shown", x = NULL, y = "Relative abundance (%)", colour = "BMI", shape = "Time") +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95", colour = "grey35"),
        axis.text.x = element_text(size = 5.5), axis.text.y = element_text(size = 5.5), axis.title.y = element_text(size = 6.5), plot.margin = margin(4, 4, 4, 4))

stem <- file.path(figures_dir, "09_9_remaining_strict_OB_CRC_direction_taxa_WK3_WK24")
w <- 170 / 25.4; h <- 82 / 25.4
ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = 600, bg = "white")
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h); print(p); dev.off()
svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(p); dev.off()
ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in", res = 600); print(p); dev.off()
