.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(svglite)
  library(ragg)
})

PROJECT <- "H:/My Drive/proj_ob_AA_CRC"
OUT <- file.path(PROJECT, "analyses/2026-10-13_SV_abundance_sensitivity")
INPUT <- file.path(OUT, "tables/candidate_abundance_sensitivity_master.tsv")
FIG_STEM <- file.path(OUT, "figures/Fig2_effect_specific_direction_counts")
DATA_OUT <- file.path(OUT, "data_for_figures/Fig2_effect_specific_direction_counts.tsv")

stopifnot(file.exists(INPUT))
master <- fread(INPUT)
required <- c(
  "BMI_SV_ID", "BMI_enriched_group", "SV_type",
  "L2_nominal_abundance_adjusted", "Target_Beta_abundance_adjusted"
)
stopifnot(all(required %in% names(master)))

detail <- master[
  L2_nominal_abundance_adjusted %in% TRUE,
  .(
    BMI_SV_ID,
    BMI_enriched_group,
    SV_type,
    Target_Beta = Target_Beta_abundance_adjusted,
    Direction = fifelse(
      Target_Beta_abundance_adjusted > 0, "Increased",
      fifelse(Target_Beta_abundance_adjusted < 0, "Decreased", "Zero")
    )
  )
]
stopifnot(nrow(detail) == 48L, !anyDuplicated(detail$BMI_SV_ID), !any(detail$Direction == "Zero"))

bmi_levels <- c("Nw", "Ow", "Ob")
sv_levels <- c("dSGV", "vSGV")
direction_levels <- c("Increased", "Decreased")

counts <- detail[, .(Count = .N), by = .(BMI_enriched_group, SV_type, Direction)]
plot_dt <- merge(
  CJ(
    BMI_enriched_group = bmi_levels,
    SV_type = sv_levels,
    Direction = direction_levels,
    unique = TRUE
  ),
  counts,
  by = c("BMI_enriched_group", "SV_type", "Direction"),
  all.x = TRUE,
  sort = FALSE
)
plot_dt[is.na(Count), Count := 0L]

overall_totals <- detail[, .(Overall_group_total = .N), by = BMI_enriched_group]
panel_totals <- detail[, .(Panel_total = .N), by = .(BMI_enriched_group, SV_type)]
plot_dt <- merge(plot_dt, overall_totals, by = "BMI_enriched_group", all.x = TRUE, sort = FALSE)
plot_dt <- merge(plot_dt, panel_totals, by = c("BMI_enriched_group", "SV_type"), all.x = TRUE, sort = FALSE)

plot_dt[, BMI_order := match(BMI_enriched_group, bmi_levels)]
plot_dt[, SV_type_order := match(SV_type, sv_levels)]
plot_dt[, Direction_order := match(Direction, direction_levels)]
setorder(plot_dt, SV_type_order, BMI_order, Direction_order)
plot_dt[, Signed_count := fifelse(Direction == "Increased", Count, -Count)]
plot_dt[, BMI_label := factor(
  BMI_enriched_group,
  levels = bmi_levels,
  labels = sprintf("%s
(total n = %d)", bmi_levels, overall_totals[match(bmi_levels, BMI_enriched_group), Overall_group_total])
)]
type_totals <- detail[, .N, by = SV_type][match(sv_levels, SV_type), N]
plot_dt[, SV_panel := factor(
  SV_type,
  levels = sv_levels,
  labels = sprintf("%s (n = %d)", sv_levels, type_totals)
)]
plot_dt[, Direction := factor(Direction, levels = direction_levels)]

source_data <- plot_dt[, .(
  BMI_enriched_group,
  BMI_order,
  BMI_label = as.character(BMI_label),
  SV_type,
  SV_type_order,
  SV_panel = as.character(SV_panel),
  Direction = as.character(Direction),
  Direction_order,
  Count,
  Signed_count,
  Panel_total,
  Overall_group_total
)]
fwrite(source_data, DATA_OUT, sep = "	", na = "NA")

palette <- c("Increased" = "#C45A54", "Decreased" = "#4B79A7")

p <- ggplot(plot_dt, aes(x = BMI_label, y = Signed_count, fill = Direction)) +
  geom_hline(yintercept = 0, linewidth = 0.45, colour = "#333333") +
  geom_col(width = 0.62, colour = "white", linewidth = 0.25) +
  geom_text(
    data = plot_dt[Count > 0],
    aes(label = Count),
    vjust = ifelse(plot_dt[Count > 0]$Signed_count > 0, -0.35, 1.35),
    size = 2.7,
    family = "Arial",
    colour = "#222222",
    show.legend = FALSE
  ) +
  facet_grid(. ~ SV_panel) +
  scale_fill_manual(
    values = palette,
    breaks = c("Increased", "Decreased"),
    labels = c("Increased in CRC", "Decreased in CRC")
  ) +
  scale_y_continuous(
    breaks = seq(-10, 20, by = 5),
    labels = function(x) abs(x),
    expand = expansion(mult = c(0.08, 0.10))
  ) +
  labs(
    x = NULL,
    y = "Number of effect-specific SVs",
    fill = NULL
  ) +
  theme_classic(base_size = 9, base_family = "Arial") +
  theme(
    axis.line.x = element_blank(),
    axis.line.y = element_line(linewidth = 0.4, colour = "#222222"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(linewidth = 0.35, colour = "#222222"),
    axis.text.x = element_text(size = 8.3, colour = "#222222", lineheight = 0.92),
    axis.text.y = element_text(size = 8, colour = "#222222"),
    axis.title.y = element_text(size = 9, margin = margin(r = 6)),
    strip.background = element_blank(),
    strip.text = element_text(size = 9, face = "bold", colour = "#222222", margin = margin(b = 4)),
    panel.spacing.x = unit(12, "pt"),
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(size = 8),
    legend.key.width = unit(12, "pt"),
    legend.key.height = unit(7, "pt"),
    plot.margin = margin(5, 7, 5, 5)
  ) +
  guides(fill = guide_legend(override.aes = list(colour = NA)))

width_mm <- 145
height_mm <- 78
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

cairo_pdf(paste0(FIG_STEM, ".pdf"), width = width_in, height = height_in, family = "Arial")
print(p)
dev.off()

svglite(paste0(FIG_STEM, ".svg"), width = width_in, height = height_in, system_fonts = list(Arial = "Arial"))
print(p)
dev.off()

agg_png(paste0(FIG_STEM, ".png"), width = width_in, height = height_in, units = "in", res = 600, background = "white")
print(p)
dev.off()

agg_tiff(paste0(FIG_STEM, ".tiff"), width = width_in, height = height_in, units = "in", res = 600, compression = "lzw", background = "white")
print(p)
dev.off()

stopifnot(all(file.exists(paste0(FIG_STEM, c(".pdf", ".svg", ".png", ".tiff")))), file.exists(DATA_OUT))
cat("Saved figure and source data for", nrow(detail), "effect-specific SVs.
")
print(source_data)