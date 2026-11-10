#!/usr/bin/env Rscript

# Targeted plotting-only redraw. Reads existing representative ROC coordinates
# and five-repeat AUC summaries; does not refit models or modify the original plot.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

ROOT <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-10-11_Ob_fixed_panel_xMarkerFinder_like'
PREFIX <- 'restricted_ObAssoc_L1L2_crossBMI_otherAlgorithms_'
curve_path <- file.path(ROOT, 'data_for_figures', paste0(PREFIX, 'Fig_transfer_ROC_representative_curves.tsv'))
auc_path <- file.path(ROOT, 'tables', paste0(PREFIX, 'AUC_summary.tsv'))

if (!file.exists(curve_path) || !file.exists(auc_path)) {
  stop('Required existing ROC coordinates or AUC summary are missing.')
}
if (!requireNamespace('svglite', quietly = TRUE) || !requireNamespace('ragg', quietly = TRUE)) {
  stop('Required export packages are missing: install.packages(c("svglite", "ragg")).')
}

curves <- fread(curve_path)
auc <- fread(auc_path)
panel_order <- c('species', 'sv', 'joint')
method_order <- c('Random forest', 'LASSO', 'Elastic Net', 'Radial SVM', 'XGBoost')
bmi_order <- c('Ob_test', 'Nw', 'Ow')
panel_names <- c(species = 'Species-only (12)', sv = 'SV-only (20)', joint = 'Species + SV (22)')
panel_colors <- c(species = '#477DB3', sv = '#4E9A91', joint = '#C44E52')
bmi_names <- c(Ob_test = 'Ob test')
method_names <- c('Random forest' = 'Random forest', 'LASSO' = 'LASSO',
                  'Elastic Net' = 'Elastic Net', 'Radial SVM' = 'Radial SVM',
                  'XGBoost' = 'XGBoost')

required_curve <- c('Panel', 'Model', 'BMI_evaluation', 'Representative_repeat',
                    'FPR', 'TPR', 'Mean_AUC', 'N')
required_auc <- c('Panel', 'Model', 'BMI_evaluation', 'Mean_AUC', 'N_repeats', 'N')
stopifnot(all(required_curve %in% names(curves)), all(required_auc %in% names(auc)))

curves <- curves[Panel %in% panel_order & Model %in% method_order & BMI_evaluation == 'Ob_test']
auc <- auc[Panel %in% panel_order & Model %in% method_order & BMI_evaluation == 'Ob_test']
setkey(auc, Panel, Model, BMI_evaluation)
setkey(curves, Panel, Model, BMI_evaluation)
stopifnot(nrow(auc) == length(panel_order) * length(method_order))
stopifnot(all(auc$N_repeats == 5L))
stopifnot(nrow(curves) > 0L, all(is.finite(curves$FPR)), all(is.finite(curves$TPR)))
stopifnot(all(curves$FPR >= 0 & curves$FPR <= 1), all(curves$TPR >= 0 & curves$TPR <= 1))

curves[, `:=`(
  Panel_label = factor(unname(panel_names[Panel]), levels = unname(panel_names[panel_order])),
  Method_label = factor(unname(method_names[Model]), levels = unname(method_names[method_order])),
  Curve_id = interaction(Panel, Model, BMI_evaluation, drop = TRUE)
)]

# Use each panel/method/BMI mean AUC as the annotation, while the drawn empirical
# ROC is the existing repeat selected nearest that five-repeat mean.
ann <- copy(auc)
ann[, `:=`(
  Panel_label = factor(unname(panel_names[Panel]), levels = unname(panel_names[panel_order])),
  Method_label = factor(unname(method_names[Model]), levels = unname(method_names[method_order])),
  Panel_short = c(species = 'Species-only', sv = 'SV-only', joint = 'Species + SV')[Panel],
  Label_y = c(species = 0.25, sv = 0.17, joint = 0.09)[Panel]
)]
ann[, Label := sprintf('%s  %.3f', Panel_short, Mean_AUC)]

g <- ggplot(curves, aes(FPR, TPR, colour = Panel_label, group = Curve_id)) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', colour = '#B8B8B8', linewidth = 0.35) +
  geom_step(direction = 'vh', linewidth = 0.55, alpha = 0.95) +
  geom_text(data = ann, aes(x = 0.97, y = Label_y, label = Label, colour = Panel_label),
            inherit.aes = FALSE, hjust = 1, size = 2.0, fontface = 'bold', show.legend = FALSE) +
  facet_wrap(~Method_label, ncol = 3) +
  scale_colour_manual(values = setNames(unname(panel_colors[panel_order]),
                                        unname(panel_names[panel_order])), drop = FALSE) +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1), expand = c(0, 0)) +
  coord_equal() +
  labs(x = 'False-positive rate', y = 'True-positive rate', colour = NULL,
       title = 'Obese CTRL vs CRC classification',
       subtitle = 'One panel per algorithm; curves use the repeat nearest mean AUC, labels show five-repeat mean AUC',
       caption = 'Ob test set: 47 CTRL and 50 CRC. Panels: Species-only (12), SV-only (20), and Species + SV (22).') +
  theme_classic(base_size = 8, base_family = 'Arial') +
  theme(
    axis.line = element_line(linewidth = 0.3, colour = '#333333'),
    axis.ticks = element_line(linewidth = 0.25, colour = '#333333'),
    axis.text = element_text(size = 6.2, colour = '#333333'),
    axis.title = element_text(size = 7.2),
    strip.background = element_rect(fill = '#F1F1F1', colour = NA),
    strip.text = element_text(face = 'bold', size = 7),
    plot.title = element_text(face = 'bold', size = 10),
    plot.subtitle = element_text(size = 7.2, colour = '#555555'),
    plot.caption = element_text(size = 6.5, colour = '#555555', hjust = 0),
    legend.position = 'bottom', legend.text = element_text(size = 7),
    legend.key.width = grid::unit(0.45, 'cm'),
    panel.spacing = grid::unit(0.35, 'cm'),
    plot.margin = margin(5, 5, 5, 5)
  )

stem <- file.path(ROOT, 'figures', paste0(PREFIX, 'Fig_Ob_only_three_panels_by_method'))
pdf_path <- paste0(stem, '.pdf')
svg_path <- paste0(stem, '.svg')
png_path <- paste0(stem, '.png')
tiff_path <- paste0(stem, '.tiff')
data_path <- file.path(ROOT, 'data_for_figures', paste0(PREFIX, 'Fig_Ob_only_three_panels_by_method.tsv'))
qc_path <- file.path(ROOT, 'tables', paste0(PREFIX, 'Fig_Ob_only_three_panels_by_method_QC.tsv'))

# Save both the exact ROC points and the exact AUC annotation positions/text.
setorder(curves, Method_label, Panel_label, FPR, TPR)
plot_points <- curves[, .(Record_type = 'ROC point', Panel, Panel_label, Model, Method_label,
                          BMI_evaluation, Representative_repeat, FPR, TPR,
                          Mean_AUC, N, Curve_id, Label = NA_character_)]
plot_labels <- ann[, .(Record_type = 'AUC annotation', Panel, Panel_label, Model, Method_label,
                       BMI_evaluation, Representative_repeat = NA_integer_,
                       FPR = 0.97, TPR = Label_y, Mean_AUC, N, Curve_id = NA_character_, Label)]
plot_data <- rbindlist(list(plot_points, plot_labels), use.names = TRUE, fill = TRUE)
fwrite(plot_data, data_path, sep = '\t', na = 'NA', quote = TRUE)

grDevices::cairo_pdf(pdf_path, width = 210 / 25.4, height = 155 / 25.4, family = 'Arial')
print(g)
dev.off()
svglite::svglite(svg_path, width = 210 / 25.4, height = 155 / 25.4)
print(g)
dev.off()
ragg::agg_png(png_path, width = 210 / 25.4, height = 155 / 25.4, units = 'in', res = 400)
print(g)
dev.off()
ragg::agg_tiff(tiff_path, width = 210 / 25.4, height = 155 / 25.4, units = 'in', res = 600,
               compression = 'lzw')
print(g)
dev.off()

endpoint_qc <- curves[, .(ok = any(FPR == 0 & TPR == 0) && any(FPR == 1 & TPR == 1)),
                      by = .(Panel, Model, BMI_evaluation)]
qc <- data.table(
  Check = c('All three panels and five algorithms are represented for Ob test only',
            'Each AUC annotation summarizes five repeats',
            'Each displayed curve is sourced from the saved representative empirical ROC coordinates',
            'Every curve starts at (0,0) and ends at (1,1)',
            'Exact plotted coordinates and AUC labels are saved',
            'PDF, SVG, PNG, and TIFF exports exist'),
  Passed = c(nrow(auc) == 15L && uniqueN(curves[, .(Panel, Model, BMI_evaluation)]) == 15L &&
               unique(curves$BMI_evaluation) == 'Ob_test',
             all(auc$N_repeats == 5L),
             all(curves$Representative_repeat %in% 1:5),
             all(endpoint_qc$ok), file.exists(data_path),
             all(file.exists(c(pdf_path, svg_path, png_path, tiff_path)))
  )
)
fwrite(qc, qc_path, sep = '\t')
stopifnot(all(qc$Passed))
print(qc)
cat('CROSS_BMI_THREE_PANELS_BY_METHOD_COMPLETE\n')
