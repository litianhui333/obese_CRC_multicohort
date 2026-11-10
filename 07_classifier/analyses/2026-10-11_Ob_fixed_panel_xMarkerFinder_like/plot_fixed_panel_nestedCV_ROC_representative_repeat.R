#!/usr/bin/env Rscript

# Targeted plotting-only redraw of the Ob fixed-panel repeated nested-CV ROC.
# Uses a representative empirical ROC staircase per panel/classifier; labels
# retain the mean AUC across all five repeats. No model is refit.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

ROOT <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-10-11_Ob_fixed_panel_xMarkerFinder_like'
PREFIX <- 'restricted_ObAssoc_L1L2_nestedCV_'
pred_path <- file.path(ROOT, 'tables', paste0(PREFIX, 'oof_predictions.tsv'))
repeat_auc_path <- file.path(ROOT, 'tables', paste0(PREFIX, 'repeat_AUC.tsv'))
summary_path <- file.path(ROOT, 'tables', paste0(PREFIX, 'summary.tsv'))
stopifnot(file.exists(pred_path), file.exists(repeat_auc_path), file.exists(summary_path))
if (!requireNamespace('svglite', quietly = TRUE) || !requireNamespace('ragg', quietly = TRUE)) {
  stop('Required export packages are missing: svglite and ragg.')
}

panel_order <- c('species', 'sv', 'joint')
panel_names <- c(species = 'Species-only (12)', sv = 'SV-only (20)',
                 joint = 'Species + SV (22)')
panel_short <- c(species = 'Species', sv = 'SV', joint = 'Species + SV')
panel_colors <- c(species = '#477DB3', sv = '#4E9A91', joint = '#C44E52')
model_order <- c('Random forest', 'LASSO', 'Elastic Net', 'Radial SVM', 'XGBoost')

pred <- fread(pred_path)
repeat_auc <- fread(repeat_auc_path)
summary <- fread(summary_path)
stopifnot(all(c('Panel', 'Repeat', 'Model', 'sample_id', 'Disease', 'score') %in% names(pred)),
          all(c('Panel', 'Repeat', 'Model', 'AUC') %in% names(repeat_auc)),
          all(c('Panel', 'Model', 'Mean_AUC', 'SD_between_repeats', 'N_repeats') %in% names(summary)),
          setequal(unique(pred$Panel), panel_order),
          setequal(unique(pred$Model), model_order),
          setequal(unique(pred$Repeat), 1:5),
          all(pred$Disease %in% c('CTRL', 'CRC')),
          all(is.finite(pred$score)))
stopifnot(nrow(repeat_auc) == length(panel_order) * length(model_order) * 5L,
          nrow(summary) == length(panel_order) * length(model_order),
          all(summary$N_repeats == 5L))

# Pick the repeat whose pooled OOF AUC is nearest to the five-repeat mean.
# This is a transparent display convention; all labels still report the mean.
choice <- merge(
  repeat_auc,
  summary[, .(Panel, Model, Mean_AUC, SD_between_repeats, N_repeats)],
  by = c('Panel', 'Model'), all = FALSE
)
choice[, Distance_to_mean := abs(AUC - Mean_AUC)]
setorder(choice, Panel, Model, Distance_to_mean, Repeat)
choice <- choice[, .SD[1L], by = .(Panel, Model)]
stopifnot(nrow(choice) == 15L)

# Empirical ROC staircase with tied prediction scores entered as complete blocks.
roc_one <- function(d) {
  d <- copy(d)
  setorderv(d, 'score', order = -1L)
  n_pos <- sum(d$Disease == 'CRC')
  n_neg <- sum(d$Disease == 'CTRL')
  stopifnot(n_pos > 0L, n_neg > 0L)
  tmp <- data.table(
    score = d$score,
    FPR = cumsum(d$Disease == 'CTRL') / n_neg,
    TPR = cumsum(d$Disease == 'CRC') / n_pos
  )
  z <- tmp[!duplicated(score, fromLast = TRUE), .(FPR, TPR)]
  out <- rbind(data.table(FPR = 0, TPR = 0), z, data.table(FPR = 1, TPR = 1))
  out[, Point_order := seq_len(.N)]
  out
}

curve_list <- vector('list', nrow(choice))
for (i in seq_len(nrow(choice))) {
  ch <- choice[i]
  d <- pred[Panel == ch$Panel & Model == ch$Model & Repeat == ch$Repeat]
  z <- roc_one(d)
  z[, `:=`(
    Panel = ch$Panel,
    Model = ch$Model,
    Representative_repeat = ch$Repeat,
    Representative_AUC = ch$AUC,
    Mean_AUC = ch$Mean_AUC,
    SD_between_repeats = ch$SD_between_repeats,
    N_repeats = ch$N_repeats,
    N = nrow(d),
    N_CTRL = sum(d$Disease == 'CTRL'),
    N_CRC = sum(d$Disease == 'CRC')
  )]
  curve_list[[i]] <- z
}
curves <- rbindlist(curve_list)
curves[, Panel_label := factor(unname(panel_names[Panel]), levels = unname(panel_names[panel_order]))]
curves[, Model_label := factor(Model, levels = model_order)]
curves[, Curve_id := interaction(Panel, Model, drop = TRUE)]
setorder(curves, Model_label, Panel_label, Point_order)

ann <- unique(curves[, .(Panel, Model, Panel_label, Model_label, Mean_AUC,
                         SD_between_repeats, Representative_repeat, Representative_AUC,
                         N_repeats, N, N_CTRL, N_CRC)])
ann[, Label_y := c(species = 0.055, sv = 0.13, joint = 0.205)[Panel]]
ann[, Label := sprintf('%s  %.3f', unname(panel_short[Panel]), Mean_AUC)]

g <- ggplot(curves, aes(FPR, TPR, colour = Panel_label, group = Curve_id)) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed',
              colour = '#B8B8B8', linewidth = 0.35) +
  geom_step(direction = 'hv', linewidth = 0.58, alpha = 0.97) +
  geom_text(data = ann, aes(x = 0.97, y = Label_y, label = Label, colour = Panel_label),
            inherit.aes = FALSE, hjust = 1, size = 2.05, fontface = 'bold',
            show.legend = FALSE) +
  facet_wrap(~Model_label, nrow = 1) +
  scale_colour_manual(values = setNames(unname(panel_colors[panel_order]),
                                        unname(panel_names[panel_order])),
                      breaks = unname(panel_names[panel_order]), name = NULL) +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title = 'Obese CTRL vs CRC: fixed-panel nested-CV ROC',
    subtitle = 'Empirical staircase from the repeat nearest the five-repeat mean; labels show mean AUC',
    x = 'False-positive rate', y = 'True-positive rate',
    caption = 'Five repeats of nested 5-fold CV; panels are fixed at Species (12), SV (20), and Species + SV (22).'
  ) +
  theme_classic(base_size = 8.5, base_family = 'Arial') +
  theme(
    axis.line = element_line(linewidth = 0.34, colour = '#333333'),
    axis.ticks = element_line(linewidth = 0.28, colour = '#333333'),
    axis.text = element_text(size = 6.5, colour = '#333333'),
    axis.title = element_text(size = 7.8),
    strip.background = element_rect(fill = '#F1F1F1', colour = NA),
    strip.text = element_text(face = 'bold', size = 7.5),
    plot.title = element_text(face = 'bold', size = 10),
    plot.subtitle = element_text(size = 7.2, colour = '#555555'),
    plot.caption = element_text(size = 6.4, colour = '#555555', hjust = 0),
    legend.position = 'bottom',
    legend.text = element_text(size = 7),
    legend.key.width = grid::unit(0.45, 'cm'),
    panel.spacing = grid::unit(0.35, 'cm'),
    plot.margin = margin(5, 5, 5, 5)
  )

stem <- file.path(ROOT, 'figures', paste0(PREFIX, 'Fig_ROC_representative_repeat'))
pdf_path <- paste0(stem, '.pdf')
svg_path <- paste0(stem, '.svg')
png_path <- paste0(stem, '.png')
tiff_path <- paste0(stem, '.tiff')
data_path <- file.path(ROOT, 'data_for_figures',
                       paste0(PREFIX, 'Fig_ROC_representative_repeat.tsv'))
qc_path <- file.path(ROOT, 'tables',
                     paste0(PREFIX, 'Fig_ROC_representative_repeat_QC.tsv'))
manifest_path <- file.path(ROOT, 'tables',
                           paste0(PREFIX, 'Fig_ROC_representative_repeat_manifest.tsv'))

# Save the exact empirical curve coordinates, annotation text, and selected repeat.
plot_points <- curves[, .(
  Record_type = 'ROC point', Panel, Panel_label, Model, Model_label,
  Representative_repeat, Representative_AUC, Mean_AUC, SD_between_repeats,
  N_repeats, N, N_CTRL, N_CRC, FPR, TPR, Point_order, Label = NA_character_,
  Label_x = NA_real_, Label_y = NA_real_
)]
plot_labels <- ann[, .(
  Record_type = 'AUC annotation', Panel, Panel_label, Model, Model_label,
  Representative_repeat, Representative_AUC, Mean_AUC, SD_between_repeats,
  N_repeats, N, N_CTRL, N_CRC, FPR = NA_real_, TPR = NA_real_,
  Point_order = NA_integer_, Label, Label_x = 0.97, Label_y
)]
plot_data <- rbindlist(list(plot_points, plot_labels), use.names = TRUE, fill = TRUE)
fwrite(plot_data, data_path, sep = '\t', na = 'NA', quote = TRUE)
fwrite(choice[, .(Panel, Model, Repeat, AUC, Mean_AUC, SD_between_repeats,
                  Distance_to_mean, N_repeats)],
       manifest_path, sep = '\t', na = 'NA', quote = TRUE)

grDevices::cairo_pdf(pdf_path, width = 210 / 25.4, height = 90 / 25.4, family = 'Arial')
print(g)
dev.off()
svglite::svglite(svg_path, width = 210 / 25.4, height = 90 / 25.4)
print(g)
dev.off()
ragg::agg_png(png_path, width = 210 / 25.4, height = 90 / 25.4, units = 'in', res = 400)
print(g)
dev.off()
ragg::agg_tiff(tiff_path, width = 210 / 25.4, height = 90 / 25.4, units = 'in', res = 600,
               compression = 'lzw')
print(g)
dev.off()

endpoint_qc <- curves[, .(
  starts_at_origin = any(FPR == 0 & TPR == 0),
  ends_at_one = any(FPR == 1 & TPR == 1),
  monotone = all(diff(FPR) >= -1e-12) && all(diff(TPR) >= -1e-12)
), by = .(Panel, Model)]
label_match <- all(abs(ann$Mean_AUC -
  summary$Mean_AUC[match(paste(ann$Panel, ann$Model),
                         paste(summary$Panel, summary$Model))]) < 1e-12)
qc <- data.table(
  Check = c(
    'All three fixed Ob panels are plotted for all five algorithms',
    'Each representative repeat is the repeat nearest its five-repeat mean AUC',
    'Displayed AUC annotations exactly match the five-repeat summary means',
    'Every empirical ROC starts at (0,0), ends at (1,1), and is monotone',
    'Plotting data and representative-repeat manifest are saved',
    'PDF, SVG, PNG, and TIFF exports exist'
  ),
  Passed = c(
    uniqueN(curves[, .(Panel, Model)]) == 15L,
    all(choice$Repeat %in% 1:5 & is.finite(choice$Distance_to_mean)),
    label_match,
    nrow(endpoint_qc) == 15L && all(endpoint_qc$starts_at_origin) &&
      all(endpoint_qc$ends_at_one) && all(endpoint_qc$monotone),
    file.exists(data_path) && file.exists(manifest_path),
    all(file.exists(c(pdf_path, svg_path, png_path, tiff_path)))
  )
)
fwrite(qc, qc_path, sep = '\t')
stopifnot(all(qc$Passed))
print(qc)
print(choice[, .(Panel, Model, Repeat, AUC, Mean_AUC)])
cat('NESTEDCV_REPRESENTATIVE_ROC_COMPLETE\n')



