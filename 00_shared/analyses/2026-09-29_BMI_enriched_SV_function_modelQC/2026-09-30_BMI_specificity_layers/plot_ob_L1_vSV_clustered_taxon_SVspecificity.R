.libPaths(c('D:/R/library', .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
})

# Targeted update: replace the stage-mean panel with audited taxon-abundance status.
# The clustered sample matrix and all unrelated figures remain unchanged.
ROOT <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-09-29_BMI_enriched_SV_function_modelQC/2026-09-30_BMI_specificity_layers'
SOURCE_STEM <- 'Heatmap_L1_Ob_vSGV_multicolor_clustered_compact'
STYLE <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else 'raw_z'
stopifnot(STYLE %in% c('raw_z', 'ctrl_centered'))
STEM <- if (STYLE == 'ctrl_centered') {
  'Heatmap_L1_Ob_vSGV_ctrl_centered_taxon_SVspecificity'
} else {
  'Heatmap_L1_Ob_vSGV_multicolor_clustered_taxon_SVspecificity'
}
matrix_suffix <- ''
if (STYLE == 'ctrl_centered') matrix_suffix <- '_unclipped'

read_plot_data <- function(suffix = '') {
  fread(file.path(ROOT, 'data_for_figures', paste0(SOURCE_STEM, suffix, '.tsv')))
}
write_plot_data <- function(x, suffix = '') {
  fwrite(x, file.path(ROOT, 'data_for_figures', paste0(STEM, suffix, '.tsv')),
         sep = '\t', na = 'NA')
}

mat_dt <- read_plot_data(matrix_suffix)
id_col <- names(mat_dt)[1]
row_ids <- mat_dt[[id_col]]
plotmat <- as.matrix(mat_dt[, -1])
storage.mode(plotmat) <- 'double'
rownames(plotmat) <- row_ids

rows <- read_plot_data('_rows')
setorder(rows, Row_order)
stopifnot(!anyNA(rows$L2_nominal), sum(rows$L2_nominal) == 11L)
rows[, SV_specificity := factor(
  ifelse(L2_nominal, 'Ob-specific SV', 'Other Ob L1 SV'),
  levels = c('Ob-specific SV', 'Other Ob L1 SV'))]
columns <- read_plot_data('_columns')
palette <- read_plot_data('_palette')
# Restore the earlier bright seven-anchor scale requested for this panel.
vlim <- max(abs(palette$Value), na.rm = TRUE)
palette <- data.table(
  Value = vlim * c(-1, -2/3, -1/3, 0, 1/3, 2/3, 1),
  Color = c('#313695', '#4575B4', '#66C2A5', '#FFFFF0', '#FEE08B',
            '#F46D43', '#A50026')
)
stopifnot(
  identical(rownames(plotmat), rows$BMI_SV_ID),
  identical(colnames(plotmat), columns$SGV_ID),
  nrow(plotmat) == 32L,
  ncol(plotmat) == 688L
)

# Remove samples that have no finite value across all 32 displayed vSVs.
# Do not replace these missing profiles with zero, because zero is a valid row-z value.
all_na_columns <- colSums(is.finite(plotmat)) == 0
removed_all_na <- copy(columns[all_na_columns])
removed_all_na[, `:=`(
  Removal_reason = 'All 32 displayed Ob L1 vSV values are NA',
  Original_column_order = Column_order
)]
plotmat <- plotmat[, !all_na_columns, drop = FALSE]
columns <- copy(columns[!all_na_columns])
columns[, Column_order := seq_len(.N)]
stopifnot(nrow(removed_all_na) == 15L, ncol(plotmat) == 673L,
          identical(colnames(plotmat), columns$SGV_ID))

# One abundance classification per SGB, using the current audited robust membership.
crosswalk <- fread(file.path(ROOT, 'tables', 'Ob_SV_taxon_abundance_crosswalk.tsv'))
taxon_status <- crosswalk[, .(
  Current_robust_in_Ob = any(Current_robust_in_Ob %in% TRUE, na.rm = TRUE),
  Current_robust_class = {
    z <- Current_robust_class[!is.na(Current_robust_class) & nzchar(Current_robust_class)]
    if ('Ob-specific' %in% z) 'Ob-specific' else if (length(z)) z[1] else NA_character_
  }
), by = SGB]
rows <- merge(rows, taxon_status, by = 'SGB', all.x = TRUE, sort = FALSE)
setorder(rows, Row_order)
rows[, Taxon_abundance_status := fcase(
  Current_robust_in_Ob %in% TRUE & Current_robust_class == 'Ob-specific', 'Ob-specific',
  Current_robust_in_Ob %in% TRUE, 'Ob-associated (non-specific)',
  default = 'Not significant in Ob'
)]
status_levels <- c('Ob-specific', 'Ob-associated (non-specific)', 'Not significant in Ob')
rows[, Taxon_abundance_status := factor(Taxon_abundance_status, levels = status_levels)]
stopifnot(!anyNA(rows$Taxon_abundance_status), identical(rows$BMI_SV_ID, rownames(plotmat)))

# Manuscript label: taxon SGB : genomic region; omit internal ObSV IDs.
rows[, Region_label := sub('^.*:', '', SV_feature)]
rows[, Region_label := gsub('_', '\u2013', Region_label, fixed = TRUE)]
rows[, Region_label := gsub(';', '; ', Region_label, fixed = TRUE)]
rows[, Row_label := paste0(SGB_label, ' ', sub('^t__', '', SGB), ' : ', Region_label)]
rows[, N_segments := lengths(strsplit(Region_label, ';', fixed = TRUE))]
rows[, Label_group := ifelse(Target_Beta > 0, 'Ob+', 'Ob-')]
rows[, Label_index := seq_len(.N), by = Label_group]
label_candidates <- copy(rows[N_segments == 1L])
label_candidates[, `:=`(L2_rank = as.integer(L2_nominal),
                        Absolute_effect = abs(Target_Beta))]
setorder(label_candidates, -L2_rank, -Absolute_effect, BMI_SV_ID)
label_ids <- character()
for (id in label_candidates$BMI_SV_ID) {
  if (length(label_ids) >= 6L) break
  candidate <- rows[BMI_SV_ID == id]
  existing <- rows[BMI_SV_ID %chin% label_ids]
  close <- existing[Label_group == candidate$Label_group &
                      abs(Label_index - candidate$Label_index) < 6L]
  if (nrow(close) == 0L) label_ids <- c(label_ids, id)
}
rows[, Selected_row_label := BMI_SV_ID %chin% label_ids]
rows[, Display_row_label := fifelse(Selected_row_label, Row_label, '')]
rows[, Label_selection_reason := fifelse(Selected_row_label,
  'Single-region label; prioritised by L2_nominal and absolute effect with spacing', '')]
stopifnot(length(label_ids) >= 5L, length(label_ids) <= 6L)

# Optional effect-centred display: preserve z-score units, but set each row
# observed Ob-CTRL mean to zero before clipping. This makes warming/cooling show
# the direction of change from the CTRL reference rather than the absolute sign
# of the original all-sample z-score.
plotmat_unclipped <- NULL
legend_title <- 'vSV row z-score (clipped)'
row_titles <- c('Ob +', 'Ob -')
figure_title <- paste0('Ob L1 vSV | ', nrow(rows), ' SVs')
n_raw_endpoint_concordant <- NA_integer_
if (STYLE == 'ctrl_centered') {
  ctrl_cols <- columns$Disease_plot == 'CTRL'
  crc_cols <- columns$Disease_plot == 'CRC'
  ctrl_mean <- rowMeans(plotmat[, ctrl_cols, drop = FALSE], na.rm = TRUE)
  crc_mean <- rowMeans(plotmat[, crc_cols, drop = FALSE], na.rm = TRUE)
  stopifnot(all(is.finite(ctrl_mean)), all(is.finite(crc_mean)))
  rows[, `:=`(
    Mean_CTRL_before_centering = ctrl_mean[BMI_SV_ID],
    Mean_CRC_before_centering = crc_mean[BMI_SV_ID]
  )]
  rows[, Raw_CTRL_to_CRC_delta := Mean_CRC_before_centering - Mean_CTRL_before_centering]
  rows[, Raw_endpoint_concordant := sign(Raw_CTRL_to_CRC_delta) == sign(Target_Beta)]
  n_raw_endpoint_concordant <- sum(rows$Raw_endpoint_concordant)
  plotmat_unclipped <- sweep(plotmat, 1, ctrl_mean, FUN = '-')
  plotmat <- pmax(pmin(plotmat_unclipped, vlim), -vlim)
  legend_title <- 'vSV z-score relative to Ob CTRL'
  row_titles <- c('CRC higher', 'CRC lower')
  figure_title <- paste0('Ob L1 vSV | CTRL-centered | ', nrow(rows), ' SVs')
}

heatmap_col <- colorRamp2(palette$Value, palette$Color)
legend_at <- seq(-vlim, vlim, length.out = 7)
status_col <- c(
  'Ob-specific' = '#7B3294',
  'Ob-associated (non-specific)' = '#008C95',
  'Not significant in Ob' = '#D9D9D9'
)
sv_specificity_col <- c('Ob-specific SV' = '#D55E00', 'Other Ob L1 SV' = '#D9D9D9')
stage_levels <- c('CTRL', 'AA', 'CRC')
stage_names <- c('CTRL', 'Adenoma', 'CRC')
column_split <- factor(columns$Disease_plot, levels = stage_levels, labels = stage_names)
row_split <- factor(ifelse(rows$Target_Beta > 0, '+', '-'),
                    levels = unique(ifelse(rows$Target_Beta > 0, '+', '-')))

top_anno <- HeatmapAnnotation(
  Stage = column_split,
  col = list(Stage = c(CTRL = '#B8B8B8', Adenoma = '#F4A261', CRC = '#B94A48')),
  annotation_name_gp = gpar(fontsize = 10),
  simple_anno_size = unit(2.5, 'mm'),
  show_legend = FALSE
)

heatmap_legend <- list(
  title = legend_title,
  at = legend_at,
  labels = round(legend_at, 2),
  direction = 'horizontal',
  legend_width = unit(62, 'mm'),
  title_gp = gpar(fontsize = 10, fontface = 'bold'),
  labels_gp = gpar(fontsize = 9)
)

hm <- Heatmap(
  plotmat,
  name = 'vSV z',
  col = heatmap_col,
  na_col = '#F0F0F0',
  column_split = column_split,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  row_split = row_split,
  row_title = row_titles,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 10),
  column_title_gp = gpar(fontsize = 10),
  column_gap = unit(2, 'mm'),
  row_gap = unit(2.5, 'mm'),
  top_annotation = top_anno,
  use_raster = FALSE,
  border = TRUE,
  rect_gp = gpar(col = NA),
  heatmap_legend_param = heatmap_legend
)

status_matrix <- matrix(
  as.character(rows$Taxon_abundance_status),
  ncol = 1,
  dimnames = list(rows$BMI_SV_ID, 'Taxon status')
)
hs <- Heatmap(
  status_matrix,
  name = 'Taxon abundance status',
  col = status_col,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cluster_row_slices = FALSE,
  row_split = row_split,
  show_row_names = FALSE,
  show_column_names = FALSE,
  width = unit(4, 'mm'),
  border = TRUE,
  heatmap_legend_param = list(
    title = 'Taxon abundance status',
    at = status_levels,
    labels = status_levels,
    nrow = 1,
    direction = 'horizontal',
    title_gp = gpar(fontsize = 10, fontface = 'bold'),
    labels_gp = gpar(fontsize = 9)
  )
)
specificity_matrix <- matrix(as.character(rows$SV_specificity), ncol = 1,
                             dimnames = list(rows$BMI_SV_ID, 'SV specificity'))
spec_hs <- Heatmap(
  specificity_matrix, name = 'SV specificity', col = sv_specificity_col,
  cluster_rows = FALSE, cluster_columns = FALSE, cluster_row_slices = FALSE,
  row_split = row_split, show_row_names = TRUE,
  row_labels = rows$Display_row_label, row_names_side = 'right',
  row_names_gp = gpar(fontsize = 6.2), row_names_max_width = unit(68, 'mm'),
  show_column_names = FALSE,
  width = unit(4, 'mm'), border = TRUE,
  heatmap_legend_param = list(title = 'SV specificity', at = names(sv_specificity_col),
    labels = names(sv_specificity_col), nrow = 1, direction = 'horizontal')
)

width <- 9.0
height <- 4.8
render <- function() {
  draw(
    hm + hs + spec_hs,
    heatmap_legend_side = 'bottom',
    merge_legends = TRUE,
    column_title = figure_title,
    column_title_gp = gpar(fontsize = 13, fontfamily = 'Arial'),
    padding = unit(c(5, 5, 5, 5), 'mm')
  )
}

path <- file.path(ROOT, 'figures', STEM)
grDevices::cairo_pdf(paste0(path, '.pdf'), width = width, height = height, family = 'Arial')
render()
dev.off()
png(paste0(path, '.png'), width = width, height = height, units = 'in', res = 300, type = 'cairo')
render()
dev.off()
svglite::svglite(paste0(path, '.svg'), width = width, height = height)
render()
dev.off()

# Save the exact data used by this revised figure. There are no stage-mean columns.
write_plot_data(data.table(BMI_SV_ID = rownames(plotmat), as.data.table(plotmat)))
if (STYLE == 'ctrl_centered') {
  write_plot_data(data.table(BMI_SV_ID = rownames(plotmat_unclipped),
                             as.data.table(plotmat_unclipped)), '_unclipped')
  write_plot_data(rows[, .(BMI_SV_ID, Target_Beta, Mean_CTRL_before_centering,
                           Mean_CRC_before_centering, Raw_CTRL_to_CRC_delta,
                           Raw_endpoint_concordant)], '_centering_audit')
}
write_plot_data(rows, '_rows')
write_plot_data(rows[, .(BMI_SV_ID, SV_specificity, L2_nominal)], '_sv_specificity')
write_plot_data(columns, '_columns')
write_plot_data(removed_all_na, '_removed_all_NA_samples')
write_plot_data(palette, '_heatmap_palette')
write_plot_data(data.table(Status = names(status_col), Color = unname(status_col)), '_taxon_status_palette')
write_plot_data(data.table(Status = names(sv_specificity_col), Color = unname(sv_specificity_col)), '_sv_specificity_palette')
write_plot_data(data.table(BMI_SV_ID = rows$BMI_SV_ID, SGB = rows$SGB,
                           Taxon_abundance_status = rows$Taxon_abundance_status,
                           Current_robust_in_Ob = rows$Current_robust_in_Ob,
                           Current_robust_class = rows$Current_robust_class),
                '_taxon_status')

qc <- rows[, .N, by = .(Taxon_abundance_status)][
  data.table(Taxon_abundance_status = factor(status_levels, levels = status_levels)),
  on = 'Taxon_abundance_status'
]
qc[is.na(N), N := 0L]
qc[, `:=`(N_SV = nrow(rows), N_samples = ncol(plotmat),
          N_removed_all_NA_samples = nrow(removed_all_na),
          Style = STYLE,
          N_raw_endpoint_concordant = n_raw_endpoint_concordant,
          Width = width, Height = height)]
fwrite(qc, file.path(ROOT, 'tables', paste0(STEM, '_QC.tsv')), sep = '\t')

append_unique <- function(path, entries) {
  old <- if (file.exists(path)) readLines(path, warn = FALSE) else character()
  new <- setdiff(entries, old)
  if (length(new)) cat(paste0(paste(new, collapse = '\n'), '\n'), file = path, append = TRUE)
}
outputs <- c(
  list.files(file.path(ROOT, 'figures'), pattern = STEM, full.names = TRUE),
  list.files(file.path(ROOT, 'data_for_figures'), pattern = STEM, full.names = TRUE),
  list.files(file.path(ROOT, 'tables'), pattern = STEM, full.names = TRUE),
  file.path(ROOT, 'plot_ob_L1_vSV_clustered_taxon_SVspecificity.R')
)
append_unique(file.path(ROOT, 'output_files.txt'), outputs)
cat('Completed taxon-status heatmap:', paste(paste(qc$Taxon_abundance_status, qc$N, sep = '='), collapse = '; '), '\n')
