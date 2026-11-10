.libPaths(c('D:/R/library', .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
})

# Final targeted display of all Ob L1 dSV labels on the cohort-level retention heatmap.
# The earlier six-label version and source plotting scripts are preserved.
ROOT <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-09-29_BMI_enriched_SV_function_modelQC/2026-09-30_BMI_specificity_layers'
STYLE <- 'vsv_matched'
STEM <- 'Heatmap_L1_Ob_dSGV_cohort_retention_allSV_labels_final'
STAGES <- c('CTRL', 'AA', 'CRC')
STAGE_LABELS <- c(CTRL = 'CTRL', AA = 'Adenoma', CRC = 'CRC')

x <- readRDS(file.path(ROOT, 'prepared_data.rds'))
l1 <- as.data.table(x$l1)
models <- as.data.table(x$models)
meta <- as.data.table(x$meta)
deletion <- x$matrices[['dSGV']]
stopifnot(all(deletion[is.finite(deletion)] %in% c(0, 1)))
retention <- 1 - deletion

# Reuse the audited L1 row order, then restrict to Ob-enriched dSVs.
source_rows <- fread(file.path(ROOT, 'data_for_figures', 'Heatmap_L1_dSGV_allBMI_rows.tsv'))
rows <- source_rows[BMI_enriched_group == 'Ob']
setorder(rows, Row_order)
rows[, Row_order := seq_len(.N)]
stopifnot(nrow(rows) == 68L, all(rows$BMI_SV_ID %in% rownames(retention)))
stopifnot(!anyNA(rows$L2_nominal), sum(rows$L2_nominal) == 24L)
rows[, SV_specificity := factor(
  ifelse(L2_nominal, 'Ob-specific SV', 'Other Ob L1 SV'),
  levels = c('Ob-specific SV', 'Other Ob L1 SV'))]

# Region-based labels: show one interval, both intervals when there are two,
# and the segment count when a feature has more than two intervals.
rows[, Region_label_full := sub('^.*:', '', SV_feature)]
rows[, Region_label_full := gsub('_', '\u2013', Region_label_full, fixed = TRUE)]
rows[, Region_label_full := gsub(';', '; ', Region_label_full, fixed = TRUE)]
compact_region <- function(z) {
  segments <- trimws(strsplit(z, ';', fixed = TRUE)[[1]])
  if (length(segments) == 1L) return(segments)
  if (length(segments) == 2L) return(paste(segments, collapse = '; '))
  paste0(length(segments), ' segments')
}
rows[, N_segments := lengths(strsplit(Region_label_full, ';', fixed = TRUE))]
rows[, Region_label := vapply(Region_label_full, compact_region, character(1))]
rows[, Row_label := paste0(SGB_label, ' ', sub('^t__', '', SGB), ' : ', Region_label)]
rows[, Direction_short := ifelse(Target_Beta > 0, 'CRC: deletion higher', 'CRC: deletion lower')]

rows[, `:=`(Selected_row_label = TRUE, Display_row_label = Row_label,
            Label_selection_reason = 'All Ob L1 dSVs labeled; one interval shown in full, two intervals listed, >2 represented by segment count')]
stopifnot(nrow(rows) == 68L, all(nzchar(rows$Display_row_label)))

# Audited abundance status of the SGB carrying each dSV, matching the vSV figure.
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
status_levels <- c('Ob-specific', 'Ob-associated (non-specific)', 'Not significant in Ob')
rows[, Taxon_abundance_status := fcase(
  Current_robust_in_Ob %in% TRUE & Current_robust_class == 'Ob-specific', 'Ob-specific',
  Current_robust_in_Ob %in% TRUE, 'Ob-associated (non-specific)',
  default = 'Not significant in Ob'
)]
rows[, Taxon_abundance_status := factor(Taxon_abundance_status, levels = status_levels)]
status_col <- c(
  'Ob-specific' = '#7B3294',
  'Ob-associated (non-specific)' = '#008C95',
  'Not significant in Ob' = '#D9D9D9'
)
sv_specificity_col <- c('Ob-specific SV' = '#D55E00', 'Other Ob L1 SV' = '#D9D9D9')

keep <- which(meta$BMI_plot == 'Ob')
submeta <- copy(meta[keep])
submeta[, Study := as.character(Study_plot)]
submeta[, Stage := as.character(Disease_plot)]
stopifnot(all(submeta$Stage %in% STAGES))

cohort_display <- function(z) {
  z <- sub('^Public_study__', '', z)
  z <- sub('^This_study_cohort[0-9]+__ONCOBIOME_', '', z)
  z <- sub('^This_study_cohort[0-9]+__', '', z)
  z <- sub('^2021_YuJ$', 'YuJ_2021', z)
  z <- sub('^2023_Yunnan$', 'Yunnan_2023', z)
  z <- sub('^2025_HK$', 'HK_2025', z)
  z <- sub('^CHN_WF-CRC$', 'WF-CRC', z)
  z
}

# Preserve first appearance in metadata, while omitting cohorts absent from a stage.
cohort_order <- unique(submeta$Study)
aggregated <- list()
stage_matrices <- list()
stage_columns <- list()

for (stage in STAGES) {
  available <- cohort_order[cohort_order %in% unique(submeta[Stage == stage, Study])]
  mat <- matrix(
    NA_real_, nrow = nrow(rows), ncol = length(available),
    dimnames = list(rows$BMI_SV_ID, available)
  )
  column_info <- list()
  for (j in seq_along(available)) {
    cohort <- available[j]
    local_idx <- which(submeta$Stage == stage & submeta$Study == cohort)
    global_idx <- keep[local_idx]
    values <- retention[rows$BMI_SV_ID, global_idx, drop = FALSE]
    denominator <- rowSums(is.finite(values))
    numerator <- rowSums(values == 1, na.rm = TRUE)
    rate <- ifelse(denominator > 0, numerator / denominator, NA_real_)
    mat[, j] <- rate
    total_n <- length(global_idx)
    display <- cohort_display(cohort)
    column_info[[j]] <- data.table(
      Stage = stage,
      Stage_label = unname(STAGE_LABELS[stage]),
      Cohort = cohort,
      Cohort_display = display,
      Total_group_N = total_n,
      Column_label = paste0(display, ' (n=', total_n, ')'),
      Column_order_within_stage = j
    )
    aggregated[[length(aggregated) + 1L]] <- data.table(
      BMI_SV_ID = rows$BMI_SV_ID,
      SGB = rows$SGB,
      Stage = stage,
      Stage_label = unname(STAGE_LABELS[stage]),
      Cohort = cohort,
      Cohort_display = display,
      Total_group_N = total_n,
      Retained_N = numerator,
      Observed_N = denominator,
      Retention_rate = rate,
      Row_order = rows$Row_order,
      Column_order_within_stage = j
    )
  }
  stage_matrices[[stage]] <- mat
  stage_columns[[stage]] <- rbindlist(column_info)
}

long <- rbindlist(aggregated)
columns <- rbindlist(stage_columns)
stopifnot(all(long$Retention_rate[is.finite(long$Retention_rate)] >= 0 &
              long$Retention_rate[is.finite(long$Retention_rate)] <= 1))

# Verify that cohort numerators/denominators reproduce the pooled stage means exactly.
weighted_check <- long[, .(
  Cohort_weighted_rate = if (sum(Observed_N) > 0) sum(Retained_N) / sum(Observed_N) else NA_real_
), by = .(BMI_SV_ID, Stage)]
pooled_check <- rbindlist(lapply(STAGES, function(stage) {
  idx <- keep[submeta$Stage == stage]
  values <- retention[rows$BMI_SV_ID, idx, drop = FALSE]
  data.table(
    BMI_SV_ID = rows$BMI_SV_ID,
    Stage = stage,
    Pooled_rate = rowMeans(values, na.rm = TRUE)
  )
}))
pooled_check[!is.finite(Pooled_rate), Pooled_rate := NA_real_]
weighted_check <- merge(weighted_check, pooled_check, by = c('BMI_SV_ID', 'Stage'))
stopifnot(isTRUE(all.equal(weighted_check$Cohort_weighted_rate,
                           weighted_check$Pooled_rate,
                           tolerance = 1e-12, check.attributes = FALSE)))

# Low retention means more deletion and is red. The optional muted style uses
# only three anchors to reduce the multicolour mosaic effect.
if (STYLE == 'muted') {
  breaks <- c(0, 0.5, 1)
  colors <- c('#B85C6B', '#E7D9D2', '#587DA3')
  missing_col <- '#F7F7F7'
} else if (STYLE == 'vsv_matched') {
  breaks <- seq(0, 1, length.out = 7)
  vsv_fun <- colorRamp2(
    c(-1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2),
    c('#313695', '#4575B4', '#66C2A5', '#FFFFF0',
      '#FEE08B', '#F46D43', '#A50026')
  )
  # Match the vSV hue family but use its inner +/-0.8 colour range because
  # nearly half of observed dSV cohort cells sit exactly at retention 0 or 1.
  colors <- vsv_fun(seq(0.8, -0.8, length.out = 7))
  missing_col <- '#F0F0F0'
} else {
  breaks <- seq(0, 1, length.out = 7)
  colors <- c('#C3445B', '#EC7868', '#F2B764', '#FFF0B0', '#70C5A8', '#6098BF', '#7C8FBD')
  missing_col <- '#F7F7F7'
}
colorfun <- colorRamp2(breaks, colors)
row_split <- factor(rows$Direction_short,
                    levels = c('CRC: deletion higher', 'CRC: deletion lower'))
row_titles <- c('CRC: deletion higher' = 'Ob +', 'CRC: deletion lower' = 'Ob -')

make_heatmap <- function(stage, show_rows = FALSE) {
  mat <- stage_matrices[[stage]]
  col_info <- stage_columns[[stage]]
  stage_name <- unname(STAGE_LABELS[stage])
  stage_anno <- HeatmapAnnotation(
    Stage = rep(stage_name, ncol(mat)),
    col = list(Stage = c(CTRL = '#B8B8B8', Adenoma = '#F4A261', CRC = '#B94A48')),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    simple_anno_size = unit(2.5, 'mm')
  )
  Heatmap(
    mat,
    name = paste0('Retention_', stage),
    col = colorfun,
    na_col = missing_col,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cluster_row_slices = FALSE,
    show_column_names = FALSE,
    column_title = stage_name,
    column_title_gp = gpar(fontsize = 10, fontface = 'bold'),
    top_annotation = stage_anno,
    show_row_names = show_rows,
    row_labels = rows$Row_label,
    row_names_side = 'right',
    row_names_gp = gpar(fontsize = 5.4),
    row_names_max_width = unit(72, 'mm'),
    row_split = row_split,
    row_title = c('Ob+', 'Ob-'),
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 7),
    row_gap = unit(2.5, 'mm'),
    width = unit(ncol(mat) * 2.1, 'mm'),
    height = unit(nrow(rows) * 1.8, 'mm'),
    border = TRUE,
    rect_gp = gpar(col = NA),
    show_heatmap_legend = stage == 'CTRL',
    heatmap_legend_param = list(
      title = 'Retention rate',
      at = c(0, 0.25, 0.5, 0.75, 1),
      labels = c('0\n(more deleted)', '0.25', '0.50', '0.75', '1\n(retained)'),
      direction = 'horizontal',
      legend_width = unit(62, 'mm')
    )
  )
}

status_matrix <- matrix(
  as.character(rows$Taxon_abundance_status),
  ncol = 1,
  dimnames = list(rows$BMI_SV_ID, 'Taxon status')
)
status_heatmap <- Heatmap(
  status_matrix,
  name = 'Taxon abundance status',
  col = status_col,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cluster_row_slices = FALSE,
  row_split = row_split,
  row_title = c('Ob+', 'Ob-'),
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 7),
  show_column_names = FALSE,
  show_row_names = FALSE,
  row_names_side = 'right',
  row_names_gp = gpar(fontsize = 4.6),
  row_names_max_width = unit(65, 'mm'),
  width = unit(4, 'mm'),
  height = unit(nrow(rows) * 1.8, 'mm'),
  border = TRUE,
  heatmap_legend_param = list(
    title = 'Taxon abundance status',
    at = status_levels,
    labels = status_levels,
    nrow = 1,
    direction = 'horizontal'
  )
)
specificity_matrix <- matrix(as.character(rows$SV_specificity), ncol = 1,
                             dimnames = list(rows$BMI_SV_ID, 'SV specificity'))
specificity_heatmap <- Heatmap(
  specificity_matrix, name = 'SV specificity', col = sv_specificity_col,
  cluster_rows = FALSE, cluster_columns = FALSE, cluster_row_slices = FALSE,
  row_split = row_split, row_title = c('Ob+', 'Ob-'), row_title_rot = 0,
  row_title_gp = gpar(fontsize = 7),
  show_column_names = FALSE, show_row_names = TRUE,
  row_labels = rows$Display_row_label, row_names_side = 'right',
  row_names_gp = gpar(fontsize = 4.6), row_names_max_width = unit(65, 'mm'),
  width = unit(4, 'mm'), border = TRUE,
  height = unit(nrow(rows) * 1.8, 'mm'),
  heatmap_legend_param = list(title = 'SV specificity', at = names(sv_specificity_col),
    labels = names(sv_specificity_col), nrow = 1, direction = 'horizontal')
)

ht <- make_heatmap('CTRL') + make_heatmap('AA') + make_heatmap('CRC') + status_heatmap + specificity_heatmap
width <- 6.4
height <- 8.0
render <- function() {
  draw(
    ht,
    heatmap_legend_side = 'bottom',
    merge_legends = TRUE,
    ht_gap = unit(3, 'mm'),
    column_title = paste0('Ob L1 dSV cohort-level retention | ', nrow(rows), ' SVs'),
    column_title_gp = gpar(fontsize = 12, fontfamily = 'Arial'),
    padding = unit(c(8, 5, 10, 18), 'mm')
  )
}

figure_path <- file.path(ROOT, 'figures', STEM)
grDevices::cairo_pdf(paste0(figure_path, '.pdf'), width = width, height = height, family = 'Arial')
render()
dev.off()
png(paste0(figure_path, '.png'), width = width, height = height,
    units = 'in', res = 300, type = 'cairo')
render()
dev.off()
svglite::svglite(paste0(figure_path, '.svg'), width = width, height = height)
render()
dev.off()

write_tsv <- function(x, suffix = '', where = 'data_for_figures') {
  fwrite(x, file.path(ROOT, where, paste0(STEM, suffix, '.tsv')), sep = '\t', na = 'NA')
}
write_tsv(long, '_long')
write_tsv(columns, '_columns')
write_tsv(rows, '_rows')
write_tsv(data.table(Value = breaks, Color = colors), '_palette')
write_tsv(data.table(Status = names(status_col), Color = unname(status_col)), '_taxon_status_palette')
write_tsv(data.table(Status = names(sv_specificity_col), Color = unname(sv_specificity_col)), '_sv_specificity_palette')
write_tsv(rows[, .(BMI_SV_ID, SGB, Taxon_abundance_status,
                   Current_robust_in_Ob, Current_robust_class)], '_taxon_status')
write_tsv(rows[, .(BMI_SV_ID, SV_specificity, L2_nominal)], '_sv_specificity')
for (stage in STAGES) {
  mat <- stage_matrices[[stage]]
  write_tsv(data.table(BMI_SV_ID = rownames(mat), as.data.table(mat)), paste0('_', stage))
}

qc <- data.table(
  N_SV = nrow(rows),
  N_deletion_higher = sum(rows$Target_Beta > 0),
  N_deletion_lower = sum(rows$Target_Beta < 0),
  N_Ob_samples = length(keep),
  N_cohorts = length(cohort_order),
  N_CTRL_cohorts = nrow(stage_columns[['CTRL']]),
  N_Adenoma_cohorts = nrow(stage_columns[['AA']]),
  N_CRC_cohorts = nrow(stage_columns[['CRC']]),
  N_cells = nrow(long),
  N_cells_without_observed_values = sum(long$Observed_N == 0),
  N_cells_with_observed_N_lt_5 = sum(long$Observed_N > 0 & long$Observed_N < 5),
  Fraction_cells_with_observed_N_lt_5 = mean(long$Observed_N > 0 & long$Observed_N < 5),
  N_Ob_specific_SVs = sum(rows$Taxon_abundance_status == 'Ob-specific'),
  N_Ob_associated_non_specific_SVs = sum(rows$Taxon_abundance_status == 'Ob-associated (non-specific)'),
  N_not_significant_in_Ob_SVs = sum(rows$Taxon_abundance_status == 'Not significant in Ob'),
  Weighted_aggregation_check = TRUE,
  Style = STYLE,
  Width = width,
  Height = height
)
write_tsv(qc, '_QC', 'tables')

append_unique <- function(path, entries) {
  old <- if (file.exists(path)) readLines(path, warn = FALSE) else character()
  new <- setdiff(entries, old)
  if (length(new)) cat(paste0(paste(new, collapse = '\n'), '\n'), file = path, append = TRUE)
}
outputs <- c(
  list.files(file.path(ROOT, 'figures'), pattern = STEM, full.names = TRUE),
  list.files(file.path(ROOT, 'data_for_figures'), pattern = STEM, full.names = TRUE),
  list.files(file.path(ROOT, 'tables'), pattern = STEM, full.names = TRUE),
  file.path(ROOT, 'plot_ob_L1_dSV_allSV_labels_final.R')
)
append_unique(file.path(ROOT, 'output_files.txt'), outputs)
cat('Completed all-label Ob L1 dSV figure:', nrow(rows), 'SVs;',
    length(cohort_order), 'cohorts;', length(keep), 'samples\n')
