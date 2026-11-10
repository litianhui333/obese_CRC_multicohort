# Microbial structural variants

Structural-variant association models, species-abundance adjustment, cross-BMI comparisons and predicted gene-function annotations.

## Code

- [00_source_PRIMARY55_inputs.R](analyses/2026-09-18_SV-reanalysis-primary55/00_source_PRIMARY55_inputs.R)
- [02_obese_differential_SV_CTRL_CRC_v2.R](analyses/2026-09-18_SV-reanalysis-primary55/02_obese_differential_SV_CTRL_CRC_v2.R)
- [00_align_allBMI_SGV_metadata.R](analyses/2026-09-19_SV-crossBMI-analysis/00_align_allBMI_SGV_metadata.R)
- [00_source_crossBMI_inputs.R](analyses/2026-09-19_SV-crossBMI-analysis/00_source_crossBMI_inputs.R)
- [01_crossBMI_CTRL_vs_CRC.R](analyses/2026-09-19_SV-crossBMI-analysis/01_crossBMI_CTRL_vs_CRC.R)
- [02_compare_original205_crossBMI_specificity.R](analyses/2026-09-19_SV-crossBMI-analysis/02_compare_original205_crossBMI_specificity.R)
- [plot_callable_SGB_overlap.R](analyses/2026-09-28_allBMI_BMIstratified_SV_landscape/plot_callable_SGB_overlap.R)
- [run.R](analyses/2026-10-13_SV_abundance_sensitivity/run.R)
- [plot_Ob_heatmaps_abundanceAdjusted_20261026.R](analyses/2026-10-13_SV_abundance_sensitivity/plot_Ob_heatmaps_abundanceAdjusted_20261026.R)
- [plot_Ob_heatmaps_specificLabels_20261026.R](analyses/2026-10-13_SV_abundance_sensitivity/plot_Ob_heatmaps_specificLabels_20261026.R)
- [plot_Ob_specific_crossBMI_dotplot.R](analyses/2026-10-13_SV_abundance_sensitivity/plot_Ob_specific_crossBMI_dotplot.R)
- [plot_effect_specific_direction_counts.R](analyses/2026-10-13_SV_abundance_sensitivity/plot_effect_specific_direction_counts.R)
- [01_classify_BMI_SV_function_categories_modelQC_speciesResolved.R](analyses/2026-09-29_BMI_enriched_SV_function_modelQC/01_classify_BMI_SV_function_categories_modelQC_speciesResolved.R)
- [02_prepare_BMI_SV_circos_inputs_modelQC_speciesResolved_withSGB.R](analyses/2026-09-29_BMI_enriched_SV_function_modelQC/02_prepare_BMI_SV_circos_inputs_modelQC_speciesResolved_withSGB.R)
- [03_plot_BMI_SV_function_circos_modelQC_speciesResolved.R](analyses/2026-09-29_BMI_enriched_SV_function_modelQC/03_plot_BMI_SV_function_circos_modelQC_speciesResolved.R)
- [run.R](analyses/2026-10-27_function_SV_integration/run.R)

## Execution notes

Review input files, project-root settings and output paths in the selected script. Required study data and intermediate results are not bundled. Scripts are separate workflows and should not be sourced consecutively.
