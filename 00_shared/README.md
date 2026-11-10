# Shared analysis utilities

Reusable functions for differential abundance, robustness checks, feature annotation and downstream plotting.

## Code

- [DA_function.R](01.scripts/DA_function.R)
- [DA_analysis_figure_utils_functions.R](01.scripts/DA_analysis_figure_utils_functions.R)
- [get_robust_features_with_coef_check.R](01.scripts/get_robust_features_with_coef_check.R)
- [alpha_analysis_functions.R](01.scripts/alpha_analysis_functions.R)
- [plot_ob_L1_dSV_allSV_labels_final.R](analyses/2026-09-29_BMI_enriched_SV_function_modelQC/2026-09-30_BMI_specificity_layers/plot_ob_L1_dSV_allSV_labels_final.R)
- [plot_ob_L1_vSV_clustered_taxon_SVspecificity.R](analyses/2026-09-29_BMI_enriched_SV_function_modelQC/2026-09-30_BMI_specificity_layers/plot_ob_L1_vSV_clustered_taxon_SVspecificity.R)
- [plot_within_study_transfer_four_panels.R](analyses/2026-10-12_Ob_fixed_joint_panel_four_comparisons/plot_within_study_transfer_four_panels.R)

## Execution notes

Review input files, project-root settings and output paths in the selected script. Required study data and intermediate results are not bundled. Scripts are separate workflows and should not be sourced consecutively.
