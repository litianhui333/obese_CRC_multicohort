# Microbiome-based CRC classification

Species-only, structural-variant-only and combined feature panels evaluated with random forest, elastic net, radial SVM and XGBoost.

## Code

- [run_ob_associated_L1L2.R](analyses/2026-10-11_Ob_fixed_panel_xMarkerFinder_like/run_ob_associated_L1L2.R)
- [evaluate_fixed_panels_nested_repeated_CV.R](analyses/2026-10-11_Ob_fixed_panel_xMarkerFinder_like/evaluate_fixed_panels_nested_repeated_CV.R)
- [plot_fixed_panel_nestedCV_ROC_representative_repeat.R](analyses/2026-10-11_Ob_fixed_panel_xMarkerFinder_like/plot_fixed_panel_nestedCV_ROC_representative_repeat.R)
- [plot_fixed_panel_nestedCV_AUC.R](analyses/2026-10-11_Ob_fixed_panel_xMarkerFinder_like/plot_fixed_panel_nestedCV_AUC.R)
- [redraw_combined_marker_importance_taxon_labels.R](analyses/2026-10-11_Ob_fixed_panel_xMarkerFinder_like/redraw_combined_marker_importance_taxon_labels.R)
- [build_ob_panel_membership_supplement.R](analyses/2026-10-11_Ob_fixed_panel_xMarkerFinder_like/build_ob_panel_membership_supplement.R)
- [evaluate_locked_ob_panel_across_bmi.R](analyses/2026-10-11_Ob_fixed_panel_xMarkerFinder_like/evaluate_locked_ob_panel_across_bmi.R)
- [evaluate_ob_panel_crossBMI_other_algorithms.R](analyses/2026-10-11_Ob_fixed_panel_xMarkerFinder_like/evaluate_ob_panel_crossBMI_other_algorithms.R)
- [plot_crossBMI_otherAlgorithms_three_panels_by_method.R](analyses/2026-10-11_Ob_fixed_panel_xMarkerFinder_like/plot_crossBMI_otherAlgorithms_three_panels_by_method.R)
- [plot_joint_crossBMI_five_methods.R](analyses/2026-10-11_Ob_fixed_panel_xMarkerFinder_like/plot_joint_crossBMI_five_methods.R)
- [run.R](analyses/2026-10-27_Ob22_nestedCV_panel_tests/run.R)
- [run_multimethod_corrected_tests.R](analyses/2026-10-27_Ob22_nestedCV_panel_tests/run_multimethod_corrected_tests.R)
- [run.R](analyses/2026-10-13_Ob_L1L2_joint_LODO/run.R)
- [run.R](analyses/2026-10-20_Ob_three_fixed_panels_LODO/run.R)
- [plot_xgboost_lodo_roc.R](analyses/2026-10-20_Ob_three_fixed_panels_LODO/plot_xgboost_lodo_roc.R)
- [run.R](analyses/2026-10-27_Ob22_STST/run.R)
- [run_multimethod.R](analyses/2026-10-27_Ob22_STST/run_multimethod.R)
- [run.R](analyses/2026-10-12_Ob_fixed_joint_panel_four_comparisons/run.R)
- [plot_results.R](analyses/2026-10-12_Ob_fixed_joint_panel_four_comparisons/plot_results.R)
- [evaluate_ctrl_crc_model_transfer.R](analyses/2026-10-12_Ob_fixed_joint_panel_four_comparisons/evaluate_ctrl_crc_model_transfer.R)
- [plot_early_crc_saved_oof.R](analyses/2026-10-12_Ob_fixed_joint_panel_four_comparisons/plot_early_crc_saved_oof.R)

## Execution notes

Review input files, project-root settings and output paths in the selected script. Required study data and intermediate results are not bundled. Scripts are separate workflows and should not be sourced consecutively.

The evaluation scripts use fixed marker panels; marker discovery is not repeated in each evaluation fold.

