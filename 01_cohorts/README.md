# Cohort preparation and quality control

Sample filtering, metadata cleaning, regional BMI classification, sample-ID alignment and clinical summaries.

## Code

- [00_sample_filtering.R](01.scripts/00_sample_filtering.R)
- [01_metadata_preprocess.R](01.scripts/01_metadata_preprocess.R)
- [02_id_mapping_and_merge.R](01.scripts/02_id_mapping_and_merge.R)
- [03_metaphlan_split.R](01.scripts/03_metaphlan_split.R)
- [04_QC_Renormalization.R](01.scripts/04_QC_Renormalization.R)
- [run.py](analyses/2026-11-02_metadata-cohort-reconciliation/run.py)
- [run_cohort_adjusted.py](analyses/2026-11-02_metadata-cohort-reconciliation/run_cohort_adjusted.py)
- [plot_clinical_baseline.R](analyses/2026-11-02_metadata-cohort-reconciliation/plot_clinical_baseline.R)

## Execution notes

Review input files, project-root settings and output paths in the selected script. Required study data and intermediate results are not bundled. Scripts are separate workflows and should not be sourced consecutively.
