# ==============================================================================
# 00_source_crossBMI_inputs.R
#
# Purpose:
#   Load and QC the finalized all-BMI SGV inputs for cross-BMI analysis.
#
# Inputs:
#   - metadata_allBMI_SGV_aligned.rds
#   - allBMI_dSGV_crossBMI_matched.rds
#   - allBMI_vSGV_crossBMI_matched.rds
#   - 19_FINAL_all_comparable_obese_to_allBMI_map.tsv
#   - 20_FINAL_unique_allBMI_SVs_for_crossBMI.tsv
#
# Output objects:
#   meta_all
#   d_all
#   v_all
#   cross_map
#   joint_features
#   ROOT
#   TABLE_DIR
#   QC_DIR
# ==============================================================================


suppressPackageStartupMessages({
  library(data.table)
})


# ==============================================================================
# 1. Paths
# ==============================================================================

PROJECT <- "H:/My Drive/proj_ob_AA_CRC"

DATA_DIR <- file.path(
  PROJECT,
  "02.processed_data",
  "05_SV_result",
  "all_BMI_crossBMI"
)

ROOT <- file.path(
  PROJECT,
  "analyses",
  "2026-09-19_SV-crossBMI-analysis"
)

TABLE_DIR <- file.path(
  ROOT,
  "tables"
)

QC_DIR <- file.path(
  ROOT,
  "qc"
)

dir.create(
  ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  QC_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 2. Input files
# ==============================================================================

META_FILE <- file.path(
  DATA_DIR,
  "metadata_allBMI_SGV_aligned.rds"
)

D_FILE <- file.path(
  DATA_DIR,
  "allBMI_dSGV_crossBMI_matched.rds"
)

V_FILE <- file.path(
  DATA_DIR,
  "allBMI_vSGV_crossBMI_matched.rds"
)

MAP_FILE <- file.path(
  DATA_DIR,
  "19_FINAL_all_comparable_obese_to_allBMI_map.tsv"
)

FEATURE_FILE <- file.path(
  DATA_DIR,
  "20_FINAL_unique_allBMI_SVs_for_crossBMI.tsv"
)


required_files <- c(
  META_FILE,
  D_FILE,
  V_FILE,
  MAP_FILE,
  FEATURE_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  cat("\nMissing files:\n")
  print(missing_files)
  
  stop(
    "Cross-BMI source inputs are incomplete."
  )
}


# ==============================================================================
# 3. Read inputs
# ==============================================================================

meta_all <- as.data.table(
  readRDS(
    META_FILE
  )
)

d_all <- as.data.table(
  readRDS(
    D_FILE
  )
)

v_all <- as.data.table(
  readRDS(
    V_FILE
  )
)

cross_map <- fread(
  MAP_FILE,
  check.names = FALSE
)

joint_features <- fread(
  FEATURE_FILE,
  check.names = FALSE
)


# ==============================================================================
# 4. Required metadata columns
# ==============================================================================

required_meta_cols <- c(
  "SGV_ID",
  "Disease",
  "Age",
  "Sex",
  "Study.name",
  "BMI_Group1"
)

missing_meta_cols <- setdiff(
  required_meta_cols,
  names(meta_all)
)

if (length(missing_meta_cols) > 0) {
  
  stop(
    "Metadata missing columns: ",
    paste(
      missing_meta_cols,
      collapse = ", "
    )
  )
}


# ==============================================================================
# 5. Matrix ID QC
# ==============================================================================

if (!"SGV_ID" %in% names(d_all)) {
  setnames(
    d_all,
    1,
    "SGV_ID"
  )
}

if (!"SGV_ID" %in% names(v_all)) {
  setnames(
    v_all,
    1,
    "SGV_ID"
  )
}


d_all[
  ,
  SGV_ID := as.character(SGV_ID)
]

v_all[
  ,
  SGV_ID := as.character(SGV_ID)
]

meta_all[
  ,
  SGV_ID := as.character(SGV_ID)
]


stopifnot(
  anyDuplicated(
    meta_all$SGV_ID
  ) == 0
)

stopifnot(
  anyDuplicated(
    d_all$SGV_ID
  ) == 0
)

stopifnot(
  anyDuplicated(
    v_all$SGV_ID
  ) == 0
)


stopifnot(
  nrow(meta_all) ==
    nrow(d_all)
)

stopifnot(
  nrow(meta_all) ==
    nrow(v_all)
)


stopifnot(
  all(
    meta_all$SGV_ID ==
      d_all$SGV_ID
  )
)

stopifnot(
  all(
    meta_all$SGV_ID ==
      v_all$SGV_ID
  )
)


# ==============================================================================
# 6. Feature QC
# ==============================================================================

d_expected <- joint_features[
  SV_type == "dSGV",
  AllBMI_SV_feature
]

v_expected <- joint_features[
  SV_type == "vSGV",
  AllBMI_SV_feature
]


stopifnot(
  length(d_expected) == 124
)

stopifnot(
  length(v_expected) == 41
)

stopifnot(
  all(
    d_expected %in%
      names(d_all)
  )
)

stopifnot(
  all(
    v_expected %in%
      names(v_all)
  )
)


extra_d <- setdiff(
  names(d_all)[-1],
  d_expected
)

extra_v <- setdiff(
  names(v_all)[-1],
  v_expected
)

if (length(extra_d) > 0) {
  
  warning(
    "Extra dSGV features detected: ",
    length(extra_d)
  )
}

if (length(extra_v) > 0) {
  
  warning(
    "Extra vSGV features detected: ",
    length(extra_v)
  )
}


# ==============================================================================
# 7. Standardize metadata for modeling
# ==============================================================================

meta_all[
  ,
  Age_model :=
    suppressWarnings(
      as.numeric(Age)
    )
]

meta_all[
  ,
  Sex_model :=
    as.character(Sex)
]

meta_all[
  ,
  Study_model :=
    as.character(Study.name)
]

meta_all[
  ,
  BMI_model :=
    as.character(BMI_Group1)
]

meta_all[
  ,
  Disease_model :=
    as.character(Disease)
]


# ==============================================================================
# 8. BMI / Disease QC
# ==============================================================================

cat(
  "\n========================================\n"
)

cat(
  "CROSS-BMI INPUT QC\n"
)

cat(
  "========================================\n"
)


cat(
  "Aligned samples:",
  nrow(meta_all),
  "\n"
)

cat(
  "dSGV features:",
  length(d_expected),
  "\n"
)

cat(
  "vSGV features:",
  length(v_expected),
  "\n"
)

cat(
  "Joint features:",
  nrow(joint_features),
  "\n"
)

cat(
  "Mapping rows:",
  nrow(cross_map),
  "\n"
)


cat(
  "\nBMI × Disease:\n"
)

print(
  table(
    meta_all$BMI_model,
    meta_all$Disease_model,
    useNA = "ifany"
  )
)


# ==============================================================================
# 9. Save source QC
# ==============================================================================

source_qc <- data.table(
  Metric = c(
    "Aligned_samples",
    "dSGV_features",
    "vSGV_features",
    "Joint_features",
    "Cross_mapping_rows"
  ),
  
  N = c(
    nrow(meta_all),
    length(d_expected),
    length(v_expected),
    nrow(joint_features),
    nrow(cross_map)
  )
)

fwrite(
  source_qc,
  file.path(
    QC_DIR,
    "00_crossBMI_input_QC.tsv"
  ),
  sep = "\t"
)


cat(
  "\nCROSS-BMI INPUTS READY.\n"
)