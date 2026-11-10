# ============================================================
# Source file: PRIMARY55 SGV analysis inputs
#
# Purpose
#   - Define paths
#   - Load PRIMARY55 obese dSGV / vSGV matrices
#   - Load feature-to-SGB mapping
#   - Load aligned obese metadata
#   - FIX sample ID using SGV_ID
#   - Perform strict ID-matching QC
#   - Standardize metadata for downstream models
#   - Cache all processed inputs as RDS
#
# Usage in downstream scripts:
#
# source(
#   "H:/My Drive/proj_ob_AA_CRC/analyses/
#    2026-09-18_SV-reanalysis-primary55/
#    00_source_PRIMARY55_inputs.R"
# )
#
# To force rebuilding cache:
# REBUILD_SV_CACHE <- TRUE
# source(...)
# ============================================================


# ============================================================
# 0. Package
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})


# ============================================================
# 1. Project paths
# ============================================================

PROJECT <- "H:/My Drive/proj_ob_AA_CRC"

SV_DATA <- file.path(
  PROJECT,
  "02.processed_data",
  "05_SV_result"
)

ROOT <- file.path(
  PROJECT,
  "analyses",
  "2026-09-18_SV-reanalysis-primary55"
)

TABLE_DIR <- file.path(
  ROOT,
  "tables"
)

QC_DIR <- file.path(
  ROOT,
  "qc"
)

CACHE_DIR <- file.path(
  ROOT,
  "cache"
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

dir.create(
  CACHE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. Input files
# ============================================================

D_FILE <- file.path(
  TABLE_DIR,
  "02_obese_dSGV_PRIMARY55.csv.gz"
)

V_FILE <- file.path(
  TABLE_DIR,
  "03_obese_vSGV_PRIMARY55.csv.gz"
)

FEATURE_MAP_FILE <- file.path(
  TABLE_DIR,
  "01_obese_PRIMARY55_all_SV_feature_map.tsv"
)

META_FILE <- file.path(
  SV_DATA,
  "metadata_obese_SGV_aligned.rds"
)

input_files <- c(
  D_FILE,
  V_FILE,
  FEATURE_MAP_FILE,
  META_FILE
)

missing_files <- input_files[
  !file.exists(input_files)
]

if (length(missing_files) > 0) {
  
  cat(
    "\nMissing required input files:\n"
  )
  
  print(missing_files)
  
  stop(
    "PRIMARY55 source setup cannot continue."
  )
}


# ============================================================
# 3. Cache
# ============================================================

CACHE_FILE <- file.path(
  CACHE_DIR,
  "PRIMARY55_analysis_inputs.rds"
)

if (
  !exists(
    "REBUILD_SV_CACHE",
    inherits = FALSE
  )
) {
  
  REBUILD_SV_CACHE <- FALSE
}


# ============================================================
# 4. Input signature
#
# Used to ensure we do not accidentally load an outdated cache
# after changing upstream files.
# ============================================================

get_input_signature <- function(files) {
  
  info <- file.info(files)
  
  data.table(
    file = normalizePath(
      files,
      winslash = "/",
      mustWork = TRUE
    ),
    size = info$size,
    mtime = as.character(
      info$mtime
    )
  )
}


current_signature <- get_input_signature(
  input_files
)


# ============================================================
# 5. Decide whether cache can be used
# ============================================================

use_cache <- FALSE

if (
  file.exists(CACHE_FILE) &&
  !REBUILD_SV_CACHE
) {
  
  cached_object <- tryCatch(
    readRDS(CACHE_FILE),
    error = function(e) NULL
  )
  
  if (
    !is.null(cached_object) &&
    "input_signature" %in% names(cached_object)
  ) {
    
    old_sig <- cached_object$input_signature
    
    same_signature <- isTRUE(
      all.equal(
        old_sig,
        current_signature,
        check.attributes = FALSE
      )
    )
    
    if (same_signature) {
      use_cache <- TRUE
    }
  }
}


# ============================================================
# 6. Load cache if valid
# ============================================================

if (use_cache) {
  
  cat(
    "\n========================================\n"
  )
  
  cat(
    "Loading cached PRIMARY55 inputs\n"
  )
  
  cat(
    "========================================\n"
  )
  
  obj <- cached_object
  
  d <- obj$d
  v <- obj$v
  
  d_id <- obj$d_id
  v_id <- obj$v_id
  
  feature_map <- obj$feature_map
  
  meta <- obj$meta
  
  meta_sgv <- obj$meta_sgv
  meta_ctrl_crc <- obj$meta_ctrl_crc
  
  id_qc <- obj$id_qc
  
  rm(
    cached_object,
    obj
  )
  
  cat(
    "Cache loaded:\n",
    CACHE_FILE,
    "\n"
  )
  
  
  # ============================================================
  # 7. Otherwise rebuild inputs
  # ============================================================
  
} else {
  
  cat(
    "\n========================================\n"
  )
  
  cat(
    "Building PRIMARY55 analysis inputs\n"
  )
  
  cat(
    "========================================\n"
  )
  
  
  # ----------------------------------------------------------
  # 7.1 Read SGV matrices
  # ----------------------------------------------------------
  
  d <- fread(
    D_FILE,
    check.names = FALSE,
    na.strings = c(
      "",
      "NA",
      "NaN",
      "nan"
    )
  )
  
  v <- fread(
    V_FILE,
    check.names = FALSE,
    na.strings = c(
      "",
      "NA",
      "NaN",
      "nan"
    )
  )
  
  
  d_id <- as.character(
    d[[1]]
  )
  
  v_id <- as.character(
    v[[1]]
  )
  
  
  # ----------------------------------------------------------
  # 7.2 Basic matrix QC
  # ----------------------------------------------------------
  
  if (anyDuplicated(d_id)) {
    stop(
      "Duplicated sample IDs in dSGV matrix."
    )
  }
  
  if (anyDuplicated(v_id)) {
    stop(
      "Duplicated sample IDs in vSGV matrix."
    )
  }
  
  
  cat(
    "\ndSGV:",
    nrow(d),
    "samples x",
    ncol(d) - 1,
    "features\n"
  )
  
  cat(
    "vSGV:",
    nrow(v),
    "samples x",
    ncol(v) - 1,
    "features\n"
  )
  
  
  # Expected current PRIMARY55 matrices
  if (
    nrow(d) != 669 ||
    nrow(v) != 669
  ) {
    
    warning(
      paste0(
        "Expected 669 obese samples, but observed ",
        nrow(d),
        " dSGV and ",
        nrow(v),
        " vSGV samples."
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # 7.3 Feature map
  # ----------------------------------------------------------
  
  feature_map <- fread(
    FEATURE_MAP_FILE,
    check.names = FALSE
  )
  
  
  required_feature_cols <- c(
    "SV_type",
    "SV_feature",
    "Reference",
    "SGB"
  )
  
  if (
    !all(
      required_feature_cols %in%
      names(feature_map)
    )
  ) {
    
    stop(
      "Feature map is missing required columns."
    )
  }
  
  
  # ----------------------------------------------------------
  # 7.4 Read aligned metadata
  # ----------------------------------------------------------
  
  meta <- as.data.table(
    readRDS(
      META_FILE
    )
  )
  
  
  required_meta_cols <- c(
    "SGV_ID",
    "Disease",
    "Age",
    "Sex",
    "Study.name"
  )
  
  missing_meta_cols <- setdiff(
    required_meta_cols,
    names(meta)
  )
  
  if (
    length(missing_meta_cols) > 0
  ) {
    
    stop(
      "Metadata missing columns: ",
      paste(
        missing_meta_cols,
        collapse = ", "
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # 7.5 Explicit ID matching QC
  #
  # IMPORTANT:
  # SGV_ID is the harmonized ID used for SGV matrices.
  # Name / Sample.Source.ID are NOT used for model joins.
  # ----------------------------------------------------------
  
  candidate_id_cols <- c(
    "SGV_ID",
    "Name",
    "Sample.Source.ID",
    "Sample.Source.ID.processed"
  )
  
  candidate_id_cols <-
    candidate_id_cols[
      candidate_id_cols %in%
        names(meta)
    ]
  
  
  id_qc <- rbindlist(
    lapply(
      candidate_id_cols,
      function(col) {
        
        x <- as.character(
          meta[[col]]
        )
        
        data.table(
          Metadata_column = col,
          
          N_metadata =
            length(x),
          
          N_unique_metadata =
            uniqueN(x),
          
          Match_dSGV =
            sum(
              d_id %in% x
            ),
          
          Match_vSGV =
            sum(
              v_id %in% x
            ),
          
          Metadata_match_dSGV =
            sum(
              x %in% d_id
            ),
          
          Metadata_match_vSGV =
            sum(
              x %in% v_id
            )
        )
      }
    )
  )
  
  
  cat(
    "\n========================================\n"
  )
  
  cat(
    "METADATA ID MATCHING QC\n"
  )
  
  cat(
    "========================================\n"
  )
  
  print(
    id_qc
  )
  
  
  fwrite(
    id_qc,
    file.path(
      QC_DIR,
      "01_metadata_SGV_ID_matching_QC.tsv"
    ),
    sep = "\t"
  )
  
  
  # ----------------------------------------------------------
  # 7.6 SGV_ID MUST be exact 1-to-1 mapping
  # ----------------------------------------------------------
  
  if (
    anyDuplicated(
      meta$SGV_ID
    )
  ) {
    
    stop(
      "Duplicated SGV_ID in metadata."
    )
  }
  
  
  if (
    !setequal(
      d_id,
      as.character(meta$SGV_ID)
    )
  ) {
    
    stop(
      paste0(
        "dSGV sample IDs are not identical ",
        "to metadata SGV_ID."
      )
    )
  }
  
  
  if (
    !setequal(
      v_id,
      as.character(meta$SGV_ID)
    )
  ) {
    
    stop(
      paste0(
        "vSGV sample IDs are not identical ",
        "to metadata SGV_ID."
      )
    )
  }
  
  
  if (
    !setequal(
      d_id,
      v_id
    )
  ) {
    
    stop(
      "dSGV and vSGV sample ID sets differ."
    )
  }
  
  
  cat(
    "\nSGV_ID exact matching confirmed:\n"
  )
  
  cat(
    "dSGV:",
    sum(
      d_id %in%
        meta$SGV_ID
    ),
    "/",
    length(d_id),
    "\n"
  )
  
  cat(
    "vSGV:",
    sum(
      v_id %in%
        meta$SGV_ID
    ),
    "/",
    length(v_id),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # 7.7 Standardized metadata
  #
  # SGV_ID is now permanently used as Sample_ID.
  # ----------------------------------------------------------
  
  meta_sgv <- meta[
    ,
    .(
      Sample_ID =
        as.character(
          SGV_ID
        ),
      
      Disease_original =
        as.character(
          Disease
        ),
      
      Age =
        suppressWarnings(
          as.numeric(
            Age
          )
        ),
      
      Sex =
        as.character(
          Sex
        ),
      
      Study =
        as.character(
          Study.name
        ),
      
      BMI =
        if (
          "BMI" %in% names(meta)
        )
          suppressWarnings(
            as.numeric(BMI)
          )
      else
        NA_real_,
      
      BMI_Group =
        if (
          "BMI_Group1" %in%
          names(meta)
        )
          as.character(
            BMI_Group1
          )
      else
        NA_character_
    )
  ]
  
  
  # ----------------------------------------------------------
  # 7.8 CTRL / CRC disease variable
  #
  # Adenoma samples remain in meta_sgv and are NOT discarded.
  # meta_ctrl_crc is specifically for current discovery model.
  # ----------------------------------------------------------
  
  disease_upper <- toupper(
    trimws(
      meta_sgv$Disease_original
    )
  )
  
  
  meta_sgv[
    ,
    Disease_CTRL_CRC :=
      fifelse(
        disease_upper %in% c(
          "CTRL",
          "CONTROL",
          "HC",
          "HEALTHY",
          "HEALTHY CONTROL"
        ),
        "CTRL",
        fifelse(
          disease_upper == "CRC",
          "CRC",
          NA_character_
        )
      )
  ]
  
  
  meta_ctrl_crc <- meta_sgv[
    Disease_CTRL_CRC %in%
      c(
        "CTRL",
        "CRC"
      )
  ]
  
  
  meta_ctrl_crc[
    ,
    Disease :=
      factor(
        Disease_CTRL_CRC,
        levels = c(
          "CTRL",
          "CRC"
        )
      )
  ]
  
  
  meta_ctrl_crc[
    ,
    Sex := factor(Sex)
  ]
  
  
  meta_ctrl_crc[
    ,
    Study := factor(Study)
  ]
  
  
  # ----------------------------------------------------------
  # 7.9 Metadata summary
  # ----------------------------------------------------------
  
  cat(
    "\n========================================\n"
  )
  
  cat(
    "STANDARDIZED METADATA\n"
  )
  
  cat(
    "========================================\n"
  )
  
  
  cat(
    "All obese SGV samples:",
    nrow(meta_sgv),
    "\n"
  )
  
  
  cat(
    "CTRL/CRC samples:",
    nrow(meta_ctrl_crc),
    "\n"
  )
  
  
  print(
    table(
      meta_ctrl_crc$Disease,
      useNA = "ifany"
    )
  )
  
  
  cat(
    "\nStudies:",
    uniqueN(
      meta_ctrl_crc$Study
    ),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # 7.10 Save standardized metadata
  # ----------------------------------------------------------
  
  fwrite(
    meta_sgv,
    file.path(
      TABLE_DIR,
      "00_metadata_obese_SGV_standardized.tsv"
    ),
    sep = "\t"
  )
  
  
  fwrite(
    meta_ctrl_crc,
    file.path(
      TABLE_DIR,
      "00_metadata_obese_CTRL_CRC_standardized.tsv"
    ),
    sep = "\t"
  )
  
  
  # ----------------------------------------------------------
  # 7.11 Save cache
  # ----------------------------------------------------------
  
  cache_object <- list(
    
    input_signature =
      current_signature,
    
    d = d,
    v = v,
    
    d_id = d_id,
    v_id = v_id,
    
    feature_map =
      feature_map,
    
    meta = meta,
    
    meta_sgv =
      meta_sgv,
    
    meta_ctrl_crc =
      meta_ctrl_crc,
    
    id_qc =
      id_qc
  )
  
  
  saveRDS(
    cache_object,
    CACHE_FILE,
    compress = TRUE
  )
  
  
  cat(
    "\nPRIMARY55 input cache written:\n",
    CACHE_FILE,
    "\n"
  )
  
  
  rm(
    cache_object
  )
}


# ============================================================
# 8. Aliases used by downstream scripts
# ============================================================

# Current primary discovery metadata
m <- copy(
  meta_ctrl_crc
)


# ============================================================
# 9. Final source report
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "PRIMARY55 INPUTS READY\n"
)

cat(
  "========================================\n"
)

cat(
  "dSGV:",
  nrow(d),
  "x",
  ncol(d) - 1,
  "\n"
)

cat(
  "vSGV:",
  nrow(v),
  "x",
  ncol(v) - 1,
  "\n"
)

cat(
  "Metadata SGV samples:",
  nrow(meta_sgv),
  "\n"
)

cat(
  "CTRL/CRC metadata samples:",
  nrow(m),
  "\n"
)

cat(
  "Sample ID key: SGV_ID\n"
)

cat(
  "Cache used:",
  use_cache,
  "\n"
)