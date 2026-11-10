# ==============================================================================
# Script: 00_align_allBMI_SGV_metadata.R
#
# Purpose:
#   Match all-BMI SGVFinder sample IDs to master metadata
#   using the SAME ID harmonization logic as previous analyses.
#
# Mapping strategy:
#   1. Direct match:
#      SGV_ID == metadata$Sample.Source.ID.processed
#
#   2. SRA bridge:
#      SGV_ID (Run)
#        -> Sample_Name
#        -> Mapped_ID defined by BioProject-specific rules
#
#   3. Create metadata_allBMI_SGV_aligned
#      ordered exactly as the SGV matrices.
#
# IMPORTANT:
#   - SGV matrix IDs are NOT overwritten.
#   - Original matrix ID is stored as SGV_ID.
# ==============================================================================


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(purrr)
  library(fs)
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

D_FILE <- file.path(
  DATA_DIR,
  "21_allBMI_dSGV_crossBMI_features.csv.gz"
)

V_FILE <- file.path(
  DATA_DIR,
  "22_allBMI_vSGV_crossBMI_features.csv.gz"
)

META_FILE <- file.path(
  PROJECT,
  "02.processed_data",
  "01_metadata_processed",
  "metadata_processed_20260208_restored.rds"
)

SRA_DIR <- file.path(
  PROJECT,
  "00.rawdata",
  "03_SRA_data"
)

OUT_DIR <- DATA_DIR

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 2. Check files
# ==============================================================================

required_files <- c(
  D_FILE,
  V_FILE,
  META_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  cat("\nMissing files:\n")
  print(missing_files)
  
  stop("Required input files are missing.")
}


# ==============================================================================
# 3. Read ONLY SGV sample-ID columns
#
# No need to load all 165 features yet.
# ==============================================================================

message(">>> Reading all-BMI SGV sample IDs...")


d_id_dt <- fread(
  D_FILE,
  select = 1,
  colClasses = "character"
)

v_id_dt <- fread(
  V_FILE,
  select = 1,
  colClasses = "character"
)


d_id_col <- names(d_id_dt)[1]
v_id_col <- names(v_id_dt)[1]


d_ids <- as.character(
  d_id_dt[[1]]
)

v_ids <- as.character(
  v_id_dt[[1]]
)


cat(
  "\ndSGV ID column:",
  d_id_col,
  "\n"
)

cat(
  "vSGV ID column:",
  v_id_col,
  "\n"
)

cat(
  "\ndSGV samples:",
  length(d_ids),
  "\n"
)

cat(
  "vSGV samples:",
  length(v_ids),
  "\n"
)

cat(
  "Unique dSGV IDs:",
  uniqueN(d_ids),
  "\n"
)

cat(
  "Unique vSGV IDs:",
  uniqueN(v_ids),
  "\n"
)


# ==============================================================================
# 4. Basic SGV matrix ID QC
# ==============================================================================

if (anyDuplicated(d_ids)) {
  stop("Duplicated sample IDs found in dSGV matrix.")
}

if (anyDuplicated(v_ids)) {
  stop("Duplicated sample IDs found in vSGV matrix.")
}


if (!setequal(d_ids, v_ids)) {
  
  cat("\nIDs only in dSGV:\n")
  print(head(setdiff(d_ids, v_ids), 30))
  
  cat("\nIDs only in vSGV:\n")
  print(head(setdiff(v_ids, d_ids), 30))
  
  stop("dSGV and vSGV sample ID sets differ.")
}


cat(
  "\ndSGV/vSGV sample sets identical: YES\n"
)


# Use dSGV order as canonical SGV order
sgv_ids <- d_ids


# ==============================================================================
# 5. Read master metadata
# ==============================================================================

message("\n>>> Reading master metadata...")


meta_df <- readRDS(
  META_FILE
)

meta_df <- as.data.table(
  meta_df
)


target_id_col <- "Sample.Source.ID.processed"


if (!target_id_col %in% names(meta_df)) {
  
  stop(
    "Metadata does not contain: ",
    target_id_col
  )
}


target_meta_ids <- as.character(
  meta_df[[target_id_col]]
)


cat(
  "Metadata rows:",
  nrow(meta_df),
  "\n"
)

cat(
  "Unique target metadata IDs:",
  uniqueN(target_meta_ids),
  "\n"
)


# Target metadata ID should preferably be unique
dup_meta <- meta_df[
  duplicated(get(target_id_col)) |
    duplicated(get(target_id_col), fromLast = TRUE)
]


if (nrow(dup_meta) > 0) {
  
  cat(
    "\nWARNING: duplicated Sample.Source.ID.processed values:",
    uniqueN(dup_meta[[target_id_col]]),
    "\n"
  )
  
  print(
    unique(
      dup_meta[[target_id_col]]
    )[1:min(
      20,
      uniqueN(
        dup_meta[[target_id_col]]
      )
    )]
  )
}


# ==============================================================================
# 6. Build SRA bridge
#
# Same logic as previous ID mapping script.
# ==============================================================================

message(
  "\n>>> Building SRA Run -> Metadata bridge..."
)


sra_files <- dir_ls(
  SRA_DIR,
  glob = "*.csv"
)


cat(
  "SRA RunInfo files:",
  length(sra_files),
  "\n"
)


read_sra_info <- function(file_path) {
  
  dt <- tryCatch(
    fread(
      file_path,
      fill = TRUE,
      showProgress = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(dt)) {
    return(NULL)
  }
  
  
  required <- c(
    "Run",
    "BioProject",
    "Sample_Name"
  )
  
  
  if (!all(required %in% names(dt))) {
    return(NULL)
  }
  
  
  dt[
    ,
    ..required
  ]
}


sra_map_raw <- map_dfr(
  sra_files,
  read_sra_info
)


sra_map_raw <- as.data.table(
  sra_map_raw
)


# Keep one row per Run, same principle as old script
sra_map_raw <- unique(
  sra_map_raw,
  by = "Run"
)


cat(
  "Unique SRA Runs:",
  nrow(sra_map_raw),
  "\n"
)


# ==============================================================================
# 7. Apply previous special BioProject naming rules
# ==============================================================================

sra_bridge <- copy(
  sra_map_raw
)


sra_bridge[
  ,
  Mapped_ID := fcase(
    
    BioProject == "PRJEB27928",
    paste0(
      Run,
      "-21-0"
    ),
    
    BioProject == "PRJEB12449",
    paste0(
      Run,
      "-27-0-0"
    ),
    
    BioProject == "PRJNA447983",
    paste0(
      "CRC_MR_",
      Run,
      "_17"
    ),
    
    BioProject == "PRJEB7774",
    paste0(
      "SID",
      Run
    ),
    
    default = Run
  )
]


# ==============================================================================
# 8. Direct matches
# ==============================================================================

message(
  "\n>>> Performing ID matching..."
)


# ==============================================================================
# Create mapping table and normalize SGV IDs for matching
# ==============================================================================

# ==============================================================================
# 8. Create ID mapping table
# ==============================================================================

message("\n>>> Performing ID matching...")


id_map <- data.table(
  SGV_ID = as.character(sgv_ids)
)


# ------------------------------------------------------------
# Normalize SGV IDs ONLY for metadata matching
#
# Example:
# C011_reformat_1_kneaddata_paired
# ->
# C011
#
# Original SGV_ID is always retained.
# ------------------------------------------------------------

id_map[
  ,
  SGV_ID_processed := sub(
    "_reformat_1_kneaddata_paired$",
    "",
    SGV_ID
  )
]


# Initialize output columns
id_map[
  ,
  `:=`(
    Target_Meta_ID = NA_character_,
    Match_method = NA_character_
  )
]


# QC: confirm processed column exists
stopifnot(
  "SGV_ID_processed" %in% names(id_map)
)


cat(
  "\nIDs modified by suffix removal:",
  sum(
    id_map$SGV_ID !=
      id_map$SGV_ID_processed
  ),
  "\n"
)


print(
  head(
    id_map[
      SGV_ID != SGV_ID_processed,
      .(
        SGV_ID,
        SGV_ID_processed
      )
    ],
    20
  )
)


# ==============================================================================
# 9. Direct matching
# ==============================================================================

# ------------------------------------------------------------
# 9.1 Original SGV ID exact match
# ------------------------------------------------------------

direct_raw <- (
  id_map$SGV_ID %in%
    target_meta_ids
)


id_map[
  direct_raw,
  `:=`(
    Target_Meta_ID = SGV_ID,
    Match_method = "Direct_raw"
  )
]


# ------------------------------------------------------------
# 9.2 Exact match after suffix removal
# ------------------------------------------------------------

direct_processed <- (
  is.na(id_map$Target_Meta_ID) &
    id_map$SGV_ID_processed %in%
    target_meta_ids
)


id_map[
  direct_processed,
  `:=`(
    Target_Meta_ID =
      SGV_ID_processed,
    
    Match_method =
      "Direct_suffix_removed"
  )
]


cat(
  "\nDirect raw matches:",
  sum(direct_raw),
  "\n"
)

cat(
  "Matches after suffix removal:",
  sum(direct_processed),
  "\n"
)


# ==============================================================================
# 10. Remaining IDs for SRA bridge
# ==============================================================================

remaining_ids <- id_map[
  is.na(Target_Meta_ID),
  SGV_ID_processed
]


cat(
  "Remaining before SRA bridge:",
  length(remaining_ids),
  "\n"
)


bridge_candidates <- sra_bridge[
  Run %in% remaining_ids
]


bridge_candidates[
  ,
  Target_Match := fcase(
    
    Sample_Name %in% target_meta_ids,
    as.character(
      Sample_Name
    ),
    
    Mapped_ID %in% target_meta_ids,
    as.character(
      Mapped_ID
    ),
    
    default = NA_character_
  )
]


bridge_candidates[
  ,
  Match_method_bridge := fcase(
    
    Sample_Name %in% target_meta_ids,
    "SRA_Sample_Name",
    
    Mapped_ID %in% target_meta_ids,
    "SRA_Mapped_ID",
    
    default = NA_character_
  )
]


bridge_matched <- bridge_candidates[
  !is.na(Target_Match)
]


# One Run should map to one metadata ID
bridge_dup <- bridge_matched[
  ,
  .(
    N_target = uniqueN(
      Target_Match
    )
  ),
  by = Run
][
  N_target > 1
]


if (nrow(bridge_dup) > 0) {
  
  cat(
    "\nERROR: some Run IDs map to multiple metadata IDs:\n"
  )
  
  print(
    bridge_dup
  )
  
  stop(
    "Ambiguous SRA bridge mapping."
  )
}


bridge_matched <- unique(
  bridge_matched[
    ,
    .(
      Run,
      Target_Match,
      Match_method_bridge
    )
  ],
  by = "Run"
)


if (nrow(bridge_matched) > 0) {
  
  idx <- match(
    bridge_matched$Run,
    id_map$SGV_ID_processed
  )
  
  
  id_map[
    idx,
    `:=`(
      Target_Meta_ID =
        bridge_matched$Target_Match,
      
      Match_method =
        bridge_matched$Match_method_bridge
    )
  ]
}


cat(
  "SRA bridge matches:",
  sum(
    id_map$Match_method %in%
      c(
        "SRA_Sample_Name",
        "SRA_Mapped_ID"
      ),
    na.rm = TRUE
  ),
  "\n"
)


# ==============================================================================
# 10. Mapping QC
# ==============================================================================

id_map[
  ,
  Matched := !is.na(
    Target_Meta_ID
  )
]


mapping_summary <- id_map[
  ,
  .N,
  by = .(
    Match_method
  )
][
  order(-N)
]


cat(
  "\n========================================\n"
)

cat(
  "ALL-BMI SGV ID MAPPING SUMMARY\n"
)

cat(
  "========================================\n"
)


print(
  mapping_summary
)


cat(
  "\nTotal SGV samples:",
  nrow(id_map),
  "\n"
)

cat(
  "Successfully matched:",
  sum(id_map$Matched),
  "\n"
)

cat(
  "Unmatched:",
  sum(!id_map$Matched),
  "\n"
)


# ==============================================================================
# 11. Check whether multiple SGV IDs map to same metadata sample
# ==============================================================================

duplicate_target <- id_map[
  Matched == TRUE,
  .N,
  by = Target_Meta_ID
][
  N > 1
]


if (nrow(duplicate_target) > 0) {
  
  cat(
    "\nWARNING: multiple SGV IDs map to the same metadata ID:\n"
  )
  
  print(
    duplicate_target
  )
}


# ==============================================================================
# 12. Save mapping QC BEFORE stopping
# ==============================================================================

fwrite(
  id_map,
  file.path(
    OUT_DIR,
    "allBMI_SGV_ID_mapping.tsv"
  ),
  sep = "\t"
)


fwrite(
  mapping_summary,
  file.path(
    OUT_DIR,
    "allBMI_SGV_ID_mapping_summary.tsv"
  ),
  sep = "\t"
)


unmatched <- id_map[
  Matched == FALSE
]


fwrite(
  unmatched,
  file.path(
    OUT_DIR,
    "allBMI_SGV_ID_unmatched.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 13. Keep matched samples only
#
# Unmatched SGV samples are retained in QC output but excluded
# from downstream cross-BMI analyses.
# ==============================================================================

if (nrow(unmatched) > 0) {
  
  cat(
    "\nNOTE:",
    nrow(unmatched),
    "SGV samples have no corresponding metadata and will be excluded.\n"
  )
  
  cat(
    "Matched samples retained:",
    sum(id_map$Matched),
    "/",
    nrow(id_map),
    "\n"
  )
}


id_map_matched <- id_map[
  Matched == TRUE
]


cat(
  "\nProceeding with",
  nrow(id_map_matched),
  "matched SGV samples.\n"
)

# ==============================================================================
# 14. Construct aligned metadata for MATCHED samples
# ==============================================================================

message(
  "\n>>> Constructing aligned metadata for matched SGV samples..."
)


meta_index <- match(
  id_map_matched$Target_Meta_ID,
  target_meta_ids
)


if (anyNA(meta_index)) {
  
  stop(
    "Internal error: matched Target_Meta_ID missing from metadata."
  )
}


final_meta <- copy(
  meta_df[
    meta_index
  ]
)


# Original SGV matrix ID
final_meta[
  ,
  SGV_ID :=
    id_map_matched$SGV_ID
]


# Keep the processed matching ID for provenance
final_meta[
  ,
  SGV_ID_processed :=
    id_map_matched$SGV_ID_processed
]


final_meta[
  ,
  SGV_match_method :=
    id_map_matched$Match_method
]


setcolorder(
  final_meta,
  c(
    "SGV_ID",
    "SGV_ID_processed",
    "SGV_match_method",
    setdiff(
      names(final_meta),
      c(
        "SGV_ID",
        "SGV_ID_processed",
        "SGV_match_method"
      )
    )
  )
)


# ==============================================================================
# 15. Final metadata QC
# ==============================================================================

stopifnot(
  nrow(final_meta) ==
    nrow(id_map_matched)
)


stopifnot(
  all(
    final_meta$SGV_ID ==
      id_map_matched$SGV_ID
  )
)


stopifnot(
  all(
    final_meta[[target_id_col]] ==
      id_map_matched$Target_Meta_ID
  )
)


stopifnot(
  anyDuplicated(final_meta$SGV_ID) == 0
)


cat(
  "\n========================================\n"
)

cat(
  "FINAL METADATA ALIGNMENT PASSED\n"
)

cat(
  "========================================\n"
)

cat(
  "Original SGV samples:",
  length(sgv_ids),
  "\n"
)

cat(
  "Matched SGV samples:",
  nrow(final_meta),
  "\n"
)

cat(
  "Excluded unmatched samples:",
  length(sgv_ids) - nrow(final_meta),
  "\n"
)


# ==============================================================================
# 16. Subset all-BMI SGV matrices to matched samples
# ==============================================================================

message(
  "\n>>> Subsetting SGV matrices to matched samples..."
)


d_all <- fread(
  D_FILE,
  check.names = FALSE,
  colClasses = list(
    character = 1
  )
)

v_all <- fread(
  V_FILE,
  check.names = FALSE,
  colClasses = list(
    character = 1
  )
)


# First column = original SGV ID
setnames(
  d_all,
  1,
  "SGV_ID"
)

setnames(
  v_all,
  1,
  "SGV_ID"
)


matched_ids <- final_meta$SGV_ID


d_matched <- d_all[
  match(
    matched_ids,
    SGV_ID
  )
]

v_matched <- v_all[
  match(
    matched_ids,
    SGV_ID
  )
]


# Strict alignment check
stopifnot(
  all(
    d_matched$SGV_ID ==
      final_meta$SGV_ID
  )
)

stopifnot(
  all(
    v_matched$SGV_ID ==
      final_meta$SGV_ID
  )
)


cat(
  "Aligned dSGV samples:",
  nrow(d_matched),
  "\n"
)

cat(
  "Aligned vSGV samples:",
  nrow(v_matched),
  "\n"
)


# ==============================================================================
# 17. Save aligned all-BMI inputs
# ==============================================================================

saveRDS(
  final_meta,
  file.path(
    OUT_DIR,
    "metadata_allBMI_SGV_aligned.rds"
  )
)


fwrite(
  final_meta,
  file.path(
    OUT_DIR,
    "metadata_allBMI_SGV_aligned.csv"
  )
)


saveRDS(
  d_matched,
  file.path(
    OUT_DIR,
    "allBMI_dSGV_crossBMI_matched.rds"
  )
)


saveRDS(
  v_matched,
  file.path(
    OUT_DIR,
    "allBMI_vSGV_crossBMI_matched.rds"
  )
)


fwrite(
  d_matched,
  file.path(
    OUT_DIR,
    "allBMI_dSGV_crossBMI_matched.csv.gz"
  )
)


fwrite(
  v_matched,
  file.path(
    OUT_DIR,
    "allBMI_vSGV_crossBMI_matched.csv.gz"
  )
)

cat(
  "\nBMI group distribution:\n"
)

print(
  table(
    final_meta$BMI_Group1,
    useNA = "ifany"
  )
)


cat(
  "\nBMI group × Disease:\n"
)

print(
  table(
    final_meta$BMI_Group1,
    final_meta$Disease,
    useNA = "ifany"
  )
)