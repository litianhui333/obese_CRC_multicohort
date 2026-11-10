# ==============================================================================
# 02_compare_original205_crossBMI_specificity.R
#
# Purpose
# -------
# Compare the original 205 obese-discovered CRC-associated SV candidates with
# their cross-BMI results in the jointly defined all-BMI SV feature space.
#
# STRICT cross-BMI pattern:
#   Ob P < 0.01
#   Nw P >= 0.05
#   Ow P >= 0.05
#
# RELAXED cross-BMI pattern:
#   Ob P < 0.05
#   Nw P >= 0.05
#   Ow P >= 0.05
#
# Supported pattern:
#   cross-BMI pattern + discovery/cross-BMI obese effect direction concordance
#
# IMPORTANT
# ---------
# These patterns support an "obesity-enriched" or "predominantly observed in
# obese individuals" interpretation. They are not, by themselves, formal
# Disease x BMI interaction tests.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ==============================================================================
# 1. Paths
# ==============================================================================

PROJECT <- "H:/My Drive/proj_ob_AA_CRC"

DISCOVERY_DIR <- file.path(
  PROJECT,
  "analyses",
  "2026-09-18_SV-reanalysis-primary55"
)

CROSS_DIR <- file.path(
  PROJECT,
  "analyses",
  "2026-09-19_SV-crossBMI-analysis"
)

DISCOVERY_FILE <- file.path(
  DISCOVERY_DIR,
  "tables",
  "12_obese_candidate205_master.tsv"
)

CROSS_FILE <- file.path(
  CROSS_DIR,
  "tables",
  "02_crossBMI_CTRL_vs_CRC_wide.tsv"
)

OUT_DIR <- file.path(
  CROSS_DIR,
  "tables"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ==============================================================================
# 2. Check inputs
# ==============================================================================

if (!file.exists(DISCOVERY_FILE)) {
  stop(
    "Cannot find original obese candidate master file:\n",
    DISCOVERY_FILE
  )
}

if (!file.exists(CROSS_FILE)) {
  stop(
    "Cannot find cross-BMI result file:\n",
    CROSS_FILE
  )
}

# ==============================================================================
# 3. Read data
# ==============================================================================

discovery <- fread(
  DISCOVERY_FILE,
  check.names = FALSE
)

cross <- fread(
  CROSS_FILE,
  check.names = FALSE
)

cat(
  "\n========================================\n",
  "INPUT SUMMARY\n",
  "========================================\n",
  sep = ""
)

cat(
  "Original obese candidates:",
  nrow(discovery),
  "\n"
)

cat(
  "Unique all-BMI SV units:",
  nrow(cross),
  "\n"
)

# ==============================================================================
# 4. Basic input QC
# ==============================================================================

required_discovery <- c(
  "Candidate_ID"
)

required_cross <- c(
  "Candidate_IDs",
  "Joint_SV_ID",
  "SGB",
  "SV_type",
  "SV_feature",
  "Beta_Ob",
  "P_Ob",
  "Beta_Ow",
  "P_Ow",
  "Beta_Nw",
  "P_Nw"
)

missing_discovery <- setdiff(
  required_discovery,
  names(discovery)
)

missing_cross <- setdiff(
  required_cross,
  names(cross)
)

if (length(missing_discovery) > 0) {
  stop(
    "Discovery file missing columns: ",
    paste(missing_discovery, collapse = ", ")
  )
}

if (length(missing_cross) > 0) {
  stop(
    "Cross-BMI file missing columns: ",
    paste(missing_cross, collapse = ", ")
  )
}

if (anyDuplicated(discovery$Candidate_ID) > 0) {
  stop(
    "Candidate_ID is not unique in the discovery master table."
  )
}

# ==============================================================================
# 5. Standardize original discovery columns
# ==============================================================================

rename_if_present <- function(dt, old, new) {
  if (
    old %in% names(dt) &&
    !(new %in% names(dt))
  ) {
    setnames(dt, old, new)
  }

  invisible(dt)
}

rename_if_present(discovery, "Beta", "Beta_discovery")
rename_if_present(discovery, "SE", "SE_discovery")
rename_if_present(discovery, "P", "P_discovery")
rename_if_present(discovery, "Q_type", "Q_type_discovery")
rename_if_present(discovery, "Model_status", "Model_status_discovery")
rename_if_present(discovery, "Singular", "Singular_discovery")
rename_if_present(discovery, "Optimizer_OK", "Optimizer_OK_discovery")

if (
  all(
    c(
      "P_discovery",
      "Q_type_discovery"
    ) %in% names(discovery)
  )
) {
  discovery[
    ,
    Discovery_candidate_criterion :=
      (!is.na(P_discovery) & P_discovery < 0.01) |
      (!is.na(Q_type_discovery) & Q_type_discovery < 0.10)
  ]
} else {
  discovery[
    ,
    Discovery_candidate_criterion := TRUE
  ]

  warning(
    "P_discovery and/or Q_type_discovery not found. ",
    "All rows are assumed to be original discovery candidates."
  )
}

# ==============================================================================
# 6. Expand Candidate_IDs from the all-BMI feature table
#
# Needed because some mappings are split/merged.
# ==============================================================================

cross[
  ,
  Cross_row_ID := .I
]

cross_long <- cross[
  ,
  {
    ids <- trimws(
      unlist(
        strsplit(
          as.character(Candidate_IDs),
          ";",
          fixed = TRUE
        )
      )
    )

    ids <- ids[
      !is.na(ids) &
      ids != ""
    ]

    .(
      Candidate_ID = ids
    )
  },
  by = names(cross)
]

cat(
  "Original candidates represented in cross-BMI map:",
  uniqueN(cross_long$Candidate_ID),
  "\n"
)

cat(
  "Expanded candidate-to-allBMI mapping rows:",
  nrow(cross_long),
  "\n"
)

# ==============================================================================
# 7. Define cross-BMI patterns
# ==============================================================================

cross_long[
  ,
  CrossBMI_valid :=
    !is.na(P_Ob) &
    !is.na(P_Ow) &
    !is.na(P_Nw)
]

cross_long[
  ,
  Cross_strict :=
    CrossBMI_valid &
    P_Ob < 0.01 &
    P_Ow >= 0.05 &
    P_Nw >= 0.05
]

cross_long[
  ,
  Cross_relaxed :=
    CrossBMI_valid &
    P_Ob < 0.05 &
    P_Ow >= 0.05 &
    P_Nw >= 0.05
]

# ==============================================================================
# 8. Merge original obese discovery with cross-BMI results
# ==============================================================================

master <- merge(
  discovery,
  cross_long,
  by = "Candidate_ID",
  all.x = TRUE,
  sort = FALSE,
  suffixes = c(
    "_discovery_table",
    "_crossBMI"
  )
)

# ==============================================================================
# 9. Candidate-level comparability
# ==============================================================================

master[
  ,
  CrossBMI_comparable :=
    !is.na(Joint_SV_ID) &
    CrossBMI_valid == TRUE
]

# ==============================================================================
# 10. Discovery vs cross-BMI obese effect concordance
# ==============================================================================

if ("Beta_discovery" %in% names(master)) {
  master[
    ,
    Discovery_crossOb_same_direction :=
      CrossBMI_comparable &
      !is.na(Beta_discovery) &
      !is.na(Beta_Ob) &
      sign(Beta_discovery) == sign(Beta_Ob)
  ]
} else {
  master[
    ,
    Discovery_crossOb_same_direction := NA
  ]

  warning(
    "Beta_discovery not available; directional concordance cannot be assessed."
  )
}

# ==============================================================================
# 11. Effect-size annotations
#
# These are annotations, not hard filters for the main strict/relaxed tables.
# ==============================================================================

master[
  ,
  Abs_Beta_Ob := abs(Beta_Ob)
]

master[
  ,
  Abs_Beta_Ow := abs(Beta_Ow)
]

master[
  ,
  Abs_Beta_Nw := abs(Beta_Nw)
]

master[
  ,
  SameDir_Nw :=
    CrossBMI_comparable &
    !is.na(Beta_Ob) &
    !is.na(Beta_Nw) &
    sign(Beta_Ob) == sign(Beta_Nw)
]

master[
  ,
  SameDir_Ow :=
    CrossBMI_comparable &
    !is.na(Beta_Ob) &
    !is.na(Beta_Ow) &
    sign(Beta_Ob) == sign(Beta_Ow)
]


master[
  ,
  Ob_stronger_than_sameDir_Nw :=
    CrossBMI_comparable &
    !is.na(Beta_Ob) &
    !is.na(Beta_Nw) &
    (
      SameDir_Nw == FALSE |
        abs(Beta_Ob) > abs(Beta_Nw)
    )
]


master[
  ,
  Ob_stronger_than_sameDir_Ow :=
    CrossBMI_comparable &
    !is.na(Beta_Ob) &
    !is.na(Beta_Ow) &
    (
      SameDir_Ow == FALSE |
        abs(Beta_Ob) > abs(Beta_Ow)
    )
]


master[
  ,
  Ob_stronger_than_sameDir_nonOb :=
    Ob_stronger_than_sameDir_Nw &
    Ob_stronger_than_sameDir_Ow
]


master[
  ,
  Any_nonOb_opposite_direction :=
    CrossBMI_comparable &
    (
      SameDir_Nw == FALSE |
        SameDir_Ow == FALSE
    )
]

master[
  ,
  Max_nonOb_absBeta :=
    fifelse(
      !is.na(Abs_Beta_Nw) &
      !is.na(Abs_Beta_Ow),
      pmax(
        Abs_Beta_Nw,
        Abs_Beta_Ow
      ),
      NA_real_
    )
]

master[
  ,
  Ob_vs_nonOb_effect_ratio :=
    fifelse(
      !is.na(Max_nonOb_absBeta) &
      Max_nonOb_absBeta > 0,
      Abs_Beta_Ob / Max_nonOb_absBeta,
      NA_real_
    )
]

master[
  ,
  Ob_Nw_same_direction :=
    CrossBMI_comparable &
    !is.na(Beta_Ob) &
    !is.na(Beta_Nw) &
    sign(Beta_Ob) == sign(Beta_Nw)
]

master[
  ,
  Ob_Ow_same_direction :=
    CrossBMI_comparable &
    !is.na(Beta_Ob) &
    !is.na(Beta_Ow) &
    sign(Beta_Ob) == sign(Beta_Ow)
]

# ==============================================================================
# 12. Define supported patterns
#
# STRICT support:
#   Ob P < 0.01; Nw/Ow P >= 0.05; obese Beta direction concordant
#
# RELAXED support:
#   Ob P < 0.05; Nw/Ow P >= 0.05; obese Beta direction concordant
# ==============================================================================

master[
  ,
  Strict_obesity_enriched_support :=
    CrossBMI_comparable &
    Cross_strict == TRUE &
    Discovery_crossOb_same_direction == TRUE
]

master[
  ,
  Relaxed_obesity_enriched_support :=
    CrossBMI_comparable &
    Cross_relaxed == TRUE &
    Discovery_crossOb_same_direction == TRUE
]

# ==============================================================================
# 13. Candidate-level summary for split mappings
#
# A split obese candidate may map to multiple all-BMI components. The full master
# table keeps each component; this table summarizes one row per Candidate_ID.
# ==============================================================================

candidate_summary <- master[
  ,
  .(
    N_crossBMI_rows = sum(!is.na(Joint_SV_ID)),

    N_valid_crossBMI_rows = sum(
      CrossBMI_comparable,
      na.rm = TRUE
    ),

    CrossBMI_comparable = any(
      CrossBMI_comparable,
      na.rm = TRUE
    ),

    Any_strict_pattern = any(
      Cross_strict == TRUE,
      na.rm = TRUE
    ),

    Any_relaxed_pattern = any(
      Cross_relaxed == TRUE,
      na.rm = TRUE
    ),

    Any_direction_concordant = any(
      Discovery_crossOb_same_direction == TRUE,
      na.rm = TRUE
    ),

    Strict_obesity_enriched_support = any(
      Strict_obesity_enriched_support == TRUE,
      na.rm = TRUE
    ),

    Relaxed_obesity_enriched_support = any(
      Relaxed_obesity_enriched_support == TRUE,
      na.rm = TRUE
    ),

    Any_Ob_stronger_than_both = any(
      Ob_stronger_than_sameDir_nonOb == TRUE,
      na.rm = TRUE
    ),

    Best_P_Ob =
      if (all(is.na(P_Ob))) {
        NA_real_
      } else {
        min(P_Ob, na.rm = TRUE)
      },

    Best_effect_ratio =
      if (all(is.na(Ob_vs_nonOb_effect_ratio))) {
        NA_real_
      } else {
        max(
          Ob_vs_nonOb_effect_ratio,
          na.rm = TRUE
        )
      }
  ),
  by = Candidate_ID
]

# Add selected original discovery information.
discovery_cols_for_summary <- intersect(
  c(
    "Candidate_ID",
    "SGB",
    "SV_type",
    "SV_feature",
    "Reference",
    "Beta_discovery",
    "SE_discovery",
    "P_discovery",
    "Q_type_discovery",
    "Model_status_discovery",
    "Singular_discovery",
    "Optimizer_OK_discovery",
    "Discovery_candidate_criterion"
  ),
  names(discovery)
)

candidate_summary <- merge(
  discovery[
    ,
    ..discovery_cols_for_summary
  ],
  candidate_summary,
  by = "Candidate_ID",
  all.x = TRUE,
  sort = FALSE
)

candidate_summary[
  is.na(CrossBMI_comparable),
  CrossBMI_comparable := FALSE
]

# ==============================================================================
# 14. Discovery-vs-crossBMI obese effect correlation
#
# Descriptive only because obese samples substantially overlap across analyses.
# ==============================================================================

effect_cor <- NA_real_

if ("Beta_discovery" %in% names(master)) {
  cor_dt <- master[
    CrossBMI_comparable == TRUE &
    !is.na(Beta_discovery) &
    !is.na(Beta_Ob)
  ]

  if (nrow(cor_dt) >= 3) {
    effect_cor <- suppressWarnings(
      cor(
        cor_dt$Beta_discovery,
        cor_dt$Beta_Ob,
        use = "complete.obs",
        method = "spearman"
      )
    )
  }
}

# ==============================================================================
# 15. Main output tables
# ==============================================================================

# A. Complete row-level master table
fwrite(
  master,
  file.path(
    OUT_DIR,
    "08_original205_crossBMI_master.tsv"
  ),
  sep = "\t"
)

# B. Candidate-level summary, one row per original Candidate_ID
fwrite(
  candidate_summary,
  file.path(
    OUT_DIR,
    "09_original205_crossBMI_candidate_summary.tsv"
  ),
  sep = "\t"
)

# C. Strict supported rows
strict_supported <- master[
  Strict_obesity_enriched_support == TRUE
][
  order(
    P_Ob,
    -Ob_vs_nonOb_effect_ratio
  )
]

fwrite(
  strict_supported,
  file.path(
    OUT_DIR,
    "10_original_obese_SVs_STRICT_crossBMI_support.tsv"
  ),
  sep = "\t"
)

# D. Relaxed supported rows
relaxed_supported <- master[
  Relaxed_obesity_enriched_support == TRUE
][
  order(
    P_Ob,
    -Ob_vs_nonOb_effect_ratio
  )
]

fwrite(
  relaxed_supported,
  file.path(
    OUT_DIR,
    "11_original_obese_SVs_RELAXED_crossBMI_support.tsv"
  ),
  sep = "\t"
)

# E. Strict supported + obese effect larger than both non-obese groups
strict_effect_enriched <- master[
  Strict_obesity_enriched_support == TRUE &
  Ob_stronger_than_sameDir_nonOb == TRUE
][
  order(
    P_Ob,
    -Ob_vs_nonOb_effect_ratio
  )
]

fwrite(
  strict_effect_enriched,
  file.path(
    OUT_DIR,
    "12_STRICT_support_Ob_effect_stronger_than_NwOw.tsv"
  ),
  sep = "\t"
)

# F. Relaxed supported + obese effect larger than both non-obese groups
relaxed_effect_enriched <- master[
  Relaxed_obesity_enriched_support == TRUE &
  Ob_stronger_than_sameDir_nonOb == TRUE
][
  order(
    P_Ob,
    -Ob_vs_nonOb_effect_ratio
  )
]

fwrite(
  relaxed_effect_enriched,
  file.path(
    OUT_DIR,
    "13_RELAXED_support_Ob_effect_stronger_than_NwOw.tsv"
  ),
  sep = "\t"
)

# G. Non-comparable original candidates
noncomparable <- candidate_summary[
  CrossBMI_comparable == FALSE
]

fwrite(
  noncomparable,
  file.path(
    OUT_DIR,
    "14_original_obese_SVs_not_comparable_crossBMI.tsv"
  ),
  sep = "\t"
)

# ==============================================================================
# 16. Overall summary
# ==============================================================================

n_comparable <- candidate_summary[
  CrossBMI_comparable == TRUE,
  uniqueN(Candidate_ID)
]

summary2 <- data.table(
  Category = c(
    "Original obese candidates",
    "Comparable original candidates in all-BMI common space",
    "Original candidates not comparable",
    "Cross-BMI strict pattern: Ob P<0.01; Nw/Ow P>=0.05",
    "Strict pattern + discovery/crossBMI obese direction concordant",
    "Cross-BMI relaxed pattern: Ob P<0.05; Nw/Ow P>=0.05",
    "Relaxed pattern + discovery/crossBMI obese direction concordant",
    "Strict support + obese effect stronger than any same-direction non-Ob effect",
    "Relaxed support + |Beta_Ob| > both |Beta_Nw| and |Beta_Ow|"
  ),

  N = c(
    uniqueN(discovery$Candidate_ID),

    n_comparable,

    candidate_summary[
      CrossBMI_comparable == FALSE,
      uniqueN(Candidate_ID)
    ],

    candidate_summary[
      Any_strict_pattern == TRUE,
      uniqueN(Candidate_ID)
    ],

    candidate_summary[
      Strict_obesity_enriched_support == TRUE,
      uniqueN(Candidate_ID)
    ],

    candidate_summary[
      Any_relaxed_pattern == TRUE,
      uniqueN(Candidate_ID)
    ],

    candidate_summary[
      Relaxed_obesity_enriched_support == TRUE,
      uniqueN(Candidate_ID)
    ],

    master[
      Strict_obesity_enriched_support == TRUE &
      Ob_stronger_than_sameDir_nonOb == TRUE,
      uniqueN(Candidate_ID)
    ],

    master[
      Relaxed_obesity_enriched_support == TRUE &
      Ob_stronger_than_sameDir_nonOb == TRUE,
      uniqueN(Candidate_ID)
    ]
  )
)

summary2[
  ,
  Percent_of_comparable :=
    fifelse(
      Category %in% c(
        "Original obese candidates",
        "Comparable original candidates in all-BMI common space",
        "Original candidates not comparable"
      ) |
      n_comparable == 0,
      NA_real_,
      100 * N / n_comparable
    )
]

fwrite(
  summary2,
  file.path(
    OUT_DIR,
    "15_original205_crossBMI_support_summary.tsv"
  ),
  sep = "\t"
)

# ==============================================================================
# 17. SV-type summary
# ==============================================================================

sv_type_col <- if (
  "SV_type" %in% names(candidate_summary)
) {
  "SV_type"
} else if (
  "SV_type_discovery_table" %in% names(candidate_summary)
) {
  "SV_type_discovery_table"
} else {
  NA_character_
}

if (!is.na(sv_type_col)) {
  sv_summary <- candidate_summary[
    ,
    .(
      N_original = .N,
      N_comparable = sum(
        CrossBMI_comparable,
        na.rm = TRUE
      ),
      N_strict_support = sum(
        Strict_obesity_enriched_support,
        na.rm = TRUE
      ),
      N_relaxed_support = sum(
        Relaxed_obesity_enriched_support,
        na.rm = TRUE
      )
    ),
    by = sv_type_col
  ]

  setnames(
    sv_summary,
    sv_type_col,
    "SV_type"
  )

  fwrite(
    sv_summary,
    file.path(
      OUT_DIR,
      "16_original205_crossBMI_support_by_SV_type.tsv"
    ),
    sep = "\t"
  )
} else {
  sv_summary <- NULL
}

# ==============================================================================
# 18. Analysis notes
# ==============================================================================

analysis_notes <- c(
  "STRICT cross-BMI pattern: Ob P<0.01 and Nw/Ow P>=0.05.",
  "RELAXED cross-BMI pattern: Ob P<0.05 and Nw/Ow P>=0.05.",
  "Supported pattern additionally requires concordant obese Beta direction between original obese-only discovery and all-BMI common-space analysis.",
  "Effect-size strengthening (|Beta_Ob| > both |Beta_Nw| and |Beta_Ow|) is retained as an annotation and secondary prioritization criterion.",
  "Candidates not comparable in the jointly defined all-BMI SV feature space are not interpreted as non-specific; they are classified as not evaluable.",
  "Cross-BMI significance patterns are descriptive evidence of obesity-enriched CRC association and are not formal Disease x BMI interaction tests.",
  paste0(
    "Spearman correlation between original obese discovery Beta and cross-BMI obese Beta = ",
    ifelse(
      is.na(effect_cor),
      "NA",
      signif(effect_cor, 4)
    )
  )
)

writeLines(
  analysis_notes,
  con = file.path(
    OUT_DIR,
    "17_crossBMI_specificity_analysis_notes.txt"
  )
)

# ==============================================================================
# 19. Console report
# ==============================================================================

cat(
  "\n\n========================================\n",
  "ORIGINAL OBESE SV vs CROSS-BMI SUMMARY\n",
  "========================================\n",
  sep = ""
)

print(
  summary2
)

if (!is.null(sv_summary)) {
  cat(
    "\nBy SV type:\n"
  )

  print(
    sv_summary
  )
}

cat(
  "\nDiscovery vs cross-BMI obese Beta Spearman rho:",
  effect_cor,
  "\n"
)

cat(
  "\nTop STRICT supported candidates:\n"
)

show_cols <- intersect(
  c(
    "Candidate_ID",
    "Joint_SV_ID",
    "SGB_crossBMI",
    "SV_type_crossBMI",
    "Beta_discovery",
    "P_discovery",
    "Q_type_discovery",
    "Beta_Nw",
    "P_Nw",
    "Beta_Ow",
    "P_Ow",
    "Beta_Ob",
    "P_Ob",
    "Ob_stronger_than_sameDir_nonOb",
    "Ob_vs_nonOb_effect_ratio",
    "Match_types",
    "Mapping_topologies"
  ),
  names(strict_supported)
)

print(
  strict_supported[
    1:min(20, .N),
    ..show_cols
  ]
)

cat(
  "\nOutput directory:\n",
  OUT_DIR,
  "\n"
)

cat(
  "\nDone.\n"
)
