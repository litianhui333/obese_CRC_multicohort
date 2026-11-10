# ==============================================================================
# 02_prepare_BMI_SV_circos_inputs.R
#
# Purpose
# -------
# Prepare six Circos/chord-diagram inputs for:
#
#   Normal-weight: dSV / vSV
#   Overweight:    dSV / vSV
#   Obese:         dSV / vSV
#
# Link weight = number of UNIQUE SVs connecting one SGB to one functional class.
#
# The script prepares four analysis subsets so the annotation does not need to
# be repeated later:
#
#   all_BMI_enriched
#   heterogeneity_supported
#   FDR_supported
#   FDR_and_heterogeneity
#
# Output
# ------
# 01_SGB_sector_labels.tsv
# 02_circos_panel_QC.tsv
# and six link tables under each subset directory.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ==============================================================================
# 1. Paths
# ==============================================================================

PROJECT <- "H:/My Drive/proj_ob_AA_CRC"

ANALYSIS_DIR <- file.path(
  PROJECT,
  "analyses",
  "2026-09-29_BMI_enriched_SV_function_modelQC_speciesResolved"
)

FUNCTION_DIR <- file.path(
  ANALYSIS_DIR,
  "01_function_classification"
)

OUT_ROOT <- file.path(
  ANALYSIS_DIR,
  "02_circos_inputs"
)

dir.create(
  OUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)

SV_FUNCTION_FILE <- file.path(
  FUNCTION_DIR,
  "03_SV_function_categories_long.tsv"
)

SGB_MAP_FILE <- file.path(
  PROJECT,
  "02.processed_data",
  "05_SV_result",
  "SGV_reference_to_SGB_final.tsv"
)

stopifnot(
  file.exists(
    SV_FUNCTION_FILE
  )
)

# ==============================================================================
# 2. Fixed order
# ==============================================================================

FUNCTION_ORDER <- c(
  "Carbohydrate metabolism",
  "Other metabolic processes",
  "Cell envelope, transport and ion homeostasis",
  "Defense, stress and protein quality control",
  "Gene regulation and signaling",
  "DNA replication, repair and chromosome dynamics",
  "Mobile genetic elements / genome plasticity",
  "Others"
)

BMI_LEVELS <- c(
  "Nw",
  "Ow",
  "Ob"
)

BMI_FILE_LABEL <- c(
  Nw = "NormalWeight",
  Ow = "Overweight",
  Ob = "Obese"
)

SV_TYPES <- c(
  "dSGV",
  "vSGV"
)

SV_FILE_LABEL <- c(
  dSGV = "dSV",
  vSGV = "vSV"
)

# ==============================================================================
# 3. Helpers
# ==============================================================================

as_bool <- function(x) {

  if (is.logical(x)) {
    return(x)
  }

  tolower(
    trimws(
      as.character(x)
    )
  ) %in% c(
    "true",
    "t",
    "1",
    "yes",
    "y"
  )
}

extract_tax_rank <- function(
    taxonomy,
    prefix
) {

  vapply(
    taxonomy,
    function(x) {

      if (
        is.na(x) ||
          x == ""
      ) {
        return(NA_character_)
      }

      pieces <- unlist(
        strsplit(
          x,
          "\\|"
        )
      )

      hit <- pieces[
        startsWith(
          pieces,
          prefix
        )
      ]

      if (length(hit) == 0) {
        return(NA_character_)
      }

      sub(
        paste0(
          "^",
          prefix
        ),
        "",
        hit[1]
      )
    },
    character(1)
  )
}

pretty_taxon <- function(x) {

  x <- gsub(
    "_",
    " ",
    x,
    fixed = TRUE
  )

  x[
    x %in% c(
      "",
      "NA",
      "unclassified",
      "unknown",
      "sp"
    )
  ] <- NA_character_

  x
}

species_with_sgb <- function(
    species,
    sgb
) {

  species <- as.character(species)
  sgb <- as.character(sgb)

  sgb_clean <- sub(
    "^t__",
    "",
    sgb
  )

  out <- ifelse(
    !is.na(species) &
      nzchar(species) &
      !is.na(sgb_clean) &
      nzchar(sgb_clean),
    paste0(
      species,
      " (",
      sgb_clean,
      ")"
    ),
    NA_character_
  )

  out
}

# ==============================================================================
# 4. Read SV-function table
# ==============================================================================

x <- fread(
  SV_FUNCTION_FILE,
  check.names = FALSE
)

required <- c(
  "BMI_SV_ID",
  "BMI_enriched_group",
  "SV_type",
  "SGB",
  "Function_category"
)

missing <- setdiff(
  required,
  names(x)
)

if (length(missing) > 0) {

  stop(
    "Missing columns: ",
    paste(
      missing,
      collapse = ", "
    )
  )
}

if (!"Heterogeneity_supported" %in% names(x)) {
  x[
    ,
    Heterogeneity_supported :=
      FALSE
  ]
}

if (!"FDR_supported" %in% names(x)) {
  x[
    ,
    FDR_supported :=
      FALSE
  ]
}

x[
  ,
  Heterogeneity_supported :=
    as_bool(
      Heterogeneity_supported
    )
]

x[
  ,
  FDR_supported :=
    as_bool(
      FDR_supported
    )
]

x <- unique(
  x[
    BMI_enriched_group %in%
      BMI_LEVELS &
      SV_type %in%
      SV_TYPES &
      Function_category %in%
      FUNCTION_ORDER
  ],
  by = c(
    "BMI_SV_ID",
    "SGB",
    "BMI_enriched_group",
    "SV_type",
    "Function_category"
  )
)

# ==============================================================================
# 5. Species-resolved SGB taxonomy labels from Original_Feature only
# ==============================================================================

sgb_ids <- sort(
  unique(
    x$SGB
  )
)

if (!file.exists(SGB_MAP_FILE)) {
  stop("Cannot find SGB mapping file: ", SGB_MAP_FILE)
}

mp <- fread(
  SGB_MAP_FILE,
  check.names = FALSE
)

if (!all(c("SGB", "Original_Feature") %in% names(mp))) {
  stop("SGB_MAP_FILE must contain SGB and Original_Feature.")
}

mp[
  ,
  Original_Feature :=
    trimws(
      as.character(
        Original_Feature
      )
    )
]

mp <- unique(
  mp[
    SGB %in% sgb_ids &
      !is.na(Original_Feature) &
      Original_Feature != "" &
      grepl("(^|\\\\|)s__", Original_Feature),
    .(
      SGB,
      Original_Feature
    )
  ],
  by = "SGB"
)

# Hard safety filter: no Excel_taxonomy / genus / raw-SGB fallback.
x <- x[
  SGB %in% mp$SGB
]

sgb_ids <- sort(
  unique(
    x$SGB
  )
)

sgb_labels <- data.table(
  SGB = sgb_ids
)

sgb_labels[
  mp,
  on = "SGB",
  Original_Feature :=
    i.Original_Feature
]

sgb_labels[
  ,
  Species :=
    pretty_taxon(
      extract_tax_rank(
        Original_Feature,
        "s__"
      )
    )
]

# After the hard Original_Feature filter, missing species labels should not occur.
sgb_labels <- sgb_labels[
  !is.na(Species) &
    Species != ""
]

sgb_labels[
  ,
  SGB_label :=
    species_with_sgb(
      Species,
      SGB
    )
]

fwrite(
  sgb_labels,
  file.path(
    OUT_ROOT,
    "01_SGB_sector_labels.tsv"
  ),
  sep = "\t"
)

x <- merge(
  x,
  sgb_labels[
    ,
    .(
      SGB,
      SGB_label,
      Species,
      Original_Feature
    )
  ],
  by = "SGB",
  all = FALSE,
  sort = FALSE
)

cat(
  "Species-resolved SGBs entering Circos:",
  uniqueN(x$SGB),
  "\n"
)

cat(
  "Species-resolved SVs entering Circos:",
  uniqueN(x$BMI_SV_ID),
  "\n"
)

# ==============================================================================
# 6. Define subsets
# ==============================================================================

subset_definitions <- list(
  all_BMI_enriched =
    rep(
      TRUE,
      nrow(x)
    ),

  heterogeneity_supported =
    x$Heterogeneity_supported,

  FDR_supported =
    x$FDR_supported,

  FDR_and_heterogeneity =
    x$FDR_supported &
    x$Heterogeneity_supported
)

qc_parts <- list()

# ==============================================================================
# 7. Generate six inputs for every subset
# ==============================================================================

for (
  subset_name in names(
    subset_definitions
  )
) {

  subset_dir <- file.path(
    OUT_ROOT,
    subset_name
  )

  dir.create(
    subset_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dat <- x[
    subset_definitions[[subset_name]]
  ]

  combined_links <- list()

  for (
    bmi in BMI_LEVELS
  ) {

    for (
      typ in SV_TYPES
    ) {

      panel <- dat[
        BMI_enriched_group == bmi &
          SV_type == typ
      ]

      links <- panel[
        ,
        .(
          N_SV =
            uniqueN(
              BMI_SV_ID
            )
        ),
        by = .(
          SGB,
          SGB_label,
          Function_category
        )
      ]

      links[
        ,
        Function_category :=
          factor(
            Function_category,
            levels = FUNCTION_ORDER
          )
      ]

      setorder(
        links,
        SGB,
        Function_category
      )

      links[
        ,
        Function_category :=
          as.character(
            Function_category
        )
      ]

      filename <- paste0(
        "circos_links_",
        BMI_FILE_LABEL[[bmi]],
        "_",
        SV_FILE_LABEL[[typ]],
        ".tsv"
      )

      fwrite(
        links,
        file.path(
          subset_dir,
          filename
        ),
        sep = "\t"
      )

      if (nrow(links) > 0) {

        combined_links[[
            paste(
              bmi,
              typ,
              sep = "_"
            )
          ]] <- copy(
          links
        )[
          ,
          `:=`(
            BMI_enriched_group =
              bmi,
            SV_type =
              typ
          )
        ]
      }

      qc_parts[[
          paste(
            subset_name,
            bmi,
            typ,
            sep = "_"
          )
        ]] <- data.table(
        Subset = subset_name,
        BMI_enriched_group = bmi,
        SV_type = typ,
        N_annotated_SV =
          uniqueN(
            panel$BMI_SV_ID
          ),
        N_SGB =
          uniqueN(
            panel$SGB
          ),
        N_function_categories =
          uniqueN(
            panel$Function_category
          ),
        N_SGB_function_links =
          nrow(
            links
          )
      )
    }
  }

  combined <- rbindlist(
    combined_links,
    fill = TRUE
  )

  fwrite(
    combined,
    file.path(
      subset_dir,
      "circos_links_ALL_6_panels.tsv"
    ),
    sep = "\t"
  )
}

qc <- rbindlist(
  qc_parts,
  fill = TRUE
)

fwrite(
  qc,
  file.path(
    OUT_ROOT,
    "02_circos_panel_QC.tsv"
  ),
  sep = "\t"
)

# ==============================================================================
# 8. Console report
# ==============================================================================

cat(
  "\n========================================\n"
)

cat(
  "CIRCOS INPUTS COMPLETE\n"
)

cat(
  "========================================\n"
)

print(
  qc[
    Subset == "all_BMI_enriched"
  ]
)

cat(
  "\nOutput root:\n",
  OUT_ROOT,
  "\n",
  sep = ""
)
