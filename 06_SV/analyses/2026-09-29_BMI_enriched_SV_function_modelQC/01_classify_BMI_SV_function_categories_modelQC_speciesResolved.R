# ==============================================================================
# 01_classify_BMI_SV_function_categories.R
#
# Purpose
# -------
# Classify eggNOG-annotated CDS/proteins from BMI-enriched CRC-associated SVs
# into the same fixed 8 functional classes used for the previous STRICT73
# analysis, then collapse protein/CDS-level assignments back to the SV level.
#
# Input
# -----
# 18_SV_CDS_eggNOG_annotation_long.tsv
# 19_ALL337_SV_function_annotation_master.tsv
#
# Output
# ------
# 01_function_category_dictionary.tsv
# 02_CDS_function_categories_long.tsv
# 03_SV_function_categories_long.tsv
# 04_SV_function_category_summary.tsv
# 05_function_classification_QC.tsv
# 06_SV_function_manual_review.tsv
#
# Notes
# -----
# - Functional unit for downstream analysis = SV, not CDS.
# - One SV can belong to multiple functional categories.
# - Only proteins with an eggNOG annotation are functionally classified.
# - The seven SVs with no CDS overlap and the eight CDS-overlapping SVs without
#   any eggNOG annotation do not enter the functional-category analysis.
# - Category rules reproduce the fixed 8-class framework used in the
#   preliminary STRICT73 curation; assignments are rule-based and auditable.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ==============================================================================
# 1. Paths
# ==============================================================================

INPUT_DIR <- "H:/My Drive/proj_ob_AA_CRC/02.processed_data/05_SV_result/SV_annotation/BMI_enriched_SV"

SGB_MAP_FILE <- "H:/My Drive/proj_ob_AA_CRC/02.processed_data/05_SV_result/SGV_reference_to_SGB_final.tsv"

ANALYSIS_DIR <- "H:/My Drive/proj_ob_AA_CRC/analyses/2026-09-29_BMI_enriched_SV_function_modelQC_speciesResolved"

OUT_DIR <- file.path(
  ANALYSIS_DIR,
  "01_function_classification"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

LONG_FILE <- file.path(
  INPUT_DIR,
  "18_SV_CDS_eggNOG_annotation_long.tsv"
)

QC_DIR <- file.path(
  "H:/My Drive/proj_ob_AA_CRC/analyses",
  "2026-09-28_groupwise_BMI_enriched_SV",
  "model_QC_filtered"
)

KEEP_FILE <- file.path(
  QC_DIR,
  "07_final_BMI_enriched_keep_IDs.tsv"
)

SV_MASTER_FILE <- file.path(
  QC_DIR,
  "08_filtered_function_annotation_master.tsv"
)

stopifnot(
  file.exists(LONG_FILE),
  file.exists(SV_MASTER_FILE),
  file.exists(KEEP_FILE),
  file.exists(SGB_MAP_FILE)
)

# ==============================================================================
# 2. Fixed functional framework
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

category_dictionary <- data.table(
  Detailed_category = c(
    "Carbohydrate metabolism",
    "Amino acid metabolism",
    "Energy metabolism",
    "Lipid metabolism",
    "Cofactor/vitamin metabolism",
    "Nucleotide metabolism",
    "Other metabolic processes",
    "Ion/metal/phosphate homeostasis",
    "Cell envelope/transport/secretion",
    "Defense/stress response",
    "Gene regulation/signaling",
    "DNA replication/repair/chromosome dynamics",
    "Mobile genetic elements/genome plasticity",
    "Protein synthesis/RNA processing",
    "Unknown"
  ),
  Function_category = c(
    "Carbohydrate metabolism",
    rep("Other metabolic processes", 6),
    rep("Cell envelope, transport and ion homeostasis", 2),
    "Defense, stress and protein quality control",
    "Gene regulation and signaling",
    "DNA replication, repair and chromosome dynamics",
    "Mobile genetic elements / genome plasticity",
    "Others",
    "Others"
  ),
  Working_definition = c(
    "Sugar uptake/utilization, glycoside hydrolases, PTS and central carbohydrate-related functions.",
    "Amino-acid biosynthesis, degradation and related metabolism.",
    "Central energy conversion, respiration and energy-conserving processes.",
    "Membrane lipid and fatty-acid related enzymes.",
    "Vitamin/cofactor biosynthesis and utilization.",
    "Purine/pyrimidine synthesis and related nucleotide metabolism.",
    "Other metabolism, including secondary-metabolite related functions not captured above.",
    "Metal, phosphate and inorganic-ion acquisition/homeostasis.",
    "Surface structures, membrane proteins, transporters, sortases and secretion machinery.",
    "Stress response, chaperones, antimicrobial/host defense and protein quality-control functions.",
    "Transcription factors, sigma factors, histidine kinases and related regulatory/signaling systems.",
    "Replication, DNA repair, recombination, topoisomerases, chromosome partition and cell-division functions.",
    "Integrases, transposases, excisionases, relaxases, reverse transcriptases and phage/mobile-element functions.",
    "Translation, ribosome and RNA-processing/modification functions.",
    "No informative assignment beyond poorly characterized/general-function annotation."
  )
)

fwrite(
  category_dictionary,
  file.path(
    OUT_DIR,
    "01_function_category_dictionary.tsv"
  ),
  sep = "\t"
)

# ==============================================================================
# 3. Helpers
# ==============================================================================

as_bool <- function(x) {

  if (is.logical(x)) {
    return(x)
  }

  y <- tolower(trimws(as.character(x)))

  y %in% c(
    "true",
    "t",
    "1",
    "yes",
    "y"
  )
}

clean_field <- function(x) {

  x <- as.character(x)

  x[
    is.na(x) |
      x %in% c(
        "",
        "-",
        "NA",
        "NaN",
        "nan"
      )
  ] <- ""

  x
}

has_keyword <- function(text, pattern) {

  grepl(
    pattern,
    text,
    ignore.case = TRUE,
    perl = TRUE
  )
}

has_cog <- function(cog, letters) {

  cog <- toupper(
    clean_field(cog)
  )

  grepl(
    paste0(
      "[",
      letters,
      "]"
    ),
    cog,
    perl = TRUE
  )
}

nonempty_annotation <- function(x) {

  clean_field(x) != ""
}

# ==============================================================================
# 4. Read data
# ==============================================================================

cds <- fread(
  LONG_FILE,
  check.names = FALSE
)

sv_master <- fread(
  SV_MASTER_FILE,
  check.names = FALSE
)

keep_ids <- fread(
  KEEP_FILE,
  check.names = FALSE
)

if (!"BMI_SV_ID" %in% names(keep_ids)) {
  stop("KEEP_FILE must contain BMI_SV_ID.")
}

KEEP_IDS <- unique(
  as.character(
    keep_ids$BMI_SV_ID
  )
)

# First restrict to model-QC-retained candidates.
cds <- cds[
  BMI_SV_ID %in% KEEP_IDS
]

sv_master <- sv_master[
  BMI_SV_ID %in% KEEP_IDS
]

cat(
  "Model-QC retained SVs before species-resolution filter:",
  uniqueN(sv_master$BMI_SV_ID),
  "\n"
)

# --------------------------------------------------------------------------
# Species-resolution filter
# --------------------------------------------------------------------------
# Keep only SGBs that map back to the MetaPhlAn4 species-level feature used in
# the manuscript (Original_Feature). We intentionally do NOT fall back to
# Excel_taxonomy, genus-only labels, or raw SGB IDs for functional figures.
#
# This is a figure/functional-interpretation filter, not a deletion from the
# statistical SV association catalog.

sgb_map <- fread(
  SGB_MAP_FILE,
  check.names = FALSE
)

if (!all(c("SGB", "Original_Feature") %in% names(sgb_map))) {
  stop("SGB_MAP_FILE must contain SGB and Original_Feature.")
}

sgb_map[
  ,
  Original_Feature :=
    trimws(
      as.character(
        Original_Feature
      )
    )
]

valid_species_map <- unique(
  sgb_map[
    !is.na(SGB) &
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

valid_sgbs <- valid_species_map$SGB

before_species_filter <- copy(sv_master)

excluded_no_original <- before_species_filter[
  !SGB %in% valid_sgbs
][
  order(
    BMI_enriched_group,
    SV_type,
    SGB,
    BMI_SV_ID
  )
]

species_filter_qc <- rbindlist(
  list(
    before_species_filter[
      ,
      .(
        Stage = "Before_species_filter",
        N_SV = uniqueN(BMI_SV_ID),
        N_SGB = uniqueN(SGB)
      ),
      by = .(
        BMI_enriched_group,
        SV_type
      )
    ],
    before_species_filter[
      SGB %in% valid_sgbs,
      .(
        Stage = "After_species_filter",
        N_SV = uniqueN(BMI_SV_ID),
        N_SGB = uniqueN(SGB)
      ),
      by = .(
        BMI_enriched_group,
        SV_type
      )
    ]
  ),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  excluded_no_original,
  file.path(
    OUT_DIR,
    "00_excluded_SVs_without_Original_Feature.tsv"
  ),
  sep = "\t"
)

fwrite(
  species_filter_qc,
  file.path(
    OUT_DIR,
    "00_species_resolution_filter_QC.tsv"
  ),
  sep = "\t"
)

sv_master <- sv_master[
  SGB %in% valid_sgbs
]

SPECIES_RESOLVED_IDS <- unique(
  as.character(
    sv_master$BMI_SV_ID
  )
)

cds <- cds[
  BMI_SV_ID %in% SPECIES_RESOLVED_IDS
]

cat(
  "Species-resolved SVs retained for function pipeline:",
  uniqueN(sv_master$BMI_SV_ID),
  "\n"
)

cat(
  "SVs excluded because SGB lacks Original_Feature:",
  uniqueN(excluded_no_original$BMI_SV_ID),
  "\n"
)

cat(
  "\nSpecies-resolution filter QC:\n"
)

print(
  species_filter_qc[
    order(
      BMI_enriched_group,
      SV_type,
      Stage
    )
  ]
)

required_cds <- c(
  "BMI_SV_ID",
  "Protein_key",
  "EggNOG_query",
  "Has_eggNOG",
  "COG_category",
  "Description",
  "Preferred_name",
  "KEGG_ko",
  "KEGG_Pathway",
  "CAZy",
  "PFAMs"
)

missing_cds <- setdiff(
  required_cds,
  names(cds)
)

if (length(missing_cds) > 0) {

  stop(
    "Missing columns in eggNOG long table: ",
    paste(
      missing_cds,
      collapse = ", "
    )
  )
}

required_sv <- c(
  "BMI_SV_ID",
  "BMI_enriched_group",
  "SGB",
  "SV_type"
)

missing_sv <- setdiff(
  required_sv,
  names(sv_master)
)

if (length(missing_sv) > 0) {

  stop(
    "Missing columns in SV master table: ",
    paste(
      missing_sv,
      collapse = ", "
    )
  )
}

cds[
  ,
  Has_eggNOG :=
    as_bool(
      Has_eggNOG
    )
]

for (col in c(
  "Description",
  "Preferred_name",
  "COG_category",
  "KEGG_ko",
  "KEGG_Pathway",
  "KEGG_Module",
  "CAZy",
  "PFAMs",
  "eggNOG_OGs"
)) {

  if (!col %in% names(cds)) {
    cds[
      ,
      (col) := ""
    ]
  }

  set(
    cds,
    j = col,
    value = clean_field(
      cds[[col]]
    )
  )
}

# Canonical SV metadata comes from the 337-SV master.
sv_meta_cols <- intersect(
  c(
    "BMI_SV_ID",
    "BMI_enriched_group",
    "SGB",
    "SV_type",
    "SV_feature",
    "Reference",
    "SV_interval",
    "Target_P",
    "Target_Q",
    "FDR_supported",
    "P_heterogeneity",
    "I2_heterogeneity",
    "Heterogeneity_supported"
  ),
  names(sv_master)
)

sv_meta <- unique(
  sv_master[
    ,
    ..sv_meta_cols
  ],
  by = "BMI_SV_ID"
)

if ("FDR_supported_modelQC" %in% names(sv_master)) {

  fdr_qc <- unique(
    sv_master[
      ,
      .(
        BMI_SV_ID,
        FDR_supported_modelQC
      )
    ],
    by = "BMI_SV_ID"
  )

  sv_meta <- merge(
    sv_meta,
    fdr_qc,
    by = "BMI_SV_ID",
    all.x = TRUE,
    sort = FALSE
  )

  sv_meta[
    ,
    FDR_supported :=
      as_bool(
        FDR_supported_modelQC
      )
  ]

} else if ("FDR_supported" %in% names(sv_meta)) {

  sv_meta[
    ,
    FDR_supported :=
      as_bool(
        FDR_supported
      )
  ]
}

if ("Heterogeneity_supported" %in% names(sv_meta)) {
  sv_meta[
    ,
    Heterogeneity_supported :=
      as_bool(
        Heterogeneity_supported
      )
  ]
}

# ==============================================================================
# 5. Build annotation text used by keyword rules
# ==============================================================================

cds[
  ,
  Annotation_text :=
    tolower(
      paste(
        Description,
        Preferred_name,
        PFAMs,
        CAZy,
        sep = " | "
      )
    )
]

# ==============================================================================
# 6. Detailed functional rules
#
# Rules combine:
#   1) COG functional classes
#   2) product/domain keywords
#   3) CAZy presence for carbohydrate-active proteins
#
# Multiple assignments are allowed.
# ==============================================================================

# ---- Carbohydrate metabolism --------------------------------------------------

flag_carb <- (
  cds$Has_eggNOG &
    (
      has_cog(
        cds$COG_category,
        "G"
      ) |
        nonempty_annotation(
          cds$CAZy
        ) |
        has_keyword(
          cds$Annotation_text,
          paste0(
            "carbohyd|",
            "glycosid|glycosyl|glucosid|galactosid|mannosid|xylosid|arabinosid|",
            "cellulas|cellulose|starch|amylas|glucan|fructos|sucros|",
            "phosphotransferase system|\\bpts\\b|sugar transport"
          )
        )
    )
)

# ---- Amino-acid metabolism ----------------------------------------------------

flag_amino <- (
  cds$Has_eggNOG &
    (
      has_cog(
        cds$COG_category,
        "E"
      ) |
        has_keyword(
          cds$Annotation_text,
          paste0(
            "amino acid|aminotransferase|transaminase|",
            "glutamate|glutamine|aspartate|asparagine|",
            "lysine|arginine|histidine|proline|leucine|isoleucine|valine|",
            "tryptophan|tyrosine|phenylalanine|methionine|cysteine|serine|threonine"
          )
        )
    )
)

# ---- Energy metabolism --------------------------------------------------------

flag_energy <- (
  cds$Has_eggNOG &
    (
      has_cog(
        cds$COG_category,
        "C"
      ) |
        has_keyword(
          cds$Annotation_text,
          paste0(
            "energy metabolism|respirat|electron transport|",
            "atp synthase|cytochrome|quinone|ferredoxin|",
            "hydrogenase|fumarate reductase|nitrate reductase"
          )
        )
    )
)

# ---- Lipid metabolism ---------------------------------------------------------

flag_lipid <- (
  cds$Has_eggNOG &
    (
      has_cog(
        cds$COG_category,
        "I"
      ) |
        has_keyword(
          cds$Annotation_text,
          paste0(
            "fatty acid|lipid|phospholipid|acyl carrier|",
            "acyl-coa|beta-oxidation|lipase"
          )
        )
    )
)

# ---- Cofactor / vitamin metabolism -------------------------------------------

flag_cofactor <- (
  cds$Has_eggNOG &
    (
      has_cog(
        cds$COG_category,
        "H"
      ) |
        has_keyword(
          cds$Annotation_text,
          paste0(
            "cofactor|vitamin|thiamin|riboflavin|flavin|",
            "folate|cobalamin|biotin|pyridox|heme biosynth|molybdopterin"
          )
        )
    )
)

# ---- Nucleotide metabolism ----------------------------------------------------

flag_nucleotide <- (
  cds$Has_eggNOG &
    (
      has_cog(
        cds$COG_category,
        "F"
      ) |
        has_keyword(
          cds$Annotation_text,
          paste0(
            "purine|pyrimidine|nucleotide metabolism|nucleoside metabolism|",
            "adenylate|guanylate|cytidylate|uridylate|thymidylate"
          )
        )
    )
)

# ---- Other metabolic processes -----------------------------------------------

flag_other_metabolic <- (
  cds$Has_eggNOG &
    has_cog(
      cds$COG_category,
      "Q"
    )
)

# ---- Ion / metal / phosphate homeostasis -------------------------------------

flag_ion <- (
  cds$Has_eggNOG &
    (
      has_cog(
        cds$COG_category,
        "P"
      ) |
        has_keyword(
          cds$Annotation_text,
          paste0(
            "ion homeostasis|metal homeostasis|phosphate transport|",
            "iron transport|ferric|ferrous|siderophore|",
            "potassium transport|sodium transport|magnesium transport|",
            "manganese transport|zinc transport|copper transport|",
            "cation transport|anion transport"
          )
        )
    )
)

# ---- Cell envelope / transport / secretion -----------------------------------

flag_envelope <- (
  cds$Has_eggNOG &
    (
      has_cog(
        cds$COG_category,
        "MUWN"
      ) |
        has_keyword(
          cds$Annotation_text,
          paste0(
            "transporter|transport system|permease|efflux|",
            "membrane protein|outer membrane|inner membrane|porin|",
            "cell wall|cell envelope|peptidoglycan|lipopolysaccharide|",
            "secretion|secretory|sortase|surface protein|",
            "abc[- ]type|abc transporter"
          )
        )
    )
)

# ---- Defense / stress / protein quality control ------------------------------

flag_defense <- (
  cds$Has_eggNOG &
    (
      has_cog(
        cds$COG_category,
        "VO"
      ) |
        has_keyword(
          cds$Annotation_text,
          paste0(
            "defense|resistance|antimicrobial|antibiotic|",
            "toxin[- ]antitoxin|restriction[- ]modification|crispr|",
            "stress response|heat shock|chaperon|protein quality|",
            "catalase|peroxidase|thioredoxin|glutaredoxin|",
            "oxidative stress|proteostasis"
          )
        )
    )
)

# ---- Gene regulation / signaling ---------------------------------------------

flag_regulation <- (
  cds$Has_eggNOG &
    (
      has_cog(
        cds$COG_category,
        "KT"
      ) |
        has_keyword(
          cds$Annotation_text,
          paste0(
            "transcription factor|transcriptional regulator|",
            "\\bregulator\\b|sigma[- ]?factor|histidine kinase|",
            "response regulator|two[- ]component|sensor kinase|",
            "helix[-_ ]turn[-_ ]helix|\\bhth\\b|",
            "\\bpas\\b|\\beal\\b|\\bggdef\\b|cyclic[- ]di[- ]gmp|",
            "signal transduction"
          )
        )
    )
)

# ---- DNA replication / repair / chromosome dynamics --------------------------

flag_dna <- (
  cds$Has_eggNOG &
    (
      has_cog(
        cds$COG_category,
        "LDB"
      ) |
        has_keyword(
          cds$Annotation_text,
          paste0(
            "dna replication|dna repair|dna polymerase|dna helicase|",
            "topoisomerase|dna gyrase|mismatch repair|excision repair|",
            "homologous recombination|chromosome partition|",
            "\\bpara\\b|\\bparb\\b|primase|cell division|septation"
          )
        )
    )
)

# ---- Mobile genetic elements / genome plasticity -----------------------------

flag_mobile <- (
  cds$Has_eggNOG &
    (
      has_cog(
        cds$COG_category,
        "X"
      ) |
        has_keyword(
          cds$Annotation_text,
          paste0(
            "phage|prophage|bacteriophage|",
            "transposase|transposition|transposable|insertion sequence|",
            "integrase|excisionase|relaxase|conjugative|conjugation|",
            "reverse transcriptase|retrotranspos|mobile element|",
            "capsid|terminase|baseplate|tape measure|tail protein|",
            "portal protein|holin|lysogen"
          )
        )
    )
)

# ---- Protein synthesis / RNA processing --------------------------------------

flag_protein_rna <- (
  cds$Has_eggNOG &
    (
      has_cog(
        cds$COG_category,
        "JA"
      ) |
        has_keyword(
          cds$Annotation_text,
          paste0(
            "ribosomal|ribosome|translation|",
            "aminoacyl[- ]trna|trna synthetase|trna modification|",
            "rrna|rna processing|rna modification|ribosome biogenesis"
          )
        )
    )
)

# ==============================================================================
# 7. Build detailed category table
# ==============================================================================

category_parts <- list()

add_category <- function(flag, detailed_category) {

  tmp <- cds[
    flag == TRUE &
      Has_eggNOG == TRUE,
    .(
      BMI_SV_ID,
      Protein_key,
      EggNOG_query,
      COG_category,
      Description,
      Preferred_name,
      KEGG_ko,
      KEGG_Pathway,
      CAZy,
      PFAMs
    )
  ]

  if (nrow(tmp) == 0) {
    return(NULL)
  }

  tmp[
    ,
    Detailed_category :=
      detailed_category
  ]

  tmp
}

category_parts[["Carbohydrate metabolism"]] <-
  add_category(
    flag_carb,
    "Carbohydrate metabolism"
  )

category_parts[["Amino acid metabolism"]] <-
  add_category(
    flag_amino,
    "Amino acid metabolism"
  )

category_parts[["Energy metabolism"]] <-
  add_category(
    flag_energy,
    "Energy metabolism"
  )

category_parts[["Lipid metabolism"]] <-
  add_category(
    flag_lipid,
    "Lipid metabolism"
  )

category_parts[["Cofactor/vitamin metabolism"]] <-
  add_category(
    flag_cofactor,
    "Cofactor/vitamin metabolism"
  )

category_parts[["Nucleotide metabolism"]] <-
  add_category(
    flag_nucleotide,
    "Nucleotide metabolism"
  )

category_parts[["Other metabolic processes"]] <-
  add_category(
    flag_other_metabolic,
    "Other metabolic processes"
  )

category_parts[["Ion/metal/phosphate homeostasis"]] <-
  add_category(
    flag_ion,
    "Ion/metal/phosphate homeostasis"
  )

category_parts[["Cell envelope/transport/secretion"]] <-
  add_category(
    flag_envelope,
    "Cell envelope/transport/secretion"
  )

category_parts[["Defense/stress response"]] <-
  add_category(
    flag_defense,
    "Defense/stress response"
  )

category_parts[["Gene regulation/signaling"]] <-
  add_category(
    flag_regulation,
    "Gene regulation/signaling"
  )

category_parts[["DNA replication/repair/chromosome dynamics"]] <-
  add_category(
    flag_dna,
    "DNA replication/repair/chromosome dynamics"
  )

category_parts[["Mobile genetic elements/genome plasticity"]] <-
  add_category(
    flag_mobile,
    "Mobile genetic elements/genome plasticity"
  )

category_parts[["Protein synthesis/RNA processing"]] <-
  add_category(
    flag_protein_rna,
    "Protein synthesis/RNA processing"
  )

cds_categories <- rbindlist(
  category_parts,
  fill = TRUE
)

# ==============================================================================
# 8. Assign remaining eggNOG-annotated proteins to Unknown -> Others
# ==============================================================================

annotated_keys <- unique(
  cds[
    Has_eggNOG == TRUE,
    .(
      BMI_SV_ID,
      Protein_key
    )
  ]
)

classified_keys <- unique(
  cds_categories[
    ,
    .(
      BMI_SV_ID,
      Protein_key
    )
  ]
)

unknown_keys <- fsetdiff(
  annotated_keys,
  classified_keys
)

if (nrow(unknown_keys) > 0) {

  unknown_rows <- merge(
    cds[
      Has_eggNOG == TRUE,
      .(
        BMI_SV_ID,
        Protein_key,
        EggNOG_query,
        COG_category,
        Description,
        Preferred_name,
        KEGG_ko,
        KEGG_Pathway,
        CAZy,
        PFAMs
      )
    ],
    unknown_keys,
    by = c(
      "BMI_SV_ID",
      "Protein_key"
    )
  )

  unknown_rows[
    ,
    Detailed_category :=
      "Unknown"
  ]

  cds_categories <- rbindlist(
    list(
      cds_categories,
      unknown_rows
    ),
    fill = TRUE
  )
}

# ==============================================================================
# 9. Map detailed categories to fixed 8 classes
# ==============================================================================

cds_categories <- merge(
  cds_categories,
  category_dictionary[
    ,
    .(
      Detailed_category,
      Function_category
    )
  ],
  by = "Detailed_category",
  all.x = TRUE,
  sort = FALSE
)

if (any(is.na(cds_categories$Function_category))) {

  stop(
    "At least one detailed category was not mapped to an 8-class category."
  )
}

cds_categories <- unique(
  cds_categories,
  by = c(
    "BMI_SV_ID",
    "Protein_key",
    "Detailed_category",
    "Function_category"
  )
)

cds_categories <- merge(
  cds_categories,
  sv_meta,
  by = "BMI_SV_ID",
  all.x = TRUE,
  sort = FALSE
)

setcolorder(
  cds_categories,
  c(
    "BMI_SV_ID",
    "BMI_enriched_group",
    "SV_type",
    "SGB",
    "Protein_key",
    "EggNOG_query",
    "Detailed_category",
    "Function_category",
    setdiff(
      names(cds_categories),
      c(
        "BMI_SV_ID",
        "BMI_enriched_group",
        "SV_type",
        "SGB",
        "Protein_key",
        "EggNOG_query",
        "Detailed_category",
        "Function_category"
      )
    )
  )
)

fwrite(
  cds_categories,
  file.path(
    OUT_DIR,
    "02_CDS_function_categories_long.tsv"
  ),
  sep = "\t"
)

# ==============================================================================
# 10. Collapse to SV x functional-category level
# ==============================================================================

sv_function <- unique(
  cds_categories[
    ,
    .(
      BMI_SV_ID,
      BMI_enriched_group,
      SV_type,
      SGB,
      Function_category
    )
  ],
  by = c(
    "BMI_SV_ID",
    "Function_category"
  )
)

sv_function <- merge(
  sv_function,
  sv_meta,
  by = c(
    "BMI_SV_ID",
    "BMI_enriched_group",
    "SV_type",
    "SGB"
  ),
  all.x = TRUE,
  sort = FALSE
)

sv_function[
  ,
  Function_category :=
    factor(
      Function_category,
      levels = FUNCTION_ORDER
    )
]

setorder(
  sv_function,
  BMI_enriched_group,
  SV_type,
  SGB,
  BMI_SV_ID,
  Function_category
)

sv_function[
  ,
  Function_category :=
    as.character(
      Function_category
    )
]

fwrite(
  sv_function,
  file.path(
    OUT_DIR,
    "03_SV_function_categories_long.tsv"
  ),
  sep = "\t"
)

# ==============================================================================
# 11. Summary by BMI group, SV type and function
# ==============================================================================

sv_function_summary <- sv_function[
  ,
  .(
    N_SV =
      uniqueN(
        BMI_SV_ID
      ),

    N_SGB =
      uniqueN(
        SGB
      )
  ),
  by = .(
    BMI_enriched_group,
    SV_type,
    Function_category
  )
]

sv_function_summary[
  ,
  Function_category :=
    factor(
      Function_category,
      levels = FUNCTION_ORDER
    )
]

setorder(
  sv_function_summary,
  BMI_enriched_group,
  SV_type,
  Function_category
)

sv_function_summary[
  ,
  Function_category :=
    as.character(
      Function_category
    )
]

fwrite(
  sv_function_summary,
  file.path(
    OUT_DIR,
    "04_SV_function_category_summary.tsv"
  ),
  sep = "\t"
)

# ==============================================================================
# 12. QC
# ==============================================================================

annotated_sv_ids <- unique(
  cds[
    Has_eggNOG == TRUE,
    BMI_SV_ID
  ]
)

categorized_sv_ids <- unique(
  sv_function$BMI_SV_ID
)

qc <- merge(
  sv_meta[
    BMI_SV_ID %in% annotated_sv_ids,
    .(
      BMI_SV_ID,
      BMI_enriched_group,
      SV_type
    )
  ],
  data.table(
    BMI_SV_ID = categorized_sv_ids,
    Categorized = TRUE
  ),
  by = "BMI_SV_ID",
  all.x = TRUE
)

qc[
  is.na(Categorized),
  Categorized := FALSE
]

qc_summary <- qc[
  ,
  .(
    N_eggNOG_annotated_SV =
      uniqueN(
        BMI_SV_ID
      ),

    N_categorized_SV =
      uniqueN(
        BMI_SV_ID[
          Categorized
        ]
      ),

    N_not_categorized =
      uniqueN(
        BMI_SV_ID[
          !Categorized
        ]
      )
  ),
  by = .(
    BMI_enriched_group,
    SV_type
  )
]

fwrite(
  qc_summary,
  file.path(
    OUT_DIR,
    "05_function_classification_QC.tsv"
  ),
  sep = "\t"
)

# ==============================================================================
# 13. Manual-review table
#
# Flag:
# - SVs assigned only to Others
# - SVs assigned to >=5 broad categories
#
# These are useful for final manuscript curation, but they are not removed.
# ==============================================================================

sv_review <- sv_function[
  ,
  .(
    N_function_categories =
      uniqueN(
        Function_category
      ),

    Function_categories =
      paste(
        sort(
          unique(
            Function_category
          )
        ),
        collapse = "; "
      )
  ),
  by = .(
    BMI_SV_ID,
    BMI_enriched_group,
    SV_type,
    SGB
  )
]

sv_review[
  ,
  Manual_review_reason :=
    fifelse(
      Function_categories == "Others",
      "Only_Others",
      fifelse(
        N_function_categories >= 5,
        "Five_or_more_categories",
        NA_character_
      )
    )
]

sv_review <- sv_review[
  !is.na(
    Manual_review_reason
  )
]

fwrite(
  sv_review,
  file.path(
    OUT_DIR,
    "06_SV_function_manual_review.tsv"
  ),
  sep = "\t"
)

# ==============================================================================
# 14. Console report
# ==============================================================================

cat(
  "\n========================================\n"
)

cat(
  "FUNCTION CLASSIFICATION COMPLETE\n"
)

cat(
  "========================================\n"
)

cat(
  "\nQC:\n"
)

print(
  qc_summary
)

cat(
  "\nSV-level function summary:\n"
)

print(
  sv_function_summary
)

cat(
  "\nOutput directory:\n",
  OUT_DIR,
  "\n",
  sep = ""
)
