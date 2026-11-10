# ==============================================================================
# 03_plot_BMI_SV_function_circos.R
#
# Purpose
# -------
# Draw six Circos/chord diagrams:
#
#   Normal-weight dSV
#   Normal-weight vSV
#   Overweight dSV
#   Overweight vSV
#   Obese dSV
#   Obese vSV
#
# Left/first sectors: SGBs
# Right/second sectors: fixed 8 functional categories
# Link width: number of unique SVs
# Link color: functional category
#
# Output
# ------
# Six individual PDF + PNG figures and one combined 3 x 2 PDF.
#
# Change SUBSET_MODE below to plot a sensitivity subset without repeating
# annotation:
#
#   all_BMI_enriched
#   heterogeneity_supported
#   FDR_supported
#   FDR_and_heterogeneity
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(circlize)
})

# ==============================================================================
# 1. User settings
# ==============================================================================

PROJECT <- "H:/My Drive/proj_ob_AA_CRC"

ANALYSIS_DIR <- file.path(
  PROJECT,
  "analyses",
  "2026-09-29_BMI_enriched_SV_function_modelQC_speciesResolved"
)

# Primary/default plot:
SUBSET_MODE <- "all_BMI_enriched"

# Options:
# SUBSET_MODE <- "heterogeneity_supported"
# SUBSET_MODE <- "FDR_supported"
# SUBSET_MODE <- "FDR_and_heterogeneity"

INPUT_DIR <- file.path(
  ANALYSIS_DIR,
  "02_circos_inputs",
  SUBSET_MODE
)

LABEL_FILE <- file.path(
  ANALYSIS_DIR,
  "02_circos_inputs",
  "01_SGB_sector_labels.tsv"
)

OUT_DIR <- file.path(
  ANALYSIS_DIR,
  "03_circos_figures",
  SUBSET_MODE
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

MAKE_COMBINED_6PANEL <- TRUE

# ==============================================================================
# 2. Fixed function order and colors
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

FUNCTION_COLORS <- c(
  "Carbohydrate metabolism" =
    "#F4A261",

  "Other metabolic processes" =
    "#D6B85A",

  "Cell envelope, transport and ion homeostasis" =
    "#55B7B4",

  "Defense, stress and protein quality control" =
    "#72A96B",

  "Gene regulation and signaling" =
    "#568FD8",

  "DNA replication, repair and chromosome dynamics" =
    "#9566C4",

  "Mobile genetic elements / genome plasticity" =
    "#D76B73",

  "Others" =
    "#B8B8B8"
)

FUNCTION_DISPLAY <- c(
  "Carbohydrate metabolism" =
    "Carbohydrate\nmetabolism",

  "Other metabolic processes" =
    "Other metabolic\nprocesses",

  "Cell envelope, transport and ion homeostasis" =
    "Cell envelope, transport\nand ion homeostasis",

  "Defense, stress and protein quality control" =
    "Defense, stress and\nprotein quality control",

  "Gene regulation and signaling" =
    "Gene regulation\nand signaling",

  "DNA replication, repair and chromosome dynamics" =
    "DNA replication, repair\nand chromosome dynamics",

  "Mobile genetic elements / genome plasticity" =
    "Mobile genetic elements /\ngenome plasticity",

  "Others" =
    "Others"
)

# ==============================================================================
# 3. Panel definitions
# ==============================================================================

PANELS <- data.table(
  BMI = c(
    "Nw",
    "Nw",
    "Ow",
    "Ow",
    "Ob",
    "Ob"
  ),
  SV_type = c(
    "dSV",
    "vSV",
    "dSV",
    "vSV",
    "dSV",
    "vSV"
  ),
  BMI_file = c(
    "NormalWeight",
    "NormalWeight",
    "Overweight",
    "Overweight",
    "Obese",
    "Obese"
  ),
  Title = c(
    "Normal-weight — dSV",
    "Normal-weight — vSV",
    "Overweight — dSV",
    "Overweight — vSV",
    "Obese — dSV",
    "Obese — vSV"
  )
)

# ==============================================================================
# 4. Read SGB labels
# ==============================================================================

if (
  file.exists(
    LABEL_FILE
  )
) {

  label_dt <- fread(
    LABEL_FILE,
    check.names = FALSE
  )

  label_map <- setNames(
    label_dt$SGB_label,
    label_dt$SGB
  )

} else {

  label_map <- character(0)
}

# ==============================================================================
# 5. Plot helper
# ==============================================================================

plot_one_circos <- function(
    links,
    panel_title,
    small_panel = FALSE
) {

  if (nrow(links) == 0) {

    plot.new()

    text(
      0.5,
      0.5,
      paste0(
        panel_title,
        "\nNo links"
      ),
      cex = 1
    )

    return(
      invisible(NULL)
    )
  }

  links <- copy(
    links
  )

  links <- links[
    N_SV > 0
  ]

  # SGB order:
  # larger total SV contribution first, then alphabetical SGB ID.
  sgb_order <- links[
    ,
    .(
      Total_link_weight =
        sum(
          N_SV
        )
    ),
    by = SGB
  ][
    order(
      -Total_link_weight,
      SGB
    ),
    SGB
  ]

  function_present <- FUNCTION_ORDER[
    FUNCTION_ORDER %in%
      links$Function_category
  ]

  sector_order <- c(
    sgb_order,
    function_present
  )

  n_sgb <- length(
    sgb_order
  )

  n_fun <- length(
    function_present
  )

  if (
    n_sgb == 0 ||
      n_fun == 0
  ) {

    plot.new()

    text(
      0.5,
      0.5,
      paste0(
        panel_title,
        "\nNo valid sectors"
      )
    )

    return(
      invisible(NULL)
    )
  }

  # Larger gap separates SGB and function blocks.
  gap_after <- c(
    rep(
      1.0,
      max(
        n_sgb - 1,
        0
      )
    ),
    8,
    rep(
      1.4,
      max(
        n_fun - 1,
        0
      )
    ),
    8
  )

  grid_col <- c(
    setNames(
      rep(
        "#333333",
        n_sgb
      ),
      sgb_order
    ),
    FUNCTION_COLORS[
      function_present
    ]
  )

  link_cols <- adjustcolor(
    FUNCTION_COLORS[
      links$Function_category
    ],
    alpha.f = 0.62
  )

  # Dynamic label size.
  sgb_cex <- if (
    small_panel
  ) {
    ifelse(
      n_sgb > 25,
      0.35,
      0.42
    )
  } else {
    ifelse(
      n_sgb > 30,
      0.48,
      ifelse(
        n_sgb > 20,
        0.55,
        0.63
      )
    )
  }

  function_cex <- if (
    small_panel
  ) {
    0.48
  } else {
    0.68
  }

  circos.clear()

  circos.par(
    start.degree = 88,
    gap.after = gap_after,
    track.margin = c(
      0.003,
      0.003
    ),
    cell.padding = c(
      0,
      0,
      0,
      0
    ),
    points.overflow.warning = FALSE
  )

  chordDiagram(
    x = links[
      ,
      .(
        SGB,
        Function_category,
        N_SV
      )
    ],
    order = sector_order,
    grid.col = grid_col,
    col = link_cols,
    transparency = 0,
    directional = 0,
    link.sort = TRUE,
    link.decreasing = TRUE,
    link.border = NA,
    
    annotationTrack = "grid",
    
    # Reduce the width of the colored outer ring
    annotationTrackHeight = mm_h(0.5),
    
    preAllocateTracks = list(
      track.height =
        ifelse(
          small_panel,
          0.12,
          0.15
        )
    )
  )

  # Outer labels.
  circos.trackPlotRegion(
    track.index = 1,
    bg.border = NA,
    panel.fun = function(
      x,
      y
    ) {

      sector <- CELL_META$sector.index

      xcenter <- CELL_META$xcenter

      ylim <- CELL_META$ylim

      if (
        sector %in%
          function_present
      ) {

        label <- FUNCTION_DISPLAY[[sector]]

        label_col <- FUNCTION_COLORS[[sector]]

        circos.text(
          x = xcenter,
          y = ylim[1],
          labels = label,
          facing = "clockwise",
          niceFacing = TRUE,
          adj = c(
            0,
            0.5
          ),
          cex = function_cex,
          font = 2,
          col = label_col
        )

      } else {

        label <- label_map[[sector]]

        if (
          is.null(label) ||
            is.na(label) ||
            label == ""
        ) {

          label <- sub(
            "^t__",
            "",
            sector
          )
        }

        circos.text(
          x = xcenter,
          y = ylim[1],
          labels = label,
          facing = "clockwise",
          niceFacing = TRUE,
          adj = c(
            0,
            0.5
          ),
          cex = sgb_cex,
          font = 3,
          col = "#222222"
        )
      }
    }
  )

  n_sv <- uniqueN(
    links[
      ,
      .(
        SGB,
        Function_category,
        N_SV
      )
    ]
  )

  title(
    main = panel_title,
    line = ifelse(
      small_panel,
      -1.5,
      -1.2
    ),
    cex.main = ifelse(
      small_panel,
      0.85,
      1.05
    ),
    font.main = 2
  )

  invisible(
    n_sv
  )
}

# ==============================================================================
# 6. Individual figures
# ==============================================================================

for (
  i in seq_len(
    nrow(
      PANELS
    )
  )
) {

  p <- PANELS[i]

  link_file <- file.path(
    INPUT_DIR,
    paste0(
      "circos_links_",
      p$BMI_file,
      "_",
      p$SV_type,
      ".tsv"
    )
  )

  if (
    !file.exists(
      link_file
    )
  ) {

    warning(
      "Missing input: ",
      link_file
    )

    next
  }

  links <- fread(
    link_file,
    check.names = FALSE
  )

  pdf_file <- file.path(
    OUT_DIR,
    paste0(
      "Circos_",
      p$BMI_file,
      "_",
      p$SV_type,
      "_",
      SUBSET_MODE,
      ".pdf"
    )
  )

  png_file <- file.path(
    OUT_DIR,
    paste0(
      "Circos_",
      p$BMI_file,
      "_",
      p$SV_type,
      "_",
      SUBSET_MODE,
      ".png"
    )
  )

  pdf(
    pdf_file,
    width = 9.5,
    height = 9.5,
    useDingbats = FALSE
  )

  par(
    mar = c(
      1,
      1,
      3,
      1
    )
  )

  plot_one_circos(
    links,
    p$Title,
    small_panel = FALSE
  )

  dev.off()


  png(
    png_file,
    width = 3000,
    height = 3000,
    res = 320
  )

  par(
    mar = c(
      1,
      1,
      3,
      1
    )
  )

  plot_one_circos(
    links,
    p$Title,
    small_panel = FALSE
  )

  dev.off()
}

# ==============================================================================
# 7. Combined 3 x 2 figure
# ==============================================================================

if (
  MAKE_COMBINED_6PANEL
) {

  combined_pdf <- file.path(
    OUT_DIR,
    paste0(
      "Circos_6panel_",
      SUBSET_MODE,
      ".pdf"
    )
  )

  pdf(
    combined_pdf,
    width = 14,
    height = 20,
    useDingbats = FALSE
  )

  par(
    mfrow = c(
      3,
      2
    ),
    mar = c(
      0.5,
      0.5,
      2.5,
      0.5
    )
  )

  for (
    i in seq_len(
      nrow(
        PANELS
      )
    )
  ) {

    p <- PANELS[i]

    link_file <- file.path(
      INPUT_DIR,
      paste0(
        "circos_links_",
        p$BMI_file,
        "_",
        p$SV_type,
        ".tsv"
      )
    )

    if (
      file.exists(
        link_file
      )
    ) {

      links <- fread(
        link_file,
        check.names = FALSE
      )

      plot_one_circos(
        links,
        p$Title,
        small_panel = TRUE
      )

    } else {

      plot.new()

      text(
        0.5,
        0.5,
        paste(
          p$Title,
          "input missing",
          sep = "\n"
        )
      )
    }
  }

  dev.off()
}

# ==============================================================================
# 8. Finish
# ==============================================================================

circos.clear()

cat(
  "\n========================================\n"
)

cat(
  "CIRCOS FIGURES COMPLETE\n"
)

cat(
  "========================================\n"
)

cat(
  "Subset: ",
  SUBSET_MODE,
  "\n",
  sep = ""
)

cat(
  "Output directory:\n",
  OUT_DIR,
  "\n",
  sep = ""
)
