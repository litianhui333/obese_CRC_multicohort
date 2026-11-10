# Windows/Google Drive compatibility launcher.
# readr removes the UTF-8 BOM; make.names reproduces base read.csv column naming.
suppressPackageStartupMessages(library(readr))

read.csv <- function(file, ...) {
    out <- as.data.frame(readr::read_csv(file, show_col_types = FALSE))
    names(out) <- make.names(names(out), unique = TRUE)
    out
}

source(
    "analyses/2026-05-20_da_feature_overall_crc_survival/rerun_optimalcutpoint_km_style.R",
    chdir = FALSE
)
