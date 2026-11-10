.libPaths(c("D:/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
})
setDTthreads(2L)

PROJECT <- "H:/My Drive/proj_ob_AA_CRC"
OUT <- file.path(PROJECT, "analyses/2026-10-13_SV_abundance_sensitivity")
OLD <- file.path(PROJECT, "analyses/2026-09-29_BMI_enriched_SV_function_modelQC/2026-09-30_BMI_specificity_layers")
RAW <- file.path(PROJECT, "02.processed_data/05_SV_result/all_BMI_crossBMI")
FILES <- c(
  candidates = file.path(OLD, "tables/L1_L2_candidate_master.tsv"),
  d = file.path(RAW, "all_BMI_dsgv.csv"),
  v = file.path(RAW, "all_BMI_vsgv.csv"),
  metadata = file.path(RAW, "metadata_allBMI_SGV_aligned.rds"),
  abundance = file.path(PROJECT, "02.processed_data/03_metaphlan_split_data/metaphlan_strain_renormalized_20260115.rds")
)
MODE <- commandArgs(TRUE)
MODE <- if (length(MODE)) MODE[1] else "manifest"
stopifnot(MODE %in% c("manifest", "pilot", "run"), all(file.exists(FILES)))
for (d in c("tables", "figures", "data_for_figures")) dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, stem, subdir = "tables") fwrite(x, file.path(OUT, subdir, paste0(stem, ".tsv")), sep = "\t", na = "NA")

hash_before <- tools::md5sum(FILES)
cand <- fread(FILES["candidates"])
stopifnot(nrow(cand) == 260L, !anyDuplicated(cand$BMI_SV_ID), all(c("BMI_SV_ID", "BMI_enriched_group", "SV_type", "SGB", "SV_feature", "L2_nominal", "L2_selected_BH") %in% names(cand)))
meta <- as.data.table(readRDS(FILES["metadata"]))
meta[, `:=`(
  SGV_ID = as.character(SGV_ID),
  Abundance_ID = as.character(Sample.Source.ID.processed),
  BMI_model = fcase(as.character(BMI_Group1) == "Normal", "Nw", as.character(BMI_Group1) == "Ow", "Ow", as.character(BMI_Group1) == "Ob", "Ob", default = NA_character_),
  Disease_model = as.character(Disease),
  Study_model = as.character(Study.name),
  Sex_model = as.character(Sex),
  Age_model = suppressWarnings(as.numeric(Age))
)]
stopifnot(!anyDuplicated(meta$SGV_ID), !anyDuplicated(meta$Abundance_ID))

ab <- readRDS(FILES["abundance"])
ab <- as.matrix(ab)
storage.mode(ab) <- "double"
ab_sgb <- sub(".*(t__SGB[0-9]+).*$", "\\1", rownames(ab))
stopifnot(!anyDuplicated(ab_sgb), all(unique(cand$SGB) %in% ab_sgb), all(is.finite(ab)))
rownames(ab) <- ab_sgb

read_sv <- function(role, type) {
  header <- names(fread(FILES[role], nrows = 0L))
  feats <- unique(cand[SV_type == type, SV_feature])
  stopifnot(all(feats %in% header))
  x <- fread(FILES[role], select = c(header[1L], feats), check.names = FALSE)
  setnames(x, 1L, "SGV_ID")
  x[, SGV_ID := as.character(SGV_ID)]
  stopifnot(!anyDuplicated(x$SGV_ID))
  x
}
dmat <- read_sv("d", "dSGV")
vmat <- read_sv("v", "vSGV")

coverage <- meta[, .(
  N_SV_metadata = .N,
  N_with_abundance = sum(Abundance_ID %in% colnames(ab)),
  N_CTRL = sum(Disease_model == "CTRL"),
  N_CRC = sum(Disease_model == "CRC"),
  N_CTRL_with_abundance = sum(Disease_model == "CTRL" & Abundance_ID %in% colnames(ab)),
  N_CRC_with_abundance = sum(Disease_model == "CRC" & Abundance_ID %in% colnames(ab))
), by = BMI_model][order(match(BMI_model, c("Nw", "Ow", "Ob")))]
manifest <- data.table(
  Metric = c("Candidate_SVs", "Candidate_SGBs", "Candidate_SGBs_matched", "SV_metadata_samples", "Abundance_samples", "SV_samples_with_abundance", "SV_samples_without_abundance", "Planned_paired_fits"),
  Value = c(nrow(cand), uniqueN(cand$SGB), sum(unique(cand$SGB) %in% rownames(ab)), nrow(meta), ncol(ab), sum(meta$Abundance_ID %in% colnames(ab)), sum(!meta$Abundance_ID %in% colnames(ab)), nrow(cand) * 3L * 2L)
)
write_tsv(manifest, "data_manifest")
write_tsv(coverage, "sample_abundance_coverage")
write_tsv(data.table(Path = unname(FILES), MD5 = unname(hash_before)), "input_hashes")
if (MODE == "manifest") {
  print(manifest)
  print(coverage)
  quit(save = "no", status = 0L)
}

capture_fit <- function(dat, type, adjusted) {
  warns <- character()
  result <- tryCatch(withCallingHandlers({
    form <- if (adjusted) Value ~ Disease + Age_z + Abundance_z + Sex + (1 | Study) else Value ~ Disease + Age_z + Sex + (1 | Study)
    if (type == "dSGV") {
      fit <- glmer(form, data = dat, family = binomial(), control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
      co <- coef(summary(fit))
      pcol <- "Pr(>|z|)"
    } else {
      fit <- lmerTest::lmer(form, data = dat, REML = FALSE, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
      co <- coef(summary(fit))
      pcol <- "Pr(>|t|)"
    }
    conv <- if (length(fit@optinfo$conv$lme4$messages)) paste(fit@optinfo$conv$lme4$messages, collapse = "; ") else ""
    abrow <- if (adjusted && "Abundance_z" %in% rownames(co)) co["Abundance_z", ] else rep(NA_real_, ncol(co))
    data.table(
      Disease_Beta = co["DiseaseCRC", "Estimate"],
      Disease_SE = co["DiseaseCRC", "Std. Error"],
      Disease_P = co["DiseaseCRC", pcol],
      Abundance_Beta = if (adjusted) abrow["Estimate"] else NA_real_,
      Abundance_SE = if (adjusted) abrow["Std. Error"] else NA_real_,
      Abundance_P = if (adjusted) abrow[pcol] else NA_real_,
      Singular = isSingular(fit, tol = 1e-4),
      Convergence_message = conv,
      LogLik = as.numeric(logLik(fit)),
      Error = ""
    )
  }, warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) data.table(
    Disease_Beta = NA_real_, Disease_SE = NA_real_, Disease_P = NA_real_,
    Abundance_Beta = NA_real_, Abundance_SE = NA_real_, Abundance_P = NA_real_,
    Singular = NA, Convergence_message = "", LogLik = NA_real_, Error = conditionMessage(e)
  ))
  result[, Warnings := paste(unique(warns), collapse = "; ")]
  result
}

fit_one <- function(z, g) {
  x <- if (z$SV_type == "dSGV") dmat else vmat
  dat <- data.table(
    SGV_ID = x$SGV_ID,
    Value = suppressWarnings(as.numeric(x[[z$SV_feature]]))
  )
  dat <- merge(dat, meta[, .(SGV_ID, Abundance_ID, BMI_model, Disease_model, Age_model, Sex_model, Study_model)], by = "SGV_ID", sort = FALSE)
  ai <- match(dat$Abundance_ID, colnames(ab))
  ri <- match(z$SGB, rownames(ab))
  dat[, Abundance_raw := ab[ri, ai]]
  dat <- dat[BMI_model == g & Disease_model %in% c("CTRL", "CRC") & is.finite(Value) & is.finite(Abundance_raw) & is.finite(Age_model) & !is.na(Sex_model) & nzchar(Sex_model) & !is.na(Study_model) & nzchar(Study_model)]
  dat[, `:=`(
    Disease = factor(Disease_model, levels = c("CTRL", "CRC")),
    Sex = factor(Sex_model),
    Study = factor(Study_model),
    Age_z = as.numeric(scale(Age_model)),
    Abundance_log1p = log1p(Abundance_raw)
  )]
  if (z$SV_type == "dSGV") dat[, Value := 1 - Value]
  abundance_sd <- sd(dat$Abundance_log1p)
  dat[, Abundance_z := if (is.finite(abundance_sd) && abundance_sd > 0) as.numeric(scale(Abundance_log1p)) else NA_real_]
  base <- data.table(
    BMI_SV_ID = z$BMI_SV_ID, BMI_enriched_group = z$BMI_enriched_group, BMI = g,
    SV_type = z$SV_type, SGB = z$SGB, SV_feature = z$SV_feature,
    N = nrow(dat), N_CTRL = sum(dat$Disease == "CTRL"), N_CRC = sum(dat$Disease == "CRC"),
    N_study = nlevels(droplevels(dat$Study)), Abundance_zero_N = sum(dat$Abundance_raw == 0),
    Abundance_raw_median = median(dat$Abundance_raw), Abundance_raw_max = max(dat$Abundance_raw),
    Abundance_log1p_SD = abundance_sd
  )
  eligible <- nrow(dat) >= 10L && min(base$N_CTRL, base$N_CRC) >= 5L && base$N_study >= 2L && is.finite(sd(dat$Age_z)) && abundance_sd > 0
  if (!eligible) {
    empty <- data.table(Disease_Beta = NA_real_, Disease_SE = NA_real_, Disease_P = NA_real_, Abundance_Beta = NA_real_, Abundance_SE = NA_real_, Abundance_P = NA_real_, Singular = NA, Convergence_message = "", LogLik = NA_real_, Error = "Not eligible after abundance matching", Warnings = "")
    return(rbind(cbind(base, Model = "matched_base", empty), cbind(base, Model = "abundance_adjusted", empty), fill = TRUE))
  }
  rbind(
    cbind(base, Model = "matched_base", capture_fit(dat, z$SV_type, FALSE)),
    cbind(base, Model = "abundance_adjusted", capture_fit(dat, z$SV_type, TRUE)),
    fill = TRUE
  )
}

targets <- if (MODE == "pilot") cand[order(BMI_enriched_group, SV_type, BMI_SV_ID)][c(1L, 2L, 90L, 91L, 160L, 161L)] else cand
targets <- targets[!is.na(BMI_SV_ID)]
started <- Sys.time()
res <- rbindlist(lapply(seq_len(nrow(targets)), function(i) {
  if (i %% 20L == 0L) message("Candidate ", i, "/", nrow(targets))
  rbindlist(lapply(c("Nw", "Ow", "Ob"), function(g) fit_one(targets[i], g)))
}), fill = TRUE)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
runtime <- data.table(Mode = MODE, Candidate_SVs = nrow(targets), Fits = nrow(res), Elapsed_seconds = elapsed, Seconds_per_fit = elapsed / nrow(res), Projected_full_seconds = elapsed / nrow(targets) * nrow(cand))
write_tsv(runtime, paste0(MODE, "_runtime"))
write_tsv(res, paste0(MODE, "_model_results_long"))
if (MODE == "pilot") {
  print(runtime)
  print(res[, .(Fits = .N, Errors = sum(nzchar(Error)), Singular = sum(Singular %in% TRUE), Warnings = sum(nzchar(Warnings))), by = Model])
  quit(save = "no", status = 0L)
}

valid <- res[, .(Valid = all(is.finite(Disease_Beta), is.finite(Disease_SE), Disease_SE > 0, is.finite(Disease_P)) && all(!nzchar(Error))), by = .(BMI_SV_ID, BMI)]
res <- merge(res, valid, by = c("BMI_SV_ID", "BMI"), all.x = TRUE, sort = FALSE)
wide <- dcast(res, BMI_SV_ID + BMI_enriched_group + BMI + SV_type + SGB + SV_feature + N + N_CTRL + N_CRC + N_study + Valid ~ Model, value.var = c("Disease_Beta", "Disease_SE", "Disease_P", "Abundance_Beta", "Abundance_SE", "Abundance_P", "Singular", "Warnings", "Error"))
wide[, `:=`(
  Disease_Beta_change = Disease_Beta_abundance_adjusted - Disease_Beta_matched_base,
  Disease_Beta_change_in_base_SE = (Disease_Beta_abundance_adjusted - Disease_Beta_matched_base) / Disease_SE_matched_base,
  Disease_abs_beta_ratio = abs(Disease_Beta_abundance_adjusted) / abs(Disease_Beta_matched_base),
  Disease_sign_changed = sign(Disease_Beta_abundance_adjusted) != sign(Disease_Beta_matched_base),
  Disease_P005_base = Disease_P_matched_base < 0.05,
  Disease_P005_adjusted = Disease_P_abundance_adjusted < 0.05
)]
write_tsv(wide, "model_comparison_by_BMI")

classify_candidates <- function(model_name) {
  bcol <- paste0("Disease_Beta_", model_name)
  scol <- paste0("Disease_SE_", model_name)
  pcol <- paste0("Disease_P_", model_name)
  rbindlist(lapply(unique(wide$BMI_SV_ID), function(id) {
    z <- wide[BMI_SV_ID == id][match(c("Nw", "Ow", "Ob"), BMI)]
    target <- match(z$BMI_enriched_group[1L], z$BMI)
    other <- setdiff(seq_len(3L), target)
    b <- z[[bcol]]; s <- z[[scol]]; p <- z[[pcol]]
    ok <- all(z$Valid) && all(is.finite(c(b, s, p))) && all(s > 0)
    if (!ok) return(data.table(BMI_SV_ID = id, Model_set = model_name, L1 = NA, L2_nominal = NA, L2_BH = NA, P_L2_joint = NA_real_, Target_Beta = NA_real_, Target_SE = NA_real_))
    effect_ok <- all(sign(b[other]) != sign(b[target]) | abs(b[other]) < abs(b[target]))
    l1 <- p[target] < 0.01 && all(p[other] >= 0.05) && effect_ok
    pdiff <- vapply(other, function(j) 2 * pnorm(abs((b[target] - b[j]) / sqrt(s[target]^2 + s[j]^2)), lower.tail = FALSE), numeric(1))
    data.table(BMI_SV_ID = id, Model_set = model_name, L1 = l1, L2_nominal = max(pdiff) < 0.05, L2_BH = FALSE, P_L2_joint = max(pdiff), Target_Beta = b[target], Target_SE = s[target])
  }))
}
cls <- rbind(classify_candidates("matched_base"), classify_candidates("abundance_adjusted"))
cls[, Q_L2_BH := p.adjust(P_L2_joint, method = "BH"), by = Model_set]
cls[, L2_BH := Q_L2_BH < 0.05]
clsw <- dcast(cls, BMI_SV_ID ~ Model_set, value.var = c("L1", "L2_nominal", "L2_BH", "P_L2_joint", "Q_L2_BH", "Target_Beta", "Target_SE"))
master <- merge(cand, clsw, by = "BMI_SV_ID", all.x = TRUE, sort = FALSE)
master[, `:=`(
  L1_adjustment_status = fcase(L1_matched_base %in% TRUE & L1_abundance_adjusted %in% TRUE, "Retained", L1_matched_base %in% TRUE & L1_abundance_adjusted %in% FALSE, "Lost after abundance adjustment", L1_matched_base %in% FALSE & L1_abundance_adjusted %in% TRUE, "Gained after abundance adjustment", default = "Absent in both"),
  L2_adjustment_status = fcase(L2_nominal_matched_base %in% TRUE & L2_nominal_abundance_adjusted %in% TRUE, "Retained", L2_nominal_matched_base %in% TRUE & L2_nominal_abundance_adjusted %in% FALSE, "Lost after abundance adjustment", L2_nominal_matched_base %in% FALSE & L2_nominal_abundance_adjusted %in% TRUE, "Gained after abundance adjustment", default = "Absent in both"),
  Target_beta_change = Target_Beta_abundance_adjusted - Target_Beta_matched_base,
  Target_abs_beta_ratio = abs(Target_Beta_abundance_adjusted) / abs(Target_Beta_matched_base),
  Target_sign_changed = sign(Target_Beta_abundance_adjusted) != sign(Target_Beta_matched_base)
)]
write_tsv(master, "candidate_abundance_sensitivity_master")

# Reproduce the previous strict target-BMI classification on the paired model sets.
# Target association uses P < 0.01; other-group association uses P < 0.05.
# Same-direction significant associations in another BMI group are classified as shared.
classify_specific_shared <- function(model_name) {
  bcol <- paste0("Disease_Beta_", model_name)
  pcol <- paste0("Disease_P_", model_name)
  rbindlist(lapply(seq_len(nrow(master)), function(i) {
    id <- master$BMI_SV_ID[i]
    g <- master$BMI_enriched_group[i]
    z <- wide[BMI_SV_ID == id][match(c("Nw", "Ow", "Ob"), BMI)]
    b <- z[[bcol]]
    p <- z[[pcol]]
    target <- match(g, z$BMI)
    other <- setdiff(seq_len(3L), target)
    sig005 <- p < 0.05
    same <- sign(b) == sign(b[target])
    effect_support <- all(sign(b[other]) != sign(b[target]) | abs(b[other]) < abs(b[target]))
    class <- if (p[target] >= 0.01) {
      "Target_not_P001"
    } else if (any(sig005[other] & same[other])) {
      "Cross_BMI_shared_same_direction"
    } else if (any(sig005[other] & !same[other])) {
      "Cross_BMI_divergent"
    } else if (effect_support) {
      "BMI_specific"
    } else {
      "Target_associated_not_effect_enriched"
    }
    data.table(
      BMI_SV_ID = id,
      BMI_enriched_group = g,
      SV_type = master$SV_type[i],
      SGB = master$SGB[i],
      Species = master$Species[i],
      Model_set = model_name,
      Specific_shared_class = class,
      Target_P = p[target],
      Target_Beta = b[target],
      P005_groups = paste(z$BMI[sig005], collapse = "+"),
      Same_direction_shared_with = paste(z$BMI[other][sig005[other] & same[other]], collapse = "+"),
      Opposite_direction_with = paste(z$BMI[other][sig005[other] & !same[other]], collapse = "+")
    )
  }))
}
specific_shared <- rbind(
  classify_specific_shared("matched_base"),
  classify_specific_shared("abundance_adjusted")
)
write_tsv(specific_shared, "specific_shared_classification")
specific_shared_summary <- specific_shared[, .(N_SV = .N), by = .(Model_set, BMI_enriched_group, SV_type, Specific_shared_class)]
setorder(specific_shared_summary, Model_set, BMI_enriched_group, SV_type, Specific_shared_class)
write_tsv(specific_shared_summary, "specific_shared_summary")

# P-value-only specificity requested for descriptive reporting. This deliberately
# ignores beta-difference/heterogeneity tests and effect-size ordering:
# target BMI P < 0.01 and both other BMI groups P > 0.05.
pvalue_only_specific <- copy(specific_shared[
  Specific_shared_class %chin% c("BMI_specific", "Target_associated_not_effect_enriched")
])
pvalue_only_specific[, Direction := fifelse(
  Target_Beta > 0, "Increased_in_CRC",
  fifelse(Target_Beta < 0, "Decreased_in_CRC", "Zero_effect")
)]
setorder(pvalue_only_specific, Model_set, BMI_enriched_group, SV_type, Direction, BMI_SV_ID)
write_tsv(pvalue_only_specific, "pvalue_only_specific_SVs")
pvalue_only_specific_summary <- pvalue_only_specific[, .(N_SV = .N),
  by = .(Model_set, BMI_enriched_group, SV_type, Direction)
]
setorder(pvalue_only_specific_summary, Model_set, BMI_enriched_group, SV_type, Direction)
write_tsv(pvalue_only_specific_summary, "pvalue_only_specific_direction_summary")
pvalue_only_specific_wide <- dcast(
  pvalue_only_specific_summary,
  Model_set + BMI_enriched_group + SV_type ~ Direction,
  value.var = "N_SV", fill = 0L
)
for (nm in c("Increased_in_CRC", "Decreased_in_CRC", "Zero_effect")) {
  if (!nm %in% names(pvalue_only_specific_wide)) pvalue_only_specific_wide[, (nm) := 0L]
}
pvalue_only_specific_wide[, Total := Increased_in_CRC + Decreased_in_CRC + Zero_effect]
setorder(pvalue_only_specific_wide, Model_set, BMI_enriched_group, SV_type)
write_tsv(pvalue_only_specific_wide, "pvalue_only_specific_summary")
specific_shared_exact <- specific_shared[
  Specific_shared_class == "Cross_BMI_shared_same_direction",
  .(N_SV = .N),
  by = .(Model_set, BMI_enriched_group, SV_type, Same_direction_shared_with)
]
setorder(specific_shared_exact, Model_set, BMI_enriched_group, SV_type, Same_direction_shared_with)
write_tsv(specific_shared_exact, "specific_shared_exact_summary")
specific_shared_wide <- dcast(
  specific_shared,
  BMI_SV_ID + BMI_enriched_group + SV_type + SGB + Species ~ Model_set,
  value.var = c("Specific_shared_class", "P005_groups", "Same_direction_shared_with", "Opposite_direction_with")
)
write_tsv(specific_shared_wide, "specific_shared_transitions")
write_tsv(
  specific_shared[Model_set == "abundance_adjusted" & Specific_shared_class == "Cross_BMI_shared_same_direction"],
  "abundance_adjusted_shared_SVs"
)

# Formal effect-difference specificity layer (L2): target-minus-each-other BMI.
# Failure to reject a difference is labelled "not distinguished", not equivalent/shared.
classify_effect_difference <- function(model_name) {
  bcol <- paste0("Disease_Beta_", model_name)
  scol <- paste0("Disease_SE_", model_name)
  rbindlist(lapply(seq_len(nrow(master)), function(i) {
    id <- master$BMI_SV_ID[i]
    g <- master$BMI_enriched_group[i]
    z <- wide[BMI_SV_ID == id][match(c("Nw", "Ow", "Ob"), BMI)]
    b <- z[[bcol]]
    s <- z[[scol]]
    target <- match(g, z$BMI)
    other <- setdiff(seq_len(3L), target)
    pdiff <- vapply(other, function(j) {
      2 * pnorm(abs((b[target] - b[j]) / sqrt(s[target]^2 + s[j]^2)), lower.tail = FALSE)
    }, numeric(1))
    names(pdiff) <- z$BMI[other]
    sig <- pdiff < 0.05
    class <- if (all(sig)) {
      "Effect_specific_both"
    } else if (sum(sig) == 1L) {
      paste0("Different_from_", names(pdiff)[sig], "; not_distinguished_from_", names(pdiff)[!sig])
    } else {
      "Not_distinguished_from_either"
    }
    data.table(
      BMI_SV_ID = id,
      BMI_enriched_group = g,
      SV_type = master$SV_type[i],
      SGB = master$SGB[i],
      Species = master$Species[i],
      Model_set = model_name,
      Effect_difference_class = class,
      Other_BMI_1 = names(pdiff)[1L],
      P_difference_1 = unname(pdiff[1L]),
      Other_BMI_2 = names(pdiff)[2L],
      P_difference_2 = unname(pdiff[2L]),
      P_L2_joint = max(pdiff)
    )
  }))
}
effect_difference <- rbind(
  classify_effect_difference("matched_base"),
  classify_effect_difference("abundance_adjusted")
)
effect_difference[, Q_L2_BH := p.adjust(P_L2_joint, method = "BH"), by = Model_set]
effect_difference[, Effect_specific_BH := Q_L2_BH < 0.05]
write_tsv(effect_difference, "effect_difference_specificity_classification")
effect_difference_summary <- effect_difference[, .(
  N_SV = .N,
  N_SGB = uniqueN(SGB),
  N_BH = sum(Effect_specific_BH)
), by = .(Model_set, BMI_enriched_group, SV_type, Effect_difference_class)]
setorder(effect_difference_summary, Model_set, BMI_enriched_group, SV_type, Effect_difference_class)
write_tsv(effect_difference_summary, "effect_difference_specificity_summary")
effect_difference_wide <- dcast(
  effect_difference,
  BMI_SV_ID + BMI_enriched_group + SV_type + SGB + Species ~ Model_set,
  value.var = c("Effect_difference_class", "P_L2_joint", "Q_L2_BH", "Effect_specific_BH")
)
write_tsv(effect_difference_wide, "effect_difference_specificity_transitions")

# Descriptive A/B decomposition of non-significant effect-difference contrasts.
# A: observed difference exceeds the prespecified moderate L3 margin but is imprecise.
# B: observed difference is within that margin. B is descriptive, not equivalence.
ab_pairwise <- rbindlist(lapply(seq_len(nrow(master)), function(i) {
  id <- master$BMI_SV_ID[i]
  g <- master$BMI_enriched_group[i]
  typ <- master$SV_type[i]
  z <- wide[BMI_SV_ID == id][match(c("Nw", "Ow", "Ob"), BMI)]
  b <- z$Disease_Beta_abundance_adjusted
  s <- z$Disease_SE_abundance_adjusted
  target <- match(g, z$BMI)
  other <- setdiff(seq_len(3L), target)
  margin <- if (typ == "dSGV") log(1.25) else 0.20 * master$Pooled_CTRL_SD[i]
  rbindlist(lapply(other, function(j) {
    delta <- b[target] - b[j]
    se_delta <- sqrt(s[target]^2 + s[j]^2)
    p_difference <- 2 * pnorm(abs(delta / se_delta), lower.tail = FALSE)
    data.table(
      BMI_SV_ID = id,
      BMI_enriched_group = g,
      Other_BMI = z$BMI[j],
      SV_type = typ,
      SGB = master$SGB[i],
      Species = master$Species[i],
      Delta = delta,
      SE_Delta = se_delta,
      Z_difference = delta / se_delta,
      P_difference = p_difference,
      Moderate_margin = margin,
      Absolute_delta_to_margin = abs(delta) / margin,
      AB_contrast_class = fcase(
        p_difference < 0.05, "Significant_difference",
        abs(delta) > margin, "A_large_point_difference_imprecise",
        default = "B_small_point_difference"
      )
    )
  }))
}))
write_tsv(ab_pairwise, "effect_difference_AB_pairwise")
ab_candidate <- rbindlist(lapply(split(ab_pairwise, ab_pairwise$BMI_SV_ID), function(z) {
  nonsig <- z[P_difference >= 0.05]
  class <- if (!nrow(nonsig)) {
    "Effect_specific_both"
  } else if (nrow(nonsig) == 1L) {
    if (nonsig$AB_contrast_class == "A_large_point_difference_imprecise") "One_nonsig_A" else "One_nonsig_B"
  } else {
    n_a <- sum(nonsig$AB_contrast_class == "A_large_point_difference_imprecise")
    if (n_a == 2L) "Both_nonsig_AA" else if (n_a == 1L) "Both_nonsig_AB" else "Both_nonsig_BB"
  }
  data.table(
    BMI_SV_ID = z$BMI_SV_ID[1L],
    BMI_enriched_group = z$BMI_enriched_group[1L],
    SV_type = z$SV_type[1L],
    SGB = z$SGB[1L],
    Species = z$Species[1L],
    Candidate_AB_class = class
  )
}))
write_tsv(ab_candidate, "effect_difference_AB_candidate")
ab_summary <- ab_candidate[, .(N_SV = .N, N_SGB = uniqueN(SGB)), by = .(BMI_enriched_group, SV_type, Candidate_AB_class)]
setorder(ab_summary, BMI_enriched_group, SV_type, Candidate_AB_class)
write_tsv(ab_summary, "effect_difference_AB_summary")

summary <- master[, .(
  N_candidates = .N,
  Original_L2_nominal = sum(L2_nominal %in% TRUE),
  Matched_base_L1 = sum(L1_matched_base %in% TRUE, na.rm = TRUE),
  Adjusted_L1 = sum(L1_abundance_adjusted %in% TRUE, na.rm = TRUE),
  Matched_base_L2 = sum(L2_nominal_matched_base %in% TRUE, na.rm = TRUE),
  Adjusted_L2 = sum(L2_nominal_abundance_adjusted %in% TRUE, na.rm = TRUE),
  Adjusted_L2_BH = sum(L2_BH_abundance_adjusted %in% TRUE, na.rm = TRUE),
  Original_L2_retained_after_adjustment = sum(L2_nominal %in% TRUE & L2_nominal_abundance_adjusted %in% TRUE, na.rm = TRUE),
  Original_L2_lost_after_adjustment = sum(L2_nominal %in% TRUE & L2_nominal_abundance_adjusted %in% FALSE, na.rm = TRUE),
  Adjusted_L2_new_vs_original = sum(L2_nominal %in% FALSE & L2_nominal_abundance_adjusted %in% TRUE, na.rm = TRUE),
  L2_retained = sum(L2_adjustment_status == "Retained", na.rm = TRUE),
  L2_lost = sum(L2_adjustment_status == "Lost after abundance adjustment", na.rm = TRUE),
  L2_gained = sum(L2_adjustment_status == "Gained after abundance adjustment", na.rm = TRUE),
  Median_target_abs_beta_ratio = median(Target_abs_beta_ratio[is.finite(Target_abs_beta_ratio)], na.rm = TRUE),
  Target_sign_changes = sum(Target_sign_changed %in% TRUE, na.rm = TRUE)
), by = .(BMI_enriched_group, SV_type)]
write_tsv(summary, "sensitivity_summary_by_BMI_SVtype")
write_tsv(master[L2_nominal %in% TRUE | L2_nominal_matched_base %in% TRUE | L2_nominal_abundance_adjusted %in% TRUE], "L2_candidate_changes")
overall <- data.table(
  Metric = c(
    "Candidate_SVs", "Original_L2_nominal", "Matched_base_L2_nominal", "Abundance_adjusted_L2_nominal",
    "Original_L2_retained_after_adjustment", "Original_L2_lost_after_adjustment", "Adjusted_L2_new_vs_original",
    "Adjusted_L2_BH", "dSGV_target_beta_Pearson", "vSGV_target_beta_Pearson",
    "dSGV_target_beta_median_abs_ratio", "vSGV_target_beta_median_abs_ratio", "Target_effect_sign_changes"
  ),
  Value = c(
    nrow(master), sum(master$L2_nominal %in% TRUE), sum(master$L2_nominal_matched_base %in% TRUE), sum(master$L2_nominal_abundance_adjusted %in% TRUE),
    sum(master$L2_nominal %in% TRUE & master$L2_nominal_abundance_adjusted %in% TRUE),
    sum(master$L2_nominal %in% TRUE & master$L2_nominal_abundance_adjusted %in% FALSE),
    sum(master$L2_nominal %in% FALSE & master$L2_nominal_abundance_adjusted %in% TRUE),
    sum(master$L2_BH_abundance_adjusted %in% TRUE),
    cor(master[SV_type == "dSGV", Target_Beta_matched_base], master[SV_type == "dSGV", Target_Beta_abundance_adjusted]),
    cor(master[SV_type == "vSGV", Target_Beta_matched_base], master[SV_type == "vSGV", Target_Beta_abundance_adjusted]),
    median(master[SV_type == "dSGV" & is.finite(Target_abs_beta_ratio), Target_abs_beta_ratio]),
    median(master[SV_type == "vSGV" & is.finite(Target_abs_beta_ratio), Target_abs_beta_ratio]),
    sum(master$Target_sign_changed %in% TRUE)
  )
)
write_tsv(overall, "sensitivity_summary_overall")

plotdat <- master[, .(BMI_SV_ID, BMI_enriched_group, SV_type, SGB, Original_L2 = L2_nominal, L2_adjustment_status, Target_Beta_matched_base, Target_Beta_abundance_adjusted, Target_abs_beta_ratio, Target_sign_changed)]
write_tsv(plotdat, "Fig1_target_effect_before_after", "data_for_figures")
cols <- c("Retained" = "#2B8CBE", "Lost after abundance adjustment" = "#D95F0E", "Gained after abundance adjustment" = "#31A354", "Absent in both" = "#BDBDBD")
draw_scatter <- function() {
  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(5, 0.7))
  op <- par(mar = c(4.2, 4.2, 2.2, 1), family = "sans")
  on.exit(par(op))
  for (typ in c("dSGV", "vSGV")) {
    z <- plotdat[SV_type == typ]
    lim <- range(c(z$Target_Beta_matched_base, z$Target_Beta_abundance_adjusted), finite = TRUE)
    plot(z$Target_Beta_matched_base, z$Target_Beta_abundance_adjusted, pch = 21, bg = cols[z$L2_adjustment_status], col = "white", cex = 0.85, xlab = "Matched-base CRC effect", ylab = "Abundance-adjusted CRC effect", main = typ, xlim = lim, ylim = lim)
    abline(0, 1, lty = 2, col = "#555555")
    abline(h = 0, v = 0, col = "#DDDDDD")
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", horiz = TRUE, legend = names(cols), pt.bg = cols, pch = 21, col = "white", bty = "n", cex = 0.78)
}
pdf(file.path(OUT, "figures/Fig1_target_effect_before_after.pdf"), width = 8.2, height = 4.5, family = "sans")
draw_scatter(); dev.off()
png(file.path(OUT, "figures/Fig1_target_effect_before_after.png"), width = 8.2, height = 4.5, units = "in", res = 300)
draw_scatter(); dev.off()

qc <- res[, .(N_fits = .N, N_errors = sum(nzchar(Error)), N_warnings = sum(nzchar(Warnings)), N_singular = sum(Singular %in% TRUE), N_valid = sum(Valid %in% TRUE)), by = Model]
write_tsv(qc, "model_QC_summary")
stopifnot(identical(unname(tools::md5sum(FILES)), unname(hash_before)))
capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))
writeLines(list.files(OUT, recursive = TRUE, full.names = TRUE), file.path(OUT, "output_files.txt"))
writeLines(c("COMPLETE", paste("Finished:", Sys.time()), "Original inputs unchanged by MD5.", "Candidate-only abundance sensitivity; no genome-wide rediscovery claim."), file.path(OUT, "COMPLETE.txt"))
print(summary)
print(qc)
