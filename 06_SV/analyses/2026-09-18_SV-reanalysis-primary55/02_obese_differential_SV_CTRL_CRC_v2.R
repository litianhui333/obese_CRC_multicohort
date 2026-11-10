# ============================================================
# Obese CTRL vs CRC differential-SV analysis
# PRIMARY55 SGBs
#
# Primary analysis:
#   All successfully fitted models with valid P values
#   are included, regardless of Model_status.
#
# Strict QC analysis:
#   Only models with Model_status == "OK"
#
# dSGV:
#   Deletion = 1 - dSGV
#   binomial GLMM:
#   Deletion ~ Disease + Age + Sex + (1 | Study)
#
# vSGV:
#   linear mixed model:
#   Value ~ Disease + Age + Sex + (1 | Study)
#
# Candidate:
#   P < 0.01 OR Q_type < 0.10
#
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
})


# ============================================================
# 1. Load standardized PRIMARY55 inputs
# ============================================================

source(
  paste0(
    "H:/My Drive/proj_ob_AA_CRC/analyses/",
    "2026-09-18_SV-reanalysis-primary55/",
    "00_source_PRIMARY55_inputs.R"
  )
)


# ============================================================
# 2. Analysis parameters
# ============================================================

MIN_N_CTRL <- 20
MIN_N_CRC <- 20
MIN_INFORMATIVE_STUDIES <- 3

MIN_DSGV_MINOR_STATE <- 5

P_CANDIDATE <- 0.01
Q_CANDIDATE <- 0.10
Q_SIGNIFICANT <- 0.05


PARAM_FILE <- file.path(
  ROOT,
  "00_analysis_parameters_CTRL_CRC.tsv"
)

params <- data.table(
  Parameter = c(
    "Population",
    "Comparison",
    "Primary_SGB_number",
    "Minimum_CTRL_N",
    "Minimum_CRC_N",
    "Minimum_informative_studies",
    "Minimum_dSGV_minor_state",
    "Candidate_P_threshold",
    "Candidate_Q_type_threshold",
    "FDR_significant_threshold",
    "Primary_model_policy",
    "Strict_QC_policy"
  ),
  
  Value = c(
    "Obese",
    "CTRL_vs_CRC",
    55,
    MIN_N_CTRL,
    MIN_N_CRC,
    MIN_INFORMATIVE_STUDIES,
    MIN_DSGV_MINOR_STATE,
    P_CANDIDATE,
    Q_CANDIDATE,
    Q_SIGNIFICANT,
    "All fitted models with valid P values",
    "Model_status == OK"
  )
)

fwrite(
  params,
  PARAM_FILE,
  sep = "\t"
)


# ============================================================
# 3. Utility: convergence information
# ============================================================

get_fit_diagnostics <- function(fit) {
  
  singular <- tryCatch(
    lme4::isSingular(
      fit,
      tol = 1e-4
    ),
    error = function(e) NA
  )
  
  
  # Optimizer convergence code
  opt_code <- tryCatch(
    fit@optinfo$conv$opt,
    error = function(e) NULL
  )
  
  
  optimizer_ok <- TRUE
  
  if (!is.null(opt_code)) {
    
    opt_values <- suppressWarnings(
      as.numeric(
        unlist(opt_code)
      )
    )
    
    opt_values <- opt_values[
      !is.na(opt_values)
    ]
    
    if (
      length(opt_values) > 0
    ) {
      
      optimizer_ok <- all(
        opt_values == 0
      )
    }
  }
  
  
  # lme4 convergence / boundary messages
  msgs <- tryCatch(
    fit@optinfo$conv$lme4$messages,
    error = function(e) NULL
  )
  
  
  if (
    is.null(msgs) ||
    length(msgs) == 0
  ) {
    
    model_message <- NA_character_
    
  } else {
    
    model_message <- paste(
      msgs,
      collapse = " | "
    )
  }
  
  
  # Strict status:
  # any lme4 message means not "OK"
  # This intentionally keeps a strict QC definition.
  if (
    optimizer_ok &&
    is.na(model_message)
  ) {
    
    model_status <- "OK"
    
  } else if (
    optimizer_ok &&
    !is.na(model_message)
  ) {
    
    model_status <- model_message
    
  } else {
    
    model_status <- paste0(
      "Optimizer_problem",
      ifelse(
        is.na(model_message),
        "",
        paste0(
          " | ",
          model_message
        )
      )
    )
  }
  
  
  list(
    Singular = singular,
    Optimizer_OK = optimizer_ok,
    Model_message = model_message,
    Model_status = model_status
  )
}


# ============================================================
# 4. Utility: prepare feature data
# ============================================================

make_feature_data <- function(
    matrix_dt,
    sample_ids,
    feature
) {
  
  x <- data.table(
    Sample_ID = sample_ids,
    Value = matrix_dt[[feature]]
  )
  
  z <- merge(
    m,
    x,
    by = "Sample_ID",
    all = FALSE,
    sort = FALSE
  )
  
  
  z <- z[
    !is.na(Value) &
      !is.na(Disease) &
      !is.na(Age) &
      !is.na(Sex) &
      !is.na(Study)
  ]
  
  
  z <- as.data.frame(z)
  
  z$Disease <- factor(
    z$Disease,
    levels = c(
      "CTRL",
      "CRC"
    )
  )
  
  z$Sex <- droplevels(
    factor(z$Sex)
  )
  
  z$Study <- droplevels(
    factor(z$Study)
  )
  
  z
}


# ============================================================
# 5. Utility: feature eligibility
# ============================================================

check_feature_eligibility <- function(
    dat,
    sv_type
) {
  
  n_ctrl <- sum(
    dat$Disease == "CTRL"
  )
  
  n_crc <- sum(
    dat$Disease == "CRC"
  )
  
  
  # Informative study:
  # at least one CTRL and one CRC
  # among non-missing samples
  informative_studies <- 0L
  
  if (
    nrow(dat) > 0
  ) {
    
    by_study <- table(
      dat$Study,
      dat$Disease
    )
    
    if (
      all(
        c(
          "CTRL",
          "CRC"
        ) %in%
        colnames(by_study)
      )
    ) {
      
      informative_studies <- sum(
        by_study[, "CTRL"] > 0 &
          by_study[, "CRC"] > 0
      )
    }
  }
  
  
  reasons <- character(0)
  
  
  if (
    n_ctrl < MIN_N_CTRL
  ) {
    reasons <- c(
      reasons,
      "N_CTRL_below_threshold"
    )
  }
  
  
  if (
    n_crc < MIN_N_CRC
  ) {
    reasons <- c(
      reasons,
      "N_CRC_below_threshold"
    )
  }
  
  
  if (
    informative_studies <
    MIN_INFORMATIVE_STUDIES
  ) {
    
    reasons <- c(
      reasons,
      "Too_few_informative_studies"
    )
  }
  
  
  if (
    sv_type == "dSGV"
  ) {
    
    value_tab <- table(
      dat$Value
    )
    
    if (
      !all(
        unique(dat$Value) %in%
        c(0, 1)
      )
    ) {
      
      reasons <- c(
        reasons,
        "dSGV_not_binary"
      )
      
    } else {
      
      n0 <- sum(
        dat$Value == 0
      )
      
      n1 <- sum(
        dat$Value == 1
      )
      
      minor_state <- min(
        n0,
        n1
      )
      
      if (
        minor_state <
        MIN_DSGV_MINOR_STATE
      ) {
        
        reasons <- c(
          reasons,
          "dSGV_minor_state_below_threshold"
        )
      }
    }
  }
  
  
  if (
    sv_type == "vSGV"
  ) {
    
    value_sd <- sd(
      dat$Value,
      na.rm = TRUE
    )
    
    if (
      is.na(value_sd) ||
      value_sd == 0
    ) {
      
      reasons <- c(
        reasons,
        "vSGV_zero_variance"
      )
    }
  }
  
  
  eligible <- length(
    reasons
  ) == 0
  
  
  list(
    Eligible = eligible,
    
    Eligibility_reason =
      ifelse(
        eligible,
        "Eligible",
        paste(
          reasons,
          collapse = ";"
        )
      ),
    
    N = nrow(dat),
    N_CTRL = n_ctrl,
    N_CRC = n_crc,
    
    N_informative_studies =
      informative_studies
  )
}


# ============================================================
# 6. Fit one dSGV
# ============================================================

fit_dSGV <- function(
    feature,
    matrix_dt,
    sample_ids,
    feature_info
) {
  
  dat <- make_feature_data(
    matrix_dt,
    sample_ids,
    feature
  )
  
  
  eli <- check_feature_eligibility(
    dat,
    "dSGV"
  )
  
  
  base_result <- data.table(
    SGB = feature_info$SGB,
    SV_type = "dSGV",
    SV_feature = feature,
    
    N = eli$N,
    N_CTRL = eli$N_CTRL,
    N_CRC = eli$N_CRC,
    
    N_informative_studies =
      eli$N_informative_studies,
    
    Eligible =
      eli$Eligible,
    
    Eligibility_reason =
      eli$Eligibility_reason
  )
  
  
  if (
    !eli$Eligible
  ) {
    
    return(
      cbind(
        base_result,
        data.table(
          CTRL_value = NA_real_,
          CRC_value = NA_real_,
          
          CTRL_deletion_rate =
            NA_real_,
          
          CRC_deletion_rate =
            NA_real_,
          
          Beta = NA_real_,
          SE = NA_real_,
          OR = NA_real_,
          P = NA_real_,
          
          Singular = NA,
          Optimizer_OK = NA,
          Model_message =
            NA_character_,
          
          Model_status =
            "Not_eligible"
        )
      )
    )
  }
  
  
  # dSGV:
  # 1 = region retained
  # 0 = region deleted
  #
  # For intuitive interpretation:
  # Deletion = 1 - dSGV
  dat$Deletion <- 1 -
    dat$Value
  
  
  ctrl_value <- mean(
    dat$Value[
      dat$Disease == "CTRL"
    ]
  )
  
  crc_value <- mean(
    dat$Value[
      dat$Disease == "CRC"
    ]
  )
  
  
  ctrl_del <- mean(
    dat$Deletion[
      dat$Disease == "CTRL"
    ]
  )
  
  crc_del <- mean(
    dat$Deletion[
      dat$Disease == "CRC"
    ]
  )
  
  
  fit <- tryCatch(
    
    glmer(
      Deletion ~
        Disease +
        Age +
        Sex +
        (1 | Study),
      
      data = dat,
      
      family = binomial(),
      
      control = glmerControl(
        optimizer = "bobyqa",
        optCtrl = list(
          maxfun = 2e5
        )
      )
    ),
    
    error = function(e) e
  )
  
  
  if (
    inherits(
      fit,
      "error"
    )
  ) {
    
    return(
      cbind(
        base_result,
        data.table(
          CTRL_value =
            ctrl_value,
          
          CRC_value =
            crc_value,
          
          CTRL_deletion_rate =
            ctrl_del,
          
          CRC_deletion_rate =
            crc_del,
          
          Beta = NA_real_,
          SE = NA_real_,
          OR = NA_real_,
          P = NA_real_,
          
          Singular = NA,
          Optimizer_OK = FALSE,
          
          Model_message =
            conditionMessage(fit),
          
          Model_status =
            paste0(
              "Model_error: ",
              conditionMessage(fit)
            )
        )
      )
    )
  }
  
  
  cf <- coef(
    summary(fit)
  )
  
  
  disease_row <- grep(
    "^DiseaseCRC$",
    rownames(cf)
  )
  
  
  if (
    length(disease_row) != 1
  ) {
    
    return(
      cbind(
        base_result,
        data.table(
          CTRL_value =
            ctrl_value,
          
          CRC_value =
            crc_value,
          
          CTRL_deletion_rate =
            ctrl_del,
          
          CRC_deletion_rate =
            crc_del,
          
          Beta = NA_real_,
          SE = NA_real_,
          OR = NA_real_,
          P = NA_real_,
          
          Singular = NA,
          Optimizer_OK = NA,
          
          Model_message =
            "Disease coefficient not identified",
          
          Model_status =
            "Coefficient_error"
        )
      )
    )
  }
  
  
  beta <- cf[
    disease_row,
    "Estimate"
  ]
  
  se <- cf[
    disease_row,
    "Std. Error"
  ]
  
  p <- cf[
    disease_row,
    "Pr(>|z|)"
  ]
  
  
  diagnostics <- get_fit_diagnostics(
    fit
  )
  
  
  cbind(
    base_result,
    
    data.table(
      CTRL_value =
        ctrl_value,
      
      CRC_value =
        crc_value,
      
      CTRL_deletion_rate =
        ctrl_del,
      
      CRC_deletion_rate =
        crc_del,
      
      Beta = beta,
      SE = se,
      
      OR = exp(beta),
      
      P = p,
      
      Singular =
        diagnostics$Singular,
      
      Optimizer_OK =
        diagnostics$Optimizer_OK,
      
      Model_message =
        diagnostics$Model_message,
      
      Model_status =
        diagnostics$Model_status
    )
  )
}


# ============================================================
# 7. Fit one vSGV
# ============================================================

fit_vSGV <- function(
    feature,
    matrix_dt,
    sample_ids,
    feature_info
) {
  
  dat <- make_feature_data(
    matrix_dt,
    sample_ids,
    feature
  )
  
  
  eli <- check_feature_eligibility(
    dat,
    "vSGV"
  )
  
  
  base_result <- data.table(
    SGB = feature_info$SGB,
    SV_type = "vSGV",
    SV_feature = feature,
    
    N = eli$N,
    N_CTRL = eli$N_CTRL,
    N_CRC = eli$N_CRC,
    
    N_informative_studies =
      eli$N_informative_studies,
    
    Eligible =
      eli$Eligible,
    
    Eligibility_reason =
      eli$Eligibility_reason
  )
  
  
  if (
    !eli$Eligible
  ) {
    
    return(
      cbind(
        base_result,
        data.table(
          CTRL_value = NA_real_,
          CRC_value = NA_real_,
          
          CTRL_deletion_rate =
            NA_real_,
          
          CRC_deletion_rate =
            NA_real_,
          
          Beta = NA_real_,
          SE = NA_real_,
          OR = NA_real_,
          P = NA_real_,
          
          Singular = NA,
          Optimizer_OK = NA,
          
          Model_message =
            NA_character_,
          
          Model_status =
            "Not_eligible"
        )
      )
    )
  }
  
  
  ctrl_value <- mean(
    dat$Value[
      dat$Disease == "CTRL"
    ],
    na.rm = TRUE
  )
  
  crc_value <- mean(
    dat$Value[
      dat$Disease == "CRC"
    ],
    na.rm = TRUE
  )
  
  
  fit <- tryCatch(
    
    lmerTest::lmer(
      Value ~
        Disease +
        Age +
        Sex +
        (1 | Study),
      
      data = dat,
      
      REML = FALSE
    ),
    
    error = function(e) e
  )
  
  
  if (
    inherits(
      fit,
      "error"
    )
  ) {
    
    return(
      cbind(
        base_result,
        data.table(
          CTRL_value =
            ctrl_value,
          
          CRC_value =
            crc_value,
          
          CTRL_deletion_rate =
            NA_real_,
          
          CRC_deletion_rate =
            NA_real_,
          
          Beta = NA_real_,
          SE = NA_real_,
          OR = NA_real_,
          P = NA_real_,
          
          Singular = NA,
          Optimizer_OK = FALSE,
          
          Model_message =
            conditionMessage(fit),
          
          Model_status =
            paste0(
              "Model_error: ",
              conditionMessage(fit)
            )
        )
      )
    )
  }
  
  
  cf <- coef(
    summary(fit)
  )
  
  
  disease_row <- grep(
    "^DiseaseCRC$",
    rownames(cf)
  )
  
  
  if (
    length(disease_row) != 1
  ) {
    
    return(
      cbind(
        base_result,
        data.table(
          CTRL_value =
            ctrl_value,
          
          CRC_value =
            crc_value,
          
          CTRL_deletion_rate =
            NA_real_,
          
          CRC_deletion_rate =
            NA_real_,
          
          Beta = NA_real_,
          SE = NA_real_,
          OR = NA_real_,
          P = NA_real_,
          
          Singular = NA,
          Optimizer_OK = NA,
          
          Model_message =
            "Disease coefficient not identified",
          
          Model_status =
            "Coefficient_error"
        )
      )
    )
  }
  
  
  beta <- cf[
    disease_row,
    "Estimate"
  ]
  
  se <- cf[
    disease_row,
    "Std. Error"
  ]
  
  p_col <- grep(
    "^Pr\\(",
    colnames(cf),
    value = TRUE
  )
  
  
  if (
    length(p_col) >= 1
  ) {
    
    p <- cf[
      disease_row,
      p_col[1]
    ]
    
  } else {
    
    p <- NA_real_
  }
  
  
  diagnostics <- get_fit_diagnostics(
    fit
  )
  
  
  cbind(
    base_result,
    
    data.table(
      CTRL_value =
        ctrl_value,
      
      CRC_value =
        crc_value,
      
      CTRL_deletion_rate =
        NA_real_,
      
      CRC_deletion_rate =
        NA_real_,
      
      Beta = beta,
      SE = se,
      
      OR = NA_real_,
      
      P = p,
      
      Singular =
        diagnostics$Singular,
      
      Optimizer_OK =
        diagnostics$Optimizer_OK,
      
      Model_message =
        diagnostics$Model_message,
      
      Model_status =
        diagnostics$Model_status
    )
  )
}


# ============================================================
# 8. Feature lists
# ============================================================

d_features <- colnames(d)[-1]
v_features <- colnames(v)[-1]


cat(
  "\ndSGV features:",
  length(d_features),
  "\n"
)

cat(
  "vSGV features:",
  length(v_features),
  "\n"
)


# ============================================================
# 9. Feature map lookup
# ============================================================

d_map <- feature_map[
  SV_type == "dSGV"
]

v_map <- feature_map[
  SV_type == "vSGV"
]


# ============================================================
# 10. Run dSGV models
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "Running dSGV models\n"
)

cat(
  "========================================\n"
)


d_results_list <- vector(
  "list",
  length(d_features)
)


for (
  i in seq_along(d_features)
) {
  
  feature <- d_features[i]
  
  
  info <- d_map[
    SV_feature == feature
  ]
  
  
  if (
    nrow(info) == 0
  ) {
    
    stop(
      "Cannot map dSGV feature: ",
      feature
    )
  }
  
  
  d_results_list[[i]] <-
    fit_dSGV(
      feature = feature,
      matrix_dt = d,
      sample_ids = d_id,
      feature_info = info[1]
    )
  
  
  if (
    i %% 100 == 0 ||
    i == length(d_features)
  ) {
    
    cat(
      "dSGV:",
      i,
      "/",
      length(d_features),
      "\n"
    )
  }
  
  
  if (
    i %% 500 == 0
  ) {
    
    checkpoint <- rbindlist(
      d_results_list[
        seq_len(i)
      ],
      fill = TRUE
    )
    
    fwrite(
      checkpoint,
      file.path(
        TABLE_DIR,
        "checkpoint_dSGV.tsv.gz"
      ),
      sep = "\t"
    )
  }
}


d_results <- rbindlist(
  d_results_list,
  fill = TRUE
)


# ============================================================
# 11. Run vSGV models
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "Running vSGV models\n"
)

cat(
  "========================================\n"
)


v_results_list <- vector(
  "list",
  length(v_features)
)


for (
  i in seq_along(v_features)
) {
  
  feature <- v_features[i]
  
  
  info <- v_map[
    SV_feature == feature
  ]
  
  
  if (
    nrow(info) == 0
  ) {
    
    stop(
      "Cannot map vSGV feature: ",
      feature
    )
  }
  
  
  v_results_list[[i]] <-
    fit_vSGV(
      feature = feature,
      matrix_dt = v,
      sample_ids = v_id,
      feature_info = info[1]
    )
  
  
  if (
    i %% 100 == 0 ||
    i == length(v_features)
  ) {
    
    cat(
      "vSGV:",
      i,
      "/",
      length(v_features),
      "\n"
    )
  }
  
  
  if (
    i %% 500 == 0
  ) {
    
    checkpoint <- rbindlist(
      v_results_list[
        seq_len(i)
      ],
      fill = TRUE
    )
    
    fwrite(
      checkpoint,
      file.path(
        TABLE_DIR,
        "checkpoint_vSGV.tsv.gz"
      ),
      sep = "\t"
    )
  }
}


v_results <- rbindlist(
  v_results_list,
  fill = TRUE
)


# ============================================================
# 12. Combine
# ============================================================

all_res <- rbindlist(
  list(
    d_results,
    v_results
  ),
  fill = TRUE
)


# ============================================================
# 13. Define TWO model sets
# ============================================================

# ------------------------------------------------------------
# Primary:
# Any eligible model that produced a valid P value
#
# This does NOT require Model_status == "OK".
# ------------------------------------------------------------

all_res[
  ,
  Primary_testable :=
    Eligible == TRUE &
    is.finite(P)
]


# ------------------------------------------------------------
# Strict QC:
# Only models with no lme4 warning/message
# ------------------------------------------------------------

all_res[
  ,
  Strict_OK :=
    Primary_testable &
    Model_status == "OK"
]


# ============================================================
# 14. Primary FDR
#
# IMPORTANT:
# This is the main result.
# Model_status is NOT used as a filter.
# ============================================================

all_res[
  ,
  `:=`(
    Q_type = NA_real_,
    Q_all = NA_real_
  )
]


for (
  typ in c(
    "dSGV",
    "vSGV"
  )
) {
  
  idx <- which(
    all_res$SV_type == typ &
      all_res$Primary_testable
  )
  
  
  all_res$Q_type[idx] <-
    p.adjust(
      all_res$P[idx],
      method = "BH"
    )
}


idx_all <- which(
  all_res$Primary_testable
)


all_res$Q_all[idx_all] <-
  p.adjust(
    all_res$P[idx_all],
    method = "BH"
  )


# ============================================================
# 15. Strict-OK FDR
#
# Separate QC/sensitivity result
# ============================================================

all_res[
  ,
  `:=`(
    Q_type_strict_OK =
      NA_real_,
    
    Q_all_strict_OK =
      NA_real_
  )
]


for (
  typ in c(
    "dSGV",
    "vSGV"
  )
) {
  
  idx <- which(
    all_res$SV_type == typ &
      all_res$Strict_OK
  )
  
  
  all_res$Q_type_strict_OK[idx] <-
    p.adjust(
      all_res$P[idx],
      method = "BH"
    )
}


idx_strict <- which(
  all_res$Strict_OK
)


all_res$Q_all_strict_OK[
  idx_strict
] <-
  p.adjust(
    all_res$P[idx_strict],
    method = "BH"
  )


# ============================================================
# 16. Candidate definitions
# ============================================================

# ------------------------------------------------------------
# PRIMARY candidate
#
# No Model_status filter
# ------------------------------------------------------------

all_res[
  ,
  Candidate :=
    Primary_testable &
    (
      P < P_CANDIDATE |
        Q_type < Q_CANDIDATE
    )
]


# ------------------------------------------------------------
# STRICT candidate
#
# Model_status must be OK
# ------------------------------------------------------------

all_res[
  ,
  Candidate_strict_OK :=
    Strict_OK &
    (
      P < P_CANDIDATE |
        Q_type_strict_OK <
        Q_CANDIDATE
    )
]


# Additional labels
all_res[
  ,
  Nominal05 :=
    Primary_testable &
    P < 0.05
]


all_res[
  ,
  FDR10 :=
    Primary_testable &
    Q_type < 0.10
]


all_res[
  ,
  FDR05 :=
    Primary_testable &
    Q_type < 0.05
]


# ============================================================
# 17. Model-category annotation
# ============================================================

all_res[
  ,
  Model_category :=
    fifelse(
      Model_status == "OK",
      "OK",
      fifelse(
        Singular == TRUE &
          Optimizer_OK == TRUE,
        "Singular_or_boundary",
        fifelse(
          Primary_testable,
          "Other_warning_but_estimated",
          "Not_testable"
        )
      )
    )
]


# ============================================================
# 18. Primary candidate table
# ============================================================

candidates <- all_res[
  Candidate == TRUE
][
  order(
    P
  )
]


# ============================================================
# 19. Strict candidate table
# ============================================================

candidates_strict <- all_res[
  Candidate_strict_OK == TRUE
][
  order(
    P
  )
]


# ============================================================
# 20. FDR-significant table
# ============================================================

fdr05 <- all_res[
  FDR05 == TRUE
][
  order(
    Q_type,
    P
  )
]


# ============================================================
# 21. Eligibility table
# ============================================================

eligibility <- all_res[
  ,
  .(
    SGB,
    SV_type,
    SV_feature,
    
    N,
    N_CTRL,
    N_CRC,
    
    N_informative_studies,
    
    Eligible,
    Eligibility_reason
  )
]


# ============================================================
# 22. Primary global summary
# ============================================================

summary_by_type <- all_res[
  ,
  .(
    N_universe = .N,
    
    N_eligible =
      sum(
        Eligible,
        na.rm = TRUE
      ),
    
    N_primary_testable =
      sum(
        Primary_testable,
        na.rm = TRUE
      ),
    
    N_strict_OK =
      sum(
        Strict_OK,
        na.rm = TRUE
      ),
    
    N_singular =
      sum(
        Singular == TRUE,
        na.rm = TRUE
      ),
    
    N_P005 =
      sum(
        Primary_testable &
          P < 0.05,
        na.rm = TRUE
      ),
    
    N_P001 =
      sum(
        Primary_testable &
          P < 0.01,
        na.rm = TRUE
      ),
    
    N_Q010 =
      sum(
        Q_type < 0.10,
        na.rm = TRUE
      ),
    
    N_Q005 =
      sum(
        Q_type < 0.05,
        na.rm = TRUE
      ),
    
    N_candidate =
      sum(
        Candidate,
        na.rm = TRUE
      ),
    
    N_candidate_strict_OK =
      sum(
        Candidate_strict_OK,
        na.rm = TRUE
      ),
    
    Min_P =
      if (
        any(
          Primary_testable
        )
      )
        min(
          P[
            Primary_testable
          ],
          na.rm = TRUE
        )
    else
      NA_real_,
    
    Min_Q_type =
      if (
        any(
          !is.na(
            Q_type
          )
        )
      )
        min(
          Q_type,
          na.rm = TRUE
        )
    else
      NA_real_
  ),
  
  by = SV_type
]


summary_all <- all_res[
  ,
  .(
    SV_type = "ALL",
    
    N_universe = .N,
    
    N_eligible =
      sum(
        Eligible,
        na.rm = TRUE
      ),
    
    N_primary_testable =
      sum(
        Primary_testable,
        na.rm = TRUE
      ),
    
    N_strict_OK =
      sum(
        Strict_OK,
        na.rm = TRUE
      ),
    
    N_singular =
      sum(
        Singular == TRUE,
        na.rm = TRUE
      ),
    
    N_P005 =
      sum(
        Primary_testable &
          P < 0.05,
        na.rm = TRUE
      ),
    
    N_P001 =
      sum(
        Primary_testable &
          P < 0.01,
        na.rm = TRUE
      ),
    
    N_Q010 =
      sum(
        Q_type < 0.10,
        na.rm = TRUE
      ),
    
    N_Q005 =
      sum(
        Q_type < 0.05,
        na.rm = TRUE
      ),
    
    N_candidate =
      sum(
        Candidate,
        na.rm = TRUE
      ),
    
    N_candidate_strict_OK =
      sum(
        Candidate_strict_OK,
        na.rm = TRUE
      ),
    
    Min_P =
      min(
        P[
          Primary_testable
        ],
        na.rm = TRUE
      ),
    
    Min_Q_type =
      min(
        Q_type,
        na.rm = TRUE
      )
  )
]


global_summary <- rbindlist(
  list(
    summary_by_type,
    summary_all
  ),
  fill = TRUE
)


# ============================================================
# 23. Model-status summary
# ============================================================

model_status_summary <- all_res[
  ,
  .N,
  by = .(
    SV_type,
    Model_category,
    Model_status
  )
][
  order(
    SV_type,
    -N
  )
]


# ============================================================
# 24. Candidate model-quality summary
# ============================================================

candidate_quality <- candidates[
  ,
  .(
    N_candidate = .N,
    
    N_OK =
      sum(
        Model_status == "OK",
        na.rm = TRUE
      ),
    
    N_singular =
      sum(
        Singular == TRUE,
        na.rm = TRUE
      ),
    
    N_optimizer_OK =
      sum(
        Optimizer_OK == TRUE,
        na.rm = TRUE
      ),
    
    N_other_warning =
      sum(
        Model_status != "OK" &
          Singular != TRUE,
        na.rm = TRUE
      ),
    
    N_SGB =
      uniqueN(
        SGB
      )
  ),
  
  by = SV_type
]


# ============================================================
# 25. Save complete results
# ============================================================

fwrite(
  d_results,
  file.path(
    TABLE_DIR,
    "04_dSGV_CTRL_vs_CRC_raw_model_results.tsv.gz"
  ),
  sep = "\t"
)


fwrite(
  v_results,
  file.path(
    TABLE_DIR,
    "05_vSGV_CTRL_vs_CRC_raw_model_results.tsv.gz"
  ),
  sep = "\t"
)


fwrite(
  all_res,
  file.path(
    TABLE_DIR,
    "06_obese_CTRL_vs_CRC_all_SV_results.tsv.gz"
  ),
  sep = "\t"
)


fwrite(
  eligibility,
  file.path(
    TABLE_DIR,
    "07_feature_eligibility.tsv"
  ),
  sep = "\t"
)


fwrite(
  all_res[
    Primary_testable &
      P < 0.05
  ][
    order(P)
  ],
  file.path(
    TABLE_DIR,
    "08_nominal_P005_SVs.tsv"
  ),
  sep = "\t"
)


# PRIMARY candidate
fwrite(
  candidates,
  file.path(
    TABLE_DIR,
    "09_candidate_SVs_PRIMARY_no_model_status_filter.tsv"
  ),
  sep = "\t"
)


# STRICT candidate
fwrite(
  candidates_strict,
  file.path(
    TABLE_DIR,
    "10_candidate_SVs_STRICT_model_status_OK.tsv"
  ),
  sep = "\t"
)


fwrite(
  fdr05,
  file.path(
    TABLE_DIR,
    "11_FDR05_SVs_PRIMARY.tsv"
  ),
  sep = "\t"
)


# ============================================================
# 26. Save QC summaries
# ============================================================

fwrite(
  global_summary,
  file.path(
    QC_DIR,
    "02_global_SV_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  model_status_summary,
  file.path(
    QC_DIR,
    "03_model_status_summary.tsv"
  ),
  sep = "\t"
)


fwrite(
  candidate_quality,
  file.path(
    QC_DIR,
    "04_PRIMARY_candidate_model_quality.tsv"
  ),
  sep = "\t"
)


capture.output(
  sessionInfo(),
  file = file.path(
    QC_DIR,
    "05_sessionInfo_CTRL_CRC.txt"
  )
)


# ============================================================
# 27. Final report
# ============================================================

cat(
  "\n\n========================================\n"
)

cat(
  "OBESE CTRL vs CRC SV ANALYSIS COMPLETE\n"
)

cat(
  "========================================\n"
)


print(
  global_summary
)


cat(
  "\nPRIMARY candidate SVs",
  "(no Model_status filter):",
  nrow(candidates),
  "\n"
)


cat(
  "STRICT candidate SVs",
  "(Model_status == OK):",
  nrow(candidates_strict),
  "\n"
)


cat(
  "\nCandidate model-quality summary:\n"
)

print(
  candidate_quality
)