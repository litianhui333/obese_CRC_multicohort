# ==============================================================================
# 01_crossBMI_CTRL_vs_CRC.R
#
# Purpose:
#   Evaluate obese-discovered SV loci in a common all-BMI feature space.
#
# Analysis:
#   165 unique all-BMI SV features
#
#   Normal:
#     CTRL vs CRC
#
#   Ow:
#     CTRL vs CRC
#
#   Ob:
#     CTRL vs CRC
#
# Models:
#
# dSGV:
#   Deletion = 1 - dSGV
#
#   glmer(
#     Deletion ~ Disease + Age + Sex + (1 | Study),
#     family = binomial()
#   )
#
# vSGV:
#   lmer(
#     Value ~ Disease + Age + Sex + (1 | Study)
#   )
#
# Primary:
#   retain every eligible model returning a finite P value.
#
# Strict QC:
#   Model_status == "OK"
#
# NOTE:
#   This is targeted cross-BMI evaluation rather than de novo discovery.
# ==============================================================================


suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
})


# ==============================================================================
# 1. Load inputs
# ==============================================================================

source(
  paste0(
    "H:/My Drive/proj_ob_AA_CRC/analyses/",
    "2026-09-19_SV-crossBMI-analysis/",
    "00_source_crossBMI_inputs.R"
  )
)


# ==============================================================================
# 2. Parameters
# ==============================================================================

BMI_GROUPS <- c(
  "Normal",
  "Ow",
  "Ob"
)

MIN_N_CTRL <- 20
MIN_N_CRC <- 20
MIN_INFORMATIVE_STUDIES <- 3

MIN_DSGV_MINOR_STATE <- 5


params <- data.table(
  Parameter = c(
    "BMI_groups",
    "Disease_comparison",
    "N_joint_SVs",
    "N_dSGV",
    "N_vSGV",
    "Minimum_CTRL_N",
    "Minimum_CRC_N",
    "Minimum_informative_studies",
    "Minimum_dSGV_minor_state",
    "Primary_policy",
    "Strict_policy"
  ),
  
  Value = c(
    paste(
      BMI_GROUPS,
      collapse = ","
    ),
    "CTRL_vs_CRC",
    165,
    124,
    41,
    MIN_N_CTRL,
    MIN_N_CRC,
    MIN_INFORMATIVE_STUDIES,
    MIN_DSGV_MINOR_STATE,
    "All eligible models with finite P",
    "Model_status == OK"
  )
)

fwrite(
  params,
  file.path(
    ROOT,
    "00_crossBMI_analysis_parameters.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 3. Fit diagnostics
# ==============================================================================

get_fit_diagnostics <- function(fit) {
  
  singular <- tryCatch(
    lme4::isSingular(
      fit,
      tol = 1e-4
    ),
    error = function(e) NA
  )
  
  
  opt_code <- tryCatch(
    fit@optinfo$conv$opt,
    error = function(e) NULL
  )
  
  
  optimizer_ok <- TRUE
  
  
  if (!is.null(opt_code)) {
    
    z <- suppressWarnings(
      as.numeric(
        unlist(opt_code)
      )
    )
    
    z <- z[
      is.finite(z)
    ]
    
    if (length(z) > 0) {
      
      optimizer_ok <- all(
        z == 0
      )
    }
  }
  
  
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
  
  
  if (
    optimizer_ok &&
    is.na(model_message)
  ) {
    
    model_status <- "OK"
    
  } else if (
    optimizer_ok
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


# ==============================================================================
# 4. Prepare one feature × BMI dataset
# ==============================================================================

make_feature_data <- function(
    matrix_dt,
    feature,
    bmi_group
) {
  
  x <- data.table(
    SGV_ID =
      matrix_dt$SGV_ID,
    
    Value =
      suppressWarnings(
        as.numeric(
          matrix_dt[[feature]]
        )
      )
  )
  
  
  m <- meta_all[
    BMI_model == bmi_group &
      Disease_model %in%
      c(
        "CTRL",
        "CRC"
      )
  ]
  
  
  dat <- merge(
    m,
    x,
    by = "SGV_ID",
    all = FALSE,
    sort = FALSE
  )
  
  
  dat <- dat[
    !is.na(Value) &
      !is.na(Age_model) &
      !is.na(Sex_model) &
      Sex_model != "" &
      !is.na(Study_model) &
      Study_model != ""
  ]
  
  
  dat[
    ,
    Disease_factor :=
      factor(
        Disease_model,
        levels = c(
          "CTRL",
          "CRC"
        )
      )
  ]
  
  
  dat[
    ,
    Sex_factor :=
      droplevels(
        factor(
          Sex_model
        )
      )
  ]
  
  
  dat[
    ,
    Study_factor :=
      droplevels(
        factor(
          Study_model
        )
      )
  ]
  
  
  dat
}


# ==============================================================================
# 5. Eligibility
# ==============================================================================

check_eligibility <- function(
    dat,
    sv_type
) {
  
  n_ctrl <- sum(
    dat$Disease_model == "CTRL"
  )
  
  n_crc <- sum(
    dat$Disease_model == "CRC"
  )
  
  
  informative_studies <- 0L
  
  
  if (nrow(dat) > 0) {
    
    tab <- table(
      dat$Study_factor,
      dat$Disease_factor
    )
    
    
    if (
      all(
        c(
          "CTRL",
          "CRC"
        ) %in%
        colnames(tab)
      )
    ) {
      
      informative_studies <- sum(
        tab[, "CTRL"] > 0 &
          tab[, "CRC"] > 0
      )
    }
  }
  
  
  reasons <- character(0)
  
  
  if (n_ctrl < MIN_N_CTRL) {
    
    reasons <- c(
      reasons,
      "N_CTRL_below_threshold"
    )
  }
  
  
  if (n_crc < MIN_N_CRC) {
    
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
  
  
  minor_state <- NA_integer_
  
  
  if (sv_type == "dSGV") {
    
    values <- unique(
      dat$Value
    )
    
    if (
      !all(
        values %in%
        c(
          0,
          1
        )
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
  
  
  if (sv_type == "vSGV") {
    
    s <- sd(
      dat$Value,
      na.rm = TRUE
    )
    
    
    if (
      !is.finite(s) ||
      s == 0
    ) {
      
      reasons <- c(
        reasons,
        "vSGV_zero_variance"
      )
    }
  }
  
  
  list(
    Eligible =
      length(reasons) == 0,
    
    Eligibility_reason =
      if (
        length(reasons) == 0
      )
        "Eligible"
    else
      paste(
        reasons,
        collapse = ";"
      ),
    
    N =
      nrow(dat),
    
    N_CTRL =
      n_ctrl,
    
    N_CRC =
      n_crc,
    
    N_informative_studies =
      informative_studies,
    
    Minor_state =
      minor_state
  )
}


# ==============================================================================
# 6. Fit one dSGV
# ==============================================================================

fit_dSGV <- function(
    feature,
    bmi_group,
    sgb,
    joint_id
) {
  
  dat <- make_feature_data(
    d_all,
    feature,
    bmi_group
  )
  
  
  eli <- check_eligibility(
    dat,
    "dSGV"
  )
  
  
  base <- data.table(
    Joint_SV_ID =
      joint_id,
    
    SGB =
      sgb,
    
    SV_type =
      "dSGV",
    
    SV_feature =
      feature,
    
    BMI_Group =
      bmi_group,
    
    N =
      eli$N,
    
    N_CTRL =
      eli$N_CTRL,
    
    N_CRC =
      eli$N_CRC,
    
    N_informative_studies =
      eli$N_informative_studies,
    
    Minor_state =
      eli$Minor_state,
    
    Eligible =
      eli$Eligible,
    
    Eligibility_reason =
      eli$Eligibility_reason
  )
  
  
  if (!eli$Eligible) {
    
    return(
      cbind(
        base,
        data.table(
          CTRL_value =
            NA_real_,
          
          CRC_value =
            NA_real_,
          
          CTRL_deletion_rate =
            NA_real_,
          
          CRC_deletion_rate =
            NA_real_,
          
          Beta =
            NA_real_,
          
          SE =
            NA_real_,
          
          OR =
            NA_real_,
          
          P =
            NA_real_,
          
          Singular =
            NA,
          
          Optimizer_OK =
            NA,
          
          Model_message =
            NA_character_,
          
          Model_status =
            "Not_eligible"
        )
      )
    )
  }
  
  
  dat[
    ,
    Deletion := 1 - Value
  ]
  
  
  ctrl_value <- mean(
    dat[
      Disease_model == "CTRL",
      Value
    ]
  )
  
  crc_value <- mean(
    dat[
      Disease_model == "CRC",
      Value
    ]
  )
  
  ctrl_del <- mean(
    dat[
      Disease_model == "CTRL",
      Deletion
    ]
  )
  
  crc_del <- mean(
    dat[
      Disease_model == "CRC",
      Deletion
    ]
  )
  
  
  fit <- tryCatch(
    
    glmer(
      Deletion ~
        Disease_factor +
        Age_model +
        Sex_factor +
        (1 | Study_factor),
      
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
        base,
        data.table(
          CTRL_value =
            ctrl_value,
          
          CRC_value =
            crc_value,
          
          CTRL_deletion_rate =
            ctrl_del,
          
          CRC_deletion_rate =
            crc_del,
          
          Beta =
            NA_real_,
          
          SE =
            NA_real_,
          
          OR =
            NA_real_,
          
          P =
            NA_real_,
          
          Singular =
            NA,
          
          Optimizer_OK =
            FALSE,
          
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
    "^Disease_factorCRC$",
    rownames(cf)
  )
  
  
  if (length(disease_row) != 1) {
    
    return(
      cbind(
        base,
        data.table(
          CTRL_value =
            ctrl_value,
          
          CRC_value =
            crc_value,
          
          CTRL_deletion_rate =
            ctrl_del,
          
          CRC_deletion_rate =
            crc_del,
          
          Beta =
            NA_real_,
          
          SE =
            NA_real_,
          
          OR =
            NA_real_,
          
          P =
            NA_real_,
          
          Singular =
            NA,
          
          Optimizer_OK =
            NA,
          
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
  
  
  diag <- get_fit_diagnostics(
    fit
  )
  
  
  cbind(
    base,
    data.table(
      CTRL_value =
        ctrl_value,
      
      CRC_value =
        crc_value,
      
      CTRL_deletion_rate =
        ctrl_del,
      
      CRC_deletion_rate =
        crc_del,
      
      Beta =
        beta,
      
      SE =
        se,
      
      OR =
        exp(beta),
      
      P =
        p,
      
      Singular =
        diag$Singular,
      
      Optimizer_OK =
        diag$Optimizer_OK,
      
      Model_message =
        diag$Model_message,
      
      Model_status =
        diag$Model_status
    )
  )
}


# ==============================================================================
# 7. Fit one vSGV
# ==============================================================================

fit_vSGV <- function(
    feature,
    bmi_group,
    sgb,
    joint_id
) {
  
  dat <- make_feature_data(
    v_all,
    feature,
    bmi_group
  )
  
  
  eli <- check_eligibility(
    dat,
    "vSGV"
  )
  
  
  base <- data.table(
    Joint_SV_ID =
      joint_id,
    
    SGB =
      sgb,
    
    SV_type =
      "vSGV",
    
    SV_feature =
      feature,
    
    BMI_Group =
      bmi_group,
    
    N =
      eli$N,
    
    N_CTRL =
      eli$N_CTRL,
    
    N_CRC =
      eli$N_CRC,
    
    N_informative_studies =
      eli$N_informative_studies,
    
    Minor_state =
      NA_integer_,
    
    Eligible =
      eli$Eligible,
    
    Eligibility_reason =
      eli$Eligibility_reason
  )
  
  
  if (!eli$Eligible) {
    
    return(
      cbind(
        base,
        data.table(
          CTRL_value =
            NA_real_,
          
          CRC_value =
            NA_real_,
          
          CTRL_deletion_rate =
            NA_real_,
          
          CRC_deletion_rate =
            NA_real_,
          
          Beta =
            NA_real_,
          
          SE =
            NA_real_,
          
          OR =
            NA_real_,
          
          P =
            NA_real_,
          
          Singular =
            NA,
          
          Optimizer_OK =
            NA,
          
          Model_message =
            NA_character_,
          
          Model_status =
            "Not_eligible"
        )
      )
    )
  }
  
  
  ctrl_value <- mean(
    dat[
      Disease_model == "CTRL",
      Value
    ],
    na.rm = TRUE
  )
  
  crc_value <- mean(
    dat[
      Disease_model == "CRC",
      Value
    ],
    na.rm = TRUE
  )
  
  
  fit <- tryCatch(
    
    lmerTest::lmer(
      Value ~
        Disease_factor +
        Age_model +
        Sex_factor +
        (1 | Study_factor),
      
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
        base,
        data.table(
          CTRL_value =
            ctrl_value,
          
          CRC_value =
            crc_value,
          
          CTRL_deletion_rate =
            NA_real_,
          
          CRC_deletion_rate =
            NA_real_,
          
          Beta =
            NA_real_,
          
          SE =
            NA_real_,
          
          OR =
            NA_real_,
          
          P =
            NA_real_,
          
          Singular =
            NA,
          
          Optimizer_OK =
            FALSE,
          
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
    "^Disease_factorCRC$",
    rownames(cf)
  )
  
  
  if (length(disease_row) != 1) {
    
    return(
      cbind(
        base,
        data.table(
          CTRL_value =
            ctrl_value,
          
          CRC_value =
            crc_value,
          
          CTRL_deletion_rate =
            NA_real_,
          
          CRC_deletion_rate =
            NA_real_,
          
          Beta =
            NA_real_,
          
          SE =
            NA_real_,
          
          OR =
            NA_real_,
          
          P =
            NA_real_,
          
          Singular =
            NA,
          
          Optimizer_OK =
            NA,
          
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
  
  
  if (length(p_col) >= 1) {
    
    p <- cf[
      disease_row,
      p_col[1]
    ]
    
  } else {
    
    p <- NA_real_
  }
  
  
  diag <- get_fit_diagnostics(
    fit
  )
  
  
  cbind(
    base,
    data.table(
      CTRL_value =
        ctrl_value,
      
      CRC_value =
        crc_value,
      
      CTRL_deletion_rate =
        NA_real_,
      
      CRC_deletion_rate =
        NA_real_,
      
      Beta =
        beta,
      
      SE =
        se,
      
      OR =
        NA_real_,
      
      P =
        p,
      
      Singular =
        diag$Singular,
      
      Optimizer_OK =
        diag$Optimizer_OK,
      
      Model_message =
        diag$Model_message,
      
      Model_status =
        diag$Model_status
    )
  )
}


# ==============================================================================
# 8. Run all models
# ==============================================================================

results <- vector(
  "list",
  nrow(joint_features) *
    length(BMI_GROUPS)
)


counter <- 0L


for (
  i in seq_len(
    nrow(joint_features)
  )
) {
  
  info <- joint_features[i]
  
  
  feature <- info$AllBMI_SV_feature
  sv_type <- info$SV_type
  sgb <- info$SGB
  joint_id <- info$Joint_SV_ID
  
  
  for (bmi in BMI_GROUPS) {
    
    counter <- counter + 1L
    
    
    if (sv_type == "dSGV") {
      
      res <- fit_dSGV(
        feature = feature,
        bmi_group = bmi,
        sgb = sgb,
        joint_id = joint_id
      )
      
    } else if (
      sv_type == "vSGV"
    ) {
      
      res <- fit_vSGV(
        feature = feature,
        bmi_group = bmi,
        sgb = sgb,
        joint_id = joint_id
      )
      
    } else {
      
      stop(
        "Unknown SV type: ",
        sv_type
      )
    }
    
    
    results[[counter]] <- res
  }
  
  
  if (
    i %% 10 == 0 ||
    i == nrow(joint_features)
  ) {
    
    cat(
      "Features completed:",
      i,
      "/",
      nrow(joint_features),
      "\n"
    )
  }
  
  
  if (
    i %% 50 == 0
  ) {
    
    tmp <- rbindlist(
      results[
        seq_len(counter)
      ],
      fill = TRUE
    )
    
    fwrite(
      tmp,
      file.path(
        TABLE_DIR,
        "checkpoint_crossBMI_models.tsv.gz"
      ),
      sep = "\t"
    )
  }
}


all_res <- rbindlist(
  results[
    seq_len(counter)
  ],
  fill = TRUE
)


# ==============================================================================
# 9. Model-status flags
# ==============================================================================

all_res[
  ,
  Primary_testable :=
    Eligible == TRUE &
    is.finite(P)
]


all_res[
  ,
  Strict_OK :=
    Primary_testable &
    Model_status == "OK"
]


all_res[
  ,
  Singular_only :=
    Primary_testable &
    Singular == TRUE &
    Optimizer_OK == TRUE
]


# ==============================================================================
# 10. Targeted BH-FDR
#
# Not used to define obesity specificity.
# Saved only as descriptive/statistical context.
# ==============================================================================

all_res[
  ,
  `:=`(
    Q_BMI_type =
      NA_real_,
    
    Q_BMI_all =
      NA_real_
  )
]


for (bmi in BMI_GROUPS) {
  
  for (
    typ in c(
      "dSGV",
      "vSGV"
    )
  ) {
    
    idx <- which(
      all_res$BMI_Group == bmi &
        all_res$SV_type == typ &
        all_res$Primary_testable
    )
    
    
    all_res$Q_BMI_type[idx] <-
      p.adjust(
        all_res$P[idx],
        method = "BH"
      )
  }
  
  
  idx <- which(
    all_res$BMI_Group == bmi &
      all_res$Primary_testable
  )
  
  
  all_res$Q_BMI_all[idx] <-
    p.adjust(
      all_res$P[idx],
      method = "BH"
    )
}


# ==============================================================================
# 11. Save long-format results
# ==============================================================================

fwrite(
  all_res,
  file.path(
    TABLE_DIR,
    "01_crossBMI_CTRL_vs_CRC_long.tsv.gz"
  ),
  sep = "\t"
)


# ==============================================================================
# 12. Model-QC summary
# ==============================================================================

model_qc <- all_res[
  ,
  .(
    N_features =
      .N,
    
    N_eligible =
      sum(
        Eligible,
        na.rm = TRUE
      ),
    
    N_testable =
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
      )
  ),
  
  by = .(
    BMI_Group,
    SV_type
  )
]


fwrite(
  model_qc,
  file.path(
    QC_DIR,
    "01_crossBMI_model_QC.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 13. Create wide table
# ==============================================================================

wide_source <- all_res[
  ,
  .(
    Joint_SV_ID,
    SGB,
    SV_type,
    SV_feature,
    BMI_Group,
    N,
    N_CTRL,
    N_CRC,
    N_informative_studies,
    Beta,
    SE,
    OR,
    P,
    Q_BMI_type,
    Singular,
    Optimizer_OK,
    Model_status,
    Primary_testable
  )
]


wide <- dcast(
  wide_source,
  Joint_SV_ID +
    SGB +
    SV_type +
    SV_feature ~
    BMI_Group,
  
  value.var = c(
    "N",
    "N_CTRL",
    "N_CRC",
    "N_informative_studies",
    "Beta",
    "SE",
    "OR",
    "P",
    "Q_BMI_type",
    "Singular",
    "Optimizer_OK",
    "Model_status",
    "Primary_testable"
  )
)


# ==============================================================================
# 14. Rename Normal → Nw for compact output
# ==============================================================================

setnames(
  wide,
  old = grep(
    "_Normal$",
    names(wide),
    value = TRUE
  ),
  new = sub(
    "_Normal$",
    "_Nw",
    grep(
      "_Normal$",
      names(wide),
      value = TRUE
    )
  )
)


# ==============================================================================
# 15. Descriptive cross-BMI pattern flags
#
# IMPORTANT:
# These are descriptive patterns, NOT formal evidence of interaction.
# ==============================================================================

wide[
  ,
  Ob_only_nominal :=
    !is.na(P_Ob) &
    !is.na(P_Nw) &
    !is.na(P_Ow) &
    P_Ob < 0.05 &
    P_Nw >= 0.05 &
    P_Ow >= 0.05
]


wide[
  ,
  Ob_stronger_absolute_effect :=
    !is.na(Beta_Ob) &
    !is.na(Beta_Nw) &
    !is.na(Beta_Ow) &
    abs(Beta_Ob) >
    abs(Beta_Nw) &
    abs(Beta_Ob) >
    abs(Beta_Ow)
]


wide[
  ,
  Ob_enriched_descriptive :=
    Ob_only_nominal &
    Ob_stronger_absolute_effect
]


wide[
  ,
  Same_direction_all_BMI :=
    !is.na(Beta_Ob) &
    !is.na(Beta_Nw) &
    !is.na(Beta_Ow) &
    sign(Beta_Ob) ==
    sign(Beta_Nw) &
    sign(Beta_Ob) ==
    sign(Beta_Ow)
]


wide[
  ,
  Max_nonOb_abs_Beta :=
    pmax(
      abs(Beta_Nw),
      abs(Beta_Ow),
      na.rm = TRUE
    )
]


wide[
  ,
  Ob_to_nonOb_effect_ratio :=
    fifelse(
      is.finite(Max_nonOb_abs_Beta) &
        Max_nonOb_abs_Beta > 0,
      abs(Beta_Ob) /
        Max_nonOb_abs_Beta,
      NA_real_
    )
]


# ==============================================================================
# 16. Add obese-candidate mapping information
# ==============================================================================

map_agg <- cross_map[
  ,
  .(
    Candidate_IDs =
      paste(
        sort(
          unique(
            Candidate_ID
          )
        ),
        collapse = ";"
      ),
    
    N_obese_candidates =
      uniqueN(
        Candidate_ID
      ),
    
    Obese_SV_features =
      paste(
        sort(
          unique(
            Obese_SV_feature
          )
        ),
        collapse = ";"
      ),
    
    Match_types =
      paste(
        sort(
          unique(
            Match_type
          )
        ),
        collapse = ";"
      ),
    
    Mapping_topologies =
      paste(
        sort(
          unique(
            Mapping_topology
          )
        ),
        collapse = ";"
      ),
    
    Min_Jaccard =
      min(
        Jaccard,
        na.rm = TRUE
      ),
    
    Max_Jaccard =
      max(
        Jaccard,
        na.rm = TRUE
      )
  ),
  
  by = .(
    SGB,
    SV_type,
    AllBMI_SV_feature
  )
]


wide <- merge(
  wide,
  map_agg,
  by.x = c(
    "SGB",
    "SV_type",
    "SV_feature"
  ),
  by.y = c(
    "SGB",
    "SV_type",
    "AllBMI_SV_feature"
  ),
  all.x = TRUE,
  sort = FALSE
)


# ==============================================================================
# 17. Order table
# ==============================================================================

setorder(
  wide,
  P_Ob
)


fwrite(
  wide,
  file.path(
    TABLE_DIR,
    "02_crossBMI_CTRL_vs_CRC_wide.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 18. Descriptive obesity-enriched subset
# ==============================================================================

ob_enriched <- wide[
  Ob_enriched_descriptive == TRUE
][
  order(
    P_Ob
  )
]


fwrite(
  ob_enriched,
  file.path(
    TABLE_DIR,
    "03_descriptive_obesity_enriched_SVs.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 19. BMI-shared directional subset
# ==============================================================================

shared_direction <- wide[
  Same_direction_all_BMI == TRUE
][
  order(
    P_Ob
  )
]


fwrite(
  shared_direction,
  file.path(
    TABLE_DIR,
    "04_same_direction_across_BMI_SVs.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 20. Pattern summary
# ==============================================================================

pattern_summary <- data.table(
  Metric = c(
    "Joint_SVs",
    "Ob_P005",
    "Nw_P005",
    "Ow_P005",
    "Ob_only_nominal",
    "Ob_stronger_absolute_effect",
    "Ob_enriched_descriptive",
    "Same_direction_all_BMI"
  ),
  
  N = c(
    nrow(wide),
    
    sum(
      wide$P_Ob < 0.05,
      na.rm = TRUE
    ),
    
    sum(
      wide$P_Nw < 0.05,
      na.rm = TRUE
    ),
    
    sum(
      wide$P_Ow < 0.05,
      na.rm = TRUE
    ),
    
    sum(
      wide$Ob_only_nominal,
      na.rm = TRUE
    ),
    
    sum(
      wide$Ob_stronger_absolute_effect,
      na.rm = TRUE
    ),
    
    sum(
      wide$Ob_enriched_descriptive,
      na.rm = TRUE
    ),
    
    sum(
      wide$Same_direction_all_BMI,
      na.rm = TRUE
    )
  )
)


fwrite(
  pattern_summary,
  file.path(
    QC_DIR,
    "02_crossBMI_pattern_summary.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 21. Session info
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    QC_DIR,
    "03_sessionInfo.txt"
  )
)


# ==============================================================================
# 22. Final report
# ==============================================================================

cat(
  "\n\n========================================\n"
)

cat(
  "CROSS-BMI SV ANALYSIS COMPLETE\n"
)

cat(
  "========================================\n"
)


cat(
  "\nModel QC:\n"
)

print(
  model_qc
)


cat(
  "\nPattern summary:\n"
)

print(
  pattern_summary
)


cat(
  "\nTop descriptive obesity-enriched SVs:\n"
)


print(
  ob_enriched[
    1:min(
      20,
      .N
    ),
    .(
      Joint_SV_ID,
      SGB,
      SV_type,
      SV_feature,
      
      Beta_Nw,
      P_Nw,
      
      Beta_Ow,
      P_Ow,
      
      Beta_Ob,
      P_Ob,
      
      Ob_to_nonOb_effect_ratio,
      
      Candidate_IDs,
      Match_types
    )
  ]
)


library(data.table)

ROOT <- paste0(
  "H:/My Drive/proj_ob_AA_CRC/analyses/",
  "2026-09-19_SV-crossBMI-analysis"
)

TABLE_DIR <- file.path(
  ROOT,
  "tables"
)

RES <- fread(
  file.path(
    TABLE_DIR,
    "02_crossBMI_CTRL_vs_CRC_wide.tsv"
  )
)

# ============================================================
# Strict cross-BMI pattern
# Ob P < 0.01
# Nw/Ow P >= 0.05
# ============================================================

strict_P001 <- RES[
  !is.na(P_Ob) &
    !is.na(P_Ow) &
    !is.na(P_Nw) &
    
    P_Ob < 0.01 &
    P_Ow >= 0.05 &
    P_Nw >= 0.05
]


# ============================================================
# Relaxed cross-BMI pattern
# Ob P < 0.05
# Nw/Ow P >= 0.05
# ============================================================

relaxed_P005 <- RES[
  !is.na(P_Ob) &
    !is.na(P_Ow) &
    !is.na(P_Nw) &
    
    P_Ob < 0.05 &
    P_Ow >= 0.05 &
    P_Nw >= 0.05
]


cat(
  "Strict P<0.01:",
  nrow(strict_P001),
  "\n"
)

cat(
  "Relaxed P<0.05:",
  nrow(relaxed_P005),
  "\n"
)

cat(
  "Additional relaxed:",
  nrow(relaxed_P005) -
    nrow(strict_P001),
  "\n"
)