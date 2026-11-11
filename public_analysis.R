# ==============================================================================
# PUBLIC REPRODUCIBILITY ANALYSIS
# Pretreatment CBC Mahalanobis distance and survival after definitive CCRT
# ==============================================================================
#
# This script reproduces the principal statistical analyses reported in the
# manuscript. It does not create manuscript tables, figures, HTML files, or RDS
# model objects. Results are written as plain CSV files to results/public_analysis.
#
# Required input files and the folder layout are described in README.md.
# ==============================================================================

# 1. Setup ---------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg)) {
  ROOT <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  ROOT <- getwd()
}

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("The survival package is required.")
}

suppressPackageStartupMessages(library(survival))

DATA_DIR <- file.path(ROOT, "data")
OUTPUT_DIR <- file.path(ROOT, "results", "public_analysis")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

N_BOOTSTRAP <- 1000L
RANDOM_SEED <- 20261106L
REPORTED_BASELINE_CUTOFF <- 3.60157833984888  # Reported as 3.602
CBC_FEATURES <- c("Hb", "MCV", "RDW", "ANC", "ALC", "AMC", "PLT", "MPV")
ENDPOINTS <- c("OS", "PFS")
TIMEPOINTS <- c("baseline", "week2", "week3", "week4")

write_result <- function(x, filename) {
  write.csv(x, file.path(OUTPUT_DIR, filename), row.names = FALSE, na = "")
}

# 2. Read and validate input data ----------------------------------------------
patient_file <- file.path(DATA_DIR, "analysis_dataset.csv")
cbc_file <- file.path(DATA_DIR, "selected_CBC_long.csv")
mean_file <- file.path(DATA_DIR, "reference_mean.csv")
covariance_file <- file.path(DATA_DIR, "reference_covariance.csv")

required_files <- c(patient_file, cbc_file, mean_file, covariance_file)
if (any(!file.exists(required_files))) {
  stop("Missing input file(s): ", paste(required_files[!file.exists(required_files)], collapse = ", "))
}

patients <- read.csv(patient_file, check.names = FALSE)
cbc_long <- read.csv(cbc_file, check.names = FALSE)
reference_mean <- read.csv(mean_file, check.names = FALSE)
reference_covariance <- read.csv(covariance_file, check.names = FALSE)

required_patient_columns <- c(
  "study_id", "age", "stage_group", "histology_group",
  "EQD2_Apoint_Gy", "treatment_duration_days",
  "os_months", "os_event", "pfs_months", "pfs_event",
  paste0("baseline_", CBC_FEATURES)
)

required_cbc_columns <- c(
  "study_id", "timepoint", "selection_rule",
  CBC_FEATURES
)

missing_patient_columns <- setdiff(required_patient_columns, names(patients))
missing_cbc_columns <- setdiff(required_cbc_columns, names(cbc_long))

if (length(missing_patient_columns)) {
  stop("Missing patient-level columns: ", paste(missing_patient_columns, collapse = ", "))
}
if (length(missing_cbc_columns)) {
  stop("Missing longitudinal CBC columns: ", paste(missing_cbc_columns, collapse = ", "))
}
if (anyDuplicated(patients$study_id)) {
  stop("study_id must be unique in analysis_dataset.csv")
}
if (!all(patients$os_event %in% 0:1) || !all(patients$pfs_event %in% 0:1)) {
  stop("os_event and pfs_event must be coded as 0 or 1.")
}

if (!identical(reference_mean$feature, CBC_FEATURES)) {
  stop("reference_mean.csv must contain the eight CBC features in the required order.")
}
if (!identical(reference_covariance$feature, CBC_FEATURES)) {
  stop("reference_covariance.csv must contain the eight CBC features in the required order.")
}

reference_covariance_matrix <- as.matrix(reference_covariance[, CBC_FEATURES])
if (min(eigen(reference_covariance_matrix, symmetric = TRUE)$values) <= 0) {
  stop("The reference covariance matrix must be positive definite.")
}

# 3. Calculate Mahalanobis distance --------------------------------------------
calculate_D <- function(x) {
  sqrt(mahalanobis(
    as.matrix(x[, CBC_FEATURES, drop = FALSE]),
    center = reference_mean$mean,
    cov = reference_covariance_matrix
  ))
}

baseline_cbc <- patients[, paste0("baseline_", CBC_FEATURES), drop = FALSE]
names(baseline_cbc) <- CBC_FEATURES
patients$D_base <- calculate_D(baseline_cbc)

if (!"baseline_NLR" %in% names(patients)) {
  patients$baseline_NLR <- patients$baseline_ANC / patients$baseline_ALC
}

cbc_primary <- cbc_long[
  cbc_long$selection_rule == "primary_latest" &
    cbc_long$timepoint %in% TIMEPOINTS,
]
cbc_primary$D <- calculate_D(cbc_primary)

for (timepoint in c("week2", "week3", "week4")) {
  current <- cbc_primary[cbc_primary$timepoint == timepoint, c("study_id", "D")]
  if (anyDuplicated(current$study_id)) {
    stop("More than one primary CBC was found for ", timepoint, ".")
  }
  patients[[paste0("D_", timepoint)]] <- current$D[match(patients$study_id, current$study_id)]
}

# Alternative pretreatment CBC used in the prespecified sensitivity analysis:
# the complete CBC with the lowest Hb during days -30 through 0.
minimum_hb_cbc <- cbc_long[
  cbc_long$selection_rule == "minimum_Hb" & cbc_long$timepoint == "baseline",
]
if (nrow(minimum_hb_cbc) != nrow(patients) || anyDuplicated(minimum_hb_cbc$study_id)) {
  stop("Exactly one minimum-Hb baseline CBC is required for every patient.")
}
minimum_hb_cbc$D <- calculate_D(minimum_hb_cbc)
patients$D_base_minHb <- minimum_hb_cbc$D[match(patients$study_id, minimum_hb_cbc$study_id)]

# Clinical variables used in every adjusted survival model.
patients$age10 <- (patients$age - mean(patients$age, na.rm = TRUE)) / 10
patients$dose10 <- (patients$EQD2_Apoint_Gy - mean(patients$EQD2_Apoint_Gy, na.rm = TRUE)) / 10
patients$duration7 <- (patients$treatment_duration_days - mean(patients$treatment_duration_days, na.rm = TRUE)) / 7
patients$stage_high <- as.integer(patients$stage_group == "IIIC2_to_IVB")
patients$non_sqcc <- as.integer(patients$histology_group == "non_sqcc")
patients$D_base_z <- as.numeric(scale(patients$D_base))

clinical_terms <- c("age10", "stage_high", "non_sqcc", "dose10", "duration7")

make_analysis_data <- function(timepoint, endpoint) {
  marker_name <- if (timepoint == "baseline") "D_base" else paste0("D_", timepoint)
  endpoint_lower <- tolower(endpoint)

  keep_columns <- c(
    "study_id", marker_name, if (timepoint != "baseline") "D_base", "D_base_z", clinical_terms,
    paste0(endpoint_lower, "_months"), paste0(endpoint_lower, "_event")
  )
  analysis_data <- patients[, keep_columns]
  analysis_data$D_current <- analysis_data[[marker_name]]
  names(analysis_data)[names(analysis_data) == paste0(endpoint_lower, "_months")] <- "time"
  names(analysis_data)[names(analysis_data) == paste0(endpoint_lower, "_event")] <- "event"

  analysis_data$D_current_z <- as.numeric(scale(analysis_data$D_current))
  complete_columns <- c("time", "event", "D_current_z", "D_base_z", clinical_terms)
  analysis_data <- analysis_data[complete.cases(analysis_data[, complete_columns]), ]
  analysis_data <- analysis_data[analysis_data$time > 0, ]
  analysis_data
}

# 4. Median follow-up using reverse Kaplan-Meier -------------------------------
reverse_km <- survfit(Surv(os_months, 1 - os_event) ~ 1, data = patients)
reverse_km_table <- summary(reverse_km)$table

cohort_summary <- data.frame(
  N = nrow(patients),
  deaths = sum(patients$os_event),
  progression_or_death_events = sum(patients$pfs_event),
  median_follow_up_months = unname(reverse_km_table["median"]),
  follow_up_95CI_low = unname(reverse_km_table["0.95LCL"]),
  follow_up_95CI_high = unname(reverse_km_table["0.95UCL"])
)
write_result(cohort_summary, "cohort_summary.csv")

# 5. Change in D during treatment ---------------------------------------------
# Each weekly value is paired with the same patient's baseline value. The three
# p values are corrected together using the Holm method, as reported in Fig. 1.
d_wide <- patients[, c("study_id", "D_base", "D_week2", "D_week3", "D_week4")]
d_change_summary <- data.frame(
  timepoint = c("Baseline", "Week 2", "Week 3", "Week 4"),
  N = vapply(d_wide[-1], function(x) sum(is.finite(x)), integer(1)),
  median = vapply(d_wide[-1], median, numeric(1), na.rm = TRUE),
  Q1 = vapply(d_wide[-1], quantile, numeric(1), probs = 0.25, na.rm = TRUE),
  Q3 = vapply(d_wide[-1], quantile, numeric(1), probs = 0.75, na.rm = TRUE)
)

d_change_tests <- do.call(rbind, lapply(c("D_week2", "D_week3", "D_week4"), function(marker) {
  complete <- is.finite(d_wide$D_base) & is.finite(d_wide[[marker]])
  test <- wilcox.test(
    d_wide$D_base[complete], d_wide[[marker]][complete],
    paired = TRUE, exact = FALSE
  )
  data.frame(
    comparison = paste("Baseline vs", sub("D_week", "Week ", marker)),
    N_pairs = sum(complete),
    raw_p = test$p.value
  )
}))
d_change_tests$Holm_adjusted_p <- p.adjust(d_change_tests$raw_p, method = "holm")

write_result(d_change_summary, "D_over_time_summary.csv")
write_result(d_change_tests, "D_over_time_paired_tests.csv")

# 6. Continuous Cox models and proportional-hazards tests ---------------------
continuous_results <- list()
model_data <- list()
result_index <- 0L

for (timepoint in TIMEPOINTS) {
  for (endpoint in ENDPOINTS) {
    analysis_data <- make_analysis_data(timepoint, endpoint)
    is_week <- timepoint != "baseline"
    model_terms <- c("D_current_z", if (is_week) "D_base_z", clinical_terms)
    model_formula <- as.formula(paste("Surv(time, event) ~", paste(model_terms, collapse = " + ")))

    fitted_model <- coxph(model_formula, data = analysis_data, ties = "efron", x = TRUE)
    model_summary <- summary(fitted_model)
    confidence_interval <- exp(confint(fitted_model)["D_current_z", ])
    ph_test <- cox.zph(fitted_model, transform = function(time) log(time), terms = FALSE)

    result_index <- result_index + 1L
    continuous_results[[result_index]] <- data.frame(
      marker = if (timepoint == "baseline") "Baseline D" else paste("Week", sub("week", "", timepoint), "D"),
      endpoint = endpoint,
      N = fitted_model$n,
      events = fitted_model$nevent,
      D_units_per_SD = sd(analysis_data$D_current),
      adjusted_HR_per_SD = exp(coef(fitted_model)["D_current_z"]),
      CI_low = confidence_interval[1],
      CI_high = confidence_interval[2],
      Wald_p = model_summary$coefficients["D_current_z", "Pr(>|z|)"],
      PH_p = ph_test$table["D_current_z", "p"]
    )

    model_data[[paste(timepoint, endpoint, sep = "_")]] <- list(
      data = analysis_data,
      full_model = fitted_model
    )
  }
}

continuous_results <- do.call(rbind, continuous_results)
write_result(continuous_results, "continuous_cox_models.csv")

# 7. Time-varying effect of baseline D -----------------------------------------
# The interaction uses log(t/12), so the main D coefficient is the HR at 12 months.
time_varying_results <- list()
time_varying_index <- 0L

for (endpoint in ENDPOINTS) {
  stored_model <- model_data[[paste("baseline", endpoint, sep = "_")]]
  analysis_data <- stored_model$data

  # Any clinical term with evidence of non-proportional hazards is also allowed
  # to vary over time, preventing its violation from being attributed to D.
  static_terms <- c("D_current_z", clinical_terms)
  ph_test <- cox.zph(stored_model$full_model, transform = function(time) log(time), terms = FALSE)
  flagged_terms <- rownames(ph_test$table)[
    rownames(ph_test$table) != "GLOBAL" & ph_test$table[, "p"] < 0.05
  ]
  interaction_terms <- unique(c("D_current_z", intersect(flagged_terms, clinical_terms)))
  model_terms <- c(static_terms, paste0("tt(", interaction_terms, ")"))
  model_formula <- as.formula(paste("Surv(time, event) ~", paste(model_terms, collapse = " + ")))

  fitted_model <- coxph(
    model_formula,
    data = analysis_data,
    ties = "efron",
    tt = function(x, time, ...) x * log(time / 12)
  )

  coefficients <- coef(fitted_model)[c("D_current_z", "tt(D_current_z)")]
  covariance <- vcov(fitted_model)[c("D_current_z", "tt(D_current_z)"), c("D_current_z", "tt(D_current_z)")]
  joint_chisq <- drop(t(coefficients) %*% solve(covariance, coefficients))

  time_varying_index <- time_varying_index + 1L
  time_varying_results[[time_varying_index]] <- data.frame(
    endpoint = endpoint,
    N = fitted_model$n,
    events = fitted_model$nevent,
    HR_at_12_months = exp(coefficients[1]),
    time_interaction_p = summary(fitted_model)$coefficients["tt(D_current_z)", "Pr(>|z|)"],
    joint_time_varying_p = pchisq(joint_chisq, df = 2, lower.tail = FALSE)
  )
}

write_result(do.call(rbind, time_varying_results), "baseline_D_time_varying_models.csv")

# Minimum-Hb baseline sensitivity analysis (Supplementary Table S5).
# HRs use the primary baseline D SD (2.085 in these data), allowing direct
# comparison with the primary baseline analysis.
minimum_hb_results <- list()
primary_D_SD <- sd(patients$D_base)

for (endpoint in ENDPOINTS) {
  endpoint_lower <- tolower(endpoint)
  analysis_data <- patients[, c(
    "D_base_minHb", clinical_terms,
    paste0(endpoint_lower, "_months"), paste0(endpoint_lower, "_event")
  )]
  names(analysis_data)[names(analysis_data) == paste0(endpoint_lower, "_months")] <- "time"
  names(analysis_data)[names(analysis_data) == paste0(endpoint_lower, "_event")] <- "event"
  analysis_data$D_current_z <-
    (analysis_data$D_base_minHb - mean(analysis_data$D_base_minHb)) / primary_D_SD
  analysis_data <- analysis_data[complete.cases(analysis_data), ]
  analysis_data <- analysis_data[analysis_data$time > 0, ]

  model_formula <- as.formula(
    paste("Surv(time, event) ~", paste(c("D_current_z", clinical_terms), collapse = " + "))
  )
  static_model <- coxph(
    model_formula, data = analysis_data, ties = "efron", x = TRUE, y = TRUE, model = TRUE
  )
  static_summary <- summary(static_model)
  static_CI <- exp(confint(static_model)["D_current_z", ])
  ph_test <- cox.zph(static_model, transform = function(time) log(time), terms = FALSE)

  set.seed(RANDOM_SEED + 50L + match(endpoint, ENDPOINTS))
  bootstrap_beta <- rep(NA_real_, N_BOOTSTRAP)
  for (bootstrap_index in seq_len(N_BOOTSTRAP)) {
    sampled_rows <- sample.int(nrow(analysis_data), nrow(analysis_data), replace = TRUE)
    bootstrap_model <- tryCatch(
      suppressWarnings(coxph(model_formula, data = analysis_data[sampled_rows, ], ties = "efron")),
      error = function(e) NULL
    )
    if (!is.null(bootstrap_model) && all(is.finite(coef(bootstrap_model)))) {
      bootstrap_beta[bootstrap_index] <- coef(bootstrap_model)["D_current_z"]
    }
  }
  if (mean(is.finite(bootstrap_beta)) < 0.95) {
    stop("More than 5% of minimum-Hb Cox bootstrap samples failed for ", endpoint, ".")
  }
  bootstrap_CI <- exp(quantile(bootstrap_beta, c(0.025, 0.975), na.rm = TRUE))

  flagged_terms <- rownames(ph_test$table)[
    rownames(ph_test$table) != "GLOBAL" & ph_test$table[, "p"] < 0.05
  ]
  interaction_terms <- unique(c("D_current_z", intersect(flagged_terms, clinical_terms)))
  time_varying_formula <- as.formula(paste(
    "Surv(time, event) ~",
    paste(c("D_current_z", clinical_terms, paste0("tt(", interaction_terms, ")")), collapse = " + ")
  ))
  time_varying_model <- coxph(
    time_varying_formula,
    data = analysis_data,
    ties = "efron",
    tt = function(x, time, ...) x * log(time / 12),
    control = coxph.control(iter.max = 100)
  )
  joint_names <- c("D_current_z", "tt(D_current_z)")
  joint_beta <- coef(time_varying_model)[joint_names]
  joint_covariance <- vcov(time_varying_model)[joint_names, joint_names, drop = FALSE]
  joint_chisq <- drop(t(joint_beta) %*% solve(joint_covariance, joint_beta))

  minimum_hb_results[[endpoint]] <- data.frame(
    endpoint = endpoint,
    N = static_model$n,
    events = static_model$nevent,
    HR_scale_primary_D_SD = primary_D_SD,
    adjusted_HR = exp(coef(static_model)["D_current_z"]),
    CI_low = static_CI[1],
    CI_high = static_CI[2],
    Wald_p = static_summary$coefficients["D_current_z", "Pr(>|z|)"],
    bootstrap_CI_low = bootstrap_CI[1],
    bootstrap_CI_high = bootstrap_CI[2],
    PH_p = ph_test$table["D_current_z", "p"],
    time_interaction_p = summary(time_varying_model)$coefficients[
      "tt(D_current_z)", "Pr(>|z|)"
    ],
    joint_time_varying_p = pchisq(joint_chisq, df = 2, lower.tail = FALSE)
  )
}

write_result(do.call(rbind, minimum_hb_results), "minimum_Hb_sensitivity_models.csv")

# 8. Added value of week-specific D beyond baseline D -------------------------
added_value_results <- list()
added_value_index <- 0L

for (timepoint in c("week2", "week3", "week4")) {
  for (endpoint in ENDPOINTS) {
    stored <- model_data[[paste(timepoint, endpoint, sep = "_")]]
    analysis_data <- stored$data
    full_model <- stored$full_model

    reduced_formula <- as.formula(
      paste("Surv(time, event) ~", paste(c("D_base_z", clinical_terms), collapse = " + "))
    )
    reduced_model <- coxph(reduced_formula, data = analysis_data, ties = "efron")
    likelihood_ratio <- 2 * (full_model$loglik[2] - reduced_model$loglik[2])

    added_value_index <- added_value_index + 1L
    added_value_results[[added_value_index]] <- data.frame(
      marker = paste("Week", sub("week", "", timepoint), "D"),
      endpoint = endpoint,
      N = full_model$n,
      events = full_model$nevent,
      correlation_with_baseline_D = cor(analysis_data$D_current, analysis_data$D_base),
      added_value_LR_p = pchisq(max(0, likelihood_ratio), df = 1, lower.tail = FALSE)
    )
  }
}

write_result(do.call(rbind, added_value_results), "weekly_D_added_value.csv")

# 9. Time-dependent cumulative/dynamic ROC analysis ---------------------------
censoring_survival <- function(time, event) {
  fit <- survfit(Surv(time, 1 - event) ~ 1)
  function(at, left_limit = FALSE) {
    index <- findInterval(at, fit$time)
    if (left_limit) {
      exact <- index > 0
      exact[exact] <- at[exact] == fit$time[index[exact]]
      index[exact] <- index[exact] - 1L
    }
    c(1, fit$surv)[index + 1L]
  }
}

roc_weights <- function(time, event, horizon) {
  censoring_km <- censoring_survival(time, event)
  cases <- event == 1 & time <= horizon
  controls <- time > horizon
  case_weights <- control_weights <- numeric(length(time))

  case_censoring_probability <- censoring_km(time[cases], left_limit = TRUE)
  horizon_censoring_probability <- censoring_km(horizon)

  if (!any(cases) || !any(controls) || any(case_censoring_probability <= 0) || horizon_censoring_probability <= 0) {
    stop("Invalid case/control sample for time-dependent ROC analysis.")
  }

  case_weights[cases] <- 1 / case_censoring_probability
  control_weights[controls] <- 1 / horizon_censoring_probability

  list(
    case = case_weights,
    control = control_weights,
    cases = sum(cases),
    controls = sum(controls),
    censored = sum(!cases & !controls)
  )
}

weighted_roc <- function(marker, weights) {
  marker_values <- sort(unique(marker), decreasing = TRUE)
  grouped_weights <- rowsum(
    cbind(weights$case, weights$control),
    match(marker, marker_values),
    reorder = TRUE
  )

  sensitivity <- c(0, cumsum(grouped_weights[, 1]) / sum(grouped_weights[, 1]))
  false_positive_rate <- c(0, cumsum(grouped_weights[, 2]) / sum(grouped_weights[, 2]))
  cutoffs <- c(Inf, (head(marker_values, -1) + tail(marker_values, -1)) / 2, -Inf)
  specificity <- 1 - false_positive_rate
  youden <- sensitivity + specificity - 1

  candidates <- which(is.finite(cutoffs) & youden == max(youden[is.finite(cutoffs)]))
  best <- candidates[which.max(specificity[candidates])]

  list(
    AUC = sum(diff(false_positive_rate) * (head(sensitivity, -1) + tail(sensitivity, -1)) / 2),
    cutoff = cutoffs[best],
    sensitivity = sensitivity[best],
    specificity = specificity[best]
  )
}

bootstrap_roc <- function(time, event, marker, horizon, seed) {
  set.seed(seed)
  auc <- cutoff <- rep(NA_real_, N_BOOTSTRAP)

  for (bootstrap_index in seq_len(N_BOOTSTRAP)) {
    sample_index <- sample.int(length(time), length(time), replace = TRUE)
    estimate <- tryCatch({
      weights <- roc_weights(time[sample_index], event[sample_index], horizon)
      weighted_roc(marker[sample_index], weights)
    }, error = function(e) NULL)

    if (!is.null(estimate)) {
      auc[bootstrap_index] <- estimate$AUC
      cutoff[bootstrap_index] <- estimate$cutoff
    }
  }

  if (mean(is.finite(auc)) < 0.95) {
    stop("More than 5% of ROC bootstrap samples failed.")
  }
  list(AUC = auc, cutoff = cutoff)
}

# Baseline D AUC at 12-60 months.
auc_horizon_results <- list()
auc_bootstrap_36 <- NULL

for (horizon in c(12, 24, 36, 48, 60)) {
  weights <- roc_weights(patients$os_months, patients$os_event, horizon)
  estimate <- weighted_roc(patients$D_base, weights)
  # The 36-month bootstrap is the primary analysis and uses RANDOM_SEED.
  # Other horizons use the same seeds as the manuscript analysis.
  bootstrap_seed <- if (horizon == 36) RANDOM_SEED else RANDOM_SEED + 500L + horizon
  bootstrap <- bootstrap_roc(
    patients$os_months, patients$os_event, patients$D_base,
    horizon, bootstrap_seed
  )

  if (horizon == 36) {
    auc_bootstrap_36 <- bootstrap
  }

  auc_horizon_results[[as.character(horizon)]] <- data.frame(
    OS_horizon_months = horizon,
    N = nrow(patients),
    deaths_by_horizon = weights$cases,
    surviving_beyond_horizon = weights$controls,
    censored_before_horizon = weights$censored,
    AUC = estimate$AUC,
    CI_low = quantile(bootstrap$AUC, 0.025, na.rm = TRUE),
    CI_high = quantile(bootstrap$AUC, 0.975, na.rm = TRUE)
  )
}

write_result(do.call(rbind, auc_horizon_results), "baseline_D_time_dependent_OS_AUC.csv")

# Baseline and week-specific D at the prespecified 36-month horizon.
roc_datasets <- list(
  Baseline_D = data.frame(
    study_id = patients$study_id,
    time = patients$os_months,
    event = patients$os_event,
    marker = patients$D_base
  )
)

for (timepoint in c("week2", "week3", "week4")) {
  current <- make_analysis_data(timepoint, "OS")
  roc_datasets[[paste0("Week_", sub("week", "", timepoint), "_D")]] <- data.frame(
    study_id = current$study_id,
    time = current$time,
    event = current$event,
    marker = current$D_current
  )
}

roc_36_results <- list()
for (marker_name in names(roc_datasets)) {
  current <- roc_datasets[[marker_name]]
  weights <- roc_weights(current$time, current$event, 36)
  estimate <- weighted_roc(current$marker, weights)
  roc_36_results[[marker_name]] <- data.frame(
    marker = marker_name,
    N = nrow(current),
    deaths_by_36_months = weights$cases,
    controls_beyond_36_months = weights$controls,
    AUC = estimate$AUC,
    Youden_cutoff = estimate$cutoff,
    sensitivity = estimate$sensitivity,
    specificity = estimate$specificity
  )
}
roc_36_results <- do.call(rbind, roc_36_results)
write_result(roc_36_results, "ROC_36month_OS_markers.csv")

# Paired patient-level bootstrap comparison of the four AUCs.
set.seed(RANDOM_SEED + 100L)
paired_auc <- matrix(
  NA_real_, nrow = N_BOOTSTRAP, ncol = length(roc_datasets),
  dimnames = list(NULL, names(roc_datasets))
)

for (bootstrap_index in seq_len(N_BOOTSTRAP)) {
  sampled_ids <- sample(patients$study_id, nrow(patients), replace = TRUE)

  for (marker_name in names(roc_datasets)) {
    current <- roc_datasets[[marker_name]]
    sample_index <- match(sampled_ids, current$study_id, nomatch = 0L)
    sample_index <- sample_index[sample_index > 0L]

    estimate <- tryCatch({
      sampled_data <- current[sample_index, ]
      weights <- roc_weights(sampled_data$time, sampled_data$event, 36)
      weighted_roc(sampled_data$marker, weights)
    }, error = function(e) NULL)

    if (!is.null(estimate)) {
      paired_auc[bootstrap_index, marker_name] <- estimate$AUC
    }
  }
}

observed_auc <- setNames(roc_36_results$AUC, roc_36_results$marker)
comparison_pairs <- combn(names(observed_auc), 2, simplify = FALSE)

auc_comparisons <- do.call(rbind, lapply(comparison_pairs, function(pair) {
  bootstrap_difference <- paired_auc[, pair[1]] - paired_auc[, pair[2]]
  raw_p <- min(1, 2 * min(
    mean(bootstrap_difference <= 0, na.rm = TRUE),
    mean(bootstrap_difference >= 0, na.rm = TRUE)
  ))

  data.frame(
    marker_1 = pair[1],
    marker_2 = pair[2],
    AUC_difference = observed_auc[pair[1]] - observed_auc[pair[2]],
    CI_low = quantile(bootstrap_difference, 0.025, na.rm = TRUE),
    CI_high = quantile(bootstrap_difference, 0.975, na.rm = TRUE),
    bootstrap_p = raw_p
  )
}))
auc_comparisons$Holm_adjusted_p <- p.adjust(auc_comparisons$bootstrap_p, method = "holm")
write_result(auc_comparisons, "ROC_36month_AUC_comparisons.csv")

# Stability of the reported baseline D cutoff.
cutoff_bootstrap <- auc_bootstrap_36$cutoff
cutoff_stability <- data.frame(
  reported_cutoff = REPORTED_BASELINE_CUTOFF,
  bootstrap_median = median(cutoff_bootstrap, na.rm = TRUE),
  bootstrap_Q1 = quantile(cutoff_bootstrap, 0.25, na.rm = TRUE),
  bootstrap_Q3 = quantile(cutoff_bootstrap, 0.75, na.rm = TRUE),
  bootstrap_2.5_percentile = quantile(cutoff_bootstrap, 0.025, na.rm = TRUE),
  bootstrap_97.5_percentile = quantile(cutoff_bootstrap, 0.975, na.rm = TRUE),
  bootstrap_SD = sd(cutoff_bootstrap, na.rm = TRUE),
  percent_within_10_percent = 100 * mean(
    abs(cutoff_bootstrap - REPORTED_BASELINE_CUTOFF) <= 0.10 * REPORTED_BASELINE_CUTOFF,
    na.rm = TRUE
  )
)
write_result(cutoff_stability, "baseline_D_cutoff_stability.csv")

# 10. Cutoff-based comparison with hemoglobin and NLR -------------------------
cutoff_model_results <- list()
cutoff_model_index <- 0L
cutoff_survival_results <- list()

cutoff_labels <- c(
  D_high = "Baseline D > 3.602",
  Hb_low = "Hemoglobin < 11 g/dL",
  NLR_high = "Neutrophil-to-lymphocyte ratio >= 3",
  age10 = "Age per 10 years",
  stage_high = "Stage IIIC2-IVB",
  non_sqcc = "Non-squamous histology",
  dose10 = "Equivalent dose in 2-Gy fractions per 10 Gy",
  duration7 = "Treatment duration per 7 days"
)

for (endpoint in ENDPOINTS) {
  endpoint_lower <- tolower(endpoint)
  analysis_data <- patients
  analysis_data$time <- analysis_data[[paste0(endpoint_lower, "_months")]]
  analysis_data$event <- analysis_data[[paste0(endpoint_lower, "_event")]]
  analysis_data$D_high <- as.integer(analysis_data$D_base > REPORTED_BASELINE_CUTOFF)
  analysis_data$Hb_low <- as.integer(analysis_data$baseline_Hb < 11)
  analysis_data$NLR_high <- as.integer(analysis_data$baseline_NLR >= 3)

  model_terms <- names(cutoff_labels)
  complete_columns <- c("time", "event", model_terms)
  analysis_data <- analysis_data[complete.cases(analysis_data[, complete_columns]), ]

  analysis_data$D_group <- factor(
    ifelse(analysis_data$D_high == 1, "High", "Low"),
    levels = c("Low", "High")
  )
  logrank_test <- survdiff(Surv(time, event) ~ D_group, data = analysis_data)
  cutoff_survival_results[[endpoint]] <- data.frame(
    endpoint = endpoint,
    low_D_N = sum(analysis_data$D_group == "Low"),
    high_D_N = sum(analysis_data$D_group == "High"),
    low_D_events = sum(analysis_data$event[analysis_data$D_group == "Low"]),
    high_D_events = sum(analysis_data$event[analysis_data$D_group == "High"]),
    log_rank_p = pchisq(logrank_test$chisq, df = 1, lower.tail = FALSE)
  )

  model_formula <- as.formula(paste("Surv(time, event) ~", paste(model_terms, collapse = " + ")))
  fitted_model <- coxph(model_formula, data = analysis_data, ties = "efron")
  model_summary <- summary(fitted_model)
  confidence_intervals <- exp(confint(fitted_model))

  for (term in model_terms) {
    cutoff_model_index <- cutoff_model_index + 1L
    cutoff_model_results[[cutoff_model_index]] <- data.frame(
      endpoint = endpoint,
      variable = cutoff_labels[term],
      adjusted_HR = exp(coef(fitted_model)[term]),
      CI_low = confidence_intervals[term, 1],
      CI_high = confidence_intervals[term, 2],
      Wald_p = model_summary$coefficients[term, "Pr(>|z|)"]
    )
  }
}

write_result(do.call(rbind, cutoff_model_results), "cutoff_multivariable_models.csv")
write_result(do.call(rbind, cutoff_survival_results), "baseline_D_cutoff_logrank.csv")

# 11. Session information ------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "sessionInfo.txt"))

cat("Public reproducibility analysis completed.\n")
cat("Results:", OUTPUT_DIR, "\n")
