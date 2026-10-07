# V7 design and artifact-recomputable metrics; no top-level fitting or writes.
# Source the frozen logistic_sglasso.R, pilot_utils.R and study_utils.R first.
# The generic simulate_logistic_sglasso_two_design_v1() is reused unchanged.

lsg_design_columns_v7 <- function() {
  c("design_id", "scenario", "n_train", "n_validation", "n_test", "groups",
    "group_size", "active_group_count", "rho_within", "rho_between",
    "prevalence", "linear_predictor_sd", "signal_pattern", "description")
}

lsg_design_checks_v7 <- function(design) {
  checks <- c(frozen_schema = is.data.frame(design) && identical(
    names(design), c(lsg_design_columns_v7(), "scenario_index")
  ))
  result <- function() data.frame(check = names(checks),
                                  passed = unname(checks))
  if (!all(checks)) return(result())
  numeric_columns <- c("n_train", "n_validation", "n_test", "groups",
                       "group_size", "active_group_count", "rho_within",
                       "rho_between", "prevalence", "linear_predictor_sd",
                       "scenario_index")
  checks["finite_numeric_design"] <- all(vapply(design[numeric_columns],
    function(x) is.numeric(x) && all(is.finite(x)), logical(1)))
  checks["eight_unique_scenarios"] <- nrow(design) == 8L &&
    !anyNA(design$scenario) && !anyDuplicated(design$scenario)
  if (!all(checks)) return(result())
  pattern <- c(rep("homogeneous", 3L), rep("weak_strong_mixed", 5L))
  rho <- c(0, 0.3, 0.9, 0, 0.3, 0.9, 0.3, 0.9)
  prevalence <- c(rep(0.5, 6L), 0.1, 0.1)
  scenario <- paste0(pattern, "_rhob_", sprintf("%.1f", rho),
                     ifelse(prevalence == 0.1, "_prev_0.1", ""))
  same <- function(x, target) isTRUE(all.equal(
    as.numeric(x), as.numeric(target), tolerance = 1e-12,
    check.attributes = FALSE
  ))
  checks["frozen_scenario_order"] <- identical(design$scenario, scenario) &&
    same(design$scenario_index, seq_len(8L))
  checks["frozen_signal_patterns"] <- identical(design$design_id, pattern) &&
    identical(design$signal_pattern, pattern)
  checks["frozen_sample_sizes"] <- all(design$n_train == 200L) &&
    all(design$n_validation == 200L) && all(design$n_test == 5000L)
  checks["frozen_group_structure"] <- all(design$groups == 200L) &&
    all(design$group_size == 3L) && all(design$active_group_count == 10L)
  checks["frozen_correlations"] <- all(abs(design$rho_within - 0.9) < 1e-12) &&
    same(design$rho_between, rho)
  checks["frozen_prevalence_and_signal_scale"] <-
    same(design$prevalence, prevalence) &&
    all(abs(design$linear_predictor_sd - 1.1) < 1e-12)
  checks["positive_definite_factor_covariance"] <-
    all(design$rho_between >= 0 &
          design$rho_between <= design$rho_within & design$rho_within < 1)
  checks["nonempty_descriptions"] <- is.character(design$description) &&
    !anyNA(design$description) && all(nzchar(trimws(design$description)))
  result()
}

lsg_read_design_v7 <- function(path) {
  design <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(names(design), lsg_design_columns_v7())) {
    stop("V7 design columns do not match the frozen schema.", call. = FALSE)
  }
  integers <- c("n_train", "n_validation", "n_test", "groups", "group_size",
                "active_group_count")
  for (name in c(integers, "rho_within", "rho_between", "prevalence",
                 "linear_predictor_sd")) {
    value <- suppressWarnings(as.numeric(design[[name]]))
    if (any(!is.finite(value)) ||
        (name %in% integers && any(value != floor(value)))) {
      stop("Invalid V7 numeric design column: ", name, call. = FALSE)
    }
    design[[name]] <- if (name %in% integers) as.integer(value) else value
  }
  design$scenario_index <- seq_len(nrow(design))
  checks <- lsg_design_checks_v7(design)
  if (!all(checks$passed)) {
    stop("V7 scientific design changed: ",
         paste(checks$check[!checks$passed], collapse = ", "), call. = FALSE)
  }
  design
}

# Counts are computed from explicit binary indicators, never coefficient size
# as a continuous score. mltools defines zero-denominator MCC as zero; keep
# an independent degeneracy flag so this convention remains visible.
lsg_binary_scores_v7 <- function(actual, predicted) {
  tp <- sum(actual & predicted)
  fp <- sum(!actual & predicted)
  tn <- sum(!actual & !predicted)
  fn <- sum(actual & !predicted)
  out <- list(tp = tp, fp = fp, tn = tn, fn = fn,
    mcc = as.numeric(mltools::mcc(TP = tp, FP = fp, TN = tn, FN = fn)),
    mcc_degenerate = any(c(tp + fp, tp + fn, tn + fp, tn + fn) == 0))
  add_ratio <- function(name, numerator, denominator, reason) {
    defined <- denominator > 0
    out[[name]] <<- if (defined) numerator / denominator else NA_real_
    out[[paste0(name, "_defined")]] <<- defined
    out[[paste0(name, "_undefined_reason")]] <<- if (defined) "" else reason
  }
  add_ratio("tpr", tp, tp + fn, "no_actual_positives")
  add_ratio("tnr", tn, tn + fp, "no_actual_negatives")
  add_ratio("precision", tp, tp + fp, "no_predicted_positives")
  add_ratio("f1", 2 * tp, 2 * tp + fp + fn,
            "no_actual_or_predicted_positives")
  out$balanced_accuracy_defined <- out$tpr_defined && out$tnr_defined
  out$balanced_accuracy <- if (out$balanced_accuracy_defined) {
    (out$tpr + out$tnr) / 2
  } else NA_real_
  out$balanced_accuracy_undefined_reason <- if (out$balanced_accuracy_defined) {
    ""
  } else "only_one_actual_class"
  # FDP = FP/max(1, discoveries), consistent with the frozen group_fdr column.
  out$fdp <- fp / max(1, tp + fp)
  out$exact_support <- fp == 0 && fn == 0
  out$selected_all <- tp + fp == length(actual)
  out$selected_none <- tp + fp == 0
  out
}

# Average precision is non-interpolated: sum over descending unique score
# thresholds of (recall increment) * precision at the complete tied-score
# block. It is not trapezoidal PR-AUC; ties cannot depend on row order.
lsg_average_precision_v7 <- function(y, probability) {
  positives <- sum(y == 1)
  if (positives == 0L) return(NA_real_)
  order_index <- order(probability, decreasing = TRUE)
  ordered_y <- y[order_index]
  block_end <- cumsum(rle(probability[order_index])$lengths)
  cumulative_tp <- cumsum(ordered_y)[block_end]
  sum(diff(c(0, cumulative_tp)) / positives * cumulative_tp / block_end)
}

lsg_metrics_v7 <- function(y, probability, coefficient, true_coefficient,
                           group, eta, true_eta, threshold = 1e-8) {
  if (!requireNamespace("mltools", quietly = TRUE)) {
    stop("V7 metrics require the installed mltools package.", call. = FALSE)
  }
  if (!is.numeric(y) && !is.logical(y)) {
    stop("y must already be numeric/logical binary 0/1.", call. = FALSE)
  }
  numeric_inputs <- list(probability, coefficient, true_coefficient, eta,
                         true_eta, threshold)
  if (any(!vapply(numeric_inputs, is.numeric, logical(1))) ||
      any(!vapply(numeric_inputs, function(x) all(is.finite(x)), logical(1)))) {
    stop("Metric numeric inputs must be finite numeric vectors.", call. = FALSE)
  }
  y <- as.numeric(y)
  probability <- as.numeric(probability)
  coefficient <- as.numeric(coefficient)
  true_coefficient <- as.numeric(true_coefficient)
  eta <- as.numeric(eta)
  true_eta <- as.numeric(true_eta)
  if (!length(y) || anyNA(y) || any(!y %in% c(0, 1)) ||
      length(probability) != length(y) || length(eta) != length(y) ||
      length(true_eta) != length(y) || any(probability < 0 | probability > 1) ||
      length(threshold) != 1L || threshold < 0 || !length(group) ||
      length(coefficient) != length(group) + 1L ||
      length(true_coefficient) != length(coefficient) || anyNA(group) ||
      any(!nzchar(as.character(group)))) {
    stop("Invalid metric dimensions, binary labels or probability range.",
         call. = FALSE)
  }
  group <- as.character(group)
  predicted_beta <- abs(coefficient[-1L]) > threshold
  active_beta <- true_coefficient[-1L] != 0
  labels <- unique(group)
  predicted_group <- vapply(labels, function(g) any(predicted_beta[group == g]),
                           logical(1))
  active_group <- vapply(labels, function(g) any(active_beta[group == g]),
                        logical(1))
  selected <- labels[predicted_group]
  active <- labels[active_group]
  base <- evaluate_binary_method(y, probability, selected, active)
  response_scores <- lsg_binary_scores_v7(y == 1, probability >= 0.5)
  predictor_scores <- lsg_binary_scores_v7(active_beta, predicted_beta)
  group_scores <- lsg_binary_scores_v7(active_group, predicted_group)
  out <- as.list(base)
  attach_scores <- function(scores, prefix) {
    for (name in names(scores)) out[[paste0(prefix, name)]] <<- scores[[name]]
  }
  attach_scores(response_scores, "y_")
  attach_scores(predictor_scores, "predictor_")
  attach_scores(group_scores, "group_")
  for (pair in list(c("tpr", "sensitivity"), c("tnr", "specificity"))) {
    for (suffix in c("", "_defined", "_undefined_reason")) {
      out[[paste0("y_", pair[2L], suffix)]] <-
        response_scores[[paste0(pair[1L], suffix)]]
    }
  }
  out$y_predicted_positive_rate <- mean(probability >= 0.5)
  out$selected_predictors <- sum(predicted_beta)
  out$group_fdr <- group_scores$fdp
  out$true_positive_groups <- group_scores$tp
  out$false_positive_groups <- group_scores$fp
  out$exact_group_support <- group_scores$exact_support
  out$auc_defined <- length(unique(y)) == 2L
  out$auc_undefined_reason <- if (out$auc_defined) "" else "only_one_response_class"
  out$average_precision <- lsg_average_precision_v7(y, probability)
  out$average_precision_defined <- any(y == 1)
  out$average_precision_undefined_reason <- if (out$average_precision_defined) {
    ""
  } else "no_response_positives"
  delta <- coefficient[-1L] - true_coefficient[-1L]
  out$coefficient_squared_l2_error <- sum(delta^2)
  out$coefficient_l2_error <- sqrt(out$coefficient_squared_l2_error)
  out$coefficient_l1_error <- sum(abs(delta))
  truth_squared_norm <- sum(true_coefficient[-1L]^2)
  out$coefficient_relative_squared_l2_error_defined <- truth_squared_norm > 0
  out$coefficient_relative_squared_l2_error <- if (truth_squared_norm > 0) {
    out$coefficient_squared_l2_error / truth_squared_norm
  } else NA_real_
  out$coefficient_relative_squared_l2_error_undefined_reason <-
    if (truth_squared_norm > 0) "" else "zero_true_coefficient_norm"
  out$intercept_absolute_error <- abs(coefficient[1L] - true_coefficient[1L])
  out$linear_predictor_rmse <- sqrt(mean((eta - true_eta)^2))
  calibration <- calibration_statistics(y, probability)
  link <- stats::qlogis(pmin(pmax(probability, 1e-8), 1 - 1e-8))
  for (name in c("intercept", "slope")) {
    reason <- if (length(unique(y)) < 2L) {
      "only_one_response_class"
    } else if (name == "slope" && length(unique(link)) < 2L) {
      "constant_predicted_log_odds"
    } else if (!is.finite(calibration[[name]])) {
      "calibration_fit_not_finite"
    } else ""
    metric <- paste0("calibration_", name)
    out[[metric]] <- if (nzchar(reason)) NA_real_ else unname(calibration[[name]])
    out[[paste0(metric, "_defined")]] <- !nzchar(reason)
    out[[paste0(metric, "_undefined_reason")]] <- reason
  }
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

lsg_metric_registry_v7 <- function() {
  make <- function(metric, family, direction) data.frame(
    metric = metric, family = family, direction = direction,
    stringsAsFactors = FALSE
  )
  rows <- list(
    make(c("log_loss", "brier", "classification_error"), "prediction", "lower"),
    make(c("auc", "average_precision"), "prediction", "higher"),
    make(c("y_mcc", "y_sensitivity", "y_specificity", "y_precision", "y_f1",
           "y_balanced_accuracy"), "classification", "higher"),
    make("calibration_intercept", "prediction", "target_0"),
    make("calibration_slope", "prediction", "target_1"),
    make(c("coefficient_squared_l2_error", "coefficient_relative_squared_l2_error",
           "coefficient_l2_error", "coefficient_l1_error", "intercept_absolute_error",
           "linear_predictor_rmse"), "coefficient_estimation", "lower"),
    make(c("selected_groups", "selected_predictors", "y_predicted_positive_rate"),
         "selection_size", "descriptive")
  )
  for (level in c("group", "predictor")) {
    family <- paste0(level, "_selection")
    rows[[length(rows) + 1L]] <- make(paste0(level, "_", c(
      "mcc", "tpr", "tnr", "precision", "f1", "exact_support"
    )), family, "higher")
    rows[[length(rows) + 1L]] <- make(paste0(level, "_fdp"), family, "lower")
    rows[[length(rows) + 1L]] <- make(paste0(level, "_", c(
      "mcc_degenerate", "selected_all", "selected_none"
    )), family, "descriptive")
  }
  rows[[length(rows) + 1L]] <- make(
    c("y_mcc_degenerate", "y_selected_all", "y_selected_none"),
    "classification", "descriptive"
  )
  do.call(rbind, rows)
}

lsg_metrics_unit_checks_v7 <- function() {
  if (!requireNamespace("mltools", quietly = TRUE)) {
    stop("Install the approved mltools dependency before running V7 checks.",
         call. = FALSE)
  }
  equal <- function(x, y) isTRUE(all.equal(x, y, tolerance = 1e-12,
                                          check.attributes = FALSE))
  actual <- c(TRUE, TRUE, FALSE, FALSE)
  perfect <- lsg_binary_scores_v7(actual, actual)
  reversed <- lsg_binary_scores_v7(actual, !actual)
  all_selected <- lsg_binary_scores_v7(actual, rep(TRUE, 4L))
  none_selected <- lsg_binary_scores_v7(actual, rep(FALSE, 4L))
  chance <- lsg_binary_scores_v7(actual, c(TRUE, FALSE, TRUE, FALSE))
  y <- rep(c(0, 1, 1, 0, 1, 0), 4L)
  probability <- rep(c(0.2, 0.3, 0.5, 0.5, 0.7, 0.8), 4L)
  coefficient <- c(0.2, 0.5, 0.5, 0, 0, 0.2, 0.2)
  true_coefficient <- c(-0.1, 0.7, 0.7, 0.3, 0.3, 0, 0)
  group <- rep(seq_len(3L), each = 2L)
  eta <- stats::qlogis(probability)
  true_eta <- eta + seq(-0.2, 0.2, length.out = length(y))
  row <- lsg_metrics_v7(y, probability, coefficient, true_coefficient, group,
                        eta, true_eta)
  empty <- lsg_metrics_v7(y, rep(0.4, length(y)), coefficient * 0,
    true_coefficient, group, rep(stats::qlogis(0.4), length(y)), true_eta)
  tied <- lsg_average_precision_v7(c(1, 0, 1, 0), c(0.9, 0.9, 0.2, 0.1))
  threshold_row <- lsg_metrics_v7(y, probability, c(0, 1e-8, 1.01e-8, 0, 0, 0, 0),
                                 true_coefficient, group, eta, true_eta)
  no_truth <- lsg_metrics_v7(y, probability, coefficient, true_coefficient * 0,
                            group, eta, true_eta)
  X <- matrix(sin(seq_len(length(y) * length(group))), nrow = length(y))
  legacy_eta <- coefficient[1L] + drop(X %*% coefficient[-1L])
  legacy_true_eta <- true_coefficient[1L] + drop(X %*% true_coefficient[-1L])
  legacy_probability <- stats::plogis(legacy_eta)
  legacy <- evaluate_study_method(y, legacy_probability, c("1", "3"),
    c("1", "2"), coefficient, true_coefficient, X)
  legacy_recomputed <- lsg_metrics_v7(y, legacy_probability, coefficient,
    true_coefficient, group, legacy_eta, legacy_true_eta)
  rejects <- function(code) inherits(try(force(code), silent = TRUE), "try-error")
  vector_agreement <- vapply(seq_len(16L) - 1L, function(mask) {
    prediction <- as.logical(intToBits(mask)[seq_len(4L)])
    count_score <- lsg_binary_scores_v7(actual, prediction)$mcc
    equal(count_score, mltools::mcc(preds = as.integer(prediction),
                                    actuals = as.integer(actual)))
  }, logical(1))
  checks <- c(
    mltools_perfect_mcc_one = perfect$mcc == 1,
    mltools_reversed_mcc_minus_one = reversed$mcc == -1,
    mltools_zero_correlation = chance$mcc == 0 && !chance$mcc_degenerate,
    mltools_all_selected_zero_flagged = all_selected$mcc == 0 &&
      all_selected$mcc_degenerate && all_selected$selected_all,
    mltools_none_selected_zero_flagged = none_selected$mcc == 0 &&
      none_selected$mcc_degenerate && none_selected$selected_none,
    mltools_count_vector_api_agree = equal(row$y_mcc, mltools::mcc(
      preds = as.integer(probability >= 0.5), actuals = as.integer(y)
    )),
    mltools_count_vector_all_binary_patterns = all(vector_agreement),
    study_all_200_groups_selected_mcc_zero =
      mltools::mcc(TP = 10, FP = 190, TN = 0, FN = 0) == 0,
    confusion_counts_sum_to_total = with(row, y_tp + y_fp + y_tn + y_fn) == length(y),
    grouped_support_counts_exact = with(row, group_tp == 1 && group_fp == 1 &&
      group_tn == 0 && group_fn == 1 && selected_groups == 2),
    equal_size_block_support_mcc_equal = equal(row$group_mcc, row$predictor_mcc),
    strict_original_scale_support_threshold = threshold_row$selected_predictors == 1,
    squared_l2_matches_l2_squared = equal(row$coefficient_squared_l2_error,
                                          row$coefficient_l2_error^2),
    coefficient_error_excludes_intercept = equal(row$coefficient_squared_l2_error,
      sum((coefficient[-1L] - true_coefficient[-1L])^2)),
    linear_predictor_rmse_uses_saved_eta = equal(row$linear_predictor_rmse,
      sqrt(mean((eta - true_eta)^2))),
    no_selected_precision_undefined = is.na(empty$predictor_precision) &&
      !empty$predictor_precision_defined &&
      empty$predictor_precision_undefined_reason == "no_predicted_positives",
    no_selected_fdp_zero = empty$group_fdp == 0 && empty$group_fdr == 0,
    constant_prediction_slope_undefined = is.na(empty$calibration_slope) &&
      empty$calibration_slope_undefined_reason == "constant_predicted_log_odds",
    zero_truth_relative_error_undefined =
      is.na(no_truth$coefficient_relative_squared_l2_error) &&
      !no_truth$coefficient_relative_squared_l2_error_defined,
    average_precision_ties_hand_calculation = equal(tied, 7 / 12),
    average_precision_tie_permutation_invariant = equal(tied,
      lsg_average_precision_v7(c(0, 1, 1, 0), c(0.9, 0.9, 0.2, 0.1))),
    average_precision_constant_equals_prevalence = equal(
      lsg_average_precision_v7(y, rep(0.5, length(y))), mean(y)),
    average_precision_no_positives_undefined = is.na(
      lsg_average_precision_v7(rep(0, 4L), c(0.9, 0.9, 0.2, 0.1))),
    average_precision_perfect_one = equal(
      lsg_average_precision_v7(c(1, 1, 0, 0), c(0.9, 0.8, 0.2, 0.1)), 1),
    frozen_legacy_metric_columns_preserved =
      equal(legacy_recomputed[names(legacy)], legacy),
    missing_outcomes_rejected = rejects(lsg_metrics_v7(c(y[-1L], NA_real_),
      probability, coefficient, true_coefficient, group, eta, true_eta)),
    invalid_probability_rejected = rejects(lsg_metrics_v7(y,
      c(probability[-1L], 1.01), coefficient, true_coefficient, group, eta, true_eta)),
    metric_registry_complete = all(lsg_metric_registry_v7()$metric %in% names(row))
  )
  data.frame(check = names(checks), passed = unname(checks), stringsAsFactors = FALSE)
}
