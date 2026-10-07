# Diagnostic-only replay of the three V7 Logistic Group Lasso failures.
#
# This module does not select a model and does not evaluate the test sample. It
# preserves the raw grpreg path before select_valid_external_grid_v2() can
# reject every point. Source the frozen V7 workflow before sourcing this file.

lsg_grpreg_diagnostic_cases_v8 <- function() {
  data.frame(
    task_id = c(126L, 134L, 150L),
    scenario = c(
      "weak_strong_mixed_rhob_0.3_prev_0.1",
      "weak_strong_mixed_rhob_0.3_prev_0.1",
      "weak_strong_mixed_rhob_0.9_prev_0.1"
    ),
    replication = c(6L, 14L, 10L),
    seed = c(730330837L, 730330845L, 730340841L),
    stringsAsFactors = FALSE
  )
}


lsg_grpreg_require_v8 <- function(root = NULL) {
  functions <- c(
    "lsg_configuration_v7", "lsg_design_for_stage_v7", "lsg_task_grid_v7",
    "lsg_with_rng_v7", "simulate_logistic_sglasso_two_design_v1",
    "lsg_data_fingerprint_v7",
    "validate_penalized_benchmark_grid_v1", "classify_grpreg_path_v2",
    "binary_log_loss"
  )
  diagnostic_environment <- environment(lsg_grpreg_require_v8)
  missing <- functions[!vapply(functions, exists, logical(1),
                               envir = diagnostic_environment,
                               mode = "function", inherits = TRUE)]
  if (length(missing)) {
    stop(
      "Source the frozen V7 workflow before the V8 grpreg diagnostic; ",
      "missing function(s): ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.null(root) &&
      (!is.character(root) || length(root) != 1L || !dir.exists(root))) {
    stop("root must be one existing project directory.", call. = FALSE)
  }
  if (!requireNamespace("grpreg", quietly = TRUE)) {
    stop("The installed grpreg package is required; no installation is attempted.",
         call. = FALSE)
  }
  invisible(TRUE)
}


lsg_grpreg_tasks_v8 <- function(root, task_ids = c(126L, 134L, 150L)) {
  lsg_grpreg_require_v8(root)
  task_ids <- as.integer(task_ids)
  frozen_cases <- lsg_grpreg_diagnostic_cases_v8()
  if (!length(task_ids) || anyNA(task_ids) || anyDuplicated(task_ids) ||
      any(!task_ids %in% frozen_cases$task_id)) {
    stop("V8 grpreg diagnostics are restricted to frozen task IDs 126, 134 and 150.",
         call. = FALSE)
  }

  configuration <- lsg_configuration_v7("pilot")
  design <- lsg_design_for_stage_v7(root, "pilot")
  tasks <- lsg_task_grid_v7(design, configuration)
  selected <- tasks[match(task_ids, tasks$task_id), , drop = FALSE]
  expected <- frozen_cases[match(task_ids, frozen_cases$task_id), , drop = FALSE]
  if (!identical(as.integer(selected$task_id), expected$task_id) ||
      !identical(as.character(selected$scenario), expected$scenario) ||
      !identical(as.integer(selected$replication), expected$replication) ||
      !identical(as.integer(selected$seed), expected$seed)) {
    stop("The frozen V7 task map no longer matches the approved diagnostic cases.",
         call. = FALSE)
  }
  list(tasks = selected, design = design, configuration = configuration)
}


lsg_capture_condition_v8 <- function(expression) {
  warnings <- character(0)
  error <- ""
  elapsed <- system.time({
    value <- tryCatch(
      withCallingHandlers(
        force(expression),
        warning = function(condition) {
          warnings <<- c(warnings, conditionMessage(condition))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(condition) {
        error <<- conditionMessage(condition)
        NULL
      }
    )
  })[["elapsed"]]
  list(
    value = value,
    warnings = unique(warnings[!is.na(warnings) & nzchar(warnings)]),
    error = error,
    elapsed_seconds = as.numeric(elapsed)
  )
}


lsg_align_numeric_v8 <- function(value, length_out) {
  out <- rep(NA_real_, length_out)
  if (length_out && length(value)) {
    take <- seq_len(min(length_out, length(value)))
    out[take] <- as.numeric(value)[take]
  }
  out
}


lsg_align_logical_v8 <- function(value, length_out) {
  out <- rep(NA, length_out)
  if (length_out && length(value)) {
    take <- seq_len(min(length_out, length(value)))
    out[take] <- as.logical(value)[take]
  }
  out
}


lsg_grpreg_usable_prefix_v8 <- function(
    iterations,
    max_iterations,
    finite_coefficient,
    finite_probability,
    finite_validation_loss,
    warnings = character(0)
) {
  length_out <- length(iterations)
  iterations <- as.numeric(iterations)
  max_iterations <- as.numeric(max_iterations)
  warnings <- unique(as.character(warnings))
  warnings <- warnings[!is.na(warnings) & nzchar(warnings)]
  if (!length_out || !is.finite(max_iterations) || max_iterations < 1 ||
      length(finite_coefficient) != length_out ||
      length(finite_probability) != length_out ||
      length(finite_validation_loss) != length_out) {
    stop("Invalid usable-prefix inputs.", call. = FALSE)
  }

  iteration_warning <- grepl(
    "Algorithm failed to converge for all values of lambda",
    warnings,
    fixed = TRUE
  )
  other_warnings <- warnings[!iteration_warning]
  valid_iterations <- is.finite(iterations) & iterations >= 0
  cumulative_iterations <- if (all(valid_iterations)) {
    cumsum(iterations)
  } else {
    rep(NA_real_, length_out)
  }
  total_iteration_limit_reached <- all(valid_iterations) &&
    sum(iterations) >= max_iterations
  warning_compatible <- !length(other_warnings) &&
    (!any(iteration_warning) || total_iteration_limit_reached)
  finite_point <- !is.na(finite_coefficient) & finite_coefficient &
    !is.na(finite_probability) & finite_probability &
    !is.na(finite_validation_loss) & finite_validation_loss
  completed_before_budget <- valid_iterations &
    is.finite(cumulative_iterations) & cumulative_iterations < max_iterations
  candidate <- warning_compatible & finite_point & completed_before_budget
  prefix <- cumprod(as.integer(candidate)) == 1L

  list(
    cumulative_iterations = cumulative_iterations,
    completed_before_total_iteration_budget = completed_before_budget,
    finite_point = finite_point,
    warning_compatible = warning_compatible,
    candidate = candidate,
    member = prefix,
    length = sum(prefix),
    total_iteration_limit_reached = total_iteration_limit_reached,
    iteration_budget_warning = any(iteration_warning),
    other_warnings = other_warnings,
    rule = paste(
      "maximal leading sequence with finite coefficients, probabilities and",
      "validation loss, and cumulative grpreg iterations strictly below max.iter;",
      "only the grpreg total-iteration warning is prefix-compatible"
    )
  )
}


lsg_grpreg_path_label_v8 <- function(
    fit_error,
    prediction_error,
    schema_compatible,
    classification,
    usable_prefix_length,
    returned_path_length
) {
  if (nzchar(fit_error)) return("fit_error")
  if (!schema_compatible) return("incompatible_raw_path_schema")
  if (nzchar(prediction_error)) return("prediction_error")
  if (is.list(classification) &&
      isTRUE(classification$path_termination_acceptable)) {
    return("legacy_path_acceptable")
  }
  if (usable_prefix_length > 0L && usable_prefix_length < returned_path_length) {
    return("legacy_path_rejected_with_finite_converged_prefix")
  }
  if (usable_prefix_length == returned_path_length && returned_path_length > 0L) {
    return("legacy_path_rejected_but_all_returned_points_prefix_usable")
  }
  "legacy_path_rejected_without_usable_prefix"
}


lsg_fit_raw_grpreg_group_lasso_v8 <- function(
    X,
    y,
    group,
    X_validation,
    y_validation,
    nlambda,
    lambda_min_ratio,
    tolerance = 1e-7,
    max_iterations = 1000000L
) {
  lsg_grpreg_require_v8()
  input <- validate_penalized_benchmark_grid_v1(
    X, y, group, X_validation, y_validation,
    alpha_grid = 1, nlambda = nlambda,
    lambda_min_ratio = lambda_min_ratio
  )
  max_iterations <- as.integer(max_iterations)
  tolerance <- as.numeric(tolerance)
  if (length(max_iterations) != 1L || is.na(max_iterations) ||
      max_iterations < 1L || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0) {
    stop("Invalid grpreg diagnostic numerical controls.", call. = FALSE)
  }

  group_factor <- factor(input$group, levels = unique(input$group))
  group_multiplier <- sqrt(as.numeric(table(group_factor)))
  fit_phase <- lsg_capture_condition_v8(grpreg::grpreg(
    X = input$X,
    y = input$y,
    group = input$group,
    penalty = "grLasso",
    family = "binomial",
    nlambda = input$nlambda,
    lambda.min = input$lambda_min_ratio,
    log.lambda = TRUE,
    alpha = 1,
    eps = tolerance,
    max.iter = max_iterations,
    dfmax = ncol(input$X),
    gmax = length(unique(input$group)),
    gamma = 3,
    group.multiplier = group_multiplier,
    warn = TRUE,
    returnX = FALSE
  ))
  fit <- fit_phase$value

  empty_summary <- function() data.frame(
    status = "fit_error",
    requested_path_length = input$nlambda,
    returned_path_length = 0L,
    coefficient_path_columns = 0L,
    iteration_vector_length = 0L,
    prediction_path_columns = 0L,
    raw_schema_compatible = FALSE,
    fit_error = fit_phase$error,
    prediction_error = "not_attempted_after_fit_error",
    fit_warning_count = length(fit_phase$warnings),
    prediction_warning_count = 0L,
    warnings = paste(fit_phase$warnings, collapse = " | "),
    iteration_sum = NA_real_,
    iteration_max = NA_real_,
    path_complete = FALSE,
    total_iteration_limit_reached = NA,
    legacy_path_termination_acceptable = FALSE,
    legacy_unexpected_warnings = "",
    raw_finite_point_count = 0L,
    legacy_numerically_eligible_count = 0L,
    usable_prefix_length = 0L,
    usable_prefix_best_lambda_index = NA_integer_,
    usable_prefix_best_validation_log_loss = NA_real_,
    fit_elapsed_seconds = fit_phase$elapsed_seconds,
    prediction_elapsed_seconds = 0,
    stringsAsFactors = FALSE
  )
  if (is.null(fit)) {
    return(list(
      schema_version = "raw_grpreg_group_lasso_diagnostic_v8",
      summary = empty_summary(),
      points = data.frame(),
      raw = list(fit = NULL, probability = NULL, iterations = numeric(0),
                 lambda = numeric(0), fit_warnings = fit_phase$warnings,
                 prediction_warnings = character(0)),
      usable_prefix_rule = NA_character_
    ))
  }

  prediction_phase <- lsg_capture_condition_v8(as.matrix(stats::predict(
    fit,
    X = input$X_validation,
    type = "response"
  )))
  probability <- prediction_phase$value
  lambda <- as.numeric(fit$lambda)
  iterations_raw <- as.numeric(fit$iter)
  coefficient <- tryCatch(as.matrix(fit$beta), error = function(e) NULL)
  path_length <- length(lambda)
  coefficient_columns <- if (is.null(coefficient)) 0L else ncol(coefficient)
  prediction_columns <- if (is.null(probability)) 0L else ncol(probability)
  coefficient_rows_ok <- !is.null(coefficient) &&
    nrow(coefficient) == ncol(input$X) + 1L
  coefficient_columns_ok <- coefficient_columns == path_length
  iteration_length_ok <- length(iterations_raw) == path_length
  prediction_shape_ok <- !is.null(probability) &&
    nrow(probability) == length(input$y_validation) &&
    prediction_columns == path_length
  schema_compatible <- path_length > 0L && coefficient_rows_ok &&
    coefficient_columns_ok && iteration_length_ok && prediction_shape_ok

  finite_coefficient <- rep(NA, path_length)
  if (coefficient_rows_ok && coefficient_columns_ok) {
    finite_coefficient <- apply(coefficient, 2L, function(value) {
      all(is.finite(value))
    })
  }
  finite_probability <- rep(NA, path_length)
  validation_loss <- rep(NA_real_, path_length)
  if (prediction_shape_ok) {
    finite_probability <- apply(probability, 2L, function(value) {
      all(is.finite(value)) && all(value >= 0 & value <= 1)
    })
    validation_loss <- vapply(seq_len(path_length), function(index) {
      if (isTRUE(finite_probability[index])) {
        binary_log_loss(input$y_validation, probability[, index])
      } else {
        NA_real_
      }
    }, numeric(1))
  }
  finite_validation_loss <- is.finite(validation_loss)
  warnings <- unique(c(fit_phase$warnings, prediction_phase$warnings))

  classification_error <- ""
  classification <- tryCatch(
    classify_grpreg_path_v2(
      penalty = "grLasso",
      iterations = iterations_raw,
      requested_path_length = input$nlambda,
      max_iterations = max_iterations,
      warnings = warnings
    ),
    error = function(condition) {
      classification_error <<- conditionMessage(condition)
      NULL
    }
  )
  prefix <- if (iteration_length_ok) {
    lsg_grpreg_usable_prefix_v8(
      iterations = iterations_raw,
      max_iterations = max_iterations,
      finite_coefficient = finite_coefficient,
      finite_probability = finite_probability,
      finite_validation_loss = finite_validation_loss,
      warnings = warnings
    )
  } else {
    list(
      cumulative_iterations = rep(NA_real_, path_length),
      completed_before_total_iteration_budget = rep(FALSE, path_length),
      finite_point = rep(FALSE, path_length),
      warning_compatible = FALSE,
      candidate = rep(FALSE, path_length),
      member = rep(FALSE, path_length),
      length = 0L,
      total_iteration_limit_reached = NA,
      iteration_budget_warning = any(grepl(
        "Algorithm failed to converge for all values of lambda",
        warnings,
        fixed = TRUE
      )),
      other_warnings = warnings,
      rule = "unavailable because fit$iter length did not match fit$lambda"
    )
  }

  legacy_point_converged <- if (is.list(classification)) {
    lsg_align_logical_v8(classification$point_converged, path_length)
  } else {
    rep(NA, path_length)
  }
  legacy_boundary_excluded <- if (is.list(classification)) {
    lsg_align_logical_v8(
      classification$returned_lower_boundary_excluded, path_length
    )
  } else {
    rep(NA, path_length)
  }
  legacy_acceptable <- is.list(classification) &&
    isTRUE(classification$path_termination_acceptable)
  legacy_eligible <- legacy_acceptable &
    !is.na(legacy_point_converged) & legacy_point_converged &
    !is.na(legacy_boundary_excluded) & !legacy_boundary_excluded &
    !is.na(finite_coefficient) & finite_coefficient &
    !is.na(finite_probability) & finite_probability &
    finite_validation_loss

  points <- data.frame(
    lambda_index = seq_len(path_length),
    lambda = lambda,
    lambda_fraction = if (path_length && is.finite(lambda[1L]) &&
                            lambda[1L] != 0) lambda / lambda[1L] else NA_real_,
    iterations = lsg_align_numeric_v8(iterations_raw, path_length),
    cumulative_iterations = prefix$cumulative_iterations,
    finite_coefficient = finite_coefficient,
    finite_probability = finite_probability,
    validation_log_loss = validation_loss,
    finite_validation_loss = finite_validation_loss,
    legacy_point_converged = legacy_point_converged,
    legacy_returned_lower_boundary_excluded = legacy_boundary_excluded,
    legacy_numerically_eligible = legacy_eligible,
    completed_before_total_iteration_budget =
      prefix$completed_before_total_iteration_budget,
    usable_prefix_candidate = prefix$candidate,
    usable_prefix_member = prefix$member,
    stringsAsFactors = FALSE
  )
  prefix_rows <- which(points$usable_prefix_member)
  best_prefix <- if (length(prefix_rows)) {
    prefix_rows[which.min(points$validation_log_loss[prefix_rows])]
  } else {
    NA_integer_
  }
  status <- lsg_grpreg_path_label_v8(
    fit_error = fit_phase$error,
    prediction_error = prediction_phase$error,
    schema_compatible = schema_compatible,
    classification = classification,
    usable_prefix_length = prefix$length,
    returned_path_length = path_length
  )
  summary <- data.frame(
    status = status,
    requested_path_length = input$nlambda,
    returned_path_length = path_length,
    coefficient_path_columns = coefficient_columns,
    iteration_vector_length = length(iterations_raw),
    prediction_path_columns = prediction_columns,
    raw_schema_compatible = schema_compatible,
    fit_error = fit_phase$error,
    prediction_error = prediction_phase$error,
    fit_warning_count = length(fit_phase$warnings),
    prediction_warning_count = length(prediction_phase$warnings),
    warnings = paste(warnings, collapse = " | "),
    iteration_sum = if (length(iterations_raw) &&
                           all(is.finite(iterations_raw))) {
      sum(iterations_raw)
    } else NA_real_,
    iteration_max = if (length(iterations_raw) &&
                           all(is.finite(iterations_raw))) {
      max(iterations_raw)
    } else NA_real_,
    path_complete = path_length == input$nlambda,
    total_iteration_limit_reached = if (is.list(classification)) {
      classification$total_iteration_limit_reached
    } else prefix$total_iteration_limit_reached,
    legacy_path_termination_acceptable = legacy_acceptable,
    legacy_unexpected_warnings = if (is.list(classification)) {
      paste(classification$unexpected_warnings, collapse = " | ")
    } else classification_error,
    raw_finite_point_count = sum(prefix$finite_point),
    legacy_numerically_eligible_count = sum(legacy_eligible),
    usable_prefix_length = prefix$length,
    usable_prefix_best_lambda_index = best_prefix,
    usable_prefix_best_validation_log_loss = if (is.na(best_prefix)) {
      NA_real_
    } else points$validation_log_loss[best_prefix],
    fit_elapsed_seconds = fit_phase$elapsed_seconds,
    prediction_elapsed_seconds = prediction_phase$elapsed_seconds,
    stringsAsFactors = FALSE
  )

  list(
    schema_version = "raw_grpreg_group_lasso_diagnostic_v8",
    summary = summary,
    points = points,
    raw = list(
      fit = fit,
      probability = probability,
      lambda = lambda,
      iterations = iterations_raw,
      group_multiplier = group_multiplier,
      fit_warnings = fit_phase$warnings,
      prediction_warnings = prediction_phase$warnings,
      path_classification = classification,
      path_classification_error = classification_error
    ),
    usable_prefix_rule = prefix$rule
  )
}


lsg_run_grpreg_task_v8 <- function(root, task_id) {
  frozen <- lsg_grpreg_tasks_v8(root, task_id)
  task <- frozen$tasks[1L, , drop = FALSE]
  scenario <- frozen$design[
    frozen$design$scenario_index == task$scenario_index,
    ,
    drop = FALSE
  ]
  generated <- lsg_with_rng_v7(
    task$seed,
    simulate_logistic_sglasso_two_design_v1(scenario, task$seed)
  )
  tuning_data <- generated[c(
    "X_train", "y_train", "X_validation", "y_validation", "group"
  )]
  rm(generated)
  data_sha256 <- list(
    training = lsg_data_fingerprint_v7(list(
      X = tuning_data$X_train, y = tuning_data$y_train,
      group = tuning_data$group
    )),
    validation = lsg_data_fingerprint_v7(list(
      X = tuning_data$X_validation, y = tuning_data$y_validation
    ))
  )
  configuration <- frozen$configuration
  diagnostic <- lsg_fit_raw_grpreg_group_lasso_v8(
    X = tuning_data$X_train,
    y = tuning_data$y_train,
    group = tuning_data$group,
    X_validation = tuning_data$X_validation,
    y_validation = tuning_data$y_validation,
    nlambda = configuration$benchmark_nlambda,
    lambda_min_ratio = configuration$benchmark_lambda_min_ratio,
    tolerance = configuration$grpreg_tolerance,
    max_iterations = configuration$grpreg_max_iterations
  )
  context <- data.frame(
    task_id = task$task_id,
    scenario_index = task$scenario_index,
    scenario = task$scenario,
    replication = task$replication,
    seed = task$seed,
    key = task$key,
    train_sample_size = nrow(tuning_data$X_train),
    validation_sample_size = nrow(tuning_data$X_validation),
    predictor_count = ncol(tuning_data$X_train),
    group_count = length(unique(tuning_data$group)),
    train_event_count = sum(tuning_data$y_train == 1),
    validation_event_count = sum(tuning_data$y_validation == 1),
    train_observed_prevalence = mean(tuning_data$y_train),
    validation_observed_prevalence = mean(tuning_data$y_validation),
    stringsAsFactors = FALSE
  )
  diagnostic$task <- task
  diagnostic$controls <- list(
    penalty = "grLasso",
    family = "binomial",
    alpha = 1,
    gamma_argument = 3,
    nlambda = configuration$benchmark_nlambda,
    lambda_min_ratio = configuration$benchmark_lambda_min_ratio,
    log_lambda = TRUE,
    tolerance = configuration$grpreg_tolerance,
    max_iterations = configuration$grpreg_max_iterations,
    dfmax = ncol(tuning_data$X_train),
    gmax = length(unique(tuning_data$group)),
    returnX = FALSE
  )
  diagnostic$data_sha256 <- data_sha256
  diagnostic$data_summary <- context
  diagnostic$summary <- cbind(context, diagnostic$summary)
  if (nrow(diagnostic$points)) {
    diagnostic$points <- cbind(
      context[rep(1L, nrow(diagnostic$points)), , drop = FALSE],
      diagnostic$points
    )
    rownames(diagnostic$points) <- NULL
  }
  diagnostic
}


lsg_run_grpreg_tasks_v8 <- function(
    root,
    task_ids = c(126L, 134L, 150L)
) {
  frozen <- lsg_grpreg_tasks_v8(root, task_ids)
  records <- lapply(frozen$tasks$task_id, function(task_id) {
    lsg_run_grpreg_task_v8(root, task_id)
  })
  names(records) <- sprintf("task_%04d", frozen$tasks$task_id)
  records
}


lsg_grpreg_diagnostic_unit_checks_v8 <- function(root) {
  frozen <- lsg_grpreg_tasks_v8(root)
  observed_cases <- frozen$tasks[c("task_id", "scenario", "replication", "seed")]
  rownames(observed_cases) <- NULL
  prefix <- lsg_grpreg_usable_prefix_v8(
    iterations = c(10, 20, 70),
    max_iterations = 100,
    finite_coefficient = rep(TRUE, 3L),
    finite_probability = rep(TRUE, 3L),
    finite_validation_loss = rep(TRUE, 3L),
    warnings = "Algorithm failed to converge for all values of lambda"
  )
  clean <- lsg_grpreg_usable_prefix_v8(
    iterations = c(1, 2, 3),
    max_iterations = 100,
    finite_coefficient = rep(TRUE, 3L),
    finite_probability = rep(TRUE, 3L),
    finite_validation_loss = rep(TRUE, 3L)
  )
  checks <- c(
    exact_three_frozen_tasks = identical(
      as.integer(frozen$tasks$task_id), c(126L, 134L, 150L)
    ),
    frozen_task_keys_and_seeds = identical(
      observed_cases, lsg_grpreg_diagnostic_cases_v8()
    ),
    v7_grpreg_controls_reused_exactly =
      frozen$configuration$benchmark_nlambda == 30L &&
      abs(frozen$configuration$benchmark_lambda_min_ratio - 0.05) < 1e-15 &&
      abs(frozen$configuration$grpreg_tolerance - 1e-7) < 1e-15 &&
      frozen$configuration$grpreg_max_iterations == 1000000L,
    budget_endpoint_excluded_from_candidate_prefix =
      identical(prefix$member, c(TRUE, TRUE, FALSE)) && prefix$length == 2L,
    clean_complete_prefix_retained =
      identical(clean$member, rep(TRUE, 3L)) && clean$length == 3L,
    no_test_metric_interface = !any(grepl(
      "X_test|y_test|test_", names(formals(lsg_fit_raw_grpreg_group_lasso_v8))
    ))
  )
  data.frame(check = names(checks), passed = unname(checks),
             stringsAsFactors = FALSE)
}


lsg_grpreg_diagnostic_smoke_v8 <- function(seed = 810260901L) {
  set.seed(as.integer(seed))
  n_train <- 60L
  n_validation <- 40L
  groups <- 4L
  group_size <- 3L
  p <- groups * group_size
  X <- matrix(stats::rnorm((n_train + n_validation) * p),
              nrow = n_train + n_validation)
  group <- rep(seq_len(groups), each = group_size)
  beta <- c(rep(0.6, group_size), rep(0, p - group_size))
  probability <- stats::plogis(-0.4 + drop(X %*% beta))
  y <- stats::rbinom(length(probability), 1, probability)
  diagnostic <- lsg_fit_raw_grpreg_group_lasso_v8(
    X = X[seq_len(n_train), , drop = FALSE],
    y = y[seq_len(n_train)],
    group = group,
    X_validation = X[n_train + seq_len(n_validation), , drop = FALSE],
    y_validation = y[n_train + seq_len(n_validation)],
    nlambda = 5L,
    lambda_min_ratio = 0.2,
    tolerance = 1e-6,
    max_iterations = 10000L
  )
  checks <- c(
    raw_fit_returned = is.list(diagnostic$raw$fit),
    exact_requested_path_recorded =
      diagnostic$summary$requested_path_length == 5L,
    point_rows_match_returned_path =
      nrow(diagnostic$points) == diagnostic$summary$returned_path_length,
    iteration_vector_preserved = identical(
      as.numeric(diagnostic$points$iterations),
      as.numeric(diagnostic$raw$iterations)
    ),
    finite_validation_loss_recomputed = all(
      is.finite(diagnostic$points$validation_log_loss[
        diagnostic$points$finite_validation_loss
      ])
    ),
    selector_not_called = !"selection" %in% names(diagnostic)
  )
  list(
    checks = data.frame(check = names(checks), passed = unname(checks),
                        stringsAsFactors = FALSE),
    diagnostic = diagnostic
  )
}
