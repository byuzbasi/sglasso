# Joint finite-path / analytic-limit tuning for the independent v7 study.
# Frozen v4/v5/v6 utilities must be sourced first. No frozen binding, option,
# package API or C++ solver is modified by this module.
#
# Selection uses validation log-loss only. Exact ties are resolved by alpha
# ascending, d ascending, analytic limit before finite within the same pair,
# then finite lambda descending. Near-ties are NOT treated as exact ties.
# This is a search of the stated candidate set, not a continuum optimum.

lsg_lambda_relative_grid_v7 <- function() {
  c(512, 256, 128, 64, 32, 16,
    exp(seq(log(8), log(0.05), length.out = 30L)))
}

lsg_finite_fitters_v7 <- function(configuration) {
  expected <- lsg_lambda_relative_grid_v7()
  if (!isTRUE(all.equal(as.numeric(configuration$lambda_relative_grid),
                        expected, tolerance = 1e-13)) ||
      !identical(as.integer(configuration$nlambda), length(expected)) ||
      !isTRUE(all.equal(as.numeric(configuration$lambda_upper_multiplier), 512)) ||
      !isTRUE(all.equal(as.numeric(configuration$lambda_min_ratio), 0.05))) {
    stop("v7 requires the frozen 36-point augmented finite lambda grid.",
         call. = FALSE)
  }
  clone <- function(original) {
    child <- new.env(parent = environment(original))
    child$logistic_sglasso_lambda_relative_grid_v5 <- local({
      frozen_grid <- expected
      function() frozen_grid
    })
    result <- original
    environment(result) <- child
    result
  }
  list(sglasso = clone(fit_extended_sglasso_validation_grid_v5),
       adelie = clone(fit_adelie_validation_grid_v5))
}

lsg_bind_tuning_rows_v7 <- function(...) {
  frames <- list(...)
  columns <- unique(unlist(lapply(frames, names), use.names = FALSE))
  frames <- lapply(frames, function(frame) {
    for (name in setdiff(columns, names(frame))) frame[[name]] <- NA
    frame[, columns, drop = FALSE]
  })
  result <- do.call(rbind, frames)
  rownames(result) <- NULL
  result
}

lsg_candidate_order_v7 <- function(tuning) {
  d <- if ("d" %in% names(tuning)) tuning$d else rep(0, nrow(tuning))
  d[is.na(d)] <- 0
  order(tuning$alpha, d, tuning$point_type != "penalty_limit",
        -tuning$lambda_relative_to_reference, tuning$lambda_index)
}

lsg_task_context_v7 <- function(frame, scenario, replication, seed) {
  frame$scenario <- rep(scenario$scenario[[1L]], nrow(frame))
  frame$replication <- rep(as.integer(replication), nrow(frame))
  frame$seed <- rep(as.integer(seed), nrow(frame))
  for (name in setdiff(names(scenario), names(frame))) {
    frame[[name]] <- rep(scenario[[name]][[1L]], nrow(frame))
  }
  frame$p <- rep(scenario$groups[[1L]] * scenario$group_size[[1L]], nrow(frame))
  frame
}

lsg_select_sglasso_candidates_v7 <- function(tuning, allowed_d) {
  allowed <- vapply(tuning$d, function(d) any(abs(d - allowed_d) <= 1e-12),
                    logical(1))
  candidates <- tuning[allowed, , drop = FALSE]
  candidates <- candidates[lsg_candidate_order_v7(candidates), , drop = FALSE]
  valid <- candidates[candidates$numerically_eligible %in% TRUE, , drop = FALSE]
  if (!nrow(valid)) {
    stop("No numerically eligible v7 SGLASSO candidate.", call. = FALSE)
  }
  selected <- valid[which.min(valid$validation_log_loss), , drop = FALSE]
  invalid_loss <- candidates$validation_log_loss[
    !candidates$numerically_eligible & is.finite(candidates$validation_log_loss)]
  best_invalid <- if (length(invalid_loss)) min(invalid_loss) else NA_real_
  finite <- candidates[candidates$point_type == "finite", , drop = FALSE]
  finite_path_keys <- paste(finite$alpha_index, finite$d_index, sep = ":")
  finite_coverage <- vapply(split(finite$numerically_eligible, finite_path_keys),
                           function(x) any(x %in% TRUE), logical(1))
  list(
    row = selected,
    validation_log_loss = selected$validation_log_loss[[1L]],
    total_candidate_points = nrow(candidates),
    numerically_eligible_points = nrow(valid),
    numerically_ineligible_points = nrow(candidates) - nrow(valid),
    nonfinite_validation_loss_points = sum(!is.finite(candidates$validation_log_loss)),
    best_invalid_validation_log_loss = best_invalid,
    invalid_candidate_competitive = length(invalid_loss) > 0L &&
      best_invalid <= selected$validation_log_loss[[1L]],
    all_finite_paths_have_valid_candidates = length(finite_coverage) > 0L &&
      all(finite_coverage),
    all_alpha_have_valid_candidates = setequal(unique(valid$alpha), unique(candidates$alpha)),
    all_d_have_valid_candidates = all(vapply(allowed_d, function(d) {
      any(abs(valid$d - d) <= 1e-12)
    }, logical(1)))
  )
}

lsg_selection_boundary_v7 <- function(selected, configuration) {
  finite <- identical(selected$point_type[[1L]], "finite")
  reference_type <- selected$lambda_reference_type[[1L]]
  extended <- reference_type %in% c("null_score_ridge_boundary", "d0_null_kkt")
  ratio <- selected$lambda_relative_to_reference[[1L]]
  upper <- finite && extended && is.finite(ratio) && abs(ratio - 512) <= 1e-10
  lower <- finite && extended && is.finite(ratio) &&
    abs(ratio - configuration$lambda_min_ratio) <= 1e-10
  # Native positive-alpha starts are their package's null-model reference,
  # not the finite 512 upper-range stress point used for shifted/ridge paths.
  list(upper = upper, lower = lower,
       native_upper = finite && !extended && is.finite(ratio) &&
         abs(ratio - 1) <= 1e-10,
       native_lower = finite && !extended && is.finite(ratio) &&
         abs(ratio - configuration$lambda_min_ratio) <= 1e-10,
       unresolved = upper || lower)
}

lsg_fit_sglasso_joint_v7 <- function(data, target_original, configuration, fitters) {
  lsg_tail_validate_data_v6(data)
  finite <- fitters$sglasso(
    data$X_train, data$y_train, data$group, data$X_validation,
    data$y_validation, target_original, configuration)
  tuning <- flatten_sglasso_validation_grid_v1(
    finite, configuration$d_grid, configuration$full_path_kkt_limit)
  tuning$point_type <- "finite"
  tuning$lambda_reference <- finite$lambda_reference[tuning$alpha_index]
  tuning$lambda_reference_type <- finite$lambda_reference_type[tuning$alpha_index]
  tuning$lambda_relative_to_reference <- tuning$lambda / tuning$lambda_reference
  tuning$lambda_relative_to_d0_kkt <-
    tuning$lambda / finite$lambda_d0_kkt_reference[tuning$alpha_index]
  tuning$finite_coefficient <- vapply(seq_len(nrow(tuning)), function(i) {
    row <- tuning[i, ]
    all(is.finite(finite$fits[[row$alpha_index]]$coefficients[
      , row$lambda_index, row$d_index]))
  }, logical(1))
  tuning$penalty_kkt <- NA_real_
  tuning$intercept_score <- vapply(seq_len(nrow(tuning)), function(i) {
    row <- tuning[i, ]
    finite$fits[[row$alpha_index]]$intercept_kkt[row$lambda_index, row$d_index]
  }, numeric(1))
  tuning$reconstruction_error <- NA_real_
  tuning$maximum_raw_objective_increase <- vapply(seq_len(nrow(tuning)), function(i) {
    row <- tuning[i, ]
    finite$fits[[row$alpha_index]]$maximum_raw_objective_increase[
      row$lambda_index, row$d_index]
  }, numeric(1))
  tuning$numerically_eligible <- tuning$numerically_eligible & tuning$finite_coefficient
  preprocess <- finite$fits[[1L]]$preprocess
  limits <- list()
  limit_rows <- list()
  limit_started <- proc.time()[["elapsed"]]
  for (ai in seq_along(finite$fits)) {
    for (di in seq_along(finite$fits[[ai]]$d)) {
      alpha <- finite$fits[[ai]]$alpha
      d <- finite$fits[[ai]]$d[di]
      endpoint <- lsg_tail_limit_fit_v6(
        data, preprocess, target_original, alpha, d, configuration)
      key <- paste(ai, di, sep = ":")
      limits[[key]] <- endpoint
      limit_rows[[key]] <- data.frame(
        alpha_index = ai, alpha = alpha, lambda_index = 0L, d_index = di,
        lambda = Inf, lambda_fraction = Inf, d = d,
        validation_log_loss = endpoint$point$validation_log_loss,
        objective = NA_real_, kkt = NA_real_,
        converged = endpoint$point$numerically_eligible, passes = NA_real_,
        selected_groups = endpoint$point$selected_groups,
        d_allowed = TRUE,
        finite_validation_loss = is.finite(endpoint$point$validation_log_loss),
        numerically_eligible = endpoint$point$numerically_eligible,
        point_type = "penalty_limit",
        lambda_reference = finite$lambda_reference[ai],
        lambda_reference_type = finite$lambda_reference_type[ai],
        lambda_relative_to_reference = Inf,
        lambda_relative_to_d0_kkt = if (alpha == 0) NA_real_ else Inf,
        finite_coefficient = all(is.finite(endpoint$coefficients)),
        penalty_kkt = endpoint$point$penalty_kkt,
        intercept_score = endpoint$point$intercept_score,
        reconstruction_error = endpoint$point$reconstruction_error,
        maximum_raw_objective_increase = NA_real_, stringsAsFactors = FALSE)
    }
  }
  limit_seconds <- proc.time()[["elapsed"]] - limit_started
  tuning <- lsg_bind_tuning_rows_v7(tuning, do.call(rbind, limit_rows))
  tuning <- tuning[lsg_candidate_order_v7(tuning), , drop = FALSE]
  rownames(tuning) <- NULL
  tuning$candidate_id <- paste("sg", tuning$alpha_index, tuning$d_index,
                               tuning$lambda_index, sep = ":")
  tuning$alpha_zero_ridge_boundary <- abs(tuning$alpha) <= 1e-12
  tuning$group_selection_capable <- !tuning$alpha_zero_ridge_boundary
  free <- lsg_select_sglasso_candidates_v7(tuning, configuration$d_grid)
  d0 <- lsg_select_sglasso_candidates_v7(tuning, 0)
  tuning$selected <- tuning$candidate_id == free$row$candidate_id[[1L]]
  tuning$selected_free_d <- tuning$selected
  tuning$selected_d0_boundary <- tuning$candidate_id == d0$row$candidate_id[[1L]]
  tuning$validation_loss_minus_selected <- tuning$validation_log_loss - free$validation_log_loss
  tuning$invalid_validation_contender <- !tuning$numerically_eligible &
    tuning$finite_validation_loss & tuning$validation_log_loss <= free$validation_log_loss
  list(finite = finite, limits = limits, tuning = tuning,
       free_selection = free, d0_selection = d0,
       limit_runtime_seconds = limit_seconds,
       runtime_seconds = finite$elapsed_seconds + limit_seconds)
}

lsg_fit_external_joint_v7 <- function(data, configuration, fitters) {
  lsg_tail_validate_data_v6(data)
  adelie <- fitters$adelie(
    data$X_train, data$y_train, data$group, data$X_validation, data$y_validation,
    alpha_grid = configuration$benchmark_alpha_grid,
    nlambda = configuration$benchmark_nlambda,
    lambda_min_ratio = configuration$benchmark_lambda_min_ratio,
    ridge_boundary_upper_multiplier = configuration$lambda_upper_multiplier,
    tolerance = configuration$adelie_tolerance,
    max_iterations = configuration$adelie_max_iterations,
    irls_tolerance = configuration$adelie_irls_tolerance,
    irls_max_iterations = configuration$adelie_irls_max_iterations)
  adelie$tuning$point_type <- "finite"
  adelie$tuning$d <- NA_real_
  adelie$tuning$finite_coefficient <- vapply(seq_len(nrow(adelie$tuning)), function(i) {
    row <- adelie$tuning[i, ]
    path <- adelie$coefficient_paths[[row$fit_index]]
    all(is.finite(c(path$intercept[row$lambda_index], path$beta[row$lambda_index, ])))
  }, logical(1))
  adelie$tuning$penalty_kkt <- NA_real_
  adelie$tuning$intercept_score <- NA_real_
  adelie$tuning$reconstruction_error <- NA_real_
  adelie$finite_tuning <- adelie$tuning
  zero_index <- which(abs(adelie$alpha_grid) <= 1e-12)
  if (length(zero_index) != 1L) stop("Adelie alpha grid must contain zero once.", call. = FALSE)
  started <- proc.time()[["elapsed"]]
  endpoint <- lsg_tail_limit_fit_v6(
    data, adelie$preprocess, numeric(ncol(data$X_train)), 0, 0, configuration)
  endpoint_seconds <- proc.time()[["elapsed"]] - started
  # Use a distinct synthetic fit ID: the finite fit's returned path metadata
  # remains untouched and its original numerical audit is still required.
  endpoint_row <- adelie$tuning[which(adelie$tuning$alpha_index == zero_index)[1L], ]
  endpoint_row$fit_index <- length(adelie$alpha_grid) + 1L
  endpoint_row$lambda_index <- 0L
  endpoint_row$lambda <- endpoint_row$lambda_fraction <- Inf
  endpoint_row$lambda_relative_to_reference <- Inf
  endpoint_row$validation_log_loss <- endpoint$point$validation_log_loss
  endpoint_row$solver_warning <- endpoint_row$unexpected_solver_warning <- ""
  endpoint_row$requested_path_length <- endpoint_row$returned_path_length <- 1L
  endpoint_row$path_complete <- endpoint_row$path_termination_acceptable <- TRUE
  endpoint_row$point_converged <- endpoint_row$numerically_eligible <- endpoint$point$numerically_eligible
  endpoint_row$finite_validation_loss <- is.finite(endpoint$point$validation_log_loss)
  endpoint_row$finite_coefficient <- all(is.finite(endpoint$coefficients))
  endpoint_row$selected <- FALSE
  endpoint_row$point_type <- "penalty_limit"
  endpoint_row$penalty_kkt <- endpoint$point$penalty_kkt
  endpoint_row$intercept_score <- endpoint$point$intercept_score
  endpoint_row$reconstruction_error <- endpoint$point$reconstruction_error
  endpoint_row$invalid_validation_contender <- FALSE
  endpoint_row$candidate_nonfinite_validation_loss <- FALSE
  endpoint_row$policy_excluded_nonfinite_validation_loss <- FALSE
  adelie$tuning <- rbind(adelie$tuning, endpoint_row)
  adelie$tuning <- adelie$tuning[lsg_candidate_order_v7(adelie$tuning), , drop = FALSE]
  adelie$tuning$selected <- FALSE
  adelie$tuning$invalid_validation_contender <- FALSE
  adelie$selection <- select_valid_external_grid_v2(
    adelie$tuning, "Logistic Group Elastic Net (adelie)")
  adelie$tuning <- mark_external_selection_v2(adelie$tuning, adelie$selection)
  adelie$endpoint <- endpoint
  adelie$endpoint_runtime_seconds <- endpoint_seconds
  adelie$runtime_seconds <- adelie$runtime_seconds + endpoint_seconds
  grpreg <- fit_grpreg_validation_paths_v2(
    data$X_train, data$y_train, data$group, data$X_validation, data$y_validation,
    nlambda = configuration$benchmark_nlambda,
    lambda_min_ratio = configuration$benchmark_lambda_min_ratio,
    tolerance = configuration$grpreg_tolerance,
    max_iterations = configuration$grpreg_max_iterations,
    penalty_specification = configuration$grpreg_penalty_specification)
  grpreg$tuning$point_type <- "finite"
  grpreg$tuning$d <- NA_real_
  grpreg$tuning$finite_coefficient <- vapply(seq_len(nrow(grpreg$tuning)), function(i) {
    row <- grpreg$tuning[i, ]
    all(is.finite(grpreg$fits[[row$fit_index]]$beta[, row$lambda_index]))
  }, logical(1))
  grpreg$tuning$penalty_kkt <- NA_real_
  grpreg$tuning$intercept_score <- NA_real_
  grpreg$tuning$reconstruction_error <- NA_real_
  tuning <- lsg_bind_tuning_rows_v7(adelie$tuning, grpreg$tuning)
  tuning$candidate_id <- paste(tuning$engine, tuning$method_path,
                               tuning$fit_index, tuning$lambda_index, sep = ":")
  list(adelie = adelie, grpreg = grpreg, tuning = tuning)
}

lsg_selected_sglasso_model_v7 <- function(grid, selection, data) {
  row <- selection$row
  ai <- row$alpha_index[[1L]]
  di <- row$d_index[[1L]]
  if (identical(row$point_type[[1L]], "penalty_limit")) {
    endpoint <- grid$limits[[paste(ai, di, sep = ":")]]
    coefficient <- as.numeric(endpoint$coefficients)
    test_eta <- endpoint$intercept_solver + drop(
      transform_lsg_newx(grid$finite$fits[[ai]]$preprocess, data$X_test) %*%
        endpoint$beta_solver)
    test_probability <- stats::plogis(test_eta)
    validation_probability <- endpoint$validation_probability
  } else {
    fit <- grid$finite$fits[[ai]]
    li <- row$lambda_index[[1L]]
    coefficient <- as.numeric(fit$coefficients[, li, di])
    test_probability <- drop(predict_logistic_sglasso(
      fit, data$X_test, type = "response", lambda_index = li, d_index = di))
    validation_probability <- drop(predict_logistic_sglasso(
      fit, data$X_validation, type = "response", lambda_index = li, d_index = di))
  }
  list(coefficient = coefficient, probability = test_probability,
       validation_probability = validation_probability,
       alpha = row$alpha[[1L]], d = row$d[[1L]], gamma = NA_real_,
       lambda = row$lambda[[1L]],
       lambda_relative = row$lambda_relative_to_reference[[1L]],
       point_type = row$point_type[[1L]], candidate_id = row$candidate_id[[1L]],
       numerical = as.list(row[c("converged", "numerically_eligible", "kkt",
         "penalty_kkt", "intercept_score", "reconstruction_error")]))
}

lsg_selected_external_model_v7 <- function(grid, method_path, data) {
  row <- grid$tuning[grid$tuning$method_path == method_path & grid$tuning$selected, ]
  if (nrow(row) != 1L) stop("Exactly one selected external candidate is required.", call. = FALSE)
  if (identical(row$engine[[1L]], "adelie")) {
    if (identical(row$point_type[[1L]], "penalty_limit")) {
      endpoint <- grid$adelie$endpoint
      coefficient <- as.numeric(endpoint$coefficients)
      probability <- stats::plogis(endpoint$intercept_solver + drop(
        transform_lsg_newx(grid$adelie$preprocess, data$X_test) %*% endpoint$beta_solver))
      validation_probability <- endpoint$validation_probability
    } else {
      selected_grid <- grid$adelie
      selected_grid$selection$row <- row
      test <- extract_adelie_group_en_solution_v2(selected_grid, data$X_test)
      validation <- extract_adelie_group_en_solution_v2(selected_grid, data$X_validation)
      coefficient <- test$coefficient
      probability <- test$probability
      validation_probability <- validation$probability
    }
  } else {
    penalty <- grid$grpreg$penalty_specification$penalty[
      grid$grpreg$penalty_specification$method_path == method_path]
    test <- extract_grpreg_group_solution_v2(grid$grpreg, penalty, data$X_test)
    validation <- extract_grpreg_group_solution_v2(grid$grpreg, penalty, data$X_validation)
    coefficient <- test$coefficient
    probability <- test$probability
    validation_probability <- validation$probability
  }
  list(coefficient = as.numeric(coefficient), probability = as.numeric(probability),
       validation_probability = as.numeric(validation_probability),
       alpha = row$alpha[[1L]], d = NA_real_, gamma = row$gamma[[1L]],
       lambda = row$lambda[[1L]],
       lambda_relative = row$lambda_relative_to_reference[[1L]],
       point_type = row$point_type[[1L]], candidate_id = row$candidate_id[[1L]],
       numerical = as.list(row[c("point_converged", "numerically_eligible",
         "path_termination_acceptable", "penalty_kkt", "intercept_score",
         "reconstruction_error", "unexpected_solver_warning")]))
}

lsg_check_selected_model_v7 <- function(model, data, configuration) {
  if (length(model$coefficient) != ncol(data$X_test) + 1L ||
      any(!is.finite(model$coefficient)) ||
      length(model$probability) != nrow(data$X_test) ||
      length(model$validation_probability) != nrow(data$X_validation) ||
      any(!is.finite(c(model$probability, model$validation_probability))) ||
      any(c(model$probability, model$validation_probability) < 0 |
            c(model$probability, model$validation_probability) > 1)) {
    stop("Invalid selected-model coefficient or probability arrays.", call. = FALSE)
  }
  model$eta <- model$coefficient[1L] + drop(data$X_test %*% model$coefficient[-1L])
  model$validation_eta <- model$coefficient[1L] +
    drop(data$X_validation %*% model$coefficient[-1L])
  model$prediction_reconstruction_error <- max(
    abs(stats::plogis(model$eta) - model$probability),
    abs(stats::plogis(model$validation_eta) - model$validation_probability))
  if (!is.finite(model$prediction_reconstruction_error) ||
      model$prediction_reconstruction_error > configuration$reconstruction_tolerance) {
    stop("Selected original-scale coefficients do not reconstruct predictions.", call. = FALSE)
  }
  model$validation_log_loss <- binary_log_loss(data$y_validation, model$validation_probability)
  model
}

lsg_run_task_v7 <- function(scenario, replication, seed, configuration) {
  if (!is.data.frame(scenario) || nrow(scenario) != 1L ||
      length(seed) != 1L || !is.finite(seed) || seed != as.integer(seed)) {
    stop("One design row and one frozen integer seed are required.", call. = FALSE)
  }
  fitters <- lsg_finite_fitters_v7(configuration)
  data <- simulate_logistic_sglasso_two_design_v1(scenario, as.integer(seed))
  # No test/truth field is permitted across the tuning interface.
  tuning_data <- data[c("X_train", "y_train", "X_validation", "y_validation", "group")]
  lsg_tail_validate_data_v6(tuning_data)
  firth <- estimate_groupwise_logistic_target(
    tuning_data$X_train, tuning_data$y_train, tuning_data$group,
    method = "firth", max_iterations = configuration$target_max_iterations,
    tolerance = configuration$target_tolerance)
  if (!isTRUE(firth$success) || firth$failed_groups != 0L ||
      any(!is.finite(firth$target_original))) {
    stop("Firth training-only target failed; no fallback or redraw is allowed.", call. = FALSE)
  }
  sglasso <- lsg_fit_sglasso_joint_v7(tuning_data, firth$target_original,
                                    configuration, fitters)
  external <- lsg_fit_external_joint_v7(tuning_data, configuration, fitters)
  method_names <- c("Logistic SGLASSO", "Logistic SGLASSO (d=0 boundary)",
                    "Logistic Group Elastic Net (adelie)",
                    "Logistic Group Lasso (grpreg)", "Logistic Group MCP (grpreg)",
                    "Logistic Group SCAD (grpreg)")
  models <- list(
    lsg_selected_sglasso_model_v7(sglasso, sglasso$free_selection, data),
    lsg_selected_sglasso_model_v7(sglasso, sglasso$d0_selection, data))
  for (method in method_names[-c(1L, 2L)]) {
    models[[length(models) + 1L]] <- lsg_selected_external_model_v7(external, method, data)
  }
  names(models) <- method_names
  selected <- list(sglasso$free_selection$row, sglasso$d0_selection$row)
  for (method in method_names[-c(1L, 2L)]) {
    selected[[length(selected) + 1L]] <- external$tuning[
      external$tuning$method_path == method & external$tuning$selected, , drop = FALSE]
  }
  names(selected) <- method_names
  true_coefficient <- c(data$intercept, data$beta)
  true_test_eta <- data$intercept + drop(data$X_test %*% data$beta)
  event_summary <- data.frame(
    train_event_count = sum(data$y_train == 1),
    validation_event_count = sum(data$y_validation == 1),
    test_event_count = sum(data$y_test == 1),
    train_observed_prevalence = mean(data$y_train),
    validation_observed_prevalence = mean(data$y_validation),
    test_observed_prevalence = mean(data$y_test))
  results <- vector("list", length(models))
  threshold <- configuration$selection_threshold
  if (is.null(threshold)) threshold <- 1e-8
  for (i in seq_along(models)) {
    model <- lsg_check_selected_model_v7(models[[i]], data, configuration)
    row <- selected[[i]]
    if (abs(model$validation_log_loss - row$validation_log_loss[[1L]]) > 1e-10) {
      stop("Stored selected probabilities do not reproduce validation loss.", call. = FALSE)
    }
    boundary <- lsg_selection_boundary_v7(row, configuration)
    metric <- lsg_metrics_v7(
      y = data$y_test, probability = model$probability,
      coefficient = model$coefficient, true_coefficient = true_coefficient,
      group = data$group, eta = model$eta, true_eta = true_test_eta,
      threshold = threshold)
    metadata <- data.frame(
      method = method_names[i], selected_alpha = model$alpha, selected_d = model$d,
      selected_gamma = model$gamma, selected_lambda = model$lambda,
      selected_lambda_reference = row$lambda_reference[[1L]],
      selected_lambda_relative_to_reference = model$lambda_relative,
      selected_lambda_reference_type = row$lambda_reference_type[[1L]],
      selected_point_type = model$point_type, selected_candidate_id = model$candidate_id,
      validation_log_loss = model$validation_log_loss,
      selected_converged = row$numerically_eligible[[1L]],
      selected_kkt = if ("kkt" %in% names(row)) row$kkt[[1L]] else NA_real_,
      selected_penalty_kkt = row$penalty_kkt[[1L]],
      selected_intercept_score = row$intercept_score[[1L]],
      selected_finite_upper_boundary = boundary$upper,
      selected_finite_lower_boundary = boundary$lower,
      selected_native_upper_boundary = boundary$native_upper,
      selected_native_lower_boundary = boundary$native_lower,
      selected_alpha_zero_ridge_boundary = abs(model$alpha) <= 1e-12,
      group_selection_capable = abs(model$alpha) > 1e-12,
      selection_unit = "predefined_group", selection_threshold = threshold,
      prediction_reconstruction_error = model$prediction_reconstruction_error,
      stringsAsFactors = FALSE)
    results[[i]] <- cbind(metadata, event_summary,
      metric[, setdiff(names(metric), c(names(metadata), names(event_summary))), drop = FALSE])
    model$numerical$selected_finite_upper_boundary <- boundary$upper
    model$numerical$selected_finite_lower_boundary <- boundary$lower
    models[[i]] <- model
  }
  external_audits <- list(external_method_audit_row_v2(
    external$adelie$tuning, external$adelie$selection, scenario, replication,
    seed, external$adelie$runtime_seconds))
  for (penalty in names(external$grpreg$selections)) {
    external_audits[[length(external_audits) + 1L]] <- external_method_audit_row_v2(
      external$grpreg$tuning, external$grpreg$selections[[penalty]], scenario,
      replication, seed, external$grpreg$elapsed_by_penalty[[penalty]])
  }
  diagnostics <- cbind(event_summary, data.frame(
    firth_failed_groups = firth$failed_groups,
    firth_separation_groups = firth$separation_groups,
    firth_runtime_seconds = firth$elapsed_seconds,
    sglasso_finite_runtime_seconds = sglasso$finite$elapsed_seconds,
    sglasso_limit_runtime_seconds = sglasso$limit_runtime_seconds,
    sglasso_runtime_seconds = sglasso$runtime_seconds,
    adelie_runtime_seconds = external$adelie$runtime_seconds,
    grpreg_runtime_seconds = external$grpreg$runtime_seconds,
    sglasso_free_minus_d0_validation_log_loss =
      sglasso$free_selection$validation_log_loss - sglasso$d0_selection$validation_log_loss,
    sglasso_invalid_candidate_competitive = sglasso$free_selection$invalid_candidate_competitive,
    sglasso_nonfinite_validation_loss_points = sglasso$free_selection$nonfinite_validation_loss_points,
    sglasso_all_alpha_have_valid_candidates = sglasso$free_selection$all_alpha_have_valid_candidates,
    sglasso_all_d_have_valid_candidates = sglasso$free_selection$all_d_have_valid_candidates,
    sglasso_all_finite_paths_have_valid_candidates =
      sglasso$free_selection$all_finite_paths_have_valid_candidates,
    maximum_prediction_reconstruction_error = max(vapply(
      models, function(x) x$prediction_reconstruction_error, numeric(1))),
    stringsAsFactors = FALSE))
  fingerprint <- lsg_data_fingerprint_v7
  artifacts <- list(
    observation_id = seq_along(data$y_test), y_test = data$y_test,
    y_train = data$y_train,
    true_test_probability = data$true_test_probability,
    true_coefficient = true_coefficient, group = data$group,
    active_groups = data$active_groups, true_test_eta = true_test_eta,
    validation_observation_id = seq_along(data$y_validation), y_validation = data$y_validation,
    models = models, firth_target_original = firth$target_original,
    training_sample_size = nrow(data$X_train),
    data_sha256 = list(
      training = fingerprint(list(X = data$X_train, y = data$y_train, group = data$group)),
      validation = fingerprint(list(X = data$X_validation, y = data$y_validation)),
      test = fingerprint(list(X = data$X_test, y = data$y_test)),
      truth = fingerprint(list(coefficient = true_coefficient, group = data$group,
                               active_groups = data$active_groups))))
  payload <- list(
    results = lsg_task_context_v7(do.call(rbind, results), scenario, replication, seed),
    sglasso_tuning = lsg_task_context_v7(sglasso$tuning, scenario, replication, seed),
    external_tuning = lsg_task_context_v7(external$tuning, scenario, replication, seed),
    diagnostics = lsg_task_context_v7(diagnostics, scenario, replication, seed),
    external_diagnostics = lsg_task_context_v7(do.call(rbind, external_audits), scenario, replication, seed),
    target_diagnostics = lsg_task_context_v7(firth$diagnostics, scenario, replication, seed),
    artifacts = artifacts)
  payload$tuning_checks <- lsg_validate_tuning_payload_v7(payload, configuration)
  payload$diagnostics$all_numerical_gates_passed <-
    all(payload$tuning_checks$passed[payload$tuning_checks$gate == "numerical"])
  payload$diagnostics$boundary_gates_passed <-
    all(payload$tuning_checks$passed[payload$tuning_checks$gate == "boundary"])
  payload
}

lsg_validate_tuning_payload_v7 <- function(payload, configuration) {
  sg <- payload$sglasso_tuning
  et <- payload$external_tuning
  results <- payload$results
  required_sg <- c(
    "candidate_id", "alpha_index", "d_index", "alpha", "d", "lambda_index",
    "lambda", "lambda_reference", "lambda_reference_type",
    "lambda_relative_to_reference", "point_type", "validation_log_loss",
    "finite_validation_loss", "finite_coefficient", "kkt", "penalty_kkt",
    "intercept_score", "reconstruction_error", "converged", "numerically_eligible",
    "selected", "selected_free_d", "selected_d0_boundary", "invalid_validation_contender")
  required_et <- c(
    "candidate_id", "engine", "method_path", "alpha", "fit_index",
    "lambda_index", "lambda", "lambda_reference", "lambda_reference_type",
    "lambda_relative_to_reference", "point_type", "validation_log_loss",
    "finite_validation_loss", "finite_coefficient", "penalty_kkt", "intercept_score",
    "reconstruction_error", "point_converged", "numerically_eligible", "selected",
    "path_complete", "path_termination_acceptable", "returned_lower_boundary_excluded",
    "unexpected_solver_warning", "requested_path_length", "returned_path_length")
  if (!is.data.frame(sg) || !is.data.frame(et) || !nrow(sg) || !nrow(et) ||
      !all(required_sg %in% names(sg)) || !all(required_et %in% names(et))) {
    stop("The v7 tuning payload is missing required raw candidate fields.", call. = FALSE)
  }
  checks <- list()
  add <- function(name, passed, gate = "numerical", detail = "") {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = name, passed = isTRUE(passed), gate = gate,
      detail = as.character(detail), stringsAsFactors = FALSE)
  }
  finite_sg <- sg$point_type == "finite"
  limit_sg <- sg$point_type == "penalty_limit"
  expected_grid <- lsg_lambda_relative_grid_v7()
  alpha_grid <- sort(unique(as.numeric(configuration$alpha_grid)))
  d_grid <- sort(unique(as.numeric(configuration$d_grid)))
  pairs <- do.call(rbind, lapply(seq_along(alpha_grid), function(ai) {
    ds <- if (abs(alpha_grid[ai] - 1) <= 1e-12) 0 else d_grid
    data.frame(alpha_index = ai, d_index = seq_along(ds), alpha = alpha_grid[ai], d = ds)
  }))
  shape <- nrow(sg) == nrow(pairs) * (length(expected_grid) + 1L) &&
    all(finite_sg | limit_sg) && !anyDuplicated(sg$candidate_id) &&
    !anyDuplicated(paste(sg$alpha_index, sg$d_index, sg$lambda_index, sep = ":"))
  finite_coverage <- logical(nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    pair <- pairs[i, ]
    x <- sg[sg$alpha_index == pair$alpha_index & sg$d_index == pair$d_index, ]
    f <- x[x$point_type == "finite", ]
    f <- f[order(f$lambda_index), ]
    endpoint <- x[x$point_type == "penalty_limit", ]
    shape <- shape && nrow(x) == length(expected_grid) + 1L &&
      all(abs(x$alpha - pair$alpha) <= 1e-12) && all(abs(x$d - pair$d) <= 1e-12) &&
      identical(as.integer(f$lambda_index), seq_along(expected_grid)) &&
      isTRUE(all.equal(as.numeric(f$lambda_relative_to_reference), expected_grid,
                       tolerance = 1e-13)) &&
      nrow(endpoint) == 1L && endpoint$lambda_index[[1L]] == 0L &&
      is.infinite(endpoint$lambda[[1L]]) && endpoint$lambda[[1L]] > 0 &&
      is.infinite(endpoint$lambda_relative_to_reference[[1L]]) &&
      endpoint$lambda_relative_to_reference[[1L]] > 0
    finite_coverage[i] <- nrow(f) > 0L && any(f$numerically_eligible %in% TRUE)
  }
  add("joint_sglasso_candidate_grid_exact", shape)
  reference_ok <- all(is.finite(sg$lambda_reference) & sg$lambda_reference > 0) &&
    all(sg$lambda_reference_type == ifelse(abs(sg$alpha) <= 1e-12,
                                          "null_score_ridge_boundary", "d0_null_kkt")) &&
    all(abs(sg$lambda[finite_sg] - sg$lambda_reference[finite_sg] *
              sg$lambda_relative_to_reference[finite_sg]) <=
          1e-12 * pmax(1, abs(sg$lambda[finite_sg])))
  add("sglasso_lambda_reference_metadata", reference_ok)
  sg_eligible <- rep(FALSE, nrow(sg))
  sg_eligible[finite_sg] <- sg$converged[finite_sg] %in% TRUE &
    sg$finite_coefficient[finite_sg] %in% TRUE &
    is.finite(sg$validation_log_loss[finite_sg]) & is.finite(sg$kkt[finite_sg]) &
    sg$kkt[finite_sg] <= configuration$full_path_kkt_limit
  sg_eligible[limit_sg] <- sg$finite_coefficient[limit_sg] %in% TRUE &
    is.finite(sg$validation_log_loss[limit_sg]) &
    is.finite(sg$penalty_kkt[limit_sg]) &
    sg$penalty_kkt[limit_sg] <= configuration$endpoint_tolerance &
    is.finite(sg$intercept_score[limit_sg]) &
    sg$intercept_score[limit_sg] <= configuration$endpoint_tolerance &
    is.finite(sg$reconstruction_error[limit_sg]) &
    sg$reconstruction_error[limit_sg] <= configuration$reconstruction_tolerance
  add("sglasso_eligibility_recomputed", !anyNA(sg$numerically_eligible) &&
        identical(as.logical(sg$numerically_eligible), sg_eligible))
  add("all_analytic_sglasso_limits_eligible", any(limit_sg) && all(sg_eligible[limit_sg]) &&
        all(sg$converged[limit_sg] %in% TRUE))
  add("all_finite_sglasso_paths_represented", all(finite_coverage))
  add("no_nonfinite_sglasso_validation_candidates",
      all(is.finite(sg$validation_log_loss)) && !anyNA(sg$finite_validation_loss) &&
        identical(as.logical(sg$finite_validation_loss), is.finite(sg$validation_log_loss)))
  recomputed_sg <- sg
  recomputed_sg$numerically_eligible <- sg_eligible
  free <- lsg_select_sglasso_candidates_v7(recomputed_sg, configuration$d_grid)
  d0 <- lsg_select_sglasso_candidates_v7(recomputed_sg, 0)
  free_row <- sg[sg$selected_free_d %in% TRUE, ]
  d0_row <- sg[sg$selected_d0_boundary %in% TRUE, ]
  selection_ok <- nrow(free_row) == 1L && nrow(d0_row) == 1L &&
    identical(as.logical(sg$selected), as.logical(sg$selected_free_d)) &&
    identical(free_row$candidate_id, free$row$candidate_id) &&
    identical(d0_row$candidate_id, d0$row$candidate_id)
  add("free_and_nested_d0_are_exact_validation_minimizers", selection_ok)
  contender <- !sg_eligible & is.finite(sg$validation_log_loss) &
    sg$validation_log_loss <= free$validation_log_loss
  add("no_competitive_ineligible_sglasso_candidate", !any(contender) &&
        !isTRUE(d0$invalid_candidate_competitive) &&
        identical(as.logical(sg$invalid_validation_contender), contender))
  add("sglasso_finite_range_resolved", selection_ok &&
        !lsg_selection_boundary_v7(free$row, configuration)$unresolved &&
        !lsg_selection_boundary_v7(d0$row, configuration)$unresolved, gate = "boundary")

  finite_et <- et$point_type == "finite"
  limit_et <- et$point_type == "penalty_limit"
  et_eligible <- rep(FALSE, nrow(et))
  et_eligible[finite_et] <- et$path_termination_acceptable[finite_et] %in% TRUE &
    et$point_converged[finite_et] %in% TRUE &
    et$returned_lower_boundary_excluded[finite_et] %in% FALSE &
    et$finite_coefficient[finite_et] %in% TRUE & is.finite(et$validation_log_loss[finite_et])
  et_eligible[limit_et] <- et$engine[limit_et] == "adelie" & abs(et$alpha[limit_et]) <= 1e-12 &
    et$finite_coefficient[limit_et] %in% TRUE & is.finite(et$validation_log_loss[limit_et]) &
    is.finite(et$penalty_kkt[limit_et]) & et$penalty_kkt[limit_et] <= configuration$endpoint_tolerance &
    is.finite(et$intercept_score[limit_et]) & et$intercept_score[limit_et] <= configuration$endpoint_tolerance &
    is.finite(et$reconstruction_error[limit_et]) &
    et$reconstruction_error[limit_et] <= configuration$reconstruction_tolerance
  add("external_eligibility_recomputed", all(finite_et | limit_et) &&
        !anyDuplicated(et$candidate_id) && !anyNA(et$numerically_eligible) &&
        identical(as.logical(et$numerically_eligible), et_eligible))
  add("single_eligible_adelie_null_limit", sum(limit_et) == 1L && all(et_eligible[limit_et]) &&
        all(is.infinite(et$lambda[limit_et]) & et$lambda[limit_et] > 0) &&
        all(is.infinite(et$lambda_relative_to_reference[limit_et]) &
              et$lambda_relative_to_reference[limit_et] > 0))
  finite_adelie <- et[finite_et & et$engine == "adelie", ]
  adelie_shape <- setequal(unique(finite_adelie$alpha), configuration$benchmark_alpha_grid)
  for (alpha in configuration$benchmark_alpha_grid) {
    x <- finite_adelie[abs(finite_adelie$alpha - alpha) <= 1e-12, ]
    x <- x[order(x$lambda_index), ]
    n <- if (alpha == 0) length(expected_grid) else configuration$benchmark_nlambda
    adelie_shape <- adelie_shape && nrow(x) == n &&
      length(unique(x$fit_index)) == 1L &&
      identical(as.integer(x$lambda_index), seq_len(n)) &&
      all(x$requested_path_length == n) && all(x$returned_path_length == n) &&
      (alpha != 0 || isTRUE(all.equal(as.numeric(x$lambda_relative_to_reference),
                                     expected_grid, tolerance = 1e-13)))
  }
  add("adelie_finite_grid_exact", adelie_shape)
  specification <- configuration$grpreg_penalty_specification
  grpreg_shape <- nrow(specification) == 3L
  for (i in seq_len(nrow(specification))) {
    definition <- specification[i, ]
    x <- et[et$method_path == definition$method_path, ]
    x <- x[order(x$lambda_index), ]
    gamma_ok <- if (is.na(definition$gamma)) all(is.na(x$gamma)) else {
      all(is.finite(x$gamma) & abs(x$gamma - definition$gamma) <= 1e-12)
    }
    grpreg_shape <- grpreg_shape && nrow(x) > 0L &&
      length(unique(x$fit_index)) == 1L && all(x$engine == "grpreg") &&
      all(x$point_type == "finite") && all(x$penalty_family == definition$penalty_family) &&
      all(x$alpha == 1) && gamma_ok &&
      identical(as.integer(x$lambda_index), seq_len(nrow(x))) &&
      all(x$requested_path_length == configuration$benchmark_nlambda) &&
      all(x$returned_path_length == nrow(x)) &&
      all(is.finite(x$lambda) & x$lambda > 0) &&
      all(x$lambda_reference_type == "native_grpreg_path_start")
  }
  add("grpreg_native_penalty_paths_unchanged", grpreg_shape)
  finite_keys <- paste(et$method_path[finite_et], et$fit_index[finite_et], sep = ":")
  finite_external_coverage <- vapply(split(et_eligible[finite_et], finite_keys), any, logical(1))
  add("all_external_finite_paths_represented", length(finite_external_coverage) > 0L &&
        all(finite_external_coverage))
  external_audit <- audit_external_nonfinite_validation_v3(et, require_selection_state = TRUE)
  add("frozen_external_nonfinite_audit", external_audit$audit_passes &&
        identical(as.logical(et$candidate_nonfinite_validation_loss),
                  external_audit$rows$candidate_nonfinite_validation_loss) &&
        identical(as.logical(et$policy_excluded_nonfinite_validation_loss),
                  external_audit$rows$policy_excluded_nonfinite_validation_loss))
  add("frozen_external_path_termination_rules", all(et$path_termination_acceptable %in% TRUE) &&
        all(et$path_complete[et$engine == "adelie" | et$penalty_family == "group_lasso"] %in% TRUE) &&
        all(!is.na(et$unexpected_solver_warning) & !nzchar(et$unexpected_solver_warning)))
  methods <- c("Logistic Group Elastic Net (adelie)", "Logistic Group Lasso (grpreg)",
               "Logistic Group MCP (grpreg)", "Logistic Group SCAD (grpreg)")
  external_selection_ok <- setequal(unique(et$method_path), methods)
  invalid_external <- FALSE
  truncation_ok <- TRUE
  external_minimizers <- list()
  for (method in methods) {
    method_rows <- et[et$method_path == method, ]
    method_rows <- method_rows[lsg_candidate_order_v7(method_rows), ]
    method_rows$numerically_eligible <- et_eligible[
      match(method_rows$candidate_id, et$candidate_id)]
    selected <- select_valid_external_grid_v2(method_rows, method)
    external_minimizers[[method]] <- selected$row
    recorded <- method_rows[method_rows$selected %in% TRUE, ]
    external_selection_ok <- external_selection_ok && nrow(recorded) == 1L &&
      identical(recorded$candidate_id, selected$row$candidate_id)
    invalid_external <- invalid_external || selected$invalid_candidate_competitive
    truncation_ok <- truncation_ok && selected$truncated_path_selection_interior
  }
  add("external_exact_validation_minimizers", external_selection_ok)
  add("no_competitive_ineligible_external_candidate", !invalid_external &&
        !any(et$invalid_validation_contender %in% TRUE))
  add("truncated_external_selection_interior", truncation_ok)
  add("adelie_finite_range_resolved", external_selection_ok &&
        !lsg_selection_boundary_v7(external_minimizers[[methods[1L]]], configuration)$unresolved,
      gate = "boundary")
  expected_methods <- c("Logistic SGLASSO", "Logistic SGLASSO (d=0 boundary)", methods)
  rows <- c(list(free$row, d0$row), external_minimizers)
  result_selection_ok <- is.data.frame(results) && nrow(results) == 6L &&
    !anyDuplicated(results$method) && setequal(results$method, expected_methods)
  if (result_selection_ok) for (i in seq_along(expected_methods)) {
    result <- results[results$method == expected_methods[i], ]
    row <- rows[[i]]
    boundary <- lsg_selection_boundary_v7(row, configuration)
    result_selection_ok <- result_selection_ok &&
      identical(result$selected_candidate_id, row$candidate_id) &&
      identical(result$selected_point_type, row$point_type) &&
      isTRUE(all.equal(result$selected_alpha, row$alpha, tolerance = 1e-13)) &&
      isTRUE(all.equal(result$selected_d, row$d, tolerance = 1e-13)) &&
      isTRUE(all.equal(result$selected_lambda, row$lambda, tolerance = 1e-13)) &&
      isTRUE(all.equal(result$selected_lambda_relative_to_reference,
                       row$lambda_relative_to_reference, tolerance = 1e-13)) &&
      identical(as.logical(result$selected_finite_upper_boundary), boundary$upper) &&
      identical(as.logical(result$selected_finite_lower_boundary), boundary$lower) &&
      isTRUE(result$selected_converged) &&
      abs(result$validation_log_loss - row$validation_log_loss) <= 1e-10
  }
  add("six_result_rows_match_raw_selected_candidates", result_selection_ok)
  target <- payload$target_diagnostics
  add("firth_all_groups_training_target_valid", is.data.frame(target) &&
        nrow(target) == length(unique(payload$artifacts$group)) &&
        !anyDuplicated(target$group) && all(target$method == "firth") &&
        all(target$success %in% TRUE) && all(target$converged %in% TRUE) &&
        all(target$finite_coefficients %in% TRUE))
  result <- do.call(rbind, checks)
  rownames(result) <- NULL
  result
}

lsg_tuning_unit_checks_v7 <- function() {
  check <- function(name, passed) data.frame(check = name, passed = isTRUE(passed),
                                             stringsAsFactors = FALSE)
  configuration <- list(lambda_relative_grid = lsg_lambda_relative_grid_v7(),
                        nlambda = 36L, lambda_upper_multiplier = 512,
                        lambda_min_ratio = 0.05)
  old_sg <- fit_extended_sglasso_validation_grid_v5
  old_ad <- fit_adelie_validation_grid_v5
  old_options <- options()
  fitters <- lsg_finite_fitters_v7(configuration)
  checks <- list(check("v7_cloned_finite_fitters_preserve_bodies",
                       identical(body(fitters$sglasso), body(old_sg)) &&
                         identical(body(fitters$adelie), body(old_ad))),
                 check("v7_no_global_frozen_binding_or_option_changes",
                       identical(old_options, options()) &&
                         identical(old_sg, fit_extended_sglasso_validation_grid_v5) &&
                         identical(old_ad, fit_adelie_validation_grid_v5)))
  checks[[length(checks) + 1L]] <- check("v7_isolated_36_point_grid",
    identical(environment(fitters$sglasso)$logistic_sglasso_lambda_relative_grid_v5(),
              lsg_lambda_relative_grid_v7()) &&
      identical(environment(fitters$adelie)$logistic_sglasso_lambda_relative_grid_v5(),
                lsg_lambda_relative_grid_v7()) &&
      identical(lsg_lambda_relative_grid_v7()[7:36],
                exp(seq(log(8), log(0.05), length.out = 30L))))
  toy <- data.frame(
    candidate_id = c("high_alpha", "finite", "endpoint", "higher_d"),
    alpha = c(0.5, 0, 0, 0), alpha_index = c(2L, 1L, 1L, 1L),
    d = c(0, 0, 0, 0.5), d_index = c(1L, 1L, 1L, 2L),
    point_type = c("penalty_limit", "finite", "penalty_limit", "finite"),
    lambda_index = c(0L, 1L, 0L, 1L),
    lambda_relative_to_reference = c(Inf, 512, Inf, 512),
    lambda_reference_type = c("d0_null_kkt", rep("null_score_ridge_boundary", 3L)),
    validation_log_loss = rep(0.6, 4L), numerically_eligible = TRUE,
    stringsAsFactors = FALSE)
  chosen <- lsg_select_sglasso_candidates_v7(toy, c(0, 0.5))
  checks[[length(checks) + 1L]] <- check("v7_exact_tie_rule_no_forced_positive_d",
    identical(chosen$row$candidate_id, "endpoint"))
  toy$validation_log_loss[2L] <- 0.6 - 1e-13
  chosen <- lsg_select_sglasso_candidates_v7(toy, c(0, 0.5))
  checks[[length(checks) + 1L]] <- check("v7_near_ties_not_coarsened",
    identical(chosen$row$candidate_id, "finite"))
  checks[[length(checks) + 1L]] <- check("v7_finite_upper_boundary_remains_a_gate",
    lsg_selection_boundary_v7(chosen$row, configuration)$unresolved &&
      !lsg_selection_boundary_v7(toy[3L, ], configuration)$unresolved)
  lower <- toy[2L, ]
  lower$lambda_relative_to_reference <- 0.05
  native <- lower
  native$lambda_reference_type <- "native_grpreg_path_start"
  checks[[length(checks) + 1L]] <- check("v7_finite_lower_gate_native_policy_unchanged",
    lsg_selection_boundary_v7(lower, configuration)$unresolved &&
      !lsg_selection_boundary_v7(native, configuration)$unresolved &&
      lsg_selection_boundary_v7(native, configuration)$native_lower)
  bad <- configuration
  bad$lambda_relative_grid[1L] <- 1024
  checks[[length(checks) + 1L]] <- check("v7_unapproved_lambda_grid_rejected",
    inherits(tryCatch(lsg_finite_fitters_v7(bad), error = identity), "error"))
  do.call(rbind, checks)
}
