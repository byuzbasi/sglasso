# Validation-tuned predefined-group logistic benchmarks.
#
# This v2 file deliberately excludes ordinary glmnet because its coordinate-
# level selection does not target the predefined all-in/all-out group support
# used by the SGLASSO simulation estimand.

audit_external_nonfinite_validation_v3 <- function(
    tuning,
    require_selection_state = FALSE
) {
  required <- c(
    "engine", "method_path", "penalty_family", "fit_index",
    "lambda_index", "validation_log_loss", "solver_warning",
    "unexpected_solver_warning", "requested_path_length",
    "returned_path_length", "path_complete",
    "saturated_path_truncation", "iteration_budget_truncation",
    "total_iteration_limit_reached", "returned_lower_boundary_excluded",
    "path_termination_acceptable", "point_converged",
    "finite_validation_loss", "numerically_eligible"
  )
  if (isTRUE(require_selection_state)) {
    required <- c(required, "selected", "invalid_validation_contender")
  }
  missing <- setdiff(required, names(tuning))
  if (length(missing) || !is.data.frame(tuning) || !nrow(tuning)) {
    stop(
      "External tuning audit is missing required nonfinite metadata: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  flag_columns <- c(
    "path_complete", "saturated_path_truncation",
    "iteration_budget_truncation", "total_iteration_limit_reached",
    "returned_lower_boundary_excluded", "path_termination_acceptable",
    "point_converged", "finite_validation_loss", "numerically_eligible"
  )
  if (isTRUE(require_selection_state)) {
    flag_columns <- c(
      flag_columns, "selected", "invalid_validation_contender"
    )
  }
  known_flags <- Reduce(
    `&`,
    lapply(flag_columns, function(column) {
      !is.na(tuning[[column]]) & tuning[[column]] %in% c(TRUE, FALSE)
    })
  )

  actual_finite <- is.finite(tuning$validation_log_loss)
  recorded_finite <- tuning$finite_validation_loss %in% TRUE
  finite_flag_agrees <- !is.na(tuning$finite_validation_loss) &
    recorded_finite == actual_finite
  boundary_excluded <-
    tuning$returned_lower_boundary_excluded %in% TRUE
  iteration_truncated <- tuning$iteration_budget_truncation %in% TRUE
  saturation_truncated <- tuning$saturated_path_truncation %in% TRUE
  path_truncated <- xor(iteration_truncated, saturation_truncated)
  selected_false <- if ("selected" %in% names(tuning)) {
    tuning$selected %in% FALSE
  } else {
    rep(TRUE, nrow(tuning))
  }
  invalid_contender_false <- if (
      "invalid_validation_contender" %in% names(tuning)) {
    tuning$invalid_validation_contender %in% FALSE
  } else {
    rep(TRUE, nrow(tuning))
  }

  fit_key_columns <- c(
    intersect(c("scenario", "replication"), names(tuning)),
    "method_path", "fit_index"
  )
  fit_key <- do.call(
    paste,
    c(unname(tuning[fit_key_columns]), list(sep = "::"))
  )
  point_key <- paste(fit_key, tuning$lambda_index, sep = "::")
  duplicate_point <- duplicated(point_key) | duplicated(point_key, fromLast = TRUE)
  exclusion_count <- ave(
    as.integer(boundary_excluded), fit_key, FUN = sum
  )
  maximum_lambda_index <- ave(
    as.numeric(tuning$lambda_index), fit_key, FUN = max
  )
  truncated_fit <- ave(
    as.integer(path_truncated), fit_key, FUN = max
  ) %in% 1L

  path_level_columns <- c(
    "engine", "penalty_family", "requested_path_length",
    "returned_path_length", "path_complete",
    "saturated_path_truncation", "iteration_budget_truncation",
    "total_iteration_limit_reached", "path_termination_acceptable",
    "solver_warning", "unexpected_solver_warning"
  )
  split_rows <- split(seq_len(nrow(tuning)), fit_key)
  consistent_fit <- vapply(split_rows, function(index) {
    all(vapply(path_level_columns, function(column) {
      length(unique(tuning[[column]][index])) == 1L
    }, logical(1)))
  }, logical(1))
  path_metadata_consistent <- unname(consistent_fit[fit_key])

  no_unexpected_warning <-
    !is.na(tuning$unexpected_solver_warning) &
    !nzchar(tuning$unexpected_solver_warning)
  recognized_saturation_warning <-
    !is.na(tuning$solver_warning) &
    grepl("Model saturated; exiting", tuning$solver_warning, fixed = TRUE)
  recognized_truncation <-
    (iteration_truncated &
       tuning$total_iteration_limit_reached %in% TRUE) |
    (saturation_truncated & recognized_saturation_warning)
  truncation_metadata_valid <-
    tuning$engine %in% "grpreg" &
    tuning$penalty_family %in% c("group_mcp", "group_scad") &
    tuning$path_complete %in% FALSE &
    tuning$path_termination_acceptable %in% TRUE &
    is.finite(tuning$requested_path_length) &
    is.finite(tuning$returned_path_length) &
    tuning$returned_path_length < tuning$requested_path_length &
    path_truncated & recognized_truncation & no_unexpected_warning
  expected_exclusion_structure <-
    (truncated_fit & exclusion_count == 1L) |
    (!truncated_fit & exclusion_count == 0L)
  qualifying_policy_exclusion <-
    boundary_excluded & truncation_metadata_valid &
    tuning$lambda_index == tuning$returned_path_length &
    tuning$lambda_index == maximum_lambda_index &
    tuning$point_converged %in% FALSE &
    tuning$numerically_eligible %in% FALSE &
    selected_false & invalid_contender_false & exclusion_count == 1L &
    !duplicate_point & finite_flag_agrees & known_flags

  candidate_nonfinite <- !actual_finite & !qualifying_policy_exclusion
  policy_excluded_nonfinite <-
    !actual_finite & qualifying_policy_exclusion
  metadata_valid <-
    !any(duplicate_point) &&
    all(finite_flag_agrees) &&
    all(known_flags) &&
    all(!(iteration_truncated & saturation_truncated)) &&
    all(path_metadata_consistent) &&
    all(expected_exclusion_structure) &&
    all(!path_truncated | truncation_metadata_valid) &&
    all(!boundary_excluded | qualifying_policy_exclusion)

  list(
    rows = data.frame(
      actual_finite_validation_loss = actual_finite,
      finite_validation_loss_flag_agrees = finite_flag_agrees,
      qualifying_policy_exclusion = qualifying_policy_exclusion,
      candidate_nonfinite_validation_loss = candidate_nonfinite,
      policy_excluded_nonfinite_validation_loss =
        policy_excluded_nonfinite,
      stringsAsFactors = FALSE
    ),
    nonfinite_validation_loss_points = sum(!actual_finite),
    candidate_nonfinite_validation_loss_points =
      sum(candidate_nonfinite),
    policy_excluded_nonfinite_validation_loss_points =
      sum(policy_excluded_nonfinite),
    metadata_valid = metadata_valid,
    audit_passes = metadata_valid && !any(candidate_nonfinite)
  )
}


select_valid_external_grid_v2 <- function(tuning, method_path) {
  rows <- tuning$method_path == method_path
  candidates <- tuning[rows, , drop = FALSE]
  if (!nrow(candidates)) {
    stop("No tuning candidates were recorded for ", method_path, ".",
         call. = FALSE)
  }
  valid <- candidates[candidates$numerically_eligible, , drop = FALSE]
  if (!nrow(valid)) {
    stop("No numerically eligible validation candidate was available for ",
         method_path, ".", call. = FALSE)
  }

  selected <- valid[which.min(valid$validation_log_loss), , drop = FALSE]
  selected_fit <- candidates[
    candidates$fit_index == selected$fit_index[[1L]], , drop = FALSE
  ]
  returned_lower_index <- max(selected_fit$lambda_index)
  returned_lower_row <- selected_fit[
    selected_fit$lambda_index == returned_lower_index, , drop = FALSE
  ][1L, , drop = FALSE]
  selected_on_returned_lower_boundary <-
    selected$lambda_index[[1L]] == returned_lower_index
  selected_fit_saturated <- any(
    selected_fit$saturated_path_truncation
  )
  selected_fit_iteration_budget_truncation <- if (
      "iteration_budget_truncation" %in% names(selected_fit)) {
    any(selected_fit$iteration_budget_truncation)
  } else {
    FALSE
  }
  returned_lower_minus_selected <-
    returned_lower_row$validation_log_loss[[1L]] -
    selected$validation_log_loss[[1L]]
  truncated_path_selection_interior <-
    !(selected_fit_saturated ||
        selected_fit_iteration_budget_truncation) ||
    !selected_on_returned_lower_boundary
  nonfinite_audit <- audit_external_nonfinite_validation_v3(candidates)
  policy_excluded <-
    nonfinite_audit$rows$qualifying_policy_exclusion
  finite_validation_loss <-
    nonfinite_audit$rows$actual_finite_validation_loss
  invalid_with_loss <- candidates[
    !candidates$numerically_eligible & finite_validation_loss &
      !policy_excluded,
    , drop = FALSE
  ]
  if (nrow(invalid_with_loss)) {
    best_invalid <- min(invalid_with_loss$validation_log_loss)
    invalid_minus_valid <- best_invalid - selected$validation_log_loss
    invalid_competitive <- best_invalid <= selected$validation_log_loss
  } else {
    best_invalid <- NA_real_
    invalid_minus_valid <- NA_real_
    invalid_competitive <- FALSE
  }

  valid_fit_indices <- unique(valid$fit_index)
  requested_fit_indices <- unique(candidates$fit_index)
  list(
    row = selected,
    total_candidate_points = nrow(candidates),
    numerically_eligible_points = nrow(valid),
    numerically_ineligible_points = nrow(candidates) - nrow(valid),
    nonfinite_validation_loss_points =
      nonfinite_audit$nonfinite_validation_loss_points,
    candidate_nonfinite_validation_loss_points =
      nonfinite_audit$candidate_nonfinite_validation_loss_points,
    policy_excluded_nonfinite_validation_loss_points = sum(
      nonfinite_audit$policy_excluded_nonfinite_validation_loss_points
    ),
    nonfinite_validation_audit_passes = nonfinite_audit$audit_passes,
    best_invalid_validation_log_loss = best_invalid,
    invalid_minus_valid_validation_log_loss = invalid_minus_valid,
    invalid_candidate_competitive = invalid_competitive,
    valid_fit_count = length(valid_fit_indices),
    requested_fit_count = length(requested_fit_indices),
    all_fit_paths_have_valid_candidates =
      setequal(valid_fit_indices, requested_fit_indices),
    selected_fit_saturated_path_truncation = selected_fit_saturated,
    selected_fit_iteration_budget_truncation =
      selected_fit_iteration_budget_truncation,
    selected_on_returned_lower_boundary =
      selected_on_returned_lower_boundary,
    returned_lower_boundary_validation_log_loss =
      returned_lower_row$validation_log_loss[[1L]],
    returned_lower_boundary_minus_selected_validation_log_loss =
      returned_lower_minus_selected,
    truncated_path_selection_interior =
      truncated_path_selection_interior
  )
}


mark_external_selection_v2 <- function(tuning, selection) {
  selected <- selection$row
  tuning$selected <- tuning$method_path == selected$method_path &
    tuning$fit_index == selected$fit_index &
    tuning$lambda_index == selected$lambda_index
  tuning$validation_loss_minus_selected <- NA_real_
  method_rows <- tuning$method_path == selected$method_path
  tuning$validation_loss_minus_selected[method_rows] <-
    tuning$validation_log_loss[method_rows] -
    selected$validation_log_loss[[1L]]
  tuning$invalid_validation_contender <- FALSE
  nonfinite_audit <- audit_external_nonfinite_validation_v3(
    tuning, require_selection_state = TRUE
  )
  policy_excluded <-
    nonfinite_audit$rows$qualifying_policy_exclusion
  tuning$candidate_nonfinite_validation_loss <-
    nonfinite_audit$rows$candidate_nonfinite_validation_loss
  tuning$policy_excluded_nonfinite_validation_loss <-
    nonfinite_audit$rows$policy_excluded_nonfinite_validation_loss
  tuning$invalid_validation_contender[method_rows] <-
    !tuning$numerically_eligible[method_rows] &
    tuning$finite_validation_loss[method_rows] &
    !policy_excluded[method_rows] &
    tuning$validation_log_loss[method_rows] <=
      selected$validation_log_loss[[1L]]
  tuning
}


classify_grpreg_path_v2 <- function(
    penalty,
    iterations,
    requested_path_length,
    max_iterations,
    warnings = character(0)
) {
  iterations <- as.numeric(iterations)
  requested_path_length <- as.integer(requested_path_length)
  max_iterations <- as.numeric(max_iterations)
  warnings <- unique(as.character(warnings))
  warnings <- warnings[!is.na(warnings) & nzchar(warnings)]
  if (!penalty %in% c("grLasso", "grMCP", "grSCAD") ||
      requested_path_length < 2L || !length(iterations) ||
      length(iterations) > requested_path_length ||
      any(!is.finite(iterations)) || any(iterations < 0) ||
      !is.finite(max_iterations) || max_iterations < 1) {
    stop("Invalid grpreg path-classification inputs.", call. = FALSE)
  }

  returned_path_length <- length(iterations)
  path_complete <- returned_path_length == requested_path_length
  saturation_warning <- grepl(
    "Model saturated; exiting", warnings, fixed = TRUE
  )
  iteration_warning <- grepl(
    "Algorithm failed to converge for all values of lambda",
    warnings,
    fixed = TRUE
  )
  total_iterations <- sum(iterations)
  if (total_iterations > max_iterations) {
    stop("grpreg iterations exceed the declared total budget.",
         call. = FALSE)
  }
  total_iteration_limit_reached <-
    total_iterations >= max_iterations
  nonconvex_penalty <- penalty %in% c("grMCP", "grSCAD")
  saturated_path_truncation <-
    nonconvex_penalty && !path_complete && any(saturation_warning)
  iteration_budget_truncation <-
    nonconvex_penalty && !path_complete && !saturated_path_truncation &&
    total_iteration_limit_reached
  recognized_warning <-
    (saturation_warning & saturated_path_truncation) |
    (iteration_warning & iteration_budget_truncation)
  unexpected_warnings <- warnings[!recognized_warning]
  path_termination_acceptable <-
    !length(unexpected_warnings) &&
    ((path_complete && !total_iteration_limit_reached) ||
       saturated_path_truncation ||
       iteration_budget_truncation)
  returned_lower_boundary_excluded <-
    rep(FALSE, returned_path_length)
  if (saturated_path_truncation || iteration_budget_truncation) {
    returned_lower_boundary_excluded[returned_path_length] <- TRUE
  }
  point_converged <- iterations < max_iterations
  if (saturated_path_truncation || iteration_budget_truncation) {
    point_converged[returned_path_length] <- FALSE
  }

  list(
    path_complete = path_complete,
    saturated_path_truncation = saturated_path_truncation,
    iteration_budget_truncation = iteration_budget_truncation,
    total_iterations = total_iterations,
    total_iteration_limit_reached = total_iteration_limit_reached,
    returned_lower_boundary_excluded =
      returned_lower_boundary_excluded,
    point_converged = point_converged,
    path_termination_acceptable = path_termination_acceptable,
    unique_warnings = warnings,
    unexpected_warnings = unexpected_warnings
  )
}


fit_adelie_validation_grid_v2 <- function(
    X,
    y,
    group,
    X_validation,
    y_validation,
    alpha_grid,
    nlambda,
    lambda_min_ratio,
    ridge_boundary_upper_multiplier = 8,
    tolerance = 1e-7,
    max_iterations = 100000L,
    irls_tolerance = 1e-7,
    irls_max_iterations = 10000L
) {
  if (!requireNamespace("adelie", quietly = TRUE)) {
    stop("The installed adelie package is required.", call. = FALSE)
  }
  input <- validate_penalized_benchmark_grid_v1(
    X, y, group, X_validation, y_validation, alpha_grid,
    nlambda, lambda_min_ratio
  )
  max_iterations <- as.integer(max_iterations)
  irls_max_iterations <- as.integer(irls_max_iterations)
  if (max_iterations < 1L || irls_max_iterations < 1L ||
      !is.finite(ridge_boundary_upper_multiplier) ||
      ridge_boundary_upper_multiplier <= 1 ||
      !is.finite(tolerance) || tolerance <= 0 ||
      !is.finite(irls_tolerance) || irls_tolerance <= 0) {
    stop("Invalid adelie numerical controls.", call. = FALSE)
  }

  preprocess <- prepare_lsg_design(input$X, input$group)
  X_validation_solver <- transform_lsg_newx(
    preprocess, input$X_validation
  )
  group_starts <- as.integer(preprocess$group_start + 1L)
  null_score_scale <- lsg_null_score_lambda_scale(preprocess, input$y)
  coefficient_method <- getS3method(
    "coef", "grpnet", envir = asNamespace("adelie")
  )
  prediction_method <- getS3method(
    "predict", "grpnet", envir = asNamespace("adelie")
  )

  fits <- vector("list", length(input$alpha_grid))
  coefficient_paths <- vector("list", length(input$alpha_grid))
  elapsed <- numeric(length(input$alpha_grid))
  warning_messages <- vector("list", length(input$alpha_grid))
  tuning <- vector("list", length(input$alpha_grid))

  for (ai in seq_along(input$alpha_grid)) {
    alpha <- input$alpha_grid[ai]
    alpha_zero <- abs(alpha) <= 1e-12
    explicit_lambda <- if (alpha_zero) {
      exp(seq(
        log(null_score_scale * ridge_boundary_upper_multiplier),
        log(null_score_scale * input$lambda_min_ratio),
        length.out = input$nlambda
      ))
    } else {
      NULL
    }
    captured_warnings <- character(0)
    elapsed[ai] <- system.time({
      fit <- withCallingHandlers(
        do.call(adelie::grpnet, list(
          X = preprocess$X,
          glm = adelie::glm.binomial(input$y),
          groups = group_starts,
          alpha = alpha,
          penalty = preprocess$group_weight,
          standardize = FALSE,
          intercept = TRUE,
          lmda_path_size = input$nlambda,
          min_ratio = input$lambda_min_ratio,
          tol = tolerance,
          max_iters = max_iterations,
          irls_tol = irls_tolerance,
          irls_max_iters = irls_max_iterations,
          screen_rule = "strong",
          early_exit = FALSE,
          check_state = TRUE,
          progress_bar = FALSE,
          n_threads = 1L,
          lambda = explicit_lambda
        )),
        warning = function(condition) {
          captured_warnings <<- c(
            captured_warnings, conditionMessage(condition)
          )
          invokeRestart("muffleWarning")
        }
      )
      coefficient_path <- coefficient_method(fit)
      lambda <- as.numeric(coefficient_path$lambda)
      probability <- as.matrix(prediction_method(
        fit,
        newx = X_validation_solver,
        lambda = lambda,
        type = "response",
        n_threads = 1L
      ))
    })[["elapsed"]]

    beta <- as.matrix(coefficient_path$betas)
    intercept <- as.numeric(coefficient_path$intercepts[, 1L])
    path_length <- length(lambda)
    if (nrow(beta) != path_length || ncol(beta) != ncol(preprocess$X) ||
        length(intercept) != path_length ||
        nrow(probability) != length(input$y_validation) ||
        ncol(probability) != path_length) {
      stop("adelie returned incompatible coefficient or prediction arrays.",
           call. = FALSE)
    }
    state_lambda <- as.numeric(fit$state$lmdas)
    lambda_aligned <- length(state_lambda) == path_length &&
      (path_length == 0L || max(abs(state_lambda - lambda)) <=
         1e-10 * max(1, max(abs(lambda))))
    path_complete <- path_length == input$nlambda && lambda_aligned
    finite_coefficient <- if (path_length) {
      apply(cbind(intercept, beta), 1L, function(value) {
        all(is.finite(value))
      })
    } else {
      logical(0)
    }
    finite_probability <- if (path_length) {
      apply(probability, 2L, function(value) all(is.finite(value)))
    } else {
      logical(0)
    }
    validation_loss <- if (path_length) {
      vapply(seq_len(path_length), function(index) {
        if (finite_probability[index]) {
          binary_log_loss(input$y_validation, probability[, index])
        } else {
          NA_real_
        }
      }, numeric(1))
    } else {
      numeric(0)
    }
    finite_validation_loss <- is.finite(validation_loss)
    warning_free <- length(captured_warnings) == 0L
    path_termination_acceptable <- path_complete & warning_free
    numerically_eligible <- path_termination_acceptable &
      finite_coefficient & finite_probability & finite_validation_loss
    lambda_reference <- if (alpha_zero) null_score_scale else lambda[1L]
    lambda_reference_type <- if (alpha_zero) {
      "null_score_ridge_boundary"
    } else {
      "native_adelie_path_start"
    }

    fits[[ai]] <- fit
    coefficient_paths[[ai]] <- list(
      beta = beta,
      intercept = intercept,
      lambda = lambda
    )
    warning_messages[[ai]] <- unique(captured_warnings)
    tuning[[ai]] <- data.frame(
      engine = "adelie",
      method_path = "Logistic Group Elastic Net (adelie)",
      penalty_family = "group_elastic_net",
      fit_index = ai,
      alpha_index = ai,
      alpha = alpha,
      gamma = NA_real_,
      lambda_index = seq_len(path_length),
      lambda = lambda,
      lambda_fraction = lambda / lambda[1L],
      lambda_reference = lambda_reference,
      lambda_reference_type = lambda_reference_type,
      lambda_relative_to_reference = lambda / lambda_reference,
      alpha_zero_ridge_boundary = alpha_zero,
      group_selection_capable = !alpha_zero,
      validation_log_loss = validation_loss,
      passes = NA_real_,
      solver_warning = paste(unique(captured_warnings), collapse = " | "),
      unexpected_solver_warning = paste(
        unique(captured_warnings), collapse = " | "
      ),
      requested_path_length = input$nlambda,
      returned_path_length = path_length,
      path_complete = path_complete,
      saturated_path_truncation = FALSE,
      iteration_budget_truncation = FALSE,
      total_iterations = NA_real_,
      total_iteration_limit_reached = FALSE,
      returned_lower_boundary_excluded = FALSE,
      path_termination_acceptable = path_termination_acceptable,
      point_converged = path_termination_acceptable,
      finite_validation_loss = finite_validation_loss,
      numerically_eligible = numerically_eligible,
      stringsAsFactors = FALSE
    )
  }

  tuning <- do.call(rbind, tuning)
  rownames(tuning) <- NULL
  method_path <- "Logistic Group Elastic Net (adelie)"
  selection <- select_valid_external_grid_v2(tuning, method_path)
  tuning <- mark_external_selection_v2(tuning, selection)

  list(
    fits = fits,
    coefficient_paths = coefficient_paths,
    preprocess = preprocess,
    alpha_grid = input$alpha_grid,
    tuning = tuning,
    selection = selection,
    elapsed_by_alpha = elapsed,
    runtime_seconds = sum(elapsed),
    warning_messages = warning_messages,
    tolerance = tolerance,
    max_iterations = max_iterations,
    irls_tolerance = irls_tolerance,
    irls_max_iterations = irls_max_iterations,
    ridge_boundary_upper_multiplier = ridge_boundary_upper_multiplier,
    null_score_scale = null_score_scale
  )
}


extract_adelie_group_en_solution_v2 <- function(grid, newx) {
  selection <- grid$selection$row[1L, , drop = FALSE]
  fit_index <- selection$fit_index[[1L]]
  lambda_index <- selection$lambda_index[[1L]]
  path <- grid$coefficient_paths[[fit_index]]
  coefficient <- recover_lsg_coefficients(
    grid$preprocess,
    path$beta[lambda_index, ],
    path$intercept[lambda_index]
  )
  prediction_method <- getS3method(
    "predict", "grpnet", envir = asNamespace("adelie")
  )
  probability <- drop(prediction_method(
    grid$fits[[fit_index]],
    newx = transform_lsg_newx(grid$preprocess, as.matrix(newx)),
    lambda = path$lambda[lambda_index],
    type = "response",
    n_threads = 1L
  ))
  list(
    coefficient = as.numeric(coefficient),
    probability = probability,
    alpha = selection$alpha[[1L]],
    gamma = NA_real_,
    penalty_family = selection$penalty_family[[1L]],
    lambda = selection$lambda[[1L]],
    lambda_fraction = selection$lambda_fraction[[1L]],
    lambda_index = lambda_index,
    path_length = nrow(path$beta),
    validation_log_loss = selection$validation_log_loss[[1L]],
    converged = isTRUE(selection$numerically_eligible[[1L]])
  )
}


default_grpreg_penalty_specification_v2 <- function() {
  data.frame(
    penalty = c("grLasso", "grMCP", "grSCAD"),
    method_path = c(
      "Logistic Group Lasso (grpreg)",
      "Logistic Group MCP (grpreg)",
      "Logistic Group SCAD (grpreg)"
    ),
    penalty_family = c("group_lasso", "group_mcp", "group_scad"),
    gamma = c(NA_real_, 3, 4),
    stringsAsFactors = FALSE
  )
}


fit_grpreg_validation_paths_v2 <- function(
    X,
    y,
    group,
    X_validation,
    y_validation,
    nlambda,
    lambda_min_ratio,
    tolerance = 1e-7,
    max_iterations = 1000000L,
    penalty_specification = default_grpreg_penalty_specification_v2()
) {
  if (!requireNamespace("grpreg", quietly = TRUE)) {
    stop("The installed grpreg package is required.", call. = FALSE)
  }
  input <- validate_penalized_benchmark_grid_v1(
    X, y, group, X_validation, y_validation,
    alpha_grid = 1, nlambda = nlambda,
    lambda_min_ratio = lambda_min_ratio
  )
  max_iterations <- as.integer(max_iterations)
  required_specification <- c(
    "penalty", "method_path", "penalty_family", "gamma"
  )
  if (!all(required_specification %in% names(penalty_specification)) ||
      anyDuplicated(penalty_specification$penalty) ||
      anyDuplicated(penalty_specification$method_path) ||
      !setequal(
        penalty_specification$penalty,
        c("grLasso", "grMCP", "grSCAD")
      ) || max_iterations < 1L || !is.finite(tolerance) ||
      tolerance <= 0) {
    stop("Invalid grpreg penalty specification or numerical controls.",
         call. = FALSE)
  }

  group_factor <- factor(input$group, levels = unique(input$group))
  group_multiplier <- sqrt(as.numeric(table(group_factor)))
  fits <- vector("list", nrow(penalty_specification))
  elapsed <- numeric(nrow(penalty_specification))
  warning_messages <- vector("list", nrow(penalty_specification))
  tuning <- vector("list", nrow(penalty_specification))

  for (pi in seq_len(nrow(penalty_specification))) {
    specification <- penalty_specification[pi, , drop = FALSE]
    gamma_for_fit <- if (specification$penalty == "grSCAD") 4 else 3
    captured_warnings <- character(0)
    elapsed[pi] <- system.time({
      fit <- withCallingHandlers(
        grpreg::grpreg(
          X = input$X,
          y = input$y,
          group = input$group,
          penalty = specification$penalty,
          family = "binomial",
          nlambda = input$nlambda,
          lambda.min = input$lambda_min_ratio,
          log.lambda = TRUE,
          alpha = 1,
          eps = tolerance,
          max.iter = max_iterations,
          dfmax = ncol(input$X),
          gmax = length(unique(input$group)),
          gamma = gamma_for_fit,
          group.multiplier = group_multiplier,
          warn = TRUE,
          returnX = FALSE
        ),
        warning = function(condition) {
          captured_warnings <<- c(
            captured_warnings, conditionMessage(condition)
          )
          invokeRestart("muffleWarning")
        }
      )
      probability <- as.matrix(stats::predict(
        fit,
        X = input$X_validation,
        type = "response"
      ))
    })[["elapsed"]]

    path_length <- length(fit$lambda)
    if (ncol(fit$beta) != path_length || length(fit$iter) != path_length ||
        nrow(probability) != length(input$y_validation) ||
        ncol(probability) != path_length) {
      stop("grpreg returned incompatible coefficient or prediction arrays.",
           call. = FALSE)
    }
    path_status <- classify_grpreg_path_v2(
      penalty = specification$penalty[[1L]],
      iterations = fit$iter,
      requested_path_length = input$nlambda,
      max_iterations = max_iterations,
      warnings = captured_warnings
    )
    finite_coefficient <- apply(fit$beta, 2L, function(value) {
      all(is.finite(value))
    })
    finite_probability <- apply(probability, 2L, function(value) {
      all(is.finite(value))
    })
    validation_loss <- vapply(seq_len(path_length), function(index) {
      if (finite_probability[index]) {
        binary_log_loss(input$y_validation, probability[, index])
      } else {
        NA_real_
      }
    }, numeric(1))
    finite_validation_loss <- is.finite(validation_loss)
    numerically_eligible <- path_status$path_termination_acceptable &
      path_status$point_converged &
      !path_status$returned_lower_boundary_excluded &
      finite_coefficient & finite_probability & finite_validation_loss

    fits[[pi]] <- fit
    warning_messages[[pi]] <- unique(captured_warnings)
    tuning[[pi]] <- data.frame(
      engine = "grpreg",
      method_path = specification$method_path,
      penalty_family = specification$penalty_family,
      fit_index = pi,
      alpha_index = 1L,
      alpha = 1,
      gamma = specification$gamma,
      lambda_index = seq_len(path_length),
      lambda = fit$lambda,
      lambda_fraction = fit$lambda / fit$lambda[1L],
      lambda_reference = fit$lambda[1L],
      lambda_reference_type = "native_grpreg_path_start",
      lambda_relative_to_reference = fit$lambda / fit$lambda[1L],
      alpha_zero_ridge_boundary = FALSE,
      group_selection_capable = TRUE,
      validation_log_loss = validation_loss,
      passes = as.numeric(fit$iter),
      solver_warning = paste(
        path_status$unique_warnings, collapse = " | "
      ),
      unexpected_solver_warning = paste(
        path_status$unexpected_warnings, collapse = " | "
      ),
      requested_path_length = input$nlambda,
      returned_path_length = path_length,
      path_complete = path_status$path_complete,
      saturated_path_truncation =
        path_status$saturated_path_truncation,
      iteration_budget_truncation =
        path_status$iteration_budget_truncation,
      total_iterations = path_status$total_iterations,
      total_iteration_limit_reached =
        path_status$total_iteration_limit_reached,
      returned_lower_boundary_excluded =
        path_status$returned_lower_boundary_excluded,
      path_termination_acceptable =
        path_status$path_termination_acceptable,
      point_converged = path_status$point_converged,
      finite_validation_loss = finite_validation_loss,
      numerically_eligible = numerically_eligible,
      stringsAsFactors = FALSE
    )
  }

  names(fits) <- penalty_specification$penalty
  names(elapsed) <- penalty_specification$penalty
  names(warning_messages) <- penalty_specification$penalty
  tuning <- do.call(rbind, tuning)
  rownames(tuning) <- NULL
  selections <- vector("list", nrow(penalty_specification))
  names(selections) <- penalty_specification$penalty
  tuning$selected <- FALSE
  tuning$validation_loss_minus_selected <- NA_real_
  tuning$invalid_validation_contender <- FALSE
  for (pi in seq_len(nrow(penalty_specification))) {
    method_path <- penalty_specification$method_path[pi]
    selections[[pi]] <- select_valid_external_grid_v2(
      tuning, method_path
    )
    marked <- mark_external_selection_v2(tuning, selections[[pi]])
    method_rows <- tuning$method_path == method_path
    tuning[method_rows, c(
      "selected", "validation_loss_minus_selected",
      "invalid_validation_contender",
      "candidate_nonfinite_validation_loss",
      "policy_excluded_nonfinite_validation_loss"
    )] <- marked[method_rows, c(
      "selected", "validation_loss_minus_selected",
      "invalid_validation_contender",
      "candidate_nonfinite_validation_loss",
      "policy_excluded_nonfinite_validation_loss"
    )]
  }

  list(
    fits = fits,
    penalty_specification = penalty_specification,
    tuning = tuning,
    selections = selections,
    elapsed_by_penalty = elapsed,
    runtime_seconds = sum(elapsed),
    warning_messages = warning_messages,
    tolerance = tolerance,
    max_iterations = max_iterations
  )
}


extract_grpreg_group_solution_v2 <- function(grid, penalty, newx) {
  if (!penalty %in% names(grid$fits)) {
    stop("Unknown grpreg penalty: ", penalty, call. = FALSE)
  }
  selection <- grid$selections[[penalty]]$row[1L, , drop = FALSE]
  fit <- grid$fits[[penalty]]
  lambda_index <- selection$lambda_index[[1L]]
  coefficient <- as.numeric(fit$beta[, lambda_index])
  probability <- drop(stats::predict(
    fit,
    X = as.matrix(newx),
    type = "response",
    which = lambda_index
  ))
  list(
    coefficient = coefficient,
    probability = probability,
    alpha = selection$alpha[[1L]],
    gamma = selection$gamma[[1L]],
    penalty_family = selection$penalty_family[[1L]],
    lambda = selection$lambda[[1L]],
    lambda_fraction = selection$lambda_fraction[[1L]],
    lambda_index = lambda_index,
    path_length = length(fit$lambda),
    validation_log_loss = selection$validation_log_loss[[1L]],
    converged = isTRUE(selection$numerically_eligible[[1L]])
  )
}


evaluate_group_benchmark_solution_v2 <- function(
    solution,
    data,
    scenario,
    replication,
    method,
    engine,
    runtime_seconds,
    runtime_scope
) {
  result <- evaluate_penalized_benchmark_v1(
    solution = solution,
    data = data,
    scenario = scenario,
    replication = replication,
    method = method,
    engine = engine,
    runtime_seconds = runtime_seconds,
    runtime_scope = runtime_scope
  )
  cbind(
    result,
    data.frame(
      penalty_family = solution$penalty_family,
      selected_gamma = solution$gamma,
      selection_unit = "predefined_group",
      stringsAsFactors = FALSE
    )
  )
}
