# Method-comparable lambda paths for the ALL BCR/ABL versus NEG V2 study.
#
# This is an analysis overlay.  It does not modify the installed sglasso
# package, the Rcpp/RcppArmadillo solver, or the frozen V19/V1 sources.

allb_common_lambda_grid_v2 <- function() {
  c(exp(seq(log(1), log(0.05), length.out = 30L)),
    0.025, 0.0125, 0.00625, 0.003125)
}

allb_sglasso_lambda_grid_v2 <- function() {
  c(8192, 4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2,
    allb_common_lambda_grid_v2())
}

allb_grid_equal_v2 <- function(x, y, tolerance = 1e-12) {
  length(x) == length(y) && isTRUE(all.equal(
    as.numeric(x), as.numeric(y), tolerance = tolerance,
    check.attributes = FALSE
  ))
}

allb_validate_lambda_grids_v2 <- function() {
  common <- allb_common_lambda_grid_v2()
  sglasso <- allb_sglasso_lambda_grid_v2()
  stopifnot(
    length(common) == 34L,
    length(sglasso) == 47L,
    all(is.finite(common)), all(common > 0),
    all(is.finite(sglasso)), all(sglasso > 0),
    all(diff(common) < 0), all(diff(sglasso) < 0),
    abs(common[[1L]] - 1) <= 1e-15,
    abs(tail(common, 1L) - 0.003125) <= 1e-15,
    all(vapply(common, function(value) {
      any(abs(sglasso - value) <= 1e-13)
    }, logical(1)))
  )
  invisible(TRUE)
}

# grpreg can stop a binomial path before the requested lower ratios when the
# fitted model becomes saturated.  V2 treats this as a solver-defined feasible
# prefix, irrespective of whether the penalty is convex or nonconvex.  The
# final returned point is conservatively excluded because it is the point at
# which the terminal condition was reported.  No missing ratio is fabricated.
allb_classify_grpreg_path_v2 <- function(
    penalty, iterations, requested_path_length, max_iterations,
    warnings = character(0)
) {
  iterations <- as.numeric(iterations)
  requested_path_length <- as.integer(requested_path_length)
  max_iterations <- as.numeric(max_iterations)
  warnings <- unique(as.character(warnings))
  warnings <- warnings[!is.na(warnings) & nzchar(warnings)]
  if (!penalty %in% c("grLasso", "grMCP", "grSCAD") ||
      length(requested_path_length) != 1L || requested_path_length < 2L ||
      !length(iterations) || length(iterations) > requested_path_length ||
      any(!is.finite(iterations)) || any(iterations < 0) ||
      length(max_iterations) != 1L || !is.finite(max_iterations) ||
      max_iterations < 1) {
    stop("Invalid ALL V2 grpreg path-classification inputs.",
      call. = FALSE)
  }

  returned <- length(iterations)
  complete <- returned == requested_path_length
  total <- sum(iterations)
  if (total > max_iterations) {
    stop("grpreg iterations exceed the declared total budget.",
      call. = FALSE)
  }
  budget_reached <- total >= max_iterations
  saturation_warning <- grepl(
    "Model saturated; exiting", warnings, fixed = TRUE
  )
  iteration_warning <- grepl(
    "Algorithm failed to converge for all values of lambda",
    warnings, fixed = TRUE
  )
  saturated <- !complete && any(saturation_warning)
  budget_truncated <- !complete && !saturated && budget_reached
  recognized <- (saturation_warning & saturated) |
    (iteration_warning & budget_truncated)
  unexpected <- warnings[!recognized]
  acceptable <- !length(unexpected) &&
    ((complete && !budget_reached) || saturated || budget_truncated)
  excluded <- rep(FALSE, returned)
  if (saturated || budget_truncated) excluded[[returned]] <- TRUE
  point_converged <- iterations < max_iterations
  if (saturated || budget_truncated) point_converged[[returned]] <- FALSE

  list(
    path_complete = complete,
    saturated_path_truncation = saturated,
    iteration_budget_truncation = budget_truncated,
    total_iterations = total,
    total_iteration_limit_reached = budget_reached,
    returned_lower_boundary_excluded = excluded,
    point_converged = point_converged,
    path_termination_acceptable = acceptable,
    safe_prefix_length = as.integer(
      if (saturated || budget_truncated) returned - 1L else returned
    ),
    unique_warnings = warnings,
    unexpected_warnings = unexpected
  )
}

allb_select_grpreg_safe_prefix_v2 <- function(tuning, method_path) {
  candidates <- tuning[tuning$method_path == method_path, , drop = FALSE]
  if (!nrow(candidates)) {
    stop("No tuning candidates were recorded for ", method_path, ".",
      call. = FALSE)
  }
  required <- c(
    "lambda_index", "lambda_relative_to_reference", "validation_log_loss",
    "finite_validation_loss", "finite_coefficient", "finite_probability",
    "path_complete", "saturated_path_truncation",
    "iteration_budget_truncation", "returned_lower_boundary_excluded",
    "path_termination_acceptable", "point_converged",
    "numerically_eligible", "unexpected_solver_warning",
    "requested_path_length", "returned_path_length",
    "requested_grid_prefix_aligned", "safe_prefix_length"
  )
  if (!all(required %in% names(candidates))) {
    stop("ALL V2 grpreg safe-prefix metadata are incomplete.",
      call. = FALSE)
  }
  candidates <- candidates[order(candidates$lambda_index), , drop = FALSE]
  n <- nrow(candidates)
  truncated <- candidates$saturated_path_truncation[[1L]] %in% TRUE ||
    candidates$iteration_budget_truncation[[1L]] %in% TRUE
  expected_exclusion <- rep(FALSE, n)
  if (truncated) expected_exclusion[[n]] <- TRUE
  eligible <- candidates$path_termination_acceptable %in% TRUE &
    candidates$point_converged %in% TRUE &
    candidates$returned_lower_boundary_excluded %in% FALSE &
    candidates$finite_validation_loss %in% TRUE &
    candidates$finite_coefficient %in% TRUE &
    candidates$finite_probability %in% TRUE
  structural <-
    identical(as.integer(candidates$lambda_index), seq_len(n)) &&
    all(candidates$requested_grid_prefix_aligned %in% TRUE) &&
    all(candidates$requested_path_length ==
      candidates$requested_path_length[[1L]]) &&
    all(candidates$returned_path_length == n) &&
    all(candidates$safe_prefix_length == sum(!expected_exclusion)) &&
    identical(
      as.logical(candidates$returned_lower_boundary_excluded),
      expected_exclusion
    ) &&
    all(!is.na(candidates$unexpected_solver_warning) &
      !nzchar(candidates$unexpected_solver_warning)) &&
    identical(as.logical(candidates$numerically_eligible), eligible) &&
    all(is.finite(candidates$validation_log_loss))
  if (!isTRUE(structural)) {
    stop("Invalid ALL V2 grpreg safe-prefix evidence for ", method_path,
      ".", call. = FALSE)
  }

  valid <- candidates[eligible, , drop = FALSE]
  if (!nrow(valid)) {
    stop("No numerically eligible validation candidate was available for ",
      method_path, ".", call. = FALSE)
  }
  selected <- valid[which.min(valid$validation_log_loss), , drop = FALSE]
  excluded_with_loss <- candidates[
    !eligible & is.finite(candidates$validation_log_loss), , drop = FALSE
  ]
  competitive <- nrow(excluded_with_loss) > 0L &&
    min(excluded_with_loss$validation_log_loss) <=
      selected$validation_log_loss[[1L]]

  list(
    row = selected,
    total_candidate_points = n,
    numerically_eligible_points = nrow(valid),
    numerically_ineligible_points = n - nrow(valid),
    nonfinite_validation_loss_points = 0L,
    candidate_nonfinite_validation_loss_points = 0L,
    policy_excluded_nonfinite_validation_loss_points = 0L,
    nonfinite_validation_audit_passes = TRUE,
    best_invalid_validation_log_loss = if (nrow(excluded_with_loss)) {
      min(excluded_with_loss$validation_log_loss)
    } else NA_real_,
    invalid_minus_valid_validation_log_loss = if (nrow(excluded_with_loss)) {
      min(excluded_with_loss$validation_log_loss) -
        selected$validation_log_loss[[1L]]
    } else NA_real_,
    invalid_candidate_competitive = competitive,
    valid_fit_count = 1L,
    requested_fit_count = 1L,
    all_fit_paths_have_valid_candidates = TRUE,
    selected_fit_saturated_path_truncation =
      candidates$saturated_path_truncation[[1L]],
    selected_fit_iteration_budget_truncation =
      candidates$iteration_budget_truncation[[1L]],
    selected_on_returned_lower_boundary = FALSE,
    returned_lower_boundary_validation_log_loss =
      candidates$validation_log_loss[[n]],
    returned_lower_boundary_minus_selected_validation_log_loss =
      candidates$validation_log_loss[[n]] -
        selected$validation_log_loss[[1L]],
    truncated_path_selection_interior = TRUE
  )
}

allb_fit_adelie_common_grid_v2 <- function(
    X, y, group, X_validation, y_validation, alpha_grid, nlambda,
    lambda_min_ratio, ridge_boundary_upper_multiplier = 8192,
    tolerance = 1e-7, max_iterations = 100000L,
    irls_tolerance = 1e-7, irls_max_iterations = 10000L,
    common_relative_grid = allb_common_lambda_grid_v2(),
    ridge_relative_grid = allb_sglasso_lambda_grid_v2()
) {
  if (!requireNamespace("adelie", quietly = TRUE)) {
    stop("The installed adelie package is required.", call. = FALSE)
  }
  allb_validate_lambda_grids_v2()
  input <- validate_penalized_benchmark_grid_v1(
    X, y, group, X_validation, y_validation, alpha_grid,
    nlambda, lambda_min_ratio
  )
  max_iterations <- as.integer(max_iterations)
  irls_max_iterations <- as.integer(irls_max_iterations)
  if (!allb_grid_equal_v2(common_relative_grid,
        allb_common_lambda_grid_v2()) ||
      !allb_grid_equal_v2(ridge_relative_grid,
        allb_sglasso_lambda_grid_v2()) ||
      input$nlambda != length(common_relative_grid) ||
      abs(input$lambda_min_ratio - tail(common_relative_grid, 1L)) > 1e-15 ||
      abs(ridge_boundary_upper_multiplier - ridge_relative_grid[[1L]]) > 1e-12 ||
      max_iterations < 1L || irls_max_iterations < 1L ||
      !is.finite(tolerance) || tolerance <= 0 ||
      !is.finite(irls_tolerance) || irls_tolerance <= 0) {
    stop("Invalid Adelie V2 common-grid controls.", call. = FALSE)
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
  base_arguments <- list(
    X = preprocess$X,
    glm = adelie::glm.binomial(input$y),
    groups = group_starts,
    penalty = preprocess$group_weight,
    standardize = FALSE,
    intercept = TRUE,
    tol = tolerance,
    max_iters = max_iterations,
    irls_tol = irls_tolerance,
    irls_max_iters = irls_max_iterations,
    screen_rule = "strong",
    early_exit = FALSE,
    check_state = TRUE,
    progress_bar = FALSE,
    n_threads = 1L
  )

  fits <- vector("list", length(input$alpha_grid))
  coefficient_paths <- vector("list", length(input$alpha_grid))
  elapsed <- numeric(length(input$alpha_grid))
  warning_messages <- vector("list", length(input$alpha_grid))
  tuning <- vector("list", length(input$alpha_grid))
  lambda_reference_values <- numeric(length(input$alpha_grid))

  for (ai in seq_along(input$alpha_grid)) {
    alpha <- input$alpha_grid[[ai]]
    alpha_zero <- abs(alpha) <= 1e-12
    captured_warnings <- character(0)
    elapsed[[ai]] <- system.time({
      if (alpha_zero) {
        lambda_reference <- null_score_scale
        relative_grid <- ridge_relative_grid
      } else {
        pilot <- withCallingHandlers(
          do.call(adelie::grpnet, c(base_arguments, list(
            alpha = alpha, lmda_path_size = 3L, min_ratio = 0.8
          ))),
          warning = function(condition) {
            captured_warnings <<- c(
              captured_warnings, conditionMessage(condition)
            )
            invokeRestart("muffleWarning")
          }
        )
        pilot_path <- coefficient_method(pilot)
        pilot_lambda <- as.numeric(pilot_path$lambda)
        if (!length(pilot_lambda) || !is.finite(pilot_lambda[[1L]]) ||
            pilot_lambda[[1L]] <= 0) {
          stop("Adelie did not return a finite native lambda reference.",
            call. = FALSE)
        }
        lambda_reference <- pilot_lambda[[1L]]
        relative_grid <- common_relative_grid
      }
      explicit_lambda <- lambda_reference * relative_grid
      fit <- withCallingHandlers(
        do.call(adelie::grpnet, c(base_arguments, list(
          alpha = alpha,
          lmda_path_size = length(explicit_lambda),
          min_ratio = tail(relative_grid, 1L),
          lambda = explicit_lambda
        ))),
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
        fit, newx = X_validation_solver, lambda = lambda,
        type = "response", n_threads = 1L
      ))
    })[["elapsed"]]

    beta <- as.matrix(coefficient_path$betas)
    intercept <- as.numeric(coefficient_path$intercepts[, 1L])
    path_length <- length(lambda)
    requested_path_length <- length(relative_grid)
    lambda_reference_values[[ai]] <- lambda_reference
    if (nrow(beta) != path_length || ncol(beta) != ncol(preprocess$X) ||
        length(intercept) != path_length ||
        nrow(probability) != length(input$y_validation) ||
        ncol(probability) != path_length) {
      stop("Adelie V2 returned incompatible coefficient arrays.",
        call. = FALSE)
    }
    state_lambda <- as.numeric(fit$state$lmdas)
    lambda_aligned <- path_length == requested_path_length &&
      allb_grid_equal_v2(lambda / lambda_reference, relative_grid) &&
      length(state_lambda) == path_length &&
      (path_length == 0L || max(abs(state_lambda - lambda)) <=
        1e-10 * max(1, max(abs(lambda))))
    finite_coefficient <- if (path_length) {
      apply(cbind(intercept, beta), 1L, function(value) {
        all(is.finite(value))
      })
    } else logical(0)
    finite_probability <- if (path_length) {
      apply(probability, 2L, function(value) all(is.finite(value)))
    } else logical(0)
    validation_loss <- if (path_length) {
      vapply(seq_len(path_length), function(index) {
        if (finite_probability[[index]]) {
          binary_log_loss(input$y_validation, probability[, index])
        } else NA_real_
      }, numeric(1))
    } else numeric(0)
    finite_validation_loss <- is.finite(validation_loss)
    warning_free <- length(captured_warnings) == 0L
    path_termination_acceptable <- lambda_aligned && warning_free
    numerically_eligible <- path_termination_acceptable &
      finite_coefficient & finite_probability & finite_validation_loss
    lambda_reference_type <- if (alpha_zero) {
      "null_score_ridge_boundary"
    } else {
      "native_adelie_path_start"
    }

    fits[[ai]] <- fit
    coefficient_paths[[ai]] <- list(
      beta = beta, intercept = intercept, lambda = lambda
    )
    warning_messages[[ai]] <- unique(captured_warnings)
    tuning[[ai]] <- data.frame(
      engine = "adelie",
      method_path = "Logistic Group Elastic Net (adelie)",
      penalty_family = "group_elastic_net",
      fit_index = ai, alpha_index = ai, alpha = alpha,
      gamma = NA_real_, lambda_index = seq_len(path_length),
      lambda = lambda, lambda_fraction = lambda / lambda[[1L]],
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
      requested_path_length = requested_path_length,
      returned_path_length = path_length,
      path_complete = lambda_aligned,
      saturated_path_truncation = FALSE,
      iteration_budget_truncation = FALSE,
      total_iterations = NA_real_,
      total_iteration_limit_reached = FALSE,
      returned_lower_boundary_excluded = FALSE,
      path_termination_acceptable = path_termination_acceptable,
      point_converged = path_termination_acceptable,
      finite_validation_loss = finite_validation_loss,
      finite_coefficient = finite_coefficient,
      finite_probability = finite_probability,
      iteration_budget = max_iterations,
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
    fits = fits, coefficient_paths = coefficient_paths,
    preprocess = preprocess, alpha_grid = input$alpha_grid,
    tuning = tuning, selection = selection,
    elapsed_by_alpha = elapsed, runtime_seconds = sum(elapsed),
    warning_messages = warning_messages, tolerance = tolerance,
    max_iterations = max_iterations, irls_tolerance = irls_tolerance,
    irls_max_iterations = irls_max_iterations,
    ridge_boundary_upper_multiplier = ridge_boundary_upper_multiplier,
    null_score_scale = null_score_scale,
    lambda_reference = lambda_reference_values,
    common_lambda_relative_grid = common_relative_grid,
    ridge_lambda_relative_grid = ridge_relative_grid
  )
}

allb_fit_grpreg_common_grid_v2 <- function(
    X, y, group, X_validation, y_validation, nlambda,
    lambda_min_ratio, tolerance = 1e-7, max_iterations = 1000000L,
    penalty_specification = default_grpreg_penalty_specification_v2(),
    common_relative_grid = allb_common_lambda_grid_v2()
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
  required <- c("penalty", "method_path", "penalty_family", "gamma")
  if (!allb_grid_equal_v2(common_relative_grid,
        allb_common_lambda_grid_v2()) ||
      input$nlambda != length(common_relative_grid) ||
      abs(input$lambda_min_ratio - tail(common_relative_grid, 1L)) > 1e-15 ||
      !all(required %in% names(penalty_specification)) ||
      anyDuplicated(penalty_specification$penalty) ||
      anyDuplicated(penalty_specification$method_path) ||
      !setequal(penalty_specification$penalty,
        c("grLasso", "grMCP", "grSCAD")) ||
      max_iterations < 1L || !is.finite(tolerance) || tolerance <= 0) {
    stop("Invalid grpreg V2 common-grid controls.", call. = FALSE)
  }

  group_factor <- factor(input$group, levels = unique(input$group))
  group_multiplier <- sqrt(as.numeric(table(group_factor)))
  fits <- vector("list", nrow(penalty_specification))
  elapsed <- numeric(nrow(penalty_specification))
  warning_messages <- vector("list", nrow(penalty_specification))
  tuning <- vector("list", nrow(penalty_specification))
  lambda_reference_values <- numeric(nrow(penalty_specification))
  base_arguments <- list(
    X = input$X, y = input$y, group = input$group,
    family = "binomial", alpha = 1, eps = tolerance,
    max.iter = max_iterations, dfmax = ncol(input$X),
    gmax = length(unique(input$group)),
    group.multiplier = group_multiplier, warn = TRUE, returnX = FALSE
  )

  for (pi in seq_len(nrow(penalty_specification))) {
    specification <- penalty_specification[pi, , drop = FALSE]
    gamma_for_fit <- if (specification$penalty == "grSCAD") 4 else 3
    captured_warnings <- character(0)
    elapsed[[pi]] <- system.time({
      pilot <- withCallingHandlers(
        do.call(grpreg::grpreg, c(base_arguments, list(
          penalty = specification$penalty,
          nlambda = 3L, lambda.min = 0.8, log.lambda = TRUE,
          gamma = gamma_for_fit
        ))),
        warning = function(condition) {
          captured_warnings <<- c(
            captured_warnings, conditionMessage(condition)
          )
          invokeRestart("muffleWarning")
        }
      )
      lambda_reference <- pilot$lambda[[1L]]
      if (!is.finite(lambda_reference) || lambda_reference <= 0) {
        stop("grpreg did not return a finite native lambda reference.",
          call. = FALSE)
      }
      explicit_lambda <- lambda_reference * common_relative_grid
      fit <- withCallingHandlers(
        do.call(grpreg::grpreg, c(base_arguments, list(
          penalty = specification$penalty,
          lambda = explicit_lambda, log.lambda = TRUE,
          gamma = gamma_for_fit
        ))),
        warning = function(condition) {
          captured_warnings <<- c(
            captured_warnings, conditionMessage(condition)
          )
          invokeRestart("muffleWarning")
        }
      )
      probability <- as.matrix(stats::predict(
        fit, X = input$X_validation, type = "response"
      ))
    })[["elapsed"]]

    lambda_reference_values[[pi]] <- lambda_reference
    path_length <- length(fit$lambda)
    requested_path_length <- length(common_relative_grid)
    returned_relative_grid <- fit$lambda / lambda_reference
    requested_grid_prefix_aligned <-
      path_length <= requested_path_length &&
      allb_grid_equal_v2(
        returned_relative_grid, head(common_relative_grid, path_length)
      )
    lambda_aligned <- path_length == requested_path_length &&
      requested_grid_prefix_aligned
    if (ncol(fit$beta) != path_length || length(fit$iter) != path_length ||
        nrow(probability) != length(input$y_validation) ||
        ncol(probability) != path_length) {
      stop("grpreg V2 returned incompatible coefficient arrays.",
        call. = FALSE)
    }
    path_status <- allb_classify_grpreg_path_v2(
      penalty = specification$penalty[[1L]],
      iterations = fit$iter,
      requested_path_length = requested_path_length,
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
      if (finite_probability[[index]]) {
        binary_log_loss(input$y_validation, probability[, index])
      } else NA_real_
    }, numeric(1))
    finite_validation_loss <- is.finite(validation_loss)
    path_termination_acceptable <-
      path_status$path_termination_acceptable &&
      requested_grid_prefix_aligned
    numerically_eligible <- path_termination_acceptable &
      path_status$point_converged &
      !path_status$returned_lower_boundary_excluded &
      finite_coefficient & finite_probability & finite_validation_loss

    fits[[pi]] <- fit
    warning_messages[[pi]] <- unique(captured_warnings)
    tuning[[pi]] <- data.frame(
      engine = "grpreg", method_path = specification$method_path,
      penalty_family = specification$penalty_family,
      fit_index = pi, alpha_index = 1L, alpha = 1,
      gamma = specification$gamma,
      lambda_index = seq_len(path_length), lambda = fit$lambda,
      lambda_fraction = returned_relative_grid,
      lambda_reference = lambda_reference,
      lambda_reference_type = "native_grpreg_path_start",
      lambda_relative_to_reference = returned_relative_grid,
      alpha_zero_ridge_boundary = FALSE,
      group_selection_capable = TRUE,
      validation_log_loss = validation_loss,
      passes = as.numeric(fit$iter),
      solver_warning = paste(unique(captured_warnings), collapse = " | "),
      unexpected_solver_warning = paste(
        path_status$unexpected_warnings, collapse = " | "
      ),
      requested_path_length = requested_path_length,
      returned_path_length = path_length,
      path_complete = path_status$path_complete && lambda_aligned,
      requested_grid_prefix_aligned = requested_grid_prefix_aligned,
      safe_prefix_length = path_status$safe_prefix_length,
      saturated_path_truncation = path_status$saturated_path_truncation,
      iteration_budget_truncation = path_status$iteration_budget_truncation,
      total_iterations = path_status$total_iterations,
      total_iteration_limit_reached =
        path_status$total_iteration_limit_reached,
      returned_lower_boundary_excluded =
        path_status$returned_lower_boundary_excluded,
      path_termination_acceptable = path_termination_acceptable,
      point_converged = path_status$point_converged,
      finite_validation_loss = finite_validation_loss,
      finite_coefficient = finite_coefficient,
      finite_probability = finite_probability,
      iteration_budget = max_iterations,
      numerically_eligible = numerically_eligible,
      stringsAsFactors = FALSE
    )
  }

  names(fits) <- penalty_specification$penalty
  names(elapsed) <- penalty_specification$penalty
  names(warning_messages) <- penalty_specification$penalty
  names(lambda_reference_values) <- penalty_specification$penalty
  tuning <- do.call(rbind, tuning)
  rownames(tuning) <- NULL
  selections <- vector("list", nrow(penalty_specification))
  names(selections) <- penalty_specification$penalty
  tuning$selected <- FALSE
  tuning$validation_loss_minus_selected <- NA_real_
  tuning$invalid_validation_contender <- FALSE
  tuning$candidate_nonfinite_validation_loss <- FALSE
  tuning$policy_excluded_nonfinite_validation_loss <- FALSE
  for (pi in seq_len(nrow(penalty_specification))) {
    method_path <- penalty_specification$method_path[[pi]]
    selections[[pi]] <- allb_select_grpreg_safe_prefix_v2(
      tuning, method_path
    )
    rows <- tuning$method_path == method_path
    selected <- selections[[pi]]$row
    tuning$selected[rows] <-
      tuning$fit_index[rows] == selected$fit_index[[1L]] &
      tuning$lambda_index[rows] == selected$lambda_index[[1L]]
    tuning$validation_loss_minus_selected[rows] <-
      tuning$validation_log_loss[rows] -
      selected$validation_log_loss[[1L]]
    tuning$invalid_validation_contender[rows] <-
      !tuning$numerically_eligible[rows] &
      tuning$finite_validation_loss[rows] &
      tuning$validation_log_loss[rows] <=
        selected$validation_log_loss[[1L]]
  }
  list(
    fits = fits, penalty_specification = penalty_specification,
    tuning = tuning, selections = selections,
    elapsed_by_penalty = elapsed, runtime_seconds = sum(elapsed),
    warning_messages = warning_messages,
    tolerance = tolerance, max_iterations = max_iterations,
    lambda_reference = lambda_reference_values,
    common_lambda_relative_grid = common_relative_grid
  )
}

allb_finite_fitters_v2 <- function(configuration) {
  expected <- allb_sglasso_lambda_grid_v2()
  common <- allb_common_lambda_grid_v2()
  if (!allb_grid_equal_v2(configuration$lambda_relative_grid, expected) ||
      !identical(as.integer(configuration$nlambda), length(expected)) ||
      !allb_grid_equal_v2(configuration$fair_common_lambda_relative_grid,
        common) ||
      configuration$benchmark_nlambda != length(common) ||
      abs(configuration$lambda_upper_multiplier - expected[[1L]]) > 1e-12 ||
      abs(configuration$lambda_min_ratio - tail(common, 1L)) > 1e-15 ||
      abs(configuration$benchmark_lambda_min_ratio -
        tail(common, 1L)) > 1e-15) {
    stop("ALL V2 requires the frozen method-comparable lambda grids.",
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
  adelie <- function(
      X, y, group, X_validation, y_validation, alpha_grid, nlambda,
      lambda_min_ratio, ridge_boundary_upper_multiplier,
      tolerance, max_iterations, irls_tolerance, irls_max_iterations
  ) {
    allb_fit_adelie_common_grid_v2(
      X, y, group, X_validation, y_validation,
      alpha_grid = alpha_grid, nlambda = nlambda,
      lambda_min_ratio = lambda_min_ratio,
      ridge_boundary_upper_multiplier = ridge_boundary_upper_multiplier,
      tolerance = tolerance, max_iterations = max_iterations,
      irls_tolerance = irls_tolerance,
      irls_max_iterations = irls_max_iterations,
      common_relative_grid = common,
      ridge_relative_grid = expected
    )
  }
  environment(adelie) <- environment()
  fitters <- list(
    sglasso = clone(fit_extended_sglasso_validation_grid_v5),
    adelie = adelie
  )
  # Preserve the frozen V19 production dispatch. The grid changes, but the
  # SGLASSO finite paths still run through the same V11.2 hybrid solver.
  hybrid <- new.env(parent = environment(fitters$sglasso))
  hybrid$fit_logistic_sglasso <- function(
      X, y, group, lambda, d, alpha, max_outer, max_inner,
      tolerance, inner_tolerance, target_original, compile, preprocess,
      use_active_set, warm_start_d, solver
  ) {
    if (!identical(solver, "hybrid_v11_2")) {
      stop("Unexpected ALL V2 SGLASSO solver dispatch.", call. = FALSE)
    }
    fit <- lsg_fit_hybrid_path_v11(
      X, y, group, lambda, d, alpha, target_original,
      preprocess = preprocess, controls = configuration$path_controls,
      compile = FALSE, use_active_set = use_active_set,
      warm_start_d = warm_start_d
    )
    fit$passes <- fit$block_sweeps + fit$apg_iterations
    fit
  }
  environment(hybrid$fit_logistic_sglasso) <- environment()
  environment(fitters$sglasso) <- hybrid
  fitters
}

allb_install_fair_tuning_v2 <- function(e) {
  required <- c(
    "lsg_finite_fitters_v7", "lsg_fit_external_joint_v7",
    "fit_grpreg_validation_paths_v2"
  )
  stopifnot(is.environment(e), all(vapply(required, exists, logical(1),
    envir = e, inherits = FALSE)))
  e$allb_v2_base <- mget(required, envir = e, inherits = FALSE)
  e$lsg_finite_fitters_v7 <- e$allb_finite_fitters_v2

  external <- e$lsg_fit_external_joint_v7
  child <- new.env(parent = environment(external))
  child$fit_grpreg_validation_paths_v2 <- function(
      X, y, group, X_validation, y_validation, nlambda,
      lambda_min_ratio, tolerance, max_iterations,
      penalty_specification
  ) {
    allb_fit_grpreg_common_grid_v2(
      X, y, group, X_validation, y_validation,
      nlambda = nlambda, lambda_min_ratio = lambda_min_ratio,
      tolerance = tolerance, max_iterations = max_iterations,
      penalty_specification = penalty_specification,
      common_relative_grid = allb_common_lambda_grid_v2()
    )
  }
  environment(child$fit_grpreg_validation_paths_v2) <- e
  environment(external) <- child
  e$lsg_fit_external_joint_v7 <- external
  invisible(e)
}
