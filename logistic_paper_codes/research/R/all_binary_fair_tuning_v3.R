# Returned-path policy for the ALL BCR/ABL versus NEG V3 study.
#
# grpreg removes coefficient columns whose iteration count is NA before it
# returns a fitted object.  Consequently, a truncated path is a shorter
# solver-defined prefix of the requested grid; its last returned point is not
# itself evidence of saturation.  V3 keeps every finite returned point whose
# iteration count is below max.iter and records, rather than rejects, selection
# on the solver-defined lower boundary.

allb_classify_grpreg_path_v3 <- function(
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
    stop("Invalid ALL V3 grpreg path-classification inputs.",
      call. = FALSE)
  }

  returned <- length(iterations)
  complete <- returned == requested_path_length
  saturation_warning <- grepl(
    "Model saturated; exiting", warnings, fixed = TRUE
  )
  iteration_warning <- grepl(
    "Algorithm failed to converge for all values of lambda",
    warnings, fixed = TRUE
  )
  recognized <- saturation_warning | iteration_warning
  unexpected <- warnings[!recognized]
  saturated <- !complete && any(saturation_warning)
  iteration_truncated <- any(iteration_warning) ||
    any(iterations >= max_iterations)
  acceptable <- !length(unexpected) &&
    (complete || saturated || iteration_truncated)
  point_converged <- iterations < max_iterations

  list(
    path_complete = complete,
    saturated_path_truncation = saturated,
    iteration_budget_truncation = iteration_truncated,
    total_iterations = sum(iterations),
    total_iteration_limit_reached = any(iterations >= max_iterations),
    returned_lower_boundary_excluded = rep(FALSE, returned),
    point_converged = point_converged,
    path_termination_acceptable = acceptable,
    safe_prefix_length = as.integer(sum(point_converged)),
    returned_path_length = as.integer(returned),
    unique_warnings = warnings,
    unexpected_warnings = unexpected
  )
}

allb_select_grpreg_returned_path_v3 <- function(tuning, method_path) {
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
    stop("ALL V3 grpreg returned-path metadata are incomplete.",
      call. = FALSE)
  }
  candidates <- candidates[order(candidates$lambda_index), , drop = FALSE]
  n <- nrow(candidates)
  truncated <- candidates$saturated_path_truncation[[1L]] %in% TRUE ||
    (!candidates$path_complete[[1L]])
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
    all(candidates$safe_prefix_length == sum(candidates$point_converged)) &&
    !any(candidates$returned_lower_boundary_excluded %in% TRUE) &&
    all(!is.na(candidates$unexpected_solver_warning) &
      !nzchar(candidates$unexpected_solver_warning)) &&
    identical(as.logical(candidates$numerically_eligible), eligible)
  if (!isTRUE(structural)) {
    stop("Invalid ALL V3 grpreg returned-path evidence for ", method_path,
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
  selected_on_solver_boundary <- isTRUE(truncated) &&
    selected$lambda_index[[1L]] == n

  list(
    row = selected,
    total_candidate_points = n,
    numerically_eligible_points = nrow(valid),
    numerically_ineligible_points = n - nrow(valid),
    nonfinite_validation_loss_points =
      sum(!candidates$finite_validation_loss),
    candidate_nonfinite_validation_loss_points =
      sum(!candidates$finite_validation_loss),
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
    selected_on_returned_lower_boundary = selected_on_solver_boundary,
    selected_on_solver_lower_boundary = selected_on_solver_boundary,
    returned_lower_boundary_validation_log_loss =
      candidates$validation_log_loss[[n]],
    returned_lower_boundary_minus_selected_validation_log_loss =
      candidates$validation_log_loss[[n]] -
        selected$validation_log_loss[[1L]],
    truncated_path_selection_interior = !selected_on_solver_boundary
  )
}

allb_fit_grpreg_common_grid_v3 <- function(...) {
  fitter <- allb_fit_grpreg_common_grid_v2
  child <- new.env(parent = environment(fitter))
  child$allb_classify_grpreg_path_v2 <- allb_classify_grpreg_path_v3
  child$allb_select_grpreg_safe_prefix_v2 <-
    allb_select_grpreg_returned_path_v3
  environment(fitter) <- child
  result <- fitter(...)
  result$returned_path_policy <-
    "all_finite_returned_points_below_max_iter_are_eligible"
  result
}

allb_install_fair_tuning_v3 <- function(e) {
  required <- c("lsg_fit_external_joint_v7", "allb_v2_base")
  stopifnot(is.environment(e), all(vapply(required, exists, logical(1),
    envir = e, inherits = FALSE)))
  external <- e$lsg_fit_external_joint_v7
  child <- new.env(parent = environment(external))
  child$fit_grpreg_validation_paths_v2 <- function(
      X, y, group, X_validation, y_validation, nlambda,
      lambda_min_ratio, tolerance, max_iterations,
      penalty_specification
  ) {
    allb_fit_grpreg_common_grid_v3(
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
