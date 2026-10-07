# V4 distinguishes a saturated grpreg path from exhaustion of max.iter,
# which grpreg applies to the *whole* path. A budget-limited returned terminal
# can be an unfinished iterate even when its individual iter is below max.iter.

allb_classify_grpreg_path_v4 <- function(
    penalty, iterations, requested_path_length, max_iterations,
    warnings = character(0)
) {
  x <- allb_classify_grpreg_path_v3(
    penalty, iterations, requested_path_length, max_iterations, warnings
  )
  total <- sum(as.numeric(iterations))
  if (total > max_iterations) stop(
    "grpreg returned iterations exceed the declared total path budget.",
    call. = FALSE
  )
  reached <- total >= max_iterations
  saturated <- isTRUE(x$saturated_path_truncation)
  budget_terminal <- reached && !saturated
  warning <- x$unique_warnings
  saturation_warning <- grepl("Model saturated; exiting", warning,
    fixed = TRUE)
  iteration_warning <- grepl(
    "Algorithm failed to converge for all values of lambda", warning,
    fixed = TRUE
  )
  recognized <- (saturation_warning & saturated) |
    (iteration_warning & budget_terminal)
  unexpected <- warning[!recognized]
  excluded <- rep(FALSE, length(iterations))
  converged <- as.numeric(iterations) < max_iterations
  if (budget_terminal) {
    excluded[[length(excluded)]] <- TRUE
    converged[[length(converged)]] <- FALSE
  }
  x$iteration_budget_truncation <- budget_terminal
  x$total_iterations <- total
  x$total_iteration_limit_reached <- reached
  x$returned_lower_boundary_excluded <- excluded
  x$point_converged <- converged
  x$path_termination_acceptable <- !length(unexpected) &&
    (isTRUE(x$path_complete) || saturated || budget_terminal)
  x$safe_prefix_length <- as.integer(sum(converged & !excluded))
  x$unexpected_warnings <- unexpected
  x
}

allb_select_grpreg_returned_path_v4 <- function(tuning, method_path) {
  candidates <- tuning[tuning$method_path == method_path, , drop = FALSE]
  if (!nrow(candidates)) stop("Missing grpreg candidate path.", call. = FALSE)
  candidates <- candidates[order(candidates$lambda_index), , drop = FALSE]
  budget <- candidates$iteration_budget_truncation[[1L]] %in% TRUE
  if (!budget) return(allb_select_grpreg_returned_path_v3(
    tuning, method_path
  ))
  n <- nrow(candidates)
  expected <- c(rep(FALSE, n - 1L), TRUE)
  structural <- n >= 2L &&
    identical(as.integer(candidates$lambda_index), seq_len(n)) &&
    all(candidates$requested_grid_prefix_aligned %in% TRUE) &&
    all(candidates$path_termination_acceptable %in% TRUE) &&
    all(candidates$iteration_budget_truncation %in% TRUE) &&
    all(candidates$returned_path_length == n) &&
    all(candidates$requested_path_length ==
      candidates$requested_path_length[[1L]]) &&
    all(candidates$safe_prefix_length == n - 1L) &&
    identical(as.logical(candidates$returned_lower_boundary_excluded),
      expected) &&
    identical(as.logical(candidates$point_converged), !expected) &&
    all(!is.na(candidates$unexpected_solver_warning) &
      !nzchar(candidates$unexpected_solver_warning))
  if (!isTRUE(structural)) stop(
    "Invalid ALL V4 grpreg budget-prefix evidence for ", method_path,
    ".", call. = FALSE
  )
  eligible <- candidates$point_converged %in% TRUE &
    candidates$returned_lower_boundary_excluded %in% FALSE &
    candidates$finite_validation_loss %in% TRUE &
    candidates$finite_coefficient %in% TRUE &
    candidates$finite_probability %in% TRUE
  if (!identical(as.logical(candidates$numerically_eligible), eligible)) {
    stop("Inconsistent ALL V4 grpreg eligibility flags.", call. = FALSE)
  }
  # The V3 selector already computes the validation audit and selection. Its
  # structural check requires all exclusion flags FALSE, so pass a private
  # copy after independently checking the exact V4 terminal-exclusion schema.
  shadow <- tuning
  rows <- shadow$method_path == method_path
  shadow$returned_lower_boundary_excluded[rows] <- FALSE
  result <- allb_select_grpreg_returned_path_v3(shadow, method_path)
  result$selected_on_returned_lower_boundary <- FALSE
  result$selected_on_solver_lower_boundary <- FALSE
  result$truncated_path_selection_interior <- TRUE
  result
}

allb_fit_grpreg_common_grid_v4 <- function(...) {
  fitter <- allb_fit_grpreg_common_grid_v2
  child <- new.env(parent = environment(fitter))
  child$allb_classify_grpreg_path_v2 <- allb_classify_grpreg_path_v4
  child$allb_select_grpreg_safe_prefix_v2 <-
    allb_select_grpreg_returned_path_v4
  environment(fitter) <- child
  result <- fitter(...)
  result$returned_path_policy <-
    "saturation_all_returned_eligible_budget_terminal_excluded"
  result
}

allb_install_fair_tuning_v4 <- function(e) {
  allb_install_fair_tuning_v3(e)
  external <- e$lsg_fit_external_joint_v7
  child <- new.env(parent = environment(external))
  child$fit_grpreg_validation_paths_v2 <- function(
      X, y, group, X_validation, y_validation, nlambda,
      lambda_min_ratio, tolerance, max_iterations,
      penalty_specification
  ) {
    allb_fit_grpreg_common_grid_v4(
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
