# V18 production repair overlay. Frozen V15 functions are never rewritten.

lsg_lambda_relative_grid_v18 <- function() {
  c(2048, 1024, 512, 256, 128, 64, 32, 16,
    exp(seq(log(8), log(0.05), length.out = 30L)))
}

lsg_configuration_v18 <- function(stage) {
  cfg <- v18_base$lsg_configuration_v7(stage)
  grid <- lsg_lambda_relative_grid_v18()
  cfg$nlambda <- length(grid)
  cfg$lambda_relative_grid <- grid
  cfg$lambda_extension_multipliers <- grid[seq_len(8L)]
  cfg$lambda_upper_multiplier <- grid[1L]
  cfg$integration_version <- "v18_fixed_2048_tail_complete_grlasso_prefix"
  cfg$group_lasso_prefix_policy <-
    "complete_path_final_total_budget_point_excluded_v18"
  cfg
}

lsg_finite_fitters_v18 <- function(configuration) {
  expected <- lsg_lambda_relative_grid_v18()
  if (!isTRUE(all.equal(as.numeric(configuration$lambda_relative_grid),
                        expected, tolerance = 1e-13)) ||
      !identical(as.integer(configuration$nlambda), length(expected)) ||
      !isTRUE(all.equal(as.numeric(configuration$lambda_upper_multiplier), 2048)) ||
      !isTRUE(all.equal(as.numeric(configuration$lambda_min_ratio), 0.05))) {
    stop("V18 requires the frozen 38-point finite lambda grid.", call. = FALSE)
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
  fitters <- list(
    sglasso = clone(fit_extended_sglasso_validation_grid_v5),
    adelie = clone(fit_adelie_validation_grid_v5)
  )

  # The frozen V5 grid builder must still dispatch each finite SGLASSO path
  # through the V12 hybrid solver. The first child environment above supplies
  # only the V18 grid; this second layer preserves the validated V12 dispatch.
  hybrid <- new.env(parent = environment(fitters$sglasso))
  hybrid$fit_logistic_sglasso <- function(X, y, group, lambda, d, alpha,
      max_outer, max_inner, tolerance, inner_tolerance, target_original,
      compile, preprocess, use_active_set, warm_start_d, solver) {
    lsg_assert_v7(
      identical(solver, "hybrid_v11_2"),
      "Unexpected V18 integration solver."
    )
    fit <- lsg_fit_hybrid_path_v11(
      X, y, group, lambda, d, alpha, target_original,
      preprocess = preprocess, controls = configuration$path_controls,
      compile = FALSE, use_active_set = use_active_set,
      warm_start_d = warm_start_d
    )
    # The legacy column records solver work, not an ABGD-only pass count.
    fit$passes <- fit$block_sweeps + fit$apg_iterations
    fit
  }
  environment(fitters$sglasso) <- hybrid
  fitters
}

lsg_selection_boundary_v18 <- function(selected, configuration) {
  finite <- identical(selected$point_type[[1L]], "finite")
  reference_type <- selected$lambda_reference_type[[1L]]
  extended <- reference_type %in% c("null_score_ridge_boundary", "d0_null_kkt")
  ratio <- selected$lambda_relative_to_reference[[1L]]
  upper <- finite && extended && is.finite(ratio) &&
    abs(ratio - configuration$lambda_upper_multiplier) <= 1e-10
  lower <- finite && extended && is.finite(ratio) &&
    abs(ratio - configuration$lambda_min_ratio) <= 1e-10
  list(
    upper = upper, lower = lower,
    native_upper = finite && !extended && is.finite(ratio) &&
      abs(ratio - 1) <= 1e-10,
    native_lower = finite && !extended && is.finite(ratio) &&
      abs(ratio - configuration$lambda_min_ratio) <= 1e-10,
    unresolved = upper || lower
  )
}

# Preserve every V9 classification except the exact grLasso case where a full
# requested path is returned and its cumulative iteration budget is first met
# at the final point. Only that final point is excluded.
lsg_classify_grpreg_path_v18 <- function(
    penalty, iterations, requested_path_length, max_iterations,
    warnings = character(0)
) {
  original <- v18_base$lsg_classify_grpreg_path_v9(
    penalty, iterations, requested_path_length, max_iterations, warnings
  )
  iterations <- as.numeric(iterations)
  complete_final_budget <- identical(penalty, "grLasso") &&
    length(iterations) == as.integer(requested_path_length) &&
    sum(iterations) == as.numeric(max_iterations) &&
    all(head(cumsum(iterations), -1L) < as.numeric(max_iterations))
  if (!complete_final_budget) return(original)

  unique_warnings <- unique(as.character(warnings))
  unique_warnings <- unique_warnings[
    !is.na(unique_warnings) & nzchar(unique_warnings)]
  iteration_warning <- grepl(
    "Algorithm failed to converge for all values of lambda",
    unique_warnings, fixed = TRUE
  )
  unexpected <- unique_warnings[!iteration_warning]
  returned <- length(iterations)
  excluded <- rep(FALSE, returned)
  excluded[returned] <- TRUE
  point_converged <- iterations < as.numeric(max_iterations)
  point_converged[returned] <- FALSE
  list(
    path_complete = TRUE,
    saturated_path_truncation = FALSE,
    iteration_budget_truncation = TRUE,
    group_lasso_safe_prefix = TRUE,
    total_iterations = sum(iterations),
    total_iteration_limit_reached = TRUE,
    returned_lower_boundary_excluded = excluded,
    point_converged = point_converged,
    path_termination_acceptable = !length(unexpected),
    usable_prefix_length = as.integer(returned - 1L),
    unique_warnings = unique_warnings,
    unexpected_warnings = unexpected
  )
}

lsg_patch_grlasso_tuning_v18 <- function(tuning, status) {
  n <- nrow(tuning)
  stopifnot(n == length(status$point_converged))
  scalar_fields <- c(
    "path_complete", "saturated_path_truncation",
    "iteration_budget_truncation", "total_iteration_limit_reached",
    "path_termination_acceptable"
  )
  for (field in scalar_fields) tuning[[field]] <- rep(status[[field]], n)
  tuning$returned_lower_boundary_excluded <-
    status$returned_lower_boundary_excluded
  tuning$point_converged <- status$point_converged
  tuning$total_iterations <- rep(status$total_iterations, n)
  tuning$unexpected_solver_warning <-
    rep(paste(status$unexpected_warnings, collapse = " | "), n)
  tuning$numerically_eligible <-
    status$path_termination_acceptable & status$point_converged &
    !status$returned_lower_boundary_excluded &
    tuning$finite_coefficient & tuning$finite_probability &
    is.finite(tuning$validation_log_loss)
  tuning
}

lsg_tuning_unit_checks_v18 <- function() {
  check <- function(name, passed) data.frame(
    check = name, passed = isTRUE(passed), stringsAsFactors = FALSE
  )
  grid <- lsg_lambda_relative_grid_v18()
  configuration <- list(
    lambda_relative_grid = grid, nlambda = length(grid),
    lambda_upper_multiplier = 2048, lambda_min_ratio = 0.05
  )
  old_sg <- fit_extended_sglasso_validation_grid_v5
  old_ad <- fit_adelie_validation_grid_v5
  old_options <- options()
  fitters <- lsg_finite_fitters_v18(configuration)
  checks <- list(
    check("v18_cloned_finite_fitters_preserve_bodies",
      identical(body(fitters$sglasso), body(old_sg)) &&
        identical(body(fitters$adelie), body(old_ad))),
    check("v18_no_global_frozen_binding_or_option_changes",
      identical(old_options, options()) &&
        identical(old_sg, fit_extended_sglasso_validation_grid_v5) &&
        identical(old_ad, fit_adelie_validation_grid_v5)),
    check("v18_isolated_38_point_grid",
      identical(get("logistic_sglasso_lambda_relative_grid_v5",
                    envir = environment(fitters$sglasso),
                    inherits = TRUE)(), grid) &&
        identical(get("logistic_sglasso_lambda_relative_grid_v5",
                      envir = environment(fitters$adelie),
                      inherits = TRUE)(), grid) &&
        identical(grid[9:38],
                  exp(seq(log(8), log(0.05), length.out = 30L)))),
    check("v18_hybrid_dispatch_layer_present",
      exists("fit_logistic_sglasso", envir = environment(fitters$sglasso),
             inherits = FALSE) &&
        grepl("lsg_fit_hybrid_path_v11", paste(deparse(body(get(
          "fit_logistic_sglasso", envir = environment(fitters$sglasso),
          inherits = FALSE))), collapse = "\n"), fixed = TRUE))
  )
  toy <- data.frame(
    candidate_id = c("upper", "old_upper", "endpoint"),
    alpha = 0, d = 0, point_type = c("finite", "finite", "penalty_limit"),
    lambda_relative_to_reference = c(2048, 512, Inf),
    lambda_reference_type = "null_score_ridge_boundary",
    validation_log_loss = c(0.5, 0.6, 0.7), numerically_eligible = TRUE,
    stringsAsFactors = FALSE
  )
  checks[[length(checks) + 1L]] <- check(
    "v18_only_2048_is_finite_upper_boundary",
    lsg_selection_boundary_v18(toy[1L, ], configuration)$unresolved &&
      !lsg_selection_boundary_v18(toy[2L, ], configuration)$unresolved &&
      !lsg_selection_boundary_v18(toy[3L, ], configuration)$unresolved
  )
  complete <- lsg_classify_grpreg_path_v18(
    "grLasso", c(100, 200, 300, 999400), 4L, 1000000L
  )
  checks[[length(checks) + 1L]] <- check(
    "v18_complete_final_budget_path_has_safe_prefix",
    complete$path_complete && complete$iteration_budget_truncation &&
      complete$group_lasso_safe_prefix &&
      complete$path_termination_acceptable &&
      identical(complete$returned_lower_boundary_excluded,
                c(FALSE, FALSE, FALSE, TRUE)) &&
      identical(complete$usable_prefix_length, 3L)
  )
  unknown <- lsg_classify_grpreg_path_v18(
    "grLasso", c(100, 200, 300, 999400), 4L, 1000000L,
    "unknown warning"
  )
  checks[[length(checks) + 1L]] <- check(
    "v18_unknown_warning_still_rejected",
    !unknown$path_termination_acceptable &&
      identical(unknown$unexpected_warnings, "unknown warning")
  )
  old_nonconvex <- v18_base$lsg_classify_grpreg_path_v9(
    "grMCP", c(10, 20, 30), 5L, 1000L, "Model saturated; exiting"
  )
  new_nonconvex <- lsg_classify_grpreg_path_v18(
    "grMCP", c(10, 20, 30), 5L, 1000L, "Model saturated; exiting"
  )
  checks[[length(checks) + 1L]] <- check(
    "v18_nonconvex_classification_unchanged",
    identical(old_nonconvex, new_nonconvex)
  )
  bad <- configuration
  bad$lambda_relative_grid[1L] <- 4096
  checks[[length(checks) + 1L]] <- check(
    "v18_unapproved_grid_rejected",
    inherits(tryCatch(lsg_finite_fitters_v18(bad), error = identity), "error")
  )
  do.call(rbind, checks)
}
