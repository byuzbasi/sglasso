# Versioned numerical repair for Logistic SGLASSO.
#
# V7/V8 functions and source artifacts remain immutable.  This module changes
# only how the already-defined convex Logistic SGLASSO objective is solved and
# how a cumulative-budget Group Lasso prefix is classified.

lsg_lambda_relative_grid_v9 <- function() {
  c(4096, 2048, 1024, 512, 256, 128, 64, 32, 16,
    exp(seq(log(8), log(0.05), length.out = 30L)))
}


lsg_compile_joint_solver_v9 <- function(root = lsg_prework_root(),
                                        rebuild = FALSE,
                                        quiet = TRUE) {
  required <- c("lsg_fit_one_joint_v9_cpp", "lsg_path_joint_v9_cpp")
  envir <- environment(lsg_compile_joint_solver_v9)
  if (!isTRUE(rebuild) && all(vapply(
    required, exists, logical(1), mode = "function", envir = envir,
    inherits = TRUE
  ))) return(invisible(TRUE))
  if (!requireNamespace("Rcpp", quietly = TRUE) ||
      !requireNamespace("RcppArmadillo", quietly = TRUE)) {
    stop("Rcpp and RcppArmadillo are required; no installation is attempted.",
         call. = FALSE)
  }
  root <- normalizePath(root, mustWork = TRUE)
  source <- file.path(root, "src", "logistic_sglasso_joint_solver_v9.cpp")
  if (!file.exists(source) || nzchar(Sys.readlink(source))) {
    stop("Missing regular V9 joint-solver source file.", call. = FALSE)
  }

  local_makevars <- file.path(root, "config", "Makevars.local")
  local_gfortran_runtime <- paste0(
    "/usr/local/gfortran/lib/gcc/",
    "aarch64-apple-darwin23/14.1.0/libemutls_w.a"
  )
  previous_makevars <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  if (file.exists(local_makevars) && file.exists(local_gfortran_runtime)) {
    Sys.setenv(R_MAKEVARS_USER = local_makevars)
  }
  on.exit({
    if (is.na(previous_makevars)) Sys.unsetenv("R_MAKEVARS_USER")
    else Sys.setenv(R_MAKEVARS_USER = previous_makevars)
  }, add = TRUE)

  Rcpp::sourceCpp(
    source, rebuild = isTRUE(rebuild), showOutput = !isTRUE(quiet),
    verbose = FALSE, env = envir
  )
  invisible(TRUE)
}


lsg_validate_joint_controls_v9 <- function(max_sweeps,
                                           kkt_tolerance,
                                           update_tolerance,
                                           intercept_tolerance,
                                           max_intercept_iterations) {
  integer_control <- function(value, name) {
    if (length(value) != 1L || is.na(value) || !is.finite(value) ||
        value < 1 || value != as.integer(value)) {
      stop(name, " must be one positive integer.", call. = FALSE)
    }
    as.integer(value)
  }
  positive_control <- function(value, name) {
    if (length(value) != 1L || is.na(value) || !is.finite(value) || value <= 0) {
      stop(name, " must be one finite positive value.", call. = FALSE)
    }
    as.numeric(value)
  }
  list(
    max_sweeps = integer_control(max_sweeps, "max_sweeps"),
    kkt_tolerance = positive_control(kkt_tolerance, "kkt_tolerance"),
    update_tolerance = positive_control(update_tolerance, "update_tolerance"),
    intercept_tolerance = positive_control(
      intercept_tolerance, "intercept_tolerance"
    ),
    max_intercept_iterations = integer_control(
      max_intercept_iterations, "max_intercept_iterations"
    )
  )
}


lsg_fit_joint_path_v9 <- function(
    X,
    y,
    group,
    lambda,
    d,
    alpha,
    target_original,
    preprocess = NULL,
    max_sweeps = 4000L,
    kkt_tolerance = 2e-6,
    update_tolerance = 1e-10,
    intercept_tolerance = 1e-12,
    max_intercept_iterations = 100L,
    use_active_set = TRUE,
    warm_start_d = TRUE,
    keep_traces = FALSE,
    lambda_order = c("decreasing", "any"),
    compile = TRUE,
    root = lsg_prework_root()
) {
  lambda_order <- match.arg(lambda_order)
  controls <- lsg_validate_joint_controls_v9(
    max_sweeps, kkt_tolerance, update_tolerance, intercept_tolerance,
    max_intercept_iterations
  )
  if (isTRUE(compile)) lsg_compile_joint_solver_v9(root)
  response <- normalize_binary_response(y)$y
  X <- as.matrix(X)
  if (nrow(X) != length(response)) {
    stop("X and y have incompatible dimensions.", call. = FALSE)
  }
  alpha <- as.numeric(alpha)
  d <- as.numeric(d)
  lambda <- as.numeric(lambda)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha < 0 || alpha > 1 ||
      !length(d) || any(!is.finite(d)) || any(d < 0 | d > 1) ||
      !length(lambda) || any(!is.finite(lambda)) || any(lambda <= 0) ||
      (lambda_order == "decreasing" && is.unsorted(-lambda, strictly = FALSE))) {
    stop("Invalid V9 alpha, d, or finite lambda path.", call. = FALSE)
  }
  if (length(use_active_set) != 1L || is.na(use_active_set) ||
      length(warm_start_d) != 1L || is.na(warm_start_d) ||
      length(keep_traces) != 1L || is.na(keep_traces)) {
    stop("V9 logical solver controls must be scalar and nonmissing.",
         call. = FALSE)
  }
  target_original <- as.numeric(target_original)
  if (length(target_original) != ncol(X) || any(!is.finite(target_original))) {
    stop("target_original must be finite and match the original predictors.",
         call. = FALSE)
  }
  if (is.null(preprocess)) {
    preprocess <- prepare_lsg_design(X, group)
  } else if (!inherits(preprocess, "lsg_preprocess") ||
             preprocess$n != nrow(X) || preprocess$p_original != ncol(X)) {
    stop("preprocess is incompatible with X.", call. = FALSE)
  }
  target <- project_lsg_original_target(preprocess, target_original)
  core <- lsg_path_joint_v9_cpp(
    preprocess$X, response, preprocess$group_start, preprocess$group_end,
    preprocess$group_weight, target, lambda, d, alpha,
    controls$max_sweeps, controls$kkt_tolerance,
    controls$update_tolerance, controls$intercept_tolerance,
    controls$max_intercept_iterations, isTRUE(use_active_set),
    isTRUE(warm_start_d), isTRUE(keep_traces)
  )

  coefficient_array <- array(
    0,
    dim = c(preprocess$p_original + 1L, length(lambda), length(d)),
    dimnames = list(c("(Intercept)", preprocess$names),
                    signif(lambda, 5), d)
  )
  for (di in seq_along(d)) {
    for (li in seq_along(lambda)) {
      coefficient_array[, li, di] <- recover_lsg_coefficients(
        preprocess, core$beta[, li, di], core$intercept[li, di]
      )
    }
  }

  selected_groups <- apply(core$beta, c(2L, 3L), function(value) {
    sum(vapply(seq_along(preprocess$group_start), function(g) {
      first <- preprocess$group_start[g] + 1L
      last <- preprocess$group_end[g] + 1L
      sqrt(sum(value[first:last]^2)) > 1e-8
    }, logical(1)))
  })
  out <- list(
    coefficients = coefficient_array,
    beta_solver = core$beta,
    intercept_solver = core$intercept,
    lambda = lambda,
    d = d,
    alpha = alpha,
    solver = core$solver,
    target = target,
    target_source = "supplied_original_coefficients",
    target_information = list(
      target = target, source = "supplied_original_coefficients",
      original_target = target_original
    ),
    preprocess = preprocess,
    y = response,
    objective = core$objective,
    kkt = core$kkt,
    intercept_kkt = core$intercept_kkt,
    group_kkt = core$group_kkt,
    converged = core$converged != 0,
    outer_iterations = core$sweeps,
    inner_iterations = matrix(0L, nrow = length(lambda), ncol = length(d)),
    passes = core$sweeps,
    sweeps = core$sweeps,
    selected_groups = selected_groups,
    termination_reason = core$termination_reason,
    intercept_iterations = core$intercept_iterations,
    backtracking_steps = core$backtracking_steps,
    kkt_scans = core$kkt_scans,
    group_updates = core$group_updates,
    active_groups = core$active_groups,
    group_curvature = core$group_curvature_upper_bound,
    maximum_raw_objective_increase = core$max_raw_objective_increase,
    maximum_accepted_objective_increase =
      core$max_accepted_objective_increase,
    objective_traces = core$objective_traces,
    use_active_set = isTRUE(use_active_set),
    warm_start_d = isTRUE(warm_start_d),
    kkt_tolerance = controls$kkt_tolerance,
    intercept_tolerance = controls$intercept_tolerance,
    call = match.call()
  )
  class(out) <- c("logistic_sglasso_prework_v9",
                  "logistic_sglasso_prework")
  out
}


lsg_fit_joint_validation_grid_v9 <- function(
    X,
    y,
    group,
    X_validation,
    y_validation,
    target_original,
    configuration
) {
  expected_grid <- lsg_lambda_relative_grid_v9()
  alpha_grid <- sort(unique(as.numeric(configuration$alpha_grid)))
  d_grid <- sort(unique(as.numeric(configuration$d_grid)))
  if (!isTRUE(all.equal(as.numeric(configuration$lambda_relative_grid),
                        expected_grid, tolerance = 1e-13)) ||
      as.integer(configuration$nlambda) != length(expected_grid) ||
      !identical(alpha_grid, seq(0, 1, by = 0.1)) ||
      !identical(d_grid, seq(0, 1, by = 0.1))) {
    stop("V9 requires the frozen 39-point lambda and 0.1 alpha/d grids.",
         call. = FALSE)
  }
  response <- normalize_binary_response(y)$y
  preprocess <- prepare_lsg_design(X, group)
  target <- project_lsg_original_target(preprocess, target_original)
  null_score_scale <- lsg_null_score_lambda_scale(preprocess, response)
  fits <- vector("list", length(alpha_grid))
  loss <- vector("list", length(alpha_grid))
  elapsed <- numeric(length(alpha_grid))
  lambda_reference <- numeric(length(alpha_grid))
  lambda_d0_kkt_reference <- rep(NA_real_, length(alpha_grid))
  lambda_reference_type <- character(length(alpha_grid))

  for (ai in seq_along(alpha_grid)) {
    alpha <- alpha_grid[ai]
    fitted_d <- if (abs(alpha - 1) <= 1e-12) 0 else d_grid
    if (abs(alpha) <= 1e-12) {
      lambda_reference[ai] <- null_score_scale
      lambda_reference_type[ai] <- "null_score_ridge_boundary"
    } else {
      diagnostic <- lsg_lambda_start_cpp(
        preprocess$X, response, preprocess$group_start,
        preprocess$group_end, preprocess$group_weight, target, alpha, 0
      )
      if (!isTRUE(diagnostic$zero_model_feasible) ||
          !is.finite(diagnostic$lambda_start) || diagnostic$lambda_start <= 0) {
        stop("Unable to construct the V9 d=0 KKT lambda reference.",
             call. = FALSE)
      }
      lambda_d0_kkt_reference[ai] <-
        diagnostic$lambda_start * (1 + 1e-8)
      lambda_reference[ai] <- lambda_d0_kkt_reference[ai]
      lambda_reference_type[ai] <- "d0_null_kkt"
    }
    lambda <- lambda_reference[ai] * expected_grid
    elapsed[ai] <- system.time({
      fits[[ai]] <- lsg_fit_joint_path_v9(
        X, y, group, lambda, fitted_d, alpha, target_original,
        preprocess = preprocess,
        max_sweeps = configuration$joint_max_sweeps,
        kkt_tolerance = configuration$joint_kkt_tolerance,
        update_tolerance = configuration$joint_update_tolerance,
        intercept_tolerance = configuration$joint_intercept_tolerance,
        max_intercept_iterations =
          configuration$joint_max_intercept_iterations,
        use_active_set = TRUE, warm_start_d = TRUE,
        keep_traces = FALSE, compile = FALSE
      )
      probability <- predict_logistic_sglasso(
        fits[[ai]], X_validation, type = "response"
      )
      loss[[ai]] <- aligned_validation_loss(y_validation, probability)
    })[["elapsed"]]
    fits[[ai]]$lambda_reference <- lambda_reference[ai]
    fits[[ai]]$lambda_reference_type <- lambda_reference_type[ai]
    fits[[ai]]$lambda_d0_kkt_reference <- lambda_d0_kkt_reference[ai]
    fits[[ai]]$lambda_relative_to_reference <- expected_grid
    fits[[ai]]$lambda_relative_to_d0_kkt <-
      fits[[ai]]$lambda / lambda_d0_kkt_reference[ai]
  }
  list(
    fits = fits, validation_loss = loss, alpha_grid = alpha_grid,
    d_grid = d_grid, elapsed_by_alpha = elapsed,
    elapsed_seconds = sum(elapsed), lambda_reference = lambda_reference,
    lambda_reference_type = lambda_reference_type,
    lambda_d0_kkt_reference = lambda_d0_kkt_reference,
    lambda_upper_multiplier = expected_grid[1L],
    lambda_min_reference_fraction = tail(expected_grid, 1L),
    lambda_relative_grid = expected_grid
  )
}


# V9 recognizes a cumulative-budget Group Lasso truncation but excludes the
# budget-ending returned point. Unknown warnings and malformed paths fail.
lsg_classify_grpreg_path_v9 <- function(
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
      length(requested_path_length) != 1L || requested_path_length < 2L ||
      !length(iterations) || length(iterations) > requested_path_length ||
      any(!is.finite(iterations)) || any(iterations < 0) ||
      length(max_iterations) != 1L || !is.finite(max_iterations) ||
      max_iterations < 1) {
    stop("Invalid V9 grpreg path-classification inputs.", call. = FALSE)
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
  nonconvex <- penalty %in% c("grMCP", "grSCAD")
  saturated <- nonconvex && !complete && any(saturation_warning)
  budget_truncated <- !complete && !saturated && budget_reached
  group_lasso_safe_prefix <- penalty == "grLasso" && budget_truncated
  recognized <- (saturation_warning & saturated) |
    (iteration_warning & budget_truncated)
  unexpected <- warnings[!recognized]
  acceptable <- !length(unexpected) &&
    ((complete && !budget_reached) || saturated || budget_truncated)
  excluded <- rep(FALSE, returned)
  if (saturated || budget_truncated) excluded[returned] <- TRUE
  point_converged <- iterations < max_iterations
  if (saturated || budget_truncated) point_converged[returned] <- FALSE
  usable_prefix <- if (saturated || budget_truncated) returned - 1L else returned
  list(
    path_complete = complete,
    saturated_path_truncation = saturated,
    iteration_budget_truncation = budget_truncated,
    group_lasso_safe_prefix = group_lasso_safe_prefix,
    total_iterations = total,
    total_iteration_limit_reached = budget_reached,
    returned_lower_boundary_excluded = excluded,
    point_converged = point_converged,
    path_termination_acceptable = acceptable,
    usable_prefix_length = as.integer(usable_prefix),
    unique_warnings = warnings,
    unexpected_warnings = unexpected
  )
}


lsg_select_grpreg_prefix_v9 <- function(validation_log_loss,
                                        finite_coefficient,
                                        finite_probability,
                                        path_status) {
  loss <- as.numeric(validation_log_loss)
  finite_coefficient <- as.logical(finite_coefficient)
  finite_probability <- as.logical(finite_probability)
  n <- length(loss)
  if (length(finite_coefficient) != n || length(finite_probability) != n ||
      length(path_status$point_converged) != n ||
      length(path_status$returned_lower_boundary_excluded) != n) {
    stop("V9 Group Lasso selection vectors have incompatible lengths.",
         call. = FALSE)
  }
  eligible <- isTRUE(path_status$path_termination_acceptable) &
    path_status$point_converged &
    !path_status$returned_lower_boundary_excluded &
    finite_coefficient & finite_probability & is.finite(loss)
  if (!any(eligible)) {
    stop("No numerically eligible V9 Group Lasso prefix point.",
         call. = FALSE)
  }
  selected <- which(eligible)[which.min(loss[eligible])]
  lower_interior <- !isTRUE(path_status$group_lasso_safe_prefix) ||
    selected < path_status$usable_prefix_length
  list(
    selected_index = as.integer(selected),
    selected_validation_log_loss = loss[selected],
    numerically_eligible = eligible,
    selected_strictly_above_usable_lower_boundary = lower_interior,
    accepted = isTRUE(lower_interior)
  )
}


lsg_portable_numeric_equal_v9 <- function(x, y, tolerance = 1e-14) {
  if (length(tolerance) != 1L || is.na(tolerance) ||
      !is.finite(tolerance) || tolerance < 0) {
    stop("Portable numeric tolerance must be one finite nonnegative value.",
         call. = FALSE)
  }
  isTRUE(all.equal(
    x, y, tolerance = tolerance, check.attributes = FALSE
  ))
}


lsg_joint_solver_unit_checks_v9 <- function(root = lsg_prework_root()) {
  compile_lsg_core(rebuild = FALSE, quiet = TRUE)
  lsg_compile_joint_solver_v9(root, rebuild = FALSE, quiet = TRUE)
  set.seed(9010901L)
  n <- 90L
  groups <- rep(seq_len(4L), each = 2L)
  X <- matrix(stats::rnorm(n * length(groups)), nrow = n)
  X[, 2L] <- 0.8 * X[, 1L] + sqrt(1 - 0.8^2) * X[, 2L]
  eta <- -1 + X %*% c(1, 1, -0.7, -0.7, 0, 0, 0, 0)
  y <- stats::rbinom(n, 1, stats::plogis(eta))
  if (length(unique(y)) != 2L) stop("Unexpected one-class V9 unit draw.")
  preprocess <- prepare_lsg_design(X, groups)
  target_original <- c(0.9, 0.9, -0.5, -0.5, rep(0, 4L))
  fit <- lsg_fit_joint_path_v9(
    X, y, groups, lambda = c(0.8, 0.3, 0.1), d = c(0, 0.5),
    alpha = 0.4, target_original = target_original,
    preprocess = preprocess, max_sweeps = 1500L,
    kkt_tolerance = 2e-6, update_tolerance = 1e-11,
    intercept_tolerance = 1e-12, keep_traces = TRUE,
    compile = FALSE, root = root
  )
  traces <- Filter(length, fit$objective_traces)
  monotone <- all(vapply(traces, function(value) {
    all(diff(value) <= 3e-12 * (1 + abs(head(value, -1L))))
  }, logical(1)))

  target <- project_lsg_original_target(preprocess, target_original)
  recomputed_objective <- matrix(NA_real_, nrow(fit$objective),
                                 ncol(fit$objective))
  recomputed_kkt <- matrix(NA_real_, nrow(fit$kkt), ncol(fit$kkt))
  for (di in seq_along(fit$d)) {
    for (li in seq_along(fit$lambda)) {
      beta <- fit$beta_solver[, li, di]
      intercept <- fit$intercept_solver[li, di]
      recomputed_objective[li, di] <- lsg_objective_cpp(
        preprocess$X, y, beta, intercept, preprocess$group_start,
        preprocess$group_end, preprocess$group_weight, target,
        fit$lambda[li], fit$alpha, fit$d[di]
      )
      recomputed_kkt[li, di] <- lsg_kkt_cpp(
        preprocess$X, y, beta, intercept, preprocess$group_start,
        preprocess$group_end, preprocess$group_weight, target,
        fit$lambda[li], fit$alpha, fit$d[di]
      )$maximum
    }
  }
  start_a <- lsg_fit_one_joint_v9_cpp(
    preprocess$X, y, preprocess$group_start, preprocess$group_end,
    preprocess$group_weight, target, 0.3, 0.4, 0.5,
    rep(0, ncol(preprocess$X)), stats::qlogis(mean(y)),
    2000L, 1e-8, 1e-12, 1e-12, 100L, TRUE, TRUE
  )
  set.seed(9010902L)
  start_b <- lsg_fit_one_joint_v9_cpp(
    preprocess$X, y, preprocess$group_start, preprocess$group_end,
    preprocess$group_weight, target, 0.3, 0.4, 0.5,
    stats::rnorm(ncol(preprocess$X), sd = 0.2), 0,
    2000L, 1e-8, 1e-12, 1e-12, 100L, TRUE, TRUE
  )
  alpha_one_0 <- lsg_fit_one_joint_v9_cpp(
    preprocess$X, y, preprocess$group_start, preprocess$group_end,
    preprocess$group_weight, target, 0.2, 1, 0,
    rep(0, ncol(preprocess$X)), stats::qlogis(mean(y)),
    2000L, 1e-8, 1e-12, 1e-12, 100L, TRUE, FALSE
  )
  alpha_one_1 <- lsg_fit_one_joint_v9_cpp(
    preprocess$X, y, preprocess$group_start, preprocess$group_end,
    preprocess$group_weight, target, 0.2, 1, 1,
    rep(0, ncol(preprocess$X)), stats::qlogis(mean(y)),
    2000L, 1e-8, 1e-12, 1e-12, 100L, TRUE, FALSE
  )

  status <- lsg_classify_grpreg_path_v9(
    "grLasso", c(12, 18, 20, 950), 5L, 1000L
  )
  prefix <- lsg_select_grpreg_prefix_v9(
    c(0.8, 0.4, 0.6, 0.7), rep(TRUE, 4L), rep(TRUE, 4L), status
  )
  boundary <- lsg_select_grpreg_prefix_v9(
    c(0.8, 0.6, 0.4, 0.7), rep(TRUE, 4L), rep(TRUE, 4L), status
  )

  checks <- c(
    joint_path_all_converged = all(fit$converged),
    joint_path_full_kkt = max(fit$kkt) <= 2e-6,
    joint_path_intercept_kkt = max(fit$intercept_kkt) <= 1e-10,
    joint_path_objective_monotone = monotone,
    joint_path_accepted_increase_is_roundoff =
      max(fit$maximum_accepted_objective_increase) <= 5e-12,
    joint_path_outputs_finite = all(is.finite(c(
      fit$beta_solver, fit$intercept_solver, fit$objective, fit$kkt
    ))),
    frozen_objective_reconstruction =
      max(abs(fit$objective - recomputed_objective)) <= 2e-12,
    frozen_kkt_reconstruction =
      max(abs(fit$kkt - recomputed_kkt)) <= 2e-12,
    independent_starts_converge = isTRUE(start_a$converged) &&
      isTRUE(start_b$converged),
    independent_starts_objective_agree =
      abs(start_a$objective - start_b$objective) <= 1e-8,
    independent_starts_prediction_agree = max(abs(
      stats::plogis(start_a$intercept + preprocess$X %*% start_a$beta) -
        stats::plogis(start_b$intercept + preprocess$X %*% start_b$beta)
    )) <= 1e-5,
    alpha_one_d_invariant = abs(
      alpha_one_0$objective - alpha_one_1$objective
    ) <= 1e-12 && max(abs(alpha_one_0$beta - alpha_one_1$beta)) <= 1e-10,
    grlasso_budget_prefix_recognized =
      isTRUE(status$group_lasso_safe_prefix) &&
      isTRUE(status$path_termination_acceptable) &&
      identical(status$usable_prefix_length, 3L) &&
      identical(status$returned_lower_boundary_excluded,
                c(FALSE, FALSE, FALSE, TRUE)),
    grlasso_interior_prefix_selected = isTRUE(prefix$accepted) &&
      identical(prefix$selected_index, 2L),
    grlasso_prefix_boundary_rejected = !isTRUE(boundary$accepted),
    portable_numeric_tolerance_accepts_last_bit =
      lsg_portable_numeric_equal_v9(1, 1 + 5e-17, 1e-14),
    portable_numeric_tolerance_rejects_material_change =
      !lsg_portable_numeric_equal_v9(1, 1 + 1e-6, 1e-14)
  )
  data.frame(check = names(checks), passed = unname(checks),
             stringsAsFactors = FALSE)
}
