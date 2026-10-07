# Hybrid numerical solver for the frozen Logistic SGLASSO objective.
#
# This research-only module does not alter the public sglasso package. It
# retains the V9 profiled block updates as a fast stage and adds a mathematically
# exact monotone proximal-gradient fallback for difficult finite-lambda points.

lsg_hybrid_defaults_v11 <- function() {
  list(
    schema_version = "logistic_sglasso_hybrid_solver_v11",
    path = list(
      block_max_sweeps = 4000L,
      block_chunk_sweeps = 100L,
      block_stall_window = 200L,
      block_stall_relative_improvement = 0.01,
      apg_max_iterations = 10000L,
      apg_kkt_check_interval = 5L,
      # V11.2: select candidates only after high-precision path fitting.
      # Scientific eligibility limits remain frozen below.
      kkt_tolerance = 1e-10,
      study_kkt_limit = 2.05e-6,
      update_tolerance = 1e-10,
      intercept_tolerance = 1e-12,
      max_intercept_iterations = 100L
    ),
    polish = list(
      block_max_sweeps = 4000L,
      block_chunk_sweeps = 100L,
      block_stall_window = 200L,
      block_stall_relative_improvement = 0.005,
      apg_max_iterations = 30000L,
      apg_kkt_check_interval = 2L,
      kkt_tolerance = 1e-11,
      study_kkt_limit = 1.05e-7,
      update_tolerance = 1e-12,
      intercept_tolerance = 1e-13,
      max_intercept_iterations = 200L
    )
  )
}


lsg_compile_hybrid_solver_v11 <- function(root = lsg_prework_root(),
                                           rebuild = FALSE,
                                           quiet = TRUE) {
  required <- c(
    "lsg_shifted_group_prox_v11_cpp",
    "lsg_fit_one_hybrid_v11_cpp",
    "lsg_path_hybrid_v11_cpp"
  )
  envir <- environment(lsg_compile_hybrid_solver_v11)
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
  source_file <- file.path(root, "src", "logistic_sglasso_hybrid_solver_v11.cpp")
  header_file <- file.path(root, "src", "logistic_sglasso_block_kernel_v11.hpp")
  if (!all(file.exists(c(source_file, header_file))) ||
      any(dir.exists(c(source_file, header_file))) ||
      any(nzchar(Sys.readlink(c(source_file, header_file))))) {
    stop("Missing regular V11 hybrid-solver source or header.", call. = FALSE)
  }

  previous_cppflags <- Sys.getenv("PKG_CPPFLAGS", unset = NA_character_)
  include_flag <- paste0("-I", shQuote(file.path(root, "src")))
  cppflags <- if (is.na(previous_cppflags) || !nzchar(previous_cppflags)) {
    include_flag
  } else paste(previous_cppflags, include_flag)
  Sys.setenv(PKG_CPPFLAGS = cppflags)

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
    if (is.na(previous_cppflags)) Sys.unsetenv("PKG_CPPFLAGS")
    else Sys.setenv(PKG_CPPFLAGS = previous_cppflags)
    if (is.na(previous_makevars)) Sys.unsetenv("R_MAKEVARS_USER")
    else Sys.setenv(R_MAKEVARS_USER = previous_makevars)
  }, add = TRUE)

  Rcpp::sourceCpp(
    source_file, rebuild = isTRUE(rebuild), showOutput = !isTRUE(quiet),
    verbose = FALSE, env = envir
  )
  invisible(TRUE)
}


lsg_validate_hybrid_controls_v11 <- function(controls) {
  required <- c(
    "block_max_sweeps", "block_chunk_sweeps", "block_stall_window",
    "block_stall_relative_improvement", "apg_max_iterations",
    "apg_kkt_check_interval", "kkt_tolerance", "update_tolerance",
    "intercept_tolerance", "max_intercept_iterations"
  )
  allowed <- c(required, "study_kkt_limit")
  if (!is.list(controls) || length(setdiff(required, names(controls))) ||
      length(setdiff(names(controls), allowed))) {
    stop("V11 numerical controls have an invalid schema.", call. = FALSE)
  }
  integer_names <- c(
    "block_max_sweeps", "block_chunk_sweeps", "block_stall_window",
    "apg_max_iterations", "apg_kkt_check_interval",
    "max_intercept_iterations"
  )
  for (name in integer_names) {
    value <- controls[[name]]
    if (length(value) != 1L || is.na(value) || !is.finite(value) ||
        value < 1 || value != as.integer(value)) {
      stop(name, " must be one positive integer.", call. = FALSE)
    }
    controls[[name]] <- as.integer(value)
  }
  positive_names <- c(
    "block_stall_relative_improvement", "kkt_tolerance",
    "update_tolerance", "intercept_tolerance"
  )
  if ("study_kkt_limit" %in% names(controls)) {
    positive_names <- c(positive_names, "study_kkt_limit")
  }
  for (name in positive_names) {
    value <- controls[[name]]
    if (length(value) != 1L || is.na(value) || !is.finite(value) || value <= 0) {
      stop(name, " must be one finite positive number.", call. = FALSE)
    }
    controls[[name]] <- as.numeric(value)
  }
  if (controls$block_stall_relative_improvement >= 1) {
    stop("block_stall_relative_improvement must be smaller than one.",
         call. = FALSE)
  }
  controls
}


lsg_fit_hybrid_path_v11 <- function(
    X,
    y,
    group,
    lambda,
    d,
    alpha,
    target_original,
    preprocess = NULL,
    controls = lsg_hybrid_defaults_v11()$path,
    use_active_set = TRUE,
    enable_fallback = TRUE,
    warm_start_d = TRUE,
    keep_traces = FALSE,
    lambda_order = c("decreasing", "any"),
    compile = TRUE,
    root = lsg_prework_root()
) {
  lambda_order <- match.arg(lambda_order)
  controls <- lsg_validate_hybrid_controls_v11(controls)
  if (isTRUE(compile)) lsg_compile_hybrid_solver_v11(root)
  response <- normalize_binary_response(y)$y
  X <- as.matrix(X)
  if (nrow(X) != length(response)) {
    stop("X and y have incompatible dimensions.", call. = FALSE)
  }
  lambda <- as.numeric(lambda)
  d <- as.numeric(d)
  alpha <- as.numeric(alpha)
  if (!length(lambda) || any(!is.finite(lambda)) || any(lambda <= 0) ||
      !length(d) || any(!is.finite(d)) || any(d < 0 | d > 1) ||
      length(alpha) != 1L || !is.finite(alpha) || alpha < 0 || alpha > 1 ||
      (lambda_order == "decreasing" && is.unsorted(-lambda))) {
    stop("Invalid V11 alpha, d, or finite lambda path.", call. = FALSE)
  }
  logical_controls <- list(
    use_active_set = use_active_set, enable_fallback = enable_fallback,
    warm_start_d = warm_start_d, keep_traces = keep_traces
  )
  valid_logical <- vapply(logical_controls, function(value) {
    is.logical(value) && length(value) == 1L && !is.na(value)
  }, logical(1))
  if (!all(valid_logical)) {
    stop("V11 logical controls must be scalar and nonmissing.", call. = FALSE)
  }
  target_original <- as.numeric(target_original)
  if (length(target_original) != ncol(X) || any(!is.finite(target_original))) {
    stop("target_original must be finite and match X.", call. = FALSE)
  }
  if (is.null(preprocess)) {
    preprocess <- prepare_lsg_design(X, group)
  } else if (!inherits(preprocess, "lsg_preprocess") ||
             preprocess$n != nrow(X) || preprocess$p_original != ncol(X)) {
    stop("preprocess is incompatible with X.", call. = FALSE)
  }
  target <- project_lsg_original_target(preprocess, target_original)
  core <- lsg_path_hybrid_v11_cpp(
    preprocess$X, response, preprocess$group_start, preprocess$group_end,
    preprocess$group_weight, target, lambda, d, alpha,
    controls$block_max_sweeps, controls$block_chunk_sweeps,
    controls$block_stall_window,
    controls$block_stall_relative_improvement,
    controls$apg_max_iterations, controls$apg_kkt_check_interval,
    controls$kkt_tolerance, controls$update_tolerance,
    controls$intercept_tolerance, controls$max_intercept_iterations,
    isTRUE(use_active_set), isTRUE(enable_fallback), isTRUE(warm_start_d),
    isTRUE(keep_traces)
  )

  coefficient_array <- array(
    0, dim = c(preprocess$p_original + 1L, length(lambda), length(d)),
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
    selected_groups = selected_groups,
    termination_reason = core$termination_reason,
    block_termination = core$block_termination,
    solver_route = core$solver_route,
    fallback_used = core$fallback_used != 0,
    block_stalled = core$block_stalled != 0,
    block_initial_kkt = core$block_initial_kkt,
    fallback_initial_kkt = core$fallback_initial_kkt,
    block_sweeps = core$block_sweeps,
    block_chunks = core$block_chunks,
    apg_iterations = core$apg_iterations,
    intercept_iterations = core$intercept_iterations,
    block_backtracking_steps = core$block_backtracking_steps,
    apg_backtracking_steps = core$apg_backtracking_steps,
    apg_restarts = core$apg_restarts,
    kkt_scans = core$kkt_scans,
    group_updates = core$group_updates,
    active_groups = core$active_groups,
    group_curvature = core$group_curvature_upper_bound,
    global_logistic_curvature = core$global_logistic_curvature,
    maximum_raw_objective_increase = core$max_raw_objective_increase,
    maximum_accepted_objective_increase =
      core$max_accepted_objective_increase,
    objective_traces = core$objective_traces,
    controls = controls,
    use_active_set = isTRUE(use_active_set),
    fallback_enabled = isTRUE(enable_fallback),
    warm_start_d = isTRUE(warm_start_d),
    call = match.call()
  )
  class(out) <- c(
    "logistic_sglasso_hybrid_prework_v11", "logistic_sglasso_prework"
  )
  out
}


lsg_shifted_group_prox_reference_v11 <- function(
    point, gradient, group_start, group_end, group_weight, target,
    lambda, alpha, d, curvature
) {
  answer <- numeric(length(point))
  for (g in seq_along(group_start)) {
    index <- seq.int(group_start[g] + 1L, group_end[g] + 1L)
    lambda1 <- lambda * alpha * group_weight[g]
    lambda2 <- lambda * (1 - alpha) * group_weight[g]
    score <- curvature * point[index] - gradient[index] +
      lambda2 * d * target[index]
    score_norm <- sqrt(sum(score^2))
    if (score_norm > lambda1 && score_norm > 0) {
      answer[index] <- (1 - lambda1 / score_norm) * score /
        (curvature + lambda2)
    }
  }
  answer
}


lsg_hybrid_solver_unit_checks_v11 <- function(root = lsg_prework_root()) {
  compile_lsg_core(rebuild = FALSE, quiet = TRUE)
  lsg_compile_hybrid_solver_v11(root, rebuild = FALSE, quiet = TRUE)
  defaults <- lsg_hybrid_defaults_v11()

  point <- c(0.4, -0.7, 0.2, 0.9, -0.1)
  gradient <- c(-0.3, 0.2, 0.5, -0.4, 0.1)
  target <- c(0.8, 0.6, -0.2, -0.3, 0.5)
  group_start <- as.integer(c(0, 2))
  group_end <- as.integer(c(1, 4))
  group_weight <- c(sqrt(2), 2.3)
  prox_errors <- unlist(lapply(c(0, 0.4, 1), function(alpha) {
    vapply(c(0, 0.5, 1), function(d) {
      observed <- lsg_shifted_group_prox_v11_cpp(
        point, gradient, group_start, group_end, group_weight, target,
        0.37, alpha, d, 1.8
      )
      reference <- lsg_shifted_group_prox_reference_v11(
        point, gradient, group_start, group_end, group_weight, target,
        0.37, alpha, d, 1.8
      )
      max(abs(observed - reference))
    }, numeric(1))
  }))

  set.seed(11011001L)
  n <- 90L
  group <- rep(seq_len(6L), each = 2L)
  X <- matrix(stats::rnorm(n * length(group)), nrow = n)
  X[, 2L] <- 0.85 * X[, 1L] + sqrt(1 - 0.85^2) * X[, 2L]
  X[, 4L] <- 0.85 * X[, 3L] + sqrt(1 - 0.85^2) * X[, 4L]
  eta <- -0.8 + X %*% c(1, 1, -0.8, -0.8, rep(0, 8L))
  y <- stats::rbinom(n, 1, stats::plogis(eta))
  if (length(unique(y)) != 2L) stop("Unexpected one-class V11 unit draw.")
  target_original <- c(0.8, 0.8, -0.6, -0.6, rep(0, 8L))
  controls <- defaults$path
  controls$block_max_sweeps <- 600L
  controls$block_chunk_sweeps <- 50L
  controls$block_stall_window <- 100L
  controls$apg_max_iterations <- 5000L
  controls$kkt_tolerance <- 1e-7
  controls$update_tolerance <- 1e-11
  fit <- lsg_fit_hybrid_path_v11(
    X, y, group, lambda = c(0.6, 0.25, 0.1), d = c(0, 0.5),
    alpha = 0.4, target_original = target_original,
    controls = controls, keep_traces = TRUE, compile = FALSE, root = root
  )
  traces <- Filter(length, fit$objective_traces)
  monotone <- all(vapply(traces, function(value) {
    all(diff(value) <= 2e-10 * (1 + abs(head(value, -1L))))
  }, logical(1)))
  reconstructed_objective <- matrix(NA_real_, nrow(fit$objective),
                                    ncol(fit$objective))
  reconstructed_kkt <- matrix(NA_real_, nrow(fit$kkt), ncol(fit$kkt))
  for (di in seq_along(fit$d)) {
    for (li in seq_along(fit$lambda)) {
      beta <- fit$beta_solver[, li, di]
      intercept <- fit$intercept_solver[li, di]
      reconstructed_objective[li, di] <- lsg_objective_cpp(
        fit$preprocess$X, y, beta, intercept, fit$preprocess$group_start,
        fit$preprocess$group_end, fit$preprocess$group_weight, fit$target,
        fit$lambda[li], fit$alpha, fit$d[di]
      )
      reconstructed_kkt[li, di] <- lsg_kkt_cpp(
        fit$preprocess$X, y, beta, intercept, fit$preprocess$group_start,
        fit$preprocess$group_end, fit$preprocess$group_weight, fit$target,
        fit$lambda[li], fit$alpha, fit$d[di]
      )$maximum
    }
  }

  preprocess <- fit$preprocess
  target_solver <- fit$target
  forced <- lsg_fit_one_hybrid_v11_cpp(
    preprocess$X, y, preprocess$group_start, preprocess$group_end,
    preprocess$group_weight, target_solver, 0.15, 0.4, 0.5,
    rep(0, ncol(preprocess$X)), stats::qlogis(mean(y)),
    1L, 1L, 1L, 0.01, 5000L, 5L, 1e-7, 1e-11, 1e-12, 100L,
    TRUE, TRUE, TRUE
  )
  deliberately_unfinished <- lsg_fit_one_hybrid_v11_cpp(
    preprocess$X, y, preprocess$group_start, preprocess$group_end,
    preprocess$group_weight, target_solver, 0.15, 0.4, 0.5,
    rep(0, ncol(preprocess$X)), stats::qlogis(mean(y)),
    1L, 1L, 1L, 0.01, 10L, 1L, 1e-12, 1e-14, 1e-13, 100L,
    TRUE, FALSE, FALSE
  )
  alpha_one_d0 <- lsg_fit_one_hybrid_v11_cpp(
    preprocess$X, y, preprocess$group_start, preprocess$group_end,
    preprocess$group_weight, target_solver, 0.2, 1, 0,
    rep(0, ncol(preprocess$X)), stats::qlogis(mean(y)),
    100L, 25L, 50L, 0.01, 5000L, 5L, 1e-8, 1e-12, 1e-12, 100L,
    TRUE, TRUE, FALSE
  )
  alpha_one_d1 <- lsg_fit_one_hybrid_v11_cpp(
    preprocess$X, y, preprocess$group_start, preprocess$group_end,
    preprocess$group_weight, target_solver, 0.2, 1, 1,
    rep(0, ncol(preprocess$X)), stats::qlogis(mean(y)),
    100L, 25L, 50L, 0.01, 5000L, 5L, 1e-8, 1e-12, 1e-12, 100L,
    TRUE, TRUE, FALSE
  )

  ridge <- lsg_fit_one_hybrid_v11_cpp(
    preprocess$X, y, preprocess$group_start, preprocess$group_end,
    preprocess$group_weight, target_solver, 0.2, 0, 0.5,
    rep(0, ncol(preprocess$X)), stats::qlogis(mean(y)),
    1L, 1L, 1L, 0.01, 5000L, 2L, 1e-9, 1e-12, 1e-13, 200L,
    TRUE, TRUE, FALSE
  )
  coordinate_weight <- rep(
    preprocess$group_weight,
    preprocess$group_end - preprocess$group_start + 1L
  )
  ridge_objective <- function(theta) {
    intercept <- theta[1L]
    beta <- theta[-1L]
    linear <- intercept + drop(preprocess$X %*% beta)
    mean(log1p(exp(-abs(linear))) + pmax(linear, 0) - y * linear) +
      0.5 * 0.2 * sum(coordinate_weight *
                        (beta - 0.5 * target_solver)^2)
  }
  ridge_gradient <- function(theta) {
    intercept <- theta[1L]
    beta <- theta[-1L]
    residual <- stats::plogis(
      intercept + drop(preprocess$X %*% beta)
    ) - y
    c(mean(residual), drop(crossprod(preprocess$X, residual)) / n +
        0.2 * coordinate_weight * (beta - 0.5 * target_solver))
  }
  reference <- stats::optim(
    c(stats::qlogis(mean(y)), rep(0, ncol(preprocess$X))),
    ridge_objective, ridge_gradient, method = "BFGS",
    control = list(maxit = 10000L, reltol = 1e-13)
  )
  ridge_prediction_error <- max(abs(
    stats::plogis(ridge$intercept + preprocess$X %*% ridge$beta) -
      stats::plogis(reference$par[1L] + preprocess$X %*%
                      reference$par[-1L])
  ))

  set.seed(11011002L)
  stress_n <- 50L
  stress_groups <- 30L
  latent <- matrix(stats::rnorm(stress_n * stress_groups), stress_n)
  stress_X <- do.call(cbind, lapply(seq_len(stress_groups), function(g) {
    sqrt(0.95) * latent[, g] +
      sqrt(0.05) * matrix(stats::rnorm(stress_n * 3L), stress_n, 3L)
  }))
  stress_eta <- -2 + stress_X[, 1L] - stress_X[, 4L]
  stress_y <- stats::rbinom(stress_n, 1, stats::plogis(stress_eta))
  if (length(unique(stress_y)) != 2L) {
    stop("Unexpected one-class V11 stress draw.", call. = FALSE)
  }
  stress_group <- rep(seq_len(stress_groups), each = 3L)
  stress_preprocess <- prepare_lsg_design(stress_X, stress_group)
  stress_target_original <- c(rep(0.4, 3L), rep(-0.4, 3L),
                              rep(0, ncol(stress_X) - 6L))
  stress_target <- project_lsg_original_target(
    stress_preprocess, stress_target_original
  )
  stress <- lsg_fit_one_hybrid_v11_cpp(
    stress_preprocess$X, stress_y, stress_preprocess$group_start,
    stress_preprocess$group_end, stress_preprocess$group_weight,
    stress_target, 0.12, 0.3, 0.7,
    rep(0, ncol(stress_preprocess$X)), stats::qlogis(mean(stress_y)),
    1L, 1L, 1L, 0.01, 8000L, 5L, 2e-6, 1e-10, 1e-12, 100L,
    TRUE, TRUE, FALSE
  )

  checks <- c(
    exact_shifted_group_prox = max(prox_errors) <= 2e-14,
    all_path_points_converged = all(fit$converged),
    all_path_kkt_within_declared_tolerance = max(fit$kkt) <= 1e-7,
    intercept_profiled = max(fit$intercept_kkt) <= 1e-10,
    accepted_objective_monotone = monotone,
    objective_reconstructed = max(abs(
      fit$objective - reconstructed_objective
    )) <= 3e-12,
    kkt_reconstructed = max(abs(fit$kkt - reconstructed_kkt)) <= 3e-12,
    forced_fallback_exercised = isTRUE(forced$fallback_used) &&
      identical(forced$solver_route,
                "profiled_block_then_monotone_apg") &&
      isTRUE(forced$converged) && forced$kkt <= 1e-7,
    unfinished_fit_is_not_false_convergence =
      !isTRUE(deliberately_unfinished$converged) &&
      identical(deliberately_unfinished$solver_route, "block_only"),
    alpha_one_is_d_invariant = isTRUE(alpha_one_d0$converged) &&
      isTRUE(alpha_one_d1$converged) &&
      abs(alpha_one_d0$objective - alpha_one_d1$objective) <= 1e-10 &&
      max(abs(alpha_one_d0$beta - alpha_one_d1$beta)) <= 1e-8,
    alpha_zero_matches_independent_smooth_reference =
      isTRUE(ridge$converged) && reference$convergence == 0L &&
      abs(ridge$objective - reference$value) <= 2e-8 &&
      ridge_prediction_error <= 2e-4,
    high_correlation_p_gt_n_fallback_converges =
      isTRUE(stress$fallback_used) && isTRUE(stress$converged) &&
      stress$kkt <= 2e-6,
    path_and_polish_tolerances_separate =
      defaults$polish$kkt_tolerance < defaults$path$kkt_tolerance &&
      defaults$polish$study_kkt_limit < defaults$path$study_kkt_limit,
    public_package_files_not_required = !any(grepl(
      "(^|/)R/(sglasso|cv\\.sglasso|helpers)\\.R$|(^|/)src/All_Functions\\.cpp$",
      c("logistic_prework/R/logistic_sglasso_hybrid_solver_v11.R",
        "logistic_prework/src/logistic_sglasso_hybrid_solver_v11.cpp")
    ))
  )
  data.frame(
    check = names(checks),
    passed = unname(vapply(checks, isTRUE, logical(1))),
    stringsAsFactors = FALSE
  )
}
