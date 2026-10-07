# Targeted numerical diagnostics for the failed V7 pilot.
#
# This module never changes the V7 estimator, tuning grid, selection rule, or
# source artifacts.  It reruns only seven frozen diagnostic cases using
# training and validation data.  Post-profiled intercepts keep beta fixed and
# are evidence for a later solver decision, not replacement fitted models.

lsg_v8_required_functions <- function(require_profiler = FALSE) {
  required <- c(
    "lsg_configuration_v7", "lsg_with_rng_v7", "lsg_data_fingerprint_v7",
    "simulate_logistic_sglasso_two_design_v1", "lsg_tail_validate_data_v6",
    "prepare_lsg_design", "project_lsg_original_target",
    "transform_lsg_newx", "recover_lsg_coefficients",
    "lsg_tail_reference_v6", "fit_logistic_sglasso",
    "predict_logistic_sglasso", "binary_log_loss", "lsg_kkt_cpp",
    "lsg_objective_cpp", "lsg_tail_limit_fit_v6",
    "lsg_select_sglasso_candidates_v7", "lsg_candidate_order_v7",
    "lsg_bind_tuning_rows_v7",
    "lsg_run_grpreg_task_v8", "select_valid_external_grid_v2"
  )
  if (isTRUE(require_profiler)) {
    required <- c(required, "lsg_profile_intercept_beta_v8_cpp",
                  "lsg_profile_intercept_offset_v8_cpp")
  }
  envir <- environment(lsg_v8_required_functions)
  missing <- required[!vapply(required, exists, logical(1), mode = "function",
                              envir = envir, inherits = TRUE)]
  if (length(missing)) {
    stop("Missing V8 diagnostic dependency: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  invisible(TRUE)
}


lsg_compile_intercept_profile_v8 <- function(root, rebuild = FALSE,
                                             quiet = TRUE) {
  compiled <- vapply(c(
    "lsg_profile_intercept_beta_v8_cpp",
    "lsg_profile_intercept_offset_v8_cpp"
  ), exists, logical(1), mode = "function",
  envir = environment(lsg_compile_intercept_profile_v8), inherits = TRUE)
  if (!isTRUE(rebuild) && all(compiled)) {
    return(invisible(TRUE))
  }
  if (!requireNamespace("Rcpp", quietly = TRUE) ||
      !requireNamespace("RcppArmadillo", quietly = TRUE)) {
    stop("Rcpp and RcppArmadillo are required; no installation is attempted.",
         call. = FALSE)
  }
  path <- file.path(root, "src", "logistic_sglasso_intercept_profile_v8.cpp")
  if (!file.exists(path) || nzchar(Sys.readlink(path))) {
    stop("Missing regular V8 intercept-profiler source file.", call. = FALSE)
  }
  # Reuse the project-local macOS toolchain shim without changing any user
  # Makevars file.  On TRUBA these files are absent and the loaded R module's
  # native compiler configuration is used unchanged.
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
    if (is.na(previous_makevars)) {
      Sys.unsetenv("R_MAKEVARS_USER")
    } else {
      Sys.setenv(R_MAKEVARS_USER = previous_makevars)
    }
  }, add = TRUE)
  Rcpp::sourceCpp(
    path, rebuild = isTRUE(rebuild), showOutput = !isTRUE(quiet),
    verbose = FALSE, env = environment(lsg_compile_intercept_profile_v8)
  )
  invisible(TRUE)
}


lsg_read_failure_cases_v8 <- function(path) {
  cases <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                           na.strings = c("", "NA"))
  required <- c(
    "case_id", "diagnostic_type", "task_id", "scenario_index", "scenario",
    "replication", "seed", "alpha", "d", "representative_rank",
    "archived_anchor_type", "archived_anchor_lambda_index",
    "archived_anchor_lambda_relative", "archived_lambda_reference",
    "archived_anchor_kkt", "archived_anchor_intercept_score",
    "archived_anchor_validation_log_loss", "source_shard_file",
    "expected_v7_status", "expected_v7_message", "description"
  )
  if (!identical(names(cases), required)) {
    stop("The V8 failure-case file has an unexpected schema.", call. = FALSE)
  }
  integer_columns <- c(
    "task_id", "scenario_index", "replication", "seed",
    "representative_rank", "archived_anchor_lambda_index"
  )
  numeric_columns <- c(
    "alpha", "d", "archived_anchor_lambda_relative",
    "archived_lambda_reference", "archived_anchor_kkt",
    "archived_anchor_intercept_score",
    "archived_anchor_validation_log_loss"
  )
  for (name in integer_columns) cases[[name]] <- as.integer(cases[[name]])
  for (name in numeric_columns) cases[[name]] <- as.numeric(cases[[name]])
  expected_types <- c(
    rep("sglasso_intercept", 3L), "sglasso_lambda_tail",
    rep("grpreg_group_lasso_raw", 3L)
  )
  if (nrow(cases) != 7L || anyDuplicated(cases$case_id) ||
      anyDuplicated(cases$task_id) || !identical(cases$diagnostic_type, expected_types) ||
      !identical(cases$task_id, c(153L, 143L, 158L, 125L, 126L, 134L, 150L)) ||
      anyNA(cases[c("case_id", "diagnostic_type", "task_id", "scenario",
                    "replication", "seed", "description")]) ||
      any(!nzchar(cases$case_id)) || any(!nzchar(cases$scenario))) {
    stop("The seven frozen V8 failure cases changed.", call. = FALSE)
  }
  numerical_cases <- cases$diagnostic_type != "grpreg_group_lasso_raw"
  if (any(!is.finite(as.matrix(cases[numerical_cases, numeric_columns]))) ||
      any(cases$alpha[numerical_cases] < 0 | cases$alpha[numerical_cases] > 1) ||
      any(cases$d[numerical_cases] < 0 | cases$d[numerical_cases] > 1)) {
    stop("Invalid frozen SGLASSO diagnostic case values.", call. = FALSE)
  }
  rownames(cases) <- NULL
  cases
}


lsg_v8_tuning_data <- function(scenario, seed) {
  generated <- lsg_with_rng_v7(
    as.integer(seed),
    simulate_logistic_sglasso_two_design_v1(scenario, as.integer(seed))
  )
  # The generator must draw the full frozen sample to preserve RNG identity;
  # test and truth fields are discarded before any diagnostic interface.
  data <- generated[c(
    "X_train", "y_train", "X_validation", "y_validation", "group"
  )]
  rm(generated)
  lsg_tail_validate_data_v6(data)
  data
}


lsg_v8_data_fingerprints <- function(data) {
  list(
    training = lsg_data_fingerprint_v7(list(
      X = data$X_train, y = data$y_train, group = data$group
    )),
    validation = lsg_data_fingerprint_v7(list(
      X = data$X_validation, y = data$y_validation
    ))
  )
}


lsg_v8_d_prefix <- function(configuration, alpha, target_d) {
  if (abs(alpha - 1) <= 1e-12) return(0)
  grid <- sort(unique(as.numeric(configuration$d_grid)))
  prefix <- grid[grid <= target_d + 1e-12]
  if (!length(prefix) || abs(tail(prefix, 1L) - target_d) > 1e-12) {
    stop("The target d is absent from the frozen V7 grid.", call. = FALSE)
  }
  prefix
}


lsg_v8_fit_path <- function(data, target_original, alpha, target_d,
                            lambda_ratios, configuration) {
  lsg_v8_required_functions(require_profiler = TRUE)
  lambda_ratios <- as.numeric(lambda_ratios)
  if (!length(lambda_ratios) || any(!is.finite(lambda_ratios)) ||
      any(lambda_ratios <= 0) || is.unsorted(-lambda_ratios, strictly = TRUE)) {
    stop("Diagnostic lambda ratios must be finite, positive and decreasing.",
         call. = FALSE)
  }
  preprocess <- prepare_lsg_design(data$X_train, data$group)
  reference <- lsg_tail_reference_v6(
    preprocess, data$y_train, target_original, alpha
  )
  d_prefix <- lsg_v8_d_prefix(configuration, alpha, target_d)
  fit <- fit_logistic_sglasso(
    data$X_train, data$y_train, data$group,
    lambda = reference * lambda_ratios,
    d = d_prefix, alpha = alpha,
    max_outer = configuration$max_passes,
    max_inner = configuration$max_inner,
    tolerance = configuration$tolerance,
    inner_tolerance = configuration$inner_tolerance,
    target_original = target_original,
    compile = FALSE, preprocess = preprocess,
    use_active_set = TRUE, warm_start_d = TRUE,
    solver = configuration$solver
  )
  d_index <- which(abs(fit$d - target_d) <= 1e-12)
  if (length(d_index) != 1L || length(fit$lambda) != length(lambda_ratios)) {
    stop("The diagnostic SGLASSO path is misaligned.", call. = FALSE)
  }
  list(fit = fit, preprocess = preprocess, lambda_reference = reference,
       lambda_ratios = lambda_ratios, d_index = d_index)
}


lsg_v8_profile_path <- function(data, fitted, target_original, alpha, target_d) {
  fit <- fitted$fit
  preprocess <- fitted$preprocess
  di <- fitted$d_index
  target_solver <- project_lsg_original_target(preprocess, target_original)
  X_validation_solver <- transform_lsg_newx(preprocess, data$X_validation)
  rows <- vector("list", length(fitted$lambda_ratios))
  for (li in seq_along(fitted$lambda_ratios)) {
    beta <- as.numeric(fit$beta_solver[, li, di])
    intercept <- as.numeric(fit$intercept_solver[li, di])
    baseline_kkt <- lsg_kkt_cpp(
      preprocess$X, data$y_train, beta, intercept,
      preprocess$group_start, preprocess$group_end, preprocess$group_weight,
      target_solver, fit$lambda[li], alpha, target_d
    )
    baseline_objective <- lsg_objective_cpp(
      preprocess$X, data$y_train, beta, intercept,
      preprocess$group_start, preprocess$group_end, preprocess$group_weight,
      target_solver, fit$lambda[li], alpha, target_d
    )
    baseline_training_probability <- stats::plogis(
      intercept + drop(preprocess$X %*% beta)
    )
    baseline_validation_probability <- stats::plogis(
      intercept + drop(X_validation_solver %*% beta)
    )
    baseline_original <- recover_lsg_coefficients(preprocess, beta, intercept)
    baseline_reconstruction_error <- max(abs(
      baseline_validation_probability - stats::plogis(
        baseline_original[1L] +
          drop(data$X_validation %*% baseline_original[-1L])
      )
    ))
    profiled <- lsg_profile_intercept_beta_v8_cpp(
      preprocess$X, data$y_train, beta, initial_intercept = intercept,
      score_tolerance = 1e-12, bracket_tolerance = 1e-12,
      max_iterations = 100L
    )
    profiled_kkt <- lsg_kkt_cpp(
      preprocess$X, data$y_train, beta, profiled$intercept,
      preprocess$group_start, preprocess$group_end, preprocess$group_weight,
      target_solver, fit$lambda[li], alpha, target_d
    )
    profiled_objective <- lsg_objective_cpp(
      preprocess$X, data$y_train, beta, profiled$intercept,
      preprocess$group_start, preprocess$group_end, preprocess$group_weight,
      target_solver, fit$lambda[li], alpha, target_d
    )
    profiled_validation_probability <- stats::plogis(
      profiled$intercept + drop(X_validation_solver %*% beta)
    )
    original <- recover_lsg_coefficients(preprocess, beta, profiled$intercept)
    reconstruction_error <- max(abs(
      profiled_validation_probability - stats::plogis(
        original[1L] + drop(data$X_validation %*% original[-1L])
      )
    ))
    rows[[li]] <- data.frame(
      lambda_index = li,
      lambda_ratio = fitted$lambda_ratios[li],
      lambda = fit$lambda[li],
      alpha = alpha,
      d = target_d,
      baseline_intercept = intercept,
      baseline_training_log_loss = binary_log_loss(
        data$y_train, baseline_training_probability
      ),
      baseline_validation_log_loss = binary_log_loss(
        data$y_validation, baseline_validation_probability
      ),
      baseline_objective = baseline_objective,
      baseline_kkt = baseline_kkt$maximum,
      baseline_intercept_kkt = baseline_kkt$intercept,
      baseline_group_kkt = max(baseline_kkt$by_group),
      baseline_kkt_decomposition_error = abs(
        baseline_kkt$maximum -
          max(baseline_kkt$intercept, max(baseline_kkt$by_group))
      ),
      baseline_converged = isTRUE(fit$converged[li, di]),
      baseline_passes = as.integer(fit$passes[li, di]),
      baseline_selected_groups = as.integer(fit$selected_groups[li, di]),
      baseline_maximum_raw_objective_increase =
        fit$maximum_raw_objective_increase[li, di],
      baseline_beta_finite = all(is.finite(beta)),
      baseline_probability_finite = all(is.finite(c(
        baseline_training_probability, baseline_validation_probability
      ))),
      baseline_reconstruction_error = baseline_reconstruction_error,
      profiled_intercept = profiled$intercept,
      profiled_intercept_shift = profiled$intercept_shift,
      profiled_training_log_loss = profiled$training_log_loss,
      profiled_validation_log_loss = binary_log_loss(
        data$y_validation, profiled_validation_probability
      ),
      profiled_objective = profiled_objective,
      profiled_objective_minus_baseline = profiled_objective - baseline_objective,
      profiled_kkt = profiled_kkt$maximum,
      profiled_intercept_kkt = profiled_kkt$intercept,
      profiled_group_kkt = max(profiled_kkt$by_group),
      profiled_kkt_decomposition_error = abs(
        profiled_kkt$maximum -
          max(profiled_kkt$intercept, max(profiled_kkt$by_group))
      ),
      profiled_internal_kkt_pass = profiled_kkt$maximum <= 2e-6,
      profiled_study_kkt_pass = profiled_kkt$maximum <= 2.05e-6,
      profile_root_converged = isTRUE(profiled$converged),
      profile_iterations = as.integer(profiled$iterations),
      profile_convergence_reason = profiled$convergence_reason,
      profile_beta_unchanged = TRUE,
      profiled_probability_finite = all(is.finite(c(
        profiled$probability, profiled_validation_probability
      ))),
      profiled_reconstruction_error = reconstruction_error,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}


lsg_run_intercept_case_v8 <- function(case, scenario, source, configuration) {
  if (!is.data.frame(case) || nrow(case) != 1L ||
      case$diagnostic_type != "sglasso_intercept") {
    stop("One frozen intercept-profile case is required.", call. = FALSE)
  }
  data <- lsg_v8_tuning_data(scenario, case$seed)
  fingerprints <- lsg_v8_data_fingerprints(data)
  if (!identical(fingerprints, source$data_sha256)) {
    stop("Regenerated training/validation data differ from the V7 shard.",
         call. = FALSE)
  }
  target <- as.numeric(source$firth_target_original)
  if (length(target) != ncol(data$X_train) || any(!is.finite(target))) {
    stop("The blinded V7 Firth target is unavailable.", call. = FALSE)
  }
  fitted <- lsg_v8_fit_path(
    data, target, case$alpha, case$d,
    configuration$lambda_relative_grid, configuration
  )
  comparison <- lsg_v8_profile_path(
    data, fitted, target, case$alpha, case$d
  )
  archived <- source$sglasso_tuning
  archived <- archived[
    archived$point_type == "finite" &
      abs(archived$alpha - case$alpha) <= 1e-12 &
      abs(archived$d - case$d) <= 1e-12,
    , drop = FALSE
  ]
  archived <- archived[order(archived$lambda_index), , drop = FALSE]
  if (nrow(archived) != nrow(comparison)) {
    stop("The archived V7 SGLASSO path is incomplete.", call. = FALSE)
  }
  comparison$source_lambda_error <- abs(comparison$lambda - archived$lambda)
  comparison$source_validation_log_loss <- archived$validation_log_loss
  comparison$source_validation_loss_error <- abs(
    comparison$baseline_validation_log_loss - archived$validation_log_loss
  )
  comparison$source_kkt <- archived$kkt
  comparison$source_kkt_error <- abs(comparison$baseline_kkt - archived$kkt)
  comparison$source_intercept_kkt <- archived$intercept_score
  comparison$source_intercept_kkt_error <- abs(
    comparison$baseline_intercept_kkt - archived$intercept_score
  )
  comparison$source_passes_match <-
    comparison$baseline_passes == as.integer(archived$passes)
  comparison$source_converged_match <-
    comparison$baseline_converged == as.logical(archived$converged)
  comparison$case_id <- case$case_id
  comparison$task_id <- case$task_id
  comparison$scenario <- case$scenario
  comparison$replication <- case$replication
  comparison$seed <- case$seed
  comparison <- comparison[c(
    "case_id", "task_id", "scenario", "replication", "seed",
    setdiff(names(comparison), c(
      "case_id", "task_id", "scenario", "replication", "seed"
    ))
  )]
  summary <- data.frame(
    case_id = case$case_id,
    task_id = case$task_id,
    alpha = case$alpha,
    d = case$d,
    points = nrow(comparison),
    source_path_zero_eligible = !any(archived$numerically_eligible %in% TRUE),
    source_path_min_kkt = min(archived$kkt),
    source_path_min_kkt_lambda_ratio =
      archived$lambda_relative_to_reference[which.min(archived$kkt)],
    lambda_reference_recomputed = fitted$lambda_reference,
    lambda_reference_archived = unique(archived$lambda_reference),
    lambda_reference_error = abs(
      fitted$lambda_reference - unique(archived$lambda_reference)
    ),
    baseline_replay_max_lambda_error = max(comparison$source_lambda_error),
    baseline_replay_max_validation_loss_error =
      max(comparison$source_validation_loss_error),
    baseline_replay_max_kkt_error = max(comparison$source_kkt_error),
    baseline_replay_max_intercept_kkt_error =
      max(comparison$source_intercept_kkt_error),
    baseline_replay_passes_exact = all(comparison$source_passes_match),
    baseline_replay_convergence_exact = all(comparison$source_converged_match),
    profile_root_all_converged = all(comparison$profile_root_converged),
    profile_max_intercept_kkt = max(comparison$profiled_intercept_kkt),
    profile_max_objective_increase =
      max(comparison$profiled_objective_minus_baseline),
    profile_study_kkt_eligible_points =
      sum(comparison$profiled_study_kkt_pass),
    profile_path_has_study_kkt_candidate =
      any(comparison$profiled_study_kkt_pass),
    stringsAsFactors = FALSE
  )
  list(
    schema_version = "sglasso_intercept_failure_diagnostic_v8",
    case = case, comparison = comparison, summary = summary,
    data_sha256 = fingerprints, test_used = FALSE,
    raw = list(lambda_reference = fitted$lambda_reference,
               d_prefix = fitted$fit$d)
  )
}


lsg_run_tail_case_v8 <- function(case, scenario, source, configuration) {
  if (!is.data.frame(case) || nrow(case) != 1L ||
      case$diagnostic_type != "sglasso_lambda_tail") {
    stop("One frozen lambda-tail case is required.", call. = FALSE)
  }
  data <- lsg_v8_tuning_data(scenario, case$seed)
  fingerprints <- lsg_v8_data_fingerprints(data)
  if (!identical(fingerprints, source$data_sha256)) {
    stop("Regenerated Task 125 data differ from the V7 shard.", call. = FALSE)
  }
  target <- as.numeric(source$firth_target_original)
  ratios <- c(4096, 2048, 1024, configuration$lambda_relative_grid)
  fitted <- lsg_v8_fit_path(
    data, target, case$alpha, case$d, ratios, configuration
  )
  finite <- lsg_v8_profile_path(data, fitted, target, case$alpha, case$d)
  # A second, exact 36-point replay starts at the original 512 boundary.  It
  # distinguishes source reproducibility from the scientifically interesting
  # 39-point extension, whose three additional warm starts may legitimately
  # change iteration counts at otherwise identical lambda values.
  replay_fitted <- lsg_v8_fit_path(
    data, target, case$alpha, case$d,
    configuration$lambda_relative_grid, configuration
  )
  v7_replay <- lsg_v8_profile_path(
    data, replay_fitted, target, case$alpha, case$d
  )
  endpoint <- lsg_tail_limit_fit_v6(
    data, fitted$preprocess, target, case$alpha, case$d, configuration
  )
  endpoint_row <- data.frame(
    case_id = case$case_id, task_id = case$task_id,
    scenario = case$scenario, replication = case$replication, seed = case$seed,
    lambda_index = 0L, lambda_ratio = Inf, lambda = Inf,
    alpha = case$alpha, d = case$d,
    baseline_intercept = endpoint$intercept_solver,
    baseline_training_log_loss = endpoint$point$training_log_loss,
    baseline_validation_log_loss = endpoint$point$validation_log_loss,
    baseline_objective = NA_real_, baseline_kkt = NA_real_,
    baseline_intercept_kkt = endpoint$point$intercept_score,
    baseline_group_kkt = endpoint$point$penalty_kkt,
    baseline_kkt_decomposition_error = 0,
    baseline_converged = endpoint$point$numerically_eligible,
    baseline_passes = NA_integer_,
    baseline_selected_groups = endpoint$point$selected_groups,
    baseline_maximum_raw_objective_increase = NA_real_,
    baseline_beta_finite = all(is.finite(endpoint$beta_solver)),
    baseline_probability_finite = all(is.finite(c(
      endpoint$training_probability, endpoint$validation_probability
    ))),
    baseline_reconstruction_error = endpoint$point$reconstruction_error,
    profiled_intercept = endpoint$intercept_solver,
    profiled_intercept_shift = 0,
    profiled_training_log_loss = endpoint$point$training_log_loss,
    profiled_validation_log_loss = endpoint$point$validation_log_loss,
    profiled_objective = NA_real_, profiled_objective_minus_baseline = NA_real_,
    profiled_kkt = max(endpoint$point$intercept_score,
                       endpoint$point$penalty_kkt),
    profiled_intercept_kkt = endpoint$point$intercept_score,
    profiled_group_kkt = endpoint$point$penalty_kkt,
    profiled_kkt_decomposition_error = 0,
    profiled_internal_kkt_pass = endpoint$point$numerically_eligible,
    profiled_study_kkt_pass = endpoint$point$numerically_eligible,
    profile_root_converged = endpoint$point$numerically_eligible,
    profile_iterations = NA_integer_, profile_convergence_reason = "analytic_limit",
    profile_beta_unchanged = TRUE,
    profiled_probability_finite = all(is.finite(c(
      endpoint$training_probability, endpoint$validation_probability
    ))),
    profiled_reconstruction_error = endpoint$point$reconstruction_error,
    point_type = "penalty_limit", stringsAsFactors = FALSE
  )
  finite$case_id <- case$case_id
  finite$task_id <- case$task_id
  finite$scenario <- case$scenario
  finite$replication <- case$replication
  finite$seed <- case$seed
  finite$point_type <- "finite"
  finite <- finite[names(endpoint_row)]
  curve <- rbind(finite, endpoint_row)

  archived <- source$sglasso_tuning
  source_selected <- archived[archived$selected_free_d %in% TRUE, , drop = FALSE]
  if (nrow(source_selected) != 1L) {
    stop("Task 125 does not have one archived free-d selection.", call. = FALSE)
  }
  old_rows <- finite[finite$lambda_ratio <= 512, , drop = FALSE]
  old_source <- archived[
    archived$point_type == "finite" &
      abs(archived$alpha - case$alpha) <= 1e-12 &
      abs(archived$d - case$d) <= 1e-12,
    , drop = FALSE
  ]
  old_source <- old_source[order(old_source$lambda_index), , drop = FALSE]
  old_rows <- old_rows[order(-old_rows$lambda_ratio), , drop = FALSE]
  old_grid_exact <- nrow(old_rows) == nrow(old_source) &&
    max(abs(old_rows$lambda_ratio - old_source$lambda_relative_to_reference)) <= 1e-12
  anchor <- finite[finite$lambda_ratio == 512, , drop = FALSE]
  anchor_error <- if (nrow(anchor) == 1L) {
    abs(anchor$baseline_validation_log_loss -
          source_selected$validation_log_loss[[1L]])
  } else Inf
  source_endpoint <- archived[
    archived$point_type == "penalty_limit" &
      abs(archived$alpha - case$alpha) <= 1e-12 &
      abs(archived$d - case$d) <= 1e-12,
    , drop = FALSE
  ]
  endpoint_loss_error <- if (nrow(source_endpoint) == 1L) {
    abs(endpoint$point$validation_log_loss -
          source_endpoint$validation_log_loss[[1L]])
  } else Inf
  source_reference <- unique(old_source$lambda_reference)
  reference_error <- if (length(source_reference) == 1L) {
    abs(fitted$lambda_reference - source_reference)
  } else Inf

  if (nrow(v7_replay) != nrow(old_source) ||
      max(abs(v7_replay$lambda_ratio -
                old_source$lambda_relative_to_reference)) > 1e-12) {
    stop("The direct Task 125 V7 replay is misaligned.", call. = FALSE)
  }
  v7_replay$source_lambda <- old_source$lambda
  v7_replay$source_lambda_error <- abs(
    v7_replay$lambda - old_source$lambda
  )
  v7_replay$source_validation_log_loss <- old_source$validation_log_loss
  v7_replay$source_validation_loss_error <- abs(
    v7_replay$baseline_validation_log_loss - old_source$validation_log_loss
  )
  v7_replay$source_kkt <- old_source$kkt
  v7_replay$source_kkt_error <- abs(
    v7_replay$baseline_kkt - old_source$kkt
  )
  v7_replay$source_intercept_kkt <- old_source$intercept_score
  v7_replay$source_intercept_kkt_error <- abs(
    v7_replay$baseline_intercept_kkt - old_source$intercept_score
  )
  v7_replay$source_passes <- as.integer(old_source$passes)
  v7_replay$source_passes_match <-
    v7_replay$baseline_passes == as.integer(old_source$passes)
  v7_replay$source_converged <- as.logical(old_source$converged)
  v7_replay$source_converged_match <-
    v7_replay$baseline_converged == as.logical(old_source$converged)
  v7_replay$case_id <- case$case_id
  v7_replay$task_id <- case$task_id
  v7_replay$scenario <- case$scenario
  v7_replay$replication <- case$replication
  v7_replay$seed <- case$seed
  v7_replay <- v7_replay[c(
    "case_id", "task_id", "scenario", "replication", "seed",
    setdiff(names(v7_replay), c(
      "case_id", "task_id", "scenario", "replication", "seed"
    ))
  )]

  # Bind the archived values only after the raw extended curve has been
  # constructed.  The three new ratios have no V7 counterpart; every one of
  # the unchanged 36 ratios and the analytic endpoint must have exactly one.
  curve$source_validation_log_loss <- NA_real_
  curve$source_kkt <- NA_real_
  curve$source_intercept_kkt <- NA_real_
  curve$source_passes <- NA_integer_
  curve$source_converged <- NA
  curve$source_validation_loss_error <- NA_real_
  curve$source_kkt_error <- NA_real_
  curve$source_intercept_kkt_error <- NA_real_
  curve$source_passes_match <- NA
  curve$source_converged_match <- NA
  old_curve_index <- vapply(
    old_source$lambda_relative_to_reference,
    function(ratio) which.min(abs(curve$lambda_ratio - ratio)),
    integer(1)
  )
  old_curve_error <- abs(
    curve$lambda_ratio[old_curve_index] -
      old_source$lambda_relative_to_reference
  )
  if (anyNA(old_curve_index) || anyDuplicated(old_curve_index) ||
      length(old_curve_index) != 36L || max(old_curve_error) > 1e-12) {
    stop("The extended Task 125 curve does not contain the frozen V7 grid.",
         call. = FALSE)
  }
  curve$source_validation_log_loss[old_curve_index] <-
    old_source$validation_log_loss
  curve$source_kkt[old_curve_index] <- old_source$kkt
  curve$source_intercept_kkt[old_curve_index] <- old_source$intercept_score
  curve$source_passes[old_curve_index] <- as.integer(old_source$passes)
  curve$source_converged[old_curve_index] <- as.logical(old_source$converged)
  curve$source_validation_loss_error[old_curve_index] <- abs(
    curve$baseline_validation_log_loss[old_curve_index] -
      old_source$validation_log_loss
  )
  curve$source_kkt_error[old_curve_index] <- abs(
    curve$baseline_kkt[old_curve_index] - old_source$kkt
  )
  curve$source_intercept_kkt_error[old_curve_index] <- abs(
    curve$baseline_intercept_kkt[old_curve_index] -
      old_source$intercept_score
  )
  curve$source_passes_match[old_curve_index] <-
    curve$baseline_passes[old_curve_index] == as.integer(old_source$passes)
  curve$source_converged_match[old_curve_index] <-
    curve$baseline_converged[old_curve_index] == as.logical(old_source$converged)
  endpoint_curve_index <- which(curve$point_type == "penalty_limit")
  if (length(endpoint_curve_index) != 1L || nrow(source_endpoint) != 1L) {
    stop("Task 125 analytic source endpoint is unavailable.", call. = FALSE)
  }
  curve$source_validation_log_loss[endpoint_curve_index] <-
    source_endpoint$validation_log_loss
  curve$source_intercept_kkt[endpoint_curve_index] <-
    source_endpoint$intercept_score
  curve$source_validation_loss_error[endpoint_curve_index] <-
    endpoint_loss_error
  curve$source_intercept_kkt_error[endpoint_curve_index] <- abs(
    curve$baseline_intercept_kkt[endpoint_curve_index] -
      source_endpoint$intercept_score
  )

  new_finite <- finite[finite$lambda_ratio > 512, , drop = FALSE]
  alpha_index <- match(case$alpha, sort(unique(configuration$alpha_grid)))
  d_index <- match(case$d, sort(unique(configuration$d_grid)))
  if (is.na(alpha_index) || is.na(d_index)) {
    stop("Task 125 alpha/d is absent from the frozen V7 grids.", call. = FALSE)
  }
  new_tuning <- data.frame(
    candidate_id = paste0("sgv8_tail:", new_finite$lambda_ratio),
    alpha_index = alpha_index, d_index = d_index,
    alpha = case$alpha, d = case$d,
    lambda_index = -seq_len(nrow(new_finite)),
    lambda = new_finite$lambda,
    lambda_relative_to_reference = new_finite$lambda_ratio,
    point_type = "finite",
    validation_log_loss = new_finite$baseline_validation_log_loss,
    numerically_eligible = new_finite$baseline_converged &
      is.finite(new_finite$baseline_kkt) &
      new_finite$baseline_kkt <= configuration$full_path_kkt_limit,
    stringsAsFactors = FALSE
  )
  combined <- lsg_bind_tuning_rows_v7(archived, new_tuning)
  selected_result <- lsg_select_sglasso_candidates_v7(
    combined, configuration$d_grid
  )
  selected <- selected_result$row
  selected_ratio <- selected$lambda_relative_to_reference[[1L]]
  range_resolved <- is.infinite(selected_ratio) ||
    !isTRUE(all.equal(selected_ratio, max(ratios), tolerance = 1e-12))
  selected_is_interior <- identical(
    as.character(selected$point_type[[1L]]), "finite"
  ) && is.finite(selected_ratio) &&
    selected_ratio > min(ratios) && selected_ratio < max(ratios)
  ordered <- combined[lsg_candidate_order_v7(combined), , drop = FALSE]
  eligible <- ordered[ordered$numerically_eligible %in% TRUE, , drop = FALSE]
  remaining <- eligible[eligible$candidate_id != selected$candidate_id[[1L]],
                        , drop = FALSE]
  next_best_loss <- if (nrow(remaining)) {
    min(remaining$validation_log_loss)
  } else NA_real_
  selection <- data.frame(
    case_id = case$case_id,
    source_selected_candidate_id = source_selected$candidate_id,
    source_selected_lambda_ratio =
      source_selected$lambda_relative_to_reference,
    source_selected_validation_log_loss = source_selected$validation_log_loss,
    diagnostic_selected_candidate_id = selected$candidate_id,
    diagnostic_selected_alpha = selected$alpha,
    diagnostic_selected_d = selected$d,
    diagnostic_selected_lambda_ratio = selected_ratio,
    diagnostic_selected_point_type = selected$point_type,
    diagnostic_selected_validation_log_loss = selected$validation_log_loss,
    diagnostic_selected_is_interior = selected_is_interior,
    diagnostic_next_best_validation_log_loss = next_best_loss,
    diagnostic_loss_gap_to_next_candidate = if (is.finite(next_best_loss)) {
      next_best_loss - selected$validation_log_loss[[1L]]
    } else NA_real_,
    diagnostic_range_resolved = range_resolved,
    largest_finite_lambda_ratio = max(ratios),
    lambda_reference_recomputed = fitted$lambda_reference,
    lambda_reference_archived = if (length(source_reference) == 1L) {
      source_reference
    } else NA_real_,
    lambda_reference_error = reference_error,
    old_v7_grid_embedded_exactly = old_grid_exact,
    ratio_512_anchor_loss_error = anchor_error,
    analytic_endpoint_loss_error = endpoint_loss_error,
    old_v7_grid_max_validation_loss_error =
      max(v7_replay$source_validation_loss_error),
    old_v7_grid_max_kkt_error =
      max(v7_replay$source_kkt_error),
    old_v7_grid_max_intercept_kkt_error =
      max(v7_replay$source_intercept_kkt_error),
    old_v7_grid_passes_exact =
      all(v7_replay$source_passes_match),
    old_v7_grid_convergence_exact =
      all(v7_replay$source_converged_match),
    extended_grid_max_validation_loss_difference_from_source =
      max(curve$source_validation_loss_error[old_curve_index]),
    extended_grid_max_kkt_difference_from_source =
      max(curve$source_kkt_error[old_curve_index]),
    stringsAsFactors = FALSE
  )
  list(
    schema_version = "sglasso_lambda_tail_failure_diagnostic_v8",
    case = case, curve = curve, v7_replay = v7_replay,
    selection = selection,
    data_sha256 = fingerprints, test_used = FALSE,
    raw = list(lambda_reference = fitted$lambda_reference,
               d_prefix = fitted$fit$d)
  )
}


lsg_grpreg_primary_reason_v8 <- function(diagnostic) {
  summary <- diagnostic$summary[1L, , drop = FALSE]
  if (nzchar(summary$fit_error)) return("fit_error")
  if (nzchar(summary$prediction_error)) return("prediction_error")
  if (!isTRUE(summary$raw_schema_compatible)) return("raw_schema_incompatible")
  if (isTRUE(summary$total_iteration_limit_reached)) {
    return("total_iteration_budget_reached")
  }
  if (!isTRUE(summary$path_complete)) return("incomplete_path")
  if (nzchar(summary$legacy_unexpected_warnings)) return("unexpected_solver_warning")
  if (summary$raw_finite_point_count < summary$returned_path_length) {
    return("nonfinite_returned_point")
  }
  if (!isTRUE(summary$legacy_path_termination_acceptable)) {
    return("legacy_path_policy_rejection")
  }
  if (summary$legacy_numerically_eligible_count == 0L) {
    return("no_point_marked_converged")
  }
  "original_failure_not_reproduced"
}


lsg_grpreg_reason_flags_v8 <- function(diagnostic) {
  summary <- diagnostic$summary[1L, , drop = FALSE]
  data.frame(
    fit_error = nzchar(summary$fit_error),
    prediction_error = nzchar(summary$prediction_error),
    raw_schema_incompatible = !isTRUE(summary$raw_schema_compatible),
    total_iteration_budget_reached =
      isTRUE(summary$total_iteration_limit_reached),
    incomplete_path = !isTRUE(summary$path_complete),
    unexpected_solver_warning =
      nzchar(summary$legacy_unexpected_warnings),
    nonfinite_returned_point =
      is.finite(summary$raw_finite_point_count) &&
      summary$raw_finite_point_count < summary$returned_path_length,
    legacy_path_policy_rejection =
      !isTRUE(summary$legacy_path_termination_acceptable),
    no_legacy_eligible_point =
      is.finite(summary$legacy_numerically_eligible_count) &&
      summary$legacy_numerically_eligible_count == 0L,
    usable_finite_prefix_available =
      is.finite(summary$usable_prefix_length) &&
      summary$usable_prefix_length > 0L,
    stringsAsFactors = FALSE
  )
}


lsg_run_grpreg_case_v8 <- function(root, case) {
  diagnostic <- lsg_run_grpreg_task_v8(root, case$task_id)
  points <- diagnostic$points
  summary <- diagnostic$summary[1L, , drop = FALSE]
  legacy_tuning <- if (nrow(points)) data.frame(
    engine = "grpreg",
    method_path = "Logistic Group Lasso (grpreg)",
    penalty_family = "group_lasso",
    fit_index = 1L,
    lambda_index = points$lambda_index,
    validation_log_loss = points$validation_log_loss,
    solver_warning = summary$warnings,
    unexpected_solver_warning = summary$legacy_unexpected_warnings,
    requested_path_length = summary$requested_path_length,
    returned_path_length = summary$returned_path_length,
    path_complete = summary$path_complete,
    saturated_path_truncation = FALSE,
    iteration_budget_truncation = FALSE,
    total_iteration_limit_reached = summary$total_iteration_limit_reached,
    returned_lower_boundary_excluded =
      points$legacy_returned_lower_boundary_excluded,
    path_termination_acceptable =
      summary$legacy_path_termination_acceptable,
    point_converged = points$legacy_point_converged,
    finite_validation_loss = points$finite_validation_loss,
    numerically_eligible = points$legacy_numerically_eligible,
    stringsAsFactors = FALSE
  ) else data.frame()
  selection_error <- ""
  selection <- if (nrow(legacy_tuning)) tryCatch(
    select_valid_external_grid_v2(
      legacy_tuning, "Logistic Group Lasso (grpreg)"
    ),
    error = function(condition) {
      selection_error <<- conditionMessage(condition)
      NULL
    }
  ) else {
    selection_error <- "No returned Group Lasso path was available for selection."
    NULL
  }
  diagnostic$legacy_tuning <- legacy_tuning
  diagnostic$selection <- selection
  diagnostic$selection_error <- selection_error
  diagnostic$raw_evidence_constructed_before_selection <- TRUE
  diagnostic$reason_flags <- lsg_grpreg_reason_flags_v8(diagnostic)
  diagnostic$summary$case_id <- case$case_id
  diagnostic$summary$primary_reason <- lsg_grpreg_primary_reason_v8(diagnostic)
  diagnostic$summary$selection_error <- selection_error
  if (nrow(diagnostic$points)) diagnostic$points$case_id <- case$case_id
  diagnostic$test_used <- FALSE
  diagnostic
}


lsg_intercept_profile_unit_checks_v8 <- function(root) {
  lsg_compile_intercept_profile_v8(root)
  set.seed(8100801L)
  X <- matrix(stats::rnorm(240), nrow = 80L, ncol = 3L)
  beta <- c(1.2, -0.7, 0.4)
  y <- c(rep(0, 68L), rep(1, 12L))
  initial <- -1
  profiled <- lsg_profile_intercept_beta_v8_cpp(
    X, y, beta, initial_intercept = initial
  )
  offset <- drop(X %*% beta)
  profiled_offset <- lsg_profile_intercept_offset_v8_cpp(
    offset, y, initial_intercept = initial
  )
  reference <- stats::uniroot(
    function(b) mean(stats::plogis(b + offset)) - mean(y),
    c(stats::qlogis(mean(y)) - max(offset),
      stats::qlogis(mean(y)) - min(offset)),
    tol = 1e-12
  )$root
  independent_score <- abs(
    mean(stats::plogis(profiled$intercept + offset)) - mean(y)
  )
  checks <- c(
    cpp_profile_converged = isTRUE(profiled$converged),
    cpp_profile_matches_uniroot = abs(profiled$intercept - reference) <= 1e-10,
    cpp_beta_and_offset_interfaces_agree =
      isTRUE(profiled_offset$converged) &&
      abs(profiled$intercept - profiled_offset$intercept) <= 1e-13 &&
      max(abs(profiled$probability - profiled_offset$probability)) <= 1e-13,
    cpp_profile_closes_independent_intercept_score =
      independent_score <= 1e-12 &&
      abs(independent_score - profiled$intercept_kkt) <= 1e-13,
    cpp_profile_does_not_increase_loss = profiled$log_loss_decrease >= -1e-14,
    cpp_profile_probability_dimensions =
      length(profiled$probability) == length(y) &&
      all(is.finite(profiled$probability))
  )
  data.frame(check = names(checks), passed = unname(checks),
             stringsAsFactors = FALSE)
}
