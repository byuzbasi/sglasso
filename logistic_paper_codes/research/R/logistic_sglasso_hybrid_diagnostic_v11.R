# Three-case, training/validation-only diagnostic for the V11 hybrid solver.
# The cases are the frozen V7/V10 SGLASSO failures sg_near,
# sg_lower_median, and sg_worst. No test outcome or truth field is exposed.

lsg_hybrid_diagnostic_defaults_v11 <- function() {
  solver <- lsg_hybrid_defaults_v11()
  list(
    schema_version = "logistic_sglasso_hybrid_diagnostic_v11",
    default_version = "logistic_sglasso_hybrid_diagnostic_v11_2",
    source_v7_version = "logistic_sglasso_prediction_selection_pilot_r20_v7",
    source_v7_signature =
      "93ce8540e29856c38b2fa0ea144a038c68e1b80ceb4ea5adca3094c34f7850ec",
    case_ids = c("sg_near", "sg_lower_median", "sg_worst"),
    configuration = list(
      lambda_relative_grid = lsg_lambda_relative_grid_v10(),
      nlambda = 63L,
      alpha_grid = seq(0, 1, by = 0.1),
      d_grid = seq(0, 1, by = 0.1),
      path_controls = solver$path,
      polish_controls = solver$polish,
      objective_reconstruction_ulp_factor = 32,
      reconstruction_tolerance = 2e-6,
      stability_objective_tolerance = 1e-7,
      stability_validation_loss_tolerance = 1e-6,
      stability_prediction_tolerance = 1e-4,
      test_sample_policy = "forbidden",
      d_path_policy = "frozen_case_d_only",
      lower_lambda_boundary_policy = "report_not_silently_extend"
    )
  )
}


lsg_hybrid_context_v11 <- function(case, frame) {
  frame$case_id <- case$case_id
  frame$task_id <- case$task_id
  frame$scenario <- case$scenario
  frame$replication <- case$replication
  frame$seed <- case$seed
  frame[c(
    "case_id", "task_id", "scenario", "replication", "seed",
    setdiff(names(frame), c(
      "case_id", "task_id", "scenario", "replication", "seed"
    ))
  )]
}


lsg_fit_one_hybrid_with_controls_v11 <- function(
    preprocess, y, target, lambda, alpha, d, beta, intercept, controls,
    keep_trace = FALSE
) {
  controls <- lsg_validate_hybrid_controls_v11(controls)
  lsg_fit_one_hybrid_v11_cpp(
    preprocess$X, y, preprocess$group_start, preprocess$group_end,
    preprocess$group_weight, target, lambda, alpha, d, beta, intercept,
    controls$block_max_sweeps, controls$block_chunk_sweeps,
    controls$block_stall_window,
    controls$block_stall_relative_improvement,
    controls$apg_max_iterations, controls$apg_kkt_check_interval,
    controls$kkt_tolerance, controls$update_tolerance,
    controls$intercept_tolerance, controls$max_intercept_iterations,
    TRUE, TRUE, isTRUE(keep_trace)
  )
}


lsg_hybrid_selected_stability_v11 <- function(
    data,
    preprocess,
    target_original,
    alpha,
    d,
    lambda,
    lambda_ratio,
    lambda_index,
    all_lambda,
    forward_beta,
    forward_intercept,
    forward_objective,
    forward_validation_probability,
    configuration
) {
  target <- project_lsg_original_target(preprocess, target_original)
  null_intercept <- stats::qlogis(mean(data$y_train))
  controls <- configuration$polish_controls
  forward <- lsg_fit_one_hybrid_with_controls_v11(
    preprocess, data$y_train, target, lambda, alpha, d,
    forward_beta, forward_intercept, controls
  )
  cold <- lsg_fit_one_hybrid_with_controls_v11(
    preprocess, data$y_train, target, lambda, alpha, d,
    rep(0, ncol(preprocess$X)), null_intercept, controls
  )

  reverse_lambda <- rev(all_lambda[seq.int(lambda_index, length(all_lambda))])
  reverse_path <- lsg_fit_hybrid_path_v11(
    data$X_train, data$y_train, data$group,
    lambda = reverse_lambda, d = d, alpha = alpha,
    target_original = target_original, preprocess = preprocess,
    controls = configuration$path_controls, use_active_set = TRUE,
    enable_fallback = TRUE, warm_start_d = FALSE, keep_traces = FALSE,
    lambda_order = "any", compile = FALSE
  )
  reverse_seed_beta <- reverse_path$beta_solver[, length(reverse_lambda), 1L]
  reverse_seed_intercept <- reverse_path$intercept_solver[
    length(reverse_lambda), 1L
  ]
  reverse <- lsg_fit_one_hybrid_with_controls_v11(
    preprocess, data$y_train, target, lambda, alpha, d,
    reverse_seed_beta, reverse_seed_intercept, controls
  )

  validation_solver <- transform_lsg_newx(preprocess, data$X_validation)
  probability <- function(fit) {
    stats::plogis(fit$intercept + drop(validation_solver %*% fit$beta))
  }
  forward_probability <- probability(forward)
  cold_probability <- probability(cold)
  reverse_probability <- probability(reverse)
  loss <- c(
    forward = binary_log_loss(data$y_validation, forward_probability),
    cold = binary_log_loss(data$y_validation, cold_probability),
    reverse = binary_log_loss(data$y_validation, reverse_probability)
  )
  objective <- c(
    forward = forward$objective,
    cold = cold$objective,
    reverse = reverse$objective
  )
  objective_error <- max(abs(outer(objective, objective, "-")))
  validation_loss_error <- max(abs(outer(loss, loss, "-")))
  prediction_error <- max(
    abs(forward_probability - cold_probability),
    abs(forward_probability - reverse_probability),
    abs(cold_probability - reverse_probability)
  )
  forward_shift <- abs(
    binary_log_loss(data$y_validation, forward_validation_probability) -
      loss[["forward"]]
  )
  checks <- c(
    all_three_starts_converged = all(c(
      forward$converged, cold$converged, reverse$converged
    )),
    all_three_starts_pass_polish_kkt = max(c(
      forward$kkt, cold$kkt, reverse$kkt
    )) <= controls$study_kkt_limit,
    objectives_agree = objective_error <=
      configuration$stability_objective_tolerance,
    validation_losses_agree = validation_loss_error <=
      configuration$stability_validation_loss_tolerance,
    validation_probabilities_agree = prediction_error <=
      configuration$stability_prediction_tolerance,
    forward_polish_is_stable =
      abs(forward$objective - forward_objective) <=
        configuration$stability_objective_tolerance &&
      forward_shift <= configuration$stability_validation_loss_tolerance &&
      max(abs(forward_probability - forward_validation_probability)) <=
        configuration$stability_prediction_tolerance
  )
  data.frame(
    selected_lambda = lambda,
    selected_lambda_ratio = lambda_ratio,
    selected_lambda_index = lambda_index,
    forward_kkt_before_polish = lsg_kkt_cpp(
      preprocess$X, data$y_train, forward_beta, forward_intercept,
      preprocess$group_start, preprocess$group_end,
      preprocess$group_weight, target, lambda, alpha, d
    )$maximum,
    forward_kkt = forward$kkt,
    cold_kkt = cold$kkt,
    reverse_kkt = reverse$kkt,
    forward_converged = isTRUE(forward$converged),
    cold_converged = isTRUE(cold$converged),
    reverse_converged = isTRUE(reverse$converged),
    forward_solver_route = forward$solver_route,
    cold_solver_route = cold$solver_route,
    reverse_solver_route = reverse$solver_route,
    maximum_objective_error = objective_error,
    maximum_validation_loss_error = validation_loss_error,
    maximum_validation_probability_error = prediction_error,
    forward_polish_objective_shift = abs(
      forward$objective - forward_objective
    ),
    forward_polish_validation_loss_shift = forward_shift,
    forward_polish_probability_shift = max(abs(
      forward_probability - forward_validation_probability
    )),
    polished_validation_log_loss = loss[["forward"]],
    accepted = all(checks),
    stringsAsFactors = FALSE
  )
}


lsg_run_hybrid_case_v11 <- function(
    case, scenario, source, source_configuration, configuration, data = NULL
) {
  required <- c(
    "lsg_v8_tuning_data", "lsg_v8_data_fingerprints",
    "lsg_tail_reference_v6", "lsg_tail_limit_fit_v6",
    "lsg_fit_hybrid_path_v11", "lsg_fit_one_hybrid_v11_cpp",
    "binary_log_loss", "lsg_kkt_cpp", "lsg_objective_cpp",
    "lsg_objective_reconstruction_audit_v10"
  )
  missing <- required[!vapply(
    required, exists, logical(1), mode = "function",
    envir = environment(lsg_run_hybrid_case_v11), inherits = TRUE
  )]
  if (length(missing)) {
    stop("Missing V11 diagnostic dependency: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  if (!is.data.frame(case) || nrow(case) != 1L ||
      !case$case_id %in% lsg_hybrid_diagnostic_defaults_v11()$case_ids) {
    stop("One frozen V11 SGLASSO case is required.", call. = FALSE)
  }
  # Optional exact input transport for local replay across R/platform versions.
  # It cannot bypass the same field and frozen fingerprint checks below.
  if (is.null(data)) data <- lsg_v8_tuning_data(scenario, case$seed)
  if (!identical(names(data), c(
      "X_train", "y_train", "X_validation", "y_validation", "group"
  ))) {
    stop("V11 diagnostic data contain a forbidden field.", call. = FALSE)
  }
  fingerprints <- lsg_v8_data_fingerprints(data)
  if (!identical(fingerprints, source$data_sha256)) {
    stop("Regenerated V11 data differ from frozen V7 evidence.",
         call. = FALSE)
  }
  target_original <- as.numeric(source$firth_target_original)
  if (length(target_original) != ncol(data$X_train) ||
      any(!is.finite(target_original))) {
    stop("The blinded V7 Firth target is unavailable for V11.",
         call. = FALSE)
  }
  preprocess <- prepare_lsg_design(data$X_train, data$group)
  target <- project_lsg_original_target(preprocess, target_original)
  lambda_reference <- lsg_tail_reference_v6(
    preprocess, data$y_train, target_original, case$alpha
  )
  ratios <- configuration$lambda_relative_grid
  lambda <- lambda_reference * ratios
  fit <- lsg_fit_hybrid_path_v11(
    data$X_train, data$y_train, data$group, lambda, case$d,
    case$alpha, target_original, preprocess = preprocess,
    controls = configuration$path_controls, use_active_set = TRUE,
    enable_fallback = TRUE, warm_start_d = FALSE, keep_traces = FALSE,
    compile = FALSE
  )
  validation_solver <- transform_lsg_newx(preprocess, data$X_validation)
  rows <- vector("list", length(ratios))
  validation_probabilities <- vector("list", length(ratios))
  for (li in seq_along(ratios)) {
    beta <- fit$beta_solver[, li, 1L]
    intercept <- fit$intercept_solver[li, 1L]
    probability <- stats::plogis(
      intercept + drop(validation_solver %*% beta)
    )
    validation_probabilities[[li]] <- probability
    original <- recover_lsg_coefficients(preprocess, beta, intercept)
    reconstruction <- max(abs(
      probability - stats::plogis(
        original[1L] + drop(data$X_validation %*% original[-1L])
      )
    ))
    objective_reference <- lsg_objective_cpp(
      preprocess$X, data$y_train, beta, intercept,
      preprocess$group_start, preprocess$group_end,
      preprocess$group_weight, target, lambda[li], case$alpha, case$d
    )
    objective_audit <- lsg_objective_reconstruction_audit_v10(
      fit$objective[li, 1L], objective_reference,
      configuration$objective_reconstruction_ulp_factor
    )
    kkt_reference <- lsg_kkt_cpp(
      preprocess$X, data$y_train, beta, intercept,
      preprocess$group_start, preprocess$group_end,
      preprocess$group_weight, target, lambda[li], case$alpha, case$d
    )
    rows[[li]] <- data.frame(
      point_type = "finite",
      lambda_index = li,
      lambda_ratio = ratios[li],
      lambda = lambda[li],
      alpha = case$alpha,
      d = case$d,
      validation_log_loss = binary_log_loss(data$y_validation, probability),
      objective = fit$objective[li, 1L],
      objective_reconstruction_error =
        objective_audit$objective_reconstruction_error,
      objective_reconstruction_ulp_units =
        objective_audit$objective_reconstruction_ulp_units,
      objective_reconstruction_pass =
        objective_audit$objective_reconstruction_pass,
      kkt = fit$kkt[li, 1L],
      kkt_reconstruction_error = abs(
        fit$kkt[li, 1L] - kkt_reference$maximum
      ),
      intercept_kkt = fit$intercept_kkt[li, 1L],
      group_kkt = fit$group_kkt[li, 1L],
      converged = isTRUE(fit$converged[li, 1L]),
      solver_route = fit$solver_route[li, 1L],
      fallback_used = isTRUE(fit$fallback_used[li, 1L]),
      block_stalled = isTRUE(fit$block_stalled[li, 1L]),
      block_termination = fit$block_termination[li, 1L],
      termination_reason = fit$termination_reason[li, 1L],
      block_sweeps = fit$block_sweeps[li, 1L],
      apg_iterations = fit$apg_iterations[li, 1L],
      selected_groups = fit$selected_groups[li, 1L],
      reconstruction_error = reconstruction,
      numerically_eligible = isTRUE(fit$converged[li, 1L]) &&
        is.finite(fit$kkt[li, 1L]) &&
        fit$kkt[li, 1L] <=
          configuration$path_controls$study_kkt_limit &&
        isTRUE(objective_audit$objective_reconstruction_pass) &&
        is.finite(reconstruction) &&
        reconstruction <= configuration$reconstruction_tolerance,
      stringsAsFactors = FALSE
    )
  }
  finite <- do.call(rbind, rows)
  endpoint <- lsg_tail_limit_fit_v6(
    data, preprocess, target_original, case$alpha, case$d,
    source_configuration
  )
  endpoint_row <- data.frame(
    point_type = "penalty_limit", lambda_index = 0L,
    lambda_ratio = Inf, lambda = Inf, alpha = case$alpha, d = case$d,
    validation_log_loss = endpoint$point$validation_log_loss,
    objective = NA_real_, objective_reconstruction_error = 0,
    objective_reconstruction_ulp_units = 0,
    objective_reconstruction_pass = TRUE, kkt = NA_real_,
    kkt_reconstruction_error = 0,
    intercept_kkt = endpoint$point$intercept_score,
    group_kkt = endpoint$point$penalty_kkt,
    converged = endpoint$point$numerically_eligible,
    solver_route = "analytic_penalty_limit", fallback_used = FALSE,
    block_stalled = FALSE, block_termination = "analytic_penalty_limit",
    termination_reason = "analytic_penalty_limit", block_sweeps = 0L,
    apg_iterations = 0L, selected_groups = endpoint$point$selected_groups,
    reconstruction_error = endpoint$point$reconstruction_error,
    numerically_eligible = endpoint$point$numerically_eligible,
    stringsAsFactors = FALSE
  )
  candidates <- rbind(endpoint_row, finite)
  eligible <- which(
    candidates$numerically_eligible &
      is.finite(candidates$validation_log_loss)
  )
  if (!length(eligible)) {
    stop("V11 hybrid path has no numerically eligible candidate.",
         call. = FALSE)
  }
  selected_row <- eligible[which.min(candidates$validation_log_loss[eligible])]
  candidates$selected <- FALSE
  candidates$selected[selected_row] <- TRUE
  selected <- candidates[selected_row, , drop = FALSE]
  invalid_loss <- candidates$validation_log_loss[
    !candidates$numerically_eligible &
      is.finite(candidates$validation_log_loss)
  ]
  invalid_competitive <- length(invalid_loss) &&
    min(invalid_loss) <= selected$validation_log_loss
  finite_selected <- selected$point_type == "finite"
  selected_lower_boundary <- finite_selected &&
    selected$lambda_index == length(ratios)

  stability <- if (finite_selected) {
    li <- selected$lambda_index
    lsg_hybrid_selected_stability_v11(
      data, preprocess, target_original, case$alpha, case$d,
      selected$lambda, selected$lambda_ratio, li, lambda,
      fit$beta_solver[, li, 1L], fit$intercept_solver[li, 1L],
      fit$objective[li, 1L], validation_probabilities[[li]], configuration
    )
  } else {
    data.frame(
      selected_lambda = Inf, selected_lambda_ratio = Inf,
      selected_lambda_index = 0L,
      forward_kkt_before_polish = NA_real_, forward_kkt = NA_real_,
      cold_kkt = NA_real_, reverse_kkt = NA_real_,
      forward_converged = TRUE, cold_converged = TRUE,
      reverse_converged = TRUE,
      forward_solver_route = "analytic_penalty_limit",
      cold_solver_route = "analytic_penalty_limit",
      reverse_solver_route = "analytic_penalty_limit",
      maximum_objective_error = 0, maximum_validation_loss_error = 0,
      maximum_validation_probability_error = 0,
      forward_polish_objective_shift = 0,
      forward_polish_validation_loss_shift = 0,
      forward_polish_probability_shift = 0,
      polished_validation_log_loss = selected$validation_log_loss,
      accepted = TRUE, stringsAsFactors = FALSE
    )
  }

  checks <- c(
    test_fields_forbidden =
      identical(configuration$test_sample_policy, "forbidden") &&
      identical(names(data), c(
        "X_train", "y_train", "X_validation", "y_validation", "group"
      )),
    frozen_training_validation_fingerprints_exact =
      identical(fingerprints, source$data_sha256),
    exact_frozen_63_point_grid = nrow(finite) == 63L &&
      isTRUE(all.equal(finite$lambda_ratio,
                       configuration$lambda_relative_grid,
                       tolerance = 1e-13)),
    all_finite_outputs = all(is.finite(c(
      finite$lambda, finite$validation_log_loss, finite$objective,
      finite$kkt, finite$intercept_kkt, finite$group_kkt,
      finite$reconstruction_error
    ))),
    all_finite_points_numerically_eligible =
      all(finite$numerically_eligible %in% TRUE),
    objective_reconstruction_within_32_ulp =
      all(finite$objective_reconstruction_pass %in% TRUE) &&
      max(finite$objective_reconstruction_ulp_units) <=
        configuration$objective_reconstruction_ulp_factor,
    kkt_reconstructed = max(finite$kkt_reconstruction_error) <= 3e-12,
    intercept_profiled = max(finite$intercept_kkt) <= 1e-10,
    prediction_reconstructed = max(finite$reconstruction_error) <=
      configuration$reconstruction_tolerance,
    selected_candidate_numerically_eligible =
      isTRUE(selected$numerically_eligible),
    no_competitive_invalid_candidate = !invalid_competitive,
    selected_solution_three_start_stable = isTRUE(stability$accepted)
  )
  outcomes <- data.frame(
    outcome = c(
      "selected_lower_lambda_boundary",
      "lambda_range_resolved_without_extension",
      "fallback_used_on_any_finite_point",
      "block_stall_detected_on_any_finite_point"
    ),
    value = c(
      selected_lower_boundary,
      !selected_lower_boundary,
      any(finite$fallback_used),
      any(finite$block_stalled)
    ),
    interpretation = c(
      "Scientific range diagnostic; not a numerical convergence failure",
      "TRUE only when the selected finite lambda is not the lower endpoint",
      "The general V11 fallback was activated on at least one path point",
      "The KKT progress-window rule activated on at least one path point"
    ),
    stringsAsFactors = FALSE
  )
  summary <- data.frame(
    alpha = case$alpha, d = case$d, finite_points = nrow(finite),
    eligible_finite_points = sum(finite$numerically_eligible),
    fallback_finite_points = sum(finite$fallback_used),
    stalled_block_points = sum(finite$block_stalled),
    selected_point_type = selected$point_type,
    selected_lambda_index = selected$lambda_index,
    selected_lambda_ratio = selected$lambda_ratio,
    selected_validation_log_loss = selected$validation_log_loss,
    polished_selected_validation_log_loss =
      stability$polished_validation_log_loss,
    selected_lower_lambda_boundary = selected_lower_boundary,
    invalid_candidate_competitive = invalid_competitive,
    maximum_finite_kkt = max(finite$kkt),
    maximum_objective_reconstruction_ulp_units =
      max(finite$objective_reconstruction_ulp_units),
    maximum_reconstruction_error = max(finite$reconstruction_error),
    stability_accepted = stability$accepted,
    stringsAsFactors = FALSE
  )
  list(
    schema_version = "logistic_sglasso_hybrid_case_v11",
    case = case,
    finite_path = lsg_hybrid_context_v11(case, finite),
    candidates = lsg_hybrid_context_v11(case, candidates),
    selected_stability = lsg_hybrid_context_v11(case, stability),
    summary = lsg_hybrid_context_v11(case, summary),
    hard_checks = data.frame(
      case_id = case$case_id, task_id = case$task_id,
      check = names(checks), passed = unname(checks),
      stringsAsFactors = FALSE
    ),
    scientific_outcomes = lsg_hybrid_context_v11(case, outcomes),
    data_sha256 = fingerprints,
    evaluation_data_used = FALSE
  )
}
