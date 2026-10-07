# Seven-case, training/validation-only validation of the V9 numerical repair.
# The cases and their V7 source evidence are inherited exactly from V8.

lsg_repair_defaults_v9 <- function() {
  list(
    schema_version = "logistic_sglasso_numerical_repair_v9",
    default_version = "logistic_sglasso_numerical_repair_diagnostic_v9",
    source_v7_version = "logistic_sglasso_prediction_selection_pilot_r20_v7",
    source_v7_signature =
      "93ce8540e29856c38b2fa0ea144a038c68e1b80ceb4ea5adca3094c34f7850ec",
    source_v8_signature =
      "7aa7219fa200d539f31c78bb7ece10e49a0571646753b70d584e2c64ff400dca",
    configuration = list(
      lambda_relative_grid = lsg_lambda_relative_grid_v9(),
      nlambda = 39L,
      alpha_grid = seq(0, 1, by = 0.1),
      d_grid = seq(0, 1, by = 0.1),
      joint_max_sweeps = 4000L,
      joint_kkt_tolerance = 2e-6,
      study_kkt_limit = 2.05e-6,
      joint_update_tolerance = 1e-10,
      joint_intercept_tolerance = 1e-12,
      joint_max_intercept_iterations = 100L,
      reconstruction_tolerance = 2e-6,
      stability_objective_tolerance = 1e-7,
      stability_validation_loss_tolerance = 1e-6,
      stability_prediction_tolerance = 1e-4,
      accepted_objective_increase_tolerance = 5e-10,
      grpreg_max_iterations = 1000000L,
      grpreg_requested_path_length = 30L,
      test_sample_policy = "forbidden"
    )
  )
}


lsg_repair_require_v9 <- function() {
  required <- c(
    "lsg_v8_tuning_data", "lsg_v8_data_fingerprints", "lsg_v8_d_prefix",
    "lsg_tail_reference_v6", "lsg_tail_limit_fit_v6",
    "lsg_fit_joint_path_v9", "lsg_fit_one_joint_v9_cpp",
    "lsg_classify_grpreg_path_v9", "lsg_select_grpreg_prefix_v9",
    "lsg_run_grpreg_task_v8", "binary_log_loss", "lsg_kkt_cpp",
    "lsg_objective_cpp"
  )
  envir <- environment(lsg_repair_require_v9)
  missing <- required[!vapply(required, exists, logical(1), mode = "function",
                              envir = envir, inherits = TRUE)]
  if (length(missing)) {
    stop("Missing V9 numerical-repair dependency: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}


lsg_repair_context_v9 <- function(case, frame) {
  frame$case_id <- case$case_id
  frame$task_id <- case$task_id
  frame$scenario <- case$scenario
  frame$replication <- case$replication
  frame$seed <- case$seed
  frame[c("case_id", "task_id", "scenario", "replication", "seed",
          setdiff(names(frame), c(
            "case_id", "task_id", "scenario", "replication", "seed"
          )))]
}


lsg_repair_selected_stability_v9 <- function(
    data,
    preprocess,
    target_original,
    alpha,
    d,
    lambda,
    lambda_ratio,
    forward_beta,
    forward_intercept,
    forward_objective,
    forward_validation_probability,
    all_lambda,
    configuration
) {
  target <- project_lsg_original_target(preprocess, target_original)
  null_intercept <- stats::qlogis(mean(data$y_train))
  cold <- lsg_fit_one_joint_v9_cpp(
    preprocess$X, data$y_train, preprocess$group_start,
    preprocess$group_end, preprocess$group_weight, target,
    lambda, alpha, d, rep(0, ncol(preprocess$X)), null_intercept,
    configuration$joint_max_sweeps,
    configuration$joint_kkt_tolerance,
    configuration$joint_update_tolerance,
    configuration$joint_intercept_tolerance,
    configuration$joint_max_intercept_iterations,
    TRUE, FALSE
  )
  reverse <- lsg_fit_joint_path_v9(
    data$X_train, data$y_train, data$group,
    lambda = rev(all_lambda), d = d, alpha = alpha,
    target_original = target_original, preprocess = preprocess,
    max_sweeps = configuration$joint_max_sweeps,
    kkt_tolerance = configuration$joint_kkt_tolerance,
    update_tolerance = configuration$joint_update_tolerance,
    intercept_tolerance = configuration$joint_intercept_tolerance,
    max_intercept_iterations =
      configuration$joint_max_intercept_iterations,
    use_active_set = TRUE, warm_start_d = FALSE, keep_traces = FALSE,
    lambda_order = "any", compile = FALSE
  )
  reverse_index <- which(abs(reverse$lambda - lambda) <=
    1e-12 * pmax(1, abs(lambda)))
  if (length(reverse_index) != 1L) {
    stop("V9 reverse path does not contain the selected lambda.",
         call. = FALSE)
  }
  X_validation_solver <- transform_lsg_newx(preprocess, data$X_validation)
  cold_probability <- stats::plogis(
    cold$intercept + drop(X_validation_solver %*% cold$beta)
  )
  reverse_probability <- stats::plogis(
    reverse$intercept_solver[reverse_index, 1L] +
      drop(X_validation_solver %*%
             reverse$beta_solver[, reverse_index, 1L])
  )
  cold_loss <- binary_log_loss(data$y_validation, cold_probability)
  reverse_loss <- binary_log_loss(data$y_validation, reverse_probability)
  objective_error <- max(abs(c(
    cold$objective - forward_objective,
    reverse$objective[reverse_index, 1L] - forward_objective
  )))
  validation_loss_error <- max(abs(c(
    cold_loss - binary_log_loss(
      data$y_validation, forward_validation_probability
    ),
    reverse_loss - binary_log_loss(
      data$y_validation, forward_validation_probability
    )
  )))
  prediction_error <- max(
    abs(cold_probability - forward_validation_probability),
    abs(reverse_probability - forward_validation_probability)
  )
  data.frame(
    selected_lambda = lambda,
    selected_lambda_ratio = lambda_ratio,
    forward_kkt = lsg_kkt_cpp(
      preprocess$X, data$y_train, forward_beta, forward_intercept,
      preprocess$group_start, preprocess$group_end,
      preprocess$group_weight, target, lambda, alpha, d
    )$maximum,
    cold_kkt = cold$kkt,
    reverse_kkt = reverse$kkt[reverse_index, 1L],
    cold_converged = isTRUE(cold$converged),
    reverse_converged = isTRUE(reverse$converged[reverse_index, 1L]),
    maximum_objective_error = objective_error,
    maximum_validation_loss_error = validation_loss_error,
    maximum_validation_probability_error = prediction_error,
    accepted = isTRUE(cold$converged) &&
      isTRUE(reverse$converged[reverse_index, 1L]) &&
      is.finite(objective_error) &&
      objective_error <= configuration$stability_objective_tolerance &&
      validation_loss_error <=
        configuration$stability_validation_loss_tolerance &&
      prediction_error <= configuration$stability_prediction_tolerance,
    stringsAsFactors = FALSE
  )
}


lsg_run_sglasso_repair_case_v9 <- function(
    case,
    scenario,
    source,
    source_configuration,
    configuration
) {
  lsg_repair_require_v9()
  if (!is.data.frame(case) || nrow(case) != 1L ||
      !case$diagnostic_type %in% c(
        "sglasso_intercept", "sglasso_lambda_tail"
      )) {
    stop("One frozen V9 SGLASSO repair case is required.", call. = FALSE)
  }
  data <- lsg_v8_tuning_data(scenario, case$seed)
  fingerprints <- lsg_v8_data_fingerprints(data)
  if (!identical(fingerprints, source$data_sha256)) {
    stop("Regenerated V9 tuning data differ from the frozen V7 evidence.",
         call. = FALSE)
  }
  target_original <- as.numeric(source$firth_target_original)
  if (length(target_original) != ncol(data$X_train) ||
      any(!is.finite(target_original))) {
    stop("The blinded V7 Firth target is unavailable for V9.",
         call. = FALSE)
  }
  preprocess <- prepare_lsg_design(data$X_train, data$group)
  lambda_reference <- lsg_tail_reference_v6(
    preprocess, data$y_train, target_original, case$alpha
  )
  ratios <- configuration$lambda_relative_grid
  lambda <- lambda_reference * ratios
  d_prefix <- lsg_v8_d_prefix(source_configuration, case$alpha, case$d)
  fit <- lsg_fit_joint_path_v9(
    data$X_train, data$y_train, data$group, lambda, d_prefix,
    case$alpha, target_original, preprocess = preprocess,
    max_sweeps = configuration$joint_max_sweeps,
    kkt_tolerance = configuration$joint_kkt_tolerance,
    update_tolerance = configuration$joint_update_tolerance,
    intercept_tolerance = configuration$joint_intercept_tolerance,
    max_intercept_iterations =
      configuration$joint_max_intercept_iterations,
    use_active_set = TRUE, warm_start_d = TRUE, keep_traces = FALSE,
    compile = FALSE
  )
  d_index <- which(abs(fit$d - case$d) <= 1e-12)
  if (length(d_index) != 1L) {
    stop("The V9 fitted d prefix does not contain the frozen target d.",
         call. = FALSE)
  }
  X_validation_solver <- transform_lsg_newx(
    preprocess, data$X_validation
  )
  target_solver <- project_lsg_original_target(
    preprocess, target_original
  )
  rows <- vector("list", length(ratios))
  validation_probabilities <- vector("list", length(ratios))
  for (li in seq_along(ratios)) {
    beta <- fit$beta_solver[, li, d_index]
    intercept <- fit$intercept_solver[li, d_index]
    probability <- stats::plogis(
      intercept + drop(X_validation_solver %*% beta)
    )
    validation_probabilities[[li]] <- probability
    original <- recover_lsg_coefficients(preprocess, beta, intercept)
    reconstruction <- max(abs(
      probability - stats::plogis(
        original[1L] + drop(data$X_validation %*% original[-1L])
      )
    ))
    frozen_kkt <- lsg_kkt_cpp(
      preprocess$X, data$y_train, beta, intercept,
      preprocess$group_start, preprocess$group_end,
      preprocess$group_weight, target_solver, lambda[li],
      case$alpha, case$d
    )
    frozen_objective <- lsg_objective_cpp(
      preprocess$X, data$y_train, beta, intercept,
      preprocess$group_start, preprocess$group_end,
      preprocess$group_weight, target_solver, lambda[li],
      case$alpha, case$d
    )
    rows[[li]] <- data.frame(
      point_type = "finite",
      lambda_index = li,
      lambda_ratio = ratios[li],
      lambda = lambda[li],
      alpha = case$alpha,
      d = case$d,
      validation_log_loss = binary_log_loss(
        data$y_validation, probability
      ),
      objective = fit$objective[li, d_index],
      objective_reconstruction_error = abs(
        fit$objective[li, d_index] - frozen_objective
      ),
      kkt = fit$kkt[li, d_index],
      kkt_reconstruction_error = abs(
        fit$kkt[li, d_index] - frozen_kkt$maximum
      ),
      intercept_kkt = fit$intercept_kkt[li, d_index],
      group_kkt = fit$group_kkt[li, d_index],
      converged = isTRUE(fit$converged[li, d_index]),
      sweeps = fit$sweeps[li, d_index],
      termination_reason = fit$termination_reason[li, d_index],
      selected_groups = fit$selected_groups[li, d_index],
      reconstruction_error = reconstruction,
      maximum_accepted_objective_increase =
        fit$maximum_accepted_objective_increase[li, d_index],
      numerically_eligible = isTRUE(fit$converged[li, d_index]) &&
        is.finite(fit$kkt[li, d_index]) &&
        fit$kkt[li, d_index] <= configuration$study_kkt_limit &&
        is.finite(reconstruction) &&
        reconstruction <= configuration$reconstruction_tolerance &&
        is.finite(binary_log_loss(data$y_validation, probability)),
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
    kkt = NA_real_, kkt_reconstruction_error = 0,
    intercept_kkt = endpoint$point$intercept_score,
    group_kkt = endpoint$point$penalty_kkt,
    converged = endpoint$point$numerically_eligible,
    sweeps = 0L, termination_reason = "analytic_penalty_limit",
    selected_groups = endpoint$point$selected_groups,
    reconstruction_error = endpoint$point$reconstruction_error,
    maximum_accepted_objective_increase = 0,
    numerically_eligible = endpoint$point$numerically_eligible,
    stringsAsFactors = FALSE
  )
  candidates <- rbind(endpoint_row, finite)
  candidates <- candidates[order(
    candidates$point_type != "penalty_limit", -candidates$lambda_ratio,
    candidates$lambda_index
  ), , drop = FALSE]
  eligible <- which(candidates$numerically_eligible %in% TRUE &
                      is.finite(candidates$validation_log_loss))
  if (!length(eligible)) {
    stop("The V9 SGLASSO repair path has no eligible candidate.",
         call. = FALSE)
  }
  selected_row <- eligible[which.min(
    candidates$validation_log_loss[eligible]
  )]
  candidates$selected <- FALSE
  candidates$selected[selected_row] <- TRUE
  selected <- candidates[selected_row, , drop = FALSE]
  invalid_loss <- candidates$validation_log_loss[
    !candidates$numerically_eligible &
      is.finite(candidates$validation_log_loss)
  ]
  invalid_competitive <- length(invalid_loss) > 0L &&
    min(invalid_loss) <= selected$validation_log_loss
  finite_selected <- selected$point_type == "finite"
  selected_interior <- !finite_selected || (
    selected$lambda_ratio > min(ratios) &&
      selected$lambda_ratio < max(ratios)
  )

  stability <- if (finite_selected) {
    li <- selected$lambda_index
    lsg_repair_selected_stability_v9(
      data, preprocess, target_original, case$alpha, case$d,
      selected$lambda, selected$lambda_ratio,
      fit$beta_solver[, li, d_index],
      fit$intercept_solver[li, d_index],
      fit$objective[li, d_index], validation_probabilities[[li]],
      lambda, configuration
    )
  } else {
    data.frame(
      selected_lambda = Inf, selected_lambda_ratio = Inf,
      forward_kkt = NA_real_, cold_kkt = NA_real_, reverse_kkt = NA_real_,
      cold_converged = TRUE, reverse_converged = TRUE,
      maximum_objective_error = 0,
      maximum_validation_loss_error = 0,
      maximum_validation_probability_error = 0,
      accepted = TRUE, stringsAsFactors = FALSE
    )
  }
  finite <- lsg_repair_context_v9(case, finite)
  candidates <- lsg_repair_context_v9(case, candidates)
  stability <- lsg_repair_context_v9(case, stability)
  summary <- lsg_repair_context_v9(case, data.frame(
    alpha = case$alpha, d = case$d,
    finite_points = nrow(finite),
    eligible_finite_points = sum(finite$numerically_eligible),
    selected_point_type = selected$point_type,
    selected_lambda_ratio = selected$lambda_ratio,
    selected_validation_log_loss = selected$validation_log_loss,
    invalid_candidate_competitive = invalid_competitive,
    selected_finite_point_interior = selected_interior,
    stability_accepted = stability$accepted,
    maximum_finite_kkt = max(finite$kkt),
    maximum_intercept_kkt = max(finite$intercept_kkt),
    maximum_reconstruction_error = max(finite$reconstruction_error),
    maximum_objective_reconstruction_error =
      max(finite$objective_reconstruction_error),
    maximum_kkt_reconstruction_error =
      max(finite$kkt_reconstruction_error),
    maximum_accepted_objective_increase =
      max(finite$maximum_accepted_objective_increase),
    stringsAsFactors = FALSE
  ))
  checks <- c(
    test_fields_absent_from_diagnostic_interface =
      identical(configuration$test_sample_policy, "forbidden") &&
      identical(names(data), c(
        "X_train", "y_train", "X_validation", "y_validation", "group"
      )),
    frozen_training_validation_fingerprints_exact =
      identical(fingerprints, source$data_sha256),
    exact_39_point_finite_grid = nrow(finite) == 39L &&
      isTRUE(all.equal(finite$lambda_ratio,
                       configuration$lambda_relative_grid,
                       tolerance = 1e-13)),
    finite_outputs = all(is.finite(c(
      finite$lambda, finite$validation_log_loss, finite$objective,
      finite$kkt, finite$intercept_kkt, finite$group_kkt,
      finite$reconstruction_error
    ))),
    frozen_objective_reconstructed =
      max(finite$objective_reconstruction_error) <= 2e-12,
    frozen_kkt_reconstructed =
      max(finite$kkt_reconstruction_error) <= 2e-12,
    formerly_zero_path_now_has_eligible_candidate =
      any(finite$numerically_eligible),
    selected_candidate_numerically_eligible =
      isTRUE(selected$numerically_eligible),
    no_competitive_invalid_candidate = !invalid_competitive,
    selected_finite_candidate_is_interior = selected_interior,
    selected_start_and_direction_stable = isTRUE(stability$accepted),
    intercept_profiled_to_declared_tolerance =
      is.finite(selected$intercept_kkt) &&
        selected$intercept_kkt <= 1e-10,
    prediction_reconstruction_passes =
      max(finite$reconstruction_error) <=
        configuration$reconstruction_tolerance,
    accepted_objective_increase_is_roundoff =
      max(finite$maximum_accepted_objective_increase) <=
        configuration$accepted_objective_increase_tolerance
  )
  list(
    schema_version = "logistic_sglasso_repair_case_v9",
    case = case, finite_path = finite, candidates = candidates,
    selected_stability = stability, summary = summary,
    data_sha256 = fingerprints, evaluation_data_used = FALSE,
    hard_checks = data.frame(
      case_id = case$case_id, task_id = case$task_id,
      check = names(checks), passed = unname(checks),
      stringsAsFactors = FALSE
    )
  )
}


lsg_run_grlasso_repair_case_v9 <- function(root, case, configuration) {
  if (!is.data.frame(case) || nrow(case) != 1L ||
      case$diagnostic_type != "grpreg_group_lasso_raw") {
    stop("One frozen V9 Group Lasso repair case is required.",
         call. = FALSE)
  }
  diagnostic <- lsg_run_grpreg_task_v8(root, case$task_id)
  points <- diagnostic$points
  summary_v8 <- diagnostic$summary[1L, , drop = FALSE]
  warnings <- unique(c(
    diagnostic$raw$fit_warnings,
    diagnostic$raw$prediction_warnings
  ))
  status <- lsg_classify_grpreg_path_v9(
    "grLasso", points$iterations,
    summary_v8$requested_path_length,
    configuration$grpreg_max_iterations, warnings
  )
  selection <- lsg_select_grpreg_prefix_v9(
    points$validation_log_loss, points$finite_coefficient,
    points$finite_probability, status
  )
  points$v9_point_converged <- status$point_converged
  points$v9_lower_boundary_excluded <-
    status$returned_lower_boundary_excluded
  points$v9_numerically_eligible <- selection$numerically_eligible
  points$v9_selected <- seq_len(nrow(points)) == selection$selected_index
  points <- lsg_repair_context_v9(case, points)
  summary <- lsg_repair_context_v9(case, data.frame(
    requested_path_length = summary_v8$requested_path_length,
    returned_path_length = summary_v8$returned_path_length,
    total_iterations = status$total_iterations,
    total_iteration_limit_reached =
      status$total_iteration_limit_reached,
    group_lasso_safe_prefix = status$group_lasso_safe_prefix,
    usable_prefix_length = status$usable_prefix_length,
    selected_lambda_index = selection$selected_index,
    selected_validation_log_loss =
      selection$selected_validation_log_loss,
    selected_strictly_above_usable_lower_boundary =
      selection$selected_strictly_above_usable_lower_boundary,
    unexpected_warning_count = length(status$unexpected_warnings),
    stringsAsFactors = FALSE
  ))
  checks <- c(
    test_fields_absent_from_diagnostic_interface =
      identical(configuration$test_sample_policy, "forbidden") &&
      !any(grepl(
        "test", c(
          names(diagnostic$points), names(diagnostic$summary),
          names(diagnostic$data_summary)
        ), ignore.case = TRUE
      )),
    raw_path_schema_compatible = isTRUE(summary_v8$raw_schema_compatible),
    raw_fit_and_prediction_succeeded =
      !nzchar(summary_v8$fit_error) && !nzchar(summary_v8$prediction_error),
    raw_returned_points_all_finite = nrow(points) > 1L &&
      all(points$finite_coefficient) && all(points$finite_probability) &&
      all(points$finite_validation_loss),
    cumulative_iteration_budget_exact =
      isTRUE(status$total_iteration_limit_reached) &&
      status$total_iterations == configuration$grpreg_max_iterations,
    safe_prefix_policy_recognized =
      isTRUE(status$group_lasso_safe_prefix) &&
      isTRUE(status$path_termination_acceptable),
    budget_ending_point_excluded =
      tail(status$returned_lower_boundary_excluded, 1L) &&
      !tail(status$point_converged, 1L),
    usable_prefix_available = status$usable_prefix_length > 1L,
    selected_prefix_candidate_eligible =
      isTRUE(selection$numerically_eligible[selection$selected_index]),
    selected_strictly_above_usable_lower_boundary =
      isTRUE(selection$selected_strictly_above_usable_lower_boundary),
    no_unexpected_solver_warning = !length(status$unexpected_warnings)
  )
  list(
    schema_version = "logistic_group_lasso_safe_prefix_case_v9",
    case = case, points = points, summary = summary,
    data_sha256 = diagnostic$data_sha256,
    controls = diagnostic$controls,
    evaluation_data_used = FALSE,
    hard_checks = data.frame(
      case_id = case$case_id, task_id = case$task_id,
      check = names(checks), passed = unname(checks),
      stringsAsFactors = FALSE
    )
  )
}


lsg_dispatch_repair_case_v9 <- function(root, audit, case, configuration) {
  scenario <- audit$specification$design[
    audit$specification$design$scenario_index == case$scenario_index,
    , drop = FALSE
  ]
  if (nrow(scenario) != 1L) {
    stop("Missing frozen scenario for V9 repair case.", call. = FALSE)
  }
  if (case$diagnostic_type %in%
      c("sglasso_intercept", "sglasso_lambda_tail")) {
    source <- audit$blinded_sources[[as.character(case$task_id)]]
    return(lsg_run_sglasso_repair_case_v9(
      case, scenario, source, audit$specification$configuration,
      configuration
    ))
  }
  lsg_run_grlasso_repair_case_v9(root, case, configuration)
}
