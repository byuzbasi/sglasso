# Logistic SGLASSO study utilities for predefined-group comparators only.

augment_sglasso_group_result_v2 <- function(result) {
  cbind(
    result,
    data.frame(
      penalty_family = "sglasso_shifted_group_elastic_net",
      selected_gamma = NA_real_,
      selection_unit = "predefined_group",
      stringsAsFactors = FALSE
    )
  )
}


external_group_tuning_frame_v2 <- function(
    adelie_grid,
    grpreg_grid,
    scenario,
    replication,
    seed
) {
  columns <- c(
    "engine", "method_path", "penalty_family", "fit_index",
    "alpha_index", "alpha", "gamma", "lambda_index", "lambda",
    "lambda_fraction", "lambda_reference", "lambda_reference_type",
    "lambda_relative_to_reference", "alpha_zero_ridge_boundary",
    "group_selection_capable", "validation_log_loss", "passes",
    "solver_warning", "unexpected_solver_warning",
    "requested_path_length", "returned_path_length", "path_complete",
    "saturated_path_truncation", "iteration_budget_truncation",
    "total_iterations", "total_iteration_limit_reached",
    "returned_lower_boundary_excluded",
    "path_termination_acceptable", "point_converged",
    "finite_validation_loss",
    "numerically_eligible", "selected",
    "validation_loss_minus_selected", "invalid_validation_contender",
    "candidate_nonfinite_validation_loss",
    "policy_excluded_nonfinite_validation_loss"
  )
  out <- rbind(
    adelie_grid$tuning[, columns, drop = FALSE],
    grpreg_grid$tuning[, columns, drop = FALSE]
  )
  out <- cbind(
    data.frame(
      scenario = scenario$scenario[[1L]],
      replication = replication,
      seed = seed,
      stringsAsFactors = FALSE
    ),
    out
  )
  rownames(out) <- NULL
  out
}


external_method_audit_row_v2 <- function(
    tuning,
    selection,
    scenario,
    replication,
    seed,
    runtime_seconds
) {
  selected <- selection$row[1L, , drop = FALSE]
  method_rows <- tuning$method_path == selected$method_path[[1L]]
  method_tuning <- tuning[method_rows, , drop = FALSE]
  warning_fit <- vapply(
    split(method_tuning$solver_warning, method_tuning$fit_index),
    function(value) any(nzchar(value[!is.na(value)])),
    logical(1)
  )
  unexpected_warning_fit <- vapply(
    split(
      method_tuning$unexpected_solver_warning,
      method_tuning$fit_index
    ),
    function(value) any(nzchar(value[!is.na(value)])),
    logical(1)
  )
  data.frame(
    scenario = scenario$scenario[[1L]],
    replication = replication,
    seed = seed,
    method = selected$method_path[[1L]],
    engine = selected$engine[[1L]],
    total_candidate_points = selection$total_candidate_points,
    numerically_eligible_points = selection$numerically_eligible_points,
    numerically_ineligible_points = selection$numerically_ineligible_points,
    nonfinite_validation_loss_points =
      selection$nonfinite_validation_loss_points,
    candidate_nonfinite_validation_loss_points =
      selection$candidate_nonfinite_validation_loss_points,
    policy_excluded_nonfinite_validation_loss_points =
      selection$policy_excluded_nonfinite_validation_loss_points,
    nonfinite_validation_audit_passes =
      selection$nonfinite_validation_audit_passes,
    best_valid_validation_log_loss = selected$validation_log_loss[[1L]],
    best_invalid_validation_log_loss =
      selection$best_invalid_validation_log_loss,
    invalid_minus_valid_validation_log_loss =
      selection$invalid_minus_valid_validation_log_loss,
    invalid_candidate_competitive =
      selection$invalid_candidate_competitive,
    valid_fit_count = selection$valid_fit_count,
    requested_fit_count = selection$requested_fit_count,
    all_fit_paths_have_valid_candidates =
      selection$all_fit_paths_have_valid_candidates,
    all_path_points_complete = all(method_tuning$path_complete),
    all_path_terminations_acceptable = all(
      method_tuning$path_termination_acceptable
    ),
    saturated_path_truncation = any(
      method_tuning$saturated_path_truncation
    ),
    iteration_budget_truncation = any(
      method_tuning$iteration_budget_truncation
    ),
    total_iterations = if (all(is.na(method_tuning$total_iterations))) {
      NA_real_
    } else {
      max(method_tuning$total_iterations, na.rm = TRUE)
    },
    total_iteration_limit_reached = any(
      method_tuning$total_iteration_limit_reached
    ),
    returned_lower_boundary_excluded = any(
      method_tuning$returned_lower_boundary_excluded
    ),
    minimum_returned_path_fraction = min(
      method_tuning$returned_path_length /
        method_tuning$requested_path_length
    ),
    path_numerical_eligibility_rate = mean(
      method_tuning$numerically_eligible
    ),
    selected_numerically_eligible =
      isTRUE(selected$numerically_eligible[[1L]]),
    selected_on_returned_lower_boundary =
      selection$selected_on_returned_lower_boundary,
    returned_lower_boundary_minus_selected_validation_log_loss =
      selection$returned_lower_boundary_minus_selected_validation_log_loss,
    truncated_path_selection_interior =
      selection$truncated_path_selection_interior,
    warning_fit_count = sum(warning_fit),
    unexpected_warning_fit_count = sum(unexpected_warning_fit),
    runtime_seconds = runtime_seconds,
    stringsAsFactors = FALSE
  )
}


run_sglasso_rhob_group_benchmark_replication_v2 <- function(
    scenario,
    replication,
    configuration
) {
  seed <- as.integer(
    configuration$seed_base +
      scenario$scenario_index[[1L]] * 100000L + replication
  )
  data <- simulate_sglasso_rhob_logistic(scenario, seed)
  firth <- estimate_groupwise_logistic_target(
    data$X_train,
    data$y_train,
    data$group,
    method = "firth",
    max_iterations = configuration$target_max_iterations,
    tolerance = configuration$target_tolerance
  )
  if (!firth$success) {
    stop(
      "Firth target failed in ", firth$failed_groups,
      " group(s); no fallback is permitted.",
      call. = FALSE
    )
  }
  fisher_started <- proc.time()[["elapsed"]]
  fisher <- estimate_null_fisher_target_original(
    data$X_train, data$y_train, data$group
  )
  fisher_seconds <- proc.time()[["elapsed"]] - fisher_started

  sglasso_grid <- fit_optimal_sglasso_validation_grid_v1(
    data$X_train,
    data$y_train,
    data$group,
    data$X_validation,
    data$y_validation,
    firth$target_original,
    configuration
  )
  sglasso_selection <- select_numerically_valid_sglasso_grid_v1(
    sglasso_grid,
    configuration$d_grid,
    configuration$full_path_kkt_limit
  )
  adelie_grid <- fit_adelie_validation_grid_v2(
    data$X_train,
    data$y_train,
    data$group,
    data$X_validation,
    data$y_validation,
    alpha_grid = configuration$benchmark_alpha_grid,
    nlambda = configuration$benchmark_nlambda,
    lambda_min_ratio = configuration$benchmark_lambda_min_ratio,
    tolerance = configuration$adelie_tolerance,
    max_iterations = configuration$adelie_max_iterations,
    irls_tolerance = configuration$adelie_irls_tolerance,
    irls_max_iterations = configuration$adelie_irls_max_iterations
  )
  grpreg_grid <- fit_grpreg_validation_paths_v2(
    data$X_train,
    data$y_train,
    data$group,
    data$X_validation,
    data$y_validation,
    nlambda = configuration$benchmark_nlambda,
    lambda_min_ratio = configuration$benchmark_lambda_min_ratio,
    tolerance = configuration$grpreg_tolerance,
    max_iterations = configuration$grpreg_max_iterations,
    penalty_specification = configuration$grpreg_penalty_specification
  )

  sglasso_result <- augment_sglasso_group_result_v2(
    evaluate_optimal_sglasso_v1(
      sglasso_grid, sglasso_selection, data, scenario, replication
    )
  )
  adelie_solution <- extract_adelie_group_en_solution_v2(
    adelie_grid, data$X_test
  )
  grpreg_solutions <- lapply(
    configuration$grpreg_penalty_specification$penalty,
    function(penalty) {
      extract_grpreg_group_solution_v2(
        grpreg_grid, penalty, data$X_test
      )
    }
  )
  names(grpreg_solutions) <-
    configuration$grpreg_penalty_specification$penalty

  external_results <- list(
    evaluate_group_benchmark_solution_v2(
      adelie_solution,
      data,
      scenario,
      replication,
      method = "Logistic Group Elastic Net (adelie)",
      engine = "adelie",
      runtime_seconds = adelie_grid$runtime_seconds,
      runtime_scope = "joint_alpha_lambda_validation_grid"
    )
  )
  for (pi in seq_len(nrow(
      configuration$grpreg_penalty_specification
  ))) {
    specification <- configuration$grpreg_penalty_specification[
      pi, , drop = FALSE
    ]
    external_results[[length(external_results) + 1L]] <-
      evaluate_group_benchmark_solution_v2(
        grpreg_solutions[[specification$penalty]],
        data,
        scenario,
        replication,
        method = specification$method_path,
        engine = "grpreg",
        runtime_seconds =
          grpreg_grid$elapsed_by_penalty[[specification$penalty]],
        runtime_scope = "single_penalty_lambda_validation_path"
      )
  }
  results <- do.call(rbind, c(list(sglasso_result), external_results))
  rownames(results) <- NULL

  sglasso_tuning <- sglasso_tuning_frame_v1(
    sglasso_grid,
    sglasso_selection,
    scenario,
    replication,
    seed,
    configuration$d_grid,
    configuration$full_path_kkt_limit
  )
  external_tuning <- external_group_tuning_frame_v2(
    adelie_grid,
    grpreg_grid,
    scenario,
    replication,
    seed
  )
  targets <- rbind(
    aligned_target_row(
      firth$target_original,
      "groupwise_firth",
      data,
      scenario,
      replication,
      firth$failed_groups,
      firth$separation_groups,
      firth$elapsed_seconds
    ),
    aligned_target_row(
      fisher$target_original,
      "fisher",
      data,
      scenario,
      replication,
      0L,
      NA_integer_,
      fisher_seconds
    )
  )
  rownames(targets) <- NULL

  external_diagnostics <- list(
    external_method_audit_row_v2(
      adelie_grid$tuning,
      adelie_grid$selection,
      scenario,
      replication,
      seed,
      adelie_grid$runtime_seconds
    )
  )
  for (pi in seq_len(nrow(
      configuration$grpreg_penalty_specification
  ))) {
    specification <- configuration$grpreg_penalty_specification[
      pi, , drop = FALSE
    ]
    external_diagnostics[[length(external_diagnostics) + 1L]] <-
      external_method_audit_row_v2(
        grpreg_grid$tuning,
        grpreg_grid$selections[[specification$penalty]],
        scenario,
        replication,
        seed,
        grpreg_grid$elapsed_by_penalty[[specification$penalty]]
      )
  }
  external_diagnostics <- do.call(rbind, external_diagnostics)
  rownames(external_diagnostics) <- NULL

  diagnostics <- data.frame(
    scenario = scenario$scenario[[1L]],
    replication = replication,
    seed = seed,
    rho_between = scenario$rho_between[[1L]],
    firth_failed_groups = firth$failed_groups,
    firth_separation_groups = firth$separation_groups,
    sglasso_runtime_seconds = sglasso_grid$elapsed_seconds,
    adelie_runtime_seconds = adelie_grid$runtime_seconds,
    grpreg_runtime_seconds = grpreg_grid$runtime_seconds,
    sglasso_path_convergence_rate = mean(sglasso_tuning$converged),
    sglasso_path_maximum_kkt = max(sglasso_tuning$kkt),
    sglasso_nonconverged_points = sum(!sglasso_tuning$converged),
    sglasso_kkt_ineligible_points = sum(
      !is.finite(sglasso_tuning$kkt) |
        sglasso_tuning$kkt > configuration$full_path_kkt_limit
    ),
    sglasso_numerically_eligible_points =
      sglasso_selection$numerically_eligible_points,
    sglasso_numerically_ineligible_points =
      sglasso_selection$numerically_ineligible_points,
    sglasso_nonfinite_validation_loss_points =
      sglasso_selection$nonfinite_validation_loss_points,
    sglasso_selected_numerically_eligible = all(
      sglasso_tuning$numerically_eligible[sglasso_tuning$selected]
    ),
    sglasso_best_valid_validation_log_loss =
      sglasso_selection$validation_log_loss,
    sglasso_best_invalid_validation_log_loss =
      sglasso_selection$best_invalid_validation_log_loss,
    sglasso_invalid_minus_valid_validation_log_loss =
      sglasso_selection$invalid_minus_valid_validation_log_loss,
    sglasso_invalid_candidate_competitive =
      sglasso_selection$invalid_candidate_competitive,
    sglasso_all_alpha_have_valid_candidates =
      sglasso_selection$all_alpha_have_valid_candidates,
    sglasso_all_d_have_valid_candidates =
      sglasso_selection$all_d_have_valid_candidates,
    sglasso_valid_alpha_count = sglasso_selection$valid_alpha_count,
    sglasso_fitted_alpha_count = sglasso_selection$fitted_alpha_count,
    sglasso_valid_d_count = sglasso_selection$valid_d_count,
    sglasso_requested_d_count = sglasso_selection$requested_d_count,
    sglasso_maximum_raw_objective_increase = max(vapply(
      sglasso_grid$fits,
      function(fit) max(fit$maximum_raw_objective_increase),
      numeric(1)
    )),
    external_selected_numerically_eligible = all(
      external_diagnostics$selected_numerically_eligible
    ),
    external_all_path_points_complete = all(
      external_diagnostics$all_path_points_complete
    ),
    external_all_path_terminations_acceptable = all(
      external_diagnostics$all_path_terminations_acceptable
    ),
    external_saturated_path_truncation_count = sum(
      external_diagnostics$saturated_path_truncation
    ),
    external_truncated_path_selections_interior = all(
      external_diagnostics$truncated_path_selection_interior
    ),
    external_unexpected_warning_fit_count = sum(
      external_diagnostics$unexpected_warning_fit_count
    ),
    external_all_fit_paths_have_valid_candidates = all(
      external_diagnostics$all_fit_paths_have_valid_candidates
    ),
    external_invalid_candidate_competitive = any(
      external_diagnostics$invalid_candidate_competitive
    ),
    external_nonfinite_validation_loss_points = sum(
      external_diagnostics$nonfinite_validation_loss_points
    ),
    external_candidate_nonfinite_validation_loss_points = sum(
      external_diagnostics$candidate_nonfinite_validation_loss_points
    ),
    external_policy_excluded_nonfinite_validation_loss_points = sum(
      external_diagnostics$policy_excluded_nonfinite_validation_loss_points
    ),
    external_nonfinite_validation_audit_passes = all(
      external_diagnostics$nonfinite_validation_audit_passes
    ),
    external_numerically_ineligible_points = sum(
      external_diagnostics$numerically_ineligible_points
    ),
    stringsAsFactors = FALSE
  )

  list(
    results = results,
    sglasso_tuning = sglasso_tuning,
    external_tuning = external_tuning,
    targets = targets,
    diagnostics = diagnostics,
    external_diagnostics = external_diagnostics,
    group_diagnostics = transform(
      firth$diagnostics,
      scenario = scenario$scenario[[1L]],
      replication = replication,
      seed = seed,
      rho_between = scenario$rho_between[[1L]]
    )
  )
}
