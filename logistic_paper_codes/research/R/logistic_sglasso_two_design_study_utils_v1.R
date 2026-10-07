# Versioned fitting utilities for the TRUBA Logistic SGLASSO study.

scenario_value_v1 <- function(scenario, name, default) {
  if (name %in% names(scenario)) scenario[[name]][[1L]] else default
}


add_two_design_metadata_v1 <- function(frame, scenario) {
  frame$design_id <- scenario_value_v1(
    scenario, "design_id", "homogeneous"
  )
  frame$signal_pattern <- scenario_value_v1(
    scenario, "signal_pattern", "homogeneous"
  )
  frame$n_train <- scenario$n_train[[1L]]
  frame$n_validation <- scenario$n_validation[[1L]]
  frame$n_test <- scenario$n_test[[1L]]
  frame$groups <- scenario$groups[[1L]]
  frame$group_size <- scenario$group_size[[1L]]
  frame$p <- scenario$groups[[1L]] * scenario$group_size[[1L]]
  frame$active_group_count <- scenario$active_group_count[[1L]]
  frame$rho_within <- scenario$rho_within[[1L]]
  frame$rho_between <- scenario$rho_between[[1L]]
  frame
}


fit_extended_sglasso_validation_grid_v1 <- function(
    X,
    y,
    group,
    X_validation,
    y_validation,
    target_original,
    configuration
) {
  alpha_grid <- sort(unique(as.numeric(configuration$alpha_grid)))
  d_grid <- sort(unique(as.numeric(configuration$d_grid)))
  if (!length(alpha_grid) || any(alpha_grid < 0 | alpha_grid > 1)) {
    stop("alpha_grid must lie in [0, 1].", call. = FALSE)
  }
  if (!length(d_grid) || any(d_grid < 0 | d_grid > 1) ||
      !any(abs(d_grid) <= 1e-12)) {
    stop("d_grid must lie in [0, 1] and contain zero.", call. = FALSE)
  }
  nlambda <- as.integer(configuration$nlambda)
  lower_fraction <- as.numeric(configuration$lambda_min_ratio)
  upper_multiplier <- as.numeric(configuration$lambda_upper_multiplier)
  if (nlambda < 2L || !is.finite(lower_fraction) ||
      lower_fraction <= 0 || lower_fraction >= 1 ||
      !is.finite(upper_multiplier) || upper_multiplier <= 1) {
    stop("Invalid extended SGLASSO lambda specification.", call. = FALSE)
  }

  response <- normalize_binary_response(y)$y
  preprocess <- prepare_lsg_design(X, group)
  target <- project_lsg_original_target(preprocess, target_original)
  fits <- vector("list", length(alpha_grid))
  loss <- vector("list", length(alpha_grid))
  elapsed <- numeric(length(alpha_grid))
  lambda_reference <- numeric(length(alpha_grid))
  lambda_d0_kkt_reference <- rep(NA_real_, length(alpha_grid))
  lambda_reference_type <- character(length(alpha_grid))
  null_score_scale <- lsg_null_score_lambda_scale(preprocess, response)

  for (ai in seq_along(alpha_grid)) {
    alpha <- alpha_grid[ai]
    fitted_d <- if (abs(alpha - 1) <= 1e-12) 0 else d_grid
    if (abs(alpha) <= 1e-12) {
      lambda_reference[ai] <- null_score_scale
      lambda_reference_type[ai] <- "null_score_ridge_boundary"
    } else {
      diagnostic <- lsg_lambda_start_cpp(
        preprocess$X,
        response,
        preprocess$group_start,
        preprocess$group_end,
        preprocess$group_weight,
        target,
        alpha,
        0
      )
      if (!isTRUE(diagnostic$zero_model_feasible) ||
          !is.finite(diagnostic$lambda_start) ||
          diagnostic$lambda_start <= 0) {
        stop("Unable to construct the d=0 KKT lambda reference.",
             call. = FALSE)
      }
      lambda_d0_kkt_reference[ai] <-
        diagnostic$lambda_start * (1 + 1e-8)
      lambda_reference[ai] <- lambda_d0_kkt_reference[ai]
      lambda_reference_type[ai] <- "d0_null_kkt"
    }
    lambda <- exp(seq(
      log(lambda_reference[ai] * upper_multiplier),
      log(lambda_reference[ai] * lower_fraction),
      length.out = nlambda
    ))
    elapsed[ai] <- system.time({
      fits[[ai]] <- fit_logistic_sglasso(
        X,
        y,
        group,
        lambda = lambda,
        d = fitted_d,
        alpha = alpha,
        max_outer = configuration$max_passes,
        max_inner = configuration$max_inner,
        tolerance = configuration$tolerance,
        inner_tolerance = configuration$inner_tolerance,
        target_original = target_original,
        compile = FALSE,
        preprocess = preprocess,
        use_active_set = TRUE,
        warm_start_d = TRUE,
        solver = configuration$solver
      )
      probability <- predict_logistic_sglasso(
        fits[[ai]], X_validation, type = "response"
      )
      loss[[ai]] <- aligned_validation_loss(y_validation, probability)
    })[["elapsed"]]
    fits[[ai]]$lambda_reference <- lambda_reference[ai]
    fits[[ai]]$lambda_reference_type <- lambda_reference_type[ai]
    fits[[ai]]$lambda_d0_kkt_reference <-
      lambda_d0_kkt_reference[ai]
    fits[[ai]]$lambda_relative_to_reference <-
      fits[[ai]]$lambda / lambda_reference[ai]
    fits[[ai]]$lambda_relative_to_d0_kkt <-
      fits[[ai]]$lambda / lambda_d0_kkt_reference[ai]
  }

  list(
    fits = fits,
    validation_loss = loss,
    alpha_grid = alpha_grid,
    d_grid = d_grid,
    elapsed_by_alpha = elapsed,
    elapsed_seconds = sum(elapsed),
    lambda_reference = lambda_reference,
    lambda_reference_type = lambda_reference_type,
    lambda_d0_kkt_reference = lambda_d0_kkt_reference,
    lambda_upper_multiplier = upper_multiplier,
    lambda_min_reference_fraction = lower_fraction
  )
}


sglasso_tuning_frame_extended_v1 <- function(
    grid,
    free_selection,
    d0_selection,
    scenario,
    replication,
    seed,
    allowed_d,
    kkt_limit
) {
  out <- sglasso_tuning_frame_v1(
    grid,
    free_selection,
    scenario,
    replication,
    seed,
    allowed_d,
    kkt_limit
  )
  out$selected_free_d <- out$selected
  out$selected_d0_boundary <-
    out$alpha_index == d0_selection$alpha_index &
    out$lambda_index == d0_selection$lambda_index &
    out$d_index == d0_selection$d_index
  reference <- grid$lambda_reference[out$alpha_index]
  kkt_reference <- grid$lambda_d0_kkt_reference[out$alpha_index]
  out$lambda_reference <- reference
  out$lambda_reference_type <-
    grid$lambda_reference_type[out$alpha_index]
  out$lambda_relative_to_reference <- out$lambda / reference
  out$lambda_relative_to_d0_kkt <- out$lambda / kkt_reference
  out$alpha_zero_ridge_boundary <- abs(out$alpha) <= 1e-12
  out$group_selection_capable <- !out$alpha_zero_ridge_boundary
  out$lambda_on_extended_upper_boundary <- out$lambda_index == 1L
  out$lambda_on_extended_lower_boundary <- vapply(
    out$alpha_index,
    function(index) length(grid$fits[[index]]$lambda),
    integer(1)
  ) == out$lambda_index
  add_two_design_metadata_v1(out, scenario)
}


run_logistic_sglasso_two_design_replication_v1 <- function(
    scenario,
    replication,
    configuration,
    data_generator
) {
  seed <- as.integer(
    configuration$seed_base +
      scenario$scenario_index[[1L]] * 100000L + replication
  )
  data <- data_generator(scenario, seed)
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

  sglasso_grid <- fit_extended_sglasso_validation_grid_v1(
    data$X_train,
    data$y_train,
    data$group,
    data$X_validation,
    data$y_validation,
    firth$target_original,
    configuration
  )
  free_selection <- select_numerically_valid_sglasso_grid_v1(
    sglasso_grid,
    configuration$d_grid,
    configuration$full_path_kkt_limit
  )
  d0_selection <- select_numerically_valid_sglasso_grid_v1(
    sglasso_grid,
    0,
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
    ridge_boundary_upper_multiplier =
      configuration$lambda_upper_multiplier,
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

  free_result <- augment_sglasso_group_result_v2(
    evaluate_optimal_sglasso_v1(
      sglasso_grid, free_selection, data, scenario, replication
    )
  )
  d0_result <- augment_sglasso_group_result_v2(
    evaluate_optimal_sglasso_v1(
      sglasso_grid, d0_selection, data, scenario, replication
    )
  )
  d0_result$method <- "Logistic SGLASSO (d=0 boundary)"
  d0_result$penalty_family <- "sglasso_d0_group_elastic_net"
  d0_result$grid_runtime_seconds <- 0
  d0_result$runtime_scope <-
    "shared_sglasso_grid_no_additional_model_fit"

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

  external_results <- list(evaluate_group_benchmark_solution_v2(
    adelie_solution,
    data,
    scenario,
    replication,
    method = "Logistic Group Elastic Net (adelie)",
    engine = "adelie",
    runtime_seconds = adelie_grid$runtime_seconds,
    runtime_scope = "joint_alpha_lambda_validation_grid"
  ))
  for (pi in seq_len(nrow(configuration$grpreg_penalty_specification))) {
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
  results <- do.call(rbind, c(
    list(free_result, d0_result), external_results
  ))
  rownames(results) <- NULL
  results$selected_lambda_relative_to_d0_kkt <- NA_real_
  results$selected_lambda_reference <- NA_real_
  results$selected_lambda_relative_to_reference <- NA_real_
  results$selected_lambda_reference_type <- NA_character_
  for (selection_name in c("free_selection", "d0_selection")) {
    selection <- get(selection_name)
    row_index <- if (selection_name == "free_selection") 1L else 2L
    results$selected_lambda_relative_to_d0_kkt[row_index] <-
      results$selected_lambda[row_index] /
      sglasso_grid$lambda_d0_kkt_reference[selection$alpha_index]
    results$selected_lambda_reference[row_index] <-
      sglasso_grid$lambda_reference[selection$alpha_index]
    results$selected_lambda_relative_to_reference[row_index] <-
      results$selected_lambda[row_index] /
      sglasso_grid$lambda_reference[selection$alpha_index]
    results$selected_lambda_reference_type[row_index] <-
      sglasso_grid$lambda_reference_type[selection$alpha_index]
  }
  adelie_row <- which(
    results$method == "Logistic Group Elastic Net (adelie)"
  )
  adelie_selected <- adelie_grid$selection$row[1L, , drop = FALSE]
  results$selected_lambda_reference[adelie_row] <-
    adelie_selected$lambda_reference[[1L]]
  results$selected_lambda_relative_to_reference[adelie_row] <-
    adelie_selected$lambda_relative_to_reference[[1L]]
  results$selected_lambda_reference_type[adelie_row] <-
    adelie_selected$lambda_reference_type[[1L]]
  results$selected_alpha_zero_ridge_boundary <-
    abs(results$selected_alpha) <= 1e-12
  results$group_selection_capable <-
    !results$selected_alpha_zero_ridge_boundary
  results <- add_two_design_metadata_v1(results, scenario)

  sglasso_tuning <- sglasso_tuning_frame_extended_v1(
    sglasso_grid,
    free_selection,
    d0_selection,
    scenario,
    replication,
    seed,
    configuration$d_grid,
    configuration$full_path_kkt_limit
  )
  external_tuning <- add_two_design_metadata_v1(
    external_group_tuning_frame_v2(
      adelie_grid, grpreg_grid, scenario, replication, seed
    ),
    scenario
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
  targets <- add_two_design_metadata_v1(targets, scenario)

  external_diagnostics <- list(external_method_audit_row_v2(
    adelie_grid$tuning,
    adelie_grid$selection,
    scenario,
    replication,
    seed,
    adelie_grid$runtime_seconds
  ))
  for (pi in seq_len(nrow(configuration$grpreg_penalty_specification))) {
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
  external_diagnostics <- add_two_design_metadata_v1(
    external_diagnostics, scenario
  )

  selected_free <- sglasso_tuning[sglasso_tuning$selected_free_d, ]
  selected_d0 <- sglasso_tuning[sglasso_tuning$selected_d0_boundary, ]
  diagnostics <- data.frame(
    scenario = scenario$scenario[[1L]],
    replication = replication,
    seed = seed,
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
      free_selection$numerically_eligible_points,
    sglasso_numerically_ineligible_points =
      free_selection$numerically_ineligible_points,
    sglasso_nonfinite_validation_loss_points =
      free_selection$nonfinite_validation_loss_points,
    sglasso_selected_numerically_eligible = all(
      selected_free$numerically_eligible
    ),
    sglasso_d0_selected_numerically_eligible = all(
      selected_d0$numerically_eligible
    ),
    sglasso_selected_alpha_lower_boundary =
      selected_free$alpha[[1L]] == min(configuration$alpha_grid),
    sglasso_selected_alpha_zero_ridge_boundary =
      selected_free$alpha_zero_ridge_boundary[[1L]],
    sglasso_selected_group_selection_capable =
      selected_free$group_selection_capable[[1L]],
    adelie_selected_alpha_zero_ridge_boundary =
      abs(adelie_grid$selection$row$alpha[[1L]]) <= 1e-12,
    adelie_selected_group_selection_capable =
      abs(adelie_grid$selection$row$alpha[[1L]]) > 1e-12,
    sglasso_selected_positive_d_upper_lambda_boundary =
      selected_free$d[[1L]] > 0 &&
        selected_free$lambda_on_extended_upper_boundary[[1L]],
    sglasso_selected_lower_lambda_boundary =
      selected_free$lambda_on_extended_lower_boundary[[1L]],
    sglasso_selected_d0_validation_log_loss =
      d0_selection$validation_log_loss,
    sglasso_free_minus_d0_validation_log_loss =
      free_selection$validation_log_loss - d0_selection$validation_log_loss,
    sglasso_best_valid_validation_log_loss =
      free_selection$validation_log_loss,
    sglasso_best_invalid_validation_log_loss =
      free_selection$best_invalid_validation_log_loss,
    sglasso_invalid_minus_valid_validation_log_loss =
      free_selection$invalid_minus_valid_validation_log_loss,
    sglasso_invalid_candidate_competitive =
      free_selection$invalid_candidate_competitive,
    sglasso_all_alpha_have_valid_candidates =
      free_selection$all_alpha_have_valid_candidates,
    sglasso_all_d_have_valid_candidates =
      free_selection$all_d_have_valid_candidates,
    sglasso_valid_alpha_count = free_selection$valid_alpha_count,
    sglasso_fitted_alpha_count = free_selection$fitted_alpha_count,
    sglasso_valid_d_count = free_selection$valid_d_count,
    sglasso_requested_d_count = free_selection$requested_d_count,
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
    external_iteration_budget_truncation_count = sum(
      external_diagnostics$iteration_budget_truncation
    ),
    external_returned_lower_boundary_exclusion_count = sum(
      external_diagnostics$returned_lower_boundary_excluded
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
  diagnostics <- add_two_design_metadata_v1(diagnostics, scenario)

  group_diagnostics <- transform(
    firth$diagnostics,
    scenario = scenario$scenario[[1L]],
    replication = replication,
    seed = seed,
    rho_between = scenario$rho_between[[1L]]
  )
  group_diagnostics <- add_two_design_metadata_v1(
    group_diagnostics, scenario
  )

  list(
    results = results,
    sglasso_tuning = sglasso_tuning,
    external_tuning = external_tuning,
    targets = targets,
    diagnostics = diagnostics,
    external_diagnostics = external_diagnostics,
    group_diagnostics = group_diagnostics
  )
}
