# Optimal Logistic SGLASSO plus penalized-logistic benchmark study utilities.
# Existing rho_b study utilities remain unchanged for source reproducibility.

fit_optimal_sglasso_validation_grid_v1 <- function(
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
  if (!length(alpha_grid) || any(alpha_grid <= 0 | alpha_grid > 1)) {
    stop("alpha_grid must lie in (0, 1].", call. = FALSE)
  }
  if (!length(d_grid) || any(d_grid < 0 | d_grid > 1) ||
      !any(abs(d_grid) <= 1e-12)) {
    stop("d_grid must lie in [0, 1] and contain zero.", call. = FALSE)
  }

  fits <- vector("list", length(alpha_grid))
  loss <- vector("list", length(alpha_grid))
  elapsed <- numeric(length(alpha_grid))
  preprocess <- prepare_lsg_design(X, group)
  for (ai in seq_along(alpha_grid)) {
    alpha <- alpha_grid[ai]
    # At alpha = 1 the shifted quadratic term vanishes, so all d values
    # describe the identical group-lasso objective. Fit that path once.
    fitted_d <- if (abs(alpha - 1) <= 1e-12) 0 else d_grid
    elapsed[ai] <- system.time({
      fits[[ai]] <- fit_logistic_sglasso(
        X,
        y,
        group,
        nlambda = configuration$nlambda,
        lambda_min_ratio = configuration$lambda_min_ratio,
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
  }

  list(
    fits = fits,
    validation_loss = loss,
    alpha_grid = alpha_grid,
    d_grid = d_grid,
    elapsed_by_alpha = elapsed,
    elapsed_seconds = sum(elapsed)
  )
}


flatten_sglasso_validation_grid_v1 <- function(
    grid,
    allowed_d,
    kkt_limit
) {
  allowed_d <- sort(unique(as.numeric(allowed_d)))
  kkt_limit <- as.numeric(kkt_limit)
  if (!length(allowed_d) || any(!is.finite(allowed_d))) {
    stop("allowed_d must contain finite values.", call. = FALSE)
  }
  if (length(kkt_limit) != 1L || !is.finite(kkt_limit) ||
      kkt_limit <= 0) {
    stop("kkt_limit must be one finite positive value.", call. = FALSE)
  }

  frames <- vector("list", length(grid$fits))
  for (ai in seq_along(grid$fits)) {
    fit <- grid$fits[[ai]]
    index <- expand.grid(
      lambda_index = seq_along(fit$lambda),
      d_index = seq_along(fit$d),
      KEEP.OUT.ATTRS = FALSE
    )
    candidate_d <- fit$d[index$d_index]
    d_allowed <- vapply(
      candidate_d,
      function(value) any(abs(value - allowed_d) <= 1e-12),
      logical(1)
    )
    validation_log_loss <- as.vector(grid$validation_loss[[ai]])
    kkt <- as.vector(fit$kkt)
    converged <- as.vector(fit$converged)
    finite_validation_loss <- is.finite(validation_log_loss)
    numerically_eligible <- d_allowed & finite_validation_loss &
      is.finite(kkt) & converged & kkt <= kkt_limit
    frames[[ai]] <- data.frame(
      alpha_index = ai,
      alpha = fit$alpha,
      lambda_index = index$lambda_index,
      d_index = index$d_index,
      lambda = fit$lambda[index$lambda_index],
      lambda_fraction = fit$lambda[index$lambda_index] / fit$lambda[1L],
      d = candidate_d,
      validation_log_loss = validation_log_loss,
      objective = as.vector(fit$objective),
      kkt = kkt,
      converged = converged,
      passes = as.vector(fit$passes),
      selected_groups = as.vector(fit$selected_groups),
      d_allowed = d_allowed,
      finite_validation_loss = finite_validation_loss,
      numerically_eligible = numerically_eligible,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, frames)
  rownames(out) <- NULL
  out
}


select_numerically_valid_sglasso_grid_v1 <- function(
    grid,
    allowed_d,
    kkt_limit
) {
  candidates <- flatten_sglasso_validation_grid_v1(
    grid, allowed_d, kkt_limit
  )
  candidates <- candidates[candidates$d_allowed, , drop = FALSE]
  valid <- candidates[candidates$numerically_eligible, , drop = FALSE]
  if (!nrow(valid)) {
    stop(
      "No Logistic SGLASSO validation candidate passed the convergence/KKT ",
      "eligibility rule.",
      call. = FALSE
    )
  }

  best_valid <- valid[which.min(valid$validation_log_loss), , drop = FALSE]
  invalid_with_loss <- candidates[
    !candidates$numerically_eligible & candidates$finite_validation_loss,
    ,
    drop = FALSE
  ]
  if (nrow(invalid_with_loss)) {
    best_invalid_validation_log_loss <- min(
      invalid_with_loss$validation_log_loss
    )
    invalid_minus_valid_validation_log_loss <-
      best_invalid_validation_log_loss - best_valid$validation_log_loss
    invalid_candidate_competitive <-
      best_invalid_validation_log_loss <= best_valid$validation_log_loss
  } else {
    best_invalid_validation_log_loss <- NA_real_
    invalid_minus_valid_validation_log_loss <- NA_real_
    invalid_candidate_competitive <- FALSE
  }

  requested_d <- sort(unique(as.numeric(allowed_d)))
  valid_d <- sort(unique(valid$d))
  list(
    alpha_index = best_valid$alpha_index[[1L]],
    lambda_index = best_valid$lambda_index[[1L]],
    d_index = best_valid$d_index[[1L]],
    validation_log_loss = best_valid$validation_log_loss[[1L]],
    numerical_eligibility_kkt_limit = as.numeric(kkt_limit),
    total_candidate_points = nrow(candidates),
    numerically_eligible_points = nrow(valid),
    numerically_ineligible_points = nrow(candidates) - nrow(valid),
    nonfinite_validation_loss_points = sum(
      !candidates$finite_validation_loss
    ),
    best_invalid_validation_log_loss =
      best_invalid_validation_log_loss,
    invalid_minus_valid_validation_log_loss =
      invalid_minus_valid_validation_log_loss,
    invalid_candidate_competitive = invalid_candidate_competitive,
    valid_alpha_count = length(unique(valid$alpha_index)),
    fitted_alpha_count = length(grid$fits),
    all_alpha_have_valid_candidates =
      length(unique(valid$alpha_index)) == length(grid$fits),
    valid_d_count = length(valid_d),
    requested_d_count = length(requested_d),
    all_d_have_valid_candidates = all(vapply(
      requested_d,
      function(value) any(abs(value - valid_d) <= 1e-12),
      logical(1)
    ))
  )
}


sglasso_tuning_frame_v1 <- function(
    grid,
    selection,
    scenario,
    replication,
    seed,
    allowed_d,
    kkt_limit
) {
  out <- flatten_sglasso_validation_grid_v1(
    grid, allowed_d, kkt_limit
  )
  out$selected <- out$alpha_index == selection$alpha_index &
    out$lambda_index == selection$lambda_index &
    out$d_index == selection$d_index
  out$validation_loss_minus_selected <-
    out$validation_log_loss - selection$validation_log_loss
  out$invalid_validation_contender <- !out$numerically_eligible &
    out$finite_validation_loss &
    out$validation_log_loss <= selection$validation_log_loss
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


evaluate_optimal_sglasso_v1 <- function(
    grid,
    selection,
    data,
    scenario,
    replication
) {
  result <- evaluate_aligned_selection(
    grid,
    selection,
    data,
    scenario,
    replication,
    method = "Logistic SGLASSO"
  )
  fit <- grid$fits[[selection$alpha_index]]
  coefficient <- as.numeric(fit$coefficients[
    , selection$lambda_index, selection$d_index
  ])
  probability <- drop(predict_logistic_sglasso(
    fit,
    data$X_test,
    type = "response",
    lambda_index = selection$lambda_index,
    d_index = selection$d_index
  ))
  reconstructed <- stats::plogis(
    coefficient[1L] + drop(data$X_test %*% coefficient[-1L])
  )
  reconstruction_error <- max(abs(reconstructed - probability))
  if (!is.finite(reconstruction_error) || reconstruction_error > 2e-6) {
    stop("Logistic SGLASSO coefficient/prediction consistency failed.",
         call. = FALSE)
  }
  cbind(
    result,
    data.frame(
      selected_lambda = fit$lambda[selection$lambda_index],
      tuning_engine = "RcppArmadillo_ABGD",
      runtime_scope = "joint_alpha_d_lambda_validation_grid",
      prediction_reconstruction_error = reconstruction_error,
      stringsAsFactors = FALSE
    )
  )
}


benchmark_tuning_frame_v1 <- function(
    grpnet_grid,
    glmnet_grid,
    scenario,
    replication,
    seed
) {
  grouped <- grpnet_grid$tuning
  grouped$selected_logistic_elastic_net <- FALSE
  grouped$method_path <- "Logistic Group Elastic Net (grpnet)"
  individual <- glmnet_grid$tuning
  individual$selected_group_elastic_net <- FALSE
  individual$selected_group_lasso <- FALSE
  individual$method_path <- "Logistic Elastic Net (glmnet)"
  columns <- c(
    "engine", "method_path", "alpha_index", "alpha", "lambda_index",
    "lambda", "lambda_fraction", "validation_log_loss", "passes",
    "converged", "selected_group_elastic_net", "selected_group_lasso",
    "selected_logistic_elastic_net"
  )
  out <- rbind(grouped[, columns], individual[, columns])
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


run_sglasso_rhob_benchmark_replication_v1 <- function(
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
  grpnet_grid <- fit_grpnet_validation_grid_v1(
    data$X_train,
    data$y_train,
    data$group,
    data$X_validation,
    data$y_validation,
    alpha_grid = configuration$benchmark_alpha_grid,
    nlambda = configuration$benchmark_nlambda,
    lambda_min_ratio = configuration$benchmark_lambda_min_ratio,
    tolerance = configuration$grpnet_tolerance,
    max_iterations = configuration$benchmark_max_iterations
  )
  glmnet_grid <- fit_glmnet_validation_grid_v1(
    data$X_train,
    data$y_train,
    data$group,
    data$X_validation,
    data$y_validation,
    alpha_grid = configuration$benchmark_alpha_grid,
    nlambda = configuration$benchmark_nlambda,
    lambda_min_ratio = configuration$benchmark_lambda_min_ratio,
    tolerance = configuration$glmnet_tolerance,
    max_iterations = configuration$benchmark_max_iterations
  )

  sglasso_result <- evaluate_optimal_sglasso_v1(
    sglasso_grid, sglasso_selection, data, scenario, replication
  )
  group_en_solution <- extract_grpnet_benchmark_solution_v1(
    grpnet_grid, grpnet_grid$selected_group_elastic_net, data$X_test
  )
  group_lasso_solution <- extract_grpnet_benchmark_solution_v1(
    grpnet_grid, grpnet_grid$selected_group_lasso, data$X_test
  )
  logistic_en_solution <- extract_glmnet_benchmark_solution_v1(
    glmnet_grid, glmnet_grid$selected, data$X_test
  )
  alpha_one <- which(abs(grpnet_grid$alpha_grid - 1) <= 1e-12)

  results <- rbind(
    sglasso_result,
    evaluate_penalized_benchmark_v1(
      group_en_solution,
      data,
      scenario,
      replication,
      method = "Logistic Group Elastic Net (grpnet)",
      engine = "grpnet",
      runtime_seconds = grpnet_grid$runtime_seconds,
      runtime_scope = "joint_alpha_lambda_validation_grid"
    ),
    evaluate_penalized_benchmark_v1(
      group_lasso_solution,
      data,
      scenario,
      replication,
      method = "Logistic Group Lasso (grpnet)",
      engine = "grpnet",
      runtime_seconds = grpnet_grid$elapsed_by_alpha[alpha_one],
      runtime_scope = "alpha_one_lambda_validation_path"
    ),
    evaluate_penalized_benchmark_v1(
      logistic_en_solution,
      data,
      scenario,
      replication,
      method = "Logistic Elastic Net (glmnet)",
      engine = "glmnet",
      runtime_seconds = glmnet_grid$runtime_seconds,
      runtime_scope = "joint_alpha_lambda_validation_grid"
    )
  )
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
  benchmark_tuning <- benchmark_tuning_frame_v1(
    grpnet_grid, glmnet_grid, scenario, replication, seed
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

  diagnostics <- data.frame(
    scenario = scenario$scenario[[1L]],
    replication = replication,
    seed = seed,
    rho_between = scenario$rho_between[[1L]],
    firth_failed_groups = firth$failed_groups,
    firth_separation_groups = firth$separation_groups,
    sglasso_runtime_seconds = sglasso_grid$elapsed_seconds,
    grpnet_runtime_seconds = grpnet_grid$runtime_seconds,
    glmnet_runtime_seconds = glmnet_grid$runtime_seconds,
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
    grpnet_path_convergence_rate = mean(grpnet_grid$tuning$converged),
    glmnet_path_convergence_rate = mean(glmnet_grid$tuning$converged),
    stringsAsFactors = FALSE
  )

  if (any(!benchmark_tuning$converged)) {
    stop("A penalized-logistic benchmark path failed to converge.",
         call. = FALSE)
  }

  list(
    results = results,
    sglasso_tuning = sglasso_tuning,
    benchmark_tuning = benchmark_tuning,
    targets = targets,
    diagnostics = diagnostics,
    group_diagnostics = transform(
      firth$diagnostics,
      scenario = scenario$scenario[[1L]],
      replication = replication,
      seed = seed,
      rho_between = scenario$rho_between[[1L]]
    )
  )
}
