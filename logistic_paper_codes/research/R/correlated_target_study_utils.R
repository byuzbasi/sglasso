# Fitting utilities for the theory-aligned correlated logistic study.

estimate_null_fisher_target_original <- function(X, y, group) {
  response <- normalize_binary_response(y)$y
  preprocess <- prepare_lsg_design(X, group)
  information <- lsg_fisher_target_cpp(
    preprocess$X,
    response,
    preprocess$group_start,
    preprocess$group_end
  )
  recovered <- recover_lsg_coefficients(
    preprocess,
    drop(information$target),
    intercept = 0
  )
  list(
    target_original = as.numeric(recovered[-1L]),
    source = "fisher_null_score"
  )
}


aligned_validation_loss <- function(y, probability) {
  dimensions <- dim(probability)
  out <- matrix(
    NA_real_,
    nrow = dimensions[2L],
    ncol = dimensions[3L]
  )
  for (di in seq_len(dimensions[3L])) {
    for (li in seq_len(dimensions[2L])) {
      out[li, di] <- binary_log_loss(y, probability[, li, di])
    }
  }
  out
}


fit_aligned_validation_grid <- function(
    X,
    y,
    group,
    X_validation,
    y_validation,
    target_original,
    configuration
) {
  fits <- vector("list", length(configuration$alpha_grid))
  loss <- vector("list", length(configuration$alpha_grid))
  elapsed <- numeric(length(configuration$alpha_grid))
  preprocess <- prepare_lsg_design(X, group)
  for (ai in seq_along(configuration$alpha_grid)) {
    elapsed[ai] <- system.time({
      fits[[ai]] <- fit_logistic_sglasso(
        X,
        y,
        group,
        nlambda = configuration$nlambda,
        lambda_min_ratio = configuration$lambda_min_ratio,
        d = configuration$d_grid,
        alpha = configuration$alpha_grid[ai],
        max_outer = configuration$max_outer,
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
        fits[[ai]],
        X_validation,
        type = "response"
      )
      loss[[ai]] <- aligned_validation_loss(y_validation, probability)
    })[["elapsed"]]
  }
  list(
    fits = fits,
    validation_loss = loss,
    alpha_grid = configuration$alpha_grid,
    d_grid = configuration$d_grid,
    elapsed_seconds = sum(elapsed)
  )
}


select_aligned_validation_grid <- function(grid, allowed_d) {
  allowed_d <- as.numeric(allowed_d)
  candidates <- list()
  for (ai in seq_along(grid$fits)) {
    fit <- grid$fits[[ai]]
    eligible <- which(vapply(
      fit$d,
      function(value) any(abs(value - allowed_d) <= 1e-12),
      logical(1)
    ))
    if (!length(eligible)) next
    restricted <- grid$validation_loss[[ai]][, eligible, drop = FALSE]
    minimum <- which(restricted == min(restricted), arr.ind = TRUE)[1L, ]
    candidates[[length(candidates) + 1L]] <- list(
      alpha_index = ai,
      lambda_index = minimum[1L],
      d_index = eligible[minimum[2L]],
      validation_log_loss = restricted[minimum[1L], minimum[2L]]
    )
  }
  if (!length(candidates)) {
    stop("No fitted point belongs to the requested d subset.", call. = FALSE)
  }
  best <- which.min(vapply(
    candidates,
    function(item) item$validation_log_loss,
    numeric(1)
  ))
  candidates[[best]]
}


evaluate_aligned_selection <- function(
    grid,
    selection,
    data,
    scenario,
    replication,
    method
) {
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
  selected_groups <- selected_groups_from_coefficients(
    coefficient[-1L],
    data$group
  )
  metrics <- evaluate_study_method(
    data$y_test,
    probability,
    selected_groups,
    data$active_groups,
    coefficient,
    c(data$intercept, data$beta),
    data$X_test
  )
  selected_d <- fit$d[selection$d_index]
  selected_lambda_fraction <-
    fit$lambda[selection$lambda_index] / fit$lambda[1L]
  cbind(
    data.frame(
      scenario = scenario$scenario[[1L]],
      replication = replication,
      seed = data$seed,
      method = method,
      selected_alpha = fit$alpha,
      selected_d = selected_d,
      selected_lambda_fraction = selected_lambda_fraction,
      validation_log_loss = selection$validation_log_loss,
      selected_kkt = fit$kkt[
        selection$lambda_index,
        selection$d_index
      ],
      selected_converged = fit$converged[
        selection$lambda_index,
        selection$d_index
      ],
      selected_lambda_on_boundary = selection$lambda_index %in%
        c(1L, length(fit$lambda)),
      selected_d_on_boundary = selected_d %in% c(0, 1),
      grid_runtime_seconds = grid$elapsed_seconds,
      stringsAsFactors = FALSE
    ),
    metrics
  )
}


aligned_target_row <- function(
    target,
    target_mode,
    data,
    scenario,
    replication,
    failed_groups,
    separation_groups,
    target_fit_seconds
) {
  quality <- direct_logistic_target_quality(
    target,
    target_mode,
    data$beta,
    data$intercept,
    data$X_test
  )
  geometry <- target_hessian_geometry(
    target,
    data$beta,
    data$intercept,
    data$X_test
  )
  cbind(
    data.frame(
      scenario = scenario$scenario[[1L]],
      replication = replication,
      seed = data$seed,
      failed_groups = failed_groups,
      separation_groups = separation_groups,
      target_fit_seconds = target_fit_seconds,
      stringsAsFactors = FALSE
    ),
    quality,
    geometry[, setdiff(names(geometry), "target_available"), drop = FALSE]
  )
}


run_aligned_target_replication <- function(
    scenario,
    replication,
    configuration
) {
  seed <- as.integer(
    configuration$seed_base +
      scenario$scenario_index[[1L]] * 100000L + replication
  )
  data <- simulate_aligned_grouped_logistic(scenario, seed)
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
    data$X_train,
    data$y_train,
    data$group
  )
  fisher_seconds <- proc.time()[["elapsed"]] - fisher_started

  grid <- fit_aligned_validation_grid(
    data$X_train,
    data$y_train,
    data$group,
    data$X_validation,
    data$y_validation,
    firth$target_original,
    configuration
  )
  selections <- list(
    "Firth SGLASSO (d >= 0)" = select_aligned_validation_grid(
      grid,
      configuration$d_grid
    ),
    "Firth SGLASSO (d > 0)" = select_aligned_validation_grid(
      grid,
      configuration$positive_d_grid
    ),
    "Group elastic net (d = 0)" = select_aligned_validation_grid(grid, 0)
  )
  results <- do.call(rbind, lapply(names(selections), function(method) {
    evaluate_aligned_selection(
      grid,
      selections[[method]],
      data,
      scenario,
      replication,
      method
    )
  }))
  rownames(results) <- NULL

  fixed_d <- do.call(rbind, lapply(configuration$d_grid, function(d_value) {
    selection <- select_aligned_validation_grid(grid, d_value)
    evaluated <- evaluate_aligned_selection(
      grid,
      selection,
      data,
      scenario,
      replication,
      paste0("Firth fixed d = ", d_value)
    )
    evaluated$fixed_d <- d_value
    evaluated
  }))
  rownames(fixed_d) <- NULL

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
    firth_failed_groups = firth$failed_groups,
    firth_separation_groups = firth$separation_groups,
    grid_runtime_seconds = grid$elapsed_seconds,
    path_convergence_rate = mean(vapply(
      grid$fits,
      function(fit) mean(fit$converged),
      numeric(1)
    )),
    path_maximum_kkt = max(vapply(
      grid$fits,
      function(fit) max(fit$kkt),
      numeric(1)
    )),
    maximum_raw_objective_increase = max(vapply(
      grid$fits,
      function(fit) max(fit$maximum_raw_objective_increase),
      numeric(1)
    )),
    stringsAsFactors = FALSE
  )

  list(
    results = results,
    fixed_d = fixed_d,
    targets = targets,
    diagnostics = diagnostics,
    group_diagnostics = transform(
      firth$diagnostics,
      scenario = scenario$scenario[[1L]],
      replication = replication,
      seed = seed
    )
  )
}
