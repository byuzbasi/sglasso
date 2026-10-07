# Validation-tuned penalized logistic benchmarks for the rho_b study.
# This file is versioned so the completed pre-benchmark study remains frozen.

validate_penalized_benchmark_grid_v1 <- function(
    X,
    y,
    group,
    X_validation,
    y_validation,
    alpha_grid,
    nlambda,
    lambda_min_ratio
) {
  X <- as.matrix(X)
  X_validation <- as.matrix(X_validation)
  storage.mode(X) <- "double"
  storage.mode(X_validation) <- "double"
  y <- normalize_binary_response(y)$y
  y_validation <- normalize_binary_response(y_validation)$y
  alpha_grid <- sort(unique(as.numeric(alpha_grid)))
  nlambda <- as.integer(nlambda)

  if (nrow(X) != length(y) || nrow(X_validation) != length(y_validation) ||
      ncol(X) != ncol(X_validation)) {
    stop("Training and validation dimensions are incompatible.", call. = FALSE)
  }
  if (length(group) != ncol(X) || anyNA(group)) {
    stop("group must contain one non-missing label per predictor.", call. = FALSE)
  }
  if (anyNA(X) || anyNA(X_validation) ||
      any(!is.finite(X)) || any(!is.finite(X_validation))) {
    stop("Benchmark design matrices must be finite and non-missing.",
         call. = FALSE)
  }
  if (!length(alpha_grid) || any(!is.finite(alpha_grid)) ||
      any(alpha_grid < 0 | alpha_grid > 1)) {
    stop("benchmark alpha_grid must lie in [0, 1].", call. = FALSE)
  }
  if (nlambda < 2L || !is.finite(lambda_min_ratio) ||
      lambda_min_ratio <= 0 || lambda_min_ratio >= 1) {
    stop("Invalid benchmark lambda grid controls.", call. = FALSE)
  }

  list(
    X = X,
    y = y,
    group = group,
    X_validation = X_validation,
    y_validation = y_validation,
    alpha_grid = alpha_grid,
    nlambda = nlambda,
    lambda_min_ratio = lambda_min_ratio
  )
}


validation_log_loss_path_v1 <- function(y, probability) {
  probability <- as.matrix(probability)
  vapply(
    seq_len(ncol(probability)),
    function(index) binary_log_loss(y, probability[, index]),
    numeric(1)
  )
}


fit_grpnet_validation_grid_v1 <- function(
    X,
    y,
    group,
    X_validation,
    y_validation,
    alpha_grid,
    nlambda,
    lambda_min_ratio,
    tolerance = 1e-6,
    max_iterations = 1000000L
) {
  if (!requireNamespace("grpnet", quietly = TRUE)) {
    stop("The installed grpnet package is required.", call. = FALSE)
  }
  input <- validate_penalized_benchmark_grid_v1(
    X, y, group, X_validation, y_validation, alpha_grid,
    nlambda, lambda_min_ratio
  )
  if (!any(abs(input$alpha_grid - 1) <= 1e-12)) {
    stop("grpnet alpha_grid must include alpha = 1 for group lasso.",
         call. = FALSE)
  }
  max_iterations <- as.integer(max_iterations)
  group_factor <- factor(input$group)
  penalty_factor <- sqrt(as.numeric(table(group_factor)))
  fits <- vector("list", length(input$alpha_grid))
  elapsed <- numeric(length(input$alpha_grid))
  tuning <- vector("list", length(input$alpha_grid))

  for (ai in seq_along(input$alpha_grid)) {
    alpha <- input$alpha_grid[ai]
    elapsed[ai] <- system.time({
      fit <- grpnet::grpnet(
        input$X,
        input$y,
        input$group,
        family = "binomial",
        alpha = alpha,
        nlambda = input$nlambda,
        lambda.min.ratio = input$lambda_min_ratio,
        penalty.factor = penalty_factor,
        penalty = "LASSO",
        standardized = FALSE,
        orthogonalized = TRUE,
        intercept = TRUE,
        thresh = tolerance,
        maxit = max_iterations,
        proglang = "Fortran",
        keep.data = FALSE
      )
      probability <- as.matrix(stats::predict(
        fit,
        newx = input$X_validation,
        s = fit$lambda,
        type = "response"
      ))
    })[["elapsed"]]
    if (ncol(probability) != length(fit$lambda) ||
        any(!is.finite(probability)) || any(!is.finite(fit$beta)) ||
        any(!is.finite(fit$a0))) {
      stop("grpnet returned an invalid prediction or coefficient path.",
           call. = FALSE)
    }
    passes <- as.numeric(fit$npasses)
    if (length(passes) == 1L) passes <- rep(passes, length(fit$lambda))
    if (length(passes) != length(fit$lambda)) {
      stop("grpnet returned incompatible iteration diagnostics.",
           call. = FALSE)
    }
    loss <- validation_log_loss_path_v1(input$y_validation, probability)
    fits[[ai]] <- fit
    tuning[[ai]] <- data.frame(
      engine = "grpnet",
      alpha_index = ai,
      alpha = alpha,
      lambda_index = seq_along(fit$lambda),
      lambda = fit$lambda,
      lambda_fraction = fit$lambda / fit$lambda[1L],
      validation_log_loss = loss,
      passes = passes,
      converged = passes < max_iterations,
      stringsAsFactors = FALSE
    )
  }
  tuning <- do.call(rbind, tuning)
  rownames(tuning) <- NULL
  selected_group_en <- which(
    tuning$validation_log_loss == min(tuning$validation_log_loss)
  )[1L]
  lasso_rows <- which(abs(tuning$alpha - 1) <= 1e-12)
  selected_group_lasso <- lasso_rows[which.min(
    tuning$validation_log_loss[lasso_rows]
  )]
  tuning$selected_group_elastic_net <- FALSE
  tuning$selected_group_lasso <- FALSE
  tuning$selected_group_elastic_net[selected_group_en] <- TRUE
  tuning$selected_group_lasso[selected_group_lasso] <- TRUE

  list(
    fits = fits,
    alpha_grid = input$alpha_grid,
    tuning = tuning,
    elapsed_by_alpha = elapsed,
    runtime_seconds = sum(elapsed),
    selected_group_elastic_net = tuning[selected_group_en, , drop = FALSE],
    selected_group_lasso = tuning[selected_group_lasso, , drop = FALSE],
    tolerance = tolerance,
    max_iterations = max_iterations
  )
}


fit_glmnet_validation_grid_v1 <- function(
    X,
    y,
    group,
    X_validation,
    y_validation,
    alpha_grid,
    nlambda,
    lambda_min_ratio,
    tolerance = 1e-7,
    max_iterations = 1000000L
) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("The installed glmnet package is required.", call. = FALSE)
  }
  input <- validate_penalized_benchmark_grid_v1(
    X, y, group, X_validation, y_validation, alpha_grid,
    nlambda, lambda_min_ratio
  )
  max_iterations <- as.integer(max_iterations)
  fits <- vector("list", length(input$alpha_grid))
  elapsed <- numeric(length(input$alpha_grid))
  tuning <- vector("list", length(input$alpha_grid))

  for (ai in seq_along(input$alpha_grid)) {
    alpha <- input$alpha_grid[ai]
    elapsed[ai] <- system.time({
      fit <- glmnet::glmnet(
        x = input$X,
        y = input$y,
        family = "binomial",
        alpha = alpha,
        nlambda = input$nlambda,
        lambda.min.ratio = input$lambda_min_ratio,
        standardize = TRUE,
        intercept = TRUE,
        control = list(
          thresh = tolerance,
          maxit = max_iterations
        )
      )
      probability <- as.matrix(stats::predict(
        fit,
        newx = input$X_validation,
        s = fit$lambda,
        type = "response"
      ))
    })[["elapsed"]]
    if (fit$jerr != 0L || ncol(probability) != length(fit$lambda) ||
        any(!is.finite(probability))) {
      stop("glmnet returned an invalid or incomplete path.", call. = FALSE)
    }
    loss <- validation_log_loss_path_v1(input$y_validation, probability)
    fits[[ai]] <- fit
    tuning[[ai]] <- data.frame(
      engine = "glmnet",
      alpha_index = ai,
      alpha = alpha,
      lambda_index = seq_along(fit$lambda),
      lambda = fit$lambda,
      lambda_fraction = fit$lambda / fit$lambda[1L],
      validation_log_loss = loss,
      passes = as.numeric(fit$npasses),
      converged = fit$jerr == 0L,
      stringsAsFactors = FALSE
    )
  }
  tuning <- do.call(rbind, tuning)
  rownames(tuning) <- NULL
  selected <- which(
    tuning$validation_log_loss == min(tuning$validation_log_loss)
  )[1L]
  tuning$selected_logistic_elastic_net <- FALSE
  tuning$selected_logistic_elastic_net[selected] <- TRUE

  list(
    fits = fits,
    alpha_grid = input$alpha_grid,
    tuning = tuning,
    elapsed_by_alpha = elapsed,
    runtime_seconds = sum(elapsed),
    selected = tuning[selected, , drop = FALSE],
    tolerance = tolerance,
    max_iterations = max_iterations
  )
}


extract_grpnet_benchmark_solution_v1 <- function(grid, selection, newx) {
  selection <- selection[1L, , drop = FALSE]
  fit <- grid$fits[[selection$alpha_index]]
  li <- selection$lambda_index
  coefficient <- c(fit$a0[li], fit$beta[, li])
  probability <- drop(stats::predict(
    fit,
    newx = as.matrix(newx),
    s = fit$lambda[li],
    type = "response"
  ))
  list(
    coefficient = as.numeric(coefficient),
    probability = probability,
    alpha = selection$alpha,
    lambda = selection$lambda,
    lambda_fraction = selection$lambda_fraction,
    lambda_index = li,
    path_length = length(fit$lambda),
    validation_log_loss = selection$validation_log_loss,
    converged = isTRUE(selection$converged)
  )
}


extract_glmnet_benchmark_solution_v1 <- function(grid, selection, newx) {
  selection <- selection[1L, , drop = FALSE]
  fit <- grid$fits[[selection$alpha_index]]
  coefficient <- as.numeric(as.matrix(stats::coef(
    fit,
    s = selection$lambda
  ))[, 1L])
  probability <- drop(stats::predict(
    fit,
    newx = as.matrix(newx),
    s = selection$lambda,
    type = "response"
  ))
  list(
    coefficient = coefficient,
    probability = probability,
    alpha = selection$alpha,
    lambda = selection$lambda,
    lambda_fraction = selection$lambda_fraction,
    lambda_index = selection$lambda_index,
    path_length = length(fit$lambda),
    validation_log_loss = selection$validation_log_loss,
    converged = isTRUE(selection$converged)
  )
}


evaluate_penalized_benchmark_v1 <- function(
    solution,
    data,
    scenario,
    replication,
    method,
    engine,
    runtime_seconds,
    runtime_scope
) {
  if (length(solution$coefficient) != ncol(data$X_test) + 1L ||
      any(!is.finite(solution$coefficient)) ||
      any(!is.finite(solution$probability))) {
    stop("Benchmark coefficients or probabilities are invalid.",
         call. = FALSE)
  }
  reconstructed <- stats::plogis(
    solution$coefficient[1L] +
      drop(data$X_test %*% solution$coefficient[-1L])
  )
  reconstruction_error <- max(abs(reconstructed - solution$probability))
  if (!is.finite(reconstruction_error) || reconstruction_error > 2e-6) {
    stop("Benchmark coefficient/prediction consistency failed for ",
         method, ".", call. = FALSE)
  }
  selected_groups <- selected_groups_from_coefficients(
    solution$coefficient[-1L], data$group
  )
  metrics <- evaluate_study_method(
    data$y_test,
    solution$probability,
    selected_groups,
    data$active_groups,
    solution$coefficient,
    c(data$intercept, data$beta),
    data$X_test
  )
  cbind(
    data.frame(
      scenario = scenario$scenario[[1L]],
      replication = replication,
      seed = data$seed,
      method = method,
      selected_alpha = solution$alpha,
      selected_d = NA_real_,
      selected_lambda_fraction = solution$lambda_fraction,
      validation_log_loss = solution$validation_log_loss,
      selected_kkt = NA_real_,
      selected_converged = solution$converged,
      selected_lambda_on_boundary = solution$lambda_index %in%
        c(1L, solution$path_length),
      selected_d_on_boundary = NA,
      grid_runtime_seconds = runtime_seconds,
      stringsAsFactors = FALSE
    ),
    metrics,
    data.frame(
      selected_lambda = solution$lambda,
      tuning_engine = engine,
      runtime_scope = runtime_scope,
      prediction_reconstruction_error = reconstruction_error,
      stringsAsFactors = FALSE
    )
  )
}
