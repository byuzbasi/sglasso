# Versioned lambda-range extension for the Logistic SGLASSO v5 audit.
#
# This overlay is sourced after the frozen v4 fitting utilities. It changes
# only the SGLASSO lambda path and the alpha = 0 Adelie path. The original
# 30 ratios from 8 down to 0.05 remain exact members of the extended grid.

logistic_sglasso_lambda_relative_grid_v5 <- function() {
  grid <- getOption("sglasso.logistic_prework.lambda_relative_grid_v5", NULL)
  grid <- as.numeric(grid)
  if (!length(grid) || any(!is.finite(grid)) || any(grid <= 0) ||
      is.unsorted(-grid, strictly = TRUE)) {
    stop("The v5 relative lambda grid is missing or invalid.", call. = FALSE)
  }
  expected_base <- exp(seq(log(8), log(0.05), length.out = 30L))
  if (length(grid) != 33L ||
      !isTRUE(all.equal(grid[seq_len(3L)], c(64, 32, 16),
                        tolerance = 1e-13)) ||
      !isTRUE(all.equal(grid[4:33], expected_base, tolerance = 1e-13))) {
    stop(
      "The v5 grid must be {64,32,16} plus the frozen 30-point v4 grid.",
      call. = FALSE
    )
  }
  grid
}


activate_logistic_sglasso_lambda_extension_v5 <- function(relative_grid) {
  options(
    sglasso.logistic_prework.lambda_relative_grid_v5 =
      as.numeric(relative_grid)
  )
  invisible(logistic_sglasso_lambda_relative_grid_v5())
}


fit_extended_sglasso_validation_grid_v5 <- function(
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
  relative_grid <- logistic_sglasso_lambda_relative_grid_v5()
  if (!isTRUE(all.equal(
        relative_grid,
        as.numeric(configuration$lambda_relative_grid),
        tolerance = 1e-13
      ))) {
    stop("The active v5 lambda grid differs from configuration.",
         call. = FALSE)
  }
  if (!length(alpha_grid) || any(alpha_grid < 0 | alpha_grid > 1)) {
    stop("alpha_grid must lie in [0, 1].", call. = FALSE)
  }
  if (!length(d_grid) || any(d_grid < 0 | d_grid > 1) ||
      !any(abs(d_grid) <= 1e-12)) {
    stop("d_grid must lie in [0, 1] and contain zero.", call. = FALSE)
  }
  if (length(relative_grid) != as.integer(configuration$nlambda) ||
      abs(relative_grid[1L] - configuration$lambda_upper_multiplier) >
        1e-12 ||
      abs(tail(relative_grid, 1L) - configuration$lambda_min_ratio) >
        1e-12) {
    stop("The v5 lambda metadata does not match the fitted grid.",
         call. = FALSE)
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
    lambda <- lambda_reference[ai] * relative_grid
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
    lambda_upper_multiplier = relative_grid[1L],
    lambda_min_reference_fraction = tail(relative_grid, 1L),
    lambda_relative_grid = relative_grid
  )
}


fit_adelie_validation_grid_v5 <- function(
    X,
    y,
    group,
    X_validation,
    y_validation,
    alpha_grid,
    nlambda,
    lambda_min_ratio,
    ridge_boundary_upper_multiplier = 64,
    tolerance = 1e-7,
    max_iterations = 100000L,
    irls_tolerance = 1e-7,
    irls_max_iterations = 10000L
) {
  if (!requireNamespace("adelie", quietly = TRUE)) {
    stop("The installed adelie package is required.", call. = FALSE)
  }
  input <- validate_penalized_benchmark_grid_v1(
    X, y, group, X_validation, y_validation, alpha_grid,
    nlambda, lambda_min_ratio
  )
  relative_grid <- logistic_sglasso_lambda_relative_grid_v5()
  max_iterations <- as.integer(max_iterations)
  irls_max_iterations <- as.integer(irls_max_iterations)
  if (max_iterations < 1L || irls_max_iterations < 1L ||
      !is.finite(ridge_boundary_upper_multiplier) ||
      ridge_boundary_upper_multiplier <= 1 ||
      abs(relative_grid[1L] - ridge_boundary_upper_multiplier) > 1e-12 ||
      abs(tail(relative_grid, 1L) - input$lambda_min_ratio) > 1e-12 ||
      !is.finite(tolerance) || tolerance <= 0 ||
      !is.finite(irls_tolerance) || irls_tolerance <= 0) {
    stop("Invalid Adelie v5 numerical or lambda-grid controls.",
         call. = FALSE)
  }

  preprocess <- prepare_lsg_design(input$X, input$group)
  X_validation_solver <- transform_lsg_newx(
    preprocess, input$X_validation
  )
  group_starts <- as.integer(preprocess$group_start + 1L)
  null_score_scale <- lsg_null_score_lambda_scale(preprocess, input$y)
  coefficient_method <- getS3method(
    "coef", "grpnet", envir = asNamespace("adelie")
  )
  prediction_method <- getS3method(
    "predict", "grpnet", envir = asNamespace("adelie")
  )

  fits <- vector("list", length(input$alpha_grid))
  coefficient_paths <- vector("list", length(input$alpha_grid))
  elapsed <- numeric(length(input$alpha_grid))
  warning_messages <- vector("list", length(input$alpha_grid))
  tuning <- vector("list", length(input$alpha_grid))

  for (ai in seq_along(input$alpha_grid)) {
    alpha <- input$alpha_grid[ai]
    alpha_zero <- abs(alpha) <= 1e-12
    explicit_lambda <- if (alpha_zero) {
      null_score_scale * relative_grid
    } else {
      NULL
    }
    requested_path_length <- if (alpha_zero) {
      length(relative_grid)
    } else {
      input$nlambda
    }
    captured_warnings <- character(0)
    elapsed[ai] <- system.time({
      fit <- withCallingHandlers(
        do.call(adelie::grpnet, list(
          X = preprocess$X,
          glm = adelie::glm.binomial(input$y),
          groups = group_starts,
          alpha = alpha,
          penalty = preprocess$group_weight,
          standardize = FALSE,
          intercept = TRUE,
          lmda_path_size = input$nlambda,
          min_ratio = input$lambda_min_ratio,
          tol = tolerance,
          max_iters = max_iterations,
          irls_tol = irls_tolerance,
          irls_max_iters = irls_max_iterations,
          screen_rule = "strong",
          early_exit = FALSE,
          check_state = TRUE,
          progress_bar = FALSE,
          n_threads = 1L,
          lambda = explicit_lambda
        )),
        warning = function(condition) {
          captured_warnings <<- c(
            captured_warnings, conditionMessage(condition)
          )
          invokeRestart("muffleWarning")
        }
      )
      coefficient_path <- coefficient_method(fit)
      lambda <- as.numeric(coefficient_path$lambda)
      probability <- as.matrix(prediction_method(
        fit,
        newx = X_validation_solver,
        lambda = lambda,
        type = "response",
        n_threads = 1L
      ))
    })[["elapsed"]]

    beta <- as.matrix(coefficient_path$betas)
    intercept <- as.numeric(coefficient_path$intercepts[, 1L])
    path_length <- length(lambda)
    if (nrow(beta) != path_length || ncol(beta) != ncol(preprocess$X) ||
        length(intercept) != path_length ||
        nrow(probability) != length(input$y_validation) ||
        ncol(probability) != path_length) {
      stop("Adelie v5 returned incompatible coefficient arrays.",
           call. = FALSE)
    }
    state_lambda <- as.numeric(fit$state$lmdas)
    lambda_aligned <- length(state_lambda) == path_length &&
      (path_length == 0L || max(abs(state_lambda - lambda)) <=
         1e-10 * max(1, max(abs(lambda))))
    path_complete <- path_length == requested_path_length && lambda_aligned
    finite_coefficient <- if (path_length) {
      apply(cbind(intercept, beta), 1L, function(value) {
        all(is.finite(value))
      })
    } else {
      logical(0)
    }
    finite_probability <- if (path_length) {
      apply(probability, 2L, function(value) all(is.finite(value)))
    } else {
      logical(0)
    }
    validation_loss <- if (path_length) {
      vapply(seq_len(path_length), function(index) {
        if (finite_probability[index]) {
          binary_log_loss(input$y_validation, probability[, index])
        } else {
          NA_real_
        }
      }, numeric(1))
    } else {
      numeric(0)
    }
    finite_validation_loss <- is.finite(validation_loss)
    warning_free <- length(captured_warnings) == 0L
    path_termination_acceptable <- path_complete & warning_free
    numerically_eligible <- path_termination_acceptable &
      finite_coefficient & finite_probability & finite_validation_loss
    lambda_reference <- if (alpha_zero) null_score_scale else lambda[1L]
    lambda_reference_type <- if (alpha_zero) {
      "null_score_ridge_boundary"
    } else {
      "native_adelie_path_start"
    }

    fits[[ai]] <- fit
    coefficient_paths[[ai]] <- list(
      beta = beta,
      intercept = intercept,
      lambda = lambda
    )
    warning_messages[[ai]] <- unique(captured_warnings)
    tuning[[ai]] <- data.frame(
      engine = "adelie",
      method_path = "Logistic Group Elastic Net (adelie)",
      penalty_family = "group_elastic_net",
      fit_index = ai,
      alpha_index = ai,
      alpha = alpha,
      gamma = NA_real_,
      lambda_index = seq_len(path_length),
      lambda = lambda,
      lambda_fraction = lambda / lambda[1L],
      lambda_reference = lambda_reference,
      lambda_reference_type = lambda_reference_type,
      lambda_relative_to_reference = lambda / lambda_reference,
      alpha_zero_ridge_boundary = alpha_zero,
      group_selection_capable = !alpha_zero,
      validation_log_loss = validation_loss,
      passes = NA_real_,
      solver_warning = paste(unique(captured_warnings), collapse = " | "),
      unexpected_solver_warning = paste(
        unique(captured_warnings), collapse = " | "
      ),
      requested_path_length = requested_path_length,
      returned_path_length = path_length,
      path_complete = path_complete,
      saturated_path_truncation = FALSE,
      iteration_budget_truncation = FALSE,
      total_iterations = NA_real_,
      total_iteration_limit_reached = FALSE,
      returned_lower_boundary_excluded = FALSE,
      path_termination_acceptable = path_termination_acceptable,
      point_converged = path_termination_acceptable,
      finite_validation_loss = finite_validation_loss,
      numerically_eligible = numerically_eligible,
      stringsAsFactors = FALSE
    )
  }

  tuning <- do.call(rbind, tuning)
  rownames(tuning) <- NULL
  method_path <- "Logistic Group Elastic Net (adelie)"
  selection <- select_valid_external_grid_v2(tuning, method_path)
  tuning <- mark_external_selection_v2(tuning, selection)

  list(
    fits = fits,
    coefficient_paths = coefficient_paths,
    preprocess = preprocess,
    alpha_grid = input$alpha_grid,
    tuning = tuning,
    selection = selection,
    elapsed_by_alpha = elapsed,
    runtime_seconds = sum(elapsed),
    warning_messages = warning_messages,
    tolerance = tolerance,
    max_iterations = max_iterations,
    irls_tolerance = irls_tolerance,
    irls_max_iterations = irls_max_iterations,
    ridge_boundary_upper_multiplier = ridge_boundary_upper_multiplier,
    null_score_scale = null_score_scale,
    lambda_relative_grid = relative_grid
  )
}


# The frozen replication driver resolves these names at runtime. Rebinding
# them after the v4 utilities are sourced keeps all non-lambda logic unchanged.
fit_extended_sglasso_validation_grid_v1 <-
  fit_extended_sglasso_validation_grid_v5
fit_adelie_validation_grid_v2 <- fit_adelie_validation_grid_v5
