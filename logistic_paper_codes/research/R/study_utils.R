# Utilities for the prespecified Logistic SGLASSO calibration and main study.

read_study_design <- function(path) {
  design <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required <- c(
    "scenario",
    "n_train",
    "n_test",
    "groups",
    "group_size",
    "active_groups",
    "rho_within",
    "rho_between",
    "prevalence",
    "linear_predictor_sd",
    "description"
  )
  missing_columns <- setdiff(required, names(design))
  if (length(missing_columns)) {
    stop(
      "Study design is missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(design) == 0L || anyDuplicated(design$scenario)) {
    stop("Study scenarios must be non-empty and unique.", call. = FALSE)
  }
  integer_columns <- c("n_train", "n_test", "groups", "group_size")
  for (column in integer_columns) {
    design[[column]] <- as.integer(design[[column]])
  }
  numeric_columns <- c(
    "rho_within",
    "rho_between",
    "prevalence",
    "linear_predictor_sd"
  )
  for (column in numeric_columns) {
    design[[column]] <- as.numeric(design[[column]])
  }
  design$active_group_list <- I(lapply(
    strsplit(design$active_groups, ";", fixed = TRUE),
    function(value) as.integer(value[nzchar(value)])
  ))

  if (any(!is.finite(as.matrix(design[, c(integer_columns, numeric_columns)])))) {
    stop("Study design contains non-finite numeric values.", call. = FALSE)
  }
  if (any(design$n_train < 20L) || any(design$n_test < 20L) ||
      any(design$groups < 2L) || any(design$group_size < 1L)) {
    stop("Study sample sizes or group dimensions are invalid.", call. = FALSE)
  }
  if (any(design$rho_within < 0 | design$rho_within >= 1) ||
      any(design$rho_between < 0 | design$rho_between >= 1)) {
    stop("Study correlations must lie in [0, 1).", call. = FALSE)
  }
  if (any(design$prevalence <= 0 | design$prevalence >= 1) ||
      any(design$linear_predictor_sd <= 0)) {
    stop("Study prevalence and signal values are invalid.", call. = FALSE)
  }
  for (i in seq_len(nrow(design))) {
    active <- design$active_group_list[[i]]
    if (!length(active) || anyDuplicated(active) ||
        any(active < 1L | active > design$groups[i])) {
      stop(
        "Invalid active groups for scenario ",
        design$scenario[i],
        ".",
        call. = FALSE
      )
    }
  }
  design
}


cv_lsg_study_grid <- function(
    X,
    y,
    group,
    fold,
    alpha_grid,
    d_grid,
    nlambda,
    lambda_min_ratio,
    max_outer = 600L,
    max_inner = 3000L,
    tolerance = 1e-7,
    inner_tolerance = 1e-9,
    target_original = NULL
) {
  alpha_grid <- sort(unique(as.numeric(alpha_grid)))
  d_grid <- sort(unique(as.numeric(d_grid)))
  if (!length(alpha_grid) || any(alpha_grid <= 0 | alpha_grid > 1)) {
    stop("alpha_grid must lie in (0, 1].", call. = FALSE)
  }
  if (!length(d_grid) || any(d_grid < 0 | d_grid > 1)) {
    stop("d_grid must lie in [0, 1].", call. = FALSE)
  }

  fits <- vector("list", length(alpha_grid))
  elapsed <- numeric(length(alpha_grid))
  minimum_loss <- numeric(length(alpha_grid))
  for (i in seq_along(alpha_grid)) {
    alpha <- alpha_grid[i]
    fitted_d <- if (abs(alpha - 1) < 1e-12) 0 else d_grid
    elapsed[i] <- system.time({
      fits[[i]] <- cv_logistic_sglasso(
        X,
        y,
        group,
        fold = fold,
        nlambda = nlambda,
        lambda_min_ratio = lambda_min_ratio,
        d = fitted_d,
        alpha = alpha,
        max_outer = max_outer,
        max_inner = max_inner,
        tolerance = tolerance,
        inner_tolerance = inner_tolerance,
        target_original = target_original,
        compile = FALSE
      )
    })[["elapsed"]]
    minimum_loss[i] <- min(fits[[i]]$mean_loss)
  }
  selected_alpha_index <- which.min(minimum_loss)
  list(
    selected = fits[[selected_alpha_index]],
    alpha_min = alpha_grid[selected_alpha_index],
    selected_alpha_index = selected_alpha_index,
    fits = fits,
    alpha_grid = alpha_grid,
    d_grid = d_grid,
    minimum_loss_by_alpha = minimum_loss,
    elapsed_by_alpha = elapsed,
    runtime_seconds = sum(elapsed)
  )
}


fit_ridge_path <- function(X, y, active_variables, ridge_lambda) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("glmnet is required for ridge comparators.", call. = FALSE)
  }
  X <- as.matrix(X)
  active_variables <- sort(unique(as.integer(active_variables)))
  ridge_lambda <- sort(unique(as.numeric(ridge_lambda)), decreasing = TRUE)
  if (!length(ridge_lambda) || any(!is.finite(ridge_lambda)) ||
      any(ridge_lambda <= 0)) {
    stop("ridge_lambda must contain finite positive values.", call. = FALSE)
  }
  if (!length(active_variables)) {
    prevalence <- min(max(mean(y), 1e-8), 1 - 1e-8)
    return(structure(
      list(
        empty = TRUE,
        intercept = stats::qlogis(prevalence),
        active_variables = integer(0),
        p_original = ncol(X),
        lambda = ridge_lambda
      ),
      class = "study_ridge_path"
    ))
  }
  fit <- glmnet::glmnet(
    x = X[, active_variables, drop = FALSE],
    y = as.numeric(y),
    family = "binomial",
    alpha = 0,
    lambda = ridge_lambda,
    standardize = TRUE,
    intercept = TRUE,
    control = list(
      thresh = 1e-9,
      maxit = 100000L
    )
  )
  structure(
    list(
      empty = FALSE,
      fit = fit,
      active_variables = active_variables,
      p_original = ncol(X),
      lambda = ridge_lambda
    ),
    class = "study_ridge_path"
  )
}


predict_study_ridge <- function(object, newx, s = object$lambda) {
  if (!inherits(object, "study_ridge_path")) {
    stop("object must be a study_ridge_path.", call. = FALSE)
  }
  newx <- as.matrix(newx)
  s <- as.numeric(s)
  if (object$empty) {
    return(matrix(
      stats::plogis(object$intercept),
      nrow = nrow(newx),
      ncol = length(s)
    ))
  }
  as.matrix(stats::predict(
    object$fit,
    newx = newx[, object$active_variables, drop = FALSE],
    type = "response",
    s = s
  ))
}


coef_study_ridge <- function(object, s) {
  if (!inherits(object, "study_ridge_path")) {
    stop("object must be a study_ridge_path.", call. = FALSE)
  }
  if (length(s) != 1L) stop("s must be scalar.", call. = FALSE)
  out <- numeric(object$p_original + 1L)
  if (object$empty) {
    out[1L] <- object$intercept
    return(out)
  }
  coefficient <- as.matrix(stats::coef(object$fit, s = s))[, 1L]
  out[1L] <- coefficient[1L]
  out[1L + object$active_variables] <- coefficient[-1L]
  out
}


fit_complete_group_lasso_path <- function(
    X,
    y,
    group,
    nlambda,
    lambda_min_ratio,
    tolerance
) {
  tolerance_sequence <- unique(c(
    tolerance,
    max(tolerance, 1e-5),
    max(tolerance, 1e-4)
  ))
  path_lengths <- integer(length(tolerance_sequence))
  for (i in seq_along(tolerance_sequence)) {
    current_tolerance <- tolerance_sequence[i]
    fit <- grpreg::grpreg(
      X,
      y,
      group,
      penalty = "grLasso",
      family = "binomial",
      nlambda = nlambda,
      lambda.min = lambda_min_ratio,
      eps = current_tolerance,
      max.iter = 100000L,
      dfmax = ncol(X),
      gmax = length(unique(group)),
      warn = FALSE
    )
    path_lengths[i] <- length(fit$lambda)
    if (path_lengths[i] == nlambda) {
      return(list(
        fit = fit,
        tolerance_used = current_tolerance,
        attempted_tolerances = tolerance_sequence[seq_len(i)],
        attempted_path_lengths = path_lengths[seq_len(i)]
      ))
    }
  }
  stop(
    "grpreg returned incomplete paths at tolerances ",
    paste(tolerance_sequence, collapse = ", "),
    " (lengths ",
    paste(path_lengths, collapse = ", "),
    ").",
    call. = FALSE
  )
}


cv_group_lasso_post_ridge <- function(
    X,
    y,
    group,
    fold,
    nlambda = 30L,
    lambda_min_ratio = 0.02,
    ridge_lambda = 10^seq(1, -4, length.out = 15L),
    tolerance = 1e-7,
    selection_threshold = 1e-8
) {
  if (!requireNamespace("grpreg", quietly = TRUE)) {
    stop("grpreg is required for the group-lasso comparator.", call. = FALSE)
  }
  X <- as.matrix(X)
  y <- normalize_binary_response(y)$y
  group_character <- as.character(group)
  group <- match(group_character, unique(group_character))
  fold_labels <- sort(unique(fold))
  ridge_lambda <- sort(unique(as.numeric(ridge_lambda)), decreasing = TRUE)
  group_loss <- matrix(
    NA_real_,
    nrow = nlambda,
    ncol = length(fold_labels)
  )
  post_loss <- array(
    NA_real_,
    dim = c(nlambda, length(ridge_lambda), length(fold_labels))
  )
  group_lasso_seconds <- 0
  ridge_seconds <- 0
  fold_tolerance_used <- numeric(length(fold_labels))

  for (fi in seq_along(fold_labels)) {
    held_out <- fold == fold_labels[fi]
    group_lasso_seconds <- group_lasso_seconds + system.time({
      fit_information <- fit_complete_group_lasso_path(
        X[!held_out, , drop = FALSE],
        y[!held_out],
        group,
        nlambda,
        lambda_min_ratio,
        tolerance
      )
      fit <- fit_information$fit
      group_probability <- as.matrix(stats::predict(
        fit,
        X[held_out, , drop = FALSE],
        type = "response"
      ))
    })[["elapsed"]]
    fold_tolerance_used[fi] <- fit_information$tolerance_used
    for (li in seq_len(nlambda)) {
      group_loss[li, fi] <- binary_log_loss(
        y[held_out],
        group_probability[, li]
      )
    }

    coefficient_path <- as.matrix(stats::coef(fit))[-1L, , drop = FALSE]
    active_sets <- lapply(
      seq_len(nlambda),
      function(li) which(abs(coefficient_path[, li]) > selection_threshold)
    )
    signatures <- vapply(
      active_sets,
      function(index) paste(index, collapse = ","),
      character(1)
    )
    for (signature in unique(signatures)) {
      lambda_indices <- which(signatures == signature)
      active <- active_sets[[lambda_indices[1L]]]
      ridge_seconds <- ridge_seconds + system.time({
        ridge_fit <- fit_ridge_path(
          X[!held_out, , drop = FALSE],
          y[!held_out],
          active,
          ridge_lambda
        )
        ridge_probability <- predict_study_ridge(
          ridge_fit,
          X[held_out, , drop = FALSE],
          s = ridge_lambda
        )
      })[["elapsed"]]
      ridge_loss <- vapply(
        seq_along(ridge_lambda),
        function(ri) binary_log_loss(
          y[held_out],
          ridge_probability[, ri]
        ),
        numeric(1)
      )
      for (li in lambda_indices) post_loss[li, , fi] <- ridge_loss
    }
  }

  mean_group_loss <- rowMeans(group_loss)
  group_index_min <- which.min(mean_group_loss)
  mean_post_loss <- apply(post_loss, c(1, 2), mean)
  post_minimum <- which(
    mean_post_loss == min(mean_post_loss),
    arr.ind = TRUE
  )[1L, ]

  group_lasso_seconds <- group_lasso_seconds + system.time({
    full_fit_information <- fit_complete_group_lasso_path(
      X,
      y,
      group,
      nlambda,
      lambda_min_ratio,
      tolerance
    )
    full_fit <- full_fit_information$fit
  })[["elapsed"]]
  full_coefficient_path <- as.matrix(stats::coef(full_fit))
  post_active <- which(
    abs(full_coefficient_path[-1L, post_minimum[1L]]) >
      selection_threshold
  )
  ridge_seconds <- ridge_seconds + system.time({
    full_post_fit <- fit_ridge_path(X, y, post_active, ridge_lambda)
  })[["elapsed"]]

  list(
    group_fit = full_fit,
    group_fold_loss = group_loss,
    group_mean_loss = mean_group_loss,
    group_index_min = group_index_min,
    group_lambda_min = full_fit$lambda[group_index_min],
    group_lambda_fraction = full_fit$lambda / full_fit$lambda[1L],
    post_fold_loss = post_loss,
    post_mean_loss = mean_post_loss,
    post_group_index_min = post_minimum[1L],
    post_ridge_index_min = post_minimum[2L],
    post_group_lambda_min = full_fit$lambda[post_minimum[1L]],
    post_ridge_lambda_min = ridge_lambda[post_minimum[2L]],
    post_active_variables = post_active,
    post_fit = full_post_fit,
    ridge_lambda = ridge_lambda,
    group_runtime_seconds = group_lasso_seconds,
    post_runtime_seconds = group_lasso_seconds + ridge_seconds,
    ridge_increment_seconds = ridge_seconds,
    grpreg_fold_tolerance_used = fold_tolerance_used,
    grpreg_full_tolerance_used = full_fit_information$tolerance_used,
    grpreg_fallback_count = sum(fold_tolerance_used > tolerance) +
      as.integer(full_fit_information$tolerance_used > tolerance),
    grpreg_maximum_tolerance_used = max(c(
      fold_tolerance_used,
      full_fit_information$tolerance_used
    ))
  )
}


predict_selected_group_lasso <- function(object, newx) {
  drop(as.matrix(stats::predict(
    object$group_fit,
    as.matrix(newx),
    type = "response"
  ))[, object$group_index_min])
}


coef_selected_group_lasso <- function(object) {
  as.matrix(stats::coef(object$group_fit))[, object$group_index_min]
}


predict_selected_post_ridge <- function(object, newx) {
  drop(predict_study_ridge(
    object$post_fit,
    newx,
    s = object$post_ridge_lambda_min
  ))
}


coef_selected_post_ridge <- function(object) {
  coef_study_ridge(object$post_fit, object$post_ridge_lambda_min)
}


cv_oracle_ridge <- function(
    X,
    y,
    active_variables,
    fold,
    ridge_lambda = 10^seq(1, -4, length.out = 15L)
) {
  X <- as.matrix(X)
  y <- normalize_binary_response(y)$y
  fold_labels <- sort(unique(fold))
  ridge_lambda <- sort(unique(as.numeric(ridge_lambda)), decreasing = TRUE)
  fold_loss <- matrix(
    NA_real_,
    nrow = length(ridge_lambda),
    ncol = length(fold_labels)
  )
  elapsed <- 0
  for (fi in seq_along(fold_labels)) {
    held_out <- fold == fold_labels[fi]
    elapsed <- elapsed + system.time({
      fit <- fit_ridge_path(
        X[!held_out, , drop = FALSE],
        y[!held_out],
        active_variables,
        ridge_lambda
      )
      probability <- predict_study_ridge(
        fit,
        X[held_out, , drop = FALSE],
        s = ridge_lambda
      )
    })[["elapsed"]]
    for (ri in seq_along(ridge_lambda)) {
      fold_loss[ri, fi] <- binary_log_loss(y[held_out], probability[, ri])
    }
  }
  mean_loss <- rowMeans(fold_loss)
  selected <- which.min(mean_loss)
  elapsed <- elapsed + system.time({
    full_fit <- fit_ridge_path(X, y, active_variables, ridge_lambda)
  })[["elapsed"]]
  list(
    fit = full_fit,
    fold_loss = fold_loss,
    mean_loss = mean_loss,
    ridge_index_min = selected,
    ridge_lambda_min = ridge_lambda[selected],
    ridge_lambda = ridge_lambda,
    runtime_seconds = elapsed
  )
}


predict_selected_oracle_ridge <- function(object, newx) {
  drop(predict_study_ridge(
    object$fit,
    newx,
    s = object$ridge_lambda_min
  ))
}


coef_selected_oracle_ridge <- function(object) {
  coef_study_ridge(object$fit, object$ridge_lambda_min)
}


selected_groups_from_coefficients <- function(
    coefficient,
    group,
    threshold = 1e-8
) {
  coefficient <- as.numeric(coefficient)
  group <- as.character(group)
  active <- tapply(abs(coefficient) > threshold, group, any)
  names(active)[active]
}


extract_lsg_solution <- function(cv_fit, newx) {
  coefficient <- cv_fit$fit$coefficients[
    ,
    cv_fit$lambda_index_min,
    cv_fit$d_index_min
  ]
  list(
    probability = predict_cv_logistic_sglasso(cv_fit, newx),
    coefficient = as.numeric(coefficient),
    selected_groups = selected_lsg_groups(cv_fit)
  )
}


calibration_statistics <- function(y, probability) {
  probability <- pmin(pmax(as.numeric(probability), 1e-8), 1 - 1e-8)
  link <- stats::qlogis(probability)
  intercept_fit <- suppressWarnings(try(
    stats::glm.fit(
      x = matrix(1, nrow = length(y), ncol = 1L),
      y = as.numeric(y),
      family = stats::binomial(),
      offset = link
    ),
    silent = TRUE
  ))
  slope_fit <- suppressWarnings(try(
    stats::glm.fit(
      x = cbind(1, link),
      y = as.numeric(y),
      family = stats::binomial()
    ),
    silent = TRUE
  ))
  intercept <- if (inherits(intercept_fit, "try-error")) {
    NA_real_
  } else {
    unname(intercept_fit$coefficients[1L])
  }
  slope <- if (inherits(slope_fit, "try-error")) {
    NA_real_
  } else {
    unname(slope_fit$coefficients[2L])
  }
  c(intercept = intercept, slope = slope)
}


evaluate_study_method <- function(
    y,
    probability,
    selected_groups,
    active_groups,
    coefficient,
    true_coefficient,
    X
) {
  coefficient <- as.numeric(coefficient)
  true_coefficient <- as.numeric(true_coefficient)
  if (length(coefficient) != length(true_coefficient) ||
      length(coefficient) != ncol(X) + 1L) {
    stop("Coefficient vectors have incompatible dimensions.", call. = FALSE)
  }
  selected_groups <- unique(as.character(selected_groups))
  active_groups <- unique(as.character(active_groups))
  base <- evaluate_binary_method(
    y,
    probability,
    selected_groups,
    active_groups
  )
  calibration <- calibration_statistics(y, probability)
  true_positive <- length(intersect(selected_groups, active_groups))
  false_positive <- length(setdiff(selected_groups, active_groups))
  exact_group_support_value <-
    true_positive == length(active_groups) && false_positive == 0L
  estimated_link <- coefficient[1L] + drop(X %*% coefficient[-1L])
  true_link <- true_coefficient[1L] + drop(X %*% true_coefficient[-1L])
  transform(
    base,
    true_positive_groups = true_positive,
    false_positive_groups = false_positive,
    # Use a uniquely named scalar: base contains a numeric selected_groups
    # column that would mask the selected_groups function argument here.
    exact_group_support = exact_group_support_value,
    coefficient_l2_error = sqrt(sum(
      (coefficient[-1L] - true_coefficient[-1L])^2
    )),
    coefficient_l1_error = sum(abs(
      coefficient[-1L] - true_coefficient[-1L]
    )),
    intercept_absolute_error = abs(coefficient[1L] - true_coefficient[1L]),
    linear_predictor_rmse = sqrt(mean((estimated_link - true_link)^2)),
    calibration_intercept = calibration[["intercept"]],
    calibration_slope = calibration[["slope"]]
  )
}


lsg_numerical_diagnostics <- function(grid) {
  cv_fit <- grid$selected
  li <- cv_fit$lambda_index_min
  di <- cv_fit$d_index_min
  data.frame(
    selected_alpha = grid$alpha_min,
    selected_d = cv_fit$d_min,
    selected_lambda_fraction = cv_fit$lambda_fraction[li],
    selected_full_kkt = cv_fit$fit$kkt[li, di],
    selected_fold_maximum_kkt = max(cv_fit$fold_kkt[li, di, ]),
    selected_full_converged = cv_fit$fit$converged[li, di],
    selected_folds_all_converged = all(cv_fit$fold_converged[li, di, ]),
    full_path_maximum_kkt = max(cv_fit$fit$kkt),
    fold_paths_maximum_kkt = max(cv_fit$fold_kkt),
    full_path_convergence_rate = mean(cv_fit$fit$converged),
    fold_paths_convergence_rate = mean(cv_fit$fold_converged),
    maximum_raw_objective_increase = max(c(
      cv_fit$fit$maximum_raw_objective_increase,
      cv_fit$fold_maximum_raw_objective_increase
    )),
    selected_lambda_on_boundary = li %in% c(1L, length(cv_fit$fit$lambda)),
    selected_d_on_boundary = cv_fit$d_min %in% c(0, 1),
    selected_d_zero_model_feasible =
      cv_fit$fit$zero_model_feasible[di],
    first_path_selected_groups = cv_fit$fit$selected_groups[1L, di],
    stringsAsFactors = FALSE
  )
}


run_study_replication <- function(scenario, replication, configuration) {
  active_groups <- scenario$active_group_list[[1L]]
  seed <- as.integer(
    configuration$seed_base +
      scenario$scenario_index * 100000L +
      replication
  )
  data <- simulate_grouped_logistic(
    n_train = scenario$n_train,
    n_test = scenario$n_test,
    groups = scenario$groups,
    group_size = scenario$group_size,
    active_groups = active_groups,
    rho_within = scenario$rho_within,
    rho_between = scenario$rho_between,
    prevalence = scenario$prevalence,
    linear_predictor_sd = scenario$linear_predictor_sd,
    seed = seed
  )
  fold <- stratified_folds(
    data$y_train,
    nfolds = configuration$nfolds,
    seed = seed + 50000L
  )

  lsg <- cv_lsg_study_grid(
    data$X_train,
    data$y_train,
    data$group,
    fold,
    alpha_grid = configuration$alpha_grid,
    d_grid = configuration$d_grid,
    nlambda = configuration$nlambda,
    lambda_min_ratio = configuration$lambda_min_ratio,
    max_outer = configuration$max_outer,
    max_inner = configuration$max_inner,
    tolerance = configuration$tolerance,
    inner_tolerance = configuration$inner_tolerance
  )
  d0 <- select_lsg_d_value(lsg, d_value = 0)
  grouped <- cv_group_lasso_post_ridge(
    data$X_train,
    data$y_train,
    data$group,
    fold,
    nlambda = configuration$nlambda,
    lambda_min_ratio = configuration$lambda_min_ratio,
    ridge_lambda = configuration$ridge_lambda,
    tolerance = configuration$grpreg_tolerance
  )
  active_variables <- which(as.character(data$group) %in% data$active_groups)
  oracle <- cv_oracle_ridge(
    data$X_train,
    data$y_train,
    active_variables,
    fold,
    ridge_lambda = configuration$ridge_lambda
  )

  lsg_solution <- extract_lsg_solution(lsg$selected, data$X_test)
  d0_solution <- extract_lsg_solution(d0$selected, data$X_test)
  group_coefficient <- coef_selected_group_lasso(grouped)
  group_solution <- list(
    probability = predict_selected_group_lasso(grouped, data$X_test),
    coefficient = group_coefficient,
    selected_groups = selected_groups_from_coefficients(
      group_coefficient[-1L],
      data$group
    )
  )
  post_coefficient <- coef_selected_post_ridge(grouped)
  post_solution <- list(
    probability = predict_selected_post_ridge(grouped, data$X_test),
    coefficient = post_coefficient,
    selected_groups = unique(as.character(
      data$group[grouped$post_active_variables]
    ))
  )
  oracle_solution <- list(
    probability = predict_selected_oracle_ridge(oracle, data$X_test),
    coefficient = coef_selected_oracle_ridge(oracle),
    selected_groups = data$active_groups
  )

  true_coefficient <- c(data$intercept, data$beta)
  solutions <- list(
    "Logistic SGLASSO" = lsg_solution,
    "Group elastic net (d=0)" = d0_solution,
    "Logistic group lasso" = group_solution,
    "Post-group-lasso ridge" = post_solution,
    "Oracle active-group ridge" = oracle_solution
  )
  runtimes <- c(
    "Logistic SGLASSO" = lsg$runtime_seconds,
    "Group elastic net (d=0)" = NA_real_,
    "Logistic group lasso" = grouped$group_runtime_seconds,
    "Post-group-lasso ridge" = grouped$post_runtime_seconds,
    "Oracle active-group ridge" = oracle$runtime_seconds
  )
  runtime_scopes <- c(
    "Logistic SGLASSO" = "full_alpha_d_cv_grid",
    "Group elastic net (d=0)" = "shared_with_sglasso_grid",
    "Logistic group lasso" = "group_lasso_cv_only",
    "Post-group-lasso ridge" = "joint_group_lasso_ridge_cv",
    "Oracle active-group ridge" = "oracle_ridge_cv_only"
  )
  result <- do.call(
    rbind,
    lapply(names(solutions), function(method) {
      solution <- solutions[[method]]
      reconstructed_probability <- stats::plogis(
        solution$coefficient[1L] +
          drop(data$X_test %*% solution$coefficient[-1L])
      )
      prediction_error <- max(abs(
        reconstructed_probability - solution$probability
      ))
      if (!is.finite(prediction_error) || prediction_error > 2e-6 ||
          any(!is.finite(solution$coefficient)) ||
          any(!is.finite(solution$probability))) {
        stop(
          "Coefficient/prediction consistency failed for ",
          method,
          ".",
          call. = FALSE
        )
      }
      metrics <- evaluate_study_method(
        data$y_test,
        solution$probability,
        solution$selected_groups,
        data$active_groups,
        solution$coefficient,
        true_coefficient,
        data$X_test
      )
      cbind(
        data.frame(
          scenario = scenario$scenario,
          replication = replication,
          seed = seed,
          method = method,
          n_train = scenario$n_train,
          p = scenario$groups * scenario$group_size,
          groups = scenario$groups,
          group_size = scenario$group_size,
          active_group_count = length(active_groups),
          target_prevalence = scenario$prevalence,
          train_prevalence = mean(data$y_train),
          test_prevalence = mean(data$y_test),
          rho_within = scenario$rho_within,
          rho_between = scenario$rho_between,
          linear_predictor_sd = scenario$linear_predictor_sd,
          runtime_seconds = unname(runtimes[[method]]),
          runtime_scope = unname(runtime_scopes[[method]]),
          stringsAsFactors = FALSE
        ),
        metrics
      )
    })
  )
  rownames(result) <- NULL

  lsg_diagnostic <- lsg_numerical_diagnostics(lsg)
  names(lsg_diagnostic) <- paste0("lsg_", names(lsg_diagnostic))
  d0_diagnostic <- lsg_numerical_diagnostics(d0)
  names(d0_diagnostic) <- paste0("d0_", names(d0_diagnostic))
  tuning <- cbind(
    data.frame(
      scenario = scenario$scenario,
      replication = replication,
      seed = seed,
      train_events = sum(data$y_train == 1),
      train_nonevents = sum(data$y_train == 0),
      group_lasso_lambda_fraction =
        grouped$group_lambda_fraction[grouped$group_index_min],
      post_group_lasso_lambda_fraction =
        grouped$group_lambda_fraction[grouped$post_group_index_min],
      post_ridge_lambda = grouped$post_ridge_lambda_min,
      oracle_ridge_lambda = oracle$ridge_lambda_min,
      stringsAsFactors = FALSE
    ),
    lsg_diagnostic,
    d0_diagnostic
  )
  tuning$grpreg_fallback_count <- grouped$grpreg_fallback_count
  tuning$grpreg_maximum_tolerance_used <-
    grouped$grpreg_maximum_tolerance_used
  list(results = result, tuning = tuning)
}


summarize_study_results <- function(results) {
  metrics <- c(
    "log_loss",
    "brier",
    "auc",
    "classification_error",
    "selected_groups",
    "group_tpr",
    "group_fdr",
    "exact_group_support",
    "coefficient_l2_error",
    "coefficient_l1_error",
    "intercept_absolute_error",
    "linear_predictor_rmse",
    "calibration_intercept",
    "calibration_slope",
    "runtime_seconds"
  )
  frames <- split(
    results,
    interaction(results$scenario, results$method, drop = TRUE)
  )
  out <- do.call(rbind, lapply(frames, function(frame) {
    row <- data.frame(
      scenario = frame$scenario[1L],
      method = frame$method[1L],
      replications = nrow(frame),
      stringsAsFactors = FALSE
    )
    for (metric in metrics) {
      value <- as.numeric(frame[[metric]])
      finite <- is.finite(value)
      row[[paste0(metric, "_n")]] <- sum(finite)
      row[[paste0(metric, "_mean")]] <- if (any(finite)) {
        mean(value[finite])
      } else {
        NA_real_
      }
      row[[paste0(metric, "_se")]] <- if (sum(finite) > 1L) {
        stats::sd(value[finite]) / sqrt(sum(finite))
      } else {
        NA_real_
      }
    }
    row
  }))
  rownames(out) <- NULL
  out
}


paired_study_log_loss <- function(results) {
  reference <- results[
    results$method == "Logistic group lasso",
    c("scenario", "replication", "log_loss")
  ]
  names(reference)[3L] <- "reference_log_loss"
  methods <- setdiff(unique(results$method), "Logistic group lasso")
  paired <- do.call(rbind, lapply(methods, function(method) {
    candidate <- results[
      results$method == method,
      c("scenario", "replication", "log_loss")
    ]
    names(candidate)[3L] <- "candidate_log_loss"
    merged <- merge(
      candidate,
      reference,
      by = c("scenario", "replication"),
      all = FALSE,
      sort = FALSE
    )
    data.frame(
      scenario = merged$scenario,
      replication = merged$replication,
      method = method,
      difference_vs_group_lasso =
        merged$candidate_log_loss - merged$reference_log_loss,
      stringsAsFactors = FALSE
    )
  }))
  rownames(paired) <- NULL
  paired
}


summarize_paired_study <- function(paired) {
  frames <- split(
    paired,
    interaction(paired$scenario, paired$method, drop = TRUE)
  )
  out <- do.call(rbind, lapply(frames, function(frame) {
    difference <- frame$difference_vs_group_lasso
    n <- length(difference)
    standard_error <- if (n > 1L) {
      stats::sd(difference) / sqrt(n)
    } else {
      NA_real_
    }
    half_width <- if (n > 1L) {
      stats::qt(0.975, df = n - 1L) * standard_error
    } else {
      NA_real_
    }
    data.frame(
      scenario = frame$scenario[1L],
      method = frame$method[1L],
      replications = n,
      mean_difference_vs_group_lasso = mean(difference),
      se_difference_vs_group_lasso = standard_error,
      lower_95 = mean(difference) - half_width,
      upper_95 = mean(difference) + half_width,
      proportion_better_than_group_lasso = mean(difference < 0),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}


atomic_save_rds <- function(object, path) {
  directory <- dirname(path)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = "checkpoint_", tmpdir = directory)
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) {
    if (!file.copy(temporary, path, overwrite = TRUE)) {
      stop("Unable to save checkpoint: ", path, call. = FALSE)
    }
    unlink(temporary)
  }
  invisible(path)
}
