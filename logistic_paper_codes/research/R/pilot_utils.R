calibrate_logistic_intercept <- function(prevalence, linear_predictor_sd) {
  if (!is.finite(prevalence) || prevalence <= 0 || prevalence >= 1) {
    stop("prevalence must lie in (0, 1).", call. = FALSE)
  }
  expectation <- function(intercept) {
    stats::integrate(
      function(z) {
        stats::plogis(intercept + linear_predictor_sd * z) *
          stats::dnorm(z)
      },
      lower = -Inf,
      upper = Inf,
      rel.tol = 1e-10
    )$value
  }
  stats::uniroot(
    function(intercept) expectation(intercept) - prevalence,
    interval = c(-20, 20),
    tol = 1e-10
  )$root
}


grouped_covariance <- function(
    groups,
    group_size,
    rho_within,
    rho_between
) {
  p <- groups * group_size
  group_index <- rep(seq_len(groups), each = group_size)
  variable_index <- rep(seq_len(group_size), times = groups)
  covariance <- matrix(0, p, p)
  for (j in seq_len(p)) {
    for (k in seq_len(p)) {
      if (j == k) {
        covariance[j, k] <- 1
      } else if (group_index[j] == group_index[k]) {
        covariance[j, k] <- rho_within
      } else {
        covariance[j, k] <-
          rho_within *
          rho_between^abs(group_index[j] - group_index[k])
      }
    }
  }
  covariance
}


simulate_grouped_logistic <- function(
    n_train = 180L,
    n_test = 600L,
    groups = 30L,
    group_size = 3L,
    active_groups = c(1L, 2L, 8L, 15L, 22L),
    rho_within = 0.8,
    rho_between = 0,
    prevalence = 0.5,
    linear_predictor_sd = 1.1,
    seed = 1L
) {
  set.seed(seed)
  n_total <- n_train + n_test
  latent_covariance <- outer(
    seq_len(groups),
    seq_len(groups),
    function(j, k) rho_between^abs(j - k)
  )
  latent <- matrix(
    stats::rnorm(n_total * groups),
    nrow = n_total,
    ncol = groups
  ) %*% chol(latent_covariance)

  p <- groups * group_size
  X <- matrix(0, nrow = n_total, ncol = p)
  group <- rep(seq_len(groups), each = group_size)
  for (g in seq_len(groups)) {
    for (k in seq_len(group_size)) {
      column <- (g - 1L) * group_size + k
      X[, column] <-
        sqrt(rho_within) * latent[, g] +
        sqrt(1 - rho_within) * stats::rnorm(n_total)
    }
  }
  colnames(X) <- paste0("G", group, "_V", rep(seq_len(group_size), groups))

  beta <- numeric(p)
  base_pattern <- rep(c(1, -0.8, 0.6), length.out = group_size)
  for (position in seq_along(active_groups)) {
    g <- active_groups[position]
    index <- which(group == g)
    beta[index] <- (-1)^(position + 1L) * base_pattern
  }
  covariance <- grouped_covariance(
    groups,
    group_size,
    rho_within,
    rho_between
  )
  current_sd <- sqrt(drop(crossprod(beta, covariance %*% beta)))
  beta <- beta * linear_predictor_sd / current_sd
  intercept <- calibrate_logistic_intercept(
    prevalence,
    linear_predictor_sd
  )
  probability <- stats::plogis(intercept + drop(X %*% beta))
  y <- stats::rbinom(n_total, 1, probability)

  list(
    X_train = X[seq_len(n_train), , drop = FALSE],
    y_train = y[seq_len(n_train)],
    X_test = X[n_train + seq_len(n_test), , drop = FALSE],
    y_test = y[n_train + seq_len(n_test)],
    group = group,
    beta = beta,
    intercept = intercept,
    active_groups = as.character(active_groups),
    parameters = list(
      n_train = n_train,
      n_test = n_test,
      groups = groups,
      group_size = group_size,
      active_groups = active_groups,
      rho_within = rho_within,
      rho_between = rho_between,
      target_prevalence = prevalence,
      realized_train_prevalence = mean(y[seq_len(n_train)]),
      realized_test_prevalence = mean(y[n_train + seq_len(n_test)]),
      linear_predictor_sd = linear_predictor_sd,
      seed = seed
    )
  )
}


cv_lsg_alpha_grid <- function(
    X,
    y,
    group,
    fold,
    alpha_grid = c(0.5, 0.8),
    d = c(0, 0.5, 1),
    nlambda = 20L,
    lambda_min_ratio = 0.04,
    tolerance = 1e-6,
    inner_tolerance = 1e-8
) {
  fits <- lapply(
    alpha_grid,
    function(alpha) {
      cv_logistic_sglasso(
        X,
        y,
        group,
        fold = fold,
        nlambda = nlambda,
        lambda_min_ratio = lambda_min_ratio,
        d = d,
        alpha = alpha,
        tolerance = tolerance,
        inner_tolerance = inner_tolerance,
        compile = FALSE
      )
    }
  )
  minima <- vapply(fits, function(fit) min(fit$mean_loss), numeric(1))
  selected <- which.min(minima)
  list(
    selected = fits[[selected]],
    alpha_min = alpha_grid[selected],
    alpha_grid = alpha_grid,
    fits = fits,
    minimum_loss_by_alpha = minima
  )
}


select_lsg_d_value <- function(alpha_grid_fit, d_value = 0) {
  candidates <- lapply(
    alpha_grid_fit$fits,
    function(cv_fit) {
      d_index <- which.min(abs(cv_fit$fit$d - d_value))
      if (abs(cv_fit$fit$d[d_index] - d_value) > 1e-12) {
        stop("Requested d value was not fitted.", call. = FALSE)
      }
      lambda_index <- which.min(cv_fit$mean_loss[, d_index])
      list(
        cv_fit = cv_fit,
        d_index = d_index,
        lambda_index = lambda_index,
        loss = cv_fit$mean_loss[lambda_index, d_index]
      )
    }
  )
  selected_alpha_index <- which.min(
    vapply(candidates, function(item) item$loss, numeric(1))
  )
  candidate <- candidates[[selected_alpha_index]]
  selected <- candidate$cv_fit
  selected$d_index_min <- candidate$d_index
  selected$lambda_index_min <- candidate$lambda_index
  selected$d_min <- selected$fit$d[candidate$d_index]
  selected$lambda_min <- selected$fit$lambda[candidate$lambda_index]
  list(
    selected = selected,
    alpha_min = alpha_grid_fit$alpha_grid[selected_alpha_index],
    cv_loss_min = candidate$loss
  )
}


cv_grpreg_manual <- function(
    X,
    y,
    group,
    fold,
    nlambda = 20L,
    lambda_min_ratio = 0.04,
    tolerance = 1e-6
) {
  if (!requireNamespace("grpreg", quietly = TRUE)) {
    stop("grpreg is required for the group-lasso comparator.", call. = FALSE)
  }
  fold_labels <- sort(unique(fold))
  fold_loss <- matrix(
    NA_real_,
    nrow = nlambda,
    ncol = length(fold_labels)
  )
  for (fold_position in seq_along(fold_labels)) {
    held_out <- fold == fold_labels[fold_position]
    fit <- grpreg::grpreg(
      X[!held_out, , drop = FALSE],
      y[!held_out],
      group,
      penalty = "grLasso",
      family = "binomial",
      nlambda = nlambda,
      lambda.min = lambda_min_ratio,
      eps = tolerance,
      max.iter = 100000,
      dfmax = ncol(X),
      gmax = length(unique(group)),
      warn = FALSE
    )
    probability <- stats::predict(
      fit,
      X[held_out, , drop = FALSE],
      type = "response"
    )
    if (ncol(probability) != nlambda) {
      stop("grpreg returned an incomplete lambda path.", call. = FALSE)
    }
    for (li in seq_len(nlambda)) {
      fold_loss[li, fold_position] <- binary_log_loss(
        y[held_out],
        probability[, li]
      )
    }
  }
  mean_loss <- rowMeans(fold_loss)
  selected_index <- which.min(mean_loss)
  full_fit <- grpreg::grpreg(
    X,
    y,
    group,
    penalty = "grLasso",
    family = "binomial",
    nlambda = nlambda,
    lambda.min = lambda_min_ratio,
    eps = tolerance,
    max.iter = 100000,
    dfmax = ncol(X),
    gmax = length(unique(group)),
    warn = FALSE
  )
  list(
    fit = full_fit,
    fold_loss = fold_loss,
    mean_loss = mean_loss,
    lambda_index_min = selected_index,
    lambda_min = full_fit$lambda[selected_index],
    lambda_fraction = full_fit$lambda / full_fit$lambda[1L]
  )
}


selected_lsg_groups <- function(cv_fit, threshold = 1e-8) {
  fit <- cv_fit$fit
  beta <- fit$beta_solver[
    ,
    cv_fit$lambda_index_min,
    cv_fit$d_index_min
  ]
  active <- vapply(
    fit$preprocess$blocks,
    function(block) {
      sqrt(sum(beta[block$solver_index]^2)) > threshold
    },
    logical(1)
  )
  fit$preprocess$group_labels[active]
}


selected_grpreg_groups <- function(cv_fit, group, threshold = 1e-8) {
  coefficient <- stats::coef(cv_fit$fit)[-1L, cv_fit$lambda_index_min]
  group_character <- as.character(group)
  active <- tapply(
    abs(coefficient) > threshold,
    group_character,
    any
  )
  names(active)[active]
}


evaluate_binary_method <- function(
    y,
    probability,
    selected_groups,
    active_groups
) {
  selected_groups <- unique(as.character(selected_groups))
  active_groups <- unique(as.character(active_groups))
  true_positive <- length(intersect(selected_groups, active_groups))
  false_positive <- length(setdiff(selected_groups, active_groups))
  data.frame(
    log_loss = binary_log_loss(y, probability),
    brier = mean((y - probability)^2),
    auc = binary_auc(y, probability),
    classification_error = mean((probability >= 0.5) != y),
    selected_groups = length(selected_groups),
    group_tpr = true_positive / length(active_groups),
    group_fdr = if (length(selected_groups) == 0L) {
      0
    } else {
      false_positive / length(selected_groups)
    },
    stringsAsFactors = FALSE
  )
}


summarize_pilot_results <- function(results) {
  metrics <- c(
    "log_loss",
    "brier",
    "auc",
    "classification_error",
    "selected_groups",
    "group_tpr",
    "group_fdr",
    "runtime_seconds"
  )
  split_results <- split(
    results,
    interaction(results$scenario, results$method, drop = TRUE)
  )
  do.call(
    rbind,
    lapply(
      split_results,
      function(frame) {
        out <- data.frame(
          scenario = frame$scenario[1L],
          method = frame$method[1L],
          replications = nrow(frame),
          stringsAsFactors = FALSE
        )
        for (metric in metrics) {
          out[[paste0(metric, "_mean")]] <- mean(frame[[metric]], na.rm = TRUE)
          out[[paste0(metric, "_se")]] <-
            stats::sd(frame[[metric]], na.rm = TRUE) / sqrt(nrow(frame))
        }
        out
      }
    )
  )
}
