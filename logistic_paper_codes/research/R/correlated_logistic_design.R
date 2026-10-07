# Theory-aligned clustered-group design for the isolated logistic prework.

parse_semicolon_numeric <- function(value, integer = FALSE) {
  out <- as.numeric(strsplit(as.character(value), ";", fixed = TRUE)[[1L]])
  if (!length(out) || any(!is.finite(out))) {
    stop("Invalid semicolon-separated numeric value.", call. = FALSE)
  }
  if (integer) as.integer(out) else out
}


read_aligned_logistic_design <- function(path) {
  design <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c(
    "scenario", "n_train", "n_validation", "n_test", "groups",
    "group_size", "cluster_size", "active_groups", "active_signs",
    "rho_within", "rho_between", "prevalence", "linear_predictor_sd",
    "description"
  )
  if (!identical(names(design), required)) {
    stop("Aligned-design columns do not match the frozen schema.", call. = FALSE)
  }
  integer_columns <- c(
    "n_train", "n_validation", "n_test", "groups", "group_size",
    "cluster_size"
  )
  for (name in integer_columns) design[[name]] <- as.integer(design[[name]])
  design$active_group_list <- lapply(
    design$active_groups,
    parse_semicolon_numeric,
    integer = TRUE
  )
  design$active_sign_list <- lapply(
    design$active_signs,
    parse_semicolon_numeric,
    integer = FALSE
  )
  design$scenario_index <- seq_len(nrow(design))

  if (anyDuplicated(design$scenario) || any(design$groups <= 0L) ||
      any(design$group_size <= 0L) || any(design$cluster_size <= 1L) ||
      any(design$groups %% design$cluster_size != 0L)) {
    stop("Invalid aligned-design dimensions.", call. = FALSE)
  }
  if (any(design$rho_within <= 0 | design$rho_within >= 1) ||
      any(design$rho_between < 0) ||
      any(design$rho_between >= design$rho_within)) {
    stop("Require 0 <= rho_between < rho_within < 1.", call. = FALSE)
  }
  if (any(design$prevalence <= 0 | design$prevalence >= 1) ||
      any(design$linear_predictor_sd <= 0)) {
    stop("Invalid prevalence or signal scale.", call. = FALSE)
  }
  for (i in seq_len(nrow(design))) {
    active <- design$active_group_list[[i]]
    signs <- design$active_sign_list[[i]]
    if (!length(active) || length(active) != length(signs) ||
        anyDuplicated(active) || any(active < 1L | active > design$groups[i]) ||
        any(!signs %in% c(-1, 1))) {
      stop("Invalid active-group specification in row ", i, ".", call. = FALSE)
    }
  }
  design
}


aligned_group_covariance <- function(
    groups,
    group_size,
    cluster_size,
    rho_within,
    rho_between
) {
  p <- groups * group_size
  group <- rep(seq_len(groups), each = group_size)
  cluster <- ceiling(seq_len(groups) / cluster_size)
  covariance <- diag(p)
  for (j in seq_len(p)) {
    for (k in seq_len(j - 1L)) {
      value <- if (group[j] == group[k]) {
        rho_within
      } else if (cluster[group[j]] == cluster[group[k]]) {
        rho_between
      } else {
        0
      }
      covariance[j, k] <- value
      covariance[k, j] <- value
    }
  }
  covariance
}


aligned_population_beta <- function(scenario) {
  groups <- scenario$groups[[1L]]
  group_size <- scenario$group_size[[1L]]
  group <- rep(seq_len(groups), each = group_size)
  active <- scenario$active_group_list[[1L]]
  signs <- scenario$active_sign_list[[1L]]
  beta <- numeric(groups * group_size)
  for (i in seq_along(active)) beta[group == active[i]] <- signs[i]
  covariance <- aligned_group_covariance(
    groups,
    group_size,
    scenario$cluster_size[[1L]],
    scenario$rho_within[[1L]],
    scenario$rho_between[[1L]]
  )
  raw_sd <- sqrt(drop(crossprod(beta, covariance %*% beta)))
  beta <- beta * scenario$linear_predictor_sd[[1L]] / raw_sd
  list(beta = beta, group = group, covariance = covariance, raw_sd = raw_sd)
}


simulate_aligned_grouped_logistic <- function(scenario, seed) {
  if (nrow(scenario) != 1L) {
    stop("scenario must contain exactly one design row.", call. = FALSE)
  }
  set.seed(as.integer(seed))
  groups <- scenario$groups[[1L]]
  group_size <- scenario$group_size[[1L]]
  cluster_size <- scenario$cluster_size[[1L]]
  cluster_count <- groups %/% cluster_size
  n_total <- scenario$n_train[[1L]] +
    scenario$n_validation[[1L]] + scenario$n_test[[1L]]
  rho_within <- scenario$rho_within[[1L]]
  rho_between <- scenario$rho_between[[1L]]

  cluster_factor <- matrix(
    stats::rnorm(n_total * cluster_count),
    nrow = n_total,
    ncol = cluster_count
  )
  group_factor <- matrix(
    stats::rnorm(n_total * groups),
    nrow = n_total,
    ncol = groups
  )
  group_cluster <- ceiling(seq_len(groups) / cluster_size)
  group <- rep(seq_len(groups), each = group_size)
  X <- matrix(0, nrow = n_total, ncol = groups * group_size)
  for (g in seq_len(groups)) {
    for (j in seq_len(group_size)) {
      column <- (g - 1L) * group_size + j
      X[, column] <-
        sqrt(rho_between) * cluster_factor[, group_cluster[g]] +
        sqrt(rho_within - rho_between) * group_factor[, g] +
        sqrt(1 - rho_within) * stats::rnorm(n_total)
    }
  }
  colnames(X) <- paste0(
    "G", group, "_V", rep(seq_len(group_size), times = groups)
  )

  population <- aligned_population_beta(scenario)
  beta <- population$beta
  intercept <- calibrate_logistic_intercept(
    scenario$prevalence[[1L]],
    scenario$linear_predictor_sd[[1L]]
  )
  probability <- stats::plogis(intercept + drop(X %*% beta))
  y <- stats::rbinom(n_total, 1, probability)
  train <- seq_len(scenario$n_train[[1L]])
  validation <- max(train) + seq_len(scenario$n_validation[[1L]])
  test <- max(validation) + seq_len(scenario$n_test[[1L]])

  list(
    X_train = X[train, , drop = FALSE],
    y_train = y[train],
    X_validation = X[validation, , drop = FALSE],
    y_validation = y[validation],
    X_test = X[test, , drop = FALSE],
    y_test = y[test],
    group = group,
    beta = beta,
    intercept = intercept,
    active_groups = as.character(scenario$active_group_list[[1L]]),
    true_test_probability = probability[test],
    covariance = population$covariance,
    seed = as.integer(seed)
  )
}


target_hessian_geometry <- function(target, beta, intercept, X) {
  target <- as.numeric(target)
  beta <- as.numeric(beta)
  X <- as.matrix(X)
  if (length(target) != length(beta) || ncol(X) != length(beta) ||
      any(!is.finite(target))) {
    return(data.frame(
      target_available = FALSE,
      hessian_inner_beta_target = NA_real_,
      hessian_target_norm_squared = NA_real_,
      line_optimal_d = NA_real_,
      line_optimal_error_ratio = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  probability <- stats::plogis(intercept + drop(X %*% beta))
  weight <- probability * (1 - probability)
  beta_link <- drop(X %*% beta)
  target_link <- drop(X %*% target)
  inner <- mean(weight * beta_link * target_link)
  target_norm_squared <- mean(weight * target_link^2)
  beta_norm_squared <- mean(weight * beta_link^2)
  d_optimal <- if (target_norm_squared > 0) {
    min(max(inner / target_norm_squared, 0), 1)
  } else {
    0
  }
  error <- beta_link - d_optimal * target_link
  data.frame(
    target_available = TRUE,
    hessian_inner_beta_target = inner,
    hessian_target_norm_squared = target_norm_squared,
    line_optimal_d = d_optimal,
    line_optimal_error_ratio = sqrt(mean(weight * error^2) / beta_norm_squared),
    stringsAsFactors = FALSE
  )
}
