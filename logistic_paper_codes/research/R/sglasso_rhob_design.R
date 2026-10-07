# SGLASSO-style exchangeable grouped-logistic design.

read_sglasso_rhob_design <- function(path) {
  design <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c(
    "scenario", "n_train", "n_validation", "n_test", "groups",
    "group_size", "active_group_count", "rho_within", "rho_between",
    "prevalence", "linear_predictor_sd", "signal_pattern", "description"
  )
  if (!identical(names(design), required)) {
    stop("SGLASSO rho_b design columns do not match the frozen schema.",
         call. = FALSE)
  }
  integer_columns <- c(
    "n_train", "n_validation", "n_test", "groups", "group_size",
    "active_group_count"
  )
  for (name in integer_columns) design[[name]] <- as.integer(design[[name]])
  design$scenario_index <- seq_len(nrow(design))
  if (anyDuplicated(design$scenario) ||
      any(design$n_train < 2L | design$n_validation < 2L | design$n_test < 2L) ||
      any(design$groups < 2L | design$group_size < 1L) ||
      any(design$active_group_count < 1L |
          design$active_group_count >= design$groups) ||
      any(design$rho_within <= 0 | design$rho_within >= 1) ||
      any(design$rho_between < 0 | design$rho_between > design$rho_within) ||
      any(design$prevalence <= 0 | design$prevalence >= 1) ||
      any(design$linear_predictor_sd <= 0) ||
      any(design$signal_pattern != "homogeneous")) {
    stop("Invalid SGLASSO rho_b design values.", call. = FALSE)
  }
  expected_rho_between <- seq(0, 0.9, by = 0.1)
  if (nrow(design) != length(expected_rho_between) ||
      !isTRUE(all.equal(sort(design$rho_between), expected_rho_between))) {
    stop("The frozen design must contain rho_between = 0, 0.1, ..., 0.9.",
         call. = FALSE)
  }
  design
}


sglasso_exchangeable_signal_variance <- function(beta, group, rho_within,
                                                  rho_between) {
  group_sum <- as.numeric(rowsum(beta, group, reorder = FALSE))
  (1 - rho_within) * sum(beta^2) +
    (rho_within - rho_between) * sum(group_sum^2) +
    rho_between * sum(beta)^2
}


simulate_sglasso_rhob_logistic <- function(scenario, seed) {
  if (nrow(scenario) != 1L) {
    stop("scenario must contain exactly one design row.", call. = FALSE)
  }
  set.seed(as.integer(seed))
  n_total <- scenario$n_train[[1L]] + scenario$n_validation[[1L]] +
    scenario$n_test[[1L]]
  groups <- scenario$groups[[1L]]
  group_size <- scenario$group_size[[1L]]
  p <- groups * group_size
  rho_within <- scenario$rho_within[[1L]]
  rho_between <- scenario$rho_between[[1L]]
  group <- rep(seq_len(groups), each = group_size)

  common_factor <- stats::rnorm(n_total)
  group_factor <- matrix(stats::rnorm(n_total * groups), nrow = n_total)
  noise <- matrix(stats::rnorm(n_total * p), nrow = n_total)
  X <- sqrt(rho_between) * common_factor +
    sqrt(rho_within - rho_between) *
      group_factor[, group, drop = FALSE] +
    sqrt(1 - rho_within) * noise
  colnames(X) <- paste0(
    "G", group, "_V", rep(seq_len(group_size), times = groups)
  )

  active_groups <- sort(sample.int(groups, scenario$active_group_count[[1L]]))
  beta <- numeric(p)
  beta[group %in% active_groups] <- 1
  raw_sd <- sqrt(sglasso_exchangeable_signal_variance(
    beta, group, rho_within, rho_between
  ))
  beta <- beta * scenario$linear_predictor_sd[[1L]] / raw_sd
  intercept <- calibrate_logistic_intercept(
    scenario$prevalence[[1L]], scenario$linear_predictor_sd[[1L]]
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
    active_groups = as.character(active_groups),
    true_test_probability = probability[test],
    seed = as.integer(seed)
  )
}
