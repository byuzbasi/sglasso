# Frozen two-design generator for the independent Logistic SGLASSO study.

read_logistic_sglasso_two_design_v1 <- function(path) {
  design <- utils::read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE
  )
  required <- c(
    "design_id", "scenario", "n_train", "n_validation", "n_test",
    "groups", "group_size", "active_group_count", "rho_within",
    "rho_between", "prevalence", "linear_predictor_sd",
    "signal_pattern", "description"
  )
  if (!identical(names(design), required)) {
    stop("Two-design columns do not match the frozen schema.", call. = FALSE)
  }
  integer_columns <- c(
    "n_train", "n_validation", "n_test", "groups", "group_size",
    "active_group_count"
  )
  numeric_columns <- c(
    "rho_within", "rho_between", "prevalence", "linear_predictor_sd"
  )
  for (name in integer_columns) design[[name]] <- as.integer(design[[name]])
  for (name in numeric_columns) design[[name]] <- as.numeric(design[[name]])
  design$scenario_index <- seq_len(nrow(design))

  numeric_matrix <- as.matrix(design[c(integer_columns, numeric_columns)])
  allowed_pattern <- c("homogeneous", "weak_strong_mixed")
  if (nrow(design) != 4L || anyDuplicated(design$scenario) ||
      any(!is.finite(numeric_matrix)) ||
      any(design$n_train < 20L | design$n_validation < 20L |
            design$n_test < 20L) ||
      any(design$groups < 2L | design$group_size < 1L) ||
      any(design$active_group_count < 1L |
            design$active_group_count >= design$groups) ||
      any(design$rho_within <= 0 | design$rho_within >= 1) ||
      any(design$rho_between < 0 |
            design$rho_between > design$rho_within) ||
      any(design$prevalence <= 0 | design$prevalence >= 1) ||
      any(design$linear_predictor_sd <= 0) ||
      any(!design$signal_pattern %in% allowed_pattern)) {
    stop("Invalid values in the frozen two-design file.", call. = FALSE)
  }
  expected_design <- c("homogeneous", "weak_strong_mixed")
  if (!setequal(unique(design$design_id), expected_design) ||
      any(vapply(split(design$rho_between, design$design_id), function(x) {
        !isTRUE(all.equal(sort(x), c(0.3, 0.9)))
      }, logical(1))) ||
      any(design$design_id != design$signal_pattern) ||
      any(design$n_train != 200L) || any(design$n_validation != 200L) ||
      any(design$n_test != 5000L) || any(design$groups != 200L) ||
      any(design$group_size != 3L) ||
      any(design$active_group_count != 10L) ||
      any(abs(design$rho_within - 0.9) > 1e-12) ||
      any(abs(design$prevalence - 0.5) > 1e-12) ||
      any(abs(design$linear_predictor_sd - 1.1) > 1e-12)) {
    stop("The two-design scientific specification has changed.",
         call. = FALSE)
  }
  design
}


logistic_sglasso_group_signal_values_v1 <- function(
    active_group_count,
    signal_pattern
) {
  active_group_count <- as.integer(active_group_count)
  if (signal_pattern == "homogeneous") {
    return(rep(1, active_group_count))
  }
  if (signal_pattern != "weak_strong_mixed") {
    stop("Unknown signal pattern: ", signal_pattern, call. = FALSE)
  }
  n_strong <- ceiling(active_group_count / 2)
  n_weak <- floor(active_group_count / 2)
  strong <- stats::runif(n_strong, min = 0.8, max = 1.2)
  weak <- stats::runif(n_weak, min = 0.15, max = 0.35)
  paired <- min(length(strong), length(weak))
  values <- as.vector(rbind(strong[seq_len(paired)], weak[seq_len(paired)]))
  if (length(values) < active_group_count) {
    values <- c(values, strong[length(strong)])
  }
  values <- values[seq_len(active_group_count)]
  values * rep(c(1, -1), length.out = active_group_count)
}


simulate_logistic_sglasso_two_design_v1 <- function(scenario, seed) {
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

  active_groups <- sort(sample.int(
    groups, scenario$active_group_count[[1L]]
  ))
  group_values <- logistic_sglasso_group_signal_values_v1(
    length(active_groups), scenario$signal_pattern[[1L]]
  )
  beta <- numeric(p)
  for (index in seq_along(active_groups)) {
    beta[group == active_groups[index]] <- group_values[index]
  }
  raw_sd <- sqrt(sglasso_exchangeable_signal_variance(
    beta, group, rho_within, rho_between
  ))
  if (!is.finite(raw_sd) || raw_sd <= 0) {
    stop("Generated coefficient vector has invalid signal variance.",
         call. = FALSE)
  }
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
    seed = as.integer(seed),
    design_id = scenario$design_id[[1L]],
    signal_pattern = scenario$signal_pattern[[1L]]
  )
}
