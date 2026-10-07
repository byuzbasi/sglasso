# Conditional tail diagnostic only. No package/core or frozen V5 function is
# replaced. The limit is the limit of the existing penalty, not a new target.

lsg_tail_defaults_v6 <- function() {
  list(
    source_version = "logistic_sglasso_two_design_r50_v5",
    ratios = c(512, 256, 128, 64),
    kkt_limit = 2.05e-6,
    anchor_loss_tolerance = 1e-5,
    lambda_reference_tolerance = 1e-10,
    reconstruction_tolerance = 2e-6,
    endpoint_tolerance = 1e-10,
    target_norm_tolerance = 1e-8,
    solver = "abgd", max_passes = 1000L, max_inner = 2500L,
    tolerance = 1e-6, inner_tolerance = 1e-8,
    target_max_iterations = 50L, target_tolerance = 1e-7,
    adelie_tolerance = 1e-7, adelie_max_iterations = 100000L,
    adelie_irls_tolerance = 1e-7, adelie_irls_max_iterations = 10000L,
    benchmark_nlambda = 30L, benchmark_lambda_min_ratio = 0.05
  )
}

lsg_tail_source_dependencies_v6 <- function(root, envir = .GlobalEnv) {
  options(sglasso.logistic_prework.root = normalizePath(root, mustWork = TRUE))
  for (name in c("logistic_sglasso.R", "pilot_utils.R",
                 "groupwise_logistic_targets.R", "sglasso_rhob_design.R",
                 "logistic_sglasso_two_design_v1.R")) {
    source(file.path(root, "R", name), local = envir)
  }
  invisible(TRUE)
}

lsg_tail_validate_data_v6 <- function(data) {
  required <- c("X_train", "y_train", "X_validation", "y_validation", "group")
  if (!is.list(data) || anyDuplicated(names(data)) ||
      !setequal(names(data), required)) {
    stop("Tail fits accept training/validation fields only; extra fields are forbidden.",
         call. = FALSE)
  }
  for (sample in c("train", "validation")) {
    X <- data[[paste0("X_", sample)]]
    y <- data[[paste0("y_", sample)]]
    if (!is.matrix(X) || any(!is.finite(X)) || anyNA(y) ||
        nrow(X) != length(y) || !all(y %in% c(0, 1)) ||
        length(unique(y)) != 2L || ncol(X) != length(data$group)) {
      stop("Invalid tail diagnostic training/validation data.", call. = FALSE)
    }
  }
  if (anyNA(data$group)) stop("Group labels must not be missing.", call. = FALSE)
  invisible(TRUE)
}

lsg_tail_regenerate_data_v6 <- function(scenario, seed) {
  # The frozen generator draws all rows jointly. Reducing n_test would change
  # training/validation RNG draws. Discard test/truth fields immediately.
  data <- simulate_logistic_sglasso_two_design_v1(scenario, seed)[
    c("X_train", "y_train", "X_validation", "y_validation", "group")
  ]
  lsg_tail_validate_data_v6(data)
  data
}

lsg_tail_penalty_limit_v6 <- function(preprocess, target_solver, alpha, d) {
  if (length(alpha) != 1L || length(d) != 1L ||
      !is.finite(alpha) || !is.finite(d) || alpha < 0 || alpha > 1 ||
      d < 0 || d > 1 || any(!is.finite(target_solver)) ||
      length(target_solver) != ncol(preprocess$X) ||
      any(!is.finite(preprocess$group_weight)) ||
      any(preprocess$group_weight <= 0)) {
    stop("Invalid penalty-limit inputs.", call. = FALSE)
  }
  beta <- numeric(length(target_solver))
  residual <- numeric(length(preprocess$group_start))
  for (g in seq_along(residual)) {
    index <- seq.int(preprocess$group_start[g] + 1L,
                     preprocess$group_end[g] + 1L)
    a <- d * target_solver[index]
    a_norm <- sqrt(sum(a^2))
    # The SAME positive weight multiplies both terms, so it cancels from
    # alpha/(1-alpha). No finite-lambda solver is called with Inf.
    if (alpha < 1 && a_norm > alpha / (1 - alpha)) {
      beta[index] <- (1 - alpha / ((1 - alpha) * a_norm)) * a
    }
    b_norm <- sqrt(sum(beta[index]^2))
    residual[g] <- preprocess$group_weight[g] * if (b_norm > 0) {
      sqrt(sum(((1 - alpha) * (beta[index] - a) +
                  alpha * beta[index] / b_norm)^2))
    } else {
      max(0, (1 - alpha) * a_norm - alpha)
    }
  }
  list(beta_solver = beta, penalty_kkt = max(residual))
}

lsg_tail_profile_intercept_v6 <- function(offset, y) {
  if (length(offset) != length(y) || any(!is.finite(offset)) ||
      !all(y %in% c(0, 1)) || length(unique(y)) != 2L) {
    stop("A finite offset and both training outcome classes are required.",
         call. = FALSE)
  }
  null <- stats::qlogis(mean(y))
  if (max(offset) == min(offset)) return(null - offset[1L])
  interval <- c(null - max(offset), null - min(offset))
  padding <- 1e-12 * max(1, abs(interval))
  stats::uniroot(
    function(intercept) mean(stats::plogis(intercept + offset)) - mean(y),
    interval = interval + c(-padding, padding), tol = 1e-12,
    maxiter = 1000L, check.conv = TRUE
  )$root
}

lsg_tail_point_v6 <- function(data, preprocess, beta, intercept, ratio,
                              lambda, eligible, finite_kkt = NA_real_,
                              penalty_kkt = NA_real_, configuration) {
  lsg_tail_validate_data_v6(data)
  beta <- as.numeric(beta)
  intercept <- as.numeric(intercept)
  if (any(!is.finite(c(beta, intercept)))) {
    stop("Nonfinite tail coefficients.", call. = FALSE)
  }
  X_validation <- transform_lsg_newx(preprocess, data$X_validation)
  training_probability <- stats::plogis(intercept + drop(preprocess$X %*% beta))
  validation_probability <- stats::plogis(intercept + drop(X_validation %*% beta))
  original <- recover_lsg_coefficients(preprocess, beta, intercept)
  reconstruction <- max(
    abs(training_probability - stats::plogis(
      original[1L] + drop(data$X_train %*% original[-1L]))),
    abs(validation_probability - stats::plogis(
      original[1L] + drop(data$X_validation %*% original[-1L])))
  )
  selected_groups <- sum(vapply(seq_along(preprocess$group_start), function(g) {
    index <- seq.int(preprocess$group_start[g] + 1L,
                     preprocess$group_end[g] + 1L)
    sqrt(sum(beta[index]^2)) > 1e-8
  }, logical(1)))
  point <- data.frame(
    lambda_relative = ratio, lambda = lambda,
    point_type = if (is.infinite(ratio)) "penalty_limit" else if (ratio == 64) {
      "v5_anchor_replay"
    } else "finite_extension",
    validation_log_loss = binary_log_loss(data$y_validation, validation_probability),
    training_log_loss = binary_log_loss(data$y_train, training_probability),
    selected_groups = selected_groups,
    finite_kkt = finite_kkt, penalty_kkt = penalty_kkt,
    intercept_score = abs(mean(training_probability - data$y_train)),
    reconstruction_error = reconstruction,
    numerically_eligible = isTRUE(eligible) &&
      all(is.finite(c(training_probability, validation_probability, reconstruction))) &&
      reconstruction <= configuration$reconstruction_tolerance,
    stringsAsFactors = FALSE
  )
  list(point = point, beta_solver = beta, intercept_solver = intercept,
       coefficients = original, training_probability = training_probability,
       validation_probability = validation_probability)
}

lsg_tail_limit_fit_v6 <- function(data, preprocess, target_original, alpha, d,
                                 configuration = lsg_tail_defaults_v6()) {
  lsg_tail_validate_data_v6(data)
  target <- project_lsg_original_target(preprocess, target_original)
  limit <- lsg_tail_penalty_limit_v6(preprocess, target, alpha, d)
  intercept <- lsg_tail_profile_intercept_v6(
    drop(preprocess$X %*% limit$beta_solver), data$y_train
  )
  fit <- lsg_tail_point_v6(
    data, preprocess, limit$beta_solver, intercept, Inf, Inf,
    eligible = limit$penalty_kkt <= configuration$endpoint_tolerance,
    penalty_kkt = limit$penalty_kkt, configuration = configuration
  )
  fit$point$numerically_eligible <- fit$point$numerically_eligible &&
    fit$point$intercept_score <= configuration$endpoint_tolerance
  fit
}

lsg_tail_reference_v6 <- function(preprocess, y, target_original, alpha) {
  if (alpha == 0) return(lsg_null_score_lambda_scale(preprocess, y))
  target <- project_lsg_original_target(preprocess, target_original)
  out <- lsg_lambda_start_cpp(
    preprocess$X, y, preprocess$group_start, preprocess$group_end,
    preprocess$group_weight, target, alpha, 0
  )
  if (!isTRUE(out$zero_model_feasible) || !is.finite(out$lambda_start) ||
      out$lambda_start <= 0) stop("Invalid reconstructed V5 lambda reference.",
                                call. = FALSE)
  out$lambda_start * (1 + 1e-8)
}

lsg_tail_sglasso_path_v6 <- function(data, preprocess, target_original, alpha, d,
                                    lambda_reference, configuration) {
  lsg_tail_validate_data_v6(data)
  if (!identical(configuration$ratios, c(512, 256, 128, 64))) {
    stop("The diagnostic finite ratios are frozen.", call. = FALSE)
  }
  fit <- fit_logistic_sglasso(
    data$X_train, data$y_train, data$group,
    lambda = lambda_reference * configuration$ratios,
    d = d, alpha = alpha, max_outer = configuration$max_passes,
    max_inner = configuration$max_inner, tolerance = configuration$tolerance,
    inner_tolerance = configuration$inner_tolerance,
    target_original = target_original, compile = FALSE, preprocess = preprocess,
    use_active_set = TRUE, warm_start_d = TRUE, solver = configuration$solver
  )
  rows <- lapply(seq_along(fit$lambda), function(i) {
    lsg_tail_point_v6(
      data, preprocess, fit$beta_solver[, i, 1L], fit$intercept_solver[i, 1L],
      configuration$ratios[i], fit$lambda[i],
      eligible = isTRUE(fit$converged[i, 1L]) && is.finite(fit$kkt[i, 1L]) &&
        fit$kkt[i, 1L] <= configuration$kkt_limit,
      finite_kkt = fit$kkt[i, 1L], configuration = configuration
    )$point
  })
  list(curves = do.call(rbind, rows))
}

lsg_tail_adelie_path_v6 <- function(data, preprocess, alpha, lambda_reference,
                                   configuration) {
  lsg_tail_validate_data_v6(data)
  if (!identical(configuration$ratios, c(512, 256, 128, 64)) || alpha != 0) {
    stop("V6 Adelie diagnostics are restricted to the frozen alpha-zero cases.",
         call. = FALSE)
  }
  lambda <- lambda_reference * configuration$ratios
  warnings <- character()
  fit <- withCallingHandlers(adelie::grpnet(
    X = preprocess$X, glm = adelie::glm.binomial(data$y_train),
    groups = as.integer(preprocess$group_start + 1L), alpha = alpha,
    penalty = preprocess$group_weight, standardize = FALSE, intercept = TRUE,
    lmda_path_size = configuration$benchmark_nlambda,
    min_ratio = configuration$benchmark_lambda_min_ratio,
    tol = configuration$adelie_tolerance,
    max_iters = configuration$adelie_max_iterations,
    irls_tol = configuration$adelie_irls_tolerance,
    irls_max_iters = configuration$adelie_irls_max_iterations,
    screen_rule = "strong", early_exit = FALSE, check_state = TRUE,
    progress_bar = FALSE, n_threads = 1L, lambda = lambda
  ), warning = function(condition) {
    warnings <<- c(warnings, conditionMessage(condition))
    invokeRestart("muffleWarning")
  })
  coefficient_method <- getS3method("coef", "grpnet", envir = asNamespace("adelie"))
  path <- coefficient_method(fit)
  beta <- as.matrix(path$betas)
  intercept <- as.numeric(path$intercepts[, 1L])
  returned <- as.numeric(path$lambda)
  state_lambda <- as.numeric(fit$state$lmdas)
  path_ok <- length(returned) == length(lambda) &&
    length(state_lambda) == length(lambda) &&
    nrow(beta) == length(lambda) && ncol(beta) == ncol(preprocess$X) &&
    length(intercept) == length(lambda) &&
    max(abs(returned - lambda)) <= 1e-10 * max(1, abs(lambda)) &&
    max(abs(state_lambda - lambda)) <= 1e-10 * max(1, abs(lambda))
  if (!path_ok) stop("Adelie diagnostic path is incomplete or misaligned.",
                     call. = FALSE)
  rows <- lapply(seq_along(lambda), function(i) {
    # Record a common-objective KKT, but preserve V5 Adelie eligibility:
    # complete warning-free path and finite predictions/coefficients.
    kkt <- lsg_kkt_cpp(
      preprocess$X, data$y_train, beta[i, ], intercept[i],
      preprocess$group_start, preprocess$group_end, preprocess$group_weight,
      rep(0, ncol(preprocess$X)), lambda[i], 0, 0
    )$maximum
    lsg_tail_point_v6(
      data, preprocess, beta[i, ], intercept[i], configuration$ratios[i], lambda[i],
      eligible = !length(warnings), finite_kkt = kkt,
      configuration = configuration
    )$point
  })
  list(curves = do.call(rbind, rows), solver_warnings = unique(warnings))
}

lsg_tail_attach_case_v6 <- function(curves, case) {
  if (nrow(case) != 1L) stop("One fixed case is required.", call. = FALSE)
  identity <- case[rep(1L, nrow(curves)),
                   c("case_id", "scenario", "replication", "method", "alpha", "d"),
                   drop = FALSE]
  rownames(identity) <- NULL
  cbind(identity, curves)
}

lsg_tail_run_task_v6 <- function(task, scenario, cases, source_shard, configuration) {
  if (nrow(task) != 1L || !nrow(cases) ||
      any(cases$task_id != task$task_id) || any(cases$seed != task$seed) ||
      any(cases$scenario != task$scenario) ||
      any(cases$replication != task$replication)) {
    stop("The diagnostic task and frozen cases differ.", call. = FALSE)
  }
  data <- lsg_tail_regenerate_data_v6(scenario, task$seed)
  fingerprint <- digest::digest(data, algo = "sha256", serialize = TRUE)
  firth <- estimate_groupwise_logistic_target(
    data$X_train, data$y_train, data$group, method = "firth",
    max_iterations = configuration$target_max_iterations,
    tolerance = configuration$target_tolerance
  )
  if (!isTRUE(firth$success) || firth$failed_groups != 0L) {
    stop("The training-only Firth target failed; no fallback is permitted.",
         call. = FALSE)
  }
  source_target <- source_shard$payload$targets
  source_target <- source_target[source_target$target_mode == "groupwise_firth", ]
  if (nrow(source_target) != 1L || !isTRUE(source_target$target_available) ||
      source_target$failed_groups != 0L || !is.finite(source_target$target_l2_norm)) {
    stop("The saved V5 Firth target audit is unavailable.", call. = FALSE)
  }
  target_norm <- sqrt(sum(firth$target_original^2))
  target_error <- abs(target_norm - source_target$target_l2_norm) /
    max(1, abs(source_target$target_l2_norm))
  preprocess <- prepare_lsg_design(data$X_train, data$group)
  curve_rows <- anchor_rows <- vector("list", nrow(cases))
  for (i in seq_len(nrow(cases))) {
    case <- cases[i, , drop = FALSE]
    adelie <- identical(case$method, "Logistic Group Elastic Net (adelie)")
    effective_d <- if (adelie) 0 else case$d
    reference <- lsg_tail_reference_v6(
      preprocess, data$y_train, firth$target_original, case$alpha
    )
    reference_error <- abs(reference - case$lambda_reference) /
      max(1, abs(case$lambda_reference))
    finite <- if (adelie) {
      lsg_tail_adelie_path_v6(data, preprocess, case$alpha, reference, configuration)
    } else {
      lsg_tail_sglasso_path_v6(
        data, preprocess, firth$target_original, case$alpha, effective_d,
        reference, configuration
      )
    }
    endpoint <- lsg_tail_limit_fit_v6(
      data, preprocess, firth$target_original, case$alpha, effective_d, configuration
    )
    curve_rows[[i]] <- lsg_tail_attach_case_v6(
      rbind(finite$curves, endpoint$point), case
    )
    replay <- finite$curves$validation_log_loss[finite$curves$lambda_relative == 64]
    error <- abs(replay - case$validation_log_loss_v5)
    anchor_rows[[i]] <- data.frame(
      case_id = case$case_id,
      lambda_reference_stored = case$lambda_reference,
      lambda_reference_recomputed = reference,
      lambda_reference_error = reference_error,
      validation_log_loss_v5 = case$validation_log_loss_v5,
      validation_log_loss_replay = replay, anchor_loss_error = error,
      anchor_within_tolerance = error <= configuration$anchor_loss_tolerance
    )
  }
  curves <- do.call(rbind, curve_rows)
  anchors <- do.call(rbind, anchor_rows)
  limits <- is.infinite(curves$lambda_relative)
  target_diagnostics <- firth$diagnostics
  target_diagnostics$target_norm_stored <- source_target$target_l2_norm
  target_diagnostics$target_norm_recomputed <- target_norm
  target_diagnostics$target_norm_relative_error <- target_error
  checks <- data.frame(
    check = c("fixed_cases_preserved", "five_points_per_case",
              "finite_points_eligible", "limit_points_eligible",
              "lambda_reference_reproduced", "v5_anchor_reproduced",
              "target_reproduced", "training_validation_only",
              "coefficient_reconstruction"),
    passed = c(
      setequal(unique(curves$case_id), cases$case_id),
      nrow(curves) == 5L * nrow(cases) && all(table(curves$case_id) == 5L),
      all(curves$numerically_eligible[!limits]),
      all(curves$numerically_eligible[limits]),
      all(anchors$lambda_reference_error <= configuration$lambda_reference_tolerance),
      all(anchors$anchor_within_tolerance),
      target_error <= configuration$target_norm_tolerance,
      setequal(names(data), c("X_train", "y_train", "X_validation", "y_validation", "group")),
      all(curves$reconstruction_error <= configuration$reconstruction_tolerance)
    ), stringsAsFactors = FALSE
  )
  list(curves = curves, anchors = anchors, target_diagnostics = target_diagnostics,
       training_validation_fingerprint = fingerprint, checks = checks, test_used = FALSE)
}
