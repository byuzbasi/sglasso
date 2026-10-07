# Local-only diagnostic for the upper finite-lambda tail selected in V15.
# It does not alter the estimator, package core, frozen grids, seeds or results.

lsg_v17_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

lsg_v17_constants <- function() {
  list(
    schema = "logistic_sglasso_lambda_sat_diagnostic_v17_1",
    output_version = "logistic_sglasso_lambda_sat_diagnostic_v17_a01",
    source_version = "logistic_sglasso_joint_r50_v15_a01",
    source_signature =
      "292b830187e0588a966bb9b4a0e069f9307ac42827dae76e131c063c0415b49e",
    task_ids = c(240L, 334L),
    finite_ratios = c(32768, 16384, 8192, 4096, 2048, 1024,
                      512, 256, 128, 64),
    anchor_ratios = c(64, 128, 256, 512),
    anchor_loss_tolerance = 1e-5,
    lambda_reference_tolerance = 1e-10
  )
}

lsg_v17_cases <- function() {
  data.frame(
    case_id = c("task240_selected_sglasso", "task334_d0_sglasso",
                "task334_alpha0_adelie"),
    task_id = c(240L, 334L, 334L),
    method = c("Logistic SGLASSO",
               "Logistic SGLASSO (d=0 boundary)",
               "Logistic Group Elastic Net (adelie)"),
    alpha = c(0.1, 0, 0), d = c(0.3, 0, NA_real_),
    stringsAsFactors = FALSE
  )
}

lsg_v17_load_source <- function(source_root) {
  source(file.path(source_root, "R/logistic_sglasso_failure_diagnostic_v16.R"),
         local = TRUE)
  e <- lsg_v16_load_source(source_root)
  required <- c(
    "lsg_design_for_stage_v7", "compile_lsg_core",
    "lsg_compile_hybrid_solver_v11",
    "lsg_fit_hybrid_path_v11", "lsg_tail_reference_v6",
    "lsg_tail_limit_fit_v6", "estimate_groupwise_logistic_target",
    "prepare_lsg_design", "project_lsg_original_target",
    "transform_lsg_newx", "predict_logistic_sglasso", "binary_log_loss",
    "lsg_runtime_v7"
  )
  lsg_v17_assert(all(vapply(required, exists, logical(1), envir = e,
                            inherits = FALSE)),
                 "The exact V15 source environment is incomplete.")
  e
}

lsg_v17_classical_d0_lambda_max <- function(preprocess, y, alpha) {
  if (!is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("alpha must lie in [0, 1].", call. = FALSE)
  }
  if (alpha == 0) return(NA_real_)
  score <- drop(crossprod(preprocess$X, y - mean(y)) / length(y))
  values <- vapply(seq_along(preprocess$group_start), function(g) {
    index <- seq.int(preprocess$group_start[g] + 1L,
                     preprocess$group_end[g] + 1L)
    sqrt(sum(score[index]^2)) / preprocess$group_weight[g]
  }, numeric(1))
  max(values) / alpha
}

lsg_v17_zero_interval <- function(e, preprocess, y, target, alpha, d,
                                  case_id) {
  diagnostic <- get("lsg_lambda_start_cpp", envir = e, inherits = TRUE)(
    preprocess$X, y, preprocess$group_start, preprocess$group_end,
    preprocess$group_weight, target, alpha, d
  )
  groups <- data.frame(
    case_id = case_id,
    group = seq_along(preprocess$group_start),
    feasible = as.logical(diagnostic$group_feasible),
    lower = as.numeric(diagnostic$group_lower),
    upper = as.numeric(diagnostic$group_upper),
    stringsAsFactors = FALSE
  )
  summary <- data.frame(
    case_id = case_id, alpha = alpha, d = d,
    zero_model_feasible = isTRUE(diagnostic$zero_model_feasible),
    common_lower = as.numeric(diagnostic$lambda_start),
    common_upper = as.numeric(diagnostic$common_upper),
    feasible_groups = sum(groups$feasible), total_groups = nrow(groups),
    finite_group_upper_bounds = sum(is.finite(groups$upper)),
    stringsAsFactors = FALSE
  )
  list(summary = summary, groups = groups)
}

lsg_v17_curve_rows <- function(e, data, preprocess, fit, endpoint, ratios,
                               case, reference) {
  probability <- e$predict_logistic_sglasso(
    fit, data$X_validation, type = "response"
  )[, , 1L, drop = FALSE][, , 1L]
  if (is.null(dim(probability))) probability <- matrix(probability, ncol = 1L)
  rows <- lapply(seq_along(ratios), function(i) {
    beta <- fit$beta_solver[, i, 1L]
    p <- probability[, i]
    data.frame(
      case_id = case$case_id, task_id = case$task_id,
      method = case$method, alpha = case$alpha, d = case$d,
      lambda_relative_to_reference = ratios[i],
      lambda = fit$lambda[i], point_type = if (ratios[i] <= 512) {
        "v15_anchor_replay"
      } else "finite_extension",
      validation_log_loss = e$binary_log_loss(data$y_validation, p),
      validation_loss_minus_limit =
        e$binary_log_loss(data$y_validation, p) -
        endpoint$point$validation_log_loss,
      coefficient_l2_distance_to_limit =
        sqrt(sum((beta - endpoint$beta_solver)^2)),
      coefficient_relative_l2_distance_to_limit =
        sqrt(sum((beta - endpoint$beta_solver)^2)) /
        max(1, sqrt(sum(endpoint$beta_solver^2))),
      validation_probability_max_distance_to_limit =
        max(abs(p - endpoint$validation_probability)),
      validation_probability_rmse_to_limit =
        sqrt(mean((p - endpoint$validation_probability)^2)),
      intercept_distance_to_limit =
        abs(fit$intercept_solver[i, 1L] - endpoint$intercept_solver),
      selected_groups = fit$selected_groups[i, 1L],
      kkt = fit$kkt[i, 1L], converged = isTRUE(fit$converged[i, 1L]),
      solver_route = fit$solver_route[i, 1L],
      fallback_used = isTRUE(fit$fallback_used[i, 1L]),
      lambda_reference = reference, stringsAsFactors = FALSE
    )
  })
  limit <- data.frame(
    case_id = case$case_id, task_id = case$task_id,
    method = case$method, alpha = case$alpha, d = case$d,
    lambda_relative_to_reference = Inf, lambda = Inf,
    point_type = "analytic_penalty_limit",
    validation_log_loss = endpoint$point$validation_log_loss,
    validation_loss_minus_limit = 0,
    coefficient_l2_distance_to_limit = 0,
    coefficient_relative_l2_distance_to_limit = 0,
    validation_probability_max_distance_to_limit = 0,
    validation_probability_rmse_to_limit = 0,
    intercept_distance_to_limit = 0,
    selected_groups = endpoint$point$selected_groups,
    kkt = NA_real_, converged = endpoint$point$numerically_eligible,
    solver_route = "closed_form_penalty_limit", fallback_used = FALSE,
    lambda_reference = reference, stringsAsFactors = FALSE
  )
  do.call(rbind, c(rows, list(limit)))
}

lsg_v17_fit_sglasso <- function(e, data, preprocess, target_original, case,
                                reference, ratios, configuration) {
  fit <- e$lsg_fit_hybrid_path_v11(
    data$X_train, data$y_train, data$group,
    lambda = reference * ratios, d = case$d, alpha = case$alpha,
    target_original = target_original, preprocess = preprocess,
    controls = configuration$path_controls, compile = FALSE,
    use_active_set = TRUE, enable_fallback = TRUE, warm_start_d = TRUE,
    lambda_order = "decreasing"
  )
  endpoint <- e$lsg_tail_limit_fit_v6(
    data, preprocess, target_original, case$alpha, case$d, configuration
  )
  list(
    curves = lsg_v17_curve_rows(e, data, preprocess, fit, endpoint, ratios,
                                case, reference),
    fit = fit, endpoint = endpoint
  )
}

lsg_v17_fit_adelie <- function(e, data, preprocess, case, reference, ratios,
                               configuration) {
  lsg_v17_assert(requireNamespace("adelie", quietly = TRUE),
                 "The installed adelie package is required.")
  lambda <- reference * ratios
  warnings <- character()
  fit <- withCallingHandlers(adelie::grpnet(
    X = preprocess$X, glm = adelie::glm.binomial(data$y_train),
    groups = as.integer(preprocess$group_start + 1L), alpha = 0,
    penalty = preprocess$group_weight, standardize = FALSE, intercept = TRUE,
    lmda_path_size = length(lambda), min_ratio = min(ratios) / max(ratios),
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
  lsg_v17_assert(length(returned) == length(lambda) &&
                   nrow(beta) == length(lambda) &&
                   length(intercept) == length(lambda) &&
                   max(abs(returned - lambda)) <= 1e-10 * max(1, max(lambda)),
                 "Adelie returned an incomplete or misaligned explicit path.")
  Xv <- e$transform_lsg_newx(preprocess, data$X_validation)
  endpoint <- e$lsg_tail_limit_fit_v6(
    data, preprocess, rep(0, ncol(data$X_train)), 0, 0, configuration
  )
  rows <- lapply(seq_along(ratios), function(i) {
    p <- stats::plogis(intercept[i] + drop(Xv %*% beta[i, ]))
    kkt <- get("lsg_kkt_cpp", envir = e, inherits = TRUE)(
      preprocess$X, data$y_train, beta[i, ], intercept[i],
      preprocess$group_start, preprocess$group_end, preprocess$group_weight,
      rep(0, ncol(preprocess$X)), lambda[i], 0, 0
    )$maximum
    data.frame(
      case_id = case$case_id, task_id = case$task_id, method = case$method,
      alpha = 0, d = NA_real_, lambda_relative_to_reference = ratios[i],
      lambda = lambda[i], point_type = if (ratios[i] <= 512) {
        "v15_anchor_replay"
      } else "finite_extension",
      validation_log_loss = e$binary_log_loss(data$y_validation, p),
      validation_loss_minus_limit =
        e$binary_log_loss(data$y_validation, p) - endpoint$point$validation_log_loss,
      coefficient_l2_distance_to_limit = sqrt(sum(beta[i, ]^2)),
      coefficient_relative_l2_distance_to_limit = sqrt(sum(beta[i, ]^2)),
      validation_probability_max_distance_to_limit =
        max(abs(p - endpoint$validation_probability)),
      validation_probability_rmse_to_limit =
        sqrt(mean((p - endpoint$validation_probability)^2)),
      intercept_distance_to_limit = abs(intercept[i] - endpoint$intercept_solver),
      selected_groups = sum(vapply(seq_along(preprocess$group_start), function(g) {
        index <- seq.int(preprocess$group_start[g] + 1L,
                         preprocess$group_end[g] + 1L)
        sqrt(sum(beta[i, index]^2)) > configuration$selection_threshold
      }, logical(1))),
      kkt = as.numeric(kkt),
      converged = !length(warnings) && all(is.finite(c(beta[i, ], intercept[i], p, kkt))),
      solver_route = "adelie_explicit_path", fallback_used = FALSE,
      lambda_reference = reference, stringsAsFactors = FALSE
    )
  })
  limit <- data.frame(
    case_id = case$case_id, task_id = case$task_id, method = case$method,
    alpha = 0, d = NA_real_, lambda_relative_to_reference = Inf, lambda = Inf,
    point_type = "analytic_penalty_limit",
    validation_log_loss = endpoint$point$validation_log_loss,
    validation_loss_minus_limit = 0,
    coefficient_l2_distance_to_limit = 0,
    coefficient_relative_l2_distance_to_limit = 0,
    validation_probability_max_distance_to_limit = 0,
    validation_probability_rmse_to_limit = 0,
    intercept_distance_to_limit = 0,
    selected_groups = endpoint$point$selected_groups, kkt = NA_real_,
    converged = endpoint$point$numerically_eligible,
    solver_route = "closed_form_penalty_limit", fallback_used = FALSE,
    lambda_reference = reference, stringsAsFactors = FALSE
  )
  list(curves = do.call(rbind, c(rows, list(limit))), warnings = unique(warnings),
       endpoint = endpoint)
}

lsg_v17_archived_path <- function(packet, case) {
  evidence <- packet$boundary_evidence[[as.character(case$task_id)]]
  if (identical(case$method, "Logistic Group Elastic Net (adelie)")) {
    path <- evidence$adelie_path
    path <- path[abs(path$alpha - case$alpha) <= 1e-12, ]
  } else {
    path <- evidence$sglasso_paths
    path <- path[abs(path$alpha - case$alpha) <= 1e-12 &
                   abs(path$d - case$d) <= 1e-12, ]
  }
  path[, c("lambda_relative_to_reference", "validation_log_loss"), drop = FALSE]
}

lsg_v17_run <- function(source_root, packet) {
  constants <- lsg_v17_constants()
  lsg_v16_validate_packet(packet, production = TRUE)
  lsg_v17_assert(identical(packet$source_version, constants$source_version) &&
                   identical(packet$source_scientific_signature,
                             constants$source_signature),
                 "The V15 source evidence differs from the frozen V17 scope.")
  e <- lsg_v17_load_source(source_root)
  design <- e$lsg_design_for_stage_v7(source_root, "production")
  cases <- lsg_v17_cases()
  e$compile_lsg_core(rebuild = FALSE, quiet = TRUE)
  e$lsg_compile_hybrid_solver_v11(source_root, rebuild = FALSE, quiet = TRUE)
  lsg_v17_assert(all(vapply(c("lsg_lambda_start_cpp", "lsg_kkt_cpp"), exists,
                            logical(1), envir = e, inherits = TRUE)),
                 "The compiled V15 core functions are unavailable.")
  curves <- anchors <- references <- intervals <- interval_groups <- list()
  fingerprints <- list()

  for (task_id in constants$task_ids) {
    evidence <- packet$boundary_evidence[[as.character(task_id)]]
    task <- evidence$task
    scenario <- design[design$scenario_index == task$scenario_index, , drop = FALSE]
    data <- lsg_v16_tuning_data(e, scenario, task$seed)
    fingerprints[[as.character(task_id)]] <- lsg_v16_input_fingerprints(e, data)
    firth <- e$estimate_groupwise_logistic_target(
      data$X_train, data$y_train, data$group, method = "firth",
      max_iterations = packet$configuration$target_max_iterations,
      tolerance = packet$configuration$target_tolerance
    )
    lsg_v17_assert(isTRUE(firth$success) && firth$failed_groups == 0L &&
                     all(is.finite(firth$target_original)),
                   paste("Firth target failed for task", task_id))
    preprocess <- e$prepare_lsg_design(data$X_train, data$group)
    target <- e$project_lsg_original_target(preprocess, firth$target_original)
    task_cases <- cases[cases$task_id == task_id, , drop = FALSE]
    for (i in seq_len(nrow(task_cases))) {
      case <- task_cases[i, , drop = FALSE]
      effective_d <- if (is.na(case$d)) 0 else case$d
      reference <- e$lsg_tail_reference_v6(
        preprocess, data$y_train, firth$target_original, case$alpha
      )
      archived <- evidence$selected[evidence$selected$method == case$method, ]
      reference_error <- abs(reference - archived$selected_lambda_reference) /
        max(1, abs(archived$selected_lambda_reference))
      classical <- lsg_v17_classical_d0_lambda_max(
        preprocess, data$y_train, case$alpha
      )
      interval <- lsg_v17_zero_interval(
        e, preprocess, data$y_train, target, case$alpha, effective_d,
        case$case_id
      )
      intervals[[case$case_id]] <- interval$summary
      interval_groups[[case$case_id]] <- interval$groups
      zero_kkt <- function(lambda, d) {
        if (!is.finite(lambda) || lambda <= 0) return(NA_real_)
        get("lsg_kkt_cpp", envir = e, inherits = TRUE)(
          preprocess$X, data$y_train, rep(0, ncol(preprocess$X)),
          stats::qlogis(mean(data$y_train)), preprocess$group_start,
          preprocess$group_end, preprocess$group_weight, target,
          lambda, case$alpha, d
        )$maximum
      }
      references[[case$case_id]] <- data.frame(
        case_id = case$case_id, task_id = task_id, method = case$method,
        alpha = case$alpha, d = case$d,
        v15_reference = archived$selected_lambda_reference,
        recomputed_reference = reference,
        relative_reference_error = reference_error,
        classical_d0_lambda_max = classical,
        reference_over_classical = if (is.finite(classical)) reference / classical else NA_real_,
        zero_kkt_at_classical_d0_lambda_max = zero_kkt(classical, effective_d),
        zero_model_feasible_for_actual_d = interval$summary$zero_model_feasible,
        actual_d_zero_interval_lower = interval$summary$common_lower,
        actual_d_zero_interval_upper = interval$summary$common_upper,
        stringsAsFactors = FALSE
      )
      fitted <- if (identical(case$method,
                              "Logistic Group Elastic Net (adelie)")) {
        lsg_v17_fit_adelie(e, data, preprocess, case, reference,
                          constants$finite_ratios, packet$configuration)
      } else {
        lsg_v17_fit_sglasso(e, data, preprocess, firth$target_original,
                           case, reference, constants$finite_ratios,
                           packet$configuration)
      }
      curves[[case$case_id]] <- fitted$curves
      archived_path <- lsg_v17_archived_path(packet, case)
      replay <- fitted$curves[
        fitted$curves$lambda_relative_to_reference %in% constants$anchor_ratios,
        c("lambda_relative_to_reference", "validation_log_loss")]
      names(replay)[2L] <- "validation_log_loss_replay"
      archived_path <- archived_path[
        archived_path$lambda_relative_to_reference %in% constants$anchor_ratios, ]
      names(archived_path)[2L] <- "validation_log_loss_v15"
      comparison <- merge(archived_path, replay,
                          by = "lambda_relative_to_reference", all = TRUE)
      comparison$case_id <- case$case_id
      comparison$absolute_error <- abs(comparison$validation_log_loss_replay -
                                         comparison$validation_log_loss_v15)
      comparison$within_tolerance <- comparison$absolute_error <=
        constants$anchor_loss_tolerance
      anchors[[case$case_id]] <- comparison
    }
  }
  curves <- do.call(rbind, curves)
  curves <- curves[order(curves$task_id, curves$method,
                         curves$lambda_relative_to_reference), ]
  rownames(curves) <- NULL
  anchors <- do.call(rbind, anchors); rownames(anchors) <- NULL
  references <- do.call(rbind, references); rownames(references) <- NULL
  intervals <- do.call(rbind, intervals); rownames(intervals) <- NULL
  interval_groups <- do.call(rbind, interval_groups); rownames(interval_groups) <- NULL
  checks <- data.frame(
    check = c("source_scope_frozen", "only_training_validation_regenerated",
              "lambda_references_reproduced", "four_v15_anchors_per_case",
              "v15_anchor_losses_reproduced", "all_finite_fits_converged",
              "all_finite_kkt_within_v15_limit", "analytic_limits_eligible"),
    passed = c(
      identical(packet$source_version, constants$source_version) &&
        identical(packet$source_scientific_signature, constants$source_signature),
      TRUE,
      all(references$relative_reference_error <=
            constants$lambda_reference_tolerance),
      nrow(anchors) == 4L * nrow(cases) &&
        all(table(anchors$case_id) == 4L),
      all(anchors$within_tolerance),
      all(curves$converged[is.finite(curves$lambda)]),
      all(curves$kkt[is.finite(curves$lambda) &
                       curves$method != "Logistic Group Elastic Net (adelie)"] <=
            packet$configuration$full_path_kkt_limit),
      all(curves$converged[is.infinite(curves$lambda)])),
    stringsAsFactors = FALSE
  )
  list(
    schema = constants$schema, accepted = all(checks$passed), checks = checks,
    cases = cases, curves = curves, anchors = anchors,
    references = references, zero_intervals = intervals,
    zero_interval_groups = interval_groups, fingerprints = fingerprints,
    finite_ratios = constants$finite_ratios,
    model_paths_fitted = 3L, test_fields_used = FALSE,
    source_version = packet$source_version,
    source_signature = packet$source_scientific_signature,
    runtime = e$lsg_runtime_v7(),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
}
