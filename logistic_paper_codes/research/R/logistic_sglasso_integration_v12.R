lsg_source_files_v12 <- function() {
  unique(c(v7_original$lsg_source_files_v7(), file.path("R", v12_extra),
    "R/logistic_sglasso_workflow_v12.R", "src/logistic_sglasso_hybrid_solver_v11.cpp",
    "src/logistic_sglasso_block_kernel_v11.hpp",
    "LOGISTIC_SGLASSO_INTEGRATION_PROTOCOL_V12.md",
    "scripts/52_run_logistic_sglasso_integration_v12.R",
    "tests/test_integration_v12.R", "tests/test_selected_replay_v12.R"))
}

lsg_configuration_v12 <- function(stage) {
  cfg <- v7_original$lsg_configuration_v7(stage)
  cfg$solver <- "hybrid_v11_2"
  cfg$path_controls <- lsg_hybrid_defaults_v11()$path
  cfg$polish_controls <- lsg_hybrid_defaults_v11()$polish
  cfg$stability_objective_tolerance <- 1e-7
  cfg$stability_validation_loss_tolerance <- 1e-6
  cfg$stability_prediction_tolerance <- 1e-4
  cfg$integration_version <- "v12_local_only"
  cfg
}

lsg_finite_fitters_v12 <- function(configuration) {
  # V7's exact 36-point check, reference calculation, preprocessing and
  # validation-loss array construction are reused without editing their bodies.
  fitters <- v7_original$lsg_finite_fitters_v7(configuration)
  child <- new.env(parent = environment(fitters$sglasso))
  child$fit_logistic_sglasso <- function(X, y, group, lambda, d, alpha,
      max_outer, max_inner, tolerance, inner_tolerance, target_original,
      compile, preprocess, use_active_set, warm_start_d, solver) {
    lsg_assert_v7(identical(solver, "hybrid_v11_2"), "Unexpected integration solver.")
    fit <- lsg_fit_hybrid_path_v11(X, y, group, lambda, d, alpha, target_original,
      preprocess = preprocess, controls = configuration$path_controls,
      compile = FALSE, use_active_set = use_active_set, warm_start_d = warm_start_d)
    # Legacy table column is work count, not an unchanged ABGD pass count.
    fit$passes <- fit$block_sweeps + fit$apg_iterations
    fit
  }
  environment(fitters$sglasso) <- child
  fitters
}

lsg_fit_sglasso_joint_v12 <- function(data, target_original, configuration, fitters) {
  result <- v7_original$lsg_fit_sglasso_joint_v7(data, target_original, configuration, fitters)
  result$stability <- list()
  selected <- unique(rbind(result$free_selection$row, result$d0_selection$row))
  for (i in seq_len(nrow(selected))) {
    row <- selected[i, , drop = FALSE]
    if (row$point_type == "penalty_limit") next
    fit <- result$finite$fits[[row$alpha_index]]
    li <- row$lambda_index; di <- row$d_index
    result$stability[[row$candidate_id]] <- lsg_hybrid_selected_stability_v11(
      data, fit$preprocess, target_original, row$alpha, row$d, row$lambda,
      row$lambda_relative_to_reference, li, fit$lambda,
      fit$beta_solver[, li, di], fit$intercept_solver[li, di], fit$objective[li, di],
      drop(predict_logistic_sglasso(fit, data$X_validation, lambda_index = li, d_index = di)),
      configuration)
  }
  # Do not replace only the winner by its polished fit: selection and reported
  # metrics continue to refer to the same 1e-10 path candidate.
  result
}

lsg_selected_sglasso_model_v12 <- function(grid, selection, data) {
  model <- v7_original$lsg_selected_sglasso_model_v7(grid, selection, data)
  model$stability <- grid$stability[[model$candidate_id]]
  model
}

lsg_stability_valid_v12 <- function(x, cfg) {
  if (!is.data.frame(x) || nrow(x) != 1L) return(FALSE)
  required <- c("forward_converged", "cold_converged", "reverse_converged", "forward_kkt",
    "cold_kkt", "reverse_kkt", "maximum_objective_error", "maximum_validation_loss_error",
    "maximum_validation_probability_error", "forward_polish_objective_shift",
    "forward_polish_validation_loss_shift", "forward_polish_probability_shift", "accepted")
  if (!all(required %in% names(x))) return(FALSE)
  numbers <- unlist(x[required[!required %in% c("forward_converged", "cold_converged",
    "reverse_converged", "accepted")]], use.names = FALSE)
  all(is.finite(numbers)) && all(numbers >= 0) && isTRUE(x$accepted) &&
    all(c(x$forward_converged, x$cold_converged, x$reverse_converged) %in% TRUE) &&
    max(x$forward_kkt, x$cold_kkt, x$reverse_kkt) <= cfg$polish_controls$kkt_tolerance &&
    max(x$maximum_objective_error, x$forward_polish_objective_shift) <= cfg$stability_objective_tolerance &&
    max(x$maximum_validation_loss_error, x$forward_polish_validation_loss_shift) <= cfg$stability_validation_loss_tolerance &&
    max(x$maximum_validation_probability_error, x$forward_polish_probability_shift) <= cfg$stability_prediction_tolerance
}

lsg_validate_tuning_payload_v12 <- function(payload, configuration) {
  checks <- v7_original$lsg_validate_tuning_payload_v7(payload, configuration)
  finite <- payload$sglasso_tuning[payload$sglasso_tuning$point_type == "finite", ]
  add <- function(name, value) {
    checks[nrow(checks) + 1L, ] <<- list(name, isTRUE(value), "numerical", "V12 integration")
  }
  add("all_finite_candidates_reach_requested_precision", nrow(finite) > 0L &&
    all(finite$converged %in% TRUE) && all(is.finite(finite$kkt)) &&
    max(finite$kkt) <= configuration$path_controls$kkt_tolerance)
  for (method in configuration$methods[1:2]) {
    model <- payload$artifacts$models[[method]]
    ok <- if (is.null(model)) FALSE else if (model$point_type == "penalty_limit") TRUE else {
      x <- model$stability
      lsg_stability_valid_v12(x, configuration) &&
        isTRUE(all.equal(x$selected_lambda, model$lambda, tolerance = 1e-13)) &&
        isTRUE(all.equal(x$selected_lambda_ratio, model$lambda_relative, tolerance = 1e-13))
    }
    add(paste0("selected_stability_", match(method, configuration$methods)), ok)
  }
  checks
}
