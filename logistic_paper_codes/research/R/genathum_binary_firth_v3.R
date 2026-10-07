# Training-only Firth target for the GenAtHum application. No frozen simulation
# fitter or package API is changed. See GENATHUM_BINARY_PROTOCOL_V1.md.

gab_firth_evaluate <- function(beta, design, y) {
  eta <- drop(design %*% beta)
  if (any(!is.finite(eta))) return(list(value = -Inf, score = NULL))
  probability <- stats::plogis(eta)
  # Stable logistic variance; no flooring, clipping, ridge, or pseudo-data.
  weight <- stats::plogis(eta) * stats::plogis(-eta)
  information <- crossprod(design, design * weight)
  factor <- tryCatch(chol(information), error = function(error) NULL)
  if (is.null(factor)) return(list(value = -Inf, score = NULL))
  value <- -sum(pmax(eta, 0) - y * eta + log1p(exp(-abs(eta)))) +
    sum(log(diag(factor)))
  inverse <- chol2inv(factor)
  leverage <- weight * rowSums((design %*% inverse) * design)
  score <- drop(crossprod(design, y - probability +
    leverage * (0.5 - probability)))
  list(value = value, score = score, information = information,
    probability = probability, weight = weight)
}

gab_firth_fit <- function(design, y, controls) {
  gab_assert(is.matrix(design) && nrow(design) == length(y) &&
    all(is.finite(design)) && all(y %in% c(0, 1)) &&
    length(unique(y)) == 2L, "Invalid Firth training design or response.")
  beta <- numeric(ncol(design))
  initial <- gab_firth_evaluate(beta, design, y)
  gab_assert(is.finite(initial$value), "Firth design is not full column rank.")
  cache_beta <- NULL; cache_value <- NULL
  evaluate <- function(b) {
    if (!identical(b, cache_beta)) {
      cache_value <<- gab_firth_evaluate(b, design, y)
      cache_beta <<- b
    }
    cache_value
  }
  objective <- function(b) -evaluate(b)$value
  gradient <- function(b) {
    z <- evaluate(b)
    gab_assert(!is.null(z$score) && all(is.finite(z$score)),
      "Firth score is not numerically representable.")
    -z$score
  }
  evaluations <- 0L; accepted <- FALSE; polish_steps <- 0L
  for (pass in seq_len(controls$target_bfgs_passes)) {
    optimum <- stats::optim(beta, objective, gradient, method = "BFGS",
      control = list(maxit = controls$target_max_iterations,
        reltol = controls$target_bfgs_reltol))
    beta <- optimum$par
    evaluations <- evaluations + unname(optimum$counts[["function"]])
    fit <- gab_firth_evaluate(beta, design, y)
    score_norm <- if (is.null(fit$score)) Inf else max(abs(fit$score))
    # A small likelihood change alone is not evidence of convergence.
    fisher_step <- if (is.finite(score_norm))
      max(abs(solve(fit$information, fit$score))) else Inf
    accepted <- optimum$convergence == 0L && all(is.finite(beta)) &&
      is.finite(fit$value) && fit$value >= initial$value - 1e-10 &&
      is.finite(score_norm) && score_norm <= controls$target_tolerance &&
      is.finite(fisher_step) && fisher_step <= controls$target_tolerance
    if (accepted) break
  }
  # Near machine precision, objective-based BFGS stopping may leave an accurate
  # score but a poorly determined coefficient direction. Refine the same score
  # equations with the numerical Jacobian of the analytic gradient.
  if (!accepted && optimum$convergence == 0L) {
    for (iteration in seq_len(controls$target_polish_steps)) {
      hessian <- stats::optimHess(beta, objective, gradient,
        control = list(ndeps = rep(controls$target_hessian_step, length(beta))))
      curvature <- tryCatch(chol(hessian), error = function(error) NULL)
      if (is.null(curvature)) break
      delta <- drop(chol2inv(curvature) %*% fit$score)
      if (any(!is.finite(delta))) break
      moved <- FALSE
      for (halving in 0:20) {
        proposed <- beta + delta * 2^(-halving)
        candidate <- gab_firth_evaluate(proposed, design, y)
        allowance <- 64 * .Machine$double.eps * max(1, abs(fit$value))
        if (is.finite(candidate$value) &&
            candidate$value >= fit$value - allowance &&
            max(abs(candidate$score)) < score_norm) {
          beta <- proposed; fit <- candidate; moved <- TRUE
          break
        }
      }
      if (!moved) break
      polish_steps <- iteration
      score_norm <- max(abs(fit$score))
      fisher_step <- max(abs(solve(fit$information, fit$score)))
      accepted <- score_norm <= controls$target_tolerance &&
        fisher_step <= controls$target_tolerance
      if (accepted) break
    }
  }
  list(coefficient = beta, success = accepted,
    log_likelihood = fit$value, initial_log_likelihood = initial$value,
    maximum_adjusted_score = score_norm, maximum_fisher_step = fisher_step,
    optimizer_code = optimum$convergence, passes = pass, polish_steps = polish_steps,
    function_evaluations = evaluations)
}

gab_estimate_training_firth_target <- function(e, X, y, group, configuration) {
  started <- proc.time()[["elapsed"]]
  preprocess <- e$prepare_lsg_design(X, group)
  target <- numeric(ncol(preprocess$X))
  diagnostics <- vector("list", length(preprocess$blocks))
  for (g in seq_along(preprocess$blocks)) {
    block <- preprocess$blocks[[g]]
    fit <- gab_firth_fit(cbind(1, preprocess$X[, block$solver_index,
      drop = FALSE]), y, configuration)
    target[block$solver_index] <- fit$coefficient[-1L]
    diagnostics[[g]] <- data.frame(
      group = block$label, original_size = block$original_size,
      solver_rank = block$rank, success = fit$success,
      maximum_adjusted_score = fit$maximum_adjusted_score,
      maximum_fisher_step = fit$maximum_fisher_step,
      optimizer_code = fit$optimizer_code, passes = fit$passes,
      polish_steps = fit$polish_steps,
      function_evaluations = fit$function_evaluations,
      log_likelihood = fit$log_likelihood,
      initial_log_likelihood = fit$initial_log_likelihood,
      stringsAsFactors = FALSE
    )
  }
  diagnostics <- do.call(rbind, diagnostics)
  original <- as.numeric(e$recover_lsg_coefficients(preprocess, target, 0)[-1L])
  replay <- e$project_lsg_original_target(preprocess, original)
  replay_error <- max(abs(replay - target))
  failed <- !diagnostics$success
  gab_assert(!any(failed), paste0(
    "Training-only Firth BFGS target failed in group(s) ",
    paste(diagnostics$group[failed], collapse = ","),
    "; maximum adjusted score = ",
    format(max(diagnostics$maximum_adjusted_score), digits = 8),
    ". No fallback is allowed."))
  gab_assert(all(is.finite(original)) &&
    replay_error <= configuration$target_replay_tolerance *
      max(1, max(abs(target))), "Firth target round-trip reconstruction failed.")
  list(target_original = original, target_solver = target, success = TRUE,
    failed_groups = 0L, source = "training_solver_coordinates_firth_bfgs_v3",
    diagnostics = diagnostics, replay_error = replay_error,
    elapsed_seconds = proc.time()[["elapsed"]] - started)
}
