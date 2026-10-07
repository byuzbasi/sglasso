# Binomial package adapter. Numerical helpers/core are validated snapshots;
# this file owns input validation, API objects, and fold-local orchestration.
.sg_assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

.sg_binary_inputs <- function(X, Y, group) {
  X <- as.matrix(X)
  .sg_assert(is.numeric(X) && length(dim(X)) == 2L && nrow(X) >= 3L &&
    ncol(X) > 0L && all(is.finite(X)), "X must be a finite numeric matrix.")
  .sg_assert(is.null(dim(Y)) && (is.numeric(Y) || is.logical(Y)) &&
    length(Y) == nrow(X) && !anyNA(Y) && all(Y %in% c(0, 1)) &&
    length(unique(Y)) == 2L,
    "Binomial Y must be a 0/1 vector containing both classes; encode labels explicitly.")
  .sg_assert(length(group) == ncol(X) && !anyNA(group),
    "group must have one nonmissing label per predictor.")
  list(X = X, Y = as.numeric(Y))
}

.sg_binomial_controls <- function(control, eps, max_iter) {
  z <- list(block_max_sweeps = 4000L, block_chunk_sweeps = 100L,
    block_stall_window = 200L, block_stall_relative_improvement = 0.01,
    apg_max_iterations = 10000L, apg_kkt_check_interval = 5L,
    kkt_tolerance = eps, update_tolerance = 1e-10,
    intercept_tolerance = 1e-12, max_intercept_iterations = 100L,
    irls_max_outer = 50L, irls_max_inner = 2000L)
  .sg_assert(length(max_iter) == 1L && is.finite(max_iter) &&
    max_iter >= 1 && max_iter == floor(max_iter), "max_iter must be a positive integer.")
  .sg_assert(is.null(control) || (is.list(control) &&
    !is.null(names(control)) && !anyDuplicated(names(control)) &&
    all(names(control) %in% names(z))), "Unknown or unnamed binomial.control entries.")
  z[names(control)] <- control
  .sg_assert(all(vapply(z, function(v) is.numeric(v) && length(v) == 1L &&
    is.finite(v) && v > 0, logical(1))), "Binomial controls must be finite positive scalars.")
  ints <- c("block_max_sweeps", "block_chunk_sweeps", "block_stall_window",
    "apg_max_iterations", "apg_kkt_check_interval", "max_intercept_iterations",
    "irls_max_outer", "irls_max_inner")
  .sg_assert(all(vapply(z[ints], function(v) v == floor(v) &&
    v <= .Machine$integer.max, logical(1))), "Iteration controls must be positive integers.")
  .sg_assert(z$block_stall_relative_improvement < 1, "Stall improvement must be below one.")
  for (key in c("block_max_sweeps", "apg_max_iterations", "irls_max_outer",
    "irls_max_inner")) z[[key]] <- min(z[[key]], max_iter)
  z
}

.sg_binomial_reference <- function(preprocess, y, alpha) {
  score <- drop(crossprod(preprocess$X, y - mean(y))) / length(y)
  q <- max(vapply(seq_along(preprocess$blocks), function(g) {
    sqrt(sum(score[preprocess$blocks[[g]]$solver_index]^2)) /
      preprocess$group_weight[g]
  }, numeric(1)))
  .sg_assert(is.finite(q) && q > 0,
    "Zero null-score scale: supply an explicit positive finite lambda grid.")
  if (alpha > 0) q / alpha else q
}

.sg_binomial_fit <- function(X, Y, group, lambda, nlambda, d, nd, alpha,
    max_iter, eps, control, relative = NULL, use_active_set = TRUE) {
  input <- .sg_binary_inputs(X, Y, group)
  X <- input$X; Y <- input$Y
  .sg_assert(is.numeric(alpha) && length(alpha) == 1L && is.finite(alpha) &&
    alpha >= 0 && alpha <= 1, "alpha must be a scalar in [0,1].")
  if (is.null(d)) {
    .sg_assert(length(nd) == 1L && is.finite(nd) && nd >= 1 && nd == floor(nd),
      "nd must be a positive integer.")
    d <- seq(0, 1, length.out = nd)
  }
  .sg_assert(is.numeric(d) && length(d) > 0L && all(is.finite(d)) &&
    all(d >= 0 & d <= 1) && !anyDuplicated(d), "d must contain unique values in [0,1].")
  d <- sort(d)
  controls <- .sg_binomial_controls(control, eps, max_iter)
  pp <- .sg_binomial_design(X, group)
  reference <- NULL
  if (is.null(lambda)) {
    reference <- .sg_binomial_reference(pp, Y, alpha)
    lambda <- reference * relative
  }
  .sg_assert(is.numeric(lambda) && length(lambda) > 0L && all(is.finite(lambda)) &&
    all(lambda > 0) && !anyDuplicated(lambda),
    "lambda must contain unique positive finite values; zero/infinite endpoints are not fitted.")
  lambda <- sort(lambda, decreasing = TRUE)
  target <- numeric(ncol(pp$X)); target_info <- list()
  if (alpha < 1 && any(d > 0)) {
    fc <- list(target_bfgs_passes = 4L, target_max_iterations = 5000L,
      target_bfgs_reltol = 1e-15, target_tolerance = 1e-6,
      target_polish_steps = 5L, target_hessian_step = 1e-4)
    for (g in seq_along(pp$blocks)) {
      idx <- pp$blocks[[g]]$solver_index
      ft <- .sg_firth_fit(cbind(1, pp$X[, idx, drop = FALSE]), Y, fc)
      .sg_assert(ft$success, paste("Firth target failed for group", pp$blocks[[g]]$label,
        "(no MLE/ridge fallback)."))
      target[idx] <- ft$coefficient[-1L]
      target_info[[pp$blocks[[g]]$label]] <- ft
    }
  }
  core <- do.call(lsg_path_hybrid_v11_cpp, c(list(X = pp$X, y = Y,
    group_start = pp$group_start, group_end = pp$group_end,
    group_weight = pp$group_weight, target = target,
    lambda = lambda, d = d, alpha = alpha), controls,
    list(reuse_apg_offsets = TRUE, use_irls = TRUE, use_active_set = use_active_set)))
  .sg_assert(all(is.finite(core$beta)) && all(is.finite(core$intercept)) &&
    all(is.finite(core$objective)) && all(is.finite(core$kkt)) &&
    all(core$converged == 1) && all(core$kkt <= controls$kkt_tolerance),
    "Binomial path did not pass convergence/KKT checks; no partial fit returned. Increase iteration budgets or inspect the grid.")
  betas <- array(0, c(ncol(X) + 1L, length(lambda), length(d)),
    dimnames = list(c("(Intercept)", pp$names), as.character(lambda), as.character(d)))
  deviance <- matrix(0, length(lambda), length(d))
  for (j in seq_along(d)) for (i in seq_along(lambda)) {
    betas[, i, j] <- .sg_binomial_recover(pp, core$beta[, i, j], core$intercept[i, j])
    eta <- drop(cbind(1, X) %*% betas[, i, j])
    deviance[i, j] <- 2 * sum(pmax(eta, 0) - Y * eta + log1p(exp(-abs(eta))))
  }
  structure(list(family = "binomial", lambda = lambda, d = d, alpha = alpha,
    betas = betas, n = nrow(X), group = group, deviance = deviance,
    df = matrix(NA_real_, length(lambda), length(d)),
    preprocess = pp, target = target, target_diagnostics = target_info,
    target_original = .sg_binomial_recover(pp, target, 0)[-1L],
    controls = controls, converged = core$converged == 1,
    kkt = core$kkt, objective = core$objective, solver = core$solver,
    solver_diagnostics = core, lambda_reference = reference,
    lambda_relative = relative, standardize = TRUE, transform = "eager",
    screen = if (use_active_set) "hybrid_active_set" else "none", selected_groups = apply(core$beta, c(2, 3),
      function(b) sum(vapply(pp$blocks, function(g)
        sqrt(sum(b[g$solver_index]^2)) > 1e-8, logical(1)))),
    dropped_constant_columns = setdiff(seq_len(ncol(X)), pp$nonconstant)),
    class = "sglasso")
}

.sg_binomial_grid <- function(nlambda, ratio) {
  .sg_assert(length(nlambda) == 1L && is.finite(nlambda) && nlambda >= 1 &&
    nlambda == floor(nlambda), "nlambda must be a positive integer.")
  .sg_assert(length(ratio) == 1L && is.finite(ratio) && ratio > 0 && ratio < 1,
    "lambda.min.ratio must be in (0,1).")
  exp(seq(0, log(ratio), length.out = nlambda))
}

.sg_binomial_options <- function(standardize, bilevel, beta_start, screen, transform) {
  .sg_assert(isTRUE(standardize) && identical(bilevel, FALSE) && is.null(beta_start),
    "Binomial fitting requires standardize=TRUE, bilevel=FALSE and beta_start=NULL.")
  .sg_assert(identical(transform, "eager"), "Binomial fits currently require transform='eager'.")
  .sg_assert(screen %in% c("SSR", "none"),
    "SSR_fast is Gaussian-only; binomial fitting uses the validated hybrid active set.")
}

.sg_cv_binomial <- function(X, Y, group, lambda, nlambda, d, nd, alpha,
    fold, nfolds, max_iter, eps, dots, screen) {
  .sg_assert(all(names(dots) %in% c("lambda.min.ratio", "binomial.control")) &&
    !anyDuplicated(names(dots)) && (length(dots) == 0L || !is.null(names(dots))),
    "Unsupported binomial CV arguments in '...'.")
  input <- .sg_binary_inputs(X, Y, group); X <- input$X; Y <- input$Y
  n <- length(Y)
  if (is.null(fold)) {
    .sg_assert(length(nfolds) == 1L && is.finite(nfolds) && nfolds >= 2 &&
      nfolds == floor(nfolds) && nfolds <= min(table(Y)),
      "nfolds must be at least two and no larger than the smaller class count.")
    fold <- integer(n)
    for (value in c(0, 1)) {
      idx <- which(Y == value)
      fold[idx] <- sample(rep(seq_len(nfolds), length.out = length(idx)))
    }
  }
  .sg_assert(length(fold) == n && !anyNA(fold), "Invalid fold assignment.")
  fold <- match(fold, unique(fold)); nfolds <- max(fold)
  .sg_assert(nfolds >= 2L && all(vapply(seq_len(nfolds), function(k)
    length(unique(Y[fold != k])) == 2L, logical(1))),
    "Each fold training sample must contain both classes.")
  ratio <- if (is.null(dots$lambda.min.ratio)) 0.005 else dots$lambda.min.ratio
  relative <- if (is.null(lambda)) .sg_binomial_grid(nlambda, ratio) else NULL
  fit_one <- function(idx) .sg_binomial_fit(X[idx, , drop = FALSE], Y[idx], group,
    lambda, nlambda, d, nd, alpha, max_iter, eps, dots$binomial.control, relative,
    use_active_set = screen != "none")
  losses <- NULL; fold_fits <- vector("list", nfolds)
  for (k in seq_len(nfolds)) {
    f <- fit_one(which(fold != k))
    if (is.null(losses)) losses <- array(NA_real_, c(n, length(f$lambda), length(f$d)))
    eta <- predict(f, X[fold == k, , drop = FALSE], type = "link", drop = FALSE)
    losses[fold == k, , ] <- pmax(eta, 0) - Y[fold == k] * eta + log1p(exp(-abs(eta)))
    fold_fits[[k]] <- list(lambda = f$lambda, center = f$preprocess$center,
      scale = f$preprocess$scale, target_original = f$target_original,
      kkt = f$kkt, training_rows = which(fold != k))
  }
  .sg_assert(all(is.finite(losses)), "Nonfinite out-of-fold log-loss; CV not selected.")
  cve <- apply(losses, c(2, 3), mean)
  # Descriptive across-fold variability; not an independent-test confidence interval.
  fold_errors <- array(0, c(dim(losses)[2:3], nfolds))
  for (k in seq_len(nfolds)) fold_errors[, , k] <-
    apply(losses[fold == k, , , drop = FALSE], c(2, 3), mean)
  cvse <- apply(fold_errors, c(1, 2), stats::sd) / sqrt(nfolds)
  best <- arrayInd(which.min(cve), dim(cve))[1, ]
  full <- fit_one(seq_len(n))
  dimnames(cve) <- dimnames(cvse) <- list(as.character(full$lambda), as.character(full$d))
  boundary <- best[1] %in% c(1L, length(full$lambda))
  if (boundary) warning("CV selected a finite lambda-grid endpoint; performance beyond the grid is unassessed.", call. = FALSE)
  structure(list(family = "binomial", lambda = full$lambda, d = full$d, alpha = alpha,
    betas = full$betas, beta_opt = full$betas[, best[1], best[2]],
    fold = fold, lambda.min = full$lambda[best[1]], d.min = full$d[best[2]],
    cve = cve, cvse = cvse, fit = full, min_ind = best,
    measure = "log-loss", fold_fits = fold_fits, fold_errors = fold_errors,
    lambda_alignment = if (is.null(lambda)) "training_fold_relative_scale" else "explicit_absolute_lambda",
    lambda_relative = relative, selected_lambda_boundary = boundary), class = "cv.sglasso")
}
