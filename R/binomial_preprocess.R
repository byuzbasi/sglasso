# Internal numerical helpers ported from the validated research implementation.
.sg_binomial_design <- function(X, group, rank_tolerance = 1e-10) {
  if (!is.matrix(X)) X <- as.matrix(X)
  storage.mode(X) <- "double"
  if (length(group) != ncol(X)) {
    stop("group must have one entry per column of X.", call. = FALSE)
  }
  if (nrow(X) < 2L || ncol(X) < 1L || anyNA(X) || any(!is.finite(X))) {
    stop("X must be a finite, non-missing numeric matrix.", call. = FALSE)
  }
  if (!is.numeric(rank_tolerance) || length(rank_tolerance) != 1L ||
      rank_tolerance <= 0) {
    stop("rank_tolerance must be positive.", call. = FALSE)
  }

  n <- nrow(X)
  p_original <- ncol(X)
  x_names <- colnames(X)
  if (is.null(x_names)) x_names <- paste0("X", seq_len(p_original))
  group_character <- as.character(group)

  center <- colMeans(X)
  centered <- sweep(X, 2L, center, "-")
  scale <- sqrt(colMeans(centered^2))
  nonconstant <- which(is.finite(scale) & scale > 1e-8)
  if (length(nonconstant) == 0L) {
    stop("All columns of X are constant.", call. = FALSE)
  }
  standardized <- sweep(
    centered[, nonconstant, drop = FALSE],
    2L,
    scale[nonconstant],
    "/"
  )

  group_nonconstant <- group_character[nonconstant]
  group_levels <- unique(group_nonconstant)
  group_index <- match(group_nonconstant, group_levels)
  order_index <- order(group_index, seq_along(group_index))
  standardized <- standardized[, order_index, drop = FALSE]
  ordered_original_index <- nonconstant[order_index]
  ordered_group <- group_index[order_index]

  transformed_blocks <- vector("list", length(group_levels))
  block_metadata <- vector("list", length(group_levels))
  retained_labels <- character(0)
  retained_original_sizes <- integer(0)
  retained_ranks <- integer(0)
  maximum_orthonormality_error <- 0

  for (g in seq_along(group_levels)) {
    columns <- which(ordered_group == g)
    if (length(columns) == 0L) next
    Xg <- standardized[, columns, drop = FALSE]
    decomposition <- svd(Xg, nu = min(dim(Xg)), nv = min(dim(Xg)))
    if (length(decomposition$d) == 0L ||
        max(decomposition$d) <= .Machine$double.eps) {
      next
    }
    keep <- which(
      decomposition$d >
        rank_tolerance * max(decomposition$d) * max(dim(Xg))
    )
    if (length(keep) == 0L) next

    transform <- sweep(
      decomposition$v[, keep, drop = FALSE],
      2L,
      sqrt(n) / decomposition$d[keep],
      "*"
    )
    transformed <- Xg %*% transform
    orthonormality_error <- max(
      abs(crossprod(transformed) / n - diag(length(keep)))
    )
    maximum_orthonormality_error <- max(
      maximum_orthonormality_error,
      orthonormality_error
    )

    transformed_blocks[[g]] <- transformed
    block_metadata[[g]] <- list(
      label = group_levels[g],
      original_index = ordered_original_index[columns],
      transform = transform,
      rank = length(keep),
      original_size = length(columns)
    )
    retained_labels <- c(retained_labels, group_levels[g])
    retained_original_sizes <- c(retained_original_sizes, length(columns))
    retained_ranks <- c(retained_ranks, length(keep))
  }

  keep_blocks <- !vapply(transformed_blocks, is.null, logical(1))
  transformed_blocks <- transformed_blocks[keep_blocks]
  block_metadata <- block_metadata[keep_blocks]
  if (length(transformed_blocks) == 0L) {
    stop("No full-rank predictor group remains after preprocessing.", call. = FALSE)
  }

  X_solver <- do.call(cbind, transformed_blocks)
  group_size <- vapply(
    block_metadata,
    function(block) block$rank,
    integer(1)
  )
  cumulative <- cumsum(group_size)
  group_start <- as.integer(c(0L, utils::head(cumulative, -1L)))
  group_end <- as.integer(cumulative - 1L)
  solver_index <- Map(
    function(first, last) seq.int(first + 1L, last + 1L),
    group_start,
    group_end
  )
  for (g in seq_along(block_metadata)) {
    block_metadata[[g]]$solver_index <- solver_index[[g]]
  }

  structure(
    list(
      X = X_solver,
      n = n,
      p_original = p_original,
      names = x_names,
      center = center,
      scale = scale,
      nonconstant = nonconstant,
      blocks = block_metadata,
      group_labels = retained_labels,
      group_original_size = retained_original_sizes,
      group_size = group_size,
      group_start = group_start,
      group_end = group_end,
      group_weight = sqrt(group_size),
      maximum_orthonormality_error = maximum_orthonormality_error,
      rank_tolerance = rank_tolerance
    ),
    class = "lsg_preprocess"
  )
}


.sg_binomial_transform <- function(preprocess, newx) {
  if (!inherits(preprocess, "lsg_preprocess")) {
    stop("preprocess must be an lsg_preprocess object.", call. = FALSE)
  }
  if (!is.matrix(newx)) newx <- as.matrix(newx)
  storage.mode(newx) <- "double"
  if (ncol(newx) != preprocess$p_original) {
    stop("newx has a different number of columns than the training X.", call. = FALSE)
  }
  if (anyNA(newx) || any(!is.finite(newx))) {
    stop("newx must be finite and non-missing.", call. = FALSE)
  }

  out <- matrix(
    0,
    nrow = nrow(newx),
    ncol = ncol(preprocess$X)
  )
  for (block in preprocess$blocks) {
    original <- block$original_index
    standardized <- sweep(
      sweep(
        newx[, original, drop = FALSE],
        2L,
        preprocess$center[original],
        "-"
      ),
      2L,
      preprocess$scale[original],
      "/"
    )
    out[, block$solver_index] <- standardized %*% block$transform
  }
  out
}


.sg_binomial_recover <- function(preprocess, beta, intercept) {
  beta <- as.numeric(beta)
  if (length(beta) != ncol(preprocess$X)) {
    stop("beta has the wrong transformed dimension.", call. = FALSE)
  }
  original_beta <- numeric(preprocess$p_original)
  for (block in preprocess$blocks) {
    standardized_beta <-
      block$transform %*% beta[block$solver_index]
    original_beta[block$original_index] <-
      drop(standardized_beta) / preprocess$scale[block$original_index]
  }
  original_intercept <- intercept -
    sum(preprocess$center * original_beta)
  stats::setNames(
    c(original_intercept, original_beta),
    c("(Intercept)", preprocess$names)
  )
}


.sg_binomial_project <- function(preprocess, target_original) {
  if (!inherits(preprocess, "lsg_preprocess")) {
    stop("preprocess must be an lsg_preprocess object.", call. = FALSE)
  }
  target_original <- as.numeric(target_original)
  if (length(target_original) != preprocess$p_original ||
      any(!is.finite(target_original))) {
    stop(
      "target_original must contain one finite slope per original predictor.",
      call. = FALSE
    )
  }

  target_solver <- numeric(ncol(preprocess$X))
  for (block in preprocess$blocks) {
    original <- block$original_index
    standardized_target <-
      target_original[original] * preprocess$scale[original]
    target_solver[block$solver_index] <- drop(qr.solve(
      block$transform,
      standardized_target,
      tol = preprocess$rank_tolerance
    ))
  }
  target_solver
}
.sg_firth_evaluate <- function(beta, design, y) {
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

.sg_firth_fit <- function(design, y, controls) {
  .sg_assert(is.matrix(design) && nrow(design) == length(y) &&
    all(is.finite(design)) && all(y %in% c(0, 1)) &&
    length(unique(y)) == 2L, "Invalid Firth training design or response.")
  beta <- numeric(ncol(design))
  initial <- .sg_firth_evaluate(beta, design, y)
  .sg_assert(is.finite(initial$value), "Firth design is not full column rank.")
  cache_beta <- NULL; cache_value <- NULL
  evaluate <- function(b) {
    if (!identical(b, cache_beta)) {
      cache_value <<- .sg_firth_evaluate(b, design, y)
      cache_beta <<- b
    }
    cache_value
  }
  objective <- function(b) -evaluate(b)$value
  gradient <- function(b) {
    z <- evaluate(b)
    .sg_assert(!is.null(z$score) && all(is.finite(z$score)),
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
    fit <- .sg_firth_evaluate(beta, design, y)
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
        candidate <- .sg_firth_evaluate(proposed, design, y)
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
