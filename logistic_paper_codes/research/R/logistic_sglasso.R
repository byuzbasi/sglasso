# Standalone R interface for the logistic SGLASSO proof of concept.
# This file intentionally does not modify or depend on the installed sglasso
# package internals.

lsg_prework_root <- function() {
  configured <- getOption("sglasso.logistic_prework.root", NULL)
  candidates <- unique(c(
    configured,
    file.path(getwd(), "logistic_prework"),
    getwd()
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  valid <- vapply(
    candidates,
    function(path) {
      file.exists(file.path(path, "src", "logistic_sglasso_core.cpp"))
    },
    logical(1)
  )
  if (!any(valid)) {
    stop(
      "Cannot locate logistic_prework. Set ",
      "options(sglasso.logistic_prework.root = '<path>').",
      call. = FALSE
    )
  }
  normalizePath(candidates[which(valid)[1]], mustWork = TRUE)
}


compile_lsg_core <- function(rebuild = FALSE, quiet = TRUE) {
  required <- c(
    "lsg_fisher_target_cpp",
    "lsg_lambda_start_cpp",
    "lsg_kkt_cpp",
    "lsg_objective_cpp",
    "lsg_path_cpp",
    "lsg_path_abgd_cpp"
  )
  if (!rebuild && all(vapply(required, exists, logical(1), mode = "function"))) {
    return(invisible(TRUE))
  }
  if (!requireNamespace("Rcpp", quietly = TRUE) ||
      !requireNamespace("RcppArmadillo", quietly = TRUE)) {
    stop("Rcpp and RcppArmadillo are required.", call. = FALSE)
  }

  root <- lsg_prework_root()
  local_makevars <- file.path(root, "config", "Makevars.local")
  local_gfortran_runtime <- paste0(
    "/usr/local/gfortran/lib/gcc/",
    "aarch64-apple-darwin23/14.1.0/libemutls_w.a"
  )
  previous_makevars <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  if (file.exists(local_makevars) && file.exists(local_gfortran_runtime)) {
    Sys.setenv(R_MAKEVARS_USER = local_makevars)
  }
  on.exit({
    if (is.na(previous_makevars)) {
      Sys.unsetenv("R_MAKEVARS_USER")
    } else {
      Sys.setenv(R_MAKEVARS_USER = previous_makevars)
    }
  }, add = TRUE)

  Rcpp::sourceCpp(
    file.path(root, "src", "logistic_sglasso_core.cpp"),
    rebuild = rebuild,
    showOutput = !quiet,
    verbose = FALSE,
    env = globalenv()
  )
  invisible(TRUE)
}


normalize_binary_response <- function(y) {
  if (is.data.frame(y) || is.matrix(y)) y <- drop(as.matrix(y))
  if (length(y) == 0L || anyNA(y)) {
    stop("y must be a non-missing binary vector.", call. = FALSE)
  }

  if (is.logical(y)) {
    labels <- c("FALSE", "TRUE")
    out <- as.numeric(y)
  } else if (is.factor(y) || is.character(y)) {
    yf <- factor(y)
    if (nlevels(yf) != 2L) {
      stop("y must contain exactly two classes.", call. = FALSE)
    }
    labels <- levels(yf)
    out <- as.numeric(yf == labels[2L])
  } else {
    values <- sort(unique(as.numeric(y)))
    if (length(values) != 2L || any(!is.finite(values))) {
      stop("y must contain exactly two finite classes.", call. = FALSE)
    }
    labels <- as.character(values)
    out <- as.numeric(as.numeric(y) == values[2L])
  }

  if (!all(out %in% c(0, 1)) || length(unique(out)) != 2L) {
    stop("Both outcome classes are required.", call. = FALSE)
  }
  list(y = out, labels = labels)
}


prepare_lsg_design <- function(X, group, rank_tolerance = 1e-10) {
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
  group_start <- as.integer(c(0L, head(cumulative, -1L)))
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


transform_lsg_newx <- function(preprocess, newx) {
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


recover_lsg_coefficients <- function(preprocess, beta, intercept) {
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


project_lsg_original_target <- function(preprocess, target_original) {
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


lsg_null_score_lambda_scale <- function(preprocess, y) {
  if (!inherits(preprocess, "lsg_preprocess")) {
    stop("preprocess must be an lsg_preprocess object.", call. = FALSE)
  }
  response <- normalize_binary_response(y)$y
  if (length(response) != preprocess$n) {
    stop("preprocess and y have incompatible dimensions.", call. = FALSE)
  }
  score <- drop(crossprod(
    preprocess$X,
    response - mean(response)
  )) / preprocess$n
  group_score <- vapply(
    seq_along(preprocess$group_start),
    function(index) {
      first <- preprocess$group_start[index] + 1L
      last <- preprocess$group_end[index] + 1L
      sqrt(sum(score[first:last]^2)) / preprocess$group_weight[index]
    },
    numeric(1)
  )
  reference <- max(group_score)
  if (!is.finite(reference) || reference <= 0) {
    stop("Unable to construct a positive null-score lambda scale.",
         call. = FALSE)
  }
  reference
}


fit_logistic_sglasso <- function(
    X,
    y,
    group,
    lambda,
    nlambda = 35L,
    lambda_min_ratio = 0.02,
    d = c(0, 0.5, 1),
    alpha = 0.8,
    rank_tolerance = 1e-10,
    max_outer = 250L,
    max_inner = 2000L,
    tolerance = 1e-7,
    inner_tolerance = 1e-9,
    keep_traces = FALSE,
    target_original = NULL,
    compile = TRUE,
    preprocess = NULL,
    use_active_set = TRUE,
    warm_start_d = TRUE,
    solver = c("abgd", "mm_reference")
) {
  if (compile) compile_lsg_core()
  solver <- match.arg(solver)
  response <- normalize_binary_response(y)
  y01 <- response$y
  if (nrow(as.matrix(X)) != length(y01)) {
    stop("X and y have incompatible dimensions.", call. = FALSE)
  }
  if (length(alpha) != 1L || !is.finite(alpha) ||
      alpha < 0 || alpha > 1) {
    stop("alpha must lie in [0, 1].", call. = FALSE)
  }
  if (length(use_active_set) != 1L || is.na(use_active_set) ||
      length(warm_start_d) != 1L || is.na(warm_start_d)) {
    stop("use_active_set and warm_start_d must be single logical values.",
         call. = FALSE)
  }
  use_active_set <- isTRUE(use_active_set)
  warm_start_d <- isTRUE(warm_start_d)
  d <- sort(unique(as.numeric(d)))
  if (length(d) == 0L || any(!is.finite(d)) || any(d < 0 | d > 1)) {
    stop("d must contain finite values in [0, 1].", call. = FALSE)
  }

  if (is.null(preprocess)) {
    preprocess <- prepare_lsg_design(
      X,
      group,
      rank_tolerance = rank_tolerance
    )
  } else if (!inherits(preprocess, "lsg_preprocess") ||
             preprocess$n != nrow(as.matrix(X)) ||
             preprocess$p_original != ncol(as.matrix(X))) {
    stop("preprocess is incompatible with X.", call. = FALSE)
  }
  if (is.null(target_original)) {
    target_information <- lsg_fisher_target_cpp(
      preprocess$X,
      y01,
      preprocess$group_start,
      preprocess$group_end
    )
    target_information$source <- "fisher_null_score"
    target <- drop(target_information$target)
  } else {
    target_original <- as.numeric(target_original)
    target <- project_lsg_original_target(preprocess, target_original)
    target_information <- list(
      target = target,
      source = "supplied_original_coefficients",
      original_target = target_original
    )
  }

  lambda_diagnostics <- lapply(
    d,
    function(d_value) {
      lsg_lambda_start_cpp(
        preprocess$X,
        y01,
        preprocess$group_start,
        preprocess$group_end,
        preprocess$group_weight,
        target,
        alpha,
        d_value
      )
    }
  )
  names(lambda_diagnostics) <- format(d, trim = TRUE)

  if (missing(lambda)) {
    if (abs(alpha) <= 1e-12) {
      stop(
        paste(
          "At alpha = 0 the null model has no finite KKT lambda",
          "boundary; supply an explicit positive lambda path."
        ),
        call. = FALSE
      )
    }
    nlambda <- as.integer(nlambda)
    if (nlambda < 2L) stop("nlambda must be at least 2.", call. = FALSE)
    if (!is.finite(lambda_min_ratio) ||
        lambda_min_ratio <= 0 || lambda_min_ratio >= 1) {
      stop("lambda_min_ratio must lie in (0, 1).", call. = FALSE)
    }
    # Use the d=0 null-model KKT boundary as a common, d-independent
    # reference scale. For d>0 the first fitted model is allowed to be
    # nonzero: a shifted quadratic target need not admit a null solution at
    # any lambda. Exact null feasibility remains available as a diagnostic.
    reference_diagnostic <- lsg_lambda_start_cpp(
      preprocess$X,
      y01,
      preprocess$group_start,
      preprocess$group_end,
      preprocess$group_weight,
      target,
      alpha,
      0
    )
    lambda_start <- reference_diagnostic$lambda_start
    lambda_upper <- reference_diagnostic$common_upper
    if (!isTRUE(reference_diagnostic$zero_model_feasible) ||
        !is.finite(lambda_start) || lambda_start <= 0) {
      stop(
        "Unable to construct the d=0 KKT reference lambda.",
        call. = FALSE
      )
    }
    lambda_start <- lambda_start * (1 + 1e-8)
    lambda <- exp(
      seq(
        log(lambda_start),
        log(lambda_start * lambda_min_ratio),
        length.out = nlambda
      )
    )
    user_lambda <- FALSE
    lambda_reference <- "d0_null_kkt"
  } else {
    lambda <- as.numeric(lambda)
    if (length(lambda) == 0L || any(!is.finite(lambda)) ||
        any(lambda <= 0)) {
      stop("lambda must contain finite positive values.", call. = FALSE)
    }
    if (is.unsorted(-lambda, strictly = FALSE)) {
      stop("lambda must be supplied in non-increasing order.", call. = FALSE)
    }
    nlambda <- length(lambda)
    user_lambda <- TRUE
    lambda_start <- lambda[1L]
    lambda_upper <- NA_real_
    lambda_reference <- "user_supplied"
  }

  if (identical(solver, "abgd")) {
    core <- lsg_path_abgd_cpp(
      preprocess$X,
      y01,
      preprocess$group_start,
      preprocess$group_end,
      preprocess$group_weight,
      target,
      lambda,
      d,
      alpha,
      as.integer(max_outer),
      tolerance,
      inner_tolerance,
      use_active_set,
      warm_start_d,
      keep_traces
    )
  } else {
    core <- lsg_path_cpp(
      preprocess$X,
      y01,
      preprocess$group_start,
      preprocess$group_end,
      preprocess$group_weight,
      target,
      lambda,
      d,
      alpha,
      as.integer(max_outer),
      as.integer(max_inner),
      tolerance,
      inner_tolerance,
      use_active_set,
      warm_start_d,
      keep_traces
    )
  }

  coefficient_array <- array(
    0,
    dim = c(preprocess$p_original + 1L, length(lambda), length(d)),
    dimnames = list(
      c("(Intercept)", preprocess$names),
      signif(lambda, 5),
      d
    )
  )
  for (di in seq_along(d)) {
    for (li in seq_along(lambda)) {
      coefficient_array[, li, di] <- recover_lsg_coefficients(
        preprocess,
        core$beta[, li, di],
        core$intercept[li, di]
      )
    }
  }

  out <- list(
    coefficients = coefficient_array,
    beta_solver = core$beta,
    intercept_solver = core$intercept,
    lambda = lambda,
    d = d,
    alpha = alpha,
    solver = solver,
    target = target,
    target_source = target_information$source,
    target_information = target_information,
    lambda_diagnostics = lambda_diagnostics,
    lambda_start = lambda_start,
    lambda_upper = lambda_upper,
    lambda_reference = lambda_reference,
    zero_model_feasible = vapply(
      lambda_diagnostics,
      function(item) item$zero_model_feasible,
      logical(1)
    ),
    user_lambda = user_lambda,
    preprocess = preprocess,
    y = y01,
    outcome_labels = response$labels,
    objective = core$objective,
    kkt = core$kkt,
    intercept_kkt = core$intercept_kkt,
    converged = core$converged != 0,
    outer_iterations = core$outer_iterations,
    inner_iterations = core$inner_iterations,
    selected_groups = core$selected_groups,
    active_set_scans = core$active_set_scans,
    passes = if (!is.null(core$passes)) core$passes else core$outer_iterations,
    kkt_scans = if (!is.null(core$kkt_scans)) core$kkt_scans else NULL,
    group_curvature = if (!is.null(core$group_curvature)) {
      core$group_curvature
    } else {
      NULL
    },
    group_updates = core$group_updates,
    use_active_set = use_active_set,
    warm_start_d = warm_start_d,
    maximum_raw_objective_increase = core$max_raw_objective_increase,
    objective_traces = core$objective_traces,
    majorization_curvature = core$majorization_curvature,
    call = match.call()
  )
  class(out) <- "logistic_sglasso_prework"
  out
}


predict_logistic_sglasso <- function(
    object,
    newx,
    type = c("response", "link"),
    lambda_index = seq_along(object$lambda),
    d_index = seq_along(object$d)
) {
  if (!inherits(object, "logistic_sglasso_prework")) {
    stop("object must be a logistic_sglasso_prework fit.", call. = FALSE)
  }
  type <- match.arg(type)
  lambda_index <- as.integer(lambda_index)
  d_index <- as.integer(d_index)
  if (any(!lambda_index %in% seq_along(object$lambda)) ||
      any(!d_index %in% seq_along(object$d))) {
    stop("Invalid lambda_index or d_index.", call. = FALSE)
  }
  transformed <- transform_lsg_newx(object$preprocess, newx)
  eta <- array(
    NA_real_,
    dim = c(nrow(transformed), length(lambda_index), length(d_index)),
    dimnames = list(
      NULL,
      signif(object$lambda[lambda_index], 5),
      object$d[d_index]
    )
  )
  for (di_position in seq_along(d_index)) {
    di <- d_index[di_position]
    eta[, , di_position] <-
      transformed %*%
        object$beta_solver[, lambda_index, di, drop = FALSE][, , 1L] +
      matrix(
        object$intercept_solver[lambda_index, di],
        nrow = nrow(transformed),
        ncol = length(lambda_index),
        byrow = TRUE
      )
  }
  if (type == "link") return(eta)
  stats::plogis(eta)
}


binary_log_loss <- function(y, probability, epsilon = 1e-12) {
  y <- as.numeric(y)
  probability <- pmin(pmax(as.numeric(probability), epsilon), 1 - epsilon)
  -mean(y * log(probability) + (1 - y) * log1p(-probability))
}


binary_auc <- function(y, probability) {
  y <- as.numeric(y)
  probability <- as.numeric(probability)
  positives <- sum(y == 1)
  negatives <- sum(y == 0)
  if (positives == 0L || negatives == 0L) return(NA_real_)
  ranks <- rank(probability, ties.method = "average")
  (sum(ranks[y == 1]) - positives * (positives + 1) / 2) /
    (positives * negatives)
}


stratified_folds <- function(y, nfolds = 5L, seed = 20260819L) {
  y <- normalize_binary_response(y)$y
  nfolds <- as.integer(nfolds)
  if (nfolds < 2L) stop("nfolds must be at least 2.", call. = FALSE)
  class_counts <- table(y)
  if (any(class_counts < nfolds)) {
    stop("Each outcome class must have at least nfolds observations.", call. = FALSE)
  }
  set.seed(seed)
  fold <- integer(length(y))
  for (class_value in c(0, 1)) {
    index <- sample(which(y == class_value))
    fold[index] <- rep(seq_len(nfolds), length.out = length(index))
  }
  fold
}


cv_logistic_sglasso <- function(
    X,
    y,
    group,
    nfolds = 5L,
    fold,
    seed = 20260819L,
    nlambda = 30L,
    lambda_min_ratio = 0.03,
    d = c(0, 0.5, 1),
    alpha = 0.8,
    rank_tolerance = 1e-10,
    max_outer = 250L,
    max_inner = 2000L,
    tolerance = 1e-7,
    inner_tolerance = 1e-9,
    target_original = NULL,
    compile = TRUE,
    solver = c("abgd", "mm_reference")
) {
  if (compile) compile_lsg_core()
  solver <- match.arg(solver)
  response <- normalize_binary_response(y)
  y01 <- response$y
  X <- as.matrix(X)
  if (missing(fold)) {
    fold <- stratified_folds(y01, nfolds = nfolds, seed = seed)
  } else {
    fold <- as.integer(fold)
    if (length(fold) != length(y01) || anyNA(fold)) {
      stop("fold must have one non-missing entry per observation.", call. = FALSE)
    }
    nfolds <- length(unique(fold))
  }
  fold_labels <- sort(unique(fold))
  d <- sort(unique(as.numeric(d)))
  fold_loss <- array(
    NA_real_,
    dim = c(nlambda, length(d), length(fold_labels)),
    dimnames = list(
      lambda_fraction = seq_len(nlambda),
      d = d,
      fold = fold_labels
    )
  )
  fold_kkt <- fold_loss
  fold_converged <- array(
    FALSE,
    dim = dim(fold_loss),
    dimnames = dimnames(fold_loss)
  )
  fold_maximum_raw_objective_increase <- numeric(length(fold_labels))
  fold_prevalence <- numeric(length(fold_labels))

  for (fold_position in seq_along(fold_labels)) {
    held_out <- fold == fold_labels[fold_position]
    fit <- fit_logistic_sglasso(
      X = X[!held_out, , drop = FALSE],
      y = y01[!held_out],
      group = group,
      nlambda = nlambda,
      lambda_min_ratio = lambda_min_ratio,
      d = d,
      alpha = alpha,
      rank_tolerance = rank_tolerance,
      max_outer = max_outer,
      max_inner = max_inner,
      tolerance = tolerance,
      inner_tolerance = inner_tolerance,
      keep_traces = FALSE,
      target_original = target_original,
      compile = FALSE,
      solver = solver
    )
    probability <- predict_logistic_sglasso(
      fit,
      X[held_out, , drop = FALSE],
      type = "response"
    )
    for (di in seq_along(d)) {
      for (li in seq_len(nlambda)) {
        fold_loss[li, di, fold_position] <- binary_log_loss(
          y01[held_out],
          probability[, li, di]
        )
      }
    }
    fold_kkt[, , fold_position] <- fit$kkt
    fold_converged[, , fold_position] <- fit$converged
    fold_maximum_raw_objective_increase[fold_position] <-
      max(fit$maximum_raw_objective_increase)
    fold_prevalence[fold_position] <- mean(y01[held_out])
  }

  mean_loss <- apply(fold_loss, c(1, 2), mean)
  se_loss <- apply(fold_loss, c(1, 2), stats::sd) /
    sqrt(length(fold_labels))
  minimum <- which(mean_loss == min(mean_loss), arr.ind = TRUE)[1L, ]

  full_fit <- fit_logistic_sglasso(
    X = X,
    y = y01,
    group = group,
    nlambda = nlambda,
    lambda_min_ratio = lambda_min_ratio,
    d = d,
    alpha = alpha,
    rank_tolerance = rank_tolerance,
    max_outer = max_outer,
    max_inner = max_inner,
    tolerance = tolerance,
    inner_tolerance = inner_tolerance,
    keep_traces = FALSE,
    target_original = target_original,
    compile = FALSE,
    solver = solver
  )

  out <- list(
    fit = full_fit,
    fold = fold,
    fold_loss = fold_loss,
    fold_kkt = fold_kkt,
    fold_converged = fold_converged,
    fold_maximum_raw_objective_increase =
      fold_maximum_raw_objective_increase,
    mean_loss = mean_loss,
    se_loss = se_loss,
    lambda_index_min = minimum[1L],
    d_index_min = minimum[2L],
    lambda_min = full_fit$lambda[minimum[1L]],
    d_min = d[minimum[2L]],
    alpha = alpha,
    solver = solver,
    lambda_fraction = full_fit$lambda / full_fit$lambda[1L],
    fold_prevalence = fold_prevalence,
    outcome_labels = response$labels,
    call = match.call()
  )
  class(out) <- "cv_logistic_sglasso_prework"
  out
}


predict_cv_logistic_sglasso <- function(object, newx, type = c("response", "link")) {
  if (!inherits(object, "cv_logistic_sglasso_prework")) {
    stop("object must be a cv_logistic_sglasso_prework object.", call. = FALSE)
  }
  type <- match.arg(type)
  drop(
    predict_logistic_sglasso(
      object$fit,
      newx,
      type = type,
      lambda_index = object$lambda_index_min,
      d_index = object$d_index_min
    )
  )
}
