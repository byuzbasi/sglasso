# Groupwise logistic targets for the isolated Logistic SGLASSO prework.

validate_groupwise_logistic_target_inputs <- function(X, y, group) {
  X <- as.matrix(X)
  y <- normalize_binary_response(y)$y
  if (nrow(X) != length(y)) {
    stop("X and y have incompatible dimensions.", call. = FALSE)
  }
  if (length(group) != ncol(X) || anyNA(group)) {
    stop("group must contain one non-missing label per predictor.", call. = FALSE)
  }
  if (length(unique(y)) != 2L) {
    stop("Both response classes are required.", call. = FALSE)
  }
  group_label <- as.character(group)
  list(
    X = X,
    y = y,
    group = group_label,
    group_levels = unique(group_label)
  )
}


capture_groupwise_fit <- function(expression) {
  warnings <- character(0)
  value <- tryCatch(
    withCallingHandlers(
      expression,
      warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) condition
  )
  list(value = value, warnings = unique(warnings))
}


estimate_groupwise_logistic_target <- function(
    X,
    y,
    group,
    method = c("mle", "firth"),
    max_iterations = 50L,
    tolerance = 1e-8
) {
  method <- match.arg(method)
  input <- validate_groupwise_logistic_target_inputs(X, y, group)
  X <- input$X
  y <- input$y
  group <- input$group
  levels <- input$group_levels
  max_iterations <- as.integer(max_iterations)
  if (max_iterations < 1L || !is.finite(tolerance) || tolerance <= 0) {
    stop("Invalid target-fitting controls.", call. = FALSE)
  }
  if (method == "firth" && !requireNamespace("logistf", quietly = TRUE)) {
    stop("The logistf package is required for the Firth target.", call. = FALSE)
  }

  target <- rep(NA_real_, ncol(X))
  diagnostics <- vector("list", length(levels))
  started <- proc.time()[["elapsed"]]

  for (position in seq_along(levels)) {
    level <- levels[position]
    index <- which(group == level)
    X_group <- X[, index, drop = FALSE]

    if (method == "mle") {
      design <- cbind("(Intercept)" = 1, X_group)
      attempt <- capture_groupwise_fit(stats::glm.fit(
        x = design,
        y = y,
        family = stats::binomial(),
        intercept = FALSE,
        control = stats::glm.control(
          epsilon = tolerance,
          maxit = max_iterations,
          trace = FALSE
        )
      ))
      fit <- attempt$value
      failed <- inherits(fit, "condition")
      coefficient <- if (failed) {
        rep(NA_real_, length(index) + 1L)
      } else {
        as.numeric(fit$coefficients)
      }
      finite <- length(coefficient) == length(index) + 1L &&
        all(is.finite(coefficient))
      converged <- !failed && isTRUE(fit$converged) && finite
      boundary <- !failed && isTRUE(fit$boundary)
      iterations <- if (failed) NA_integer_ else as.integer(fit$iter)
      message <- if (failed) conditionMessage(fit) else {
        paste(attempt$warnings, collapse = " | ")
      }
    } else {
      frame <- data.frame(response = y, X_group, check.names = FALSE)
      names(frame) <- c("response", paste0("x", seq_along(index)))
      control <- logistf::logistf.control(
        maxit = max_iterations,
        lconv = tolerance,
        gconv = tolerance,
        xconv = tolerance
      )
      attempt <- capture_groupwise_fit(logistf::logistf(
        response ~ .,
        data = frame,
        pl = FALSE,
        firth = TRUE,
        control = control,
        model = FALSE
      ))
      fit <- attempt$value
      failed <- inherits(fit, "condition")
      coefficient <- if (failed) {
        rep(NA_real_, length(index) + 1L)
      } else {
        as.numeric(fit$coefficients)
      }
      finite <- length(coefficient) == length(index) + 1L &&
        all(is.finite(coefficient))
      convergence_measure <- if (failed || is.null(fit$conv)) {
        rep(Inf, 3L)
      } else {
        abs(as.numeric(fit$conv))
      }
      converged <- !failed && finite &&
        all(is.finite(convergence_measure)) &&
        all(convergence_measure <= tolerance)
      boundary <- FALSE
      iterations <- if (failed || is.null(fit$iter)) {
        NA_integer_
      } else {
        as.integer(fit$iter[["full"]])
      }
      message <- if (failed) conditionMessage(fit) else {
        paste(attempt$warnings, collapse = " | ")
      }
    }

    slopes <- coefficient[-1L]
    if (converged) target[index] <- slopes
    maximum_absolute_slope <- if (any(is.finite(slopes))) {
      max(abs(slopes[is.finite(slopes)]))
    } else {
      NA_real_
    }
    warning_text <- paste(attempt$warnings, collapse = " | ")
    separation_suspected <- boundary ||
      grepl(
        "0 or 1 occurred|did not converge|separation",
        warning_text,
        ignore.case = TRUE
      ) ||
      (is.finite(maximum_absolute_slope) && maximum_absolute_slope > 20)

    diagnostics[[position]] <- data.frame(
      method = method,
      group = level,
      group_size = length(index),
      success = converged,
      converged = converged,
      finite_coefficients = finite,
      boundary = boundary,
      separation_suspected = separation_suspected,
      iterations = iterations,
      maximum_absolute_slope = maximum_absolute_slope,
      message = message,
      stringsAsFactors = FALSE
    )
  }

  diagnostics <- do.call(rbind, diagnostics)
  rownames(diagnostics) <- NULL
  success <- all(diagnostics$success) && all(is.finite(target))
  list(
    target_original = target,
    method = method,
    source = paste0("groupwise_logistic_", method),
    success = success,
    failed_groups = sum(!diagnostics$success),
    separation_groups = sum(diagnostics$separation_suspected),
    diagnostics = diagnostics,
    elapsed_seconds = proc.time()[["elapsed"]] - started
  )
}


direct_logistic_target_quality <- function(
    target,
    target_mode,
    true_beta,
    true_intercept,
    X
) {
  target <- as.numeric(target)
  true_beta <- as.numeric(true_beta)
  X <- as.matrix(X)
  valid <- length(target) == length(true_beta) && all(is.finite(target))
  if (!valid) {
    return(data.frame(
      target_mode = target_mode,
      target_available = FALSE,
      target_l2_norm = NA_real_,
      coefficient_l2_error = NA_real_,
      cosine_alignment = NA_real_,
      linear_predictor_rmse = NA_real_,
      relative_linear_predictor_rmse = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  difference <- target - true_beta
  target_link_error <- drop(X %*% difference)
  zero_link_error <- drop(X %*% true_beta)
  denominator <- sqrt(sum(target^2) * sum(true_beta^2))
  data.frame(
    target_mode = target_mode,
    target_available = TRUE,
    target_l2_norm = sqrt(sum(target^2)),
    coefficient_l2_error = sqrt(sum(difference^2)),
    cosine_alignment = if (denominator > 0) {
      sum(target * true_beta) / denominator
    } else {
      NA_real_
    },
    linear_predictor_rmse = sqrt(mean(target_link_error^2)),
    relative_linear_predictor_rmse =
      sqrt(mean(target_link_error^2)) / sqrt(mean(zero_link_error^2)),
    stringsAsFactors = FALSE
  )
}


logistic_target_on_original_scale <- function(fit) {
  if (!inherits(fit, "logistic_sglasso_prework")) {
    stop("fit must be a logistic_sglasso_prework object.", call. = FALSE)
  }
  recovered <- recover_lsg_coefficients(
    fit$preprocess,
    fit$target,
    intercept = 0
  )
  as.numeric(recovered[-1L])
}
