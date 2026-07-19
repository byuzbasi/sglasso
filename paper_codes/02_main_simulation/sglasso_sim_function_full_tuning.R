############################################################
# SGLASSO reviewer-revision simulation functions
#
# Purpose
# -------
# This file provides a clean simulation framework for the revised
# SGLASSO manuscript.  The code is written to directly address the
# reviewers' requests:
#
#   1. Remove the old gglasso baseline.
#   2. Add ADELIE as a recent competitor.
#   3. Remove old gliunet/gllreg naming and use SGLASSO explicitly.
#   4. Add heterogeneous coefficient settings:
#        - homogeneous
#        - mixed_signs
#        - weak_strong_mixed
#   5. Use a clear train-validation-test workflow:
#        - fit on training data,
#        - select tuning parameters using a chosen criterion,
#        - evaluate final performance on the independent test data.
#   6. Add additional evaluation metrics requested by the reviewer:
#        - test prediction MSE,
#        - MSE_y,
#        - MSE_beta,
#        - PVE,
#        - FDR,
#        - selected groups,
#        - selected features,
#        - runtime,
#        - selected alpha/lambda/d.
#   7. Add tuning criteria for sensitivity analysis:
#        - Min_val: validation MSE
#        - Risk: oracle risk, used only to check consistency with the
#          previous manuscript simulations
#        - AIC
#        - BIC
#        - EBIC
#        - GCV
#
# Important design note
# ---------------------
# The main method-comparison tables should use the SAME tuning criterion
# for all methods, typically Min_val.  The Risk criterion uses the true
# beta and should only be used for diagnostic comparison with the older
# simulation code.  AIC/BIC/EBIC/GCV are intended for the reviewer's
# requested tuning-criterion validation analysis.
############################################################


############################################################
# 0. Required packages
############################################################

required_packages <- c(
  "simstudy", "data.table", "grpnet", "sglasso", "grpreg", "adelie",
  "caret", "dplyr", "MASS", "Metrics", "matrixcalc", "readr",
  "foreach", "doParallel", "doRNG", "ggplot2", "tibble"
)

github_packages <- list(
  sglasso = "byuzbasi/sglasso"
  #grpnet  = "YaohuiZeng/grpnet"
)

install_and_load_packages <- function(pkgs = required_packages,
                                      install_missing = TRUE,
                                      lib = Sys.getenv("R_LIBS_USER")) {

  if (nzchar(lib)) {
    dir.create(lib, recursive = TRUE, showWarnings = FALSE)
    .libPaths(c(lib, .libPaths()))
  } else {
    lib <- .libPaths()[1]
  }

  options(repos = c(CRAN = "https://cloud.r-project.org"))

  for (pkg in pkgs) {

    if (!requireNamespace(pkg, quietly = TRUE)) {

      if (!install_missing) {
        stop("Package not installed: ", pkg)
      }

      message("Installing missing package: ", pkg)

      if (pkg %in% names(github_packages)) {

        if (!requireNamespace("remotes", quietly = TRUE)) {
          install.packages("remotes", lib = lib)
        }

        remotes::install_github(
          github_packages[[pkg]],
          lib = lib,
          upgrade = "never",
          force = FALSE
        )

      } else {

        install.packages(pkg, lib = lib)
      }
    }

    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }

  invisible(TRUE)
}


############################################################
# 1. ADELIE helper functions
############################################################
# ADELIE expects group information as the starting column index of each
# group, not as a repeated group label vector.  For example, if
# group = c(1,1,1,2,2,2), the starts are c(1,4).

make_adelie_group_starts <- function(group) {
  as.integer(tapply(seq_along(group), group, min))
}

# In some ADELIE versions, coef(fit) may fail because of a family-object
# issue.  The following functions extract coefficient paths using the
# internal predict.grpnet method.  The returned matrix has dimension
# (p + 1) x nlambda, where the first row is the intercept.

adelie_coef_object <- function(fit) {
  adelie:::predict.grpnet(fit, type = "coef")
}

adelie_coef_matrix <- function(fit) {
  pc <- adelie_coef_object(fit)
  B <- as.matrix(pc$betas)
  intercepts <- as.numeric(pc$intercepts)
  rbind(intercepts, t(B))
}

adelie_lambda <- function(fit) {
  pc <- adelie_coef_object(fit)
  as.numeric(pc$lambda)
}


############################################################
# 2. Heterogeneous coefficient generation
############################################################
# This function generates the group-level signal strengths for the active
# groups.  Each active group receives one scalar signal value, and all
# features inside that group receive this same value.
#
# homogeneous:
#   beta values are all equal to eff_nonzero.
#
# mixed_signs:
#   beta values have equal magnitude but alternating signs.
#
# weak_strong_mixed:
#   active groups contain both strong and weak signals, with alternating
#   signs.  This setting addresses the reviewer request for a more
#   heterogeneous and realistic coefficient structure.

make_beta_group_values <- function(
    strong_J,
    eff_nonzero = 1,
    signal_pattern = c("homogeneous", "mixed_signs", "weak_strong_mixed")
) {

  signal_pattern <- match.arg(signal_pattern)

  if (signal_pattern == "homogeneous") {
    return(rep(eff_nonzero, strong_J))
  }

  if (signal_pattern == "mixed_signs") {
    vals <- rep(eff_nonzero, strong_J)
    signs <- rep(c(1, -1), length.out = strong_J)
    return(vals * signs)
  }

  if (signal_pattern == "weak_strong_mixed") {

    n_strong <- ceiling(strong_J / 2)
    n_weak   <- floor(strong_J / 2)

    strong_vals <- runif(
      n_strong,
      min = 0.8 * eff_nonzero,
      max = 1.2 * eff_nonzero
    )

    weak_vals <- runif(
      n_weak,
      min = 0.15 * eff_nonzero,
      max = 0.35 * eff_nonzero
    )

    vals <- c(rbind(
      strong_vals[seq_len(min(length(strong_vals), length(weak_vals)))],
      weak_vals
    ))

    if (length(vals) < strong_J) {
      vals <- c(vals, strong_vals[length(strong_vals)])
    }

    vals <- vals[seq_len(strong_J)]
    signs <- rep(c(1, -1), length.out = strong_J)

    return(vals * signs)
  }
}


############################################################
# 3. Data-generation functions
############################################################

simsham_corrmat_toeplitz <- function(p, rho) {
  if (rho == 0) return(diag(p))
  rho^abs(outer(seq_len(p), seq_len(p), "-"))
}

block_sim_data <- function(n_train = 100,
                           n_val = 100,
                           n_test = 400,
                           pj = 5,
                           J = 20,
                           strong_J = 4,
                           rho_w = 0.9,
                           rho_b = 0.5,
                           eff_nonzero = 1,
                           corrmat_type = c("Exchangeable", "Decay", "Toeplitz"),
                           snr = 1.528,
                           nonzero_id = NULL,
                           signal_pattern = c("homogeneous",
                                              "mixed_signs",
                                              "weak_strong_mixed")) {

  corrmat_type <- match.arg(corrmat_type)
  signal_pattern <- match.arg(signal_pattern)

  p <- pj * J
  if (is.null(nonzero_id)) nonzero_id <- sample(J, strong_J)

  X_train <- matrix(rnorm(n_train * p), n_train, p)
  X_val   <- matrix(rnorm(n_val * p), n_val, p)
  X_test  <- matrix(rnorm(n_test * p), n_test, p)

  if (corrmat_type == "Decay") {
    corrmat <- simstudy::blockDecayMat(
      ninds = pj, nperiods = J,
      rho_w = rho_w, r = rho_b,
      pattern = "xsection"
    )
  } else if (corrmat_type == "Exchangeable") {
    corrmat <- simstudy::blockExchangeMat(
      ninds = pj, nperiods = J,
      rho_w = rho_w, rho_b = rho_b,
      pattern = "xsection"
    )
  } else {
    corrmat <- simsham_corrmat_toeplitz(p, rho_b)
  }

  sv <- svd(corrmat)
  corrmat_half <- sv$u %*% sqrt(diag(pmax(sv$d, 0))) %*% t(sv$v)

  X_train <- X_train %*% corrmat_half
  X_val   <- X_val %*% corrmat_half
  X_test  <- X_test %*% corrmat_half

  group <- rep(seq_len(J), each = pj)

  beta_values <- make_beta_group_values(
    strong_J = strong_J,
    eff_nonzero = eff_nonzero,
    signal_pattern = signal_pattern
  )

  beta <- rep(0, p)
  for (k in seq_len(strong_J)) {
    idx <- ((nonzero_id[k] - 1) * pj + 1):(nonzero_id[k] * pj)
    beta[idx] <- beta_values[k]
  }

  risk_null <- as.numeric(t(beta) %*% corrmat %*% beta)
  sigma <- sqrt(risk_null / snr)
  error_null <- risk_null + sigma^2

  y_train <- as.numeric(X_train %*% beta + rnorm(n_train) * sigma)
  y_val   <- as.numeric(X_val   %*% beta + rnorm(n_val)   * sigma)
  y_test  <- as.numeric(X_test  %*% beta + rnorm(n_test)  * sigma)

  true_groups <- sort(nonzero_id)

  list(
    X_train = X_train,
    X_val = X_val,
    X_test = X_test,
    y_train = y_train,
    y_val = y_val,
    y_test = y_test,
    beta = beta,
    group = group,
    true_groups = true_groups,
    corrmat = corrmat,
    sigma = sigma,
    risk_null = risk_null,
    error_null = error_null,
    n_train = n_train,
    n_val = n_val,
    n_test = n_test,
    p = p,
    pj = pj,
    J = J,
    strong_J = strong_J,
    rho_w = rho_w,
    rho_b = rho_b,
    snr = snr,
    signal_pattern = signal_pattern
  )
}

block_sim_data_exchangeable_fast <- function(n_train = 100,
                                             n_val = 100,
                                             n_test = 400,
                                             pj = 5,
                                             J = 20,
                                             strong_J = 4,
                                             rho_w = 0.9,
                                             rho_b = 0.5,
                                             eff_nonzero = 1,
                                             snr = 1.528,
                                             nonzero_id = NULL,
                                             signal_pattern = c("homogeneous",
                                                                "mixed_signs",
                                                                "weak_strong_mixed")) {
  signal_pattern <- match.arg(signal_pattern)

  if (rho_b < 0 || rho_w < rho_b || rho_w > 1) {
    stop("Fast exchangeable generator requires 0 <= rho_b <= rho_w <= 1.")
  }

  p <- pj * J
  if (is.null(nonzero_id)) nonzero_id <- sample(J, strong_J)

  make_X <- function(n) {
    if (n == 0) return(matrix(numeric(0), nrow = 0, ncol = p))
    common_factor <- rnorm(n)
    group_factor <- matrix(rnorm(n * J), nrow = n, ncol = J)
    noise <- matrix(rnorm(n * p), nrow = n, ncol = p)
    sqrt(rho_b) * common_factor +
      sqrt(rho_w - rho_b) * group_factor[, rep(seq_len(J), each = pj), drop = FALSE] +
      sqrt(1 - rho_w) * noise
  }

  X_train <- make_X(n_train)
  X_val <- make_X(n_val)
  X_test <- make_X(n_test)

  group <- rep(seq_len(J), each = pj)

  beta_values <- make_beta_group_values(
    strong_J = strong_J,
    eff_nonzero = eff_nonzero,
    signal_pattern = signal_pattern
  )

  beta <- rep(0, p)
  for (k in seq_len(strong_J)) {
    idx <- ((nonzero_id[k] - 1) * pj + 1):(nonzero_id[k] * pj)
    beta[idx] <- beta_values[k]
  }

  group_sum <- as.numeric(rowsum(beta, group, reorder = FALSE))
  risk_null <- (1 - rho_w) * sum(beta^2) +
    (rho_w - rho_b) * sum(group_sum^2) +
    rho_b * sum(beta)^2
  sigma <- sqrt(risk_null / snr)
  error_null <- risk_null + sigma^2

  y_train <- as.numeric(X_train %*% beta + rnorm(n_train) * sigma)
  y_val <- as.numeric(X_val %*% beta + rnorm(n_val) * sigma)
  y_test <- as.numeric(X_test %*% beta + rnorm(n_test) * sigma)

  true_groups <- sort(nonzero_id)

  list(
    X_train = X_train,
    X_val = X_val,
    X_test = X_test,
    y_train = y_train,
    y_val = y_val,
    y_test = y_test,
    beta = beta,
    group = group,
    true_groups = true_groups,
    corrmat = NULL,
    sigma = sigma,
    risk_null = risk_null,
    error_null = error_null,
    n_train = n_train,
    n_val = n_val,
    n_test = n_test,
    p = p,
    pj = pj,
    J = J,
    strong_J = strong_J,
    rho_w = rho_w,
    rho_b = rho_b,
    snr = snr,
    signal_pattern = signal_pattern
  )
}


############################################################
# 4. Train-only standardization
############################################################
# If standardize = TRUE, only training-data means and standard deviations
# are used.  The same transformation is then applied to validation and
# test data.  This avoids information leakage.

standardize_by_train <- function(X_train, X_val, X_test, y_train,
                                 standardize = TRUE) {

  if (!standardize) {
    return(list(
      X_train = X_train,
      X_val = X_val,
      X_test = X_test,
      y_train = y_train,
      x_center = rep(0, ncol(X_train)),
      x_scale = rep(1, ncol(X_train)),
      y_center = 0,
      standardized = FALSE
    ))
  }

  x_center <- colMeans(X_train)
  x_scale <- apply(X_train, 2, sd)
  x_scale[is.na(x_scale) | x_scale == 0] <- 1
  y_center <- mean(y_train)

  list(
    X_train = scale(X_train, center = x_center, scale = x_scale),
    X_val   = scale(X_val,   center = x_center, scale = x_scale),
    X_test  = scale(X_test,  center = x_center, scale = x_scale),
    y_train = y_train - y_center,
    x_center = x_center,
    x_scale = x_scale,
    y_center = y_center,
    standardized = TRUE
  )
}

make_intercept_design <- function(X) cbind(1, X)

standardized_coef_to_original <- function(B, x_center, x_scale, y_center) {
  if (is.null(dim(B))) B <- matrix(B, ncol = 1)
  B <- as.matrix(B)

  beta_std <- B[-1, , drop = FALSE]
  beta_orig <- sweep(beta_std, 1, x_scale, "/")
  intercept_orig <- B[1, ] + y_center -
    as.numeric(crossprod(x_center, beta_orig))

  rbind(intercept_orig, beta_orig)
}

true_coef_to_standardized <- function(beta_true, x_center, x_scale, y_center) {
  beta_std <- beta_true * x_scale
  intercept_std <- sum(x_center * beta_true) - y_center
  c(intercept_std, beta_std)
}

exchangeable_block_quadratic <- function(delta, group, rho_w, rho_b) {
  group_sum <- as.numeric(rowsum(delta, group, reorder = FALSE))
  (1 - rho_w) * sum(delta^2) +
    (rho_w - rho_b) * sum(group_sum^2) +
    rho_b * sum(delta)^2
}

call_sglasso <- function(...,
                         standardize = TRUE,
                         screen = NULL,
                         transform = NULL) {
  args <- list(...)
  sglasso_args <- names(formals(sglasso::sglasso))
  if ("standardize" %in% sglasso_args) {
    args$standardize <- standardize
  }
  if (!is.null(screen) && "screen" %in% sglasso_args) {
    args$screen <- screen
  }
  if (!is.null(transform) && "transform" %in% sglasso_args) {
    args$transform <- transform
  }
  do.call(sglasso::sglasso, args)
}

parallel_lapply <- function(X, FUN, cores = 1, ...) {
  if (cores <= 1 || length(X) <= 1) {
    return(lapply(X, FUN, ...))
  }

  if (.Platform$OS.type == "windows") {
    warning("parallel_lapply falls back to serial execution on Windows.",
            call. = FALSE)
    return(lapply(X, FUN, ...))
  }

  parallel::mclapply(
    X,
    FUN,
    ...,
    mc.cores = min(as.integer(cores), length(X)),
    mc.preschedule = FALSE
  )
}


############################################################
# 5. Generic helper functions
############################################################

coef_to_matrix <- function(obj) {
  cf <- coef(obj)
  if (is.matrix(cf)) return(cf)
  if (inherits(cf, "dgCMatrix")) return(as.matrix(cf))
  if (is.list(cf)) {
    if (!is.null(cf$beta)) return(as.matrix(cf$beta))
    if (!is.null(cf$coef)) return(as.matrix(cf$coef))
  }
  as.matrix(cf)
}

safe_lambda <- function(lambda_vec, lambda_id) {
  if (is.null(lambda_vec)) return(NA_real_)
  if (length(lambda_vec) < lambda_id) return(NA_real_)
  as.numeric(lambda_vec[lambda_id])
}

selected_groups_from_beta <- function(beta, group, tol = 1e-8) {
  sort(unique(group[abs(beta) > tol]))
}

selection_metrics <- function(selected_groups, true_groups, J) {
  selected_groups <- sort(unique(selected_groups))
  true_groups <- sort(unique(true_groups))

  TP <- length(intersect(selected_groups, true_groups))
  FP <- length(setdiff(selected_groups, true_groups))
  FN <- length(setdiff(true_groups, selected_groups))
  TN <- J - TP - FP - FN

  FDR <- ifelse(TP + FP == 0, 0, FP / (TP + FP))
  TPR <- ifelse(TP + FN == 0, NA_real_, TP / (TP + FN))
  TNR <- ifelse(TN + FP == 0, NA_real_, TN / (TN + FP))

  c(TP = TP, FP = FP, TN = TN, FN = FN,
    FDR = FDR, TPR = TPR, TNR = TNR)
}


############################################################
# 6. Tuning criteria
############################################################
# B is a coefficient-path matrix of dimension (p + 1) x nlambda.
# The first row is the intercept.
#
# Min_val:
#   validation MSE. This is the main practical tuning rule.
#
# Risk:
#   oracle risk using true beta and Sigma. This should only be used for
#   checking consistency with the previous simulations because it uses
#   information unavailable in real data.
#
# AIC/BIC/EBIC/GCV:
#   information-criterion and generalized cross-validation alternatives.
#   These are used for the reviewer-requested tuning sensitivity analysis.

validation_mse_path <- function(B, X_val_std, y_val_centered) {
  pred <- make_intercept_design(X_val_std) %*% B
  colMeans((y_val_centered - pred)^2)
}

compute_tuning_score <- function(
    B,
    dat,
    std_dat,
    criterion = c("Min_val", "Val_Risk", "Risk", "AIC", "BIC", "EBIC", "GCV"),
    ebic_gamma = 0.5
) {

  criterion <- match.arg(criterion)

  beta_path <- B[-1, , drop = FALSE]

  pred_val <- make_intercept_design(std_dat$X_val) %*% B
  y_val_centered <- dat$y_val - std_dat$y_center

  val_mse <- colMeans((y_val_centered - pred_val)^2)

  if (criterion == "Min_val") {
    return(val_mse)
  }

  if (criterion == "Val_Risk") {
    true_val <- as.numeric(dat$X_val %*% dat$beta - std_dat$y_center)

    val_risk <- colMeans(
      sweep(pred_val, 1, true_val, "-")^2
    )

    return(val_risk)
  }

  if (criterion == "Risk") {
    B_orig <- standardized_coef_to_original(
      B,
      std_dat$x_center,
      std_dat$x_scale,
      std_dat$y_center
    )
    beta_path_orig <- B_orig[-1, , drop = FALSE]
    delta <- sweep(beta_path_orig, 1, dat$beta, "-")
    risk <- if (!is.null(dat$corrmat)) {
      diag(t(delta) %*% dat$corrmat %*% delta)
    } else {
      apply(delta, 2, exchangeable_block_quadratic,
            group = dat$group, rho_w = dat$rho_w, rho_b = dat$rho_b)
    }
    return(risk)
  }

  pred_train <- make_intercept_design(std_dat$X_train) %*% B
  rss <- colSums((std_dat$y_train - pred_train)^2)

  n <- nrow(std_dat$X_train)
  J <- dat$J

  df_features <- colSums(abs(beta_path) > 1e-8) + 1

  df_groups <- apply(beta_path, 2, function(b) {
    length(unique(dat$group[abs(b) > 1e-8]))
  })

  if (criterion == "AIC") {
    return(n * log(pmax(rss / n, .Machine$double.eps)) + 2 * df_features)
  }

  if (criterion == "BIC") {
    return(n * log(pmax(rss / n, .Machine$double.eps)) + log(n) * df_features)
  }

  if (criterion == "EBIC") {
    ebic_penalty <- 2 * ebic_gamma * lchoose(J, pmin(df_groups, J))

    return(
      n * log(pmax(rss / n, .Machine$double.eps)) +
        log(n) * df_features +
        ebic_penalty
    )
  }

  if (criterion == "GCV") {
    return((rss / n) / pmax((1 - df_features / n)^2, .Machine$double.eps))
  }
}


############################################################
# 7. Evaluation
############################################################

empty_method_result <- function(method,
                                message = NA_character_,
                                tuning_criterion = NA_character_) {
  data.frame(
    rep = NA_integer_,
    method = method,
    time = NA_real_,
    iter = NA_real_,
    risk = NA_real_,
    error_test = NA_real_,
    pve = NA_real_,
    rte_bayes = NA_real_,
    rte_null = NA_real_,
    pe_mse = NA_real_,
    pe_mae = NA_real_,
    pe_rmse = NA_real_,
    MSE_y = NA_real_,
    MSE_beta = NA_real_,
    MSE_beta_cov = NA_real_,
    selected_groups = NA_real_,
    selected_features = NA_real_,
    sparsity = NA_real_,
    TP = NA_real_,
    FP = NA_real_,
    TN = NA_real_,
    FN = NA_real_,
    FDR = NA_real_,
    TPR = NA_real_,
    TNR = NA_real_,
    alpha = NA_real_,
    lambda = NA_real_,
    d = NA_real_,
    tuning_criterion = tuning_criterion,
    converged = FALSE,
    message = message,
    stringsAsFactors = FALSE
  )
}

evaluate_beta <- function(method,
                          beta_hat,
                          beta_true,
                          X_test_orig,
                          X_test_std,
                          y_test,
                          y_center,
                          x_center,
                          x_scale,
                          group,
                          true_groups,
                          corrmat,
                          COV_X_test,
                          sigma2,
                          error_null,
                          risk_null,
                          alpha = NA_real_,
                          lambda = NA_real_,
                          d = NA_real_,
                          tuning_criterion = NA_character_,
                          time = NA_real_,
                          iter = NA_real_,
                          message = NA_character_) {

  if (is.null(dim(beta_hat))) beta_hat <- matrix(beta_hat, ncol = 1)
  beta_hat <- as.matrix(beta_hat)

  beta_orig <- standardized_coef_to_original(
    beta_hat,
    x_center,
    x_scale,
    y_center
  )[, 1]

  beta_no_int <- beta_orig[-1]
  delta <- beta_no_int - beta_true

  risk <- if (!is.null(corrmat)) {
    as.numeric(t(delta) %*% corrmat %*% delta)
  } else {
    exchangeable_block_quadratic(
      delta = delta,
      group = group,
      rho_w = attr(beta_true, "rho_w"),
      rho_b = attr(beta_true, "rho_b")
    )
  }
  error_test <- risk + sigma2
  pve <- 1 - error_test / error_null
  rte_bayes <- error_test / sigma2
  rte_null <- risk / risk_null

  pred_test <- as.numeric(make_intercept_design(X_test_orig) %*% beta_orig)
  pe_mse <- mean((y_test - pred_test)^2)
  pe_mae <- mean(abs(y_test - pred_test))
  pe_rmse <- sqrt(pe_mse)

  true_mean_test <- as.numeric(X_test_orig %*% beta_true)
  MSE_y <- mean((true_mean_test - pred_test)^2)
  MSE_beta <- sum(delta^2)
  MSE_beta_cov <- if (!is.null(COV_X_test)) {
    as.numeric(t(delta) %*% COV_X_test %*% delta)
  } else {
    risk
  }

  sel_groups <- selected_groups_from_beta(beta_no_int, group)
  sel_metrics <- selection_metrics(sel_groups, true_groups, length(unique(group)))

  data.frame(
    rep = NA_integer_,
    method = method,
    time = as.numeric(time),
    iter = as.numeric(iter),
    risk = risk,
    error_test = error_test,
    pve = pve,
    rte_bayes = rte_bayes,
    rte_null = rte_null,
    pe_mse = pe_mse,
    pe_mae = pe_mae,
    pe_rmse = pe_rmse,
    MSE_y = MSE_y,
    MSE_beta = MSE_beta,
    MSE_beta_cov = MSE_beta_cov,
    selected_groups = length(sel_groups),
    selected_features = sum(abs(beta_no_int) > 1e-8),
    sparsity = mean(abs(beta_no_int) <= 1e-8),
    TP = sel_metrics["TP"],
    FP = sel_metrics["FP"],
    TN = sel_metrics["TN"],
    FN = sel_metrics["FN"],
    FDR = sel_metrics["FDR"],
    TPR = sel_metrics["TPR"],
    TNR = sel_metrics["TNR"],
    alpha = alpha,
    lambda = lambda,
    d = d,
    tuning_criterion = tuning_criterion,
    converged = TRUE,
    message = message,
    stringsAsFactors = FALSE
  )
}


############################################################
# 8. Method-specific fitting functions
############################################################

fit_oracle <- function(dat, std_dat, COV_X_test,
                       tuning_criterion = "Oracle") {
  beta_hat <- true_coef_to_standardized(
    dat$beta,
    std_dat$x_center,
    std_dat$x_scale,
    std_dat$y_center
  )

  evaluate_beta(
    method = "ORACLE",
    beta_hat = beta_hat,
    beta_true = dat$beta,
    X_test_orig = dat$X_test,
    X_test_std = std_dat$X_test,
    y_test = dat$y_test,
    y_center = std_dat$y_center,
    x_center = std_dat$x_center,
    x_scale = std_dat$x_scale,
    group = dat$group,
    true_groups = dat$true_groups,
    corrmat = dat$corrmat,
    COV_X_test = COV_X_test,
    sigma2 = dat$sigma^2,
    error_null = dat$error_null,
    risk_null = dat$risk_null,
    tuning_criterion = tuning_criterion,
    time = 0,
    iter = 0,
    message = "True active groups known."
  )
}

fit_grpreg_method <- function(method_name,
                              penalty,
                              dat,
                              std_dat,
                              COV_X_test,
                              nlambda,
                              eps,
                              maxit,
                              tuning_criterion = "Min_val",
                              ebic_gamma = 0.5,
                              lambda_min_ratio = 0.01) {
  out <- tryCatch({

    t0 <- proc.time()

    fit <- grpreg::grpreg(
      std_dat$X_train,
      std_dat$y_train,
      group = dat$group,
      penalty = penalty,
      max.iter = maxit,
      eps = eps,
      lambda.min = lambda_min_ratio,
      nlambda = nlambda
    )

    elapsed <- (proc.time() - t0)[3]

    B <- as.matrix(coef(fit))

    score <- compute_tuning_score(
      B = B,
      dat = dat,
      std_dat = std_dat,
      criterion = tuning_criterion,
      ebic_gamma = ebic_gamma
    )

    id <- which.min(score)

    evaluate_beta(
      method = method_name,
      beta_hat = B[, id],
      beta_true = dat$beta,
      X_test_orig = dat$X_test,
      X_test_std = std_dat$X_test,
      y_test = dat$y_test,
      y_center = std_dat$y_center,
      x_center = std_dat$x_center,
      x_scale = std_dat$x_scale,
      group = dat$group,
      true_groups = dat$true_groups,
      corrmat = dat$corrmat,
      COV_X_test = COV_X_test,
      sigma2 = dat$sigma^2,
      error_null = dat$error_null,
      risk_null = dat$risk_null,
      alpha = 1,
      lambda = safe_lambda(fit$lambda, id),
      tuning_criterion = tuning_criterion,
      time = elapsed,
      iter = sum(fit$iter, na.rm = TRUE)
    )
  }, error = function(e) {
    empty_method_result(method_name, conditionMessage(e), tuning_criterion)
  })

  out
}

fit_adelie <- function(dat,
                       std_dat,
                       COV_X_test,
                       alpha_seq,
                       nlambda,
                       eps,
                       maxit,
                       n_threads = 1,
                       alpha_cores = 1,
                       tuning_criterion = "Min_val",
                       ebic_gamma = 0.5,
                       lambda_min_ratio = 0.01) {
  out <- tryCatch({

    group_starts <- make_adelie_group_starts(dat$group)

    t0 <- proc.time()

    fits <- parallel_lapply(alpha_seq, function(a) {
      adelie::grpnet(
        X = std_dat$X_train,
        glm = adelie::glm.gaussian(std_dat$y_train),
        groups = group_starts,
        alpha = a,
        standardize = FALSE,
        intercept = TRUE,
        lmda_path_size = nlambda,
        min_ratio = lambda_min_ratio,
        tol = eps,
        max_iters = as.integer(maxit),
        screen_rule = "strong",
        progress_bar = FALSE,
        n_threads = n_threads
      )
    }, cores = alpha_cores)

    elapsed <- (proc.time() - t0)[3]

    best <- list(
      score = Inf,
      alpha_id = NA_integer_,
      lambda_id = NA_integer_,
      B = NULL,
      lambda = NA_real_
    )

    for (ai in seq_along(fits)) {

      B <- adelie_coef_matrix(fits[[ai]])
      lambda_vec <- adelie_lambda(fits[[ai]])

      score <- compute_tuning_score(
        B = B,
        dat = dat,
        std_dat = std_dat,
        criterion = tuning_criterion,
        ebic_gamma = ebic_gamma
      )

      li <- which.min(score)

      if (score[li] < best$score) {
        best <- list(
          score = score[li],
          alpha_id = ai,
          lambda_id = li,
          B = B[, li],
          lambda = safe_lambda(lambda_vec, li)
        )
      }
    }

    evaluate_beta(
      method = "ADELIE",
      beta_hat = best$B,
      beta_true = dat$beta,
      X_test_orig = dat$X_test,
      X_test_std = std_dat$X_test,
      y_test = dat$y_test,
      y_center = std_dat$y_center,
      x_center = std_dat$x_center,
      x_scale = std_dat$x_scale,
      group = dat$group,
      true_groups = dat$true_groups,
      corrmat = dat$corrmat,
      COV_X_test = COV_X_test,
      sigma2 = dat$sigma^2,
      error_null = dat$error_null,
      risk_null = dat$risk_null,
      alpha = alpha_seq[best$alpha_id],
      lambda = best$lambda,
      tuning_criterion = tuning_criterion,
      time = elapsed / length(alpha_seq),
      iter = NA_real_
    )
  }, error = function(e) {
    empty_method_result("ADELIE", conditionMessage(e), tuning_criterion)
  })

  out
}

fit_genet <- function(dat,
                      std_dat,
                      COV_X_test,
                      alpha_seq,
                      nlambda,
                      eps,
                      maxit,
                      alpha_cores = 1,
                      tuning_criterion = "Min_val",
                      ebic_gamma = 0.5,
                      lambda_min_ratio = 0.01) {
  out <- tryCatch({

    t0 <- proc.time()

    fits <- parallel_lapply(alpha_seq, function(a) {
      grpnet::grpnet(
        std_dat$X_train,
        std_dat$y_train,
        dat$group,
        nlambda = nlambda,
        alpha = a,
        thresh = eps,
        maxit = maxit,
        lambda.min.ratio = lambda_min_ratio
      )
    }, cores = alpha_cores)

    elapsed <- (proc.time() - t0)[3]

    coef_grpnet <- getS3method(
      "coef",
      "grpnet",
      envir = asNamespace("grpnet")
    )

    best <- list(
      score = Inf,
      alpha_id = NA_integer_,
      lambda_id = NA_integer_,
      B = NULL,
      fit = NULL
    )

    for (ai in seq_along(fits)) {

      B <- as.matrix(coef_grpnet(fits[[ai]]))

      score <- compute_tuning_score(
        B = B,
        dat = dat,
        std_dat = std_dat,
        criterion = tuning_criterion,
        ebic_gamma = ebic_gamma
      )

      score[!is.finite(score)] <- Inf

      if (all(is.infinite(score))) {
        next
      }

      li <- which.min(score)

      if (score[li] < best$score) {
        best <- list(
          score = score[li],
          alpha_id = ai,
          lambda_id = li,
          B = B[, li],
          fit = fits[[ai]]
        )
      }
    }

    if (is.null(best$B)) {
      stop("GENET failed: all tuning scores are non-finite.")
    }

    evaluate_beta(
      method = "GENET",
      beta_hat = best$B,
      beta_true = dat$beta,
      X_test_orig = dat$X_test,
      X_test_std = std_dat$X_test,
      y_test = dat$y_test,
      y_center = std_dat$y_center,
      x_center = std_dat$x_center,
      x_scale = std_dat$x_scale,
      group = dat$group,
      true_groups = dat$true_groups,
      corrmat = dat$corrmat,
      COV_X_test = COV_X_test,
      sigma2 = dat$sigma^2,
      error_null = dat$error_null,
      risk_null = dat$risk_null,
      alpha = alpha_seq[best$alpha_id],
      lambda = safe_lambda(best$fit$lambda, best$lambda_id),
      tuning_criterion = tuning_criterion,
      time = elapsed / length(alpha_seq),
      iter = mean(sapply(fits, function(z) {
        sum(z$npasses, na.rm = TRUE)
      }), na.rm = TRUE)
    )
  }, error = function(e) {
    empty_method_result(
      "GENET",
      paste0("GENET failed: ", conditionMessage(e)),
      tuning_criterion
    )
  })

  out
}

fit_sglasso <- function(dat,
                        std_dat,
                        COV_X_test,
                        alpha_seq,
                        d_seq,
                        nlambda,
                        eps,
                        maxit,
                        alpha_cores = 1,
                        tuning_criterion = "Min_val",
                        ebic_gamma = 0.5,
                        lambda_min_ratio = 0.01) {
  out <- tryCatch({

    t0 <- proc.time()

    fits <- parallel_lapply(alpha_seq, function(a) {
      call_sglasso(
        X = std_dat$X_train,
        Y = std_dat$y_train,
        group = dat$group,
        nlambda = nlambda,
        alpha = a,
        d = d_seq,
        eps = eps,
        max_iter = maxit,
        standardize = !isTRUE(std_dat$standardized),
        lambda.min.ratio = lambda_min_ratio
      )
    }, cores = alpha_cores)

    elapsed <- (proc.time() - t0)[3]

    best <- list(
      score = Inf,
      alpha_id = NA_integer_,
      lambda_id = NA_integer_,
      d_id = NA_integer_,
      B = NULL,
      fit = NULL
    )

    for (ai in seq_along(fits)) {
      for (di in seq_along(d_seq)) {

        B <- fits[[ai]]$betas[, , di, drop = FALSE][, , 1]

        score <- compute_tuning_score(
          B = B,
          dat = dat,
          std_dat = std_dat,
          criterion = tuning_criterion,
          ebic_gamma = ebic_gamma
        )

        li <- which.min(score)

        if (score[li] < best$score) {
          best <- list(
            score = score[li],
            alpha_id = ai,
            lambda_id = li,
            d_id = di,
            B = B[, li],
            fit = fits[[ai]]
          )
        }
      }
    }

    iter_val <- NA_real_
    if (!is.null(fits[[1]]$iter)) {
      iter_val <- mean(sapply(fits, function(z) {
        sum(z$iter, na.rm = TRUE)
      }), na.rm = TRUE)
    }

    evaluate_beta(
      method = "SGLASSO",
      beta_hat = best$B,
      beta_true = dat$beta,
      X_test_orig = dat$X_test,
      X_test_std = std_dat$X_test,
      y_test = dat$y_test,
      y_center = std_dat$y_center,
      x_center = std_dat$x_center,
      x_scale = std_dat$x_scale,
      group = dat$group,
      true_groups = dat$true_groups,
      corrmat = dat$corrmat,
      COV_X_test = COV_X_test,
      sigma2 = dat$sigma^2,
      error_null = dat$error_null,
      risk_null = dat$risk_null,
      alpha = alpha_seq[best$alpha_id],
      lambda = safe_lambda(best$fit$lambda, best$lambda_id),
      d = d_seq[best$d_id],
      tuning_criterion = tuning_criterion,
      time = elapsed / (length(alpha_seq) * length(d_seq)),
      iter = iter_val
    )
  }, error = function(e) {
    empty_method_result("SGLASSO", conditionMessage(e), tuning_criterion)
  })

  out
}


############################################################
# 9. One replicate and simulation wrapper
############################################################

fit_one_sim_replicate <- function(rep_id = 1,
                                  n_train = 100,
                                  n_val = 100,
                                  n_test = 100,
                                  pj = 5,
                                  J = 20,
                                  strong_J = 4,
                                  rho_w = 0.9,
                                  rho_b = 0.5,
                                  eff_nonzero = 1,
                                  corrmat_type = "Exchangeable",
                                  snr = 1.528,
                                  signal_pattern = "homogeneous",
                                  alpha_seq = seq(0.1, 0.9, by = 0.1),
                                  d_seq = seq(0, 1, by = 0.1),
                                  nlambda = 100,
                                  standardize = TRUE,
                                  tuning_criterion = "Min_val",
                                  ebic_gamma = 0.5,
                                  lambda_min_ratio = 0.01,
                                  alpha_cores = 1,
                                  eps = 1e-4,
                                  maxit = 1e5,
                                  fast_exchangeable_data = TRUE,
                                  seed = NULL) {

  if (!is.null(seed)) set.seed(seed + rep_id)

  nonzero_id <- sample(J, strong_J)

  dat <- if (isTRUE(fast_exchangeable_data) &&
             identical(corrmat_type, "Exchangeable")) {
    block_sim_data_exchangeable_fast(
      n_train = n_train,
      n_val = n_val,
      n_test = n_test,
      pj = pj,
      J = J,
      strong_J = strong_J,
      rho_w = rho_w,
      rho_b = rho_b,
      eff_nonzero = eff_nonzero,
      snr = snr,
      nonzero_id = nonzero_id,
      signal_pattern = signal_pattern
    )
  } else {
    block_sim_data(
      n_train = n_train,
      n_val = n_val,
      n_test = n_test,
      pj = pj,
      J = J,
      strong_J = strong_J,
      rho_w = rho_w,
      rho_b = rho_b,
      eff_nonzero = eff_nonzero,
      corrmat_type = corrmat_type,
      snr = snr,
      nonzero_id = nonzero_id,
      signal_pattern = signal_pattern
    )
  }
  attr(dat$beta, "rho_w") <- dat$rho_w
  attr(dat$beta, "rho_b") <- dat$rho_b

  std_dat <- standardize_by_train(
    dat$X_train,
    dat$X_val,
    dat$X_test,
    dat$y_train,
    standardize
  )

  COV_X_test <- if (isTRUE(fast_exchangeable_data) &&
                    identical(corrmat_type, "Exchangeable")) {
    NULL
  } else {
    cov(dat$X_test)
  }

  res <- dplyr::bind_rows(
    fit_oracle(dat, std_dat, COV_X_test,
               tuning_criterion = "Oracle"),

    fit_grpreg_method("GLASSO", "grLasso",
                      dat, std_dat, COV_X_test,
                      nlambda, eps, maxit,
                      tuning_criterion, ebic_gamma, lambda_min_ratio),

    fit_adelie(dat, std_dat, COV_X_test,
               alpha_seq, nlambda, eps, maxit,
               alpha_cores = alpha_cores,
               tuning_criterion = tuning_criterion,
               ebic_gamma = ebic_gamma,
               lambda_min_ratio = lambda_min_ratio),

    fit_genet(dat, std_dat, COV_X_test,
              alpha_seq, nlambda, eps, maxit,
              alpha_cores = alpha_cores,
              tuning_criterion = tuning_criterion,
              ebic_gamma = ebic_gamma,
              lambda_min_ratio = lambda_min_ratio),

    fit_sglasso(dat, std_dat, COV_X_test,
                alpha_seq, d_seq, nlambda, eps, maxit,
                alpha_cores = alpha_cores,
                tuning_criterion = tuning_criterion,
                ebic_gamma = ebic_gamma,
                lambda_min_ratio = lambda_min_ratio),

    fit_grpreg_method("GSCAD", "grSCAD",
                      dat, std_dat, COV_X_test,
                      nlambda, eps, maxit,
                      tuning_criterion, ebic_gamma, lambda_min_ratio),

    fit_grpreg_method("GMCP", "grMCP",
                      dat, std_dat, COV_X_test,
                      nlambda, eps, maxit,
                      tuning_criterion, ebic_gamma, lambda_min_ratio)
  )

  res$rep <- rep_id
  res$setting <- NA_character_
  res$n_train <- n_train
  res$n_val <- n_val
  res$n_test <- n_test
  res$p <- pj * J
  res$pj <- pj
  res$J <- J
  res$strong_J <- strong_J
  res$rho_w <- rho_w
  res$rho_b <- rho_b
  res$snr <- snr
  res$corrmat_type <- corrmat_type
  res$signal_pattern <- signal_pattern
  res$standardize <- standardize
  res$fast_exchangeable_data <- isTRUE(fast_exchangeable_data)

  res
}

simulation_sglasso <- function(
    repeatnum = 10,
    seed = 2026,
    cores = 1,
    snr = 1.528,
    snr_grid = NULL,
    rho_b = 0.5,
    rho_b_grid = NULL,
    signal_pattern = "homogeneous",
    signal_pattern_grid = NULL,
    ...) {

  ############################################################
  # This wrapper runs the simulation over:
  #   - simulation repetitions,
  #   - SNR values,
  #   - between-group correlation values rho_b.
  #
  # If snr_grid is NULL, the function uses the single value snr.
  # If rho_b_grid is NULL, the function uses the single value rho_b.
  #
  # Examples:
  #
  # Single SNR and single rho_b:
  # simulation_sglasso(repeatnum = 10, snr = 1.528, rho_b = 0.5)
  #
  # Multiple SNR values:
  # simulation_sglasso(repeatnum = 10, snr_grid = get_snr_grid(), rho_b = 0.5)
  #
  # Multiple SNR and rho_b values:
  # simulation_sglasso(
  #   repeatnum = 10,
  #   snr_grid = get_snr_grid(),
  #   rho_b_grid = seq(0, 0.9, by = 0.1)
  # )
  ############################################################

  if (is.null(snr_grid)) {
    snr_grid <- snr
  }

  if (is.null(rho_b_grid)) {
    rho_b_grid <- rho_b
  }

  if (is.null(signal_pattern_grid)) {
    signal_pattern_grid <- signal_pattern
  }

  snr_grid <- as.numeric(snr_grid)
  rho_b_grid <- as.numeric(rho_b_grid)

  ############################################################
  # Create one task for every combination of:
  #   SNR x rho_b x signal_pattern_val x replicate.
  #
  # This avoids nested foreach loops, which are not supported by doRNG.
  ############################################################

  task_grid <- expand.grid(
    signal_pattern_val = signal_pattern_grid,
    snr_val = snr_grid,
    rho_b_val = rho_b_grid,
    rep_id = seq_len(repeatnum),
    stringsAsFactors = FALSE
  )

  ############################################################
  # Parallel version
  ############################################################

  if (cores > 1) {

    doParallel::registerDoParallel(cores = cores)
    doRNG::registerDoRNG(seed)

    out <- foreach::foreach(
      ii = seq_len(nrow(task_grid)),
      .combine = dplyr::bind_rows,
      .packages = required_packages
    ) %dopar% {

      fit_one_sim_replicate(
        rep_id = task_grid$rep_id[ii],
        seed = seed,
        snr = task_grid$snr_val[ii],
        rho_b = task_grid$rho_b_val[ii],
        signal_pattern = task_grid$signal_pattern_val[ii],
        ...
      )
    }

    foreach::registerDoSEQ()

  } else {

    ############################################################
    # Sequential version
    ############################################################

    out <- dplyr::bind_rows(lapply(seq_len(nrow(task_grid)), function(ii) {

      fit_one_sim_replicate(
        rep_id = task_grid$rep_id[ii],
        seed = seed,
        snr = task_grid$snr_val[ii],
        rho_b = task_grid$rho_b_val[ii],
        signal_pattern = task_grid$signal_pattern_val[ii],
        ...
      )

    }))
  }

  rownames(out) <- NULL

  out
}

############################################################
# 10. Paper settings and short tests
############################################################

get_paper_settings <- function() {
  data.frame(
    setting = c("LD-4", "LD-8", "LD-12", "HD-10", "HD-20", "HD-40"),

    n_train = c(500, 500, 500, 100, 100, 100),
    n_val   = c(200, 200, 200, 100, 100, 100),
    n_test  = c(200, 200, 200, 100, 100, 100),

    p        = c(100, 100, 100, 1000, 1000, 1000),
    strong_J = c(4, 8, 12, 10, 20, 40),
    J        = c(20, 20, 20, 200, 200, 200),
    pj       = c(5, 5, 5, 5, 5, 5),

    dimensionality = c("LD", "LD", "LD", "HD", "HD", "HD"),
    sparsity_level = c("Sparse", "Moderate", "Dense",
                       "Sparse", "Moderate", "Dense")
  )
}

get_snr_grid <- function() {
  exp(seq(log(0.05), log(6), length.out = 8))
}

run_short_sim_test <- function(tuning_criterion = "Min_val") {
  simulation_sglasso(
    repeatnum = 1,
    seed = 2026,
    cores = 1,
    n_train = 100,
    n_val = 50,
    n_test = 50,
    pj = 5,
    J = 20,
    strong_J = 4,
    rho_w = 0.9,
    rho_b = 0.5,
    eff_nonzero = 1,
    corrmat_type = "Exchangeable",
    snr = 1.528,
    signal_pattern = "weak_strong_mixed",
    alpha_seq = c(0.1, 0.5),
    d_seq = c(0.1, 0.5),
    nlambda = 20,
    standardize = TRUE,
    tuning_criterion = tuning_criterion,
    eps = 1e-4,
    maxit = 1e5
  )
}


############################################################
# 11. Summary tables
############################################################

make_pm <- function(mean_value, sd_value, digits = 3) {
  paste0(
    formatC(mean_value, digits = digits, format = "f"),
    " (",
    formatC(sd_value, digits = digits, format = "f"),
    ")"
  )
}

summarise_simulation_results <- function(res) {
  res %>%
    dplyr::group_by(setting, signal_pattern, rho_b, snr,
                    tuning_criterion, method) %>%
    dplyr::summarise(
      n_total = dplyr::n(),
      n_ok = sum(converged, na.rm = TRUE),
      time_mean = mean(time, na.rm = TRUE),
      time_sd = sd(time, na.rm = TRUE),
      risk_mean = mean(risk, na.rm = TRUE),
      risk_sd = sd(risk, na.rm = TRUE),
      pe_mse_mean = mean(pe_mse, na.rm = TRUE),
      pe_mse_sd = sd(pe_mse, na.rm = TRUE),
      MSE_y_mean = mean(MSE_y, na.rm = TRUE),
      MSE_y_sd = sd(MSE_y, na.rm = TRUE),
      MSE_beta_mean = mean(MSE_beta, na.rm = TRUE),
      MSE_beta_sd = sd(MSE_beta, na.rm = TRUE),
      pve_mean = mean(pve, na.rm = TRUE),
      pve_sd = sd(pve, na.rm = TRUE),
      FDR_mean = mean(FDR, na.rm = TRUE),
      FDR_sd = sd(FDR, na.rm = TRUE),
      selected_groups_mean = mean(selected_groups, na.rm = TRUE),
      selected_groups_sd = sd(selected_groups, na.rm = TRUE),
      selected_features_mean = mean(selected_features, na.rm = TRUE),
      selected_features_sd = sd(selected_features, na.rm = TRUE),
      alpha_median = median(alpha, na.rm = TRUE),
      lambda_median = median(lambda, na.rm = TRUE),
      d_median = median(d, na.rm = TRUE),
      .groups = "drop"
    )
}

make_revision_paper_table <- function(summary_results) {
  summary_results %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      Success = paste0(n_ok, "/", n_total),
      Runtime = make_pm(time_mean, time_sd, 3),
      Risk = make_pm(risk_mean, risk_sd, 4),
      Test_MSE = make_pm(pe_mse_mean, pe_mse_sd, 4),
      MSE_y = make_pm(MSE_y_mean, MSE_y_sd, 4),
      MSE_beta = make_pm(MSE_beta_mean, MSE_beta_sd, 4),
      PVE = make_pm(pve_mean, pve_sd, 4),
      FDR = make_pm(FDR_mean, FDR_sd, 4),
      Selected_Groups = make_pm(selected_groups_mean,
                                selected_groups_sd, 2),
      Selected_Features = make_pm(selected_features_mean,
                                  selected_features_sd, 2)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      setting, signal_pattern, rho_b, snr,
      tuning_criterion, method,
      Success, Runtime, Risk, Test_MSE,
      MSE_y, MSE_beta, PVE, FDR,
      Selected_Groups, Selected_Features,
      alpha_median, lambda_median, d_median
    )
}


############################################################
# 12. Scalability experiment
############################################################
# This function addresses the reviewer request for runtime comparisons
# across increasing p and J.  It records only runtimes and is separated
# from the main performance simulation for efficiency.

run_scalability_experiment <- function(p_grid = c(1000, 5000, 10000),
                                       J_grid = c(200, 500, 1000),
                                       pj = 5,
                                       n_train = 100,
                                       n_val = 0,
                                       n_test = 0,
                                       strong_frac = 0.10,
                                       rho_w = 0.9,
                                       rho_b = 0.5,
                                       snr = 1.528,
                                       repeatnum = 5,
                                       seed = 2026,
                                       alpha_fixed = 0.5,
                                       grpreg_alpha_fixed = 1,
                                       d_fixed = 0.5,
                                       nlambda = 50,
                                       standardize = FALSE,
                                       lambda_min_ratio = 0.01,
                                       eps = 1e-4,
                                       maxit = 1e5,
                                       signal_pattern = "weak_strong_mixed",
                                       checkpoint_dir = NULL,
                                       cores = 1,
                                       fast_exchangeable_data = TRUE) {

  set.seed(seed)
  cores <- max(1L, as.integer(cores))
  if (!is.null(checkpoint_dir)) {
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  }

  progress_msg <- function(...) {
    cat(
      sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      sprintf(...),
      "\n",
      sep = ""
    )
    flush.console()
  }

  make_task_key <- function(p_value, J_value, rep_id) {
    paste(p_value, J_value, rep_id, sep = "__")
  }

  save_scalability_checkpoint <- function(task_result,
                                          axis_name,
                                          p_value,
                                          J_value,
                                          rep_id,
                                          task_id,
                                          reused_duplicate = FALSE) {
    if (is.null(checkpoint_dir)) {
      return(task_result)
    }

    task_result$task_id <- task_id
    task_result$reused_duplicate <- reused_duplicate
    task_result$checkpoint_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

    checkpoint_file <- file.path(
      checkpoint_dir,
      sprintf(
        "task_%03d_axis_%s_p_%d_J_%d_rep_%03d.rds",
        task_id,
        axis_name,
        p_value,
        J_value,
        rep_id
      )
    )
    saveRDS(task_result, checkpoint_file)
    progress_msg("Saved scalability checkpoint: %s", checkpoint_file)
    invisible(task_result)
  }

  run_one_runtime <- function(axis_name, p_value, J_value, rep_id) {

    strong_J <- max(1, round(strong_frac * J_value))

    dat <- if (isTRUE(fast_exchangeable_data)) {
      block_sim_data_exchangeable_fast(
        n_train = n_train,
        n_val = n_val,
        n_test = n_test,
        pj = pj,
        J = J_value,
        strong_J = strong_J,
        rho_w = rho_w,
        rho_b = rho_b,
        snr = snr,
        nonzero_id = sample(J_value, strong_J),
        signal_pattern = signal_pattern
      )
    } else {
      block_sim_data(
        n_train = n_train,
        n_val = n_val,
        n_test = n_test,
        pj = pj,
        J = J_value,
        strong_J = strong_J,
        rho_w = rho_w,
        rho_b = rho_b,
        snr = snr,
        nonzero_id = sample(J_value, strong_J),
        signal_pattern = signal_pattern
      )
    }

    std_dat <- standardize_by_train(
      dat$X_train,
      dat$X_val,
      dat$X_test,
      dat$y_train,
      standardize
    )

    time_method <- function(expr, denom = 1) {
      t0 <- proc.time()
      err <- NULL
      ok <- tryCatch(
        {
          force(expr)
          TRUE
        },
        error = function(e) {
          err <<- conditionMessage(e)
          FALSE
        }
      )
      list(
        runtime = as.numeric((proc.time() - t0)[3]) / denom,
        success = ok,
        error = if (ok) NA_character_ else err
      )
    }

    t_glasso <- time_method(grpreg::grpreg(
      std_dat$X_train, std_dat$y_train,
      group = dat$group,
      penalty = "grLasso",
      alpha = grpreg_alpha_fixed,
      max.iter = maxit,
      eps = eps,
      lambda.min = lambda_min_ratio,
      nlambda = nlambda
    ))

    t_adelie <- time_method(adelie::grpnet(
      X = std_dat$X_train,
      glm = adelie::glm.gaussian(std_dat$y_train),
      groups = make_adelie_group_starts(dat$group),
      alpha = alpha_fixed,
      standardize = !isTRUE(std_dat$standardized),
      intercept = TRUE,
      lmda_path_size = nlambda,
      min_ratio = lambda_min_ratio,
      tol = eps,
      max_iters = as.integer(maxit),
      screen_rule = "strong",
      progress_bar = FALSE
    ))

    t_genet <- time_method(grpnet::grpnet(
      std_dat$X_train,
      std_dat$y_train,
      dat$group,
      nlambda = nlambda,
      alpha = alpha_fixed,
      thresh = eps,
      maxit = maxit,
      lambda.min.ratio = lambda_min_ratio
    ))

    t_sglasso <- time_method(call_sglasso(
      X = std_dat$X_train,
      Y = std_dat$y_train,
      group = dat$group,
      nlambda = nlambda,
      alpha = alpha_fixed,
      d = d_fixed,
      eps = eps,
      max_iter = maxit,
      standardize = !isTRUE(std_dat$standardized),
      screen = "SSR_fast",
      transform = "lazy",
      lambda.min.ratio = lambda_min_ratio
    ))

    t_gscad <- time_method(grpreg::grpreg(
      std_dat$X_train, std_dat$y_train,
      group = dat$group,
      penalty = "grSCAD",
      alpha = grpreg_alpha_fixed,
      max.iter = maxit,
      eps = eps,
      lambda.min = lambda_min_ratio,
      nlambda = nlambda
    ))

    t_gmcp <- time_method(grpreg::grpreg(
      std_dat$X_train, std_dat$y_train,
      group = dat$group,
      penalty = "grMCP",
      alpha = grpreg_alpha_fixed,
      max.iter = maxit,
      eps = eps,
      lambda.min = lambda_min_ratio,
      nlambda = nlambda
    ))

    data.frame(
      axis = axis_name,
      p = p_value,
      J = J_value,
      rep_id = rep_id,
      method = c("GLASSO", "ADELIE", "GENET",
                 "SGLASSO", "GSCAD", "GMCP"),
      alpha = c(grpreg_alpha_fixed, alpha_fixed, alpha_fixed,
                alpha_fixed, grpreg_alpha_fixed, grpreg_alpha_fixed),
      d = c(NA_real_, NA_real_, NA_real_,
            d_fixed, NA_real_, NA_real_),
      n_paths = 1,
      runtime = c(t_glasso$runtime, t_adelie$runtime, t_genet$runtime,
                  t_sglasso$runtime, t_gscad$runtime, t_gmcp$runtime),
      success = c(t_glasso$success, t_adelie$success, t_genet$success,
                  t_sglasso$success, t_gscad$success, t_gmcp$success),
      error = c(t_glasso$error, t_adelie$error, t_genet$error,
                t_sglasso$error, t_gscad$error, t_gmcp$error)
    )
  }

  p_tasks <- expand.grid(
    p = p_grid,
    rep_id = seq_len(repeatnum),
    KEEP.OUT.ATTRS = FALSE
  )
  p_tasks$J <- as.integer(p_tasks$p / pj)
  p_tasks$axis <- "p"

  j_tasks <- expand.grid(
    J = J_grid,
    rep_id = seq_len(repeatnum),
    KEEP.OUT.ATTRS = FALSE
  )
  j_tasks$p <- as.integer(j_tasks$J * pj)
  j_tasks$axis <- "J"

  all_tasks <- rbind(
    p_tasks[, c("axis", "p", "J", "rep_id")],
    j_tasks[, c("axis", "p", "J", "rep_id")]
  )
  all_tasks$task_id <- seq_len(nrow(all_tasks))
  all_tasks$key <- mapply(
    make_task_key,
    all_tasks$p,
    all_tasks$J,
    all_tasks$rep_id,
    USE.NAMES = FALSE
  )

  fit_tasks <- all_tasks[!duplicated(all_tasks$key), c("key", "p", "J", "rep_id", "task_id")]
  row.names(fit_tasks) <- NULL

  progress_msg(
    "Scalability task plan: output_tasks=%d | unique_fit_tasks=%d | cores=%d",
    nrow(all_tasks),
    nrow(fit_tasks),
    cores
  )

  run_fit_task <- function(ii) {
    fit_task <- fit_tasks[ii, , drop = FALSE]
    task_rows <- all_tasks[all_tasks$key == fit_task$key[[1]], , drop = FALSE]
    task_start <- proc.time()[3]

    progress_msg(
      "[%d/%d] Scalability fit started | p=%d | J=%d | rep=%d",
      ii,
      nrow(fit_tasks),
      fit_task$p[[1]],
      fit_task$J[[1]],
      fit_task$rep_id[[1]]
    )

    set.seed(seed + fit_task$task_id[[1]])
    fit_result <- run_one_runtime(
      axis_name = task_rows$axis[[1]],
      p_value = fit_task$p[[1]],
      J_value = fit_task$J[[1]],
      rep_id = fit_task$rep_id[[1]]
    )

    output <- lapply(seq_len(nrow(task_rows)), function(jj) {
      task_result <- fit_result
      task_result$axis <- task_rows$axis[[jj]]
      save_scalability_checkpoint(
        task_result,
        axis_name = task_rows$axis[[jj]],
        p_value = task_rows$p[[jj]],
        J_value = task_rows$J[[jj]],
        rep_id = task_rows$rep_id[[jj]],
        task_id = task_rows$task_id[[jj]],
        reused_duplicate = jj > 1L
      )
    })

    progress_msg(
      "[%d/%d] Scalability fit finished | p=%d | J=%d | rep=%d | output_tasks=%d | elapsed=%.2f sec",
      ii,
      nrow(fit_tasks),
      fit_task$p[[1]],
      fit_task$J[[1]],
      fit_task$rep_id[[1]],
      nrow(task_rows),
      proc.time()[3] - task_start
    )

    dplyr::bind_rows(output)
  }

  if (cores > 1L && nrow(fit_tasks) > 1L) {
    worker_count <- min(cores, nrow(fit_tasks))
    progress_msg("Starting parallel scalability fits with %d workers", worker_count)
    doParallel::registerDoParallel(cores = worker_count)
    doRNG::registerDoRNG(seed)

    out <- foreach::foreach(
      ii = seq_len(nrow(fit_tasks)),
      .combine = dplyr::bind_rows,
      .packages = required_packages,
      .export = c(
        "fit_tasks",
        "all_tasks",
        "seed",
        "run_fit_task",
        "run_one_runtime",
        "save_scalability_checkpoint",
        "progress_msg",
        "block_sim_data",
        "standardize_by_train",
        "make_adelie_group_starts"
      )
    ) %dopar% {
      run_fit_task(ii)
    }

    foreach::registerDoSEQ()
  } else {
    progress_msg("Starting serial scalability fits")
    out <- dplyr::bind_rows(
      lapply(seq_len(nrow(fit_tasks)), run_fit_task)
    )
  }

  out
}

summarise_scalability <- function(scalability_results) {
  scalability_results %>%
    dplyr::group_by(axis, p, J, method, alpha, d) %>%
    dplyr::summarise(
      n_ok = sum(success, na.rm = TRUE),
      n_total = dplyr::n(),
      median_runtime = median(runtime[success], na.rm = TRUE),
      mean_runtime = mean(runtime[success], na.rm = TRUE),
      sd_runtime = sd(runtime[success], na.rm = TRUE),
      .groups = "drop"
    )
}
