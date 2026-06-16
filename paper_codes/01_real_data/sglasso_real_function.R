# ============================================================
# Real-data benchmark for grouped penalized regression methods
# TRUBA/server version
# ============================================================

# Tuning parameters are selected by cross-validation.
# Final runtime, iteration count, and prediction metrics are computed
# from a non-cross-validated refit using the selected tuning parameters.
#
# This design avoids inflating computational metrics with hyperparameter
# search overhead and enables a fairer comparison of solver efficiency.

# ------------------------------------------------------------
# 1. Package management
# ------------------------------------------------------------

# grpreg     : group lasso, group SCAD, group MCP
# grpnet     : group elastic net benchmark method
# adelie     : efficient group elastic net solver
# sglasso    : scaled group lasso package, installed once before the job
# foreach    : parallel loop infrastructure
# doParallel : parallel backend
# doRNG      : reproducible parallel RNG
# tibble     : clean result tables
# remotes    : install sglasso from GitHub when needed

user_lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_lib, .libPaths()))
}
install_lib <- if (nzchar(user_lib)) user_lib else NULL

cran_packages <- c(
  "grpreg",
  "grpnet",
  "adelie",
  "foreach",
  "doParallel",
  "doRNG",
  "tibble",
  "remotes"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, lib = install_lib, repos = "https://cloud.r-project.org")
  }
}

if (!requireNamespace("sglasso", quietly = TRUE)) {
  remotes::install_github(
    "byuzbasi/sglasso",
    lib = install_lib,
    upgrade = "never",
    dependencies = TRUE
  )
}

if (!requireNamespace("sglasso", quietly = TRUE)) {
  stop("Package 'sglasso' could not be installed from GitHub.")
}

library(grpnet)
library(sglasso)
library(foreach)
library(doParallel)
library(doRNG)
library(tibble)

# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------

safe_num <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  as.numeric(x)[1]
}

safe_sum <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  sum(x, na.rm = TRUE)
}


# Count selected groups/features from fitted coefficient vector.
# Intercept, if present, is removed automatically when coef length is p + 1.
extract_beta <- function(fit, p, lambda = NULL, coef_fun = stats::coef) {
  
  # Direct extraction for grpnet package objects
  if (!is.null(fit$beta)) {
    b <- fit$beta
    b <- as.vector(b)
    
    if (length(b) >= p) {
      return(tail(as.numeric(b), p))
    }
  }
  
  # Direct extraction for sglasso package objects
  if (!is.null(fit$betas)) {
    b <- fit$betas
    b <- as.vector(b)
    
    if (length(b) >= p) {
      return(tail(as.numeric(b), p))
    }
  }
  
  # General fallback
  b <- tryCatch({
    if (!is.null(lambda)) {
      coef_fun(fit, lambda = lambda)
    } else {
      coef_fun(fit)
    }
  }, error = function(e) NULL)
  
  if (is.null(b)) {
    return(rep(NA_real_, p))
  }
  
  if (is.list(b)) {
    if (!is.null(b$beta)) {
      b <- b$beta
    } else if (!is.null(b$coef)) {
      b <- b$coef
    } else {
      b <- unlist(b, recursive = TRUE, use.names = FALSE)
    }
  }
  
  if (inherits(b, "dgCMatrix") || inherits(b, "matrix")) {
    b <- as.vector(b)
  }
  
  b <- suppressWarnings(as.numeric(b))
  b <- b[is.finite(b)]
  
  if (length(b) == p + 1) {
    b <- b[-1]
  }
  
  if (length(b) > p) {
    b <- tail(b, p)
  }
  
  if (length(b) < p) {
    return(rep(NA_real_, p))
  }
  
  b
}
selection_summary <- function(beta, group, tol = 1e-6) {
  if (all(is.na(beta))) {
    return(list(
      selected_features = NA_real_,
      selected_groups = NA_real_,
      sparsity = NA_real_
    ))
  }

  active_feature <- abs(beta) > tol
  selected_features <- sum(active_feature, na.rm = TRUE)
  selected_groups <- sum(tapply(active_feature, group, any), na.rm = TRUE)

  list(
    selected_features = selected_features,
    selected_groups = selected_groups,
    sparsity = 1 - selected_features / length(beta)
  )
}

# ------------------------------------------------------------
# 3. Main benchmark function
# ------------------------------------------------------------

real_slasso <- function(
    X,
    y,
    group,
    n = nrow(X),
    nrep = 100,
    ncores = {
      cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
      if (is.na(cpus) || cpus < 1L) 1L else cpus
    },
    seed = 2025,
    train_ratio = 0.70,
    nfolds = 5,
    maxit = 1e8,
    eps = 1e-5,
    lambda_min_ratio = 0.05,
    nlambda = 20,
    alpha_seq = seq(0.1, 0.9, length.out = 5),
    nd = 11
) {
  
  X <- as.matrix(X)
  y <- as.numeric(y)
  group_membership <- as.numeric(group)
  
  make_adelie_group_starts <- function(group) {
    as.integer(tapply(seq_along(group), group, min))
  }
  
  adelie_groups <- make_adelie_group_starts(group_membership)
  
  if (nrow(X) != length(y)) {
    stop("nrow(X) must be equal to length(y).")
  }
  
  if (ncol(X) != length(group_membership)) {
    stop("ncol(X) must be equal to length(group).")
  }
  
  if (n != nrow(X)) {
    stop("n must be equal to nrow(X).")
  }
  
  ntrain <- round(n * train_ratio)
  
  methods <- c(
    "GLASSO",
    "ADELIE",
    "GENET",
    "SGLASSO",
    "GSCAD",
    "GMCP"
  )
  
  nmethods <- length(methods)
  
  cl <- parallel::makeCluster(ncores)
  
  user_lib <- Sys.getenv("R_LIBS_USER")
  parallel::clusterExport(cl, "user_lib", envir = environment())
  
  parallel::clusterEvalQ(cl, {
    if (nzchar(user_lib)) {
      .libPaths(c(user_lib, .libPaths()))
    }
    
    library(grpreg)
    library(grpnet)
    library(adelie)
    library(sglasso)
    library(foreach)
    library(doParallel)
    library(doRNG)
    library(tibble)
    
    NULL
  })
  
  doParallel::registerDoParallel(cl)
  doRNG::registerDoRNG(seed)
  
  cat("\nParallel setup\n")
  cat("Requested cores   :", ncores, "\n")
  cat("Available cores   :", parallel::detectCores(), "\n")
  cat("Registered workers:", foreach::getDoParWorkers(), "\n\n")
  
  on.exit({
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
    cat("Parallel cluster stopped.\n")
  }, add = TRUE)
  
  results_from_foreach <- foreach(
    i = seq_len(nrep),
    .combine = rbind,
    .packages = c(
      "grpreg",
      "sglasso",
      "foreach",
      "doParallel",
      "doRNG"
    ),
    .export = c("safe_sum", "safe_num", "extract_beta", "selection_summary")
  ) %dopar% {
    
    results <- data.frame(
      rep = rep(i, nmethods),
      method = methods,
      time = NA_real_,
      iter = NA_real_,
      mse = NA_real_,
      mae = NA_real_,
      alpha = NA_real_,
      lambda = NA_real_,
      d = NA_real_,
      selected_groups = NA_real_,
      selected_features = NA_real_,
      sparsity = NA_real_
    )
    
    trainid <- sample.int(n, size = ntrain)
    testid <- setdiff(seq_len(n), trainid)
    foldid <- sample(rep(seq_len(nfolds), length.out = ntrain))
    
    # --------------------------------------------------------
    # 1. GLASSO: tuning and fit by cv.grpreg
    # --------------------------------------------------------
    
    counter <- 1
    
    # ---- CV for lambda selection ----
    
    cvmod <- grpreg::cv.grpreg(
      X = X[trainid, ],
      y = y[trainid],
      group = group_membership,
      fold = foldid,
      alpha = 1,
      penalty = "grLasso",
      nlambda = nlambda,
      lambda.min = lambda_min_ratio,
      eps = eps,
      max.iter = maxit,
      dfmax = ncol(X),
      gmax = length(unique(group_membership))
    )
    
    best_lambda <- cvmod$lambda.min
    
    # ---- Final fit without CV ----
    
    tic <- proc.time()
    
    fit <- grpreg::grpreg(
      X = X[trainid, ],
      y = y[trainid],
      group = group_membership,
      alpha = 1,
      penalty = "grLasso",
      lambda = best_lambda,
      eps = eps,
      max.iter = maxit,
      dfmax = ncol(X),
      gmax = length(unique(group_membership))
    )
    
    toc <- proc.time() - tic
    
    pred <- predict(fit, X[testid, ])
    
    results$time[counter] <- toc[3]
    results$iter[counter] <- safe_sum(fit$iter)
    results$mse[counter] <- mean((y[testid] - pred)^2)
    results$mae[counter] <- mean(abs(y[testid] - pred))
    results$alpha[counter] <- 1
    results$lambda[counter] <- best_lambda
    beta_hat <- extract_beta(fit, p = ncol(X))
    sel <- selection_summary(beta_hat, group_membership)
    results$selected_groups[counter] <- sel$selected_groups
    results$selected_features[counter] <- sel$selected_features
    results$sparsity[counter] <- sel$sparsity
    
    # --------------------------------------------------------
    # 2. ADELIE: select alpha by CV, then refit selected alpha
    # --------------------------------------------------------
    
    counter <- 2
    
    x_train <- X[trainid, ]
    y_train <- y[trainid]
    x_test  <- X[testid, ]
    
    # ---- CV for alpha/lambda selection ----
    
    adelie_fits <- lapply(alpha_seq, function(a) {
      
      adelie::cv.grpnet(
        X = x_train,
        glm = adelie::glm.gaussian(y = y_train),
        groups = adelie_groups,
        foldid = foldid,
        alpha = a,
        min_ratio = lambda_min_ratio,
        lmda_path_size = nlambda,
        n_threads = 1,
        tol = eps,
        max_iters = maxit,
        irls_tol = eps,
        irls_max_iters = maxit,
        newton_tol = eps,
        newton_max_iters = maxit,
        early_exit = FALSE
      )
      
    })
    
    adelie_cvm <- sapply(adelie_fits, function(fit) {
      min(fit$cvm, na.rm = TRUE)
    })
    
    best_ind <- which.min(adelie_cvm)
    
    best_alpha <- alpha_seq[best_ind]
    best_lambda <- adelie_fits[[best_ind]]$lambda.min
    
    if (!is.finite(best_lambda)) {
      lambda_path <- adelie_fits[[best_ind]]$lambda
      
      best_lambda <- min(
        lambda_path[is.finite(lambda_path) & lambda_path > 0],
        na.rm = TRUE
      )
    }
    
    
    # ---- Final fit without CV ----
    
    tic <- proc.time()
    
    fit <- adelie::grpnet(
      X = x_train,
      glm = adelie::glm.gaussian(y = y_train),
      groups = adelie_groups,
      alpha = best_alpha,
      lambda = best_lambda,
      n_threads = 1,
      tol = eps,
      max_iters = maxit,
      irls_tol = eps,
      irls_max_iters = maxit,
      newton_tol = eps,
      newton_max_iters = maxit,
      early_exit = FALSE
    )
    
    toc <- proc.time() - tic
    
    
    
    pred <- getS3method(
      f = "predict",
      class = "grpnet",
      envir = asNamespace("adelie")
    )(
      fit,
      newx = x_test,
      lambda = best_lambda
    )
    
    results$time[counter] <- toc[3]
    results$iter[counter] <- NA_real_
    results$mse[counter] <- mean((y[testid] - pred)^2)
    results$mae[counter] <- mean(abs(y[testid] - pred))
    results$alpha[counter] <- best_alpha
    results$lambda[counter] <- best_lambda
    beta_hat <- extract_beta(
      fit,
      p = ncol(X),
      lambda = best_lambda,
      coef_fun = getS3method("coef", "grpnet", envir = asNamespace("adelie"))
    )
    sel <- selection_summary(beta_hat, group_membership)
    results$selected_groups[counter] <- sel$selected_groups
    results$selected_features[counter] <- sel$selected_features
    results$sparsity[counter] <- sel$sparsity
    
    # --------------------------------------------------------
    # 3. GENET: select alpha by CV, then refit selected alpha
    # --------------------------------------------------------
    
    counter <- 3
    
    # ---- CV for alpha/lambda selection ----
    
    grpnet_fits <- lapply(alpha_seq, function(a) {
      
      grpnet::cv.grpnet(
        x = X[trainid, ],
        y = y[trainid],
        group = group_membership,
        foldid = foldid,
        type.measure = "mse",
        alpha = a,
        thresh = eps,
        maxit = maxit,
        lambda.min.ratio = lambda_min_ratio,
        verbose = FALSE
      )
      
    })
    
    grpnet_cvm <- sapply(grpnet_fits, function(fit) {
      min(fit$cvm, na.rm = TRUE)
    })
    
    best_ind <- which.min(grpnet_cvm)
    
    best_alpha <- alpha_seq[best_ind]
    best_lambda <- grpnet_fits[[best_ind]]$lambda.min
    
    # ---- Final fit without CV ----
    
    tic <- proc.time()
    
    fit <- grpnet::grpnet(
      x = X[trainid, ],
      y = y[trainid],
      group = group_membership,
      alpha = best_alpha,
      lambda = best_lambda,
      thresh = eps,
      maxit = maxit
    )
    
    toc <- proc.time() - tic
    
    pred <- getS3method(
      f = "predict",
      class = "grpnet",
      envir = asNamespace("grpnet")
    )(
      fit,
      newx = X[testid, ]
    )
    
    results$time[counter] <- toc[3]
    results$iter[counter] <- safe_sum(fit$npasses)
    results$mse[counter] <- mean((y[testid] - pred)^2)
    results$mae[counter] <- mean(abs(y[testid] - pred))
    results$alpha[counter] <- best_alpha
    results$lambda[counter] <- best_lambda
    beta_hat <- extract_beta(
      fit,
      p = ncol(X),
      coef_fun = getS3method("coef", "grpnet", envir = asNamespace("grpnet"))
    )
    sel <- selection_summary(beta_hat, group_membership)
    results$selected_groups[counter] <- sel$selected_groups
    results$selected_features[counter] <- sel$selected_features
    results$sparsity[counter] <- sel$sparsity
    
    # --------------------------------------------------------
    # 4. SGLASSO: select alpha and d by CV, then final fit
    # --------------------------------------------------------
    
    counter <- 4
    
    sglasso_fits <- lapply(alpha_seq, function(a) {
      sglasso::cv.sglasso(
        X = X[trainid, ],
        Y = y[trainid],
        group = group_membership,
        fold = foldid,
        nd = nd,
        lambda.min.ratio = lambda_min_ratio,
        alpha = a,
        eps = eps,
        max_iter = maxit
      )
    })
    
    sglasso_cve <- sapply(sglasso_fits, function(fit) {
      min(fit$cve, na.rm = TRUE)
    })
    
    best_ind <- which.min(sglasso_cve)
    best_alpha <- alpha_seq[best_ind]
    best_d <- sglasso_fits[[best_ind]]$d.min
    best_lambda <- sglasso_fits[[best_ind]]$lambda.min
    
    tic <- proc.time()
    
    fit <- sglasso::sglasso(
      X = X[trainid, ],
      Y = y[trainid],
      group = group_membership,
      d = best_d,
      alpha = best_alpha,
      lambda = best_lambda,
      eps = eps,
      max_iter = maxit
    )
    
    toc <- proc.time() - tic
    
    pred <- as.numeric(
      predict(
        fit,
        newx = X[testid, ],
        type = "response"
      )
    )
    
    
    results$time[counter] <- toc[3]
    results$iter[counter] <- safe_sum(fit$iter)
    results$mse[counter] <- mean((y[testid] - pred)^2)
    results$mae[counter] <- mean(abs(y[testid] - pred))
    results$alpha[counter] <- best_alpha
    results$lambda[counter] <- best_lambda
    results$d[counter] <- safe_num(best_d)
    beta_hat <- extract_beta(fit, p = ncol(X))
    sel <- selection_summary(beta_hat, group_membership)
    results$selected_groups[counter] <- sel$selected_groups
    results$selected_features[counter] <- sel$selected_features
    results$sparsity[counter] <- sel$sparsity
    
    # --------------------------------------------------------
    # 5. GSCAD: tuning and fit by cv.grpreg
    # --------------------------------------------------------
    
    counter <- 5
    
    cvmod <- grpreg::cv.grpreg(
      X = X[trainid, ],
      y = y[trainid],
      group = group_membership,
      fold = foldid,
      alpha = 1,
      penalty = "grSCAD",
      nlambda = nlambda,
      lambda.min = lambda_min_ratio,
      eps = eps,
      max.iter = maxit,
      dfmax = ncol(X),
      gmax = length(unique(group_membership))
    )
    
    best_lambda <- cvmod$lambda.min
    
    # ---- Final fit without CV ----
    
    tic <- proc.time()
    
    fit <- grpreg::grpreg(
      X = X[trainid, ],
      y = y[trainid],
      group = group_membership,
      alpha = 1,
      penalty = "grSCAD",
      lambda = best_lambda,
      eps = eps,
      max.iter = maxit,
      dfmax = ncol(X),
      gmax = length(unique(group_membership))
    )
    
    toc <- proc.time() - tic
    
    pred <- predict(fit, X[testid, ])
    
    results$time[counter] <- toc[3]
    results$iter[counter] <- safe_sum(fit$iter)
    results$mse[counter] <- mean((y[testid] - pred)^2)
    results$mae[counter] <- mean(abs(y[testid] - pred))
    results$alpha[counter] <- 1
    results$lambda[counter] <- safe_num(cvmod$lambda.min)
    beta_hat <- extract_beta(fit, p = ncol(X))
    sel <- selection_summary(beta_hat, group_membership)
    results$selected_groups[counter] <- sel$selected_groups
    results$selected_features[counter] <- sel$selected_features
    results$sparsity[counter] <- sel$sparsity
    
    # --------------------------------------------------------
    # 6. GMCP: tuning and fit by cv.grpreg
    # --------------------------------------------------------
    
    counter <- 6
    cvmod <- grpreg::cv.grpreg(
      X = X[trainid, ],
      y = y[trainid],
      group = group_membership,
      fold = foldid,
      alpha = 1,
      penalty = "grMCP",
      nlambda = nlambda,
      lambda.min = lambda_min_ratio,
      eps = eps,
      max.iter = maxit,
      dfmax = ncol(X),
      gmax = length(unique(group_membership))
    )
    
    best_lambda <- cvmod$lambda.min
    
    # ---- Final fit without CV ----
    
    tic <- proc.time()
    
    fit <- grpreg::grpreg(
      X = X[trainid, ],
      y = y[trainid],
      group = group_membership,
      alpha = 1,
      penalty = "grMCP",
      lambda = best_lambda,
      eps = eps,
      max.iter = maxit,
      dfmax = ncol(X),
      gmax = length(unique(group_membership))
    )
    
    toc <- proc.time() - tic
    
    pred <- predict(fit, X[testid, ])
    
    results$time[counter] <- toc[3]
    results$iter[counter] <- safe_sum(fit$iter)
    results$mse[counter] <- mean((y[testid] - pred)^2)
    results$mae[counter] <- mean(abs(y[testid] - pred))
    results$alpha[counter] <- 1
    results$lambda[counter] <- best_lambda
    beta_hat <- extract_beta(fit, p = ncol(X))
    sel <- selection_summary(beta_hat, group_membership)
    results$selected_groups[counter] <- sel$selected_groups
    results$selected_features[counter] <- sel$selected_features
    results$sparsity[counter] <- sel$sparsity
    
    results$method <- factor(results$method, levels = methods)
    
    results
  }
  
  return(results_from_foreach)
}
