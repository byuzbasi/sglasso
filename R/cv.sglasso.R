#' Cross-validation for scaled group lasso
#'
#' @description
#' Performs K-fold cross-validation for the scaled group lasso over a grid of
#' lambda and d values.
#'
#' @param X Design matrix.
#' @param Y Response vector.
#' @param group Group membership vector.
#' @param lambda Optional lambda sequence.
#' @param nlambda Number of lambda values.
#' @param d Optional d sequence.
#' @param nd Number of d values.
#' @param alpha Mixing parameter.
#' @param nfolds Number of cross-validation folds.
#' @param fold Optional fold assignment vector.
#' @param beta_start Optional initial coefficient vector.
#' @param family Model family. Currently "gaussian".
#' @param bilevel Logical; whether bilevel selection is used.
#' @param max_iter Maximum number of iterations.
#' @param eps Convergence tolerance.
#' @param standardize Logical flag indicating whether \code{X} should be
#' centered and scaled internally before fitting.
#' @param screen Screening rule passed to \code{\link{sglasso}}.
#' @param ... Additional arguments.
#'
#' @return An object of class \code{cv.sglasso}.
#'
#' @seealso \code{\link{sglasso}}
#'
#' @examples
#' data(Birthwt, package = "grpreg")
#' X <- Birthwt$X
#' group <- Birthwt$group
#' Y <- Birthwt$bwt
#' CVsglasso_fit <- cv.sglasso(X, Y, group)
#'
#' @export
cv.sglasso <- function(
    X,
    Y,
    group = 1:ncol(X),
    lambda,
    nlambda = 100,
    d,
    nd = 11,
    alpha = 0.5,
    nfolds = 10,
    fold,
    beta_start = NULL,
    family = "gaussian",
    bilevel = FALSE,
    max_iter = 1e8,
    eps = 1e-4,
    standardize = TRUE,
    screen = c("SSR", "none", "SSR_fast"),
    ...
) {
  screen <- screen[1L]
  screen <- match.arg(screen, choices = c("SSR", "none", "SSR_fast"))
  
  fit_sglasso <- sglasso(
    X = X,
    Y = Y,
    group = group,
    lambda = lambda,
    nlambda = nlambda,
    d = d,
    nd = nd,
    alpha = alpha,
    beta_start = beta_start,
    family = family,
    bilevel = bilevel,
    max_iter = max_iter,
    eps = eps,
    standardize = standardize,
    screen = screen
  )
  
  XG <- fit_sglasso$Xs
  group <- XG$g
  X <- XG$X
  Y <- fit_sglasso$Ys
  
  n <- nrow(X)
  lambda <- fit_sglasso$lambda
  d <- fit_sglasso$d
  betas <- coef(fit_sglasso, drop = FALSE)
  alpha <- fit_sglasso$alpha
  
  nlambda_eff <- length(lambda)
  nd_eff <- length(d)
  
  if (missing(fold)) {
    fold <- sample(rep(seq_len(nfolds), length.out = n))
  } else {
    nfolds <- max(fold)
  }
  
  errors <- array(
    NA_real_,
    dim = c(nlambda_eff, nd_eff, length(unique(fold)))
  )
  
  for (j in unique(fold)) {
    
    Xtrain <- X[fold != j, , drop = FALSE]
    Ytrain <- Y[fold != j]
    
    Xtest <- X[fold == j, , drop = FALSE]
    Ytest <- Y[fold == j]
    
    result_j <- sglasso(
      X = Xtrain,
      Y = Ytrain,
      group = group,
      alpha = alpha,
      lambda = lambda,
      nlambda = nlambda_eff,
      d = d,
      nd = nd_eff,
      beta_start = NULL,
      family = family,
      bilevel = bilevel,
      max_iter = max_iter,
      eps = eps,
      standardize = standardize,
      screen = screen
    )
    
    y_hat_array <- array(
      NA_real_,
      dim = c(nrow(Xtest), nlambda_eff, nd_eff)
    )
    
    result_betas <- coef(result_j, drop = FALSE)
    y_hat_array[] <- apply(
      result_betas,
      3,
      function(beta_mat) {
        cbind(1, Xtest) %*% beta_mat
      }
    )
    
    y_obs_array <- array(
      Ytest,
      dim = dim(y_hat_array)
    )
    
    errors[, , j] <- apply(
      (y_obs_array - y_hat_array)^2,
      c(2, 3),
      mean
    )
  }
  
  cve <- apply(errors, c(1, 2), mean, na.rm = TRUE)
  cvse <- apply(errors, c(1, 2), stats::sd, na.rm = TRUE) /
    sqrt(length(unique(fold)))
  
  colnames(cve) <- colnames(cvse) <- d
  rownames(cve) <- rownames(cvse) <- round(lambda, 5)
  
  min_ind <- which(cve == min(cve, na.rm = TRUE), arr.ind = TRUE)
  min_ind <- min_ind[1, ]
  
  out <- list(
    lambda = lambda,
    d = d,
    alpha = alpha,
    betas = betas,
    beta_opt = betas[, min_ind[1], min_ind[2]],
    fold = fold,
    lambda.min = lambda[min_ind[1]],
    d.min = d[min_ind[2]],
    cve = cve,
    cvse = cvse,
    fit = fit_sglasso,
    min_ind = min_ind
  )
  
  class(out) <- "cv.sglasso"
  out
}
