#' Extract coefficients from a fitted sglasso object
#'
#' @param object A fitted \code{sglasso} object.
#' @param lambda Lambda values for coefficient extraction.
#' @param d Scaling parameter values.
#' @param which Indices of lambda values.
#' @param drop Logical; should singleton dimensions be dropped?
#' @param ... Additional arguments for compatibility.
#'
#' @return Estimated coefficient array or matrix.
#'
#' @export
transform_sglasso_betas <- function(object) {
  if (!is.null(object$betas)) {
    return(object$betas)
  }
  if (is.null(object$raw_betas)) {
    stop("No coefficient path is stored in this sglasso object.", call. = FALSE)
  }

  Xtilde <- object$Xs
  lambda <- object$lambda
  d <- object$d
  raw <- object$raw_betas
  p <- dim(raw)[1]
  P <- length(object$group)
  nlambda <- length(lambda)
  nd <- length(d)

  beta_mat <- array(
    NA_real_,
    dim = c(P + 1, nlambda, nd),
    dimnames = list(c("intercept", paste("X", seq_len(P), sep = "")),
                    round(lambda, 4), d)
  )
  for (a in seq_len(nd)) {
    b <- rbind(object$y_mean, matrix(raw[, , a], nrow = p))
    bu <- unorthogonalize(b, Xtilde$X, Xtilde$g)
    beta_mat[, , a] <- unstandardize(bu, Xtilde)
  }
  beta_mat
}

#' @export
coef.sglasso <- function(object, lambda, d, which=1:length(object$lambda), drop=TRUE, ...) {
  if (missing(lambda)) lambda <- object$lambda[which]
  if (missing(d)) d <- object$d
  if (any(lambda > max(object$lambda) | lambda < min(object$lambda))) stop('lambda must lie within the range of the fitted coefficient path', call.=FALSE)
  if (any(d > max(object$d) | d < min(object$d))) stop('d must lie within the range of the fitted coefficient path', call.=FALSE)
  if (length(object$lambda) == 1L) {
    rl <- rep(1L, length(lambda))
  } else {
    ind_l <- approx(rev(object$lambda), rev(seq_along(object$lambda)), lambda)$y
    rl <- pmin(pmax(ceiling(ind_l), 1), length(object$lambda))
  }
  if (length(object$d) == 1L) {
    rd <- rep(1L, length(d))
  } else {
    ind_d <- approx(object$d, seq_along(object$d), d)$y
    rd <- pmin(pmax(ceiling(ind_d), 1), length(object$d))
  }
  beta_path <- transform_sglasso_betas(object)
  if (length(dim(beta_path)) == 3) {
    beta <- beta_path[, rl, rd, drop=FALSE]
  }else beta <- beta_path[, , rd, drop=FALSE]
  if (drop) return(drop(beta)) else return(beta)
}
