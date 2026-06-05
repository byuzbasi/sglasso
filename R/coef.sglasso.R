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
coef.sglasso <- function(object, lambda, d, which=1:length(object$lambda), drop=TRUE, ...) {
  if (missing(lambda)) lambda <- object$lambda[which]
  if (missing(d)) d <- object$d
  if (any(lambda > max(object$lambda) | lambda < min(object$lambda))) stop('lambda must lie within the range of the fitted coefficient path', call.=FALSE)
  if (any(d > max(object$d) | d < min(object$d))) stop('d must lie within the range of the fitted coefficient path', call.=FALSE)
  ind_l <- approx(rev(object$lambda), rev(seq_along(object$lambda)), lambda)$y
  rl <- pmin(pmax(ceiling(ind_l), 1), length(object$lambda))
  ind_d <- approx(object$d, seq_along(object$d), d)$y
  rd <- pmin(pmax(ceiling(ind_d), 1), length(object$d))
  if (length(dim(object$betas)) == 3) {
    beta <- object$betas[, rl, rd, drop=FALSE]
  }else beta <- object$betas[, , rd, drop=FALSE]
  if (drop) return(drop(beta)) else return(beta)
}
