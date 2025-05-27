#' Extract coefficients from a sglasso object
#'
#' @method coef sglasso
#' @export
coef.sglasso <- function(object, lambda, d, which=1:length(object$lambda), drop=TRUE) {
  if (any(lambda > max(object$lambda) | lambda < min(object$lambda))) stop('lambda must lie within the range of the fitted coefficient path', call.=FALSE)
  if (any(d > max(object$d) | d < min(object$d))) stop('lambda must lie within the range of the fitted coefficient path', call.=FALSE)
  ind_l <- approx(object$lambda, seq(object$lambda), lambda)$y
  ll <- floor(ind_l)
  rl <- ceiling(ind_l)
  wl <- ind_l %% 1
  ind_d <- approx(object$d, seq(object$d), d)$y
  ld <- floor(ind_d)
  rd <- ceiling(ind_d)
  wd <- ind_d %% 1
  if (length(dim(object$betas)) == 3) {
    beta <- object$betas[, rl, rd, drop=FALSE]
  }else beta <- object$betas[, , rd, drop=FALSE]
  if (drop) return(drop(beta)) else return(beta)
}  