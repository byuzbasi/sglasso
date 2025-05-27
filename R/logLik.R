#' logLik method for sglasso
#' @method logLik sglasso
#' @export
logLik.sglasso <- function(object) {
  n <- as.integer(object$n)
  df <- object$df
  RSS <- object$deviance
  l <- -n/2 * (log(2*pi) + log(RSS) - log(n)) - n/2
  df <- df + 1
  structure(l, df=df, nobs=n, class='logLik')
}