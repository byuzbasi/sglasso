#' Log-likelihood for sglasso
#'
#' @param object A fitted \code{sglasso} object.
#' @param ... Additional arguments.
#'
#' @return Log-likelihood value.
#'
#' @export
logLik.sglasso <- function(object, ...) {
  if (identical(object$family, "binomial")) {
    return(structure(-object$deviance / 2, df = object$df,
      nobs = object$n, class = "logLik"))
  }
  n <- as.integer(object$n)
  df <- object$df
  RSS <- object$deviance
  l <- -n/2 * (log(2*pi) + log(RSS) - log(n)) - n/2
  df <- df + 1
  structure(l, df=df, nobs=n, class='logLik')
}
