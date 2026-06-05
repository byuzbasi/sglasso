#' @method predict cv.sglasso
#' @export
predict.cv.sglasso <- function(object, newx = NULL, lambda, d, which=1:length(object$lambda),
                               s=c("ALL","opt"),
                               type=c("response", "coefficients", "vars", "groups"), ...) {
  # cv.sglasso stores the selected model separately, but all prediction modes
  # are implemented by predict.sglasso on the underlying fit.
  type <- match.arg(type)
  s <- match.arg(s)
  if(missing(lambda)){lambda = object$lambda}
  if(missing(d)){d = object$d}
  return(predict(object$fit, newx=newx, lambda=lambda, d=d, which=which, s=s, opt_beta = object$beta_opt, type=type, ...))
}


