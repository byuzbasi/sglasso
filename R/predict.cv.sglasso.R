#' @method predict cv.sglasso
#' @export
predict.cv.sglasso <- function(object, newx, lambda, d, which=1:length(object$lambda),
                               s=c("ALL","opt"),
                               type=c("response", "coefficients", "vars", "groups"), ...) {
  type <- match.arg(type)
  if(missing(lambda)){lambda = object$lambda}
  if(missing(d)){d = object$d}
  return(predict(object$fit, newx=newx, lambda=lambda, d=d, which=which, s=s, opt_beta = object$beta_opt, type=type, ...))
}





