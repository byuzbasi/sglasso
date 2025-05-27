#' @method coef cv.sglasso
#' @export
coef.cv.sglasso <- function(object, 
                            lambda=c("lambda.1se","lambda.min"), 
                            d=c("d.min"),
                            which=1:length(object$lambda), ...) {
  lmin = match.arg(lambda, c("lambda.min", "lambda.1se"))
  dmin = d
  predict(object$fit, lambda=object[[lmin]], d=object[[dmin]], which=which, type = "coefficients", ...)
}
