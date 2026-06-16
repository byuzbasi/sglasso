#' @method coef cv.sglasso
#' @export
coef.cv.sglasso <- function(object, 
                            lambda=c("lambda.min","lambda.1se"),
                            d=c("d.min"),
                            which=1:length(object$lambda), ...) {
  # Convert the cv-selected labels to numeric tuning values, then delegate to
  # coef.sglasso on the fitted path object.
  lmin = match.arg(lambda, c("lambda.min", "lambda.1se"))
  dmin = match.arg(d, c("d.min"))
  if (is.null(object[[lmin]])) {
    stop(lmin, " is not available in this cv.sglasso object", call. = FALSE)
  }
  coef(object$fit, lambda=object[[lmin]], d=object[[dmin]], which=which, ...)
}
