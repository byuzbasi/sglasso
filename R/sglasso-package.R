#' sglasso
#'
#' Scaled group lasso methods.
#'
#' @useDynLib sglasso, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom graphics abline arrows axis mtext points
#' @importFrom stats AIC BIC approx coef glm logLik model.matrix predict relevel sd
#' @importFrom Matrix bdiag
"_PACKAGE"