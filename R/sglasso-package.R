#' sglasso: Scaled Group Lasso
#'
#' Tools for scaled group lasso estimation and cross-validation.
#'
#' @keywords internal
#' @useDynLib sglasso, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom graphics abline arrows axis mtext points
#' @importFrom stats AIC BIC approx coef glm logLik model.matrix predict relevel
#' @importFrom Matrix bdiag
"_PACKAGE"