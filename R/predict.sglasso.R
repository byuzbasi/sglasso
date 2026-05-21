#' Predict method for sglasso objects
#'
#' @param object Fitted sglasso object.
#' @param newx New design matrix.
#' @param type Prediction type.
#' @param lambda Lambda values.
#' @param d Scaling parameter values.
#' @param s Selection mode.
#' @param opt_beta Optional coefficient vector.
#' @param which Indices of lambda values.
#' @param drop Logical; should dimensions be dropped?
#' @param ... Additional arguments.
#' @method predict sglasso
#' @export
predict.sglasso <- function(
    object,
    newx,
    type = c("response", "coefficients", "vars", "groups"),
    lambda,
    d,
    s = "ALL",
    opt_beta = NULL,
    which = 1:length(object$lambda),
    drop = TRUE,
    ...
) {
  
  type <- match.arg(type)
  
  if (missing(lambda)) lambda <- object$lambda
  if (missing(d)) d <- object$d
  
  newx <- as.matrix(newx)
  
  # ------------------------------------------------------------
  # Special case: final fit with a single lambda and single d
  # ------------------------------------------------------------
  
  if (length(object$lambda) == 1 && length(object$d) == 1) {
    
    beta <- object$beta
    
    if (length(dim(beta)) == 3) {
      beta <- beta[, 1, 1]
    }
    
    beta <- as.numeric(beta)
    
    if (type == "coefficients") {
      return(beta)
    }
    
    if (type == "response") {
      if (length(beta) == ncol(newx) + 1) {
        return(as.numeric(cbind(1, newx) %*% beta))
      } else {
        return(as.numeric(newx %*% beta))
      }
    }
  }
  
  # ------------------------------------------------------------
  # General path prediction
  # ------------------------------------------------------------
  
  beta <- coef(
    object,
    lambda = lambda,
    d = d,
    which = which,
    drop = FALSE
  )
  
  if (type == "coefficients") {
    if (s == "ALL") {
      if (drop) return(drop(beta)) else return(beta)
    } else {
      return(opt_beta)
    }
  }
  
  if (type == "vars") {
    return(drop(apply(beta[-1, , , drop = TRUE] != 0, 2, which)))
  }
  
  if (type == "groups") {
    return(drop(apply(
      beta[-1, , , drop = TRUE] != 0,
      2,
      function(x) unique(object$group[x])
    )))
  }
  
  if (type == "response") {
    if (s == "ALL") {
      y_hat <- array(
        NA,
        dim = c(nrow(newx), length(lambda), length(d))
      )
      
      y_hat[] <- apply(
        beta,
        3,
        function(x) cbind(1, newx) %*% x
      )
      
      return(y_hat)
    } else {
      return(as.numeric(cbind(1, newx) %*% opt_beta))
    }
  }
}