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
    newx = NULL,
    type = c("response", "coefficients", "vars", "groups"),
    lambda,
    d,
    s = c("ALL", "opt"),
    opt_beta = NULL,
    which = 1:length(object$lambda),
    drop = TRUE,
    ...
) {
  
  type <- match.arg(type)
  s <- match.arg(s)
  
  if (missing(lambda)) lambda <- object$lambda
  if (missing(d)) d <- object$d

  if (type == "response" && is.null(newx)) {
    stop("newx is required when type = 'response'", call. = FALSE)
  }

  if (!is.null(newx)) newx <- as.matrix(newx)
  
  # ------------------------------------------------------------
  # Special case: final fit with a single lambda and single d
  # ------------------------------------------------------------
  
  if (length(object$lambda) == 1 && length(object$d) == 1) {
    
    beta <- coef(object, drop = FALSE)
    
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

    if (type == "vars") {
      return(which(beta[-1] != 0))
    }

    if (type == "groups") {
      return(unique(object$group[beta[-1] != 0]))
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
      if (is.null(opt_beta)) {
        stop("opt_beta is required when s = 'opt'", call. = FALSE)
      }
      return(opt_beta)
    }
  }
  
  if (type == "vars") {
    if (s == "opt") {
      if (is.null(opt_beta)) {
        stop("opt_beta is required when s = 'opt'", call. = FALSE)
      }
      return(which(opt_beta[-1] != 0))
    }

    out <- apply(
      beta[-1, , , drop = FALSE] != 0,
      c(2, 3),
      which
    )
    if (drop) return(drop(out)) else return(out)
  }
  
  if (type == "groups") {
    if (s == "opt") {
      if (is.null(opt_beta)) {
        stop("opt_beta is required when s = 'opt'", call. = FALSE)
      }
      return(unique(object$group[opt_beta[-1] != 0]))
    }

    out <- apply(
      beta[-1, , , drop = FALSE] != 0,
      c(2, 3),
      function(x) unique(object$group[x])
    )
    if (drop) return(drop(out)) else return(out)
  }
  
  if (type == "response") {
    if (s == "ALL") {
      y_hat <- array(
        NA_real_,
        dim = c(nrow(newx), length(lambda), length(d))
      )

      x_int <- cbind(1, newx)
      for (di in seq_along(d)) {
        y_hat[, , di] <- x_int %*% beta[, , di, drop = FALSE][, , 1]
      }

      return(y_hat)
    } else {
      if (is.null(opt_beta)) {
        stop("opt_beta is required when s = 'opt'", call. = FALSE)
      }
      return(as.numeric(cbind(1, newx) %*% opt_beta))
    }
  }
}
