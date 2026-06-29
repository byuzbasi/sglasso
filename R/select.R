#' Select tuning parameters
#'
#' @param obj An object.
#' @param ... Additional arguments.
#'
#' @export
select <- function(obj, ...) {
  UseMethod("select")
}


#' Select tuning parameters for sglasso
#'
#' @param obj A fitted \code{sglasso} object.
#' @param criterion Selection criterion.
#' @param ebic_gamma EBIC gamma parameter. Default is 0.5.
#' @param ebic_level EBIC penalty level. Use \code{"group"} for group-level
#'   EBIC or \code{"feature"} for feature-level EBIC.
#' @param tol Tolerance used to determine nonzero coefficients.
#' @param ... Additional arguments.
#'
#' @export
select.sglasso <- function(obj,
                           criterion = c("AIC", "BIC", "EBIC", "GCV", "AICc"),
                           ebic_gamma = 0.5,
                           ebic_level = c("group", "feature"),
                           tol = 1e-8,
                           ...) {
  
  criterion <- match.arg(criterion)
  ebic_level <- match.arg(ebic_level)
  
  length_lambda <- length(obj$lambda)
  length_d <- length(obj$d)
  
  ll <- logLik(obj)
  df <- as.numeric(attr(ll, "df"))
  ll_num <- as.numeric(ll)
  
  if (length(ll_num) != length_lambda * length_d) {
    stop("logLik length does not match lambda-by-d grid.")
  }
  
  if (length(df) != length_lambda * length_d) {
    stop("df length does not match lambda-by-d grid.")
  }
  
  beta_all <- coef(obj, drop = FALSE)
  d_obj <- dim(beta_all)
  p <- d_obj[1] - 1
  
  IC <- switch(
    criterion,
    
    AIC = {
      matrix(
        AIC(ll),
        nrow = length_lambda,
        ncol = length_d,
        byrow = FALSE
      )
    },
    
    BIC = {
      matrix(
        BIC(ll),
        nrow = length_lambda,
        ncol = length_d,
        byrow = FALSE
      )
    },
    
    GCV = {
      gcv_vec <- (-2 * ll_num / obj$n) /
        pmax((1 - df / obj$n)^2, .Machine$double.eps)
      
      matrix(
        gcv_vec,
        nrow = length_lambda,
        ncol = length_d,
        byrow = FALSE
      )
    },
    
    AICc = {
      aicc_vec <- as.numeric(AIC(ll)) +
        2 * df * (df + 1) /
        pmax(obj$n - df - 1, .Machine$double.eps)
      
      matrix(
        aicc_vec,
        nrow = length_lambda,
        ncol = length_d,
        byrow = FALSE
      )
    },
    
    EBIC = {
      
      base_bic <- as.numeric(BIC(ll))
      
      if (ebic_level == "feature") {
        
        j <- pmax(df - 2, 0)
        j <- pmin(j, p)
        
        ebic_penalty <- 2 * ebic_gamma * (
          lgamma(p + 1) -
            lgamma(j + 1) -
            lgamma(p - j + 1)
        )
        
      } else {
        
        if (is.null(obj$group)) {
          stop("Group-level EBIC requires obj$group.")
        }
        
        group <- obj$group
        J <- length(unique(group))
        
        beta_path <- beta_all[-1, , , drop = FALSE]
        
        selected_groups <- apply(
          beta_path,
          c(2, 3),
          function(b) {
            length(unique(group[abs(b) > tol]))
          }
        )
        
        g <- pmin(selected_groups, J)
        
        ebic_penalty <- 2 * ebic_gamma * (
          lgamma(J + 1) -
            lgamma(g + 1) -
            lgamma(J - g + 1)
        )
      }
      
      ebic_vec <- base_bic + as.numeric(ebic_penalty)
      
      matrix(
        ebic_vec,
        nrow = length_lambda,
        ncol = length_d,
        byrow = FALSE
      )
    }
  )
  
  IC[!is.finite(IC)] <- Inf
  
  min_ind <- which(IC == min(IC, na.rm = TRUE), arr.ind = TRUE)
  
  if (nrow(min_ind) > 1) {
    min_ind <- min_ind[1, , drop = FALSE]
  }
  
  list(
    beta = beta_all[, min_ind[1], min_ind[2]],
    lambda = obj$lambda[min_ind[1]],
    d = obj$d[min_ind[2]],
    df = obj$df[min_ind[1], min_ind[2]],
    IC = IC,
    criterion = criterion,
    ebic_gamma = ebic_gamma,
    ebic_level = ebic_level
  )
}
