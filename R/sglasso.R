#' @title Extract coefficients from a cv. sglasso object.
#' @param object Fitted sglasso object.
#' @param ... Additional arguments for compatibility.

#' @return A vector of coefficients
#' 
#' @description 
#' Extract coefficients from a sglasso object.
#' 
#' @seealso 
#' \code{\link{sglasso}}
#' 
#' @examples 
#' data(Birthwt,package = "grpreg")
#' X <- Birthwt$X
#' group <- Birthwt$group
#' Y <- Birthwt$bwt
#' sglasso(X,Y,group)
#' @export
sglasso <- function(X ,Y, group=1:ncol(X), 
                    lambda, nlambda = 50,
                    d, nd = 11,  
                    alpha = 0.5, 
                    beta_start= NULL,
                    family= "gaussian", bilevel = F, max_iter=1e+08, eps = 1e-04,
                    lambda.min = 0.005,
                    dfmax=p, gmax=length(unique(group))){
  ###
  ##############################################################################
  if (nlambda < 2) 
    stop("nlambda must be at least 2", call. = FALSE)
  if (nd < 2) 
    stop("nd must be at least 2", call. = FALSE)
  if (alpha > 1 | alpha < 0) 
    stop("alpha must be in [0, 1]", call. = FALSE)
  if (nrow(X) != length(Y)) 
    stop("X and Y do not have the same number of observations", call. = FALSE)
  
  ##############################################################################
  P <- ncol(X)
  group <- factor(group)
  K <- table(group)
  group.multiplier <- sqrt(K)
  ##############################################################################
  Ytilde <- newY(Y, family)
  Xtilde <- newXG(X, group, group.multiplier, attr(Ytilde, "m"), bilevel)
  n <- length(Ytilde)
  p <- ncol(Xtilde$X)
  K <- as.integer(table(Xtilde$g))
  K0 <- as.integer(if (min(Xtilde$g) == 0) K[1] else 0)
  K1 <- as.integer(if (min(Xtilde$g) == 0) cumsum(K) else c(0,cumsum(K)))
  group.multiplier <- sqrt(K)
  
  
  ###
  if (missing(lambda)) {
    lambda_values <- lambda_sglasso(Xtilde$X,Ytilde,group,alpha,lambda.min,nlambda)
    lambda <- lambda_values[[1]]
    lambda_max <- lambda_values[[2]]
    user.lambda <- FALSE
  }
  else {
    lambda_max <- -1
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }
  ###
  if (missing(d)) {
    d = seq(0, 1, length = nd)
    d_max <- d[nd]
    user.d <- FALSE
  }
  else {
    d_max <- -1
    nd <- length(d)
    user.d <- TRUE
  }
  
  
  ####################################################################################
  
  #[ToDo]  Check for starting point beta_start. If none supplied, initialize with a vector of zeros. If supplied, check for compatibility with Xtilde in terms of p
  if(length(beta_start) == 0){
    beta_start = rep(0, p)
  }else{
    if(length(beta_start) != p){
      stop('The number of features in beta_init and X cannot be different')
    }
  }
  
  
  
  # [ToDo] Fit sglasso on a sequence of values using gd_sglasso_ssr 
  # (make sure the parameters carry over)
  sglasso_fit = gd_sglasso_ssr(Xtilde$X, Ytilde,lambda,lambda_max,d, 
                               alpha, K, K1, K0, group.multiplier, beta_start, max_iter, eps,
                               dfmax,gmax)
  #min_ind <- which(lasso_fit$fmin_vec == min(lasso_fit$fmin_vec), arr.ind = TRUE)
  
  
  beta_mat = array(NA,dim=c(P+1,nlambda,nd), 
                   dimnames = list(c("intercept",paste("X", 1:P, sep="")),round(lambda,4),d)) 
  ## add intercept term
  for (a in 1:nd) {
    b <- rbind(mean(Y), matrix(sglasso_fit$beta_mat[,,a], nrow = p))
    bu <- unorthogonalize(b, Xtilde$X, Xtilde$g)
    beta_mat[,,a] <- unstandardize(bu, Xtilde)
  }
  
  
  
  # [ToDo] Perform back scaling and centering to get original intercept and coefficient vector for each lambda
  #b <- rbind(mean(Y), matrix(lasso_fit$beta_mat[,,], nrow = p))
  #bu <- unorthogonalize(b, Xtilde$X, Xtilde$g)
  #beta <- unstandardize(bu, Xtilde)
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  obj = list(lambda = lambda, d = d, alpha=alpha,
             deviance = sglasso_fit$loss,
             n=n,df=sglasso_fit$df,
             Xs = Xtilde, Ys=Ytilde,
             iter=sglasso_fit$iter,
             betas = beta_mat,
             group=group)
  class(obj) = "sglasso"
  return(obj)
}

################################################################
################################################################
################################################################