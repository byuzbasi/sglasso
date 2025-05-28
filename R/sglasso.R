#' @title Fit a scaled group lasso regression path
#' 
#' @description 
#' Computes a scaled group lasso regularized linear models. 
#' 
#' 
#' @param X The design matrix, without an intercept.  \code{sglasso}
#' standardizes the data and includes an intercept by default.
#' @param y The response vector
#' @param group A vector describing the grouping of the coefficients.
#' @param family is "gaussian", not other option at this moment, depending on the response.
#' @param bilevel bi-level selection is not supported at this moment.
#' @param nlambda The number of \code{lambda} values.  Default is 50
#' @param lambda A user supplied sequence of \code{lambda} values.  Typically,
#' this is left unspecified, and the function automatically computes a grid of
#' lambda values that ranges uniformly on the log scale over the relevant range
#' of lambda values.
#' @param nd The number of \code{d} values.  Default is 11.
#' @param d The scale parameter between 0 and 1.
#' @param alpha Elastic Net tuning constant: the value must be between 0 and 1. Default is 0.5.
#' @param lambda.min The smallest value for \code{lambda}, as a fraction of
#' \code{lambda.max}.  Default is .0005.
#' @param eps Convergence threshhold.  The algorithm iterates until the BCD
#' for the change in linear predictors for each coefficient is less than
#' \code{eps}.  Default is \code{1e-4}.  
#' @param max.iter Maximum number of iterations (total across entire path).
#' Default is 1e+08.
#' @param dfmax Limit on the number of parameters allowed to be nonzero.  If
#' this limit is exceeded, the algorithm will exit early from the
#' regularization path.
#' @param gmax Limit on the number of groups allowed to have nonzero elements.
#' If this limit is exceeded, the algorithm will exit early from the
#' regularization path.
#' 
#' @return An object with S3 class \code{"sglasso"} containing:
#' \describe{
#' \item{beta}{The fitted matrix of coefficients.  The number of rows is equal
#' to the number of coefficients, and the number of columns is equal to
#' \code{nlambda}.}
#' \item{family}{Same as above.}
#' \item{group}{Same as above.}
#' \item{lambda}{The sequence of \code{lambda} values in the path.}
#' \item{alpha}{Same as above.}
#' \item{deviance}{A vector containing the deviance of the fitted model at each value of `lambda`.}
#' \item{n}{Number of observations.}
#' \item{penalty}{Same as above.}
#' \item{df}{A vector of length `nlambda` containing estimates of effective number of model parameters all the points along the regularization path.  For details on how this is calculated, see Breheny and Huang (2009).}
#' \item{iter}{A vector of length `nlambda` containing the number of iterations until convergence at each value of `lambda`.}
#' \item{group.multiplier}{A named vector containing the multiplicative constant applied to each group's penalty.}
#' }
#' 
#' 
#' @author Bahadir Yuzbasi
#' 
#' @seealso 
#' \code{\link{cv.sglasso}}
#' 
#' @examples 
#' data(GenAtHum,package = "sglasso")
#' X <- GenAtHum$X
#' y <- GenAtHum$y
#' group <- GenAtHum$group
#' n =  nrow(X)
#' p =  ncol(X)
#' set.seed(123)
#' fit <- sglasso(X,y,group, nlambda = 20, nd = 3)
#' select(fit,"EBIC")
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