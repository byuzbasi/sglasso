#' @title Extract coefficients from a cv.sglasso object.
#' @param object Fitted cv.sglasso object.
#' @param ... Additional arguments for compatibility.

#' @return An object of class cv.sglasso, a list with entries
#' 
#' @description 
#' Extract coefficients from a cv.sglasso object.
#' 
#' @seealso 
#' \code{\link{sglasso}}
#' 
#' @examples 
#' data(Birthwt,package = "grpreg")
#' X <- Birthwt$X
#' group <- Birthwt$group
#' Y <- Birthwt$bwt
#' CVsglasso_fit <- cv.sglasso(X,Y,group)
#' @export
cv.sglasso <- function(X ,Y, group=1:ncol(X), 
                       lambda, nlambda = 100,
                       d, nd = 11,  
                       alpha = 0.5,
                       nfolds=10,
                       fold,
                       beta_start= NULL,
                       family= "gaussian", bilevel = F, max_iter=1e+08, eps = 1e-04){
  # [ToDo] Fit Lasso on original data using fitLASSO
  fit_sglasso = sglasso(X ,Y, group=group, nlambda=nlambda, lambda=lambda, d=d,
                        nd=nd, alpha=alpha, beta_start=beta_start, family=family, 
                        bilevel = bilevel, max_iter=max_iter, eps=eps)
  XG       = fit_sglasso$Xs
  group    = XG$g
  X        = XG$X
  n        = nrow(X)
  Y        = fit_sglasso$Ys
  lambda   = fit_sglasso$lambda
  d        = fit_sglasso$d
  betas    = fit_sglasso$betas
  alpha    = fit_sglasso$alpha
  
  
  
  # [ToDo] If fold is NULL, split the data randomly into k folds. If fold is not NULL, split the data according to supplied fold.
  if (missing(fold)) {
    fold <- sample(rep(seq(nfolds), length = n)) 
    } else nfolds <- max(fold)

  
  
  
  # errors container of size k by len(lambda_sec)
  #errors = matrix(NA, length(unique(fold)), length(lambda_seq))
  errors = array(NA, dim=c(length(lambda), length(d),length(unique(fold))))
  
  # [ToDo] Calculate LASSO on each fold using fitLASSO, and perform any additional calculations needed for CV(lambda) and SE_CV(lambda
  for (j in unique(fold)) {
    # Xtrain and Ytrain for each fold
    Xtrain = X[fold != j, ]
    Ytrain = Y[fold != j]
    
    # Xtest and Ytest for each fold
    Xtest = X[fold == j, ]
    Ytest = Y[fold == j]
    # Fit Model j
    result_j = sglasso(Xtrain, Ytrain, group, alpha=alpha,
                       lambda=lambda, nlambda=length(lambda),
                       d=d, nd=length(d),beta_start= NULL,family= "gaussian", 
                       bilevel=bilevel, max_iter=max_iter,eps = eps)
    
    # Vectorize the calculation of losses per fold. store the results in a corresponding row in the errors_all matrix
    
    y_hat_matrix_1 =  matrix(rep(Ytest, each = length(result_j$lambda)), length(Ytest), byrow = TRUE)
    y_hat_matrix =  array(NA,dim=c(dim(Xtest)[1], nlambda,nd))
    for (a in 1:nd) {y_hat_matrix[,,a] <- y_hat_matrix_1}
    
    
    # subtract the intercept from the matrix
    #loss_intercept = sweep(y_hat_matrix, 2, result_beta0)
    
    # calculate X'beta in matrix form
    loss_xbeta <- array(NA, dim=c(dim(Xtest)[1], nlambda,nd))
    loss_xbeta[] <- apply(result_j$betas, 3, function(x) cbind(1,Xtest)%*%x)
    errors[,,j] = apply((y_hat_matrix-loss_xbeta)^2, 3, colMeans)
    
  }
  
  
  
  # calculate cvlambda
  #cvm = colMeans(errors)
  cve = apply(errors,c(1,2),mean)
  cvse <- apply(errors, c(1,2), sd)/sqrt(n)
  colnames(cve) <- colnames(cvse) <- d
  rownames(cve) <- rownames(cvse) <- round(lambda,5)
  min_ind <- which(cve == min(cve), arr.ind = TRUE)
  min_ind <- min_ind[1,]
  
  
  # Return output
  # Output from fitLASSO on the whole data
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  # fold - used splitting into folds from 1 to k (either as supplied or as generated in the beginning)
  # lambda_min - selected lambda based on minimal rule
  # lambda_1se - selected lambda based on 1SE rule
  # cvm - values of CV(lambda) for each lambda
  # cvse - values of SE_CV(lambda) for each lambda
  #return(list(lambda_seq = lambda_seq, beta_mat = rbind(beta0_vec,beta_mat), beta0_vec = beta0_vec, fold = fold, lambda_min = lambda_min, lambda_1se = lambda_1se, cvm = cvm, cvse = cvse))
  return(structure(list(lambda    = lambda, 
              d         = d,
              alpha     = alpha,
              betas     = betas, 
              beta_opt  = betas[,min_ind[1],min_ind[2]], 
              fold      = fold, 
              lambda.min= lambda[min_ind[1]],
              d.min     = d[min_ind[2]],
              cve       = cve,
              cvse      = cvse,
              fit       = fit_sglasso,
              min_ind   = min_ind),
              class = "cv.sglasso")
         )
}
