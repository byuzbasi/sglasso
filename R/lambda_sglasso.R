lambda_sglasso <- function(X, y, group, alpha, lambda.min, nlambda) {
  n <- length(y)
  K <- table(group)
  K1 <- c(0, cumsum(K))
  fit <- glm(y~1, family="gaussian")
  ## Determine lambda.max
  r <- fit$residuals
  zmax <- lambda_max_c(X,r,K,K1)/n
  if(alpha == 0) alpha=alpha+0.001
  lambda.max <- zmax/alpha
  lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), length=nlambda))
  return(list(lambda,lambda.max))
}
