lambda_sglasso <- function(X, y, group, alpha, lambda.min, nlambda) {
  n <- length(y)
  # The C++ lambda_max routine expects contiguous group sizes and cumulative
  # group boundaries on the transformed design used by the solver.
  K <- table(group)
  K1 <- c(0, cumsum(K))

  # Start the path at lambda.max, where the null Gaussian model is active.
  r <- y - mean(y)
  zmax <- lambda_max_c(X,r,K,K1)/n

  # alpha=0 is a ridge-only edge case; use a tiny positive value so the
  # logarithmic lambda path is still well-defined for the sparse penalty part.
  if(alpha == 0) alpha=alpha+0.001
  lambda.max <- zmax/alpha
  lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), length=nlambda))
  return(list(lambda,lambda.max))
}
