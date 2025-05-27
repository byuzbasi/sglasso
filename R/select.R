#' @rdname select
#' @export
select <- function(obj,...) UseMethod("select")


#' Extract coefficients from a sglasso object
#'
#' @method select sglasso
#' @export
select.sglasso <- function(obj, criterion=c("BIC","AIC","GCV","AICc","EBIC")) {
  criterion <- match.arg(criterion)
  length_lambda <- length(obj$lambda)
  length_d <- length(obj$d)
  ll <- logLik(obj)
  df <- as.double(attr(ll,"df"))
  d <- dim(obj$beta)
  p <- d[1] - 1
  j <-  df - 2
  IC <- switch(criterion,
               AIC = AIC(ll),
               BIC = BIC(ll),
               GCV = matrix((1/obj$n) * (-2) * as.double(ll) / (1-df/obj$n)^2,length_lambda,length_d,byrow = F),
               AICc = AIC(ll) + 2*df*(df+1)/(obj$n-df-1),
               EBIC = BIC(ll) + 2*(lgamma(p+1) - lgamma(j+1) - lgamma(p-j+1)))
  min_ind <- which(IC == min(IC), arr.ind = TRUE)
  if(dim(min_ind)[1] > 1) min_ind = min_ind[1,]
  return(list(beta=obj$beta[,min_ind[1],min_ind[2]],
              lambda=obj$lambda[min_ind[1]],
              d = obj$d[min_ind[2]],
              df=obj$df[min_ind[1],min_ind[2]],
              IC=IC))
}