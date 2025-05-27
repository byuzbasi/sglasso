#' @method predict sglasso
#' @export
predict.sglasso <- function(object, newx, type=c("response", "coefficients", "vars", "groups"), 
                            lambda, d, s, opt_beta, which=1:length(object$lambda), drop=T,...) {
  if(missing(lambda)) lambda= object$lambda
  if(missing(d)) d= object$d
  beta <- coef(object, lambda=lambda, d=d, which=which, drop=F)
  if (type=="coefficients") {
    if(s=="ALL"){if (drop) return(drop(beta)) else return(beta)
      } else return(opt_beta)
  } 
  if (type=="vars") return(drop(apply(beta[-1, , , drop=T]!=0, 2, FUN=which)))
  if (type=="groups") return(drop(apply(beta[-1, , , drop=T]!=0, 2, function(x) unique(object$group[x]))))
  if(type=="response") {
    if(s=="ALL"){
      y_hat =  array(NA,dim=c(dim(newx)[1], length(lambda),length(d)))
      y_hat[] = apply(beta, 3, function(x) cbind(1,newx)%*%x)  
      return(y_hat)  
    } else return(cbind(1,newx)%*%opt_beta)
  }
}








