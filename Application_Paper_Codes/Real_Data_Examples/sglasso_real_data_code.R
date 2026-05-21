chooseCRANmirror(ind=11)

packages <- c("simstudy","data.table","grpnet","sglasso","grpreg","gglasso",
              "caret","dplyr","MASS","glmnet","Metrics","matrixcalc","readr",
              "foreach","doParallel","doFuture","doRNG","mltools","splines",
              "limma","umap")



# go over each package individually 
for (package in packages) { 
  # check whether the package is installed
  path <- system.file(package = package)
  # if it is not installed...
  if (path == "") {
    install.packages(package)
  } else { 
    cat(sprintf("%s is already installed\n", package))
  }
}

# after checking that all the required packages are installed, load 
# the packages
for (package in packages) { 
  library(package, character.only = TRUE) 
}


real_sglasso <- function(X,y,group,n,nrep){
### cv design
ntrain <- round(n * 0.7)
ntest <- n - ntrain
methods <- c("GLASSO", "GENET", "SGLASSO", "GSCAD", "GMCP")
nmethods <- length(methods)

### to hold results
repeatnum_parallel = 1
results <- data.frame(rep = rep(1:repeatnum_parallel, each = nmethods),
                      method = rep(methods, times = repeatnum_parallel),
                      time = NA,
                      iter = NA,
                      mse = NA,
                      mae = NA)
results$method <- factor(results$method, levels = methods)

maxit=1e8
eps = 1e-5
alpha_seq <- seq(0.1,0.9,length=5)
l_alphas <- length(alpha_seq)


# Detect number of available cores in your computer/server
# detectCores()

# Decide how many core will you use for the application
number_cores <- 5

# pleminaries for foreach
registerDoParallel(cores = number_cores)# Shows the number of Parallel Workers to be used
registerDoRNG(2025) # This ensures you get the same results with the paper.
# foreach computing
results_from_foreach <- 
  foreach(icount(nrep), .combine=rbind, .packages = c("grpnet","sglasso","grpreg")) %dopar% {
    
  ## define training id
  trainid <- sample.int(n, size = ntrain)
  
  ## define foldid
  k_fold <- 5 
  set_folds <- sample(rep(1:k_fold, length.out = ntrain))
  
  ## grpreg Lasso
  counter <-  1
  tic <- proc.time()
  cvmod  <- cv.grpreg(X[trainid,], y[trainid], group=group,
                      fold = set_folds, alpha = 1,
                      penalty="grLasso",max.iter=maxit,eps=eps,lambda.min = 0.05)
  toc <- proc.time() - tic
  results$time[counter] <- toc[3]
  results$iter[counter] <- sum(cvmod$fit$iter)
  pred <- predict(cvmod, X = X[-trainid,])
  results$mse[counter] <- mean((y[-trainid] - pred)^2)
  results$mae[counter] <- mean(abs(y[-trainid] - pred))
  
  
  ## grpnet
  counter <- 2
  tic <- proc.time()
  fits.cvmod <- lapply(1:l_alphas, function(ind){
    cv.grpnet(X[trainid,], y[trainid], group = group, foldid = set_folds,type.measure = "mse",
              alpha = alpha_seq[ind], thresh=eps,maxit=maxit,lambda.min.ratio = 0.05)
  })
  toc <- proc.time() - tic
  results$time[counter] <- toc[3]/l_alphas
  cvs.mod <- sapply(fits.cvmod, function(x) {min(x$cvm)})
  cvmod   <- fits.cvmod[[which.min(cvs.mod)]]
  results$iter[counter] <- sum(cvmod$grpnet.fit$npasses)
  pred <- predict(cvmod, newx = X[-trainid,])
  results$mse[counter] <- mean((y[-trainid] - pred)^2)
  results$mae[counter] <- mean(abs(y[-trainid] - pred))
  
  # sglasso
  counter <- 3
  tic <- proc.time()
  fits.cvmod <- lapply(1:l_alphas, function(ind){
    cv.sglasso(X[trainid,], y[trainid], group = group, fold = set_folds, nd=5,
               alpha = alpha_seq[ind], eps=eps,max_iter=maxit)
  })
  toc <- proc.time() - tic
  results$time[counter] <- toc[3]/(l_alphas*5)
  cvs.mod <- sapply(fits.cvmod, function(x) {min(x$cve)})
  cvmod   <- fits.cvmod[[which.min(cvs.mod)]]
  results$iter[counter] <- sum(cvmod$fit$iter)/5
  pred <- predict(cvmod, newx=X[-trainid,], s="opt")
  results$mse[counter] <- mean((y[-trainid] - pred)^2)
  results$mae[counter] <- mean(abs(y[-trainid] - pred))
  
  
  ## grpreg SCAD
  counter <- 4
  tic <- proc.time()
  cvmod  <- cv.grpreg(X[trainid,], y[trainid], group=group,
                      fold = set_folds, alpha = 1,
                      penalty="grSCAD",max.iter=maxit,eps=eps,lambda.min = 0.05)
  toc <- proc.time() - tic
  results$time[counter] <- toc[3]
  results$iter[counter] <- sum(cvmod$fit$iter)
  pred <- predict(cvmod, X = X[-trainid,])
  results$mse[counter] <- mean((y[-trainid] - pred)^2)
  results$mae[counter] <- mean(abs(y[-trainid] - pred))
  
  
  ## grpreg MCP
  counter <- 5
  tic <- proc.time()
  cvmod  <- cv.grpreg(X[trainid,], y[trainid], group=group,
                      fold = set_folds, alpha = 1,
                      penalty="grMCP",max.iter=maxit,eps=eps,lambda.min = 0.05)
  toc <- proc.time() - tic
  results$time[counter] <- toc[3]
  results$iter[counter] <- sum(cvmod$fit$iter)
  pred <- predict(cvmod, X = X[-trainid,])
  results$mse[counter] <- mean((y[-trainid] - pred)^2)
  results$mae[counter] <- mean(abs(y[-trainid] - pred))
  return(results)
  }
# When it'is done, clean up the cluster
stopImplicitCluster()
return(results_from_foreach)
}
