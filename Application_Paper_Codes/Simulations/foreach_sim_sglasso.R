# Select a CRAN Mirror
chooseCRANmirror(ind=11)

# The packages that is needed
packages <- c("simstudy","data.table","grpnet","sglasso","grpreg","gglasso",
              "caret","dplyr","MASS","glmnet","Metrics","matrixcalc","readr",
              "foreach","doParallel","doFuture","doRNG","mltools")


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
#######################################################
### functions used in the simulation ###
#######################################################
simsham_corrmat_toeplitz <- function (p, rho) 
{
  if (rho == 0) {
    corrmat_identity(p)
  }
  else {
    rho^abs(outer(1:p, 1:p, "-"))
  }
}
###########
block_sim_data <- function (n.train=100,
                            n.val = 100,
                            n.test = 400,
                            pj = 3,
                            J  = 10,
                            strong_J = 3,
                            rho_w = 0.7,
                            rho_b = 0.2,
                            eff.nonzero = 1, 
                            corrmat_type = c("Exchangeable", "Decay") ,
                            snr = .5,
                            nonzero.id)
{
  p = pj * J
  X_train <- matrix(rnorm(n.train * p), n.train, p)
  X_val   <- matrix(rnorm(n.val * p), n.val, p)
  X_test  <- matrix(rnorm(n.test * p), n.test, p)
  #
  if(corrmat_type == c("Decay")){
    corrmat <- simstudy::blockDecayMat(ninds = pj , nperiods = J, rho_w = rho_w, r = rho_b, pattern = "xsection")
  }else if (corrmat_type == c("Exchangeable")){
    corrmat <- simstudy::blockExchangeMat(ninds = pj , nperiods = J, rho_w = rho_w, rho_b = rho_b, pattern = "xsection")  
  } else{
    corrmat <- simsham_corrmat_toeplitz(p, rho_b)
  }
  obj <- svd(corrmat)
  corrmat_half <- obj$u %*% (sqrt(diag(obj$d))) %*% t(obj$v)
  X_train <- X_train %*% corrmat_half
  X_val   <- X_val %*% corrmat_half
  X_test  <- X_test %*% corrmat_half
  # beta
  #beta.nonzero <- runif(strong_J, -eff.nonzero, eff.nonzero)
  beta.nonzero <- rep(eff.nonzero,strong_J)
  beta <- rep(0, p)
  #nonzero.id <- sample(J, strong_J) # random true group indices
  for (k in 1:strong_J) {
    beta[((nonzero.id[k]-1)*pj+1):(nonzero.id[k]*pj)] <- rep(beta.nonzero[k], pj)
  }
  #
  risk_null <- as.numeric(t(beta) %*% corrmat %*% beta)
  sigma <- sqrt(risk_null/snr)
  error_null <- risk_null + sigma^2
  #
  y_train <- as.numeric(X_train %*% beta + rnorm(n.train) * sigma)
  y_val   <- as.numeric(X_val %*% beta + rnorm(n.val) * sigma)
  y_test  <- as.numeric(X_test %*% beta + rnorm(n.test) * sigma)
  group <- rep(1:J,each=pj)
  true_groups_indx  <- sort(nonzero.id)
  beta_true_groups <- rep(0,J); 
  beta_true_groups[true_groups_indx] = 1; 
  beta_true_groups <- as.factor(beta_true_groups)
  # Return results
  list(X_train = X_train, 
       X_val = X_val, 
       X_test = X_test, 
       y_train = y_train, 
       y_val = y_val,
       y_test = y_test, 
       beta=beta,
       group=group,
       true_groups_indx=true_groups_indx,
       corrmat=corrmat,sigma=sigma,
       risk_null=risk_null, error_null=error_null,
       beta_true_groups=beta_true_groups)
}

###############################################
MSE.funct_y <- function(beta.est, newx, newy, mean_y.train, mean_X.train, beta.true)
{
    fit = newx %*% beta.true - mean_y.train - cbind(1,scale(newx, mean_X.train,F)) %*% beta.est
  return(drop( t (fit) %*% fit / length (newy)))
}


#######################################################
#######################################################
#######################################################
#### simulation function
#######################################################
# ntrain the number of train observations
# nvalidate the number of validation set observations
# ntest the number of test set observations
# n  the total number of observations
# pj the number of elements of each group
# J the number of groups
# strong_J the number of strong (non-zero) groups
# rho_w correlation of within groups
# rho_b correlation of between groups
# eff.nonzero strength of non-zero coefficients
# corrmat_type the type correlation matrix of groups from c("Exchangeable", "Decay")
# tun_par the optimal selection criteria "CP" or Min value from validation set
# snr signal to noise ratio
# repeatnum the number of simulations
# nfolds the fold of cross-validations
# dimensionality is only for paper simulation, does not affect simulation
# alpha_seq elastic net parameter set
# d_seq scale parameter set
# eps convergence parameter
# maxit maximum number of iterations
# set_seed the reproducibility
simulation_sglasso <- function(ntrain,nvalidate,ntest,pj,J,strong_J,rho_w,rho_b,eff.nonzero,corrmat_type,tun_par,snr,
                              repeatnum,set_seed,nlambda=100,
                              alpha_seq = seq(0.1,0.9,length=9),
                              d_seq= seq(0,1,length=11),
                              dimensionality=c("low"),
                              standardize = FALSE,
                              eps = 1e-5,
                              maxit=1e8){
  
  set.seed(set_seed)
  nonzero.id <- sample(J, strong_J) # random true group indices
  n <- ntrain+nvalidate+ntest
  # methods
  methods <- c("oracle","grpreg_Lasso", "grpnet", "sglasso", "grpreg_SCAD", "grpreg_MCP")
  nmethods <- length(methods)
  repeatnum_parallel = 1
  # initialize performance data frame
  results <- data.frame(rep = rep(1:repeatnum_parallel, each = nmethods),
                        dimensionality=dimensionality,
                        method = rep(methods, times = repeatnum_parallel),
                        corrmat_type=corrmat_type,
                        tun_par=tun_par,
                        standardize=standardize,
                        time=NA,
                        iter=NA,
                        n = n,
                        n.train=ntrain,
                        n.val = nvalidate,
                        n.test = ntest,
                        p = pj * J,
                        pj = pj,
                        J = J,
                        strong_J=strong_J,
                        rho_w=rho_w,
                        rho_b = rho_b,
                        eff.nonzero=eff.nonzero,
                        snr=snr,
                        prop = NA,
                        nzs = NA,
                        precision=NA,
                        MSE_y=NA)
  results$method <- factor(results$method, levels = methods)
  # pleminaries for foreach
  registerDoParallel(cores=2) # Shows the number of Parallel Workers to be used, you can change accordingly
  registerDoRNG(set_seed)
  # foreach computing
  results_from_foreach <- 
    foreach(icount(repeatnum), .combine=rbind) %dopar% {
    
    # generate data via block_sim_data
    data_r <- block_sim_data(ntrain,nvalidate,ntest, pj, J, strong_J,rho_w,rho_b,
                             eff.nonzero, corrmat_type, snr,nonzero.id)
    
    # Assign results of block_sim_data
    beta.true <- data_r$beta                       # true beta
    group <- data_r$group                          # group index
    true_groups_indx = data_r$true_groups_indx     # non zero groups
    beta_true_groups <- data_r$beta_true_groups    # non zero group index of beta
    # Covariance Matrix, variance, error and risk of null
    Big_sigma = data_r$corrmat                     # Covariance matrix
    Small_sigma2 = data_r$sigma^2                  # varince of sigma
    error_null = data_r$error_null                 
    risk_null = data_r$risk_null
    
    
    # Assign train/val and test data of X
    x.train     <- data_r$X_train
    x.validate  <- data_r$X_val
    x.test      <- data_r$X_test
    
    COV_X_test = cov(x.test)
    
    # Assign train/val and test data of y
    y.train     <- data_r$y_train
    y.validate  <- data_r$y_val
    y.test      <- data_r$y_test
    
    
    
    if (standardize == TRUE) {
      mean_X.train <- apply (x.train, 2, mean)
      mean_y.train <- mean (y.train)
      x.train <- scale (x.train, mean_X.train, FALSE)
      y.train <- y.train - mean_y.train 
      }else {mean_X.train <- rep (0, length = ncol (x.train)) 
      mean_y.train <- 0}

    
    ### ORACLE
    counter <- 1
    results$time[counter] = 0
    results$iter[counter] = 0
    betahat  = beta.true
    delta    = 0
    risk = 0
    err.test <- risk + Small_sigma2
    results$prop[counter]     = (1 - err.test / error_null)
    results$nzs[counter] = strong_J
    beta.est_gruops <- beta_true_groups
    cm_est <- confusionMatrix(factor(beta.est_gruops,levels=c('0','1')),beta_true_groups,mode = "everything")
    results$precision[counter] = cm_est$byClass[5]
    results$MSE_y[counter] = MSE.funct_y(c(0,betahat), x.test, y.test, mean_y.train, mean_X.train, beta.true)
    
    
    
    ## grpreg
    counter <- 2
    time_start <- proc.time()
    model_grpreg  <- grpreg(x.train, y.train, group=group, penalty="grLasso",max.iter=maxit,eps=eps,lambda.min = 0.05)
    time_end <- proc.time() - time_start
    results$time[counter] = time_end[3]
    results$iter[counter] = sum(model_grpreg$iter)
    error_train <- apply((y.train-predict(model_grpreg,x.train))^2,2,sum)
    error_val   <- apply((y.validate-predict(model_grpreg,x.validate))^2,2,sum)
    betahat_c= as.matrix(coef(model_grpreg))
    betahat  = betahat_c[-1,]
    delta    = sweep(betahat,1,beta.true,"-")
    risk = diag(t(delta) %*% Big_sigma %*% delta)
    err.test <- risk + Small_sigma2
    val_risk   <- apply((sweep(cbind(1,scale(x.validate, mean_X.train,F)) %*% betahat_c,1,x.validate %*% beta.true - mean_y.train,"-"))^2,2,mean) 
    ############################################################
    # Select optimal lambda       ##############################
    model_df = model_grpreg$df
    if(tun_par == "CP"){
      rss.big = rev(error_val)[1]
      df.big  = NROW(x.train) - rev(model_df)[1]
      sigma2  = rss.big/df.big
      Cp_val =  error_val/sigma2 - NROW(x.train) + 2 * model_df
      tun.val  = which.min(Cp_val) # minimum lambda index
    }else if (tun_par == "Min_val") {
      tun.val  = which.min(error_val) # minimum lambda index  
    } else tun.val = which.min(val_risk)
    results$tun_par[counter] = tun_par
    ############################################################
    results$prop[counter]     = (1 - err.test / error_null)[tun.val]
    nzs      = sapply(apply(betahat!=0, 2, function(x) unique(group[x])), length)
    if(length(nzs) > nlambda) {nzs = rep(J,nlambda) }
    results$nzs[counter] = nzs[tun.val]
    selected_groups <-  c(apply(as.matrix(betahat[,tun.val]) !=0, 2, function(x) unique(group[x])))
    beta.est_gruops <- rep(0,J); beta.est_gruops[selected_groups] = 1
    cm_est <- confusionMatrix(factor(beta.est_gruops,levels=c('0','1')),beta_true_groups,mode = "everything")
    results$precision[counter] = cm_est$byClass[5]
    results$MSE_y[counter] = MSE.funct_y(betahat_c[,tun.val], x.test, y.test, mean_y.train, mean_X.train, beta.true)
    
    
    
    
    
    # grpnet
    counter <- 3
    l_alphas <- length(alpha_seq)
    time_start <- proc.time()
    fits.grpnet <- lapply(1:l_alphas, function(ind){
      grpnet(x.train, y.train, group, nlambda=nlambda,type.measure = "mse",
             alpha = alpha_seq[ind], thresh=eps,maxit=maxit,lambda.min.ratio = 0.05)
    })
    time_end <- proc.time() - time_start
    results$time[counter] = time_end[3]/l_alphas
    results$iter[counter] = floor(mean(sapply(1:l_alphas, function(x) sum(fits.grpnet[[x]]$npasses))))/l_alphas
    error_val_mat <- sapply(1:l_alphas, function(ind){
      apply((sweep(cbind(1,scale(x.validate, mean_X.train,F)) %*% coef(fits.grpnet[[ind]]),1,x.validate %*% beta.true - mean_y.train,"-"))^2,2,mean)
    })
    min_ind_grpnet  <- which(error_val_mat == min(error_val_mat), arr.ind = TRUE)
    fits.grpnet_opt <- fits.grpnet[[min_ind_grpnet[2]]]
    error_train <- apply((y.train-predict(fits.grpnet_opt,x.train))^2,2,sum)
    error_val   <- apply((y.validate-predict(fits.grpnet_opt,x.validate))^2,2,sum)
    betahat_c= as.matrix(coef(fits.grpnet_opt))
    betahat  = betahat_c[-1,]
    delta    = sweep(betahat,1,beta.true,"-")
    risk = diag(t(delta) %*% Big_sigma %*% delta)
    err.test <- risk + Small_sigma2
    val_risk   <- apply((sweep(cbind(1,scale(x.validate, mean_X.train,F)) %*% betahat_c,1,x.validate %*% beta.true - mean_y.train,"-"))^2,2,mean) 
    ############################################################
    # Select optimal lambda       ##############################
    model_df = model_grpreg$df
    if(tun_par == "CP"){
      rss.big = rev(error_val)[1]
      df.big  = NROW(x.train) - rev(model_df)[1]
      sigma2  = rss.big/df.big
      Cp_val =  error_val/sigma2 - NROW(x.train) + 2 * model_df
      tun.val  = which.min(Cp_val) # minimum lambda index
    }else if (tun_par == "Min_val") {
      tun.val  = which.min(error_val) # minimum lambda index  
    } else tun.val = which.min(val_risk)
    results$tun_par[counter] = tun_par
    ############################################################
    results$prop[counter]     = (1 - err.test / error_null)[tun.val]
    nzs      = sapply(apply(betahat!=0, 2, function(x) unique(group[x])), length)
    if(length(nzs) > nlambda) {nzs = rep(J,nlambda) }
    results$nzs[counter] = nzs[tun.val]
    selected_groups <-  c(apply(as.matrix(betahat[,tun.val]) !=0, 2, function(x) unique(group[x])))
    beta.est_gruops <- rep(0,J); beta.est_gruops[selected_groups] = 1
    cm_est <- confusionMatrix(factor(beta.est_gruops,levels=c('0','1')),beta_true_groups,mode = "everything")
    results$precision[counter] = cm_est$byClass[5]
    results$MSE_y[counter] = MSE.funct_y(betahat_c[,tun.val], x.test, y.test, mean_y.train, mean_X.train, beta.true)
    
    
    ## sglasso
    counter <- 4
    l_ds <- 11
    error_val_array <- array(NA,dim=c(nlambda,l_alphas,l_ds))
    time_start <- proc.time()
    fits.sglasso <- lapply(1:l_alphas, function(ind){
      sglasso(x.train, y.train, group, nlambda=nlambda,alpha = alpha_seq[ind], nd=l_ds,
              eps=eps,max_iter=maxit, lambda.min = 0.05)
    })
    time_end <- proc.time() - time_start
    results$time[counter] = time_end[3]/(l_alphas*l_ds)
    results$iter[counter] = floor(sum(sapply(1:l_alphas, function(x) sum(fits.sglasso[[x]]$iter)))/(l_alphas*l_ds))
    for (ind in 1:l_alphas){
      for(di in 1:l_ds){
        #error_val_array[,ind,di]   <- apply((y.validate-cbind(1,x.validate) %*% fits.sglasso[[ind]]$betas[,,di])^2,2,sum)
        error_val_array[,ind,di]   <- apply((sweep(cbind(1,scale(x.validate, mean_X.train,F)) %*% fits.sglasso[[ind]]$betas[,,di],1,x.validate %*% beta.true - mean_y.train,"-"))^2,2,mean) 
      }
    }
    min_ind_sglasso <- which(error_val_array == min(error_val_array), arr.ind = TRUE)
    fit.sglasso   <- fits.sglasso[[min_ind_sglasso[2]]]
    error_train <- apply((y.train-cbind(1,x.train) %*% fit.sglasso$betas[,,min_ind_sglasso[3]])^2,2,sum)
    error_val   <- apply((y.validate-cbind(1,x.validate) %*% fit.sglasso$betas[,,min_ind_sglasso[3]])^2,2,sum)
    betahat_c= fit.sglasso$betas[,,min_ind_sglasso[3]]
    betahat  = betahat_c[-1,]
    delta    = sweep(betahat,1,beta.true,"-")
    risk = diag(t(delta) %*% Big_sigma %*% delta)
    err.test <- risk + Small_sigma2
    val_risk   <- apply((sweep(cbind(1,scale(x.validate, mean_X.train,F)) %*% betahat_c,1,x.validate %*% beta.true - mean_y.train,"-"))^2,2,mean) 
    ############################################################
    # Select optimal lambda       ##############################
    model_df = model_grpreg$df
    if(tun_par == "CP"){
      rss.big = rev(error_val)[1]
      df.big  = NROW(x.train) - rev(model_df)[1]
      sigma2  = rss.big/df.big
      Cp_val =  error_val/sigma2 - NROW(x.train) + 2 * model_df
      tun.val  = which.min(Cp_val) # minimum lambda index
    }else if (tun_par == "Min_val") {
      tun.val  = which.min(error_val) # minimum lambda index  
    } else tun.val = which.min(val_risk)
    results$tun_par[counter] = tun_par
    ############################################################
    results$prop[counter]     = (1 - err.test / error_null)[tun.val]
    nzs      = sapply(apply(betahat!=0, 2, function(x) unique(group[x])), length)
    if(length(nzs) > nlambda) {nzs = rep(J,nlambda) }
    results$nzs[counter] = nzs[tun.val]
    selected_groups <-  c(apply(as.matrix(betahat[,tun.val]) !=0, 2, function(x) unique(group[x])))
    beta.est_gruops <- rep(0,J); beta.est_gruops[selected_groups] = 1
    cm_est <- confusionMatrix(factor(beta.est_gruops,levels=c('0','1')),beta_true_groups,mode = "everything")
    results$precision[counter] = cm_est$byClass[5]
    results$MSE_y[counter] = MSE.funct_y(betahat_c[,tun.val], x.test, y.test, mean_y.train, mean_X.train, beta.true)
    
    
    ## grpreg SCAD
    counter <- 5
    time_start <- proc.time()
    model_grpreg  <- grpreg(x.train, y.train, group=group, penalty="grSCAD",max.iter=maxit,eps=eps,lambda.min = 0.05)
    time_end <- proc.time() - time_start
    results$time[counter] = time_end[3]
    results$iter[counter] = sum(model_grpreg$iter)
    error_train <- apply((y.train-predict(model_grpreg,x.train))^2,2,sum)
    error_val   <- apply((y.validate-predict(model_grpreg,x.validate))^2,2,sum)
    betahat_c= as.matrix(coef(model_grpreg))
    betahat  = betahat_c[-1,]
    delta    = sweep(betahat,1,beta.true,"-")
    risk = diag(t(delta) %*% Big_sigma %*% delta)
    err.test <- risk + Small_sigma2
    val_risk   <- apply((sweep(cbind(1,scale(x.validate, mean_X.train,F)) %*% betahat_c,1,x.validate %*% beta.true - mean_y.train,"-"))^2,2,mean) 
    ############################################################
    # Select optimal lambda       ##############################
    model_df = model_grpreg$df
    if(tun_par == "CP"){
      rss.big = rev(error_val)[1]
      df.big  = NROW(x.train) - rev(model_df)[1]
      sigma2  = rss.big/df.big
      Cp_val =  error_val/sigma2 - NROW(x.train) + 2 * model_df
      tun.val  = which.min(Cp_val) # minimum lambda index
    }else if (tun_par == "Min_val") {
      tun.val  = which.min(error_val) # minimum lambda index  
    } else tun.val = which.min(val_risk)
    results$tun_par[counter] = tun_par
    ############################################################
    results$prop[counter]     = (1 - err.test / error_null)[tun.val]
    nzs      = sapply(apply(betahat!=0, 2, function(x) unique(group[x])), length)
    if(length(nzs) > nlambda) {nzs = rep(J,nlambda) }
    results$nzs[counter] = nzs[tun.val]
    selected_groups <-  c(apply(as.matrix(betahat[,tun.val]) !=0, 2, function(x) unique(group[x])))
    beta.est_gruops <- rep(0,J); beta.est_gruops[selected_groups] = 1
    cm_est <- confusionMatrix(factor(beta.est_gruops,levels=c('0','1')),beta_true_groups,mode = "everything")
    results$precision[counter] = cm_est$byClass[5]
    results$MSE_y[counter] = MSE.funct_y(betahat_c[,tun.val], x.test, y.test, mean_y.train, mean_X.train, beta.true)
    
    
    
    ## grpreg SCAD
    counter <- 6
    time_start <- proc.time()
    model_grpreg  <- grpreg(x.train, y.train, group=group, penalty="grMCP",max.iter=maxit,eps=eps,lambda.min = 0.05)
    time_end <- proc.time() - time_start
    results$time[counter] = time_end[3]
    results$iter[counter] = sum(model_grpreg$iter)
    error_train <- apply((y.train-predict(model_grpreg,x.train))^2,2,sum)
    error_val   <- apply((y.validate-predict(model_grpreg,x.validate))^2,2,sum)
    betahat_c= as.matrix(coef(model_grpreg))
    betahat  = betahat_c[-1,]
    delta    = sweep(betahat,1,beta.true,"-")
    risk = diag(t(delta) %*% Big_sigma %*% delta)
    err.test <- risk + Small_sigma2
    val_risk   <- apply((sweep(cbind(1,scale(x.validate, mean_X.train,F)) %*% betahat_c,1,x.validate %*% beta.true - mean_y.train,"-"))^2,2,mean) 
    ############################################################
    # Select optimal lambda       ##############################
    model_df = model_grpreg$df
    if(tun_par == "CP"){
      rss.big = rev(error_val)[1]
      df.big  = NROW(x.train) - rev(model_df)[1]
      sigma2  = rss.big/df.big
      Cp_val =  error_val/sigma2 - NROW(x.train) + 2 * model_df
      tun.val  = which.min(Cp_val) # minimum lambda index
    }else if (tun_par == "Min_val") {
      tun.val  = which.min(error_val) # minimum lambda index  
    } else tun.val = which.min(val_risk)
    results$tun_par[counter] = tun_par
    ############################################################
    results$prop[counter]     = (1 - err.test / error_null)[tun.val]
    nzs      = sapply(apply(betahat!=0, 2, function(x) unique(group[x])), length)
    if(length(nzs) > nlambda) {nzs = rep(J,nlambda) }
    results$nzs[counter] = nzs[tun.val]
    selected_groups <-  c(apply(as.matrix(betahat[,tun.val]) !=0, 2, function(x) unique(group[x])))
    beta.est_gruops <- rep(0,J); beta.est_gruops[selected_groups] = 1
    cm_est <- confusionMatrix(factor(beta.est_gruops,levels=c('0','1')),beta_true_groups,mode = "everything")
    results$precision[counter] = cm_est$byClass[5]
    results$MSE_y[counter] = MSE.funct_y(betahat_c[,tun.val], x.test, y.test, mean_y.train, mean_X.train, beta.true)
    #### 
    return(results)
    }
  # When it'is done, clean up the cluster
  stopImplicitCluster()
  return(results_from_foreach)
}