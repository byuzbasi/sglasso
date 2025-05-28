##########################################################
library("gglasso")
data(bardet, package = "gglasso")
group <- rep(1:20,each=5)
y <- bardet$y
X <- bardet$x
# n and p
n <- nrow(X)
p <- ncol(X) 
#length(group) == p # check
#################


replication_number = c(5) # In the paper, replication_number = c(100)
Sim_Results <- tibble()
  tryCatch({
    # Get Each Results
    real_Results = real_sglasso(X,y,group,n,replication_number)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})


# To see the results
# real_Results
### To save 
#readr::write_rds(real_Results, "Bardet_5folds_100rep.rds", compress = "gz")
