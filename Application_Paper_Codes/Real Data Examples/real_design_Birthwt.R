##########################################################
library("grpreg")
data(Birthwt,package = "grpreg")
X <- Birthwt$X
y <- Birthwt$bwt
K <- as.integer(table(Birthwt$group))
group <- rep(1:length(K),K)
# n and p
n <- nrow(X)
p <- ncol(X)
#length(group) == p # check
#################


replication_number = c(5) # In the paper, replication_number = c(100)
Sim_Results <- tibble()
  tryCatch({
    # Get Each Results
    Sim_Results = real_sglasso(X,y,group,n,replication_number)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

# To see the results
# Sim_Results
#readr::write_rds(Sim_Results, "Birthwt_5folds_100rep.rds", compress = "gz")
