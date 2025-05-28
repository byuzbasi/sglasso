##########################################################
data(GenAtHum,package = "sglasso")
X <- GenAtHum$X
y <- GenAtHum$y
group <- GenAtHum$group
n =  nrow(X)
p =  ncol(X)
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
#readr::write_rds(real_Results, "gene_atlas_plots_real_5Folds.rds", compress = "gz")