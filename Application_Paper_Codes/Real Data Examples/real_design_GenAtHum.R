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


# To see all results
Sim_Results

# To see results in the paper
options(pillar.sigfig = 4)
Sim_Results %>% 
  group_by(method) %>% 
  mutate(mse_se = sd(mse)/sqrt(replication_number)) %>%
  dplyr::summarise(across(everything(), median),.groups = 'drop',na.rm = TRUE) %>%
  mutate(PE = paste0( round(mse,3), "(", round(mse_se,3), ")")) %>%
  mutate(RPE = mse/mse[3]) %>%
  dplyr::select(c("method","time","iter","PE","RPE"))%>% 
  dplyr::slice(c(3, 1, 2, 4, 5))