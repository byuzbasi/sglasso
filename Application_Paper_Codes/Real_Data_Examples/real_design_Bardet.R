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


################################################################################
# This may take some time depending on the replication number.
replication_number = c(5) # In the paper, replication_number = c(100)
Sim_Results <- tibble()
tryCatch({
  # Get Each Results
  Sim_Results = real_sglasso(X,y,group,n,replication_number)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
################################################################################


################################################################################
# To see all results
Sim_Results
################################################################################


################################################################################
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
################################################################################