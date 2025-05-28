library(dplyr)
#
alpha_seq = seq(0.1,.9,length=9)
d_seq= seq(0,1,length=11)
# Parameter setting in our simulation study
dimension_param <- dplyr::tibble(
  expand.grid(
  n.train = c(500),
  n.validate = c(200),
  n.nest = c(200),
  p = c(100),
  J = c(20),
  Strong_J = c(4,8,16),
  dimensionality = c("Low")
  )
)

# Calculate Pj, the number of group elements
dimension_param <- dimension_param %>%
  dplyr::mutate(Pj = p / J)


# strength of beta coefficients
Eff.nonzero = c(1)
# signal to noise
SNR = exp(seq(log(0.05),log(6),length=8))
# tune parameter
tun_par = c("Risk")
std = c("TRUE")


# Correlation type
cor_param <- dplyr::as_tibble(
  expand.grid(
    Corrmat_type = c("Exchangeable"),
    Rho_Within = c(0.9),
    Rho_Between = seq(0,0.9,length=10),
    SNR = SNR,
    Eff.nonzero=Eff.nonzero,
    tun_par=tun_par,
    std=std
  )
)
# combine the all cor_param
cor_param <- cor_param %>% dplyr::mutate(across(where(is.factor), as.character))
# combine all
sim_param <- merge(dimension_param, cor_param)
# use sapply for simulation loop
parameters = sapply(sim_param,list)






length(parameters$n.train)