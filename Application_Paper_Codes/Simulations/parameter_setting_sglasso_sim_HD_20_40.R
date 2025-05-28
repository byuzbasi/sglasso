library(dplyr)
#
alpha_seq = seq(0.1,0.9,length=9)
d_seq= seq(0,1,length=11)

# Parameter setting in our simulation study
dimension_param_1 <- dplyr::tibble(
  expand.grid(
    n.train = c(100),
    n.validate = c(100),
    n.nest = c(100),
    p = c(1000),
    J = c(200),
    dimensionality = c("High")
  )
)



# 
dimension_param_s1 <- dimension_param_1 %>%
  dplyr::mutate(Strong_J = floor(J*c(.1))) %>%
  dplyr::mutate(Pj = p / J)
dimension_param_s2 <- dimension_param_1 %>%
  dplyr::mutate(Strong_J = floor(J*c(.2))) %>%
  dplyr::mutate(Pj = p / J)


dimension_param <- dplyr::bind_rows(dimension_param_s1,dimension_param_s2)
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
    Rho_Between = seq(0,.9,.1),
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

