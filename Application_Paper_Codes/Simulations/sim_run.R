replication_number = c(50)  # The number of the replications for HD cases
#replication_number = c(500) # The number of the replications for LD cases
set_seed = 2025             # For reproducible results
Sim_Results <- tibble()
# START the Simulations
for (i in 1:NROW(sim_param)) {
  tryCatch({
    # Get Each Results
    get_resutls = simulation_gllreg(ntrain = parameters$n.train[i],
                                    nvalidate = parameters$n.validate[i],
                                    ntest= parameters$n.nest[i], 
                                    pj = parameters$Pj[i], 
                                    J  = parameters$J[i], 
                                    strong_J = parameters$Strong_J[i],
                                    rho_w = parameters$Rho_Within[i],
                                    rho_b = parameters$Rho_Between[i], 
                                    eff.nonzero = parameters$Eff.nonzero[i], 
                                    corrmat_type = parameters$Corrmat_type[i], 
                                    snr = parameters$SNR[i],
                                    standardize = parameters$std[i],
                                    tun_par = parameters$tun_par[i],
                                    dimensionality=parameters$dimensionality[i],
                                    d_seq=d_seq,
                                    alpha_seq=alpha_seq,
                                    repeatnum=replication_number,
                                    set_seed=set_seed)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  # Storage all results here
  Sim_Results = bind_rows(Sim_Results,get_resutls)
  }

### to see the results run 
Sim_Results