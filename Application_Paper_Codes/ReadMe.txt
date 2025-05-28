This package has two folders:

Simulations: this contains scripts for all simulations.
Inside of this folder:

1- foreach_sim_sglasso.R supports a function (block_sim_data) generation block design matrix with correlated and a function (simulation_sglasso) for simulation function with parallel computing.  
2-  parameter_setting_sglasso_sim_HD_20_40.R supports the simulation setting for HD cases
3-  parameter_setting_sglasso_sim_LD_4_8_16.R supports the simulation setting for LD cases
4-   After running step 1 and step 2 (or 3), simply run the script "sim_run.R" to get results of either HD or LD cases.


Real Data Examples: this contains scripts for all real data examples.
Inside of this folder:

1- sglasso_real_data_code.R script is a general function for running real data examples.
2- After run step 1, simply run the script  real_design_Birthwt.R to get results of Birthwt data example
3- After run step 1, simply run the script  real_design_Bardet.R to get results of Bardet data example
4- After run step 1, simply run the script  real_design_GenAtHum.R to get results of GenAtHum data example
 

