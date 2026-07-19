project_dir <- normalizePath(Sys.getenv("SGLASSO_PROJECT_DIR", unset = getwd()))
setwd(project_dir)

source("paper_codes/04_scalability/sglasso_sim_function_full_tuning.R")
install_and_load_packages(install_missing = FALSE)

res_scalability <- run_scalability_experiment(
  repeatnum = 50,
  checkpoint_dir = "results/scalability_low_to_high/checkpoints",
  cores = 1
)

summary_scalability <- summarise_scalability(res_scalability)

print(summary_scalability, n = Inf, width = Inf)

dir.create("results/scalability_low_to_high", recursive = TRUE, showWarnings = FALSE)

readr::write_csv(
  res_scalability,
  "results/scalability_low_to_high/scalability_raw.csv"
)

readr::write_csv(
  summary_scalability,
  "results/scalability_low_to_high/scalability_summary.csv"
)

saveRDS(
  res_scalability,
  "results/scalability_low_to_high/scalability_raw.rds"
)

saveRDS(
  summary_scalability,
  "results/scalability_low_to_high/scalability_summary.rds"
)
