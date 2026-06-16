############################################################
# Fixed-tuning scalability analysis for reviewer response.
# This compares one regularization path per method with externally scaled data.
############################################################

script_file <- tryCatch(
  normalizePath(sys.frame(1)$ofile),
  error = function(e) NA_character_
)

if (is.na(script_file)) {
  cmd_args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args_full, value = TRUE)
  script_file <- if (length(file_arg) > 0L) {
    normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE)
  } else {
    normalizePath("04_scalability/run_scalability_analysis.R",
                  mustWork = FALSE)
  }
}

script_dir <- dirname(script_file)
source_file <- file.path(script_dir, "sglasso_sim_function_full_tuning.R")
if (!file.exists(source_file)) {
  stop("Cannot find sglasso_sim_function_full_tuning.R near ", script_file)
}
source(source_file)

install_and_load_packages(
  install_missing = FALSE
)

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}

parse_num_grid <- function(x) {
  as.numeric(strsplit(x, ",", fixed = TRUE)[[1]])
}

log_msg <- function(...) {
  cat(
    sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf(...),
    "\n"
  )
}

log_msg("Starting scalability experiment")

OUTDIR <- get_arg("outdir", "results/scalability")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

p_grid <- parse_num_grid(
  get_arg("p_grid", "1000,5000,10000,15000")
)
J_grid <- parse_num_grid(
  get_arg("J_grid", "200,500,1000,1500")
)
repeatnum <- as.integer(get_arg("nrep", "5"))
seed <- as.integer(get_arg("seed", "2026"))
alpha_fixed <- as.numeric(get_arg("alpha", "0.5"))
d_fixed <- as.numeric(get_arg("d", "0"))
nlambda <- as.integer(get_arg("nlambda", "50"))

log_msg(
  "Settings: p_grid=%s; J_grid=%s; nrep=%d; alpha=%.3f; d=%.3f; nlambda=%d",
  paste(p_grid, collapse = ","),
  paste(J_grid, collapse = ","),
  repeatnum,
  alpha_fixed,
  d_fixed,
  nlambda
)

res_scalability <- run_scalability_experiment(
  p_grid = p_grid,
  J_grid = J_grid,
  pj = 5,
  
  n_train = 100,
  n_val = 0,
  n_test = 0,
  
  strong_frac = 0.10,
  
  rho_w = 0.9,
  rho_b = 0.5,
  snr = 1.528,
  
  repeatnum = repeatnum,
  seed = seed,
  
  alpha_fixed = alpha_fixed,
  d_fixed = d_fixed,
  
  nlambda = nlambda,
  standardize = FALSE,
  
  eps = 1e-5,
  maxit = 1e5,
  
  signal_pattern = "weak_strong_mixed"
)

summary_scalability <- summarise_scalability(
  res_scalability
)

readr::write_csv(
  res_scalability,
  file.path(OUTDIR, "scalability_raw.csv")
)

readr::write_csv(
  summary_scalability,
  file.path(OUTDIR, "scalability_summary.csv")
)

saveRDS(
  res_scalability,
  file.path(OUTDIR, "scalability_raw.rds")
)

saveRDS(
  summary_scalability,
  file.path(OUTDIR, "scalability_summary.rds")
)

print(summary_scalability)

log_msg("Scalability experiment completed")
