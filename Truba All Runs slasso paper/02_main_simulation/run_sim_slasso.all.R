############################################################
# Run SGLASSO simulation scenarios for manuscript revision
#
# Main aim:
# - Use Min_val as the main tuning criterion because it mimics
#   the validation/CV strategy used in the real-data analysis.
# - Include both unfavorable and favorable regimes for SGLASSO:
#   low/moderate/high between-group correlation and low/high SNR.
# - Save both raw replicate-level results and summary tables.
############################################################

options(stringsAsFactors = FALSE)

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      sprintf(...), "\n", sep = "")
  flush.console()
}

############################################################
# 1. Source functions and load packages
############################################################

cmd_args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args_full, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = TRUE))
} else {
  getwd()
}

source_candidates <- c(
  file.path(script_dir, "sglasso_sim_function_full_tuning.R"),
  file.path(getwd(), "02_main_simulation", "sglasso_sim_function_full_tuning.R"),
  file.path(getwd(), "sglasso_sim_function_full_tuning.R")
)
source_file <- source_candidates[file.exists(source_candidates)][1L]
if (is.na(source_file)) {
  stop("Cannot find sglasso_sim_function_full_tuning.R near the run script or current directory.")
}

source(source_file)
install_and_load_packages(install_missing = FALSE)

############################################################
# 2. Output folders
############################################################

outdir <- Sys.getenv("OUTDIR", unset = "results/main_simulation")
rds_dir <- file.path(outdir, "rds")
tables_dir <- file.path(outdir, "tables")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

############################################################
# 3. Simulation design
############################################################

# For final server run, increase repeatnum if desired.
repeatnum <- as.integer(Sys.getenv("REPEATNUM", unset = "100"))
cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "55"))
cores <- max(1L, cores)

seed <- 2026

# Main manuscript scenarios:
# LD-4: sparse, often unfavorable for SGLASSO
# LD-8: moderate sparsity, useful transition case
# LD-16: denser active group structure, more favorable for SGLASSO
paper_settings <- get_paper_settings()

# Between-group correlations: include moderate and high regimes.
rho_b_grid <- seq(0, 0.9, by = 0.1)

# Low and high SNR settings. The high-SNR settings showed the clearest
# contrast between prediction-oriented and sparsity-oriented behavior.
snr_grid <- get_snr_grid()

# Include a simple homogeneous signal and a more realistic heterogeneous signal.
signal_pattern_grid <- c("homogeneous","mixed_signs","weak_strong_mixed")

# Main tuning criterion. Min_val is the main choice because it corresponds
# to validation/CV-type tuning. Val_Risk can be added only as a diagnostic.
tuning_criterion_grid <- c("Min_val")
#tuning_criterion_grid <- c("Min_val", "Val_Risk")

alpha_seq <- seq(0.1, 0.9, by = 0.1)
d_seq <- seq(0, 1, by = 0.1)
nlambda <- 100

readr::write_csv(paper_settings, file.path(tables_dir, "simulation_paper_settings.csv"))
readr::write_csv(
  expand.grid(
    setting = paper_settings$setting,
    rho_b = rho_b_grid,
    snr = snr_grid,
    signal_pattern = signal_pattern_grid,
    tuning_criterion = tuning_criterion_grid,
    stringsAsFactors = FALSE
  ),
  file.path(tables_dir, "simulation_scenario_grid.csv")
)

log_msg("Simulation script started")
log_msg("Function source     : %s", source_file)
log_msg("Output directory    : %s", outdir)
log_msg("repeatnum          : %d", repeatnum)
log_msg("cores              : %d", cores)
log_msg("settings           : %s", paste(paper_settings$setting, collapse = ", "))
log_msg("rho_b_grid         : %s", paste(rho_b_grid, collapse = ", "))
log_msg("snr_grid           : %s", paste(round(snr_grid, 3), collapse = ", "))
log_msg("signal patterns    : %s", paste(signal_pattern_grid, collapse = ", "))
log_msg("tuning criteria    : %s", paste(tuning_criterion_grid, collapse = ", "))

safe_num <- function(x, digits = 3) {
  gsub("\\.", "p", format(round(x, digits), nsmall = digits, trim = TRUE))
}

############################################################
# 4. Run all scenarios
############################################################

res_list <- list()
job_id <- 1L
n_jobs <- nrow(paper_settings) * length(tuning_criterion_grid)

for (i in seq_len(nrow(paper_settings))) {
  
  st <- paper_settings[i, ]
  
  for (criterion in tuning_criterion_grid) {
    
    log_msg(
      "[%d/%d] Running %s | criterion=%s | p=%d | J=%d | strong_J=%d",
      job_id,
      n_jobs,
      st$setting,
      criterion,
      st$p,
      st$J,
      st$strong_J
    )
    
    res <- simulation_sglasso(
      repeatnum = repeatnum,
      seed = seed,
      cores = cores,
      
      n_train = st$n_train,
      n_val = st$n_val,
      n_test = st$n_test,
      pj = st$pj,
      J = st$J,
      strong_J = st$strong_J,
      
      rho_w = 0.9,
      rho_b_grid = rho_b_grid,
      
      eff_nonzero = 1,
      corrmat_type = "Exchangeable",
      
      snr_grid = snr_grid,
      signal_pattern_grid = signal_pattern_grid,
      
      alpha_seq = alpha_seq,
      d_seq = d_seq,
      nlambda = nlambda,
      
      standardize = TRUE,
      tuning_criterion = criterion,
      eps = 1e-5,
      maxit = 1e5
    )
    
    res$setting <- st$setting
    res$dimensionality <- st$dimensionality
    res$tuning_criterion <- criterion
    
    if ("sparsity_level" %in% names(st)) {
      res$sparsity_level <- st$sparsity_level
    }
    
    fname_base <- sprintf(
      "%s__%s__all_scenarios",
      st$setting,
      criterion
    )
    
    saveRDS(
      res,
      file = file.path(rds_dir, paste0(fname_base, ".rds"))
    )
    
    readr::write_csv(
      res,
      file.path(tables_dir, paste0(fname_base, ".csv"))
    )
    
    res_list[[length(res_list) + 1L]] <- res
    job_id <- job_id + 1L
  }
}

all_res <- dplyr::bind_rows(res_list)

############################################################
# 5. Save combined raw results
############################################################

saveRDS(all_res, file.path(rds_dir, "all_sim_slasso_results.rds"))
readr::write_csv(all_res, file.path(tables_dir, "all_sim_slasso_results.csv"))

############################################################
# 6. Summary tables
############################################################

paper_summary_snr_rhob <- all_res %>%
  dplyr::filter(converged == TRUE) %>%
  dplyr::group_by(setting, dimensionality, signal_pattern,
                  tuning_criterion, rho_b, snr, method) %>%
  dplyr::summarise(
    risk = mean(risk, na.rm = TRUE),
    pe_mse = mean(pe_mse, na.rm = TRUE),
    MSE_y = mean(MSE_y, na.rm = TRUE),
    MSE_beta = mean(MSE_beta, na.rm = TRUE),
    pve = mean(pve, na.rm = TRUE),
    FDR = mean(FDR, na.rm = TRUE),
    selected_groups = mean(selected_groups, na.rm = TRUE),
    selected_features = mean(selected_features, na.rm = TRUE),
    time = mean(time, na.rm = TRUE),
    n_ok = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    snr_label = round(snr, 3),
    rho_b_label = round(rho_b, 1)
  )

paper_summary_setting <- all_res %>%
  dplyr::filter(converged == TRUE) %>%
  dplyr::group_by(setting, dimensionality, signal_pattern,
                  tuning_criterion, method) %>%
  dplyr::summarise(
    risk = mean(risk, na.rm = TRUE),
    pe_mse = mean(pe_mse, na.rm = TRUE),
    MSE_y = mean(MSE_y, na.rm = TRUE),
    MSE_beta = mean(MSE_beta, na.rm = TRUE),
    pve = mean(pve, na.rm = TRUE),
    FDR = mean(FDR, na.rm = TRUE),
    selected_groups = mean(selected_groups, na.rm = TRUE),
    selected_features = mean(selected_features, na.rm = TRUE),
    time = mean(time, na.rm = TRUE),
    n_ok = dplyr::n(),
    .groups = "drop"
  )

paper_summary_overall <- all_res %>%
  dplyr::filter(converged == TRUE) %>%
  dplyr::group_by(tuning_criterion, method) %>%
  dplyr::summarise(
    risk = mean(risk, na.rm = TRUE),
    pe_mse = mean(pe_mse, na.rm = TRUE),
    MSE_y = mean(MSE_y, na.rm = TRUE),
    MSE_beta = mean(MSE_beta, na.rm = TRUE),
    pve = mean(pve, na.rm = TRUE),
    FDR = mean(FDR, na.rm = TRUE),
    selected_groups = mean(selected_groups, na.rm = TRUE),
    selected_features = mean(selected_features, na.rm = TRUE),
    time = mean(time, na.rm = TRUE),
    n_ok = dplyr::n(),
    .groups = "drop"
  )

saveRDS(paper_summary_snr_rhob, file.path(rds_dir, "paper_summary_snr_rhob.rds"))
saveRDS(paper_summary_setting, file.path(rds_dir, "paper_summary_setting.rds"))
saveRDS(paper_summary_overall, file.path(rds_dir, "paper_summary_overall.rds"))

readr::write_csv(paper_summary_snr_rhob, file.path(tables_dir, "paper_summary_snr_rhob.csv"))
readr::write_csv(paper_summary_setting, file.path(tables_dir, "paper_summary_setting.csv"))
readr::write_csv(paper_summary_overall, file.path(tables_dir, "paper_summary_overall.csv"))

log_msg("Summary by setting/rho_b/SNR:")
print(paper_summary_snr_rhob, n = Inf)

log_msg("Overall summary:")
print(paper_summary_overall, n = Inf)

log_msg("Simulation script finished successfully")
