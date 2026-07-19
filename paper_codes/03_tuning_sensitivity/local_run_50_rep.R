############################################################
# Local 50-replication run for SGLASSO tuning sensitivity
############################################################
#
# Purpose:
#   Run the reviewer-requested tuning-criterion sensitivity analysis
#   locally and save all raw outputs, summaries, and reviewer-ready
#   tables under:
#
#     results/tuning_sensitivity
#
# Scenarios are defined in run_sglasso_tuning_sensitivity.R:
#   - settings: LD-8 and HD-10
#   - rho_b: 0.3 and 0.9
#   - SNR: middle paper-grid value
#   - signal pattern: weak_strong_mixed
#   - criteria: Validation, AIC, BIC, EBIC, GCV
#
# Current computational settings:
#   - 50 replications
#   - fast exchangeable data generation
#   - lambda.min.ratio = 0.01
#   - eps = 1e-4
#   - nlambda = 100
#   - SGLASSO screen = SSR_fast
#   - SGLASSO transform = lazy
#
# Note:
#   This script intentionally writes to the same final paper-results
#   folder used by the server scripts, so rerunning it will refresh
#   those result files.

############################################################
# Paths and run settings
############################################################

project_dir <- normalizePath(Sys.getenv("SGLASSO_PROJECT_DIR", unset = getwd()))
setwd(project_dir)

script_file <- "paper_codes/03_tuning_sensitivity/run_sglasso_tuning_sensitivity.R"
source_file <- "paper_codes/03_tuning_sensitivity/sglasso_sim_function_full_tuning.R"
outdir <- "results/tuning_sensitivity"

nrep <- 50L
cores <- 1L
fast_data <- TRUE
lambda_min_ratio <- 0.01
eps <- 1e-4
maxit <- 1e5
nlambda <- 100L
sglasso_screen <- "SSR_fast"
sglasso_transform <- "lazy"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

############################################################
# Load local helper functions and packages
############################################################

source(source_file)
install_and_load_packages(install_missing = FALSE)

############################################################
# Run the tuning sensitivity script
############################################################

cmd_args <- c(
  "--vanilla",
  script_file,
  paste0("--nrep=", nrep),
  paste0("--cores=", cores),
  paste0("--outdir=", outdir),
  paste0("--fast_data=", fast_data),
  paste0("--lambda_min_ratio=", lambda_min_ratio),
  paste0("--eps=", eps),
  paste0("--maxit=", maxit),
  paste0("--nlambda=", nlambda),
  paste0("--sglasso_screen=", sglasso_screen),
  paste0("--sglasso_transform=", sglasso_transform)
)

cat("Running tuning sensitivity analysis with command:\n")
cat("Rscript", paste(cmd_args, collapse = " "), "\n\n")

status <- system2("Rscript", args = cmd_args)
if (!identical(status, 0L)) {
  stop("Tuning sensitivity run failed with exit status: ", status)
}

############################################################
# Read raw results and regenerate summaries
############################################################

raw_file <- file.path(outdir, "sglasso_tuning_all_raw.rds")
if (!file.exists(raw_file)) {
  stop("Raw tuning result file was not created: ", raw_file)
}

all_tune <- readRDS(raw_file)

summary_tune_setting <- all_tune %>%
  dplyr::filter(converged == TRUE) %>%
  dplyr::group_by(
    setting,
    dimensionality,
    sparsity_level,
    tuning_criterion
  ) %>%
  dplyr::summarise(
    risk = mean(risk, na.rm = TRUE),
    pe_mse = mean(pe_mse, na.rm = TRUE),
    MSE_y = mean(MSE_y, na.rm = TRUE),
    MSE_beta = mean(MSE_beta, na.rm = TRUE),
    pve = mean(pve, na.rm = TRUE),
    FDR = mean(FDR, na.rm = TRUE),
    selected_groups = mean(selected_groups, na.rm = TRUE),
    selected_features = mean(selected_features, na.rm = TRUE),
    alpha = median(alpha, na.rm = TRUE),
    lambda = median(lambda, na.rm = TRUE),
    d = median(d, na.rm = TRUE),
    n_ok = dplyr::n(),
    .groups = "drop"
  )

summary_tune_scenario <- all_tune %>%
  dplyr::filter(converged == TRUE) %>%
  dplyr::group_by(
    setting,
    tuning_criterion,
    rho_b,
    snr,
    signal_pattern
  ) %>%
  dplyr::summarise(
    risk = mean(risk, na.rm = TRUE),
    pe_mse = mean(pe_mse, na.rm = TRUE),
    MSE_y = mean(MSE_y, na.rm = TRUE),
    MSE_beta = mean(MSE_beta, na.rm = TRUE),
    pve = mean(pve, na.rm = TRUE),
    FDR = mean(FDR, na.rm = TRUE),
    selected_groups = mean(selected_groups, na.rm = TRUE),
    selected_features = mean(selected_features, na.rm = TRUE),
    alpha = median(alpha, na.rm = TRUE),
    lambda = median(lambda, na.rm = TRUE),
    d = median(d, na.rm = TRUE),
    n_ok = dplyr::n(),
    .groups = "drop"
  )

summary_tune_overall <- all_tune %>%
  dplyr::filter(converged == TRUE) %>%
  dplyr::group_by(tuning_criterion) %>%
  dplyr::summarise(
    risk = mean(risk, na.rm = TRUE),
    pe_mse = mean(pe_mse, na.rm = TRUE),
    MSE_y = mean(MSE_y, na.rm = TRUE),
    MSE_beta = mean(MSE_beta, na.rm = TRUE),
    pve = mean(pve, na.rm = TRUE),
    FDR = mean(FDR, na.rm = TRUE),
    selected_groups = mean(selected_groups, na.rm = TRUE),
    selected_features = mean(selected_features, na.rm = TRUE),
    alpha = median(alpha, na.rm = TRUE),
    lambda = median(lambda, na.rm = TRUE),
    d = median(d, na.rm = TRUE),
    n_ok = dplyr::n(),
    .groups = "drop"
  )

############################################################
# Save summaries
############################################################

readr::write_csv(
  summary_tune_setting,
  file.path(outdir, "sglasso_tuning_summary_by_setting.csv")
)
readr::write_csv(
  summary_tune_scenario,
  file.path(outdir, "sglasso_tuning_summary_by_scenario.csv")
)
readr::write_csv(
  summary_tune_overall,
  file.path(outdir, "sglasso_tuning_summary_overall.csv")
)

saveRDS(
  summary_tune_setting,
  file.path(outdir, "sglasso_tuning_summary_by_setting.rds")
)
saveRDS(
  summary_tune_scenario,
  file.path(outdir, "sglasso_tuning_summary_by_scenario.rds")
)
saveRDS(
  summary_tune_overall,
  file.path(outdir, "sglasso_tuning_summary_overall.rds")
)

############################################################
# Reviewer-ready compact table
############################################################

tab_reviewer <- summary_tune_scenario %>%
  dplyr::mutate(
    tuning_criterion = dplyr::recode(
      tuning_criterion,
      "Min_val" = "Validation"
    ),
    tuning_criterion = factor(
      tuning_criterion,
      levels = c("Validation", "AIC", "BIC", "EBIC", "GCV")
    ),
    MSE_y = round(MSE_y, 3),
    MSE_beta = round(MSE_beta, 3),
    pve = round(pve, 3),
    FDR = round(FDR, 3),
    selected_groups = round(selected_groups, 1),
    selected_features = round(selected_features, 1),
    alpha = round(alpha, 2),
    lambda = round(lambda, 3),
    d = round(d, 2)
  ) %>%
  dplyr::select(
    setting,
    rho_b,
    tuning_criterion,
    MSE_y,
    MSE_beta,
    pve,
    FDR,
    selected_groups,
    selected_features,
    alpha,
    lambda,
    d
  ) %>%
  dplyr::arrange(setting, rho_b, tuning_criterion)

colnames(tab_reviewer) <- c(
  "Setting",
  "$\\rho_b$",
  "Criterion",
  "$\\mathrm{MSE}_y$",
  "$\\mathrm{MSE}_{\\beta}$",
  "PVE",
  "FDR",
  "$\\widehat{s}$",
  "$\\widehat{p}$",
  "$\\alpha$",
  "$\\lambda$",
  "$d$"
)

readr::write_csv(
  tab_reviewer,
  file.path(outdir, "summary_tune_scenario_table.csv")
)

if (requireNamespace("kableExtra", quietly = TRUE)) {
  latex_tab <- kableExtra::kbl(
    tab_reviewer,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    align = "c",
    caption = paste0(
      "\\label{tab:tuning_sensitivity} ",
      "Sensitivity of SGLASSO performance to alternative tuning criteria."
    )
  ) %>%
    kableExtra::collapse_rows(
      columns = c(1, 2),
      valign = "middle",
      latex_hline = "full"
    ) %>%
    kableExtra::footnote(
      general = paste(
        "$\\\\widehat{s}$ is the average number of selected groups.",
        "$\\\\widehat{p}$ is the average number of selected predictors.",
        "Validation denotes tuning based on the validation-set prediction error."
      ),
      threeparttable = TRUE,
      escape = FALSE
    ) %>%
    kableExtra::kable_styling(
      latex_options = c("hold_position", "scale_down")
    )

  cat(
    latex_tab,
    file = file.path(outdir, "summary_tune_scenario_table.tex")
  )
}

############################################################
# Quick checks printed to console
############################################################

cat("\nCompleted local 50-rep tuning sensitivity run.\n")
cat("Results directory:", normalizePath(outdir), "\n\n")

cat("Main output files:\n")
print(list.files(outdir, full.names = TRUE))

cat("\nSummary by setting:\n")
print(summary_tune_setting)

cat("\nOverall summary:\n")
print(summary_tune_overall)
