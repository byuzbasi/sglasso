############################################################
# SGLASSO tuning-criterion sensitivity analysis
# Only SGLASSO is fitted
# Tuning criteria are selected using sglasso::select()
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
    normalizePath("03_tuning_sensitivity/run_sglasso_tuning_sensitivity.R",
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

log_msg <- function(...) {
  cat(
    sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf(...),
    "\n"
  )
  flush.console()
}

############################################################
# Settings
############################################################

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}

repeatnum <- as.integer(get_arg("nrep", "5"))
seed <- as.integer(get_arg("seed", "2026"))
OUTDIR <- get_arg("outdir", "results/tuning_sensitivity")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

slurm_cpus <- suppressWarnings(
  as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
)
default_cores <- if (!is.na(slurm_cpus) && slurm_cpus > 1L) {
  slurm_cpus - 1L
} else {
  1L
}
cores <- as.integer(get_arg("cores", as.character(default_cores)))
cores <- max(1L, cores)

# Small test settings first.
# For final run, use: paper_settings <- get_paper_settings()
paper_settings <- get_paper_settings()[c(2, 4), ]

rho_b_grid <- c(0.3, 0.9)

snr_grid <- c(
  get_snr_grid()[5]
)

signal_pattern_grid <- c(
  "weak_strong_mixed"
)

tuning_criterion_grid <- c(
  "Min_val",
  "AIC",
  "BIC",
  "EBIC",
  "GCV"
)

alpha_seq <- seq(0.1, 0.9, by = 0.1)
d_seq <- seq(0, 1, by = 0.1)
nlambda <- 100

############################################################
# Task grid
############################################################

task_grid <- expand.grid(
  setting_id = seq_len(nrow(paper_settings)),
  tuning_criterion = tuning_criterion_grid,
  rho_b = rho_b_grid,
  snr = snr_grid,
  signal_pattern = signal_pattern_grid,
  rep_id = seq_len(repeatnum),
  stringsAsFactors = FALSE
)

log_msg("Total SGLASSO tuning tasks: %d", nrow(task_grid))
log_msg("Settings: nrep=%d, seed=%d, cores=%d", repeatnum, seed, cores)

############################################################
# One SGLASSO-only task
############################################################

run_one_sglasso_tuning_task <- function(ii) {
  
  ##########################################################
  # Extract one task safely
  ##########################################################
  
  task <- task_grid[ii, , drop = FALSE]
  st <- paper_settings[task$setting_id[[1]], ]
  
  criterion_i <- as.character(task$tuning_criterion[[1]])
  rho_b_i <- as.numeric(task$rho_b[[1]])
  snr_i <- as.numeric(task$snr[[1]])
  signal_i <- as.character(task$signal_pattern[[1]])
  rep_i <- as.integer(task$rep_id[[1]])
  
  set.seed(seed + ii)
  
  ##########################################################
  # Generate data
  ##########################################################
  
  nonzero_id <- sample(st$J, st$strong_J)
  
  dat <- block_sim_data(
    n_train = st$n_train,
    n_val = st$n_val,
    n_test = st$n_test,
    pj = st$pj,
    J = st$J,
    strong_J = st$strong_J,
    rho_w = 0.9,
    rho_b = rho_b_i,
    eff_nonzero = 1,
    corrmat_type = "Exchangeable",
    snr = snr_i,
    nonzero_id = nonzero_id,
    signal_pattern = signal_i
  )
  
  std_dat <- standardize_by_train(
    dat$X_train,
    dat$X_val,
    dat$X_test,
    dat$y_train,
    standardize = TRUE
  )
  
  COV_X_test <- cov(std_dat$X_test)
  
  ##########################################################
  # Fit SGLASSO for each alpha separately
  # sglasso() expects a single alpha value, not alpha_seq.
  ##########################################################
  
  fits <- lapply(alpha_seq, function(a) {
    sglasso::sglasso(
      X = std_dat$X_train,
      Y = std_dat$y_train,
      group = dat$group,
      nlambda = nlambda,
      alpha = a,
      d = d_seq,
      eps = 1e-5,
      max_iter = 1e5,
      lambda.min = 0.05
    )
  })
  
  ##########################################################
  # Case 1: Min_val
  # Select by validation MSE using compute_tuning_score()
  ##########################################################
  
  if (identical(criterion_i, "Min_val")) {
    
    best <- list(
      score = Inf,
      alpha_id = NA_integer_,
      lambda_id = NA_integer_,
      d_id = NA_integer_,
      beta = NULL,
      fit = NULL
    )
    
    for (ai in seq_along(fits)) {
      for (di in seq_along(d_seq)) {
        
        B <- fits[[ai]]$betas[, , di, drop = FALSE][, , 1]
        
        score <- compute_tuning_score(
          B = B,
          dat = dat,
          std_dat = std_dat,
          criterion = "Min_val"
        )
        
        score[!is.finite(score)] <- Inf
        
        if (all(is.infinite(score))) {
          next
        }
        
        li <- which.min(score)
        
        if (score[li] < best$score) {
          best <- list(
            score = score[li],
            alpha_id = ai,
            lambda_id = li,
            d_id = di,
            beta = B[, li],
            fit = fits[[ai]]
          )
        }
      }
    }
    
    if (is.null(best$beta)) {
      
      res <- empty_method_result(
        method = "SGLASSO",
        message = "Min_val failed: all tuning scores are non-finite.",
        tuning_criterion = criterion_i
      )
      
    } else {
      
      res <- evaluate_beta(
        method = "SGLASSO",
        beta_hat = best$beta,
        beta_true = dat$beta,
        X_test_std = std_dat$X_test,
        y_test = dat$y_test,
        y_center = std_dat$y_center,
        group = dat$group,
        true_groups = dat$true_groups,
        corrmat = dat$corrmat,
        COV_X_test = COV_X_test,
        sigma2 = dat$sigma^2,
        error_null = dat$error_null,
        risk_null = dat$risk_null,
        alpha = alpha_seq[best$alpha_id],
        lambda = safe_lambda(best$fit$lambda, best$lambda_id),
        d = d_seq[best$d_id],
        tuning_criterion = criterion_i,
        time = NA_real_,
        iter = NA_real_
      )
    }
    
    ##########################################################
    # Case 2: AIC/BIC/EBIC/GCV
    # Select by package-level sglasso::select()
    ##########################################################
    
  } else {
    
    best <- list(
      score = Inf,
      alpha_id = NA_integer_,
      beta = NULL,
      lambda = NA_real_,
      d = NA_real_
    )
    
    for (ai in seq_along(fits)) {
      
      sel <- tryCatch(
        {
          sglasso::select(
            fits[[ai]],
            criterion = criterion_i,
            ebic_gamma = 0.5,
            ebic_level = "group"
          )
        },
        error = function(e) e
      )
      
      if (inherits(sel, "error")) {
        next
      }
      
      ic_min <- min(sel$IC, na.rm = TRUE)
      
      if (is.finite(ic_min) && ic_min < best$score) {
        best <- list(
          score = ic_min,
          alpha_id = ai,
          beta = sel$beta,
          lambda = sel$lambda,
          d = sel$d
        )
      }
    }
    
    if (is.null(best$beta)) {
      
      res <- empty_method_result(
        method = "SGLASSO",
        message = paste0(criterion_i, " failed for all alpha values."),
        tuning_criterion = criterion_i
      )
      
    } else {
      
      res <- evaluate_beta(
        method = "SGLASSO",
        beta_hat = best$beta,
        beta_true = dat$beta,
        X_test_std = std_dat$X_test,
        y_test = dat$y_test,
        y_center = std_dat$y_center,
        group = dat$group,
        true_groups = dat$true_groups,
        corrmat = dat$corrmat,
        COV_X_test = COV_X_test,
        sigma2 = dat$sigma^2,
        error_null = dat$error_null,
        risk_null = dat$risk_null,
        alpha = alpha_seq[best$alpha_id],
        lambda = best$lambda,
        d = best$d,
        tuning_criterion = criterion_i,
        time = NA_real_,
        iter = NA_real_
      )
    }
  }
  
  ##########################################################
  # Attach scenario information
  ##########################################################
  
  res$rep <- rep_i
  res$setting <- st$setting
  res$dimensionality <- st$dimensionality
  res$sparsity_level <- st$sparsity_level
  
  res$n_train <- st$n_train
  res$n_val <- st$n_val
  res$n_test <- st$n_test
  res$p <- st$p
  res$pj <- st$pj
  res$J <- st$J
  res$strong_J <- st$strong_J
  
  res$rho_w <- 0.9
  res$rho_b <- rho_b_i
  res$snr <- snr_i
  res$signal_pattern <- signal_i
  res$corrmat_type <- "Exchangeable"
  res$standardize <- TRUE
  
  res
}

############################################################
# Run all tasks
############################################################

if (cores > 1) {
  
  doParallel::registerDoParallel(cores = cores)
  doRNG::registerDoRNG(seed)
  
  all_tune <- foreach::foreach(
    ii = seq_len(nrow(task_grid)),
    .combine = dplyr::bind_rows,
    .packages = required_packages,
    .export = c(
      "task_grid",
      "paper_settings",
      "seed",
      "alpha_seq",
      "d_seq",
      "nlambda",
      "run_one_sglasso_tuning_task"
    )
  ) %dopar% {
    run_one_sglasso_tuning_task(ii)
  }
  
  foreach::registerDoSEQ()
  
} else {
  
  all_tune <- dplyr::bind_rows(
    lapply(seq_len(nrow(task_grid)), run_one_sglasso_tuning_task)
  )
}

############################################################
# Save raw results
############################################################

saveRDS(
  all_tune,
  file.path(OUTDIR, "sglasso_tuning_all_raw.rds")
)

readr::write_csv(
  all_tune,
  file.path(OUTDIR, "sglasso_tuning_all_raw.csv")
)

############################################################
# Summaries
############################################################

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
  file.path(OUTDIR, "sglasso_tuning_summary_by_setting.csv")
)

readr::write_csv(
  summary_tune_scenario,
  file.path(OUTDIR, "sglasso_tuning_summary_by_scenario.csv")
)

readr::write_csv(
  summary_tune_overall,
  file.path(OUTDIR, "sglasso_tuning_summary_overall.csv")
)

saveRDS(
  summary_tune_setting,
  file.path(OUTDIR, "sglasso_tuning_summary_by_setting.rds")
)

saveRDS(
  summary_tune_scenario,
  file.path(OUTDIR, "sglasso_tuning_summary_by_scenario.rds")
)

saveRDS(
  summary_tune_overall,
  file.path(OUTDIR, "sglasso_tuning_summary_overall.rds")
)

############################################################
# Compact LaTeX table
############################################################

if (requireNamespace("kableExtra", quietly = TRUE)) {
  
  tab_tuning <- summary_tune_setting %>%
    dplyr::select(
      setting,
      tuning_criterion,
      MSE_y,
      MSE_beta,
      pve,
      FDR,
      selected_groups,
      selected_features
    ) %>%
    dplyr::arrange(
      setting,
      tuning_criterion
    )
  
  latex_table <- kableExtra::kbl(
    tab_tuning,
    format = "latex",
    digits = 3,
    booktabs = TRUE,
    longtable = FALSE,
    caption =
      "Sensitivity of SGLASSO to alternative tuning criteria."
  ) %>%
    kableExtra::kable_styling(
      latex_options = c(
        "hold_position",
        "scale_down"
      )
    )
  
  cat(
    latex_table,
    file = file.path(OUTDIR, "sglasso_tuning_table.tex")
  )
}

############################################################
# Print
############################################################

print(summary_tune_setting)
print(summary_tune_overall)

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
  file.path(OUTDIR, "summary_tune_scenario_table.csv")
)

if (requireNamespace("kableExtra", quietly = TRUE)) {
  
  latex_tab <- kableExtra::kbl(
    tab_reviewer,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    align = "c",
    caption = "\\label{tab:tuning_sensitivity} Sensitivity of SGLASSO performance to alternative tuning criteria."
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
  
  cat(latex_tab)
  
  cat(
    latex_tab,
    file = file.path(OUTDIR, "summary_tune_scenario_table.tex")
  )
}

log_msg("SGLASSO tuning sensitivity analysis completed")
