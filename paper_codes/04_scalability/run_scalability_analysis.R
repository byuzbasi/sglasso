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
    normalizePath("paper_codes/04_scalability/run_scalability_analysis.R",
                  mustWork = FALSE)
  }
}

r_library_dir <- Sys.getenv(
  "R_LIBRARY_DIR",
  unset = Sys.getenv("R_LIBS_USER", unset = "")
)
if (nzchar(r_library_dir)) {
  dir.create(r_library_dir, recursive = TRUE, showWarnings = FALSE)
  .libPaths(unique(c(r_library_dir, .libPaths())))
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

parse_logical <- function(x) {
  tolower(x) %in% c("true", "t", "1", "yes", "y")
}

parse_optional_integer <- function(x) {
  if (tolower(x) %in% c("", "na", "nan", "null", "none")) return(NA_integer_)
  as.integer(x)
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
CHECKPOINT_DIR <- file.path(OUTDIR, "checkpoints")
dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

p_grid <- parse_num_grid(
  get_arg("p_grid", "3000,10000,20000,50000")
)
J_grid <- parse_num_grid(
  get_arg("J_grid", "1000,2000,5000,10000")
)
repeatnum <- as.integer(get_arg("nrep", "5"))
seed <- as.integer(get_arg("seed", "2026"))
pj <- as.integer(get_arg("pj", "3"))
strong_frac <- as.numeric(get_arg("strong_frac", "0.10"))
active_groups <- parse_optional_integer(get_arg("active_groups", "10"))
alpha_fixed <- as.numeric(get_arg("alpha", "0.5"))
grpreg_alpha_fixed <- as.numeric(get_arg("grpreg_alpha", "1"))
d_fixed <- as.numeric(get_arg("d", "0.5"))
nlambda <- as.integer(get_arg("nlambda", "30"))
standardize <- parse_logical(get_arg("standardize", "FALSE"))
lambda_min_ratio <- as.numeric(get_arg("lambda_min_ratio", "0.01"))
eps <- as.numeric(get_arg("eps", "1e-4"))
maxit <- as.numeric(get_arg("maxit", "1e5"))
sglasso_screen <- get_arg("sglasso_screen", "SSR_fast")
sglasso_transform <- get_arg("sglasso_transform", "lazy")
adelie_early_exit <- parse_logical(get_arg("adelie_early_exit", "FALSE"))
fast_data <- parse_logical(get_arg("fast_data", "TRUE"))
slurm_cpus <- suppressWarnings(
  as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
)
default_cores <- if (!is.na(slurm_cpus) && slurm_cpus > 1L) {
  slurm_cpus
} else {
  1L
}
cores <- as.integer(get_arg("cores", as.character(default_cores)))
cores <- max(1L, cores)

log_msg(
  "Settings: p_grid=%s; J_grid=%s; pj=%d; active_groups=%s; strong_frac=%.3f; nrep=%d; alpha=%.3f; grpreg_alpha=%.3f; d=%.3f; nlambda=%d; lambda_min_ratio=%.3f; eps=%g; standardize=%s; sglasso_screen=%s; sglasso_transform=%s; adelie_early_exit=%s; signal_pattern=%s; cores=%d; fast_data=%s",
  paste(p_grid, collapse = ","),
  paste(J_grid, collapse = ","),
  pj,
  ifelse(is.na(active_groups), "NA", as.character(active_groups)),
  strong_frac,
  repeatnum,
  alpha_fixed,
  grpreg_alpha_fixed,
  d_fixed,
  nlambda,
  lambda_min_ratio,
  eps,
  standardize,
  sglasso_screen,
  sglasso_transform,
  adelie_early_exit,
  "weak_strong_mixed",
  cores,
  fast_data
)
log_msg("Output directory: %s", OUTDIR)
log_msg("Task checkpoint directory: %s", CHECKPOINT_DIR)
log_msg(
  "Total scalability tasks: %d",
  (length(p_grid) + length(J_grid)) * repeatnum
)

scale_start <- proc.time()[3]
log_msg("Starting run_scalability_experiment()")
res_scalability <- run_scalability_experiment(
  p_grid = p_grid,
  J_grid = J_grid,
  pj = pj,
  
  n_train = 100,
  n_val = 0,
  n_test = 0,
  
  strong_frac = strong_frac,
  active_groups = active_groups,
  
  rho_w = 0.9,
  rho_b = 0.5,
  snr = 1.528,
  
  repeatnum = repeatnum,
  seed = seed,
  
  alpha_fixed = alpha_fixed,
  grpreg_alpha_fixed = grpreg_alpha_fixed,
  d_fixed = d_fixed,
  
  nlambda = nlambda,
  standardize = standardize,
  lambda_min_ratio = lambda_min_ratio,
  
  eps = eps,
  maxit = maxit,
  
  signal_pattern = "weak_strong_mixed",
  sglasso_screen = sglasso_screen,
  sglasso_transform = sglasso_transform,
  adelie_early_exit = adelie_early_exit,
  checkpoint_dir = CHECKPOINT_DIR,
  cores = cores,
  fast_exchangeable_data = fast_data
)
log_msg(
  "Finished run_scalability_experiment() | rows=%d | elapsed=%.2f min",
  nrow(res_scalability),
  (proc.time()[3] - scale_start) / 60
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
log_msg("Scalability outputs saved in %s", OUTDIR)
log_msg(
  "Scalability task checkpoints saved: %d file(s)",
  length(list.files(CHECKPOINT_DIR, pattern = "\\.rds$", full.names = TRUE))
)

print(summary_scalability)

log_msg("Scalability experiment completed")
