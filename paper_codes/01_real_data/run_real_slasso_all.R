# ============================================================
# Run all real-data sglasso benchmarks on server
# Datasets:
#   1) bardet  from gglasso
#   2) Birthwt from grpreg
#   3) GenAtHum from sglasso
#
# Outputs are saved automatically after each dataset.
# ============================================================

options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

user_lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_lib, .libPaths()))
}

# -----------------------------
# 0. User/server options
# -----------------------------

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}

NREP <- as.integer(get_arg("nrep", "100"))
SEED <- as.integer(get_arg("seed", "2025"))
OUTDIR <- get_arg("outdir", "results/real_data")

slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA)
if (!is.na(slurm_cpus) && nzchar(slurm_cpus)) {
  NCORES <- max(1L, as.integer(slurm_cpus) - 1L)
} else {
  NCORES <- max(1L, parallel::detectCores() - 1L)
}

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat("\n============================================================\n")
cat("Real-data grouped penalized regression benchmark\n")
cat("Started:", format(Sys.time()), "\n")
cat("nrep   :", NREP, "\n")
cat("seed   :", SEED, "\n")
cat("ncores :", NCORES, "\n")
cat("outdir :", OUTDIR, "\n")
cat("============================================================\n\n")

# -----------------------------
# 1. Source benchmark function
# -----------------------------

cmd_args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args_all, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
} else {
  getwd()
}

function_candidates <- c(
  file.path(script_dir, "sglasso_real_function.R"),
  file.path(getwd(), "paper_codes", "01_real_data", "sglasso_real_function.R"),
  file.path(getwd(), "01_real_data", "sglasso_real_function.R"),
  file.path(getwd(), "sglasso_real_function.R")
)
function_file <- function_candidates[file.exists(function_candidates)][1L]
if (is.na(function_file)) {
  stop("Cannot find sglasso_real_function.R near the run script or current directory.")
}
source(function_file)

if (!exists("real_slasso")) {
  stop("Function real_slasso() was not found after sourcing ", function_file)
}

# -----------------------------
# 2. Packages and helpers
# -----------------------------

needed_pkgs <- c(
  "tibble",
  "gglasso",
  "grpreg",
  "grpnet",
  "adelie",
  "foreach",
  "doParallel",
  "doRNG",
  "sglasso"
)
for (pkg in needed_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package not available: ", pkg,
         ". Install it before running the server job.")
  }
}

summarise_results <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(data.frame())

  numeric_cols <- c(
    "time",
    "iter",
    "mse",
    "mae",
    "alpha",
    "lambda",
    "d",
    "selected_groups",
    "selected_features",
    "sparsity"
  )
  out <- do.call(
    rbind,
    lapply(split(df, df$method), function(z) {
      row <- data.frame(
        method = as.character(z$method[1]),
        nrep_success = length(unique(z$rep)),
        stringsAsFactors = FALSE
      )

      for (cc in numeric_cols) {
        if (cc %in% names(z)) {
          row[[paste0(cc, "_mean")]] <- mean(z[[cc]], na.rm = TRUE)
          row[[paste0(cc, "_sd")]] <- stats::sd(z[[cc]], na.rm = TRUE)
          row[[paste0(cc, "_median")]] <- stats::median(z[[cc]], na.rm = TRUE)
        }
      }
      row
    })
  )

  rownames(out) <- NULL
  out
}

save_one_dataset <- function(dataset_name, X, y, group, nrep = NREP, seed = SEED) {
  cat("\n------------------------------------------------------------\n")
  cat("Dataset:", dataset_name, "\n")
  cat("n      :", nrow(X), "\n")
  cat("p      :", ncol(X), "\n")
  cat("groups :", length(unique(group)), "\n")
  cat("start  :", format(Sys.time()), "\n")
  cat("------------------------------------------------------------\n")

  dataset_dir <- file.path(OUTDIR, dataset_name)
  dir.create(dataset_dir, showWarnings = FALSE, recursive = TRUE)

  meta <- list(
    dataset = dataset_name,
    n = nrow(X),
    p = ncol(X),
    ngroups = length(unique(group)),
    group_table = table(group),
    nrep = nrep,
    seed = seed,
    ncores = NCORES,
    start_time = Sys.time()
  )

  saveRDS(meta, file = file.path(dataset_dir, paste0(dataset_name, "_meta.rds")))

  result <- tryCatch(
    {
      real_slasso(
        X = X,
        y = y,
        group = group,
        n = nrow(X),
        nrep = nrep,
        ncores = NCORES,
        seed = seed
      )
    },
    error = function(e) {
      cat("ERROR in", dataset_name, ":", conditionMessage(e), "\n")
      err <- list(
        dataset = dataset_name,
        message = conditionMessage(e),
        call = conditionCall(e),
        time = Sys.time()
      )
      saveRDS(err, file = file.path(dataset_dir, paste0(dataset_name, "_ERROR.rds")))
      return(NULL)
    }
  )

  if (is.null(result)) {
    cat("Dataset failed:", dataset_name, "\n")
    return(invisible(NULL))
  }

  meta$end_time <- Sys.time()
  meta$elapsed_minutes <- as.numeric(difftime(meta$end_time, meta$start_time, units = "mins"))

  summary_result <- summarise_results(result)

  rds_file <- file.path(dataset_dir, paste0(dataset_name, "_raw_results.rds"))
  csv_file <- file.path(dataset_dir, paste0(dataset_name, "_raw_results.csv"))
  summary_rds_file <- file.path(dataset_dir, paste0(dataset_name, "_summary.rds"))
  summary_csv_file <- file.path(dataset_dir, paste0(dataset_name, "_summary.csv"))
  meta_file <- file.path(dataset_dir, paste0(dataset_name, "_meta_final.rds"))

  saveRDS(result, rds_file)
  write.csv(result, csv_file, row.names = FALSE)

  saveRDS(summary_result, summary_rds_file)
  write.csv(summary_result, summary_csv_file, row.names = FALSE)

  saveRDS(meta, meta_file)

  cat("Saved raw RDS    :", rds_file, "\n")
  cat("Saved raw CSV    :", csv_file, "\n")
  cat("Saved summary CSV:", summary_csv_file, "\n")
  cat("elapsed minutes  :", round(meta$elapsed_minutes, 3), "\n")
  cat("end              :", format(Sys.time()), "\n")

  invisible(result)
}

# -----------------------------
# 3. Build datasets
# -----------------------------

datasets <- list()

# 3.1 Bardet from gglasso
data(bardet, package = "gglasso")
X_bardet <- as.matrix(bardet$x)
y_bardet <- as.numeric(bardet$y)
group_bardet <- rep(1:20, each = 5)

if (ncol(X_bardet) != length(group_bardet)) {
  stop("Bardet group length does not match ncol(X).")
}

datasets$bardet <- list(
  X = X_bardet,
  y = y_bardet,
  group = group_bardet
)

# 3.2 Birthwt from grpreg
data(Birthwt, package = "grpreg")
X_birthwt <- as.matrix(Birthwt$X)
y_birthwt <- as.numeric(Birthwt$bwt)
K_birthwt <- as.integer(table(Birthwt$group))
group_birthwt <- rep(seq_along(K_birthwt), K_birthwt)

if (ncol(X_birthwt) != length(group_birthwt)) {
  stop("Birthwt group length does not match ncol(X).")
}

datasets$Birthwt <- list(
  X = X_birthwt,
  y = y_birthwt,
  group = group_birthwt
)

# 3.3 GenAtHum from sglasso
data(GenAtHum, package = "sglasso")
X_gen <- as.matrix(GenAtHum$X)
y_gen <- as.numeric(GenAtHum$y)
group_gen <- as.numeric(GenAtHum$group)

if (ncol(X_gen) != length(group_gen)) {
  stop("GenAtHum group length does not match ncol(X).")
}

datasets$GenAtHum <- list(
  X = X_gen,
  y = y_gen,
  group = group_gen
)

# -----------------------------
# 4. Run all datasets
# -----------------------------

all_status <- data.frame(
  dataset = character(),
  status = character(),
  n = integer(),
  p = integer(),
  ngroups = integer(),
  stringsAsFactors = FALSE
)

for (nm in names(datasets)) {
  dat <- datasets[[nm]]

  ans <- save_one_dataset(
    dataset_name = nm,
    X = dat$X,
    y = dat$y,
    group = dat$group,
    nrep = NREP,
    seed = SEED
  )

  all_status <- rbind(
    all_status,
    data.frame(
      dataset = nm,
      status = if (is.null(ans)) "failed" else "completed",
      n = nrow(dat$X),
      p = ncol(dat$X),
      ngroups = length(unique(dat$group)),
      stringsAsFactors = FALSE
    )
  )

  write.csv(all_status, file.path(OUTDIR, "run_status.csv"), row.names = FALSE)
  saveRDS(all_status, file.path(OUTDIR, "run_status.rds"))
}

# -----------------------------
# 5. Save combined outputs
# -----------------------------

raw_files <- file.path(OUTDIR, names(datasets), paste0(names(datasets), "_raw_results.rds"))
raw_files <- raw_files[file.exists(raw_files)]

if (length(raw_files) > 0) {
  combined <- do.call(
    rbind,
    lapply(raw_files, function(f) {
      x <- readRDS(f)
      x$dataset <- basename(dirname(f))
      x
    })
  )

  combined <- combined[, c("dataset", setdiff(names(combined), "dataset"))]

  saveRDS(combined, file.path(OUTDIR, "ALL_raw_results.rds"))
  write.csv(combined, file.path(OUTDIR, "ALL_raw_results.csv"), row.names = FALSE)

  combined_summary <- do.call(
    rbind,
    lapply(split(combined, combined$dataset), function(z) {
      ss <- summarise_results(z)
      ss$dataset <- z$dataset[1]
      ss[, c("dataset", setdiff(names(ss), "dataset"))]
    })
  )

  saveRDS(combined_summary, file.path(OUTDIR, "ALL_summary.rds"))
  write.csv(combined_summary, file.path(OUTDIR, "ALL_summary.csv"), row.names = FALSE)
}

sink_file <- file.path(OUTDIR, "sessionInfo.txt")
writeLines(capture.output(sessionInfo()), sink_file)

cat("\n============================================================\n")
cat("All datasets finished:", format(Sys.time()), "\n")
cat("Outputs saved in:", OUTDIR, "\n")
cat("============================================================\n")
