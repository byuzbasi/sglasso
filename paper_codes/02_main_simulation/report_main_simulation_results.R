############################################################
# Main simulation result report
############################################################
#
# This script reads the saved main-simulation outputs and creates
# compact diagnostic/report tables for manuscript revision.
#
# Default input:
#   results/main_simulation_rep50_fixed
#
# Default output:
#   results/main_simulation_rep50_fixed/report_tables
#
# Usage from R:
#   source("paper_codes/02_main_simulation/report_main_simulation_results.R")
#
# Usage from terminal:
#   Rscript --vanilla paper_codes/02_main_simulation/report_main_simulation_results.R \
#     --results-dir=results/main_simulation_rep50_fixed

cmd_args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), cmd_args, value = TRUE)
  if (length(hit) == 0L) {
    return(default)
  }
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

parse_bool <- function(x) {
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

results_dir <- get_arg("results-dir", "results/main_simulation_rep50_fixed")
tables_dir <- file.path(results_dir, "tables")
out_dir <- get_arg("out-dir", file.path(results_dir, "report_tables"))
top_n <- as.integer(get_arg("top-n", "2"))
expected_patterns <- strsplit(
  get_arg("expected-patterns", "homogeneous,weak_strong_mixed"),
  ",",
  fixed = TRUE
)[[1L]]
include_oracle_in_ranks <- parse_bool(get_arg("include-oracle-in-ranks", "false"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_required_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path, call. = FALSE)
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

write_report_csv <- function(x, filename) {
  path <- file.path(out_dir, filename)
  utils::write.csv(x, path, row.names = FALSE)
  path
}

rank_metric <- function(x, decreasing = FALSE) {
  if (decreasing) {
    rank(-x, ties.method = "min", na.last = "keep")
  } else {
    rank(x, ties.method = "min", na.last = "keep")
  }
}

safe_divide <- function(num, den) {
  out <- rep(NA_real_, length(num))
  ok <- !is.na(num) & !is.na(den) & den != 0
  out[ok] <- num[ok] / den[ok]
  out
}

cat("Main simulation report\n")
cat("Results directory:", results_dir, "\n")
cat("Tables directory :", tables_dir, "\n")
cat("Output directory :", out_dir, "\n\n")

all_results <- read_required_csv(file.path(tables_dir, "all_sim_slasso_results.csv"))
summary_overall <- read_required_csv(file.path(tables_dir, "paper_summary_overall.csv"))
summary_setting <- read_required_csv(file.path(tables_dir, "paper_summary_setting.csv"))

scenario_grid_path <- file.path(tables_dir, "simulation_scenario_grid.csv")
scenario_grid <- if (file.exists(scenario_grid_path)) {
  utils::read.csv(scenario_grid_path, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  NULL
}

############################################################
# 1. Output integrity and signal-pattern audit
############################################################

realized_patterns <- as.data.frame(table(all_results$signal_pattern), stringsAsFactors = FALSE)
names(realized_patterns) <- c("signal_pattern", "realized_rows")

expected_pattern_df <- data.frame(
  signal_pattern = expected_patterns,
  expected = TRUE,
  stringsAsFactors = FALSE
)

pattern_audit <- merge(
  expected_pattern_df,
  realized_patterns,
  by = "signal_pattern",
  all = TRUE
)
pattern_audit$expected[is.na(pattern_audit$expected)] <- FALSE
pattern_audit$realized_rows[is.na(pattern_audit$realized_rows)] <- 0L
pattern_audit$realized <- pattern_audit$realized_rows > 0L
pattern_audit$status <- ifelse(
  pattern_audit$expected & pattern_audit$realized,
  "present",
  ifelse(pattern_audit$expected & !pattern_audit$realized, "missing", "unexpected")
)

if (!is.null(scenario_grid) && "signal_pattern" %in% names(scenario_grid)) {
  planned <- as.data.frame(table(scenario_grid$signal_pattern), stringsAsFactors = FALSE)
  names(planned) <- c("signal_pattern", "planned_grid_rows")
  pattern_audit <- merge(pattern_audit, planned, by = "signal_pattern", all = TRUE)
  pattern_audit$planned_grid_rows[is.na(pattern_audit$planned_grid_rows)] <- 0L
}

key_metrics <- intersect(
  c("risk", "pe_mse", "MSE_y", "MSE_beta", "pve", "FDR",
    "selected_groups", "selected_features", "time"),
  names(all_results)
)

integrity_summary <- data.frame(
  item = c(
    "n_rows",
    "n_columns",
    "n_settings",
    "n_methods",
    "n_signal_patterns_realized",
    "n_missing_expected_signal_patterns",
    "n_key_metric_na"
  ),
  value = c(
    nrow(all_results),
    ncol(all_results),
    length(unique(all_results$setting)),
    length(unique(all_results$method)),
    length(unique(all_results$signal_pattern)),
    sum(pattern_audit$status == "missing"),
    sum(is.na(all_results[, key_metrics, drop = FALSE]))
  )
)

############################################################
# 2. Overall method rankings
############################################################

rank_pool <- summary_overall
if (!include_oracle_in_ranks && "method" %in% names(rank_pool)) {
  rank_pool <- rank_pool[rank_pool$method != "ORACLE", , drop = FALSE]
}

overall_ranking <- rank_pool
overall_ranking$rank_MSE_y <- rank_metric(overall_ranking$MSE_y)
overall_ranking$rank_MSE_beta <- rank_metric(overall_ranking$MSE_beta)
overall_ranking$rank_pve <- rank_metric(overall_ranking$pve, decreasing = TRUE)
overall_ranking$rank_FDR <- rank_metric(overall_ranking$FDR)
overall_ranking <- overall_ranking[order(overall_ranking$rank_MSE_y, overall_ranking$method), ]

############################################################
# 3. SGLASSO ranks by setting and best methods by setting
############################################################

setting_pool <- summary_setting
if (!include_oracle_in_ranks && "method" %in% names(setting_pool)) {
  setting_pool <- setting_pool[setting_pool$method != "ORACLE", , drop = FALSE]
}

setting_split <- split(setting_pool, setting_pool$setting)
setting_ranked <- do.call(rbind, lapply(setting_split, function(d) {
  d$rank_MSE_y <- rank_metric(d$MSE_y)
  d$rank_MSE_beta <- rank_metric(d$MSE_beta)
  d$rank_pve <- rank_metric(d$pve, decreasing = TRUE)
  d$rank_FDR <- rank_metric(d$FDR)
  d
}))
rownames(setting_ranked) <- NULL

sglasso_setting_ranks <- setting_ranked[setting_ranked$method == "SGLASSO", , drop = FALSE]
sglasso_setting_ranks <- sglasso_setting_ranks[order(sglasso_setting_ranks$setting), ]

best_methods_by_setting <- do.call(rbind, lapply(split(setting_ranked, setting_ranked$setting), function(d) {
  d <- d[order(d$MSE_y, d$method), ]
  d[seq_len(min(top_n, nrow(d))), , drop = FALSE]
}))
rownames(best_methods_by_setting) <- NULL

############################################################
# 4. SGLASSO comparisons against competitors
############################################################

sg_overall <- summary_overall[summary_overall$method == "SGLASSO", , drop = FALSE]
competitors <- summary_overall[summary_overall$method != "SGLASSO", , drop = FALSE]
if (!include_oracle_in_ranks) {
  competitors <- competitors[competitors$method != "ORACLE", , drop = FALSE]
}

sglasso_vs_competitors_overall <- competitors
sglasso_vs_competitors_overall$MSE_y_over_SGLASSO <- safe_divide(
  sglasso_vs_competitors_overall$MSE_y,
  sg_overall$MSE_y
)
sglasso_vs_competitors_overall$SGLASSO_over_method_MSE_y <- safe_divide(
  sg_overall$MSE_y,
  sglasso_vs_competitors_overall$MSE_y
)
sglasso_vs_competitors_overall$time_over_SGLASSO <- safe_divide(
  sglasso_vs_competitors_overall$time,
  sg_overall$time
)
sglasso_vs_competitors_overall$SGLASSO_over_method_time <- safe_divide(
  sg_overall$time,
  sglasso_vs_competitors_overall$time
)
sglasso_vs_competitors_overall <- sglasso_vs_competitors_overall[
  order(sglasso_vs_competitors_overall$MSE_y),
]

############################################################
# 5. Save report tables
############################################################

written <- c(
  write_report_csv(integrity_summary, "integrity_summary.csv"),
  write_report_csv(pattern_audit, "signal_pattern_audit.csv"),
  write_report_csv(overall_ranking, "overall_method_ranking.csv"),
  write_report_csv(sglasso_setting_ranks, "sglasso_setting_ranks.csv"),
  write_report_csv(best_methods_by_setting, "best_methods_by_setting.csv"),
  write_report_csv(sglasso_vs_competitors_overall, "sglasso_vs_competitors_overall.csv")
)

############################################################
# 6. Console report
############################################################

cat("Integrity summary:\n")
print(integrity_summary, row.names = FALSE)

cat("\nSignal-pattern audit:\n")
print(pattern_audit, row.names = FALSE)

cat("\nOverall ranking by MSE_y", if (!include_oracle_in_ranks) "(ORACLE excluded)" else "", ":\n")
print(
  overall_ranking[, intersect(
    c("method", "MSE_y", "rank_MSE_y", "MSE_beta", "rank_MSE_beta",
      "pve", "rank_pve", "FDR", "rank_FDR", "selected_groups",
      "selected_features", "time"),
    names(overall_ranking)
  )],
  row.names = FALSE
)

cat("\nSGLASSO ranks by setting:\n")
print(
  sglasso_setting_ranks[, intersect(
    c("setting", "signal_pattern", "MSE_y", "rank_MSE_y", "MSE_beta",
      "rank_MSE_beta", "pve", "rank_pve", "FDR", "rank_FDR",
      "selected_groups", "selected_features", "time"),
    names(sglasso_setting_ranks)
  )],
  row.names = FALSE
)

cat("\nWritten report files:\n")
cat(paste0(" - ", written, collapse = "\n"), "\n")

if (any(pattern_audit$status == "missing")) {
  cat("\nWARNING: At least one expected signal pattern is missing from the realized results.\n")
  cat("Missing patterns:\n")
  print(pattern_audit$signal_pattern[pattern_audit$status == "missing"])
}
