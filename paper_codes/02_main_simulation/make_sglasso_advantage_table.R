############################################################
# Advantageous regimes table for SGLASSO
############################################################
#
# This script summarizes the regimes where SGLASSO is most competitive
# in the main simulation.  It ranks methods within each
# setting x signal_pattern x rho_b x SNR cell, excluding ORACLE, and
# reports the proportion of cells where SGLASSO is among the top two
# methods for prediction error.
#
# Default input:
#   results/main_simulation_rep50_fixed/tables/paper_summary_snr_rhob.csv
#
# Default output:
#   results/main_simulation_rep50_fixed/report_tables/
#
# Usage from R:
#   source("paper_codes/02_main_simulation/make_sglasso_advantage_table.R")
#
# Usage from terminal:
#   Rscript --vanilla paper_codes/02_main_simulation/make_sglasso_advantage_table.R \
#     --results-dir=results/main_simulation_rep50_fixed

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

cmd_args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), cmd_args, value = TRUE)
  if (length(hit) == 0L) {
    return(default)
  }
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

results_dir <- get_arg(
  "results-dir",
  Sys.getenv("MAIN_SIM_OUTDIR", unset = "results/main_simulation_rep50_fixed")
)
summary_file <- file.path(results_dir, "tables", "paper_summary_snr_rhob.csv")
out_dir <- get_arg("out-dir", file.path(results_dir, "report_tables"))

if (!file.exists(summary_file)) {
  stop("Cannot find summary file: ", summary_file, call. = FALSE)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

format_pct <- function(x) {
  sprintf("%.1f", 100 * x)
}

format_num <- function(x, digits = 3) {
  sprintf(paste0("%.", digits, "f"), x)
}

latex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([&_#%])", "\\\\\\1", x, perl = TRUE)
  x
}

make_condition_text <- function(type, lower = NULL, upper = NULL, values = NULL) {
  if (!is.null(values)) {
    return(paste(values, collapse = ", "))
  }
  if (!is.null(lower) && !is.null(upper)) {
    return(sprintf("%s in [%.3g, %.3g]", type, lower, upper))
  }
  if (!is.null(lower)) {
    return(sprintf("%s >= %.3g", type, lower))
  }
  if (!is.null(upper)) {
    return(sprintf("%s <= %.3g", type, upper))
  }
  "All grid values"
}

summary_dat <- readr::read_csv(summary_file, show_col_types = FALSE) %>%
  filter(method != "ORACLE")

ranked <- summary_dat %>%
  group_by(setting, signal_pattern, rho_b, snr) %>%
  mutate(
    rank_MSE_y = min_rank(MSE_y),
    rank_MSE_beta = min_rank(MSE_beta),
    rank_pve = min_rank(desc(pve)),
    rank_FDR = min_rank(FDR)
  ) %>%
  ungroup()

sglasso_cells <- ranked %>%
  filter(method == "SGLASSO")

############################################################
# Regimes selected for the manuscript table
############################################################

regimes <- tibble::tribble(
  ~regime_id, ~design, ~signal_pattern, ~condition_label, ~condition_tex, ~filter_expr, ~interpretation,
  "LD-Hom-rho-high",
  "LD",
  "Homogeneous",
  make_condition_text("rho_b", lower = 0.7),
  "$\\rho_b \\ge 0.7$",
  "dimensionality == 'LD' & signal_pattern == 'homogeneous' & rho_b >= 0.7",
  "Coherent grouped signals with strong between-group correlation.",

  "HD-Hom-snr-mid-high",
  "HD",
  "Homogeneous",
  make_condition_text("SNR", lower = 0.196),
  "$\\mathrm{SNR} \\ge 0.196$",
  "dimensionality == 'HD' & signal_pattern == 'homogeneous' & snr >= 0.196",
  "High-dimensional coherent signals once the signal strength is moderate.",

  "HD-WSM-snr-mid-high",
  "HD",
  "Weak/strong mixed",
  make_condition_text("SNR", lower = 0.389, upper = 3.028),
  "$0.389 \\le \\mathrm{SNR} \\le 3.028$",
  "dimensionality == 'HD' & signal_pattern == 'weak_strong_mixed' & snr >= 0.389 & snr <= 3.028",
  "High-dimensional settings with heterogeneous signal magnitudes.",

  "HD-WSM-rho-mid-high",
  "HD",
  "Weak/strong mixed",
  make_condition_text("rho_b", lower = 0.5),
  "$\\rho_b \\ge 0.5$",
  "dimensionality == 'HD' & signal_pattern == 'weak_strong_mixed' & rho_b >= 0.5",
  "Correlated high-dimensional groups with weak and strong active effects.",

  "LD-WSM-rho-high",
  "LD",
  "Weak/strong mixed",
  "rho_b = 0.9",
  "$\\rho_b = 0.9$",
  "dimensionality == 'LD' & signal_pattern == 'weak_strong_mixed' & abs(rho_b - 0.9) < 1e-12",
  "A low-dimensional stress case where SGLASSO remains competitive under strong correlation."
)

summarise_regime <- function(filter_expr) {
  dat <- sglasso_cells %>%
    filter(eval(parse(text = filter_expr)))

  if (nrow(dat) == 0L) {
    return(tibble::tibble(
      n_cells = 0L,
      top2_MSE_y_rate = NA_real_,
      top3_MSE_y_rate = NA_real_,
      top2_pve_rate = NA_real_,
      top2_MSE_beta_rate = NA_real_,
      mean_pve = NA_real_,
      mean_FDR = NA_real_,
      mean_selected_groups = NA_real_,
      mean_selected_features = NA_real_
    ))
  }

  tibble::tibble(
    n_cells = nrow(dat),
    top2_MSE_y_rate = mean(dat$rank_MSE_y <= 2, na.rm = TRUE),
    top3_MSE_y_rate = mean(dat$rank_MSE_y <= 3, na.rm = TRUE),
    top2_pve_rate = mean(dat$rank_pve <= 2, na.rm = TRUE),
    top2_MSE_beta_rate = mean(dat$rank_MSE_beta <= 2, na.rm = TRUE),
    mean_pve = mean(dat$pve, na.rm = TRUE),
    mean_FDR = mean(dat$FDR, na.rm = TRUE),
    mean_selected_groups = mean(dat$selected_groups, na.rm = TRUE),
    mean_selected_features = mean(dat$selected_features, na.rm = TRUE)
  )
}

advantage_table <- regimes %>%
  rowwise() %>%
  mutate(summary = list(summarise_regime(filter_expr))) %>%
  tidyr::unnest(summary) %>%
  ungroup() %>%
  transmute(
    regime_id,
    Design = design,
    `Signal pattern` = signal_pattern,
    `Favorable regime` = condition_label,
    `Favorable regime (LaTeX)` = condition_tex,
    `Cells` = n_cells,
    `Top-2 MSE_y (%)` = format_pct(top2_MSE_y_rate),
    `Top-3 MSE_y (%)` = format_pct(top3_MSE_y_rate),
    `Top-2 PVE (%)` = format_pct(top2_pve_rate),
    `Top-2 MSE_beta (%)` = format_pct(top2_MSE_beta_rate),
    `Mean PVE` = format_num(mean_pve),
    `Mean FDR` = format_num(mean_FDR),
    `Mean selected groups` = format_num(mean_selected_groups, 1),
    Interpretation = interpretation
  )

csv_file <- file.path(out_dir, "sglasso_advantage_regimes.csv")
tex_file <- file.path(out_dir, "sglasso_advantage_regimes.tex")

readr::write_csv(advantage_table, csv_file)

latex_table <- advantage_table %>%
  select(
    Design,
    `Signal pattern`,
    `Favorable regime (LaTeX)`,
    `Top-2 MSE_y (%)`,
    `Top-2 PVE (%)`,
    `Mean PVE`,
    `Mean FDR`,
    `Mean selected groups`
  )

latex_lines <- c(
  "\\begin{table}[!ht]",
  "\\centering",
  "\\caption{Advantageous simulation regimes for SGLASSO.}",
  "\\label{tab:sglasso_advantage_regimes}",
  "\\begin{tabular}{lllccccc}",
  "\\toprule",
  "Design & Signal pattern & Favorable regime & Top-2 $\\mathrm{MSE}_y$ (\\%) & Top-2 PVE (\\%) & Mean PVE & Mean FDR & $\\widehat{s}$ \\\\",
  "\\midrule"
)

for (ii in seq_len(nrow(latex_table))) {
  row <- latex_table[ii, ]
  latex_lines <- c(
    latex_lines,
    paste(
      latex_escape(as.character(row$Design)),
      latex_escape(as.character(row$`Signal pattern`)),
      as.character(row$`Favorable regime (LaTeX)`),
      latex_escape(as.character(row$`Top-2 MSE_y (%)`)),
      latex_escape(as.character(row$`Top-2 PVE (%)`)),
      latex_escape(as.character(row$`Mean PVE`)),
      latex_escape(as.character(row$`Mean FDR`)),
      latex_escape(as.character(row$`Mean selected groups`)),
      sep = " & "
    ) |>
      paste0(" \\\\")
  )
}

latex_lines <- c(
  latex_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{flushleft}",
  "\\footnotesize Note: Methods are ranked within each setting, signal pattern, $\\rho_b$, and SNR cell, excluding ORACLE. Top-2 rates report the percentage of cells in which SGLASSO is among the two best methods. $\\widehat{s}$ denotes the mean number of selected groups.",
  "\\end{flushleft}",
  "\\end{table}"
)

writeLines(latex_lines, tex_file)

cat("Advantage table written:\n")
cat(" - ", csv_file, "\n", sep = "")
cat(" - ", tex_file, "\n", sep = "")
cat("\nPreview:\n")
print(advantage_table, n = Inf)
