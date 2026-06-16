############################################################
# Plot SGLASSO simulation results
#
# Input:
#   results/main_simulation/rds/*.rds
#
# Output:
#   results/main_simulation/plots/*.pdf
#   results/main_simulation/plots/*.eps
#
# Main figures:
#   1. HD generalization risk vs rho_b
#   2. HD PVE vs rho_b
#   3. HD FDR vs rho_b
#   4. HD selected groups vs SNR
#   5. LD generalization risk vs rho_b
#   6. LD PVE vs rho_b
#   7. LD FDR vs rho_b
#   8. LD selected groups vs SNR
############################################################

options(stringsAsFactors = FALSE)

############################################################
# 1. Packages
############################################################

req_pkgs <- c(
  "dplyr",
  "purrr",
  "readr",
  "ggplot2",
  "stringr",
  "forcats",
  "scales"
)

missing_pkgs <- req_pkgs[
  !sapply(req_pkgs, requireNamespace, quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nPlease install them before running this script."
  )
}

library(dplyr)
library(purrr)
library(readr)
library(ggplot2)
library(stringr)
library(forcats)
library(scales)

############################################################
# 2. Helper functions
############################################################

log_msg <- function(...) {
  cat(
    sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf(...),
    "\n",
    sep = ""
  )
  flush.console()
}

save_plot_both <- function(plot_obj,
                           filename,
                           width = 8.5,
                           height = 5.5,
                           out_dir = "results/main_simulation/plots") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  pdf_file <- file.path(out_dir, paste0(filename, ".pdf"))
  eps_file <- file.path(out_dir, paste0(filename, ".eps"))
  
  ggsave(
    filename = pdf_file,
    plot = plot_obj,
    width = width,
    height = height,
    device = cairo_pdf
  )
  
  ggsave(
    filename = eps_file,
    plot = plot_obj,
    width = width,
    height = height,
    device = cairo_ps
  )
  
  log_msg("Saved: %s", pdf_file)
  log_msg("Saved: %s", eps_file)
}

clean_method_names <- function(x) {
  dplyr::recode(
    as.character(x),
    "SGLASSO" = "SGLASSO",
    "GLASSO" = "GLASSO",
    "GENET" = "GENET",
    "GSCAD" = "GSCAD",
    "GMCP" = "GMCP",
    "ORACLE" = "ORACLE",
    .default = as.character(x)
  )
}

make_common_theme <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      strip.background = element_rect(fill = "grey90", colour = "grey40"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = base_size - 1),
      axis.text.y = element_text(size = base_size - 1),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
}

############################################################
# 3. Read RDS results
############################################################

main_outdir <- Sys.getenv("MAIN_SIM_OUTDIR", unset = "results/main_simulation")
rds_dir <- file.path(main_outdir, "rds")
tables_dir <- file.path(main_outdir, "tables")

if (!dir.exists(rds_dir)) {
  stop("Cannot find directory: ", rds_dir)
}

rds_files <- list.files(
  rds_dir,
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(rds_files) == 0) {
  stop("No .rds files found in: ", rds_dir)
}

log_msg("Found %d RDS files.", length(rds_files))

# Prefer combined all-results file if available.
combined_file <- file.path(rds_dir, "all_sim_slasso_results.rds")

if (file.exists(combined_file)) {
  log_msg("Reading combined file: %s", combined_file)
  all_res <- readRDS(combined_file)
} else {
  log_msg("Combined file not found. Reading and binding all RDS files.")
  all_res <- purrr::map_dfr(rds_files, readRDS)
}

############################################################
# 4. Basic cleaning
############################################################

required_cols <- c(
  "setting",
  "dimensionality",
  "method",
  "rho_b",
  "snr",
  "signal_pattern",
  "converged",
  "risk",
  "MSE_y",
  "MSE_beta",
  "pve",
  "FDR",
  "selected_groups",
  "selected_features"
)

missing_cols <- setdiff(required_cols, names(all_res))

if (length(missing_cols) > 0) {
  stop(
    "The following required columns are missing from results: ",
    paste(missing_cols, collapse = ", ")
  )
}

all_res <- all_res %>%
  mutate(
    method = clean_method_names(method),
    method = factor(
      method,
      levels = c("SGLASSO", "ADELIE", "GENET", "GLASSO", "GSCAD", "GMCP", "ORACLE")
    ),
    setting = factor(
      setting,
      levels = c("LD-4", "LD-8", "LD-16", "HD-20", "HD-40")
    ),
    dimensionality = factor(dimensionality, levels = c("LD", "HD")),
    rho_b = as.numeric(rho_b),
    snr = as.numeric(snr),
    snr_label = paste0("SNR = ", round(snr, 3)),
    rho_b_label = paste0(expression(rho[b]), " = ", round(rho_b, 1)),
    signal_pattern = factor(
      signal_pattern,
      levels = c("homogeneous", "mixed_signs", "weak_strong_mixed"),
      labels = c("Homogeneous", "Mixed signs", "Weak/strong mixed")
    )
  )

all_res <- all_res %>%
  filter(
    converged == TRUE,
    !is.na(method),
    method != "ORACLE"
  )

log_msg("Rows after filtering converged non-oracle methods: %d", nrow(all_res))

############################################################
# 5. Summary tables
############################################################

summary_snr_rhob <- all_res %>%
  group_by(
    setting,
    dimensionality,
    signal_pattern,
    rho_b,
    snr,
    snr_label,
    method
  ) %>%
  summarise(
    risk = mean(risk, na.rm = TRUE),
    MSE_y = mean(MSE_y, na.rm = TRUE),
    MSE_beta = mean(MSE_beta, na.rm = TRUE),
    pve = mean(pve, na.rm = TRUE),
    FDR = mean(FDR, na.rm = TRUE),
    selected_groups = mean(selected_groups, na.rm = TRUE),
    selected_features = mean(selected_features, na.rm = TRUE),
    n_ok = dplyr::n(),
    .groups = "drop"
  )

summary_setting <- all_res %>%
  group_by(
    setting,
    dimensionality,
    signal_pattern,
    method
  ) %>%
  summarise(
    risk = mean(risk, na.rm = TRUE),
    MSE_y = mean(MSE_y, na.rm = TRUE),
    MSE_beta = mean(MSE_beta, na.rm = TRUE),
    pve = mean(pve, na.rm = TRUE),
    FDR = mean(FDR, na.rm = TRUE),
    selected_groups = mean(selected_groups, na.rm = TRUE),
    selected_features = mean(selected_features, na.rm = TRUE),
    n_ok = dplyr::n(),
    .groups = "drop"
  )

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(
  summary_snr_rhob,
  file.path(tables_dir, "plot_summary_by_snr_rhob.csv")
)

readr::write_csv(
  summary_setting,
  file.path(tables_dir, "plot_summary_by_setting.csv")
)

############################################################
# 6. Plot design choices
############################################################

# Representative SNR levels used in the manuscript figures.
# If exact values are not available, the closest values are used.
snr_targets <- c(0.05, 1.528, 6)

snr_keep <- sort(unique(summary_snr_rhob$snr))[
  sapply(
    snr_targets,
    function(x) which.min(abs(sort(unique(summary_snr_rhob$snr)) - x))
  )
]

snr_keep <- unique(snr_keep)

log_msg(
  "Using representative SNR values: %s",
  paste(round(snr_keep, 3), collapse = ", ")
)

# Representative rho values for SNR-trend figures.
rho_targets <- c(0.2, 0.5, 0.8)

rho_keep <- sort(unique(summary_snr_rhob$rho_b))[
  sapply(
    rho_targets,
    function(x) which.min(abs(sort(unique(summary_snr_rhob$rho_b)) - x))
  )
]

rho_keep <- unique(rho_keep)

log_msg(
  "Using representative rho_b values: %s",
  paste(round(rho_keep, 1), collapse = ", ")
)

############################################################
# 7. High-dimensional figures
############################################################

hd_data <- summary_snr_rhob %>%
  filter(
    dimensionality == "HD",
    snr %in% snr_keep
  )

fig_hd_risk <- ggplot(
  hd_data,
  aes(
    x = rho_b,
    y = FDR,
    color = method,
    linetype = method,
    group = method
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_wrap(setting ~ snr_label, scales = "free") +
  scale_x_continuous(
    breaks = seq(0, 0.9, by = 0.3),
    limits = c(0, 0.9)
  ) +
  labs(
    title = "High-dimensional simulation: generalization risk",
    x = expression(rho[b]),
    y = "Generalization risk"
  ) +
  make_common_theme()
fig_hd_risk

save_plot_both(
  fig_hd_risk,
  "plot_HD_gen_risk",
  width = 10,
  height = 6
)

fig_hd_pve <- ggplot(
  hd_data,
  aes(
    x = rho_b,
    y = pve,
    color = method,
    linetype = method,
    group = method
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_grid(setting ~ snr_label) +
  scale_x_continuous(
    breaks = seq(0, 0.9, by = 0.3),
    limits = c(0, 0.9)
  ) +
  labs(
    title = "High-dimensional simulation: PVE",
    x = expression(rho[b]),
    y = "PVE"
  ) +
  make_common_theme()

save_plot_both(
  fig_hd_pve,
  "plot_HD_pve",
  width = 10,
  height = 6
)

fig_hd_fdr <- ggplot(
  hd_data,
  aes(
    x = rho_b,
    y = FDR,
    color = method,
    linetype = method,
    group = method
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_grid(setting ~ snr_label) +
  scale_x_continuous(
    breaks = seq(0, 0.9, by = 0.3),
    limits = c(0, 0.9)
  ) +
  scale_y_continuous(limits = c(0, 1), oob = scales::squish) +
  labs(
    title = "High-dimensional simulation: false discovery rate",
    x = expression(rho[b]),
    y = "FDR"
  ) +
  make_common_theme()

save_plot_both(
  fig_hd_fdr,
  "plot_HD_fdr",
  width = 10,
  height = 6
)

# Selected groups vs SNR for representative rho_b values.
hd_snr_data <- summary_snr_rhob %>%
  filter(
    dimensionality == "HD",
    rho_b %in% rho_keep
  )

fig_hd_selected_groups <- ggplot(
  hd_snr_data,
  aes(
    x = snr,
    y = selected_groups,
    color = method,
    linetype = method,
    group = method
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_grid(setting ~ paste0("rho_b = ", rho_b), scales = "free_y") +
  scale_x_log10(
    breaks = sort(unique(hd_snr_data$snr)),
    labels = round(sort(unique(hd_snr_data$snr)), 3)
  ) +
  labs(
    title = "High-dimensional simulation: selected groups",
    x = "SNR",
    y = expression(hat(s))
  ) +
  make_common_theme()

save_plot_both(
  fig_hd_selected_groups,
  "plot_HD_selected_groups",
  width = 10,
  height = 6
)

############################################################
# 8. Low-dimensional figures
############################################################

ld_data <- summary_snr_rhob %>%
  filter(
    dimensionality == "LD",
    snr %in% snr_keep
  )

fig_ld_risk <- ggplot(
  ld_data,
  aes(
    x = rho_b,
    y = risk,
    color = method,
    linetype = method,
    group = method
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_grid(setting ~ snr_label, scales = "free_y") +
  scale_x_continuous(
    breaks = seq(0, 0.9, by = 0.3),
    limits = c(0, 0.9)
  ) +
  labs(
    title = "Low-dimensional simulation: generalization risk",
    x = expression(rho[b]),
    y = "Generalization risk"
  ) +
  make_common_theme()

save_plot_both(
  fig_ld_risk,
  "plot_LD_gen_risk",
  width = 10,
  height = 7
)

fig_ld_pve <- ggplot(
  ld_data,
  aes(
    x = rho_b,
    y = pve,
    color = method,
    linetype = method,
    group = method
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_grid(setting ~ snr_label) +
  scale_x_continuous(
    breaks = seq(0, 0.9, by = 0.3),
    limits = c(0, 0.9)
  ) +
  labs(
    title = "Low-dimensional simulation: PVE",
    x = expression(rho[b]),
    y = "PVE"
  ) +
  make_common_theme()

save_plot_both(
  fig_ld_pve,
  "plot_LD_pve",
  width = 10,
  height = 7
)

fig_ld_fdr <- ggplot(
  ld_data,
  aes(
    x = rho_b,
    y = FDR,
    color = method,
    linetype = method,
    group = method
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_grid(setting ~ snr_label) +
  scale_x_continuous(
    breaks = seq(0, 0.9, by = 0.3),
    limits = c(0, 0.9)
  ) +
  scale_y_continuous(limits = c(0, 1), oob = scales::squish) +
  labs(
    title = "Low-dimensional simulation: false discovery rate",
    x = expression(rho[b]),
    y = "FDR"
  ) +
  make_common_theme()

save_plot_both(
  fig_ld_fdr,
  "plot_LD_fdr",
  width = 10,
  height = 7
)

ld_snr_data <- summary_snr_rhob %>%
  filter(
    dimensionality == "LD",
    rho_b %in% rho_keep
  )

fig_ld_selected_groups <- ggplot(
  ld_snr_data,
  aes(
    x = snr,
    y = selected_groups,
    color = method,
    linetype = method,
    group = method
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_grid(setting ~ paste0("rho_b = ", rho_b), scales = "free_y") +
  scale_x_log10(
    breaks = sort(unique(ld_snr_data$snr)),
    labels = round(sort(unique(ld_snr_data$snr)), 3)
  ) +
  labs(
    title = "Low-dimensional simulation: selected groups",
    x = "SNR",
    y = expression(hat(s))
  ) +
  make_common_theme()

save_plot_both(
  fig_ld_selected_groups,
  "plot_LD_selected_groups",
  width = 10,
  height = 7
)

############################################################
# 9. Combined 3-panel figures for manuscript
############################################################

if (requireNamespace("patchwork", quietly = TRUE)) {
  
  library(patchwork)
  
  fig_hd_three <- fig_hd_risk / fig_hd_pve / fig_hd_fdr +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  save_plot_both(
    fig_hd_three,
    "plot_HD_risk_pve_fdr_combined",
    width = 11,
    height = 13
  )
  
  fig_ld_three <- fig_ld_risk / fig_ld_pve / fig_ld_fdr +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  save_plot_both(
    fig_ld_three,
    "plot_LD_risk_pve_fdr_combined",
    width = 11,
    height = 14
  )
  
} else {
  log_msg(
    "Package 'patchwork' is not installed. Combined figures were skipped."
  )
}

############################################################
# 10. Optional: tuning sensitivity plots if available
############################################################

tune_file <- "results/tuning_sensitivity/sglasso_tuning_all_raw.rds"

if (file.exists(tune_file)) {
  
  log_msg("Tuning sensitivity file found: %s", tune_file)
  
  tune_res <- readRDS(tune_file)
  
  tune_summary <- tune_res %>%
    filter(converged == TRUE) %>%
    mutate(
      tuning_criterion = dplyr::recode(
        tuning_criterion,
        "Min_val" = "Validation",
        .default = tuning_criterion
      ),
      tuning_criterion = factor(
        tuning_criterion,
        levels = c("Validation", "AIC", "BIC", "EBIC", "GCV", "AICc")
      )
    ) %>%
    group_by(
      setting,
      rho_b,
      snr,
      signal_pattern,
      tuning_criterion
    ) %>%
    summarise(
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
  
  readr::write_csv(
    tune_summary,
    file.path(tables_dir, "tuning_sensitivity_plot_summary.csv")
  )
  
  fig_tune_msey <- ggplot(
    tune_summary,
    aes(
      x = tuning_criterion,
      y = MSE_y,
      fill = tuning_criterion
    )
  ) +
    geom_col(show.legend = FALSE) +
    facet_grid(setting ~ paste0("rho_b = ", rho_b), scales = "free_y") +
    labs(
      title = "Sensitivity to tuning criteria",
      x = NULL,
      y = expression(MSE[y])
    ) +
    make_common_theme() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  save_plot_both(
    fig_tune_msey,
    "plot_tuning_sensitivity_MSE_y",
    width = 9,
    height = 6
  )
  
  fig_tune_model_size <- ggplot(
    tune_summary,
    aes(
      x = tuning_criterion,
      y = selected_groups,
      fill = tuning_criterion
    )
  ) +
    geom_col(show.legend = FALSE) +
    facet_grid(setting ~ paste0("rho_b = ", rho_b), scales = "free_y") +
    labs(
      title = "Tuning criteria and selected model size",
      x = NULL,
      y = expression(hat(s))
    ) +
    make_common_theme() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  save_plot_both(
    fig_tune_model_size,
    "plot_tuning_sensitivity_selected_groups",
    width = 9,
    height = 6
  )
  
} else {
  log_msg("No tuning sensitivity file found. Skipping tuning plots.")
}

############################################################
# 11. Done
############################################################

log_msg("All plots completed successfully.")
