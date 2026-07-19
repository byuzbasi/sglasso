############################################################
# Paper plots for the main simulation
############################################################
#
# This script reads the summary tables created by
# run_sim_slasso.all.R and produces manuscript-ready figures.
#
# Default input:
#   results/main_simulation_rep50_fixed/tables/paper_summary_snr_rhob.csv
#
# Default output:
#   results/main_simulation_rep50_fixed/figures
#
# Usage from R:
#   source("paper_codes/02_main_simulation/plots_paper.R")
#
# Usage from terminal:
#   Rscript --vanilla paper_codes/02_main_simulation/plots_paper.R \
#     --results-dir=results/main_simulation_rep50_fixed

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
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

split_arg <- function(x) {
  x <- trimws(x)
  if (!nzchar(x)) {
    return(character())
  }
  trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
}

results_dir <- get_arg(
  "results-dir",
  Sys.getenv("MAIN_SIM_OUTDIR", unset = "results/main_simulation_rep50_fixed")
)
tables_dir <- file.path(results_dir, "tables")
fig_dir <- get_arg("fig-dir", file.path(results_dir, "figures"))
summary_file <- file.path(tables_dir, "paper_summary_snr_rhob.csv")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(summary_file)) {
  stop("Cannot find summary file: ", summary_file, call. = FALSE)
}

summary_dat <- readr::read_csv(summary_file, show_col_types = FALSE)

############################################################
# Method styling
############################################################

# Keep these identical to paper_codes/04_scalability/plot_scalability_runtime.R.
method_levels <- c("ADELIE", "GLASSO", "GENET", "SGLASSO", "GSCAD", "GMCP")

method_colors <- c(
  ADELIE  = "#4D4D4D",
  GLASSO  = "#2F8F46",
  GENET   = "#4AA6B8",
  SGLASSO = "#1F5AB8",
  GSCAD   = "#FF5733",
  GMCP    = "#FF8C28"
)

method_linetypes <- c(
  ADELIE  = "solid",
  GLASSO  = "longdash",
  GENET   = "dotted",
  SGLASSO = "solid",
  GSCAD   = "longdash",
  GMCP    = "dotdash"
)

method_linewidths <- c(
  ADELIE  = 0.95,
  GLASSO  = 1.15,
  GENET   = 1.15,
  SGLASSO = 1.30,
  GSCAD   = 1.15,
  GMCP    = 1.15
)

############################################################
# Plot settings
############################################################

selected_methods <- split_arg(get_arg(
  "methods",
  paste(method_levels, collapse = ",")
))

selected_settings <- split_arg(get_arg(
  "settings",
  paste(c("LD-4", "LD-8", "LD-12", "HD-10", "HD-20", "HD-40"), collapse = ",")
))

selected_snr <- as.numeric(split_arg(get_arg("snr-values", "1.528")))
selected_rho_for_snr_plots <- as.numeric(split_arg(get_arg("rho-values", "0.5")))

selected_patterns <- split_arg(get_arg(
  "signal-patterns",
  "homogeneous,weak_strong_mixed"
))

format_snr_label <- function(x) {
  format(round(as.numeric(x), 3), trim = TRUE, scientific = FALSE)
}

clean_pattern_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

format_pattern_label <- function(x) {
  dplyr::recode(
    x,
    homogeneous = "Homogeneous",
    weak_strong_mixed = "Weak/strong mixed",
    mixed_signs = "Mixed signs",
    .default = x
  )
}

plot_dat <- summary_dat %>%
  filter(
    method %in% selected_methods,
    setting %in% selected_settings,
    signal_pattern %in% selected_patterns
  ) %>%
  mutate(
    method = factor(method, levels = method_levels),
    setting = factor(setting, levels = selected_settings),
    signal_pattern_label = factor(
      format_pattern_label(signal_pattern),
      levels = format_pattern_label(selected_patterns)
    ),
    snr_label = factor(format_snr_label(snr),
                       levels = format_snr_label(sort(unique(snr)))),
    rho_b_label = factor(format(rho_b, trim = TRUE, scientific = FALSE),
                         levels = format(sort(unique(rho_b)), trim = TRUE,
                                         scientific = FALSE))
  ) %>%
  arrange(signal_pattern, setting, snr, rho_b, method)

if (nrow(plot_dat) == 0L) {
  stop("No rows remain after filtering methods/settings/signal patterns.", call. = FALSE)
}

theme_paper <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      strip.background = element_rect(fill = "grey82", color = "grey35",
                                      linewidth = 0.55),
      strip.text = element_text(face = "bold", size = base_size),
      panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.55),
      panel.grid.major.x = element_line(color = "grey84", linewidth = 0.48),
      panel.grid.major.y = element_line(color = "grey86", linewidth = 0.48),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_line(color = "grey94", linewidth = 0.22),
      axis.text = element_text(face = "bold", color = "grey30", size = base_size - 3),
      axis.title = element_text(face = "bold", size = base_size + 1),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = element_text(face = "bold", size = base_size - 2),
      legend.key.width = unit(1.45, "cm"),
      plot.margin = margin(8, 12, 8, 8)
    )
}

add_method_scales <- function(p) {
  p +
    scale_color_manual(values = method_colors, drop = TRUE) +
    scale_linetype_manual(values = method_linetypes, drop = TRUE) +
    scale_linewidth_manual(values = method_linewidths, drop = TRUE) +
    guides(
      color = guide_legend(nrow = 1, byrow = TRUE),
      linetype = guide_legend(nrow = 1, byrow = TRUE),
      linewidth = "none"
    )
}

save_plot_all_formats <- function(plot, basename, width = 16, height = 10) {
  pdf_file <- file.path(fig_dir, paste0(basename, ".pdf"))
  eps_file <- file.path(fig_dir, paste0(basename, ".eps"))
  png_file <- file.path(fig_dir, paste0(basename, ".png"))

  ggsave(pdf_file, plot, width = width, height = height, units = "in")
  ggsave(eps_file, plot, width = width, height = height, units = "in",
         device = "eps")
  ggsave(png_file, plot, width = width, height = height, units = "in",
         dpi = 300)

  c(pdf_file, eps_file, png_file)
}

make_rhob_plot <- function(data, metric, y_label, selected_snr_values,
                           free_y = TRUE) {
  metric_sym <- rlang::sym(metric)

  data <- data %>%
    filter(round(snr, 3) %in% round(selected_snr_values, 3)) %>%
    mutate(
      facet_label = factor(
        paste(signal_pattern_label, setting_label, sep = "\n"),
        levels = as.vector(t(outer(
          levels(signal_pattern_label),
          levels(setting_label),
          paste,
          sep = "\n"
        )))
      )
    )

  p <- ggplot(
    data,
    aes(
      x = rho_b,
      y = !!metric_sym,
      color = method,
      linetype = method,
      linewidth = method,
      group = method
    )
  ) +
    geom_line(lineend = "round") +
    facet_wrap(
      ~facet_label,
      scales = if (free_y) "free" else "fixed",
      ncol = length(levels(data$setting_label))
    ) +
    scale_x_continuous(
      breaks = seq(0, 0.9, by = 0.1),
      minor_breaks = NULL,
      expand = expansion(mult = c(0.015, 0.02))
    ) +
    labs(
      x = expression(rho[b]),
      y = y_label,
      color = NULL,
      linetype = NULL,
      linewidth = NULL
    ) +
    theme_paper()

  add_method_scales(p)
}

make_snr_plot <- function(data, metric, y_label, selected_rho_values,
                          free_y = TRUE, add_pve_reference = FALSE) {
  metric_sym <- rlang::sym(metric)

  data <- data %>%
    filter(round(rho_b, 3) %in% round(selected_rho_values, 3)) %>%
    mutate(
      rho_b_label = factor(
        paste0("rho_b = ", format(rho_b, trim = TRUE, scientific = FALSE)),
        levels = paste0(
          "rho_b = ",
          format(selected_rho_values, trim = TRUE, scientific = FALSE)
        )
      ),
      facet_label = if (length(selected_rho_values) == 1L) {
        factor(
          paste(signal_pattern_label, setting, sep = "\n"),
          levels = as.vector(t(outer(
            levels(signal_pattern_label),
            levels(setting),
            paste,
            sep = "\n"
          )))
        )
      } else {
        factor(
          paste(signal_pattern_label, setting, rho_b_label, sep = "\n"),
          levels = as.vector(t(outer(
            levels(signal_pattern_label),
            as.vector(t(outer(
              levels(setting),
              paste0(
                "rho_b = ",
                format(selected_rho_values, trim = TRUE, scientific = FALSE)
              ),
              paste,
              sep = "\n"
            ))),
            paste,
            sep = "\n"
          )))
        )
      }
    )

  p <- ggplot(
    data,
    aes(
      x = snr,
      y = !!metric_sym,
      color = method,
      linetype = method,
      linewidth = method,
      group = method
    )
  ) +
    geom_line(lineend = "round") +
    facet_wrap(
      ~facet_label,
      scales = if (free_y) "free" else "fixed",
      ncol = if (length(selected_rho_values) == 1L) {
        length(levels(data$setting))
      } else {
        length(selected_rho_values)
      }
    ) +
    scale_x_continuous(
      trans = "log10",
      breaks = sort(unique(data$snr)),
      labels = format_snr_label,
      minor_breaks = NULL,
      expand = expansion(mult = c(0.025, 0.025))
    ) +
    labs(
      x = "Signal-to-noise ratio",
      y = y_label,
      color = NULL,
      linetype = NULL,
      linewidth = NULL
    ) +
    theme_paper() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))

  if (isTRUE(add_pve_reference)) {
    ref_dat <- data %>%
      distinct(facet_label, snr) %>%
      mutate(theoretical_pve = snr / (1 + snr))

    p <- p +
      geom_line(
        data = ref_dat,
        aes(x = snr, y = theoretical_pve, group = facet_label),
        inherit.aes = FALSE,
        color = "grey20",
        linetype = "dashed",
        linewidth = 0.55
      )
  }

  add_method_scales(p)
}

############################################################
# Generate plots
############################################################

written <- character()

dimension_groups <- list(
  LD = c("LD-4", "LD-8", "LD-12"),
  HD = c("HD-10", "HD-20", "HD-40")
)

metric_specs <- list(
  list(metric = "MSE_y", label = expression(MSE[y]), name = "mse_y",
       free_y = TRUE),
  list(metric = "MSE_beta", label = expression(MSE[beta]), name = "mse_beta",
       free_y = TRUE),
  list(metric = "pve", label = "Proportion of variance explained",
       name = "pve", free_y = TRUE),
  list(metric = "FDR", label = "False discovery rate",
       name = "fdr", free_y = TRUE),
  list(metric = "selected_groups", label = "Selected groups",
       name = "selected_groups", free_y = TRUE),
  list(metric = "selected_features", label = "Selected predictors",
       name = "selected_features", free_y = TRUE)
)

for (dimension_name in names(dimension_groups)) {
  dimension_dat <- plot_dat %>%
    filter(setting %in% dimension_groups[[dimension_name]]) %>%
    mutate(
      setting = factor(setting, levels = dimension_groups[[dimension_name]]),
      setting_label = factor(
        as.character(setting),
        levels = dimension_groups[[dimension_name]]
      )
    )

  if (nrow(dimension_dat) == 0L) {
    next
  }

  for (spec in metric_specs) {
    rhob_plot <- make_rhob_plot(
      dimension_dat,
      metric = spec$metric,
      y_label = spec$label,
      selected_snr_values = selected_snr,
      free_y = spec$free_y
    )

    written <- c(
      written,
      save_plot_all_formats(
        rhob_plot,
        paste0("main_sim_", dimension_name, "_", spec$name, "_vs_rhob"),
        width = 14,
        height = 8
      )
    )
  }

  snr_pve_plot <- make_snr_plot(
    dimension_dat,
    metric = "pve",
    y_label = "Proportion of variance explained",
    selected_rho_values = selected_rho_for_snr_plots,
    free_y = TRUE,
    add_pve_reference = TRUE
  )

  written <- c(
    written,
    save_plot_all_formats(
      snr_pve_plot,
      paste0("main_sim_", dimension_name, "_pve_vs_snr"),
      width = 14,
      height = 8
    )
  )

  snr_plot <- make_snr_plot(
    dimension_dat,
    metric = "selected_groups",
    y_label = "Selected groups",
    selected_rho_values = selected_rho_for_snr_plots,
    free_y = TRUE
  )

  written <- c(
    written,
    save_plot_all_formats(
      snr_plot,
      paste0("main_sim_", dimension_name, "_selected_groups_vs_snr"),
      width = 14,
      height = 8
    )
  )
}

cat("Read summary file:\n")
cat(" - ", summary_file, "\n", sep = "")
cat("Methods included:\n")
cat(" - ", paste(selected_methods, collapse = ", "), "\n", sep = "")
cat("Signal patterns included:\n")
cat(" - ", paste(selected_patterns, collapse = ", "), "\n", sep = "")
cat("Figures written:\n")
cat(paste0(" - ", written, collapse = "\n"), "\n")
