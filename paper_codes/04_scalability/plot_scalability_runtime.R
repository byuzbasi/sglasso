############################################################
# Plot scalability runtime curves.
############################################################

options(stringsAsFactors = FALSE)

library(readr)
library(dplyr)
library(ggplot2)

base_dir <- "results/scalability_low_to_high"
fig_dir <- file.path(base_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

summary_file <- file.path(base_dir, "scalability_summary.csv")
summary_scalability <- readr::read_csv(summary_file, show_col_types = FALSE)

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

format_big_number <- function(x) {
  ifelse(
    x >= 1000,
    paste0(format(x / 1000, trim = TRUE, scientific = FALSE), "k"),
    format(x, trim = TRUE, scientific = FALSE)
  )
}

plot_dat <- summary_scalability %>%
  mutate(
    method = factor(method, levels = method_levels),
    axis_label = ifelse(axis == "p", "Runtime vs p", "Runtime vs J"),
    x_value = ifelse(axis == "p", p, J)
  ) %>%
  arrange(axis, method, x_value)

axis_ranges <- plot_dat %>%
  group_by(axis_label) %>%
  summarise(
    x_min = min(x_value),
    x_max = max(x_value),
    .groups = "drop"
  )

major_x_breaks <- sort(unique(c(
  axis_ranges$x_min,
  axis_ranges$x_max,
  3000, 6000, 10000, 15000, 30000, 50000,
  1000, 3000, 6000, 10000, 15000
)))

make_runtime_plot <- function(data,
                              include_adelie = TRUE,
                              log_y = TRUE) {
  if (!include_adelie) {
    data <- dplyr::filter(data, method != "ADELIE")
  }

  y_label <- if (log_y) {
    "Median runtime (seconds, log scale)"
  } else {
    "Median runtime (seconds)"
  }

  p <- ggplot(
    data,
    aes(
      x = x_value,
      y = median_runtime,
      color = method,
      linetype = method,
      linewidth = method,
      group = method
    )
  ) +
    geom_vline(
      data = axis_ranges,
      aes(xintercept = x_min),
      inherit.aes = FALSE,
      color = "grey72",
      linewidth = 0.45
    ) +
    geom_vline(
      data = axis_ranges,
      aes(xintercept = x_max),
      inherit.aes = FALSE,
      color = "grey72",
      linewidth = 0.45
    ) +
    geom_line(lineend = "round") +
    facet_wrap(~axis_label, scales = "free_x", nrow = 1) +
    scale_color_manual(values = method_colors, drop = TRUE) +
    scale_linetype_manual(values = method_linetypes, drop = TRUE) +
    scale_linewidth_manual(values = method_linewidths, drop = TRUE) +
    scale_x_continuous(
      breaks = major_x_breaks,
      labels = format_big_number,
      minor_breaks = NULL,
      expand = expansion(mult = c(0.025, 0.025))
    ) +
    labs(
      x = NULL,
      y = y_label,
      color = NULL,
      linetype = NULL,
      linewidth = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(
      strip.background = element_rect(fill = "grey82", color = "grey35",
                                      linewidth = 0.55),
      strip.text = element_text(face = "bold", size = 14),
      panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.55),
      panel.grid.major.x = element_line(color = "grey84", linewidth = 0.48),
      panel.grid.major.y = element_line(color = "grey86", linewidth = 0.48),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_line(color = "grey94", linewidth = 0.22),
      axis.text = element_text(face = "bold", color = "grey30", size = 10),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, margin = margin(t = 3)),
      axis.title.y = element_text(face = "bold", size = 15),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = element_text(face = "bold", size = 12),
      legend.key.width = unit(1.5, "cm"),
      plot.margin = margin(8, 12, 8, 8)
    ) +
    guides(
      color = guide_legend(nrow = 1, byrow = TRUE),
      linetype = guide_legend(nrow = 1, byrow = TRUE),
      linewidth = "none"
    )

  if (log_y) {
    p <- p +
      scale_y_log10(
        breaks = c(0.01, 0.03, 0.1, 0.3, 1, 3, 10),
        labels = function(x) format(x, trim = TRUE, scientific = FALSE)
      )
  }

  p
}

plot_all_log <- make_runtime_plot(
  plot_dat,
  include_adelie = TRUE,
  log_y = TRUE
)

plot_all_linear <- make_runtime_plot(
  plot_dat,
  include_adelie = TRUE,
  log_y = FALSE
)

plot_no_adelie_log <- make_runtime_plot(
  plot_dat,
  include_adelie = FALSE,
  log_y = TRUE
)

plot_no_adelie_linear <- make_runtime_plot(
  plot_dat,
  include_adelie = FALSE,
  log_y = FALSE
)

ggsave(
  file.path(fig_dir, "runtime_scalability_curve_all_methods_log.pdf"),
  plot_all_log,
  width = 10.5,
  height = 4.8
)
ggsave(
  file.path(fig_dir, "runtime_scalability_curve_all_methods_log.eps"),
  plot_all_log,
  width = 10.5,
  height = 4.8,
  device = cairo_ps
)
ggsave(
  file.path(fig_dir, "runtime_scalability_curve_all_methods_log.png"),
  plot_all_log,
  width = 10.5,
  height = 4.8,
  dpi = 300
)

ggsave(
  file.path(fig_dir, "runtime_scalability_curve_all_methods_linear.pdf"),
  plot_all_linear,
  width = 10.5,
  height = 4.8
)
ggsave(
  file.path(fig_dir, "runtime_scalability_curve_all_methods_linear.eps"),
  plot_all_linear,
  width = 10.5,
  height = 4.8,
  device = cairo_ps
)
ggsave(
  file.path(fig_dir, "runtime_scalability_curve_all_methods_linear.png"),
  plot_all_linear,
  width = 10.5,
  height = 4.8,
  dpi = 300
)

ggsave(
  file.path(fig_dir, "runtime_scalability_curve_no_adelie_log.pdf"),
  plot_no_adelie_log,
  width = 10.5,
  height = 4.8
)
ggsave(
  file.path(fig_dir, "runtime_scalability_curve_no_adelie_log.eps"),
  plot_no_adelie_log,
  width = 10.5,
  height = 4.8,
  device = cairo_ps
)
ggsave(
  file.path(fig_dir, "runtime_scalability_curve_no_adelie_log.png"),
  plot_no_adelie_log,
  width = 10.5,
  height = 4.8,
  dpi = 300
)

ggsave(
  file.path(fig_dir, "runtime_scalability_curve_no_adelie_linear.pdf"),
  plot_no_adelie_linear,
  width = 10.5,
  height = 4.8
)
ggsave(
  file.path(fig_dir, "runtime_scalability_curve_no_adelie_linear.eps"),
  plot_no_adelie_linear,
  width = 10.5,
  height = 4.8,
  device = cairo_ps
)
ggsave(
  file.path(fig_dir, "runtime_scalability_curve_no_adelie_linear.png"),
  plot_no_adelie_linear,
  width = 10.5,
  height = 4.8,
  dpi = 300
)

largest_summary <- summary_scalability %>%
  group_by(axis) %>%
  filter(if (unique(axis) == "p") p == max(p) else J == max(J)) %>%
  ungroup() %>%
  arrange(axis, median_runtime)

readr::write_csv(
  largest_summary,
  file.path(base_dir, "largest_scenario_runtime_summary.csv")
)

print(plot_all_log)
message("Saved scalability figures in: ", normalizePath(fig_dir))
