############################################################
# Real-data figures for the SGLASSO paper.
############################################################

options(stringsAsFactors = FALSE)

library(dplyr)
library(ggplot2)

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  key <- paste0("--", name, "=")
  hit <- grep(paste0("^", key), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(key, "", hit[[1]])
}

base_dir <- get_arg("input-dir", Sys.getenv("REAL_DATA_DIR", "results/real_data"))
fig_dir <- get_arg("fig-dir", file.path(base_dir, "figures"))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

raw_file <- file.path(base_dir, "ALL_raw_results.rds")
summary_file <- file.path(base_dir, "ALL_summary.rds")

if (!file.exists(raw_file)) {
  stop("Cannot find raw results file: ", raw_file, call. = FALSE)
}
if (!file.exists(summary_file)) {
  stop("Cannot find summary results file: ", summary_file, call. = FALSE)
}

raw_results <- readRDS(raw_file)
summary_results <- readRDS(summary_file)

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

method_shapes <- c(
  ADELIE  = 21,
  GLASSO  = 22,
  GENET   = 23,
  SGLASSO = 24,
  GSCAD   = 25,
  GMCP    = 21
)

theme_real <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      strip.background = element_rect(fill = "grey82", color = "grey35",
                                      linewidth = 0.55),
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "grey35", fill = NA, linewidth = 0.55),
      panel.grid.major = element_line(color = "grey86", linewidth = 0.45),
      panel.grid.minor = element_line(color = "grey94", linewidth = 0.22),
      axis.text = element_text(face = "bold", color = "grey30"),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(face = "bold"),
      legend.key.width = unit(1.25, "cm"),
      plot.margin = margin(8, 12, 8, 8)
    )
}

save_plot <- function(plot, stem, width = 9, height = 4.8) {
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot,
         width = width, height = height)
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot,
         width = width, height = height, dpi = 300)
  ggsave(file.path(fig_dir, paste0(stem, ".eps")), plot,
         width = width, height = height, device = cairo_ps)
}

summary_plot_data <- summary_results %>%
  filter(method %in% method_levels) %>%
  mutate(
    method = factor(method, levels = method_levels),
    dataset = factor(dataset, levels = unique(dataset))
  )

write.csv(
  summary_plot_data,
  file.path(fig_dir, "real_data_summary_for_figures.csv"),
  row.names = FALSE
)

sglasso_ref <- raw_results %>%
  filter(method == "SGLASSO") %>%
  select(dataset, rep, sglasso_mse = mse)

tradeoff_data <- raw_results %>%
  filter(method %in% method_levels) %>%
  left_join(sglasso_ref, by = c("dataset", "rep")) %>%
  mutate(
    rel_mse = mse / sglasso_mse,
    method = factor(method, levels = method_levels),
    dataset = factor(dataset, levels = unique(dataset))
  ) %>%
  group_by(dataset, method) %>%
  summarise(
    rel_mse_median = median(rel_mse, na.rm = TRUE),
    rel_mse_mean = mean(rel_mse, na.rm = TRUE),
    sparsity_median = median(sparsity, na.rm = TRUE),
    selected_groups_median = median(selected_groups, na.rm = TRUE),
    selected_features_median = median(selected_features, na.rm = TRUE),
    time_median = median(time, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  tradeoff_data,
  file.path(fig_dir, "real_data_tradeoff_for_figures.csv"),
  row.names = FALSE
)

p_tradeoff <- ggplot(
  tradeoff_data,
  aes(
    x = sparsity_median,
    y = rel_mse_median,
    color = method,
    shape = method,
    size = time_median
  )
) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey45") +
  geom_point(fill = "white", stroke = 1.0, alpha = 0.9) +
  geom_text(
    aes(label = method),
    vjust = -0.75,
    size = 3,
    color = "grey20",
    show.legend = FALSE
  ) +
  facet_wrap(~dataset, scales = "free_x", nrow = 1) +
  scale_color_manual(values = method_colors, drop = FALSE) +
  scale_shape_manual(values = method_shapes, drop = FALSE) +
  labs(
    x = "Median sparsity ratio",
    y = "Median relative MSE",
    color = NULL,
    shape = NULL,
    size = "Median time (s)"
  ) +
  theme_real() +
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE),
    shape = guide_legend(nrow = 1, byrow = TRUE)
  )

save_plot(p_tradeoff, "real_data_mse_sparsity_tradeoff", width = 10, height = 4.8)

make_metric_plot <- function(data, metric, y_label, stem, log_y = FALSE) {
  p <- ggplot(
    data,
    aes(
      x = method,
      y = .data[[metric]],
      color = method,
      shape = method
    )
  ) +
    geom_point(size = 2.8, fill = "white", stroke = 1.0) +
    facet_wrap(~dataset, scales = "free_y", nrow = 1) +
    scale_color_manual(values = method_colors, drop = FALSE) +
    scale_shape_manual(values = method_shapes, drop = FALSE) +
    labs(
      x = NULL,
      y = y_label,
      color = NULL,
      shape = NULL
    ) +
    theme_real() +
    guides(
      color = guide_legend(nrow = 1, byrow = TRUE),
      shape = guide_legend(nrow = 1, byrow = TRUE)
    )

  if (log_y) {
    p <- p + scale_y_log10(labels = function(x) format(x, trim = TRUE))
  }

  save_plot(p, stem, width = 10, height = 4.8)
  invisible(p)
}

p_mse <- make_metric_plot(
  summary_plot_data,
  "mse_median",
  "Median test MSE",
  "real_data_mse_median"
)

p_mae <- make_metric_plot(
  summary_plot_data,
  "mae_median",
  "Median test MAE",
  "real_data_mae_median"
)

p_selected_groups <- make_metric_plot(
  summary_plot_data,
  "selected_groups_median",
  "Median selected groups",
  "real_data_selected_groups_median"
)

p_selected_features <- make_metric_plot(
  summary_plot_data,
  "selected_features_median",
  "Median selected predictors",
  "real_data_selected_features_median"
)

p_time <- make_metric_plot(
  summary_plot_data,
  "time_median",
  "Median runtime (seconds)",
  "real_data_runtime_median",
  log_y = TRUE
)

message("Real-data figures saved in: ", normalizePath(fig_dir, mustWork = FALSE))
