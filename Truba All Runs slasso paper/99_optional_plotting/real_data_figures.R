# ============================================================
# Real-data figure: Relative MSE vs sparsity
# ============================================================

library(dplyr)
library(ggplot2)

Sim_Results <- readRDS("results/real_data/ALL_raw_results.rds")

Figure_Data <- Sim_Results |>
  group_by(dataset, rep) |>
  mutate(
    sglasso_mse = mse[method == "SGLASSO"][1],
    rel_mse = mse / sglasso_mse
  ) |>
  ungroup() |>
  group_by(dataset, method) |>
  summarise(
    rel_mse_med = median(rel_mse, na.rm = TRUE),
    sparsity_med = median(sparsity, na.rm = TRUE),
    time_med = median(time, na.rm = TRUE),
    .groups = "drop"
  )

p_tradeoff <- ggplot(
  Figure_Data,
  aes(
    x = sparsity_med,
    y = rel_mse_med,
    label = method
  )
) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(aes(size = time_med), alpha = 0.8) +
  geom_text(
    vjust = -0.7,
    size = 3
  ) +
  facet_wrap(~ dataset, scales = "free_x") +
  labs(
    x = "Median sparsity ratio",
    y = "Relative MSE",
    size = "Median time (s)",
    caption = "Relative MSE is computed with respect to SGLASSO within each dataset and replication."
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

p_tradeoff

ggsave(
  filename = "Figure_realdata_mse_sparsity_tradeoff.pdf",
  plot = p_tradeoff,
  width = 9,
  height = 4.5
)

ggsave(
  filename = "Figure_realdata_mse_sparsity_tradeoff.png",
  plot = p_tradeoff,
  width = 9,
  height = 4.5,
  dpi = 300
)


postscript(
  file = "Figure_realdata_mse_sparsity_tradeoff.eps",
  width = 9,
  height = 4.5,
  horizontal = FALSE,
  onefile = FALSE,
  paper = "special"
)

print(p_tradeoff)

dev.off()
