# ============================================================
# Real-data benchmark summary tables
# Creates two publication-ready LaTeX tables:
#
# Table 3: Prediction and computational performance
# Table 4: Model sparsity and selection complexity
# ============================================================

library(dplyr)
library(tibble)
library(kableExtra)

# ------------------------------------------------------------
# 1. Read raw benchmark results
# ------------------------------------------------------------

Sim_Results <- readRDS("results/real_data/ALL_raw_results.rds")

# Choose summary type:
summary_type <- "median"
# summary_type <- "mean"

if (!summary_type %in% c("mean", "median")) {
  stop("summary_type must be either 'mean' or 'median'.")
}

main_fun <- switch(
  summary_type,
  mean = mean,
  median = median
)

spread_fun <- switch(
  summary_type,
  mean = function(x) stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))),
  median = function(x) stats::IQR(x, na.rm = TRUE)
)

spread_label <- switch(
  summary_type,
  mean = "standard error",
  median = "IQR"
)

# ------------------------------------------------------------
# 2. Dataset information
# ------------------------------------------------------------

Dataset_Info <- tibble(
  dataset = c("Birthwt", "GenAtHum", "bardet"),
  n = c(189, 158, 120),
  p = c(16, 2045, 100),
  groups_total = c(8, 79, 20)
)

# ------------------------------------------------------------
# 3. Compute benchmark summaries
# ------------------------------------------------------------

Benchmark_Summary <- Sim_Results |>
  mutate(
    lambda_clean = if_else(
      is.finite(lambda) & lambda > 0,
      lambda,
      NA_real_
    ),
    rmse_raw = sqrt(mse)
  ) |>
  group_by(dataset, rep) |>
  mutate(
    sglasso_rmse = rmse_raw[method == "SGLASSO"][1],
    rmse = rmse_raw / sglasso_rmse
  ) |>
  ungroup() |>
  group_by(dataset, method) |>
  summarise(
    stat_rmse = main_fun(rmse, na.rm = TRUE),
    spread_rmse = spread_fun(rmse),
    
    stat_mse = main_fun(mse, na.rm = TRUE),
    spread_mse = spread_fun(mse),
    
    stat_time = main_fun(time, na.rm = TRUE),
    spread_time = spread_fun(time),
    
    stat_iter = main_fun(iter, na.rm = TRUE),
    spread_iter = spread_fun(iter),
    
    stat_groups = main_fun(selected_groups, na.rm = TRUE),
    spread_groups = spread_fun(selected_groups),
    
    stat_features = main_fun(selected_features, na.rm = TRUE),
    spread_features = spread_fun(selected_features),
    
    stat_sparsity = main_fun(sparsity, na.rm = TRUE),
    spread_sparsity = spread_fun(sparsity),
    
    stat_alpha = main_fun(alpha, na.rm = TRUE),
    
    stat_lambda = main_fun(lambda_clean, na.rm = TRUE),
    spread_lambda = spread_fun(lambda_clean),
    
    stat_d = main_fun(d, na.rm = TRUE),
    
    nrep_success = n_distinct(rep),
    
    .groups = "drop"
  )

# ============================================================
# 4. Table 3: Prediction and computational performance
# ============================================================

Prediction_Table <- Benchmark_Summary |>
  mutate(
    RMSE = sprintf("%.4f", stat_rmse),
    
    MSE = sprintf("%.3f (%.3f)", stat_mse, spread_mse),
    
    `Time (s)` = sprintf("%.3f (%.3f)", stat_time, spread_time),
    
    Iter = if_else(
      is.na(stat_iter) | is.nan(stat_iter),
      "-",
      sprintf("%.1f (%.1f)", stat_iter, spread_iter)
    ),
    
    Alpha = if_else(
      is.na(stat_alpha) | is.nan(stat_alpha),
      "-",
      sprintf("%.2f", stat_alpha)
    ),
    
    Lambda = if_else(
      is.na(stat_lambda) | is.nan(stat_lambda),
      "-",
      sprintf("%.3f (%.3f)", stat_lambda, spread_lambda)
    ),
    
    d = if_else(
      is.na(stat_d) | is.nan(stat_d),
      "-",
      sprintf("%.2f", stat_d)
    )
  ) |>
  select(
    Dataset = dataset,
    Method = method,
    RMSE,
    MSE,
    `Time (s)`,
    Iter,
    Alpha,
    Lambda,
    d
  ) |>
  arrange(
    Dataset,
    desc(Method == "SGLASSO"),
    RMSE
  )

write.csv(
  Prediction_Table,
  "Prediction_Table_By_Dataset.csv",
  row.names = FALSE
)

Latex_Prediction_Table <- Prediction_Table

colnames(Latex_Prediction_Table) <- c(
  "Dataset",
  "Method",
  "RMSE",
  "MSE",
  "Time (s)",
  "Iter",
  "$\\alpha$",
  "$\\lambda$",
  "$d$"
)

latex_prediction <- Latex_Prediction_Table |>
  kbl(
    format = "latex",
    booktabs = TRUE,
    longtable = FALSE,
    align = "llccccccc",
    caption = "Prediction and computational benchmark comparison across real-data examples.",
    label = "tab:real_prediction_benchmark",
    escape = FALSE,
    linesep = ""
  ) |>
  kable_styling(
    latex_options = "hold_position",
    font_size = 9,
    full_width = FALSE,
    position = "center"
  ) |>
  collapse_rows(
    columns = 1,
    valign = "top"
  ) |>
  row_spec(0, bold = TRUE) |>
  footnote(
    general = paste(
      "Values are reported as", summary_type,
      "(", spread_label, ").",
      "RMSE is computed relative to SGLASSO within each dataset and replication.",
      "MSE and runtime are on their original scales.",
      "Runtime refers to the final non-cross-validated refit after tuning.",
      "$d$ is reported only for SGLASSO."
    ),
    general_title = "Note:",
    footnote_as_chunk = FALSE,
    threeparttable = TRUE,
    escape = FALSE
  )

cat(latex_prediction, file = "Prediction_Table_By_Dataset.tex")

# ============================================================
# 5. Table 4: Sparsity and model complexity
# ============================================================

Sparsity_Table <- Benchmark_Summary |>
  left_join(
    Dataset_Info,
    by = "dataset"
  ) |>
  mutate(
    `Sel. Groups` = if_else(
      is.na(stat_groups) | is.nan(stat_groups),
      "-",
      sprintf("%.1f (%.1f)", stat_groups, spread_groups)
    ),
    
    `Sel. Features` = if_else(
      is.na(stat_features) | is.nan(stat_features),
      "-",
      sprintf("%.1f (%.1f)", stat_features, spread_features)
    ),
    
    Sparsity = if_else(
      is.na(stat_sparsity) | is.nan(stat_sparsity),
      "-",
      sprintf("%.3f (%.3f)", stat_sparsity, spread_sparsity)
    )
  ) |>
  select(
    Dataset = dataset,
    n,
    p,
    Groups = groups_total,
    Method = method,
    `Sel. Groups`,
    `Sel. Features`,
    Sparsity
  ) |>
  arrange(
    Dataset,
    desc(Method == "SGLASSO")
  )

write.csv(
  Sparsity_Table,
  "Sparsity_Table_By_Dataset.csv",
  row.names = FALSE
)

Latex_Sparsity_Table <- Sparsity_Table

colnames(Latex_Sparsity_Table) <- c(
  "Dataset",
  "$n$",
  "$p$",
  "Groups",
  "Method",
  "Sel. Groups",
  "Sel. Features",
  "Sparsity"
)

latex_sparsity <- Latex_Sparsity_Table |>
  kbl(
    format = "latex",
    booktabs = TRUE,
    longtable = FALSE,
    align = "lccclccc",
    caption = "Model sparsity and selection complexity across real-data examples.",
    label = "tab:real_sparsity_benchmark",
    escape = FALSE,
    linesep = ""
  ) |>
  kable_styling(
    latex_options = "hold_position",
    font_size = 9,
    full_width = FALSE,
    position = "center"
  ) |>
  collapse_rows(
    columns = 1:4,
    valign = "top"
  ) |>
  row_spec(0, bold = TRUE) |>
  footnote(
    general = paste(
      "Values are reported as", summary_type,
      "(", spread_label, ").",
      "Selected groups and selected features are computed from the final refitted model.",
      "Sparsity is defined as one minus the proportion of selected features."
    ),
    general_title = "Note:",
    footnote_as_chunk = FALSE,
    threeparttable = TRUE,
    escape = FALSE
  )

cat(latex_sparsity, file = "Sparsity_Table_By_Dataset.tex")

# ------------------------------------------------------------
# 6. Print tables in R console
# ------------------------------------------------------------

Prediction_Table
Sparsity_Table

cat("Saved files:\n")
cat("- Prediction_Table_By_Dataset.csv\n")
cat("- Prediction_Table_By_Dataset.tex\n")
cat("- Sparsity_Table_By_Dataset.csv\n")
cat("- Sparsity_Table_By_Dataset.tex\n")
