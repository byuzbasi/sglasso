# ============================================================
# Benchmark summary tables by dataset
# Options:
#   summary_type = "mean"   -> mean (standard error)
#   summary_type = "median" -> median (IQR)
# ============================================================

library(dplyr)

Sim_Results <- readRDS("results/real_data/ALL_raw_results.rds")

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
# 1. Compute benchmark summaries separately for each dataset
# ------------------------------------------------------------

Benchmark_Summary <- Sim_Results |>
  dplyr::mutate(
    lambda_clean = dplyr::if_else(
      is.finite(lambda) & lambda > 0,
      lambda,
      NA_real_
    ),
    rmse_raw = sqrt(mse)
  ) |>
  dplyr::group_by(dataset, rep) |>
  dplyr::mutate(
    sglasso_rmse = rmse_raw[method == "SGLASSO"][1],
    rmse = rmse_raw / sglasso_rmse
  ) |>
  dplyr::ungroup() |>
  dplyr::group_by(dataset, method) |>
  dplyr::summarise(
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
    
    nrep_success = dplyr::n_distinct(rep),
    
    .groups = "drop"
  )

# ------------------------------------------------------------
# 2. Publication-ready combined table
# ------------------------------------------------------------

Benchmark_Table_By_Dataset <- Benchmark_Summary |>
  dplyr::mutate(
    RMSE = sprintf("%.4f", stat_rmse),
    
    MSE = sprintf("%.3f (%.3f)", stat_mse, spread_mse),
    
    `Time (s)` = sprintf("%.3f (%.3f)", stat_time, spread_time),
    
    Iter = dplyr::if_else(
      is.na(stat_iter) | is.nan(stat_iter),
      "-",
      sprintf("%.1f (%.1f)", stat_iter, spread_iter)
    ),
    
    `Sel. Groups` = dplyr::if_else(
      is.na(stat_groups) | is.nan(stat_groups),
      "-",
      sprintf("%.1f (%.1f)", stat_groups, spread_groups)
    ),
    
    `Sel. Features` = dplyr::if_else(
      is.na(stat_features) | is.nan(stat_features),
      "-",
      sprintf("%.1f (%.1f)", stat_features, spread_features)
    ),
    
    Sparsity = dplyr::if_else(
      is.na(stat_sparsity) | is.nan(stat_sparsity),
      "-",
      sprintf("%.3f (%.3f)", stat_sparsity, spread_sparsity)
    ),
    
    Alpha = dplyr::if_else(
      is.na(stat_alpha) | is.nan(stat_alpha),
      "-",
      sprintf("%.2f", stat_alpha)
    ),
    
    Lambda = dplyr::if_else(
      is.na(stat_lambda) | is.nan(stat_lambda),
      "-",
      sprintf("%.3f (%.3f)", stat_lambda, spread_lambda)
    ),
    
    d = dplyr::if_else(
      is.na(stat_d) | is.nan(stat_d),
      "-",
      sprintf("%.2f", stat_d)
    )
  ) |>
  dplyr::select(
    Dataset = dataset,
    Method = method,
    RMSE,
    MSE,
    `Time (s)`,
    Iter,
    `Sel. Groups`,
    `Sel. Features`,
    Sparsity,
    Alpha,
    Lambda,
    d
  ) |>
  dplyr::arrange(
    Dataset,
    dplyr::desc(Method == "SGLASSO"),
    RMSE
  )

Benchmark_Table_By_Dataset

# ------------------------------------------------------------
# 3. Separate tables for each dataset
# ------------------------------------------------------------

Benchmark_Tables_Separate <- split(
  Benchmark_Table_By_Dataset,
  Benchmark_Table_By_Dataset$Dataset
)

Benchmark_Tables_Separate$bardet
Benchmark_Tables_Separate$Birthwt
Benchmark_Tables_Separate$GenAtHum

# ------------------------------------------------------------
# 4. Save CSV files
# ------------------------------------------------------------

write.csv(
  Benchmark_Table_By_Dataset,
  "Benchmark_Table_By_Dataset.csv",
  row.names = FALSE
)

for (nm in names(Benchmark_Tables_Separate)) {
  write.csv(
    Benchmark_Tables_Separate[[nm]],
    paste0("Benchmark_Table_", nm, ".csv"),
    row.names = FALSE
  )
}

# ============================================================
# 5. Export combined table to LaTeX using kableExtra
# ============================================================

if (!requireNamespace("kableExtra", quietly = TRUE)) {
  install.packages("kableExtra")
}

library(kableExtra)

Latex_Table_By_Dataset <- Benchmark_Table_By_Dataset

colnames(Latex_Table_By_Dataset) <- c(
  "Dataset",
  "Method",
  "RMSE",
  "MSE",
  "Time (s)",
  "Iter",
  "Sel. Groups",
  "Sel. Features",
  "Sparsity",
  "$\\alpha$",
  "$\\lambda$",
  "$d$"
)

latex_table <- Latex_Table_By_Dataset |>
  kbl(
    format = "latex",
    booktabs = TRUE,
    longtable = FALSE,
    align = "llcccccccccc",
    caption = paste(
      "Benchmark comparison of grouped penalized regression methods",
      "across real-data examples."
    ),
    label = "tab:benchmark_real_by_dataset",
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
      "RMSE is relative to SGLASSO within each dataset and replication.",
      "MSE and runtime are on their original scales.",
      "Runtime refers to the final refit after tuning.",
      "Sparsity is defined as one minus the proportion of selected features. ",
      "$d$ is reported only for SGLASSO."
    ),
    general_title = "Note:",
    footnote_as_chunk = FALSE,
    threeparttable = TRUE,
    escape = FALSE
  )


cat(latex_table, file = "Benchmark_Table_By_Dataset.tex")
cat("LaTeX table saved as Benchmark_Table_By_Dataset.tex\n")

# ============================================================
# 6. Export separate LaTeX tables for each dataset
# ============================================================

for (nm in names(Benchmark_Tables_Separate)) {
  
  tab_nm <- Benchmark_Tables_Separate[[nm]] |>
    dplyr::select(-Dataset)
  
  colnames(tab_nm) <- c(
    "Method",
    "RMSE",
    "MSE",
    "Time (s)",
    "Iter",
    "Sel. Groups",
    "Sel. Features",
    "Sparsity",
    "$\\alpha$",
    "$\\lambda$",
    "$d$"
  )
  
  latex_nm <- tab_nm |>
    kbl(
      format = "latex",
      booktabs = TRUE,
      longtable = FALSE,
      align = "lcccccccccc",
      caption = paste0(
        "Benchmark comparison of grouped penalized regression methods for the ",
        nm,
        " dataset."
      ),
      label = paste0("tab:benchmark_real_", nm),
      escape = FALSE,
      linesep = ""
    ) |>
    kable_styling(
      latex_options = "hold_position",
      font_size = 9,
      full_width = FALSE,
      position = "center"
    ) |>
    row_spec(0, bold = TRUE) |>
    footnote(
      general = paste(
        "Values are reported as", summary_type,
        "(", spread_label, ").",
        "RMSE is computed relative to SGLASSO within each replication.",
        "Runtime is measured in elapsed seconds for the final non-cross-validated refit.",
        "Sparsity is defined as one minus the proportion of selected features.",
        "$d$ is reported only for SGLASSO."
      ),
      general_title = "Note:",
      footnote_as_chunk = FALSE,
      threeparttable = TRUE,
      escape = FALSE
    )
  
  cat(
    latex_nm,
    file = paste0("Benchmark_Table_", nm, ".tex")
  )
}

cat("Separate LaTeX tables saved for each dataset.\n")
