############################################################
# Real-data tables for the SGLASSO paper.
############################################################

options(stringsAsFactors = FALSE)

library(dplyr)
library(knitr)
library(kableExtra)

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  key <- paste0("--", name, "=")
  hit <- grep(paste0("^", key), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(key, "", hit[[1]])
}

base_dir <- get_arg("input-dir", Sys.getenv("REAL_DATA_DIR", "results/real_data"))
table_dir <- get_arg("table-dir", file.path(base_dir, "tables"))
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

summary_file <- file.path(base_dir, "ALL_summary.rds")
if (!file.exists(summary_file)) {
  stop("Cannot find summary results file: ", summary_file, call. = FALSE)
}

summary_results <- readRDS(summary_file)

method_levels <- c("ADELIE", "GLASSO", "GENET", "SGLASSO", "GSCAD", "GMCP")

fmt_num <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    "--",
    formatC(x, digits = digits, format = "f")
  )
}

fmt_time <- function(x) {
  ifelse(
    is.na(x),
    "--",
    ifelse(x < 0.01, formatC(x, digits = 4, format = "f"),
           formatC(x, digits = 3, format = "f"))
  )
}

table_data <- summary_results %>%
  filter(method %in% method_levels) %>%
  mutate(
    dataset = factor(dataset, levels = unique(summary_results$dataset)),
    method = factor(method, levels = method_levels)
  ) %>%
  group_by(dataset) %>%
  mutate(
    sglasso_mse = mse_median[method == "SGLASSO"][1],
    rel_mse_sglasso = mse_median / sglasso_mse
  ) %>%
  ungroup() %>%
  arrange(dataset, method) %>%
  transmute(
    dataset,
    method = as.character(method),
    nrep_success,
    mse_median,
    mse_sd,
    mse_se = mse_sd / sqrt(nrep_success),
    rel_mse_sglasso,
    selected_groups_median,
    selected_features_median,
    time_median,
    alpha_median,
    d_median
  )

write.csv(table_data, file.path(table_dir, "real_data_summary_table.csv"),
          row.names = FALSE)

latex_rows <- table_data %>%
  mutate(
    mse_value_tex = fmt_num(mse_median, 4),
    mse_se_tex = fmt_num(mse_se, 4),
    mse_tex = paste0(mse_value_tex, " (", mse_se_tex, ")"),
    rel_mse_tex = fmt_num(rel_mse_sglasso, 3),
    groups_tex = fmt_num(selected_groups_median, 1),
    features_tex = fmt_num(selected_features_median, 1),
    time_tex = fmt_time(time_median),
    alpha_tex = fmt_num(alpha_median, 1),
    d_tex = fmt_num(d_median, 1)
  )

make_kable_data <- function(data) {
  data %>%
    transmute(
      Method = method,
      `MSE (SE)` = mse_tex,
      `Rel. MSE` = rel_mse_tex,
      `$\\widehat{s}$` = groups_tex,
      `$\\widehat{p}$` = features_tex,
      Time = time_tex,
      `$\\alpha$` = alpha_tex,
      `$d$` = d_tex
    )
}

make_real_data_table <- function(data, dataset_name) {
  dataset_display <- dplyr::case_when(
    dataset_name == "bardet" ~ "Bardet",
    TRUE ~ dataset_name
  )
  caption_text <- paste0(
    "Real-data prediction and selection performance for the ",
    dataset_display,
    " dataset."
  )
  label_text <- paste0("real_data_", gsub("[^A-Za-z0-9]+", "_", tolower(dataset_name)))

  make_kable_data(data) %>%
    knitr::kable(
      format = "latex",
      booktabs = TRUE,
      escape = FALSE,
      align = c("l", rep("r", 8)),
      caption = caption_text,
      label = label_text
    ) %>%
    kableExtra::kable_styling(
      latex_options = c("scale_down"),
      position = "center",
      font_size = 9
    )
}

for (dataset_name in levels(latex_rows$dataset)) {
  dataset_rows <- dplyr::filter(latex_rows, dataset == dataset_name)
  dataset_stub <- gsub("[^A-Za-z0-9]+", "_", tolower(dataset_name))
  write.csv(
    dplyr::filter(table_data, as.character(dataset) == dataset_name),
    file.path(table_dir, paste0("real_data_", dataset_stub, "_table.csv")),
    row.names = FALSE
  )
  writeLines(
    as.character(make_real_data_table(dataset_rows, dataset_name)),
    file.path(table_dir, paste0("real_data_", dataset_stub, "_table.tex"))
  )
}

message("Real-data table files saved in: ", normalizePath(table_dir, mustWork = FALSE))
