# Research-only integration. Every overlay lives in a private R environment;
# neither frozen V7 functions nor the installed package are rebound globally.
lsg_load_v12 <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  e <- new.env(parent = globalenv())
  source(file.path(root, "R/logistic_sglasso_workflow_v7.R"), local = e)
  e$lsg_source_v7(root, e)
  extra <- c("logistic_sglasso_numerical_repair_v9.R",
             "logistic_sglasso_numerical_repair_v10.R",
             "logistic_sglasso_hybrid_solver_v11.R",
             "logistic_sglasso_hybrid_diagnostic_v11.R",
             "logistic_sglasso_integration_v12.R", "logistic_sglasso_progress_v12.R",
             "logistic_sglasso_dispatch_v12.R")
  for (f in extra) source(file.path(root, "R", f), local = e)
  originals <- c("lsg_configuration_v7", "lsg_source_files_v7", "lsg_required_packages_v7",
    "lsg_finite_fitters_v7", "lsg_fit_sglasso_joint_v7", "lsg_selected_sglasso_model_v7",
    "lsg_validate_tuning_payload_v7", "lsg_manifest_files_v7")
  e$v7_original <- mget(originals, envir = e)
  e$v12_extra <- extra
  e$lsg_configuration_v7 <- e$lsg_configuration_v12
  e$lsg_source_files_v7 <- e$lsg_source_files_v12
  e$lsg_required_packages_v7 <- function() c(e$v7_original$lsg_required_packages_v7(), "jsonlite")
  e$lsg_finite_fitters_v7 <- e$lsg_finite_fitters_v12
  e$lsg_fit_sglasso_joint_v7 <- e$lsg_fit_sglasso_joint_v12
  e$lsg_selected_sglasso_model_v7 <- e$lsg_selected_sglasso_model_v12
  e$lsg_validate_tuning_payload_v7 <- e$lsg_validate_tuning_payload_v12
  e$lsg_manifest_files_v7 <- function(output, version) {
    setdiff(e$v7_original$lsg_manifest_files_v7(output, version), c("progress.json", "progress.tsv"))
  }
  e
}
