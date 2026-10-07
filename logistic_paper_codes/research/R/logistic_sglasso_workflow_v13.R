# V13 research-only continuation of the locally validated V12 workflow.
lsg_load_v13 <- function(root) {
  scope <- new.env(parent = globalenv())
  source(file.path(root, "R/logistic_sglasso_workflow_v12.R"), local = scope)
  e <- scope$lsg_load_v12(root)
  e$v13_original <- mget(c("classify_grpreg_path_v2", "select_valid_external_grid_v2",
    "audit_external_nonfinite_validation_v3", "fit_grpreg_validation_paths_v2"), envir = e)
  source(file.path(root, "R/logistic_sglasso_group_prefix_v13.R"), local = e)
  e$lsg_install_prefix_adapters_v13(e)
  e$lsg_configuration_v7 <- function(stage) {
    cfg <- e$lsg_configuration_v12(stage)
    cfg$integration_version <- "v13_local_only_group_lasso_prefix"
    cfg$group_lasso_prefix_policy <- "finite_budget_prefix_strict_interior_v13"
    cfg
  }
  e$lsg_source_files_v7 <- function() unique(c(e$lsg_source_files_v12(),
    "R/logistic_sglasso_group_prefix_v13.R", "R/logistic_sglasso_workflow_v13.R",
    "scripts/53_run_logistic_sglasso_integration_v13.R",
    "scripts/54_export_logistic_sglasso_tail125_inputs_v13.R",
    "scripts/51_export_logistic_sglasso_replay_inputs_v11.R",
    "tests/test_group_prefix_v13.R", "tests/test_tail125_export_v13.R",
    "tests/test_integration_outputs_v13.R",
    "LOGISTIC_SGLASSO_PREFIX_PROTOCOL_V13.md"))
  e
}
