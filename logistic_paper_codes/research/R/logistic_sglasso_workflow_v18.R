# Clean V18 production overlay: fixed 2048/1024 tail and complete-path
# Group Lasso safe-prefix repair. V15 remains immutable.
lsg_load_v18 <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  scope <- new.env(parent = globalenv())
  source(file.path(root, "R/logistic_sglasso_workflow_v15.R"), local = scope)
  e <- scope$lsg_load_v15(root)
  e$v18_base <- mget(c(
    "lsg_configuration_v7", "lsg_source_files_v7", "lsg_finite_fitters_v7",
    "lsg_selection_boundary_v7", "lsg_tuning_unit_checks_v7",
    "lsg_classify_grpreg_path_v9", "lsg_classify_grpreg_path_v13",
    "lsg_run_v15"
  ), envir = e)
  source(file.path(root, "R/logistic_sglasso_repair_v18.R"), local = e)
  source(file.path(root, "R/logistic_sglasso_release_v18.R"), local = e)

  e$lsg_lambda_relative_grid_v7 <- e$lsg_lambda_relative_grid_v18
  e$lsg_configuration_v7 <- e$lsg_configuration_v18
  e$lsg_finite_fitters_v7 <- e$lsg_finite_fitters_v18
  e$lsg_selection_boundary_v7 <- e$lsg_selection_boundary_v18
  e$lsg_tuning_unit_checks_v7 <- e$lsg_tuning_unit_checks_v18
  e$lsg_classify_grpreg_path_v9 <- e$lsg_classify_grpreg_path_v18
  e$lsg_source_files_v7 <- function() unique(c(
    e$v18_base$lsg_source_files_v7(),
    "R/logistic_sglasso_repair_v18.R",
    "R/logistic_sglasso_release_v18.R",
    "R/logistic_sglasso_workflow_v18.R",
    "scripts/65_run_logistic_sglasso_production_v18.R",
    "scripts/66_validate_logistic_sglasso_production_v18.R",
    "tests/test_production_v18.R",
    "LOGISTIC_SGLASSO_PRODUCTION_PROTOCOL_V18.md",
    "truba/run_logistic_sglasso_production_v18.slurm",
    "truba/build_logistic_sglasso_production_bundle_v18.R",
    "truba/README_TRUBA_LOGISTIC_SGLASSO_V18.md"
  ))

  # Run the unchanged dispatcher in a child environment where only the
  # operational preflight call resolves to the V18 release guard.
  run <- e$v18_base$lsg_run_v15
  child <- new.env(parent = environment(run))
  child$lsg_operational_preflight_v15 <- e$lsg_operational_preflight_v18
  environment(run) <- child
  e$lsg_run_v18 <- run
  e
}
