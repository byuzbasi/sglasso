# V19 production overlay: two additional upper-tail ratios and a schema-
# preserving complete-final-budget Group Lasso audit repair.
lsg_load_v19 <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  scope <- new.env(parent = globalenv())
  source(file.path(root, "R/logistic_sglasso_workflow_v18.R"), local = scope)
  e <- scope$lsg_load_v18(root)
  e$v19_base <- mget(c(
    "lsg_lambda_relative_grid_v7", "lsg_configuration_v7",
    "lsg_source_files_v7", "lsg_finite_fitters_v7",
    "lsg_selection_boundary_v7", "lsg_tuning_unit_checks_v7",
    "audit_external_nonfinite_validation_v3", "v13_extended_audit",
    "select_valid_external_grid_v2", "mark_external_selection_v2",
    "lsg_validate_tuning_payload_v7", "lsg_run_v18"
  ), envir = e)
  source(file.path(root, "R/logistic_sglasso_repair_v19.R"), local = e)
  source(file.path(root, "R/logistic_sglasso_tail_diagnostic_v19.R"), local = e)
  source(file.path(root, "R/logistic_sglasso_release_v19.R"), local = e)

  e$lsg_lambda_relative_grid_v7 <- e$lsg_lambda_relative_grid_v19
  e$lsg_configuration_v7 <- e$lsg_configuration_v19
  e$lsg_finite_fitters_v7 <- e$lsg_finite_fitters_v19
  e$lsg_selection_boundary_v7 <- e$lsg_selection_boundary_v19
  e$lsg_tuning_unit_checks_v7 <- e$lsg_tuning_unit_checks_v19
  e$audit_external_nonfinite_validation_v3 <-
    e$lsg_build_external_audit_v19()
  e$lsg_source_files_v7 <- function() unique(c(
    e$v19_base$lsg_source_files_v7(),
    "R/logistic_sglasso_repair_v19.R",
    "R/logistic_sglasso_tail_diagnostic_v19.R",
    "R/logistic_sglasso_release_v19.R",
    "R/logistic_sglasso_workflow_v19.R",
    "config/logistic_sglasso_task234_v18_anchors_v19.csv",
    "scripts/67_run_logistic_sglasso_production_v19.R",
    "scripts/68_validate_logistic_sglasso_production_v19.R",
    "scripts/69_run_logistic_sglasso_tail_diagnostic_v19.R",
    "tests/test_production_v19.R",
    "LOGISTIC_SGLASSO_PRODUCTION_PROTOCOL_V19.md",
    "truba/run_logistic_sglasso_production_v19.slurm",
    "truba/build_logistic_sglasso_production_bundle_v19.R",
    "truba/README_TRUBA_LOGISTIC_SGLASSO_V19.md"
  ))

  # Reuse the frozen dispatcher; only its release preflight resolves to V19.
  run <- e$v19_base$lsg_run_v18
  child <- new.env(parent = environment(run))
  child$lsg_operational_preflight_v15 <- e$lsg_operational_preflight_v19
  environment(run) <- child
  e$lsg_run_v19 <- run
  e
}
