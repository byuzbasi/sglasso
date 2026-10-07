# Production orchestration only: the frozen V13 numerical methods are reused.
lsg_load_v15 <- function(root) {
  scope <- new.env(parent = globalenv())
  source(file.path(root, "R/logistic_sglasso_workflow_v13.R"), local = scope)
  e <- scope$lsg_load_v13(root)
  e$v15_base <- mget(c("lsg_configuration_v7", "lsg_task_grid_v7", "lsg_source_files_v7",
    "lsg_design_for_stage_v7", "lsg_run_v12", "lsg_progress_new_v12",
    "lsg_progress_emit_v12", "lsg_read_shard_v7", "lsg_collect_v7"), envir = e)
  for (f in c("logistic_sglasso_progress_v15.R", "logistic_sglasso_release_v15.R"))
    source(file.path(root, "R", f), local = e)
  e$v15_progress <- NULL
  e$lsg_configuration_v7 <- function(stage) {
    stopifnot(stage %in% c("smoke", "pilot", "production"))
    cfg <- e$v15_base$lsg_configuration_v7(if (stage == "smoke") "smoke" else "pilot")
    cfg$stage <- if (stage == "smoke") "smoke" else "production"
    cfg$replications_per_scenario <- if (stage == "smoke") 2L else 50L
    cfg$alpha_grid <- cfg$d_grid <- cfg$benchmark_alpha_grid <- seq(0, 1, .1)
    cfg$integration_version <- "v15_joint_r50_production"
    cfg
  }
  e$lsg_design_for_stage_v7 <- function(root, stage)
    e$v15_base$lsg_design_for_stage_v7(root, if (stage == "smoke") "smoke" else "pilot")
  e$lsg_task_grid_v7 <- function(design, configuration) {
    tasks <- e$v15_base$lsg_task_grid_v7(design, configuration)
    if (configuration$stage == "production") {
      frozen <- utils::read.csv(file.path(root, "config/logistic_sglasso_production_tasks_v15.csv"),
        stringsAsFactors = FALSE)
      e$lsg_assert_v7(identical(tasks, frozen), "Frozen 400-task/seed map changed.")
    }
    tasks
  }
  e$lsg_source_files_v7 <- function() unique(c(e$v15_base$lsg_source_files_v7(),
    "R/logistic_sglasso_workflow_v15.R", "R/logistic_sglasso_progress_v15.R",
    "R/logistic_sglasso_release_v15.R", "config/logistic_sglasso_production_tasks_v15.csv",
    "scripts/59_run_logistic_sglasso_production_v15.R",
    "scripts/60_validate_logistic_sglasso_production_v15.R", "tests/test_production_v15.R",
    "LOGISTIC_SGLASSO_PRODUCTION_PROTOCOL_V15.md",
    "truba/run_logistic_sglasso_production_v15.slurm",
    "truba/build_logistic_sglasso_production_bundle_v15.R",
    "truba/README_TRUBA_LOGISTIC_SGLASSO_V15.md"))
  e$lsg_progress_new_v12 <- e$lsg_progress_new_v15
  e$lsg_progress_emit_v12 <- e$lsg_progress_emit_v15
  e$lsg_progress_eta_v12 <- e$lsg_progress_eta_v15
  e$lsg_read_shard_v7 <- e$lsg_read_shard_progress_v15
  run <- e$lsg_replace_expression_v13(e$v15_base$lsg_run_v12,
    quote(stage <- "smoke"), quote(stage <- match.arg(stage, c("smoke", "production"))))
  run <- e$lsg_replace_expression_v13(run,
    quote(lsg_assert_v7(cores %in% 1:2 && is.finite(max_seconds) && max_seconds > 0 && max_seconds <= 120,
      "Local integration is restricted to at most two workers and a 120-second soft budget.")),
    quote(lsg_operational_preflight_v15(root, stage, cores, max_seconds, task_limit, interval)))
  run <- e$lsg_replace_expression_v13(run,
    quote(old <- lsg_collect_v7(output, spec, deep = TRUE, require_all = FALSE)),
    quote(old <- lsg_scan_existing_v15(output, spec, cores, interval)))
  formals(run) <- c(formals(run), alist(stage = "smoke"))
  e$lsg_run_v15 <- run
  e
}
