# Frozen small pilot; numerical fitting and scientific checks remain V13.
lsg_load_v14 <- function(root) {
  scope <- new.env(parent = globalenv())
  source(file.path(root, "R/logistic_sglasso_workflow_v13.R"), local = scope)
  e <- scope$lsg_load_v13(root)
  e$v14_base <- mget(c("lsg_configuration_v7", "lsg_task_grid_v7",
    "lsg_source_files_v7", "lsg_run_v12"), envir = e)
  source(file.path(root, "R/logistic_sglasso_release_v14.R"), local = e)
  e$lsg_configuration_v7 <- function(stage) {
    cfg <- e$v14_base$lsg_configuration_v7(stage)
    cfg$replications_per_scenario <- if (stage == "pilot") 3L else 2L
    cfg$alpha_grid <- cfg$d_grid <- cfg$benchmark_alpha_grid <- seq(0, 1, .1)
    cfg$integration_version <- "v14_joint_small_pilot"
    cfg
  }
  e$lsg_task_grid_v7 <- function(design, configuration) {
    tasks <- e$v14_base$lsg_task_grid_v7(design, configuration)
    if (configuration$stage == "pilot") {
      tasks$source_v7_task_id <- as.integer((tasks$scenario_index - 1L) * 20L + tasks$replication)
      frozen <- utils::read.csv(file.path(root, "config/logistic_sglasso_pilot_tasks_v14.csv"),
        stringsAsFactors = FALSE)
      e$lsg_assert_v7(identical(tasks, frozen), "Frozen 24-task/seed map differs from configuration.")
    }
    tasks
  }
  e$lsg_source_files_v7 <- function() unique(c(e$v14_base$lsg_source_files_v7(),
    "R/logistic_sglasso_workflow_v14.R", "R/logistic_sglasso_release_v14.R",
    "config/logistic_sglasso_pilot_tasks_v14.csv",
    "scripts/57_run_logistic_sglasso_pilot_v14.R",
    "scripts/58_validate_logistic_sglasso_pilot_v14.R",
    "tests/test_pilot_v14.R", "LOGISTIC_SGLASSO_PILOT_PROTOCOL_V14.md",
    "truba/run_logistic_sglasso_pilot_v14.slurm",
    "truba/build_logistic_sglasso_pilot_bundle_v14.R",
    "truba/README_TRUBA_LOGISTIC_SGLASSO_V14.md"))
  # Reuse the entire validated dispatcher, worker, checkpoint and finalizer.
  # Only the entry-stage and operational budget guards change. Exact single
  # AST matches fail closed if the frozen implementation has changed.
  run <- e$lsg_replace_expression_v13(e$v14_base$lsg_run_v12,
    quote(stage <- "smoke"), quote(stage <- match.arg(stage, c("smoke", "pilot"))))
  run <- e$lsg_replace_expression_v13(run,
    quote(lsg_assert_v7(cores %in% 1:2 && is.finite(max_seconds) && max_seconds > 0 && max_seconds <= 120,
      "Local integration is restricted to at most two workers and a 120-second soft budget.")),
    quote(lsg_operational_preflight_v14(root, stage, cores, max_seconds, task_limit, interval)))
  formals(run) <- c(formals(run), alist(stage = "smoke"))
  e$lsg_run_v14 <- run
  e
}
