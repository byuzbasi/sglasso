lsg_operational_preflight_v15 <- function(root, stage, cores, max_seconds, task_limit, interval) {
  lsg_assert_v7(stage %in% c("smoke", "production"), "Invalid stage.")
  lsg_assert_v7(length(cores) == 1L && is.finite(cores) && cores == floor(cores) &&
    cores >= 1 && cores <= if (stage == "smoke") 2 else 56, "Invalid worker count.")
  lsg_assert_v7(length(max_seconds) == 1L && is.finite(max_seconds) && max_seconds > 0 &&
    max_seconds <= if (stage == "smoke") 120 else 16200, "Invalid soft dispatch budget.")
  lsg_assert_v7(length(task_limit) == 1L && !is.na(task_limit) && task_limit >= 0 &&
    (is.infinite(task_limit) || task_limit == floor(task_limit)), "Invalid task limit.")
  lsg_assert_v7(length(interval) == 1L && is.finite(interval) && interval > 0 &&
    interval <= 60, "Heartbeat interval must be positive and at most 60 seconds.")
  if (stage == "production") {
    lsg_assert_v7(nzchar(Sys.getenv("SLURM_JOB_ID")),
      "Production execution is reserved for the user-submitted TRUBA job.")
    lsg_assert_v7(identical(Sys.getenv("SLURM_CPUS_PER_TASK"), "56"), "Production requires 56 allocated CPUs.")
    lsg_validate_local_release_v15(root)
  }
  invisible(TRUE)
}

lsg_production_identity_v15 <- function(root) {
  cfg <- lsg_configuration_v7("production")
  design <- lsg_design_for_stage_v7(root, "production")
  list(configuration = cfg, design = design, tasks = lsg_task_grid_v7(design, cfg))
}

lsg_validate_local_release_v15 <- function(root, quiet = FALSE) {
  path <- lsg_path_v7(root, "release/v15_a01/LOCAL_VALIDATED_v15.rds")
  lsg_assert_v7(file.exists(path), "Missing local V15 release receipt; no production computation allowed.")
  r <- readRDS(path)
  lsg_assert_v7(identical(r$schema, "local_validated_joint_production_v15") &&
    isTRUE(r$accepted) && nrow(r$checks) > 0 && all(r$checks$passed %in% TRUE), "Invalid release receipt.")
  lsg_assert_v7(identical(r$production_identity, lsg_production_identity_v15(root)),
    "Released production design, grid, controls or seeds changed.")
  lsg_assert_v7(identical(r$sources, lsg_inventory_v7(root, lsg_source_files_v7())),
    "Source hashes changed; create a new locally validated release.")
  lsg_assert_v7(nrow(r$evidence) > 0 && identical(r$evidence, lsg_inventory_v7(root, r$evidence$file)),
    "Local evidence size/checksum mismatch.")
  for (item in c("unit", "outputs")) {
    x <- readRDS(lsg_path_v7(root, r$test_receipts[[item]]))
    lsg_assert_v7(identical(x$sources, r$sources) && all(x$checks$passed %in% TRUE),
      paste("Invalid source-bound test receipt:", item))
  }
  smoke <- readRDS(lsg_path_v7(root, paste0("outputs/study/", r$smoke_version, "/study_specification.rds")))
  lsg_validate_spec_v7(smoke, root, "smoke", r$smoke_version)
  lsg_assert_v7(identical(smoke$scientific_signature, r$smoke_signature) && nrow(smoke$tasks) == 2L,
    "Released full-grid smoke identity mismatch.")
  runtime <- lsg_runtime_v7()
  if (!quiet) {
    cat("Local frozen release verified: TRUE; preflight model fits: 0\n")
    cat("Local smoke runtime:", r$runtime$r_version, r$runtime$platform, "\n")
    cat("Execution runtime:", runtime$r_version, runtime$platform, "\n")
    print(runtime$package_versions)
    cat("Library paths:\n", paste(.libPaths(), collapse = "\n"), "\n")
    cat("Frozen production: 8 scenarios x 50 = 400 tasks; 2400 selected method rows.\n")
  }
  invisible(r)
}
