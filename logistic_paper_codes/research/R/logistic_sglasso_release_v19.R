lsg_operational_preflight_v19 <- function(root, stage, cores, max_seconds,
                                          task_limit, interval) {
  lsg_assert_v7(stage %in% c("smoke", "production"), "Invalid V19 stage.")
  lsg_assert_v7(length(cores) == 1L && is.finite(cores) &&
    cores == floor(cores) && cores >= 1 &&
    cores <= if (stage == "smoke") 2 else 56,
    "Invalid V19 worker count.")
  lsg_assert_v7(length(max_seconds) == 1L && is.finite(max_seconds) &&
    max_seconds > 0 && max_seconds <= if (stage == "smoke") 300 else 18000,
    "Invalid V19 soft dispatch budget.")
  lsg_assert_v7(length(task_limit) == 1L && !is.na(task_limit) &&
    task_limit >= 0 && (is.infinite(task_limit) || task_limit == floor(task_limit)),
    "Invalid V19 task limit.")
  lsg_assert_v7(length(interval) == 1L && is.finite(interval) &&
    interval > 0 && interval <= 60,
    "Heartbeat interval must be positive and at most 60 seconds.")
  if (stage == "production") {
    lsg_assert_v7(nzchar(Sys.getenv("SLURM_JOB_ID")),
      "Production execution is reserved for the user-submitted TRUBA job.")
    lsg_assert_v7(identical(Sys.getenv("SLURM_CPUS_PER_TASK"), "56"),
      "Production requires 56 allocated CPUs.")
    lsg_validate_local_release_v19(root)
  }
  invisible(TRUE)
}

lsg_production_identity_v19 <- function(root) {
  cfg <- lsg_configuration_v7("production")
  design <- lsg_design_for_stage_v7(root, "production")
  list(
    configuration = cfg,
    design = design,
    tasks = lsg_task_grid_v7(design, cfg)
  )
}

lsg_validate_local_release_v19 <- function(root, quiet = FALSE) {
  path <- lsg_path_v7(root, "release/v19_a01/LOCAL_VALIDATED_v19.rds")
  lsg_assert_v7(file.exists(path),
    "Missing local V19 release receipt; no production computation allowed.")
  receipt <- readRDS(path)
  lsg_assert_v7(
    identical(receipt$schema, "local_validated_joint_production_v19") &&
      isTRUE(receipt$accepted) && nrow(receipt$checks) > 0L &&
      all(receipt$checks$passed %in% TRUE),
    "Invalid V19 release receipt."
  )
  lsg_assert_v7(
    identical(receipt$production_identity, lsg_production_identity_v19(root)),
    "Released V19 production design, grid, controls or seeds changed."
  )
  lsg_assert_v7(
    identical(receipt$sources, lsg_inventory_v7(root, lsg_source_files_v7())),
    "V19 source hashes changed; create a new locally validated release."
  )
  lsg_assert_v7(
    nrow(receipt$evidence) > 0L &&
      identical(receipt$evidence,
                lsg_inventory_v7(root, receipt$evidence$file)),
    "V19 local evidence size/checksum mismatch."
  )
  for (item in c("unit", "outputs")) {
    x <- readRDS(lsg_path_v7(root, receipt$test_receipts[[item]]))
    lsg_assert_v7(
      identical(x$sources, receipt$sources) &&
        all(x$checks$passed %in% TRUE),
      paste("Invalid V19 source-bound test receipt:", item)
    )
  }
  tail <- readRDS(lsg_path_v7(root, receipt$tail_diagnostic_result))
  lsg_assert_v7(
    isTRUE(tail$accepted) && all(tail$checks$passed %in% TRUE) &&
      identical(tail$sources, receipt$sources) &&
      identical(tail$test_fields_used, FALSE),
    "Invalid source-bound V19 task-234 tail diagnostic."
  )
  smoke <- readRDS(lsg_path_v7(root, paste0(
    "outputs/study/", receipt$smoke_version, "/study_specification.rds"
  )))
  lsg_validate_spec_v7(smoke, root, "smoke", receipt$smoke_version)
  lsg_assert_v7(
    identical(smoke$scientific_signature, receipt$smoke_signature) &&
      nrow(smoke$tasks) == 2L,
    "Released V19 full-grid smoke identity mismatch."
  )
  if (!quiet) {
    runtime <- lsg_runtime_v7()
    cat("Local frozen V19 release verified: TRUE; preflight model fits: 0\n")
    cat("Execution runtime:", runtime$r_version, runtime$platform, "\n")
    print(runtime$package_versions)
    cat("Library paths:\n", paste(.libPaths(), collapse = "\n"), "\n")
    cat("Frozen production: 8 scenarios x 50 = 400 tasks; 2400 selected rows.\n")
    cat("Finite lambda ratios:", lsg_lambda_relative_grid_v19(), "plus Inf\n")
  }
  invisible(receipt)
}
