lsg_operational_preflight_v18 <- function(root, stage, cores, max_seconds,
                                          task_limit, interval) {
  lsg_assert_v7(stage %in% c("smoke", "production"), "Invalid V18 stage.")
  lsg_assert_v7(length(cores) == 1L && is.finite(cores) && cores == floor(cores) &&
    cores >= 1 && cores <= if (stage == "smoke") 2 else 56,
    "Invalid V18 worker count.")
  lsg_assert_v7(length(max_seconds) == 1L && is.finite(max_seconds) &&
    max_seconds > 0 && max_seconds <= if (stage == "smoke") 300 else 18000,
    "Invalid V18 soft dispatch budget.")
  lsg_assert_v7(length(task_limit) == 1L && !is.na(task_limit) && task_limit >= 0 &&
    (is.infinite(task_limit) || task_limit == floor(task_limit)),
    "Invalid V18 task limit.")
  lsg_assert_v7(length(interval) == 1L && is.finite(interval) && interval > 0 &&
    interval <= 60, "Heartbeat interval must be positive and at most 60 seconds.")
  if (stage == "production") {
    lsg_assert_v7(nzchar(Sys.getenv("SLURM_JOB_ID")),
      "Production execution is reserved for the user-submitted TRUBA job.")
    lsg_assert_v7(identical(Sys.getenv("SLURM_CPUS_PER_TASK"), "56"),
      "Production requires 56 allocated CPUs.")
    lsg_validate_local_release_v18(root)
  }
  invisible(TRUE)
}

lsg_production_identity_v18 <- function(root) {
  cfg <- lsg_configuration_v7("production")
  design <- lsg_design_for_stage_v7(root, "production")
  list(configuration = cfg, design = design,
       tasks = lsg_task_grid_v7(design, cfg))
}

lsg_validate_local_release_v18 <- function(root, quiet = FALSE) {
  path <- lsg_path_v7(root, "release/v18_a01/LOCAL_VALIDATED_v18.rds")
  lsg_assert_v7(file.exists(path),
    "Missing local V18 release receipt; no production computation allowed.")
  receipt <- readRDS(path)
  lsg_assert_v7(
    identical(receipt$schema, "local_validated_joint_production_v18") &&
      isTRUE(receipt$accepted) && nrow(receipt$checks) > 0L &&
      all(receipt$checks$passed %in% TRUE),
    "Invalid V18 release receipt."
  )
  lsg_assert_v7(
    identical(receipt$production_identity, lsg_production_identity_v18(root)),
    "Released V18 production design, grid, controls or seeds changed."
  )
  lsg_assert_v7(
    identical(receipt$sources, lsg_inventory_v7(root, lsg_source_files_v7())),
    "V18 source hashes changed; create a new locally validated release."
  )
  lsg_assert_v7(
    nrow(receipt$evidence) > 0L &&
      identical(receipt$evidence,
                lsg_inventory_v7(root, receipt$evidence$file)),
    "V18 local evidence size/checksum mismatch."
  )
  for (item in c("unit", "outputs")) {
    x <- readRDS(lsg_path_v7(root, receipt$test_receipts[[item]]))
    lsg_assert_v7(identical(x$sources, receipt$sources) &&
      all(x$checks$passed %in% TRUE),
      paste("Invalid V18 source-bound test receipt:", item))
  }
  smoke <- readRDS(lsg_path_v7(root, paste0(
    "outputs/study/", receipt$smoke_version, "/study_specification.rds")))
  lsg_validate_spec_v7(smoke, root, "smoke", receipt$smoke_version)
  lsg_assert_v7(
    identical(smoke$scientific_signature, receipt$smoke_signature) &&
      nrow(smoke$tasks) == 2L,
    "Released V18 full-grid smoke identity mismatch."
  )
  if (!quiet) {
    runtime <- lsg_runtime_v7()
    cat("Local frozen V18 release verified: TRUE; preflight model fits: 0\n")
    cat("Execution runtime:", runtime$r_version, runtime$platform, "\n")
    print(runtime$package_versions)
    cat("Library paths:\n", paste(.libPaths(), collapse = "\n"), "\n")
    cat("Frozen production: 8 scenarios x 50 = 400 tasks; 2400 selected rows.\n")
    cat("Finite lambda ratios:", lsg_lambda_relative_grid_v18(), "plus Inf\n")
  }
  invisible(receipt)
}
