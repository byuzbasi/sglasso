# Portable local-validation evidence; remote preflight never fits a model.
lsg_operational_preflight_v14 <- function(root, stage, cores, max_seconds, task_limit, interval) {
  lsg_assert_v7(length(cores) == 1L && is.finite(cores) && cores == floor(cores) &&
    cores >= 1 && cores <= if (stage == "smoke") 2 else 24, "Invalid worker count.")
  lsg_assert_v7(length(max_seconds) == 1L && is.finite(max_seconds) && max_seconds > 0 &&
    max_seconds <= if (stage == "smoke") 120 else 6600, "Invalid soft dispatch budget.")
  lsg_assert_v7(length(task_limit) == 1L && !is.na(task_limit) && task_limit >= 0 &&
    (is.infinite(task_limit) || task_limit == floor(task_limit)), "Invalid task limit.")
  lsg_assert_v7(length(interval) == 1L && is.finite(interval) && interval > 0 &&
    interval <= 60, "Heartbeat interval must be positive and at most 60 seconds.")
  if (stage == "pilot") lsg_validate_local_release_v14(root)
  invisible(TRUE)
}

lsg_pilot_identity_v14 <- function(root) {
  cfg <- lsg_configuration_v7("pilot")
  design <- lsg_design_for_stage_v7(root, "pilot")
  list(configuration = cfg, design = design, tasks = lsg_task_grid_v7(design, cfg))
}

lsg_validate_local_release_v14 <- function(root, quiet = FALSE) {
  path <- lsg_path_v7(root, "release/LOCAL_VALIDATED_v14.rds")
  lsg_assert_v7(file.exists(path), "Missing local release receipt; run local validation before transfer.")
  r <- readRDS(path)
  lsg_assert_v7(identical(r$schema, "local_validated_joint_pilot_v14") &&
    isTRUE(r$accepted) && is.data.frame(r$checks) && nrow(r$checks) > 0 &&
    all(r$checks$passed %in% TRUE), "Invalid/failed local release receipt.")
  lsg_assert_v7(identical(r$pilot_identity, lsg_pilot_identity_v14(root)),
    "Released pilot design, grid, controls or seeds changed.")
  lsg_assert_v7(identical(r$sources, lsg_inventory_v7(root, lsg_source_files_v7())),
    "Released source hashes changed; revalidate locally under a new release.")
  lsg_assert_v7(is.data.frame(r$evidence) && nrow(r$evidence) > 0 &&
    identical(r$evidence, lsg_inventory_v7(root, r$evidence$file)),
    "Local validation evidence checksum/size mismatch.")
  for (item in c("unit", "outputs")) {
    x <- readRDS(lsg_path_v7(root, r$test_receipts[[item]]))
    lsg_assert_v7(identical(x$sources, r$sources) && all(x$checks$passed %in% TRUE),
      paste("Test receipt source/acceptance mismatch:", item))
  }
  smoke <- readRDS(lsg_path_v7(root, paste0("outputs/study/", r$smoke_version,
    "/study_specification.rds")))
  lsg_validate_spec_v7(smoke, root, "smoke", r$smoke_version)
  lsg_assert_v7(identical(smoke$scientific_signature, r$smoke_signature) &&
    nrow(smoke$tasks) == 2L && identical(smoke$configuration$alpha_grid, seq(0, 1, .1)) &&
    identical(smoke$configuration$d_grid, seq(0, 1, .1)), "Released full-grid smoke identity mismatch.")
  # Different OS/R/BLAS is expected on TRUBA, not a reason to re-fit the smoke.
  # Runtime is recorded and becomes exact-match mandatory on each pilot resume.
  runtime <- lsg_runtime_v7()
  if (!quiet) {
    cat("Local frozen release verified: TRUE; preflight model fits: 0\n")
    cat("Local smoke runtime:", r$runtime$r_version, r$runtime$platform, "\n")
    cat("Execution runtime:", runtime$r_version, runtime$platform, "\n")
    print(runtime$package_versions)
    cat("Library paths:\n", paste(.libPaths(), collapse = "\n"), "\n")
    cat("Frozen pilot: 8 scenarios x 3 = 24 tasks; 144 selected method rows.\n")
  }
  invisible(r)
}
