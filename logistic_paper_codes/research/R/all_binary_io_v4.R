# Versioned V4 identity and checkpoint I/O; V3 artifacts are never reused.

allb_make_spec_v4 <- function(e, root, stage, version) {
  inherited <- allb_clone_with_bindings_v2(allb_make_spec_v3, list(
    allb_configuration_v3 = allb_configuration_v4,
    allb_source_files_v3 = allb_source_files_v4
  ))
  spec <- inherited(e, root, stage, version)
  spec$schema_version <- "all_binary_study_v4"
  spec$scientific_signature <- NULL
  identity <- spec[!names(spec) %in% c("runtime", "created_utc")]
  spec$scientific_signature <- allb_hash(e, identity)
  spec
}

allb_validate_spec_v4 <- function(e, spec, root, stage, version) {
  current <- allb_make_spec_v4(e, root, stage, version)
  fields <- c(
    "schema_version", "version", "stage", "configuration", "tasks",
    "data_identity", "source_manifest", "scientific_signature", "runtime"
  )
  allb_assert(identical(spec[fields], current[fields]),
    "ALL V4 scientific identity, sources, data, runtime, or seeds changed.")
  invisible(TRUE)
}

allb_validate_release_v4 <- function(e, root, version) {
  path <- file.path(root, "release", "all_binary_v4",
    "LOCAL_VALIDATED_all_binary_v4.rds")
  allb_assert(file.exists(path),
    "Missing locally validated ALL Binary V4 release receipt.")
  receipt <- readRDS(path)
  current <- allb_make_spec_v4(e, root, "production", version)
  fields <- c("scientific_signature", "source_manifest", "data_identity",
    "configuration", "tasks")
  allb_assert(
    identical(receipt$schema_version, "all_binary_portable_release_v4") &&
      isTRUE(receipt$accepted) && all(receipt$checks$passed) &&
      identical(receipt$production_version, version) &&
      identical(receipt$identity[fields], current[fields]),
    "ALL V4 local release does not match production identity."
  )
  invisible(receipt)
}

allb_valid_shard_v4 <- function(path, task, signature, configuration) {
  if (!file.exists(path) || dir.exists(path)) return(FALSE)
  tryCatch({
    shard <- readRDS(path)
    checks <- allb_validate_payload_v4(shard$payload, task, configuration)
    identical(shard$schema_version, "all_binary_shard_v4") &&
      identical(shard$scientific_signature, signature) &&
      identical(shard$task, task) && all(checks$passed)
  }, error = function(error) FALSE)
}

allb_progress_snapshot_v4 <- function(spec, out, status, started, phase,
                                       invocation_id,
                                       completed_this_invocation = 0L) {
  inherited <- allb_clone_with_bindings_v2(allb_progress_snapshot, list(
    allb_valid_shard = allb_valid_shard_v4
  ))
  result <- inherited(spec, out, status, started, phase, invocation_id,
    completed_this_invocation)
  result$schema_version <- "all_binary_progress_v4"
  result$fail_fast <- TRUE
  result
}

allb_worker_v4 <- function(task_index, spec, data, e, out) {
  task <- spec$tasks[task_index, , drop = FALSE]
  started <- proc.time()[["elapsed"]]
  tryCatch({
    payload <- allb_run_task_v4(e, data, task, spec$configuration)
    checks <- allb_validate_payload_v4(payload, task, spec$configuration)
    allb_assert(all(checks$passed), paste(
      "Task checks failed:", paste(checks$check[!checks$passed],
        collapse = ", ")
    ))
    shard <- list(
      schema_version = "all_binary_shard_v4",
      scientific_signature = spec$scientific_signature,
      task = task, checks = checks, payload = payload,
      runtime_seconds = proc.time()[["elapsed"]] - started,
      completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
    allb_atomic_rds(shard, file.path(out, "shards", task$shard_file))
    list(task_id = task$task_id[[1L]], status = "passed", message = "")
  }, error = function(error) list(
    task_id = task$task_id[[1L]], status = "failed",
    message = conditionMessage(error)
  ))
}

allb_run_pending_v4 <- function(spec, data, e, out, cores, interval) {
  allb_run_pending_fail_fast_v3(
    spec, data, e, out, cores, interval,
    valid_shard_fun = allb_valid_shard_v4,
    worker_fun = allb_worker_v4,
    progress_snapshot_fun = allb_progress_snapshot_v4,
    write_progress_fun = allb_write_progress_v3
  )
}

allb_finalize_v4 <- function(e, spec, out) {
  inherited <- allb_clone_with_bindings_v2(allb_finalize, list(
    allb_valid_shard = allb_valid_shard_v4,
    allb_final_names = allb_final_names_v3
  ))
  inherited(e, spec, out)
}

allb_verify_outputs_v4 <- function(spec, out) {
  inherited <- allb_clone_with_bindings_v2(allb_verify_outputs, list(
    allb_valid_shard = allb_valid_shard_v4,
    allb_final_names = allb_final_names_v3
  ))
  inherited(spec, out)
}

allb_run_v4 <- function(root, stage, version, cores = 1L, interval = 30) {
  e <- allb_load_environment_v4(root)
  if (stage == "production") {
    allb_assert(nzchar(Sys.getenv("SLURM_JOB_ID")),
      "Production is reserved for the user-submitted TRUBA job.")
    allb_assert(identical(Sys.getenv("SLURM_CPUS_PER_TASK"), "56"),
      "TRUBA Hamsi production requires exactly 56 allocated CPUs.")
    allb_validate_release_v4(e, root, version)
  }
  proposed <- allb_make_spec_v4(e, root, stage, version)
  out <- allb_output_directory(root, version)
  dir.create(file.path(out, "shards"), recursive = TRUE,
    showWarnings = FALSE)
  specification <- file.path(out, "study_specification.rds")
  if (file.exists(specification)) {
    spec <- readRDS(specification)
    allb_validate_spec_v4(e, spec, root, stage, version)
  } else {
    spec <- proposed
    allb_atomic_rds(spec, specification)
  }
  cat("Stage:", stage, "\nVersion:", version,
    "\nScientific signature:", spec$scientific_signature,
    "\nWork units:", nrow(spec$tasks),
    "\nRequested workers:", cores, "\nFail-fast: TRUE\n")
  print(spec$runtime)
  data <- allb_load_data(root, spec$configuration)
  status <- allb_run_pending_v4(spec, data, e, out, cores, interval)
  if (any(status == "failed")) stop(
    "An ALL V4 task failed; valid shards remain resumable.", call. = FALSE
  )
  marker_path <- file.path(out, paste0("COMPLETED_", spec$version, ".rds"))
  if (file.exists(marker_path)) {
    marker <- readRDS(marker_path)
    allb_assert(
      identical(marker$scientific_signature, spec$scientific_signature) &&
        marker$completed_tasks == nrow(spec$tasks) &&
        file.exists(file.path(out, "MANIFEST.csv")),
      "Existing ALL V4 completion evidence is invalid."
    )
  } else {
    marker <- allb_finalize_v4(e, spec, out)
  }
  allb_verify_outputs_v4(spec, out)
  progress <- jsonlite::read_json(file.path(out, "progress.json"),
    simplifyVector = TRUE)
  progress$phase <- "complete"
  progress$eta_seconds <- 0
  progress$eta_status <- "complete_verified"
  progress$heartbeat_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  progress$estimated_completion_utc <- progress$heartbeat_utc
  allb_write_progress_v3(progress, out)
  cat("Complete: TRUE\nOutput directory:", normalizePath(out),
    "\nCompleted tasks:", marker$completed_tasks, "\n")
  invisible(list(spec = spec, output = out, marker = marker))
}
