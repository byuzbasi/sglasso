# Versioned V5 local-smoke identity and checkpoint I/O. Production is disabled.

allb_make_spec_v5 <- function(e, root, stage, version) {
  inherited <- allb_clone_with_bindings_v2(allb_make_spec_v3, list(
    allb_configuration_v3 = allb_configuration_v5,
    allb_source_files_v3 = allb_source_files_v5
  ))
  spec <- inherited(e, root, stage, version)
  spec$schema_version <- "all_binary_study_v5"
  spec$scientific_signature <- NULL
  identity <- spec[!names(spec) %in% c("runtime", "created_utc")]
  spec$scientific_signature <- allb_hash(e, identity)
  spec
}

allb_validate_spec_v5 <- function(e, spec, root, stage, version) {
  current <- allb_make_spec_v5(e, root, stage, version)
  fields <- c(
    "schema_version", "version", "stage", "configuration", "tasks",
    "data_identity", "source_manifest", "scientific_signature", "runtime"
  )
  allb_assert(identical(spec[fields], current[fields]),
    "ALL V5 scientific identity, sources, data, runtime, or seeds changed.")
  invisible(TRUE)
}

allb_valid_shard_v5 <- function(path, task, signature, configuration) {
  if (!file.exists(path) || dir.exists(path)) return(FALSE)
  tryCatch({
    shard <- readRDS(path)
    checks <- allb_validate_payload_v5(shard$payload, task, configuration)
    identical(shard$schema_version, "all_binary_shard_v5") &&
      identical(shard$scientific_signature, signature) &&
      identical(shard$task, task) && all(checks$passed)
  }, error = function(error) FALSE)
}

allb_progress_snapshot_v5 <- function(spec, out, status, started, phase,
                                       invocation_id,
                                       completed_this_invocation = 0L) {
  inherited <- allb_clone_with_bindings_v2(allb_progress_snapshot, list(
    allb_valid_shard = allb_valid_shard_v5
  ))
  result <- inherited(spec, out, status, started, phase, invocation_id,
    completed_this_invocation)
  result$schema_version <- "all_binary_progress_v5"
  result$fail_fast <- TRUE
  result
}

allb_worker_v5 <- function(task_index, spec, data, e, out) {
  task <- spec$tasks[task_index, , drop = FALSE]
  started <- proc.time()[["elapsed"]]
  tryCatch({
    payload <- allb_run_task_v5(e, data, task, spec$configuration)
    checks <- allb_validate_payload_v5(payload, task, spec$configuration)
    allb_assert(all(checks$passed), paste(
      "Task checks failed:", paste(checks$check[!checks$passed],
        collapse = ", ")
    ))
    shard <- list(
      schema_version = "all_binary_shard_v5",
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

allb_run_pending_v5 <- function(spec, data, e, out, cores, interval) {
  allb_run_pending_fail_fast_v3(
    spec, data, e, out, cores, interval,
    valid_shard_fun = allb_valid_shard_v5,
    worker_fun = allb_worker_v5,
    progress_snapshot_fun = allb_progress_snapshot_v5,
    write_progress_fun = allb_write_progress_v3
  )
}

allb_finalize_v5 <- function(e, spec, out) {
  inherited <- allb_clone_with_bindings_v2(allb_finalize, list(
    allb_valid_shard = allb_valid_shard_v5,
    allb_final_names = allb_final_names_v3
  ))
  inherited(e, spec, out)
}

allb_verify_outputs_v5 <- function(spec, out) {
  inherited <- allb_clone_with_bindings_v2(allb_verify_outputs, list(
    allb_valid_shard = allb_valid_shard_v5,
    allb_final_names = allb_final_names_v3
  ))
  inherited(spec, out)
}

allb_run_v5 <- function(root, stage, version, cores = 1L, interval = 30) {
  allb_assert(identical(stage, "smoke"),
    "ALL V5 is a local validation candidate; production is disabled.")
  e <- allb_load_environment_v5(root)
  proposed <- allb_make_spec_v5(e, root, stage, version)
  out <- allb_output_directory(root, version)
  dir.create(file.path(out, "shards"), recursive = TRUE,
    showWarnings = FALSE)
  specification <- file.path(out, "study_specification.rds")
  if (file.exists(specification)) {
    spec <- readRDS(specification)
    allb_validate_spec_v5(e, spec, root, stage, version)
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
  status <- allb_run_pending_v5(spec, data, e, out, cores, interval)
  if (any(status == "failed")) stop(
    "An ALL V5 task failed; valid shards remain resumable.", call. = FALSE
  )
  marker_path <- file.path(out, paste0("COMPLETED_", spec$version, ".rds"))
  if (file.exists(marker_path)) {
    marker <- readRDS(marker_path)
    allb_assert(
      identical(marker$scientific_signature, spec$scientific_signature) &&
        marker$completed_tasks == nrow(spec$tasks) &&
        file.exists(file.path(out, "MANIFEST.csv")),
      "Existing ALL V5 completion evidence is invalid."
    )
  } else {
    marker <- allb_finalize_v5(e, spec, out)
  }
  allb_verify_outputs_v5(spec, out)
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
