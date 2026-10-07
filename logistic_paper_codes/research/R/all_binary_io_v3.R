# Checkpoint, fail-fast execution, release, and finalization for ALL Binary V3.

allb_make_spec_v3 <- function(e, root, stage, version) {
  configuration <- allb_configuration_v3(e, stage)
  data <- allb_load_data(root, configuration)
  tasks <- allb_tasks(configuration)
  source_manifest <- allb_inventory(e, root,
    unique(c(e$lsg_source_files_v7(), allb_source_files_v3())))
  identity <- list(
    schema_version = "all_binary_study_v3",
    version = version, stage = stage, configuration = configuration,
    tasks = tasks,
    data_identity = list(
      source_path_basename = basename(data$source_path),
      source_sha256 = data$source_sha256,
      dimensions = dim(data$X),
      sample_id_sha256 = allb_hash(e, data$sample_id),
      response_sha256 = allb_hash(e, data$y),
      predictor_id_sha256 = allb_hash(e, colnames(data$X)),
      predictor_group_sha256 = allb_hash(e, data$group),
      group_name_sha256 = allb_hash(e, data$group_name),
      preprocessing = data$preprocessing,
      provenance = data$provenance
    ),
    source_manifest = source_manifest
  )
  identity$scientific_signature <- allb_hash(e, identity)
  c(identity, list(
    runtime = allb_runtime(),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  ))
}

allb_validate_spec_v3 <- function(e, spec, root, stage, version) {
  current <- allb_make_spec_v3(e, root, stage, version)
  fields <- c(
    "schema_version", "version", "stage", "configuration", "tasks",
    "data_identity", "source_manifest", "scientific_signature"
  )
  allb_assert(identical(spec[fields], current[fields]),
    "ALL V3 scientific identity, data, sources, grids, or seeds changed.")
  allb_assert(identical(spec$runtime, current$runtime),
    "Runtime/package identity changed; resume under the original environment.")
  invisible(TRUE)
}

allb_validate_release_v3 <- function(e, root, version) {
  path <- file.path(root, "release", "all_binary_v3",
    "LOCAL_VALIDATED_all_binary_v3.rds")
  allb_assert(file.exists(path),
    "Missing local ALL Binary V3 release receipt; production is blocked.")
  receipt <- readRDS(path)
  current <- allb_make_spec_v3(e, root, "production", version)
  allb_assert(
    identical(receipt$schema_version, "all_binary_portable_release_v3") &&
      isTRUE(receipt$accepted) && all(receipt$checks$passed) &&
      identical(receipt$production_version, version) &&
      identical(receipt$source_manifest, current$source_manifest) &&
      identical(receipt$data_identity, current$data_identity) &&
      identical(receipt$production_configuration, current$configuration) &&
      identical(receipt$production_tasks, current$tasks),
    paste(
      "ALL Binary V3 local release, processed data, source,",
      "returned-path policy, common grid, or seeds changed."
    )
  )
  invisible(receipt)
}

allb_valid_shard_v3 <- function(path, task, signature, configuration) {
  if (!file.exists(path) || dir.exists(path)) return(FALSE)
  tryCatch({
    x <- readRDS(path)
    checks <- allb_validate_payload_v3(x$payload, task, configuration)
    identical(x$schema_version, "all_binary_shard_v3") &&
      identical(x$scientific_signature, signature) &&
      identical(x$task, task) && all(checks$passed)
  }, error = function(error) FALSE)
}

allb_progress_snapshot_v3 <- function(spec, out, status, started, phase,
                                       invocation_id,
                                       completed_this_invocation = 0L) {
  inherited <- allb_clone_with_bindings_v2(allb_progress_snapshot, list(
    allb_valid_shard = allb_valid_shard_v3
  ))
  result <- inherited(
    spec, out, status, started, phase, invocation_id,
    completed_this_invocation
  )
  result$schema_version <- "all_binary_progress_v3"
  result$fail_fast <- TRUE
  result
}

allb_write_progress_v3 <- function(snapshot, out) {
  allb_write_progress(snapshot, out)
}

allb_worker_v3 <- function(task_index, spec, data, e, out) {
  task <- spec$tasks[task_index, , drop = FALSE]
  started <- proc.time()[["elapsed"]]
  tryCatch({
    payload <- allb_run_task_v3(e, data, task, spec$configuration)
    checks <- allb_validate_payload_v3(payload, task, spec$configuration)
    allb_assert(all(checks$passed), paste(
      "Task checks failed:",
      paste(checks$check[!checks$passed], collapse = ", ")
    ))
    shard <- list(
      schema_version = "all_binary_shard_v3",
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

allb_run_pending_fail_fast_v3 <- function(
    spec, data, e, out, cores, interval,
    valid_shard_fun = allb_valid_shard_v3,
    worker_fun = allb_worker_v3,
    progress_snapshot_fun = allb_progress_snapshot_v3,
    write_progress_fun = allb_write_progress_v3
) {
  valid <- vapply(seq_len(nrow(spec$tasks)), function(i) {
    valid_shard_fun(
      file.path(out, "shards", spec$tasks$shard_file[i]),
      spec$tasks[i, , drop = FALSE], spec$scientific_signature,
      spec$configuration
    )
  }, logical(1))
  pending <- which(!valid)
  status <- rep("pending", nrow(spec$tasks))
  status[valid] <- "passed"
  cores <- min(as.integer(cores), length(pending))
  started <- proc.time()[["elapsed"]]
  invocation <- paste0(
    format(Sys.time(), "%Y%m%dT%H%M%S"), "_", Sys.getpid()
  )
  if (!length(pending)) {
    write_progress_fun(progress_snapshot_fun(
      spec, out, status, started, "validating_existing_outputs",
      invocation, 0L
    ), out)
    return(status)
  }

  jobs <- list()
  outcomes <- list()
  completed_now <- 0L
  next_task <- 1L
  last_update <- -Inf
  failure_detected <- FALSE
  launch <- function(index) {
    status[index] <<- "running"
    job <- parallel::mcparallel(
      worker_fun(index, spec, data, e, out), silent = TRUE
    )
    jobs[[as.character(job$pid)]] <<- list(
      process = job, task_index = index
    )
  }
  while (next_task <= length(pending) && length(jobs) < cores) {
    launch(pending[[next_task]])
    next_task <- next_task + 1L
  }

  repeat {
    collected <- if (length(jobs)) parallel::mccollect(
      lapply(jobs, `[[`, "process"), wait = FALSE
    ) else NULL
    if (length(collected)) {
      for (pid in names(collected)) {
        info <- jobs[[pid]]
        result <- collected[[pid]]
        outcomes[[as.character(info$task_index)]] <- result
        if (inherits(result, "try-error") || is.null(result$status)) {
          status[info$task_index] <- "failed"
          failure_detected <- TRUE
        } else {
          status[info$task_index] <- result$status
          if (identical(result$status, "passed")) {
            completed_now <- completed_now + 1L
          }
          if (identical(result$status, "failed")) {
            failure_detected <- TRUE
            cat(
              "Task ", result$task_id, " failed: ", result$message,
              "\n", sep = ""
            )
          }
        }
        jobs[[pid]] <- NULL
      }
    }

    if (failure_detected && length(jobs)) {
      mckill <- get("mckill", envir = asNamespace("parallel"))
      for (pid in names(jobs)) {
        info <- jobs[[pid]]
        try(mckill(info$process, signal = 15L), silent = TRUE)
        status[info$task_index] <- "pending"
        outcomes[[as.character(info$task_index)]] <- list(
          task_id = spec$tasks$task_id[[info$task_index]],
          status = "aborted_after_peer_failure",
          message = "Terminated by V3 fail-fast policy."
        )
      }
      try(parallel::mccollect(
        lapply(jobs, `[[`, "process"), wait = TRUE
      ), silent = TRUE)
      jobs <- list()
    }

    if (!failure_detected) {
      while (next_task <= length(pending) && length(jobs) < cores) {
        launch(pending[[next_task]])
        next_task <- next_task + 1L
      }
    }

    elapsed <- proc.time()[["elapsed"]] - started
    finished <- if (failure_detected) !length(jobs) else
      !length(jobs) && next_task > length(pending)
    if (elapsed - last_update >= interval || finished) {
      phase <- if (finished) {
        if (failure_detected || any(status == "failed")) {
          "failed"
        } else "fitting_complete"
      } else "fitting"
      write_progress_fun(progress_snapshot_fun(
        spec, out, status, started, phase, invocation, completed_now
      ), out)
      last_update <- elapsed
    }
    if (finished) break
    Sys.sleep(1)
  }

  allb_atomic_rds(list(
    invocation_id = invocation, status = status, outcomes = outcomes,
    fail_fast_triggered = failure_detected,
    scientific_signature = spec$scientific_signature
  ), file.path(out, "attempts", paste0(invocation, ".rds")))
  status
}

allb_run_pending_v3 <- function(spec, data, e, out, cores, interval) {
  allb_run_pending_fail_fast_v3(spec, data, e, out, cores, interval)
}

allb_final_names_v3 <- function() c(
  allb_final_names(), "fair_tuning_audit"
)

allb_finalize_v3 <- function(e, spec, out) {
  inherited <- allb_clone_with_bindings_v2(allb_finalize, list(
    allb_valid_shard = allb_valid_shard_v3,
    allb_final_names = allb_final_names_v3
  ))
  inherited(e, spec, out)
}

allb_verify_outputs_v3 <- function(spec, out) {
  inherited <- allb_clone_with_bindings_v2(allb_verify_outputs, list(
    allb_valid_shard = allb_valid_shard_v3,
    allb_final_names = allb_final_names_v3
  ))
  inherited(spec, out)
}

allb_run_v3 <- function(root, stage, version, cores = 1L, interval = 30) {
  e <- allb_load_environment_v3(root)
  if (stage == "production") {
    allb_assert(nzchar(Sys.getenv("SLURM_JOB_ID")),
      "Production is reserved for the user-submitted TRUBA job.")
    allb_assert(identical(Sys.getenv("SLURM_CPUS_PER_TASK"), "56"),
      "TRUBA Hamsi production requires exactly 56 allocated CPUs.")
    allb_validate_release_v3(e, root, version)
  }
  proposed <- allb_make_spec_v3(e, root, stage, version)
  out <- allb_output_directory(root, version)
  dir.create(file.path(out, "shards"), recursive = TRUE,
    showWarnings = FALSE)
  specification <- file.path(out, "study_specification.rds")
  if (file.exists(specification)) {
    spec <- readRDS(specification)
    allb_validate_spec_v3(e, spec, root, stage, version)
  } else {
    spec <- proposed
    allb_atomic_rds(spec, specification)
  }
  cat(
    "Stage:", stage, "\nVersion:", version,
    "\nScientific signature:", spec$scientific_signature,
    "\nWork units:", nrow(spec$tasks),
    "\nRequested workers:", cores,
    "\nCommon lambda ratios:",
    paste(spec$configuration$fair_tuning$common_relative_grid,
      collapse = ","),
    "\nFail-fast: TRUE\n"
  )
  print(spec$runtime)
  data <- allb_load_data(root, spec$configuration)
  status <- allb_run_pending_v3(spec, data, e, out, cores, interval)
  if (any(status == "failed")) stop(
    paste(
      "An ALL V3 task failed; no new tasks were launched after the first",
      "failure and valid signature-matching shards remain resumable."
    ), call. = FALSE
  )
  marker_path <- file.path(out, paste0("COMPLETED_", spec$version, ".rds"))
  if (file.exists(marker_path)) {
    marker <- readRDS(marker_path)
    allb_assert(
      identical(marker$scientific_signature, spec$scientific_signature) &&
        marker$completed_tasks == nrow(spec$tasks) &&
        file.exists(file.path(out, "MANIFEST.csv")),
      "Existing ALL V3 completion evidence is invalid."
    )
  } else {
    marker <- allb_finalize_v3(e, spec, out)
  }
  allb_verify_outputs_v3(spec, out)
  progress <- jsonlite::read_json(file.path(out, "progress.json"),
    simplifyVector = TRUE)
  progress$phase <- "complete"
  progress$eta_seconds <- 0
  progress$eta_status <- "complete_verified"
  progress$heartbeat_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  progress$estimated_completion_utc <- progress$heartbeat_utc
  allb_write_progress_v3(progress, out)
  cat(
    "Complete: TRUE\nOutput directory:", normalizePath(out),
    "\nCompleted tasks:", marker$completed_tasks, "\n"
  )
  invisible(list(spec = spec, output = out, marker = marker))
}
