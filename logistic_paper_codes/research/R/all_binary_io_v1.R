# Checkpoint, progress, release, and finalization support for ALL Binary V1.

allb_runtime <- function() gab_runtime()
allb_hash <- function(e, object) e$lsg_hash_v7(object)
allb_inventory <- function(e, root, files) e$lsg_inventory_v7(root, files)
allb_atomic_rds <- gab_atomic_rds
allb_atomic_json <- gab_atomic_json

allb_output_directory <- function(root, version) {
  allb_assert(grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", version),
    "Unsafe ALL binary version name.")
  file.path(root, "outputs", "study", version)
}

allb_make_spec <- function(e, root, stage, version) {
  configuration <- allb_configuration(e, stage)
  data <- allb_load_data(root, configuration)
  tasks <- allb_tasks(configuration)
  source_manifest <- allb_inventory(e, root,
    unique(c(e$lsg_source_files_v7(), allb_source_files())))
  identity <- list(
    schema_version = "all_binary_study_v1",
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
  c(identity, list(runtime = allb_runtime(),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)))
}

allb_validate_spec <- function(e, spec, root, stage, version) {
  current <- allb_make_spec(e, root, stage, version)
  fields <- c("schema_version", "version", "stage", "configuration", "tasks",
    "data_identity", "source_manifest", "scientific_signature")
  allb_assert(identical(spec[fields], current[fields]),
    "ALL scientific identity, data, sources, grids, or seeds changed.")
  allb_assert(identical(spec$runtime, current$runtime),
    "Runtime/package identity changed; resume under the original environment.")
  invisible(TRUE)
}

allb_validate_release <- function(e, root, version) {
  path <- file.path(root, "release", "all_binary_v1",
    "LOCAL_VALIDATED_all_binary_v1.rds")
  allb_assert(file.exists(path),
    "Missing local ALL binary release receipt; production is blocked.")
  receipt <- readRDS(path)
  current <- allb_make_spec(e, root, "production", version)
  allb_assert(identical(receipt$schema_version,
      "all_binary_portable_release_v1") && isTRUE(receipt$accepted) &&
    all(receipt$checks$passed) && identical(receipt$production_version, version) &&
    identical(receipt$source_manifest, current$source_manifest) &&
    identical(receipt$data_identity, current$data_identity) &&
    identical(receipt$production_configuration, current$configuration) &&
    identical(receipt$production_tasks, current$tasks),
    "ALL binary local release, processed data, source, grid, or seeds changed.")
  invisible(receipt)
}

allb_valid_shard <- function(path, task, signature, configuration) {
  if (!file.exists(path) || dir.exists(path)) return(FALSE)
  tryCatch({
    x <- readRDS(path)
    checks <- allb_validate_payload(x$payload, task, configuration)
    identical(x$schema_version, "all_binary_shard_v1") &&
      identical(x$scientific_signature, signature) &&
      identical(x$task, task) && all(checks$passed)
  }, error = function(error) FALSE)
}

allb_progress_snapshot <- function(spec, out, status, started, phase,
                                    invocation_id,
                                    completed_this_invocation = 0L) {
  valid <- vapply(seq_len(nrow(spec$tasks)), function(i) {
    allb_valid_shard(file.path(out, "shards", spec$tasks$shard_file[i]),
      spec$tasks[i, , drop = FALSE], spec$scientific_signature,
      spec$configuration)
  }, logical(1))
  elapsed <- proc.time()[["elapsed"]] - started
  completed <- sum(valid)
  running <- sum(status == "running")
  failed <- sum(status == "failed")
  pending <- nrow(spec$tasks) - completed - running - failed
  rate <- if (elapsed > 0) completed_this_invocation / elapsed else NA_real_
  eta <- if (is.finite(rate) && rate > 0 && pending + running > 0) {
    (pending + running) / rate
  } else NA_real_
  now <- Sys.time()
  remaining_slurm <- suppressWarnings(as.numeric(Sys.getenv(
    "SLURM_JOB_END_TIME", unset = NA_character_))) - as.numeric(now)
  list(
    schema_version = "all_binary_progress_v1",
    run_id = spec$version,
    scientific_signature = spec$scientific_signature,
    slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = "local"),
    invocation_id = invocation_id, phase = phase,
    work_unit = "outer_replication_all_six_methods_validated",
    total = nrow(spec$tasks), completed = completed, running = running,
    failed = failed, pending = pending,
    percent_completed = 100 * completed / nrow(spec$tasks),
    elapsed_seconds = elapsed,
    completed_this_invocation = completed_this_invocation,
    throughput_tasks_per_second = if (is.finite(rate)) rate else NULL,
    eta_seconds = if (is.finite(eta)) eta else NULL,
    estimated_completion_utc = if (is.finite(eta))
      format(now + eta, tz = "UTC", usetz = TRUE) else NULL,
    eta_status = if (is.finite(eta))
      "measured_current_invocation" else "estimating_or_unknown",
    eta_scope = "fitting_and_task_validation_only_excludes_final_aggregation",
    slurm_remaining_seconds = if (is.finite(remaining_slurm))
      remaining_slurm else NULL,
    heartbeat_utc = format(now, tz = "UTC", usetz = TRUE),
    most_recently_completed_work_unit = if (any(valid))
      tail(spec$tasks$key[valid], 1L) else NULL
  )
}

allb_write_progress <- function(snapshot, out) {
  allb_atomic_json(snapshot, file.path(out, "progress.json"))
  row <- data.frame(
    timestamp_utc = snapshot$heartbeat_utc,
    phase = snapshot$phase, total = snapshot$total,
    completed = snapshot$completed, running = snapshot$running,
    failed = snapshot$failed, pending = snapshot$pending,
    percent_completed = snapshot$percent_completed,
    elapsed_seconds = snapshot$elapsed_seconds,
    throughput_tasks_per_second = snapshot$throughput_tasks_per_second %||%
      NA_real_,
    eta_seconds = snapshot$eta_seconds %||% NA_real_,
    slurm_remaining_seconds = snapshot$slurm_remaining_seconds %||% NA_real_,
    stringsAsFactors = FALSE
  )
  path <- file.path(out, "progress.tsv")
  utils::write.table(row, path, sep = "\t", row.names = FALSE,
    col.names = !file.exists(path), quote = FALSE, append = file.exists(path))
  cat(sprintf(
    "[%s] phase=%s completed=%d/%d running=%d failed=%d pending=%d (%.1f%%); ETA=%s; SLURM remaining=%s\n",
    snapshot$heartbeat_utc, snapshot$phase, snapshot$completed, snapshot$total,
    snapshot$running, snapshot$failed, snapshot$pending,
    snapshot$percent_completed,
    if (is.null(snapshot$eta_seconds)) "estimating" else
      paste0(round(snapshot$eta_seconds), "s"),
    if (is.null(snapshot$slurm_remaining_seconds)) "unavailable" else
      paste0(round(snapshot$slurm_remaining_seconds), "s")
  ))
  flush.console()
  invisible(snapshot)
}

allb_worker <- function(task_index, spec, data, e, out) {
  task <- spec$tasks[task_index, , drop = FALSE]
  started <- proc.time()[["elapsed"]]
  tryCatch({
    payload <- allb_run_task(e, data, task, spec$configuration)
    checks <- allb_validate_payload(payload, task, spec$configuration)
    allb_assert(all(checks$passed), paste(
      "Task checks failed:", paste(checks$check[!checks$passed], collapse = ", ")
    ))
    shard <- list(
      schema_version = "all_binary_shard_v1",
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

allb_run_pending <- function(spec, data, e, out, cores, interval) {
  valid <- vapply(seq_len(nrow(spec$tasks)), function(i) {
    allb_valid_shard(file.path(out, "shards", spec$tasks$shard_file[i]),
      spec$tasks[i, , drop = FALSE], spec$scientific_signature,
      spec$configuration)
  }, logical(1))
  pending <- which(!valid)
  status <- rep("pending", nrow(spec$tasks)); status[valid] <- "passed"
  cores <- min(as.integer(cores), length(pending))
  started <- proc.time()[["elapsed"]]
  invocation <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "_", Sys.getpid())
  if (!length(pending)) {
    allb_write_progress(allb_progress_snapshot(spec, out, status, started,
      "validating_existing_outputs", invocation, 0L), out)
    return(status)
  }
  jobs <- list(); outcomes <- list(); completed_now <- 0L
  next_task <- 1L; last_update <- -Inf
  launch <- function(index) {
    status[index] <<- "running"
    job <- parallel::mcparallel(allb_worker(index, spec, data, e, out),
      silent = TRUE)
    jobs[[as.character(job$pid)]] <<- list(process = job, task_index = index)
  }
  while (next_task <= length(pending) && length(jobs) < cores) {
    launch(pending[next_task]); next_task <- next_task + 1L
  }
  repeat {
    collected <- if (length(jobs)) parallel::mccollect(
      lapply(jobs, `[[`, "process"), wait = FALSE
    ) else NULL
    if (length(collected)) {
      for (pid in names(collected)) {
        info <- jobs[[pid]]; result <- collected[[pid]]
        outcomes[[as.character(info$task_index)]] <- result
        if (inherits(result, "try-error") || is.null(result$status)) {
          status[info$task_index] <- "failed"
        } else {
          status[info$task_index] <- result$status
          if (identical(result$status, "passed"))
            completed_now <- completed_now + 1L
          if (identical(result$status, "failed")) cat(
            "Task ", result$task_id, " failed: ", result$message, "\n", sep = ""
          )
        }
        jobs[[pid]] <- NULL
      }
      while (next_task <= length(pending) && length(jobs) < cores) {
        launch(pending[next_task]); next_task <- next_task + 1L
      }
    }
    elapsed <- proc.time()[["elapsed"]] - started
    if (elapsed - last_update >= interval || !length(jobs)) {
      phase <- if (!length(jobs) && next_task > length(pending)) {
        if (any(status == "failed")) "failed" else "fitting_complete"
      } else "fitting"
      allb_write_progress(allb_progress_snapshot(
        spec, out, status, started, phase, invocation, completed_now
      ), out)
      last_update <- elapsed
    }
    if (!length(jobs) && next_task > length(pending)) break
    Sys.sleep(1)
  }
  allb_atomic_rds(list(
    invocation_id = invocation, status = status, outcomes = outcomes,
    scientific_signature = spec$scientific_signature
  ), file.path(out, "attempts", paste0(invocation, ".rds")))
  status
}

allb_final_names <- function() c(
  "results", "predictions", "selected_groups", "inner_candidates",
  "inner_folds", "selected_rows", "firth_refit_diagnostics"
)

allb_finalize <- function(e, spec, out) {
  paths <- file.path(out, "shards", spec$tasks$shard_file)
  valid <- vapply(seq_along(paths), function(i) allb_valid_shard(
    paths[i], spec$tasks[i, , drop = FALSE], spec$scientific_signature,
    spec$configuration
  ), logical(1))
  allb_assert(all(valid), paste("Valid shards:", sum(valid), "/", length(valid)))
  shards <- lapply(paths, readRDS)
  bind <- function(name) do.call(rbind,
    lapply(shards, function(x) x$payload[[name]]))
  final <- file.path(out, "final")
  dir.create(final, recursive = TRUE, showWarnings = FALSE)
  objects <- stats::setNames(lapply(allb_final_names(), bind), allb_final_names())
  for (name in names(objects)) {
    allb_atomic_rds(objects[[name]], file.path(final, paste0(name, ".rds")),
      identical_ok = TRUE)
    csv <- file.path(final, paste0(name, ".csv"))
    allb_assert(!file.exists(csv), paste("Refusing to overwrite:", csv))
    utils::write.csv(objects[[name]], csv, row.names = FALSE)
  }
  manifest_files <- c(file.path(out, "study_specification.rds"),
    paths, list.files(final, full.names = TRUE))
  manifest_files <- manifest_files[file.info(manifest_files)$isdir %in% FALSE]
  manifest <- data.frame(
    file = substring(manifest_files, nchar(out) + 2L),
    bytes = file.info(manifest_files)$size,
    sha256 = vapply(manifest_files, digest::digest, character(1),
      file = TRUE, algo = "sha256"), stringsAsFactors = FALSE
  )
  manifest <- manifest[order(manifest$file), ]
  manifest_path <- file.path(out, "MANIFEST.csv")
  allb_assert(!file.exists(manifest_path), "Refusing to overwrite MANIFEST.csv.")
  utils::write.csv(manifest, manifest_path, row.names = FALSE)
  marker <- list(
    version = spec$version, scientific_signature = spec$scientific_signature,
    completed_tasks = nrow(spec$tasks), selected_method_rows = nrow(objects$results),
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  allb_atomic_rds(marker,
    file.path(out, paste0("COMPLETED_", spec$version, ".rds")))
  invisible(marker)
}

allb_verify_outputs <- function(spec, out) {
  marker_path <- file.path(out, paste0("COMPLETED_", spec$version, ".rds"))
  manifest_path <- file.path(out, "MANIFEST.csv")
  allb_assert(file.exists(marker_path) && file.exists(manifest_path),
    "Completion marker or output manifest is missing.")
  marker <- readRDS(marker_path)
  allb_assert(identical(marker$scientific_signature, spec$scientific_signature) &&
    marker$completed_tasks == nrow(spec$tasks) &&
    marker$selected_method_rows == 6L * nrow(spec$tasks),
    "Completion counts or scientific signature do not match.")
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  final_files <- unlist(lapply(allb_final_names(), function(name) {
    file.path("final", paste0(name, c(".rds", ".csv")))
  }), use.names = FALSE)
  expected <- c("study_specification.rds",
    file.path("shards", spec$tasks$shard_file), final_files)
  allb_assert(identical(names(manifest), c("file", "bytes", "sha256")) &&
    !anyDuplicated(manifest$file) && setequal(manifest$file, expected),
    "ALL output manifest membership differs from the frozen schema.")
  paths <- file.path(out, manifest$file)
  allb_assert(all(file.exists(paths)) &&
    all(file.info(paths)$size == manifest$bytes) &&
    identical(unname(vapply(paths, digest::digest, character(1),
      file = TRUE, algo = "sha256")), manifest$sha256),
    "An ALL output size or SHA-256 checksum does not match.")
  valid <- vapply(seq_len(nrow(spec$tasks)), function(i) {
    allb_valid_shard(file.path(out, "shards", spec$tasks$shard_file[i]),
      spec$tasks[i, , drop = FALSE], spec$scientific_signature,
      spec$configuration)
  }, logical(1))
  allb_assert(all(valid), "One or more ALL shards fail identity/numerical checks.")
  shards <- lapply(file.path(out, "shards", spec$tasks$shard_file), readRDS)
  for (name in allb_final_names()) {
    combined <- do.call(rbind, lapply(shards, function(x) x$payload[[name]]))
    allb_assert(identical(combined,
      readRDS(file.path(out, "final", paste0(name, ".rds")))),
      paste("ALL final output/shard mismatch:", name))
  }
  invisible(TRUE)
}

allb_run <- function(root, stage, version, cores = 1L, interval = 30) {
  e <- allb_load_environment(root)
  if (stage == "production") {
    allb_assert(nzchar(Sys.getenv("SLURM_JOB_ID")),
      "Production is reserved for the user-submitted TRUBA job.")
    allb_assert(identical(Sys.getenv("SLURM_CPUS_PER_TASK"), "56"),
      "TRUBA Hamsi production requires exactly 56 allocated CPUs.")
    allb_validate_release(e, root, version)
  }
  proposed <- allb_make_spec(e, root, stage, version)
  out <- allb_output_directory(root, version)
  dir.create(file.path(out, "shards"), recursive = TRUE, showWarnings = FALSE)
  specification <- file.path(out, "study_specification.rds")
  if (file.exists(specification)) {
    spec <- readRDS(specification)
    allb_validate_spec(e, spec, root, stage, version)
  } else {
    spec <- proposed
    allb_atomic_rds(spec, specification)
  }
  cat("Stage:", stage, "\nVersion:", version,
      "\nScientific signature:", spec$scientific_signature,
      "\nWork units:", nrow(spec$tasks),
      "\nRequested workers:", cores, "\n")
  print(spec$runtime)
  data <- allb_load_data(root, spec$configuration)
  status <- allb_run_pending(spec, data, e, out, cores, interval)
  if (any(status == "failed")) stop(
    "One or more ALL tasks failed; valid signature-matching shards are resumable.",
    call. = FALSE
  )
  marker_path <- file.path(out, paste0("COMPLETED_", spec$version, ".rds"))
  if (file.exists(marker_path)) {
    marker <- readRDS(marker_path)
    allb_assert(identical(marker$scientific_signature,
      spec$scientific_signature) && marker$completed_tasks == nrow(spec$tasks) &&
      file.exists(file.path(out, "MANIFEST.csv")),
      "Existing ALL completion evidence is invalid.")
  } else {
    marker <- allb_finalize(e, spec, out)
  }
  allb_verify_outputs(spec, out)
  progress <- jsonlite::read_json(file.path(out, "progress.json"),
    simplifyVector = TRUE)
  progress$phase <- "complete"
  progress$eta_seconds <- 0
  progress$eta_status <- "complete_verified"
  progress$heartbeat_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  progress$estimated_completion_utc <- progress$heartbeat_utc
  allb_write_progress(progress, out)
  cat("Complete: TRUE\nOutput directory:", normalizePath(out),
      "\nCompleted tasks:", marker$completed_tasks, "\n")
  invisible(list(spec = spec, output = out, marker = marker))
}
