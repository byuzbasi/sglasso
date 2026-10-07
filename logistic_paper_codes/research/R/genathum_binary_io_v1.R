# Checkpoint, progress, validation, and finalization for GenAtHum binary V1.

gab_runtime <- function() {
  packages <- c("Rcpp", "RcppArmadillo", "digest", "adelie", "grpreg",
    "logistf", "mltools", "jsonlite")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  gab_assert(!length(missing), paste(
    "Missing installed packages (no automatic installation):",
    paste(missing, collapse = ", ")
  ))
  list(
    r_version = R.version.string, platform = R.version$platform,
    package_versions = stats::setNames(vapply(packages, function(package) {
      as.character(utils::packageVersion(package))
    }, character(1)), packages),
    library_paths = .libPaths(),
    rng_kind = c("Mersenne-Twister", "Inversion", "Rejection")
  )
}

gab_hash <- function(e, object) e$lsg_hash_v7(object)

gab_inventory <- function(e, root, files) e$lsg_inventory_v7(root, files)

gab_atomic_rds <- function(object, path, identical_ok = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) {
    if (identical_ok && identical(readRDS(path), object)) return(invisible(path))
    stop("Refusing to overwrite existing artifact: ", path, call. = FALSE)
  }
  temporary <- tempfile(pattern = paste0(".", basename(path), "."),
    tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, version = 3)
  gab_assert(file.rename(temporary, path), paste("Atomic rename failed:", path))
  invisible(path)
}

gab_atomic_json <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".progress_", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  jsonlite::write_json(object, temporary, auto_unbox = TRUE, pretty = FALSE,
    null = "null", na = "null")
  gab_assert(file.rename(temporary, path), "Atomic progress.json rename failed.")
  invisible(path)
}

gab_make_spec <- function(e, project_root, root, stage, version) {
  configuration <- gab_configuration(e, root, stage)
  data <- gab_load_data(project_root, configuration)
  tasks <- gab_tasks(configuration)
  source_manifest <- gab_inventory(e, root,
    unique(c(e$lsg_source_files_v7(), gab_source_files())))
  identity <- list(
    schema_version = "genathum_binary_study_v1", version = version,
    stage = stage, configuration = configuration, tasks = tasks,
    data_identity = list(
      source_path_basename = basename(data$source_path),
      source_sha256 = data$source_sha256,
      dimensions = dim(data$X),
      sample_id_sha256 = gab_hash(e, data$sample_id),
      sample_block_sha256 = gab_hash(e, data$sample_block),
      predictor_group_sha256 = gab_hash(e, data$group)
    ),
    source_manifest = source_manifest
  )
  identity$scientific_signature <- gab_hash(e, identity)
  c(identity, list(runtime = gab_runtime(),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)))
}

gab_validate_spec <- function(e, spec, project_root, root, stage, version) {
  current <- gab_make_spec(e, project_root, root, stage, version)
  fields <- c("schema_version", "version", "stage", "configuration", "tasks",
    "data_identity", "source_manifest", "scientific_signature")
  gab_assert(identical(spec[fields], current[fields]),
    "GenAtHum scientific identity, data, source, grid, or seeds changed.")
  gab_assert(identical(spec$runtime, current$runtime),
    "Runtime/package identity changed; resume under the original environment.")
  invisible(TRUE)
}

gab_output_directory <- function(root, version) {
  gab_assert(grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", version),
    "Unsafe GenAtHum version name.")
  file.path(root, "outputs", "study", version)
}

gab_validate_release <- function(e, project_root, root, version) {
  path <- file.path(root, "release", "genathum_binary_v3",
    "LOCAL_VALIDATED_genathum_binary_v3.rds")
  gab_assert(file.exists(path),
    "Missing local GenAtHum release receipt; production is blocked.")
  receipt <- readRDS(path)
  current <- gab_make_spec(e, project_root, root, "production", version)
  configuration_equal <- isTRUE(all.equal(
    receipt$production_configuration,
    current$configuration,
    tolerance = 1e-13,
    check.attributes = FALSE
  ))
  tasks_equal <- isTRUE(all.equal(
    receipt$production_tasks,
    current$tasks,
    tolerance = 0,
    check.attributes = FALSE
  ))
  gab_assert(identical(receipt$schema_version,
      "genathum_binary_portable_release_v3") &&
    isTRUE(receipt$accepted) && all(receipt$checks$passed) &&
    identical(receipt$production_version, version) &&
    identical(receipt$source_manifest, current$source_manifest) &&
    identical(receipt$data_identity, current$data_identity) &&
    configuration_equal && tasks_equal,
    paste(
      "Portable release receipt mismatch:",
      paste(c(
        if (!identical(receipt$source_manifest, current$source_manifest))
          "source_manifest",
        if (!identical(receipt$data_identity, current$data_identity))
          "data_identity",
        if (!configuration_equal) "production_configuration",
        if (!tasks_equal) "production_tasks"
      ), collapse = ", ")
    ))
  invisible(receipt)
}

gab_valid_shard <- function(path, task, signature, configuration) {
  if (!file.exists(path) || dir.exists(path)) return(FALSE)
  tryCatch({
    x <- readRDS(path)
    checks <- gab_validate_payload(x$payload, task, configuration)
    identical(x$schema_version, "genathum_binary_shard_v1") &&
      identical(x$scientific_signature, signature) &&
      identical(x$task, task) && all(checks$passed)
  }, error = function(error) FALSE)
}

gab_progress_snapshot <- function(spec, out, status, started, phase,
                                  invocation_id, completed_this_invocation = 0L) {
  valid <- vapply(seq_len(nrow(spec$tasks)), function(i) {
    gab_valid_shard(file.path(out, "shards", spec$tasks$shard_file[i]),
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
    schema_version = "genathum_binary_progress_v1",
    run_id = spec$version, scientific_signature = spec$scientific_signature,
    slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = "local"),
    invocation_id = invocation_id, phase = phase,
    work_unit = "endpoint_outer_replication_all_six_methods",
    total = nrow(spec$tasks), completed = completed, running = running,
    failed = failed, pending = pending,
    percent_completed = 100 * completed / nrow(spec$tasks),
    elapsed_seconds = elapsed,
    completed_this_invocation = completed_this_invocation,
    throughput_tasks_per_second = if (is.finite(rate)) rate else NULL,
    eta_seconds = if (is.finite(eta)) eta else NULL,
    estimated_completion_utc = if (is.finite(eta))
      format(now + eta, tz = "UTC", usetz = TRUE) else NULL,
    eta_status = if (is.finite(eta)) "measured_current_invocation" else "estimating_or_unknown",
    eta_scope = "task_fitting_only_excludes_final_aggregation",
    slurm_remaining_seconds = if (is.finite(remaining_slurm)) remaining_slurm else NULL,
    heartbeat_utc = format(now, tz = "UTC", usetz = TRUE),
    most_recently_completed_work_unit = if (any(valid))
      tail(spec$tasks$key[valid], 1L) else NULL
  )
}

gab_write_progress <- function(snapshot, out) {
  gab_atomic_json(snapshot, file.path(out, "progress.json"))
  row <- data.frame(
    timestamp_utc = snapshot$heartbeat_utc, phase = snapshot$phase,
    total = snapshot$total, completed = snapshot$completed,
    running = snapshot$running, failed = snapshot$failed,
    pending = snapshot$pending,
    percent_completed = snapshot$percent_completed,
    elapsed_seconds = snapshot$elapsed_seconds,
    throughput_tasks_per_second = snapshot$throughput_tasks_per_second %||% NA_real_,
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

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

gab_worker <- function(task_index, spec, data, e, out) {
  task <- spec$tasks[task_index, , drop = FALSE]
  outcome <- spec$configuration$outcomes[
    spec$configuration$outcomes$outcome_id == task$outcome_id, , drop = FALSE
  ]
  started <- proc.time()[["elapsed"]]
  tryCatch({
    payload <- gab_run_task(e, data, outcome, task, spec$configuration)
    checks <- gab_validate_payload(payload, task, spec$configuration)
    gab_assert(all(checks$passed), paste(
      "Task checks failed:", paste(checks$check[!checks$passed], collapse = ", ")
    ))
    shard <- list(
      schema_version = "genathum_binary_shard_v1",
      scientific_signature = spec$scientific_signature,
      task = task, checks = checks, payload = payload,
      runtime_seconds = proc.time()[["elapsed"]] - started,
      completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
    gab_atomic_rds(shard, file.path(out, "shards", task$shard_file))
    list(task_id = task$task_id[[1L]], status = "passed", message = "")
  }, error = function(error) list(
    task_id = task$task_id[[1L]], status = "failed",
    message = conditionMessage(error)
  ))
}

gab_run_pending <- function(spec, data, e, out, cores, interval) {
  valid <- vapply(seq_len(nrow(spec$tasks)), function(i) {
    gab_valid_shard(file.path(out, "shards", spec$tasks$shard_file[i]),
      spec$tasks[i, , drop = FALSE], spec$scientific_signature,
      spec$configuration)
  }, logical(1))
  pending <- which(!valid)
  status <- rep("pending", nrow(spec$tasks)); status[valid] <- "passed"
  cores <- min(as.integer(cores), length(pending))
  started <- proc.time()[["elapsed"]]
  invocation <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "_", Sys.getpid())
  if (!length(pending)) {
    gab_write_progress(gab_progress_snapshot(spec, out, status, started,
      "validating_existing_outputs", invocation, 0L), out)
    return(status)
  }
  jobs <- list(); completed_now <- 0L; next_task <- 1L; last_update <- -Inf
  outcomes <- list()
  launch <- function(index) {
    status[index] <<- "running"
    job <- parallel::mcparallel(gab_worker(index, spec, data, e, out),
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
          if (identical(result$status, "passed")) completed_now <- completed_now + 1L
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
      gab_write_progress(gab_progress_snapshot(
        spec, out, status, started, phase, invocation, completed_now
      ), out)
      last_update <- elapsed
    }
    if (!length(jobs) && next_task > length(pending)) break
    Sys.sleep(1)
  }
  gab_atomic_rds(list(invocation_id = invocation, status = status,
    outcomes = outcomes, scientific_signature = spec$scientific_signature),
    file.path(out, "attempts", paste0(invocation, ".rds")))
  status
}

gab_verify_outputs <- function(spec, out) {
  marker_path <- file.path(out, paste0("COMPLETED_", spec$version, ".rds"))
  manifest_path <- file.path(out, "MANIFEST.csv")
  gab_assert(file.exists(marker_path) && file.exists(manifest_path),
    "Completion marker or output manifest is missing.")
  marker <- readRDS(marker_path)
  gab_assert(identical(marker$scientific_signature, spec$scientific_signature) &&
    marker$completed_tasks == nrow(spec$tasks) &&
    marker$selected_method_rows == 6L * nrow(spec$tasks),
    "Completion counts or scientific signature do not match.")
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  expected <- c("study_specification.rds",
    file.path("shards", spec$tasks$shard_file),
    paste0("final/", rep(c("results", "predictions", "selected_groups",
      "inner_candidates", "inner_folds", "selected_rows",
      "firth_refit_diagnostics"), each = 2L), rep(c(".rds", ".csv"), 7L)))
  gab_assert(identical(names(manifest), c("file", "bytes", "sha256")) &&
    !anyDuplicated(manifest$file) && setequal(manifest$file, expected),
    "Output manifest membership differs from the frozen schema.")
  paths <- file.path(out, manifest$file)
  gab_assert(all(file.exists(paths)) &&
    all(file.info(paths)$size == manifest$bytes) &&
    identical(unname(vapply(paths, digest::digest, character(1),
      file = TRUE, algo = "sha256")), manifest$sha256),
    "An output size or SHA-256 checksum does not match.")
  valid <- vapply(seq_len(nrow(spec$tasks)), function(i) {
    gab_valid_shard(file.path(out, "shards", spec$tasks$shard_file[i]),
      spec$tasks[i, , drop = FALSE], spec$scientific_signature, spec$configuration)
  }, logical(1))
  gab_assert(all(valid), "One or more final shards fail identity/numerical checks.")
  shards <- lapply(file.path(out, "shards", spec$tasks$shard_file), readRDS)
  for (name in c("results", "predictions", "selected_groups", "inner_candidates",
                "inner_folds", "selected_rows", "firth_refit_diagnostics")) {
    combined <- do.call(rbind, lapply(shards, function(x) x$payload[[name]]))
    gab_assert(identical(combined, readRDS(file.path(out, "final",
      paste0(name, ".rds")))), paste("Final output/shard mismatch:", name))
  }
  invisible(TRUE)
}

gab_finalize <- function(e, spec, out) {
  paths <- file.path(out, "shards", spec$tasks$shard_file)
  valid <- vapply(seq_along(paths), function(i) gab_valid_shard(
    paths[i], spec$tasks[i, , drop = FALSE], spec$scientific_signature,
    spec$configuration
  ), logical(1))
  gab_assert(all(valid), paste("Valid shards:", sum(valid), "/", length(valid)))
  shards <- lapply(paths, readRDS)
  bind <- function(name) do.call(rbind, lapply(shards, function(x) x$payload[[name]]))
  final <- file.path(out, "final")
  dir.create(final, recursive = TRUE, showWarnings = FALSE)
  objects <- list(
    results = bind("results"), predictions = bind("predictions"),
    selected_groups = bind("selected_groups"),
    inner_candidates = bind("inner_candidates"),
    inner_folds = bind("inner_folds"),
    selected_rows = bind("selected_rows"),
    firth_refit_diagnostics = bind("firth_refit_diagnostics")
  )
  for (name in names(objects)) {
    gab_atomic_rds(objects[[name]], file.path(final, paste0(name, ".rds")),
      identical_ok = TRUE)
    csv <- file.path(final, paste0(name, ".csv"))
    if (!file.exists(csv)) utils::write.csv(objects[[name]], csv, row.names = FALSE)
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
  utils::write.csv(manifest, file.path(out, "MANIFEST.csv"), row.names = FALSE)
  marker <- list(
    version = spec$version, scientific_signature = spec$scientific_signature,
    completed_tasks = nrow(spec$tasks), selected_method_rows = nrow(objects$results),
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  gab_atomic_rds(marker, file.path(out, paste0("COMPLETED_", spec$version, ".rds")),
    identical_ok = TRUE)
  invisible(marker)
}

gab_run <- function(project_root, root, stage, version, cores = 1L,
                    interval = 30) {
  e <- gab_load_environment(root)
  if (stage == "production") {
    gab_assert(nzchar(Sys.getenv("SLURM_JOB_ID")),
      "Production is reserved for the user-submitted TRUBA job.")
    gab_assert(identical(Sys.getenv("SLURM_CPUS_PER_TASK"), "56"),
      "TRUBA Hamsi production requires exactly 56 allocated CPUs.")
    gab_validate_release(e, project_root, root, version)
  }
  proposed <- gab_make_spec(e, project_root, root, stage, version)
  out <- gab_output_directory(root, version)
  dir.create(file.path(out, "shards"), recursive = TRUE, showWarnings = FALSE)
  specification <- file.path(out, "study_specification.rds")
  if (file.exists(specification)) {
    spec <- readRDS(specification)
    gab_validate_spec(e, spec, project_root, root, stage, version)
  } else {
    spec <- proposed
    gab_atomic_rds(spec, specification)
  }
  cat("Stage:", stage, "\nVersion:", version,
      "\nScientific signature:", spec$scientific_signature,
      "\nWork units:", nrow(spec$tasks),
      "\nRequested workers:", cores, "\n")
  print(spec$runtime)
  data <- gab_load_data(project_root, spec$configuration)
  status <- gab_run_pending(spec, data, e, out, cores, interval)
  if (any(status == "failed")) stop(
    "One or more tasks failed; valid signature-matching shards are resumable.",
    call. = FALSE
  )
  marker_path <- file.path(out, paste0("COMPLETED_", spec$version, ".rds"))
  if (file.exists(marker_path)) {
    marker <- readRDS(marker_path)
    gab_assert(identical(marker$scientific_signature,
      spec$scientific_signature) && marker$completed_tasks == nrow(spec$tasks) &&
      file.exists(file.path(out, "MANIFEST.csv")),
      "Existing completion evidence is invalid.")
  } else {
    marker <- gab_finalize(e, spec, out)
  }
  gab_verify_outputs(spec, out)
  progress <- jsonlite::read_json(file.path(out, "progress.json"),
    simplifyVector = TRUE)
  progress$phase <- "complete"
  progress$eta_seconds <- 0
  progress$eta_status <- "complete_verified"
  progress$heartbeat_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  progress$estimated_completion_utc <- progress$heartbeat_utc
  gab_write_progress(progress, out)
  cat("Complete: TRUE\nOutput directory:", normalizePath(out),
      "\nCompleted tasks:", marker$completed_tasks, "\n")
  invisible(list(spec = spec, output = out, marker = marker))
}
