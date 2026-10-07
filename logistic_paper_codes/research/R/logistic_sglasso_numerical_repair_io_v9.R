# Immutable seven-shard I/O and verification for the V9 numerical repair.

lsg_repair_source_files_v9 <- function() {
  unique(c(
    lsg_v8_io_source_files(),
    file.path("R", c(
      "logistic_sglasso_numerical_repair_v9.R",
      "logistic_sglasso_numerical_repair_diagnostic_v9.R",
      "logistic_sglasso_numerical_repair_io_v9.R"
    )),
    "src/logistic_sglasso_joint_solver_v9.cpp",
    "config/logistic_sglasso_v7_failure_cases_v8.csv",
    "LOGISTIC_SGLASSO_NUMERICAL_REPAIR_PROTOCOL_V9.md",
    "scripts/45_validate_logistic_sglasso_numerical_repair_v9.R",
    "scripts/46_run_logistic_sglasso_numerical_repair_v9.R",
    "truba/run_logistic_sglasso_numerical_repair_smoke_v9.slurm",
    "truba/run_logistic_sglasso_numerical_repair_diagnostic_v9.slurm"
  ))
}


lsg_repair_safe_version_v9 <- function(version) {
  if (length(version) != 1L || is.na(version) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9._-]{2,100}$", version)) {
    stop("Unsafe V9 output version.", call. = FALSE)
  }
  version
}


lsg_repair_output_v9 <- function(root, version) {
  file.path(normalizePath(root, mustWork = TRUE), "outputs", "study",
            lsg_repair_safe_version_v9(version))
}


lsg_repair_task_grid_v9 <- function(cases) {
  data.frame(
    diagnostic_id = seq_len(nrow(cases)),
    case_id = cases$case_id,
    task_id = cases$task_id,
    diagnostic_type = cases$diagnostic_type,
    shard_file = sprintf("repair_shard_%02d_%s.rds", seq_len(nrow(cases)),
                         cases$case_id),
    stringsAsFactors = FALSE
  )
}


lsg_repair_spec_signature_v9 <- function(spec) {
  lsg_hash_v7(spec[c(
    "schema_version", "version", "source_v7_version",
    "source_v7_signature", "source_v8_signature", "configuration",
    "cases", "tasks", "source_input_manifest", "code_manifest", "runtime"
  )])
}


lsg_repair_make_spec_v9 <- function(root, audit, version) {
  defaults <- lsg_repair_defaults_v9()
  spec <- list(
    schema_version = defaults$schema_version,
    version = lsg_repair_safe_version_v9(version),
    source_v7_version = defaults$source_v7_version,
    source_v7_signature = defaults$source_v7_signature,
    source_v8_signature = defaults$source_v8_signature,
    configuration = defaults$configuration,
    cases = audit$cases,
    tasks = lsg_repair_task_grid_v9(audit$cases),
    source_input_manifest = audit$input_manifest,
    code_manifest = lsg_inventory_v7(root, lsg_repair_source_files_v9()),
    runtime = lsg_runtime_v7(),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  spec$scientific_signature <- lsg_repair_spec_signature_v9(spec)
  spec
}


lsg_repair_validate_spec_v9 <- function(spec, root, audit = NULL,
                                        require_runtime = TRUE) {
  defaults <- lsg_repair_defaults_v9()
  exact_identity <- identical(spec$schema_version, defaults$schema_version) &&
    identical(spec$source_v7_version, defaults$source_v7_version) &&
    identical(spec$source_v7_signature, defaults$source_v7_signature) &&
    identical(spec$source_v8_signature, defaults$source_v8_signature) &&
    identical(spec$scientific_signature,
              lsg_repair_spec_signature_v9(spec))
  if (!isTRUE(exact_identity)) {
    stop("V9 numerical-repair specification identity/signature mismatch.",
         call. = FALSE)
  }
  current_code <- lsg_inventory_v7(root, lsg_repair_source_files_v9())
  if (!identical(spec$code_manifest, current_code) ||
      !lsg_portable_numeric_equal_v9(
        spec$configuration, defaults$configuration, 1e-14
      )) {
    stop("V9 numerical-repair code or configuration changed.",
         call. = FALSE)
  }
  if (!is.null(audit)) {
    if (!lsg_portable_numeric_equal_v9(spec$cases, audit$cases, 1e-14) ||
        !identical(spec$tasks, lsg_repair_task_grid_v9(audit$cases)) ||
        !identical(spec$source_input_manifest, audit$input_manifest)) {
      stop("V9 frozen cases or V7 source evidence changed.", call. = FALSE)
    }
  }
  if (isTRUE(require_runtime) && !identical(spec$runtime, lsg_runtime_v7())) {
    stop("V9 run/resume requires the exact recorded runtime.",
         call. = FALSE)
  }
  invisible(TRUE)
}


lsg_repair_shard_paths_v9 <- function(output, task) {
  relative <- file.path("shards", task$shard_file)
  c(shard = lsg_path_v7(output, relative),
    receipt = lsg_path_v7(output, paste0(relative, ".receipt.rds")))
}


lsg_repair_shard_receipt_v9 <- function(shard, path) {
  list(
    schema_version = "logistic_sglasso_numerical_repair_shard_receipt_v9",
    scientific_signature = shard$scientific_signature,
    diagnostic_id = shard$task$diagnostic_id,
    case_id = shard$task$case_id,
    task_id = shard$task$task_id,
    shard_bytes = as.numeric(file.info(path)$size),
    shard_sha256 = lsg_file_hash_v7(path)
  )
}


lsg_repair_publish_shard_v9 <- function(output, shard) {
  paths <- lsg_repair_shard_paths_v9(output, shard$task)
  lsg_atomic_v7(shard, output,
                file.path("shards", shard$task$shard_file), "rds")
  receipt <- lsg_repair_shard_receipt_v9(shard, paths[["shard"]])
  lsg_atomic_v7(receipt, output,
                file.path("shards", paste0(
                  shard$task$shard_file, ".receipt.rds"
                )), "rds")
  invisible(TRUE)
}


lsg_repair_read_shard_v9 <- function(output, spec, task) {
  paths <- lsg_repair_shard_paths_v9(output, task)
  if (!all(file.exists(paths)) || any(dir.exists(paths)) ||
      any(nzchar(Sys.readlink(paths)))) {
    stop("Missing or unsafe V9 shard/receipt pair: ", task$case_id,
         call. = FALSE)
  }
  shard <- readRDS(paths[["shard"]])
  receipt <- readRDS(paths[["receipt"]])
  expected_receipt <- lsg_repair_shard_receipt_v9(
    shard, paths[["shard"]]
  )
  if (!identical(shard$schema_version,
                 "logistic_sglasso_numerical_repair_shard_v9") ||
      !identical(shard$scientific_signature, spec$scientific_signature) ||
      !identical(shard$task, task) ||
      !identical(receipt, expected_receipt) ||
      !is.data.frame(shard$payload$hard_checks) ||
      anyDuplicated(shard$payload$hard_checks$check)) {
    stop("V9 shard identity, receipt, or hard-check schema failed: ",
         task$case_id, call. = FALSE)
  }
  shard
}


lsg_repair_collect_v9 <- function(output, spec) {
  tasks <- spec$tasks
  present <- logical(nrow(tasks))
  shards <- vector("list", nrow(tasks))
  for (i in seq_len(nrow(tasks))) {
    paths <- lsg_repair_shard_paths_v9(output, tasks[i, , drop = FALSE])
    if (xor(file.exists(paths[["shard"]]), file.exists(paths[["receipt"]]))) {
      stop("Asymmetric V9 shard/receipt publication: ", tasks$case_id[i],
           call. = FALSE)
    }
    present[i] <- all(file.exists(paths))
    if (present[i]) shards[[i]] <- lsg_repair_read_shard_v9(
      output, spec, tasks[i, , drop = FALSE]
    )
  }
  names(shards) <- tasks$case_id
  list(present = present, shards = shards)
}


lsg_repair_bind_v9 <- function(items) {
  items <- Filter(function(value) is.data.frame(value) && nrow(value), items)
  if (!length(items)) return(data.frame())
  columns <- unique(unlist(lapply(items, names), use.names = FALSE))
  items <- lapply(items, function(value) {
    for (name in setdiff(columns, names(value))) value[[name]] <- NA
    value[columns]
  })
  out <- do.call(rbind, items)
  rownames(out) <- NULL
  out
}


lsg_repair_final_tables_v9 <- function(shards) {
  sglasso <- Filter(function(shard) {
    identical(shard$payload$schema_version,
              "logistic_sglasso_repair_case_v9")
  }, shards)
  grlasso <- Filter(function(shard) {
    identical(shard$payload$schema_version,
              "logistic_group_lasso_safe_prefix_case_v9")
  }, shards)
  list(
    sglasso_finite_paths = lsg_repair_bind_v9(lapply(
      sglasso, function(shard) shard$payload$finite_path
    )),
    sglasso_candidates = lsg_repair_bind_v9(lapply(
      sglasso, function(shard) shard$payload$candidates
    )),
    sglasso_stability = lsg_repair_bind_v9(lapply(
      sglasso, function(shard) shard$payload$selected_stability
    )),
    sglasso_summary = lsg_repair_bind_v9(lapply(
      sglasso, function(shard) shard$payload$summary
    )),
    grlasso_points = lsg_repair_bind_v9(lapply(
      grlasso, function(shard) shard$payload$points
    )),
    grlasso_summary = lsg_repair_bind_v9(lapply(
      grlasso, function(shard) shard$payload$summary
    )),
    hard_checks = lsg_repair_bind_v9(lapply(
      shards, function(shard) shard$payload$hard_checks
    ))
  )
}


lsg_repair_integrity_checks_v9 <- function(tables, spec) {
  checks <- c(
    seven_frozen_cases_present =
      length(unique(tables$hard_checks$case_id)) == 7L,
    four_sglasso_paths_present =
      nrow(tables$sglasso_summary) == 4L &&
      length(unique(tables$sglasso_summary$case_id)) == 4L,
    exact_156_sglasso_finite_points =
      nrow(tables$sglasso_finite_paths) == 4L * 39L,
    four_sglasso_stability_rows =
      nrow(tables$sglasso_stability) == 4L,
    three_grlasso_paths_present =
      nrow(tables$grlasso_summary) == 3L &&
      length(unique(tables$grlasso_summary$case_id)) == 3L,
    all_hard_checks_pass = nrow(tables$hard_checks) > 0L &&
      all(tables$hard_checks$passed %in% TRUE),
    all_sglasso_repairs_have_eligible_candidates =
      all(tables$sglasso_summary$eligible_finite_points > 0L),
    no_competitive_invalid_sglasso_candidate =
      all(!tables$sglasso_summary$invalid_candidate_competitive),
    all_selected_sglasso_candidates_interior =
      all(tables$sglasso_summary$selected_finite_point_interior),
    all_selected_sglasso_candidates_start_stable =
      all(tables$sglasso_summary$stability_accepted),
    all_grlasso_safe_prefixes_accepted =
      all(tables$grlasso_summary$group_lasso_safe_prefix) &&
      all(tables$grlasso_summary$
          selected_strictly_above_usable_lower_boundary),
    case_inventory_exact = setequal(
      unique(tables$hard_checks$case_id), spec$cases$case_id
    )
  )
  data.frame(check = names(checks),
             passed = unname(vapply(checks, isTRUE, logical(1))),
             stringsAsFactors = FALSE)
}


lsg_repair_final_names_v9 <- function(version) {
  c(
    sglasso_finite_paths = "final/sglasso_finite_paths.csv",
    sglasso_candidates = "final/sglasso_candidates.csv",
    sglasso_stability = "final/sglasso_selected_stability.csv",
    sglasso_summary = "final/sglasso_repair_summary.csv",
    grlasso_points = "final/grlasso_safe_prefix_points.csv",
    grlasso_summary = "final/grlasso_safe_prefix_summary.csv",
    hard_checks = "final/hard_checks.csv",
    integrity_checks = "final/final_integrity_checks.csv",
    record = paste0("final/diagnostic_record_", version, ".rds"),
    acceptance = paste0("final/ACCEPTED_", version, ".rds")
  )
}


lsg_repair_inventory_absolute_v9 <- function(base, files) {
  files <- sort(unique(files), method = "radix")
  paths <- file.path(base, files)
  if (!all(file.exists(paths)) || any(dir.exists(paths)) ||
      any(nzchar(Sys.readlink(paths)))) {
    stop("A V9 manifest candidate is absent, a directory, or a symlink.",
         call. = FALSE)
  }
  data.frame(
    file = files,
    bytes = unname(as.numeric(file.info(paths)$size)),
    sha256 = unname(vapply(paths, lsg_file_hash_v7, character(1))),
    stringsAsFactors = FALSE
  )
}


lsg_repair_verify_inventory_v9 <- function(base, manifest) {
  if (!is.data.frame(manifest) ||
      !identical(base::names(manifest), c("file", "bytes", "sha256")) ||
      anyDuplicated(manifest$file) || anyNA(manifest$file) ||
      anyNA(manifest$bytes) || anyNA(manifest$sha256)) {
    stop("V9 output-manifest schema or unique file identity failed.",
         call. = FALSE)
  }
  expected <- lsg_repair_inventory_absolute_v9(base, manifest$file)
  exact_identity <- identical(
    as.character(manifest$file), as.character(expected$file)
  ) && identical(
    as.character(manifest$sha256), as.character(expected$sha256)
  ) && identical(
    as.numeric(manifest$bytes), as.numeric(expected$bytes)
  )
  if (!isTRUE(exact_identity)) {
    stop("V9 output-manifest bytes or SHA-256 verification failed.",
         call. = FALSE)
  }
  invisible(TRUE)
}


lsg_repair_finalize_v9 <- function(output, spec, shards) {
  tables <- lsg_repair_final_tables_v9(shards)
  integrity <- lsg_repair_integrity_checks_v9(tables, spec)
  if (!all(integrity$passed %in% TRUE)) {
    print(integrity[!integrity$passed, ], row.names = FALSE)
    stop("V9 final scientific/numerical gates failed.", call. = FALSE)
  }
  final_paths <- lsg_repair_final_names_v9(spec$version)
  for (name in setdiff(base::names(final_paths), c(
    "integrity_checks", "record", "acceptance"
  ))) {
    lsg_atomic_v7(tables[[name]], output, final_paths[[name]], "csv", TRUE)
  }
  lsg_atomic_v7(
    integrity, output, final_paths[["integrity_checks"]], "csv", TRUE
  )
  record <- list(
    schema_version = "logistic_sglasso_numerical_repair_record_v9",
    version = spec$version,
    scientific_signature = spec$scientific_signature,
    tables = tables,
    integrity_checks = integrity
  )
  lsg_atomic_v7(record, output, final_paths[["record"]], "rds", TRUE)
  acceptance <- list(
    schema_version = "logistic_sglasso_numerical_repair_acceptance_v9",
    version = spec$version,
    scientific_signature = spec$scientific_signature,
    completed_cases = nrow(spec$tasks),
    hard_checks = nrow(tables$hard_checks),
    integrity_checks = nrow(integrity),
    all_hard_gates_passed = all(tables$hard_checks$passed),
    all_integrity_gates_passed = all(integrity$passed),
    accepted = TRUE
  )
  lsg_atomic_v7(
    acceptance, output, final_paths[["acceptance"]], "rds", TRUE
  )

  manifest_name <- paste0("output_manifest_", spec$version, ".csv")
  marker_name <- paste0("COMPLETED_", spec$version, ".txt")
  files <- list.files(output, recursive = TRUE, all.files = TRUE,
                      include.dirs = FALSE, no.. = TRUE)
  files <- files[!startsWith(files, ".run_lock/") &
                   !files %in% c(manifest_name, marker_name)]
  manifest <- lsg_repair_inventory_absolute_v9(output, files)
  lsg_atomic_v7(manifest, output, manifest_name, "csv", TRUE)
  marker <- c(
    paste("version", spec$version),
    paste("scientific_signature", spec$scientific_signature),
    paste("output_manifest_sha256", lsg_file_hash_v7(
      file.path(output, manifest_name)
    )),
    paste("completed_cases", nrow(spec$tasks)),
    "accepted TRUE"
  )
  lsg_atomic_v7(marker, output, marker_name, "text", TRUE)
  invisible(list(tables = tables, integrity = integrity,
                 acceptance = acceptance, manifest = manifest))
}


lsg_repair_initialize_v9 <- function(root, source_root, version) {
  audit <- lsg_v8_io_audit_source(
    root, source_root = source_root, require_runtime = TRUE
  )
  defaults <- lsg_repair_defaults_v9()
  if (!identical(audit$specification$scientific_signature,
                 defaults$source_v7_signature)) {
    stop("The audited V7 source signature is not the frozen V9 source.",
         call. = FALSE)
  }
  output <- lsg_repair_output_v9(root, version)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  spec_path <- file.path(output, "study_specification.rds")
  spec <- if (file.exists(spec_path)) readRDS(spec_path) else {
    value <- lsg_repair_make_spec_v9(root, audit, version)
    lsg_atomic_v7(value, output, "study_specification.rds", "rds")
    value
  }
  lsg_repair_validate_spec_v9(spec, root, audit, require_runtime = TRUE)
  lsg_atomic_v7(spec$cases, output, "diagnostic_case_inventory.csv", "csv", TRUE)
  lsg_atomic_v7(spec$tasks, output, "diagnostic_task_grid.csv", "csv", TRUE)
  lsg_atomic_v7(spec$source_input_manifest, output,
                "source_input_manifest.csv", "csv", TRUE)
  lsg_atomic_v7(spec$code_manifest, output, "code_manifest.csv", "csv", TRUE)
  lsg_atomic_v7(spec$runtime, output, "runtime.rds", "rds", TRUE)
  list(output = output, audit = audit, specification = spec)
}


lsg_run_numerical_repair_v9 <- function(
    root,
    source_root = NULL,
    version = lsg_repair_defaults_v9()$default_version,
    cores = 1L,
    max_seconds = Inf,
    max_tasks = Inf
) {
  root <- normalizePath(root, mustWork = TRUE)
  version <- lsg_repair_safe_version_v9(version)
  integer_argument <- function(value, name, infinite = FALSE) {
    if (isTRUE(infinite) && identical(as.numeric(value), Inf)) return(Inf)
    numeric <- as.numeric(value)
    if (length(numeric) != 1L || is.na(numeric) || !is.finite(numeric) ||
        numeric < 1 || numeric != as.integer(numeric)) {
      stop(name, " must be a positive integer", if (infinite) " or Inf" else "",
           ".", call. = FALSE)
    }
    as.integer(numeric)
  }
  cores <- integer_argument(cores, "cores")
  max_tasks <- integer_argument(max_tasks, "max_tasks", TRUE)
  max_seconds <- as.numeric(max_seconds)
  if (length(max_seconds) != 1L || is.na(max_seconds) || max_seconds <= 0) {
    stop("max_seconds must be positive.", call. = FALSE)
  }
  if (nzchar(Sys.getenv("SLURM_CPUS_PER_TASK")) &&
      cores > as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))) {
    stop("V9 workers exceed allocated SLURM CPUs.", call. = FALSE)
  }
  if (.Platform$OS.type != "unix" && cores > 1L) {
    stop("Parallel V9 diagnostics require Unix fork workers.",
         call. = FALSE)
  }
  output <- lsg_repair_output_v9(root, version)
  complete <- file.path(output, paste0("COMPLETED_", version, ".txt"))
  if (file.exists(complete)) {
    verified <- lsg_verify_numerical_repair_v9(
      root, version = version, source_root = source_root,
      require_source = TRUE, quiet = TRUE
    )
    cat("Valid completed V9 repair diagnostic reused; no files changed.\n")
    return(invisible(c(verified, list(output = output, complete = TRUE))))
  }
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  lock <- file.path(output, ".run_lock")
  if (!dir.create(lock, showWarnings = FALSE)) {
    stop("V9 run lock exists; inspect its owner before any removal: ", lock,
         call. = FALSE)
  }
  owner <- list(pid = Sys.getpid(), host = Sys.info()[["nodename"]],
                slurm_job = Sys.getenv("SLURM_JOB_ID"),
                started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  on.exit(lsg_release_lock_v7(lock, owner), add = TRUE)
  lsg_atomic_v7(owner, output, ".run_lock/owner.rds", "rds")

  initialized <- lsg_repair_initialize_v9(root, source_root, version)
  audit <- initialized$audit
  spec <- initialized$specification
  collected <- lsg_repair_collect_v9(output, spec)
  pending <- which(!collected$present)
  if (length(pending)) {
    compile_lsg_core(quiet = TRUE)
    lsg_compile_joint_solver_v9(root, quiet = TRUE)
  }
  limit <- min(length(pending), as.integer(min(max_tasks, .Machine$integer.max)))
  queue <- if (limit) head(pending, limit) else integer()
  started <- proc.time()[["elapsed"]]
  run_one <- function(i) {
    task <- spec$tasks[i, , drop = FALSE]
    case <- spec$cases[spec$cases$case_id == task$case_id, , drop = FALSE]
    task_started <- proc.time()[["elapsed"]]
    tryCatch({
      payload <- lsg_dispatch_repair_case_v9(
        root, audit, case, spec$configuration
      )
      shard <- list(
        schema_version = "logistic_sglasso_numerical_repair_shard_v9",
        scientific_signature = spec$scientific_signature,
        task = task,
        payload = payload,
        runtime_seconds = proc.time()[["elapsed"]] - task_started
      )
      lsg_repair_publish_shard_v9(output, shard)
      list(diagnostic_id = task$diagnostic_id, case_id = task$case_id,
           status = if (all(payload$hard_checks$passed %in% TRUE))
             "passed" else "gate_failed",
           message = paste(payload$hard_checks$check[
             !(payload$hard_checks$passed %in% TRUE)
           ], collapse = ";"))
    }, error = function(condition) {
      list(diagnostic_id = task$diagnostic_id, case_id = task$case_id,
           status = "error", message = conditionMessage(condition))
    })
  }
  results <- list()
  unstarted <- setdiff(pending, queue)
  while (length(queue) &&
         proc.time()[["elapsed"]] - started < max_seconds) {
    wave <- head(queue, min(cores, 7L))
    queue <- queue[-seq_along(wave)]
    outcome <- if (length(wave) > 1L && cores > 1L) {
      parallel::mclapply(
        wave, run_one, mc.cores = min(cores, length(wave)),
        mc.preschedule = FALSE, mc.set.seed = FALSE
      )
    } else {
      lapply(wave, run_one)
    }
    results <- c(results, outcome)
  }
  unstarted <- c(unstarted, queue)
  attempt <- list(
    schema_version = "logistic_sglasso_numerical_repair_attempt_v9",
    scientific_signature = spec$scientific_signature,
    started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    requested_cores = cores, pending_before = pending,
    attempted = as.integer(setdiff(pending, unstarted)),
    unstarted = as.integer(unstarted), results = results,
    slurm_job = Sys.getenv("SLURM_JOB_ID")
  )
  attempt_name <- paste0(
    "attempts/attempt_", format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"),
    "_", Sys.getpid(), ".rds"
  )
  lsg_atomic_v7(attempt, output, attempt_name, "rds")

  collected <- lsg_repair_collect_v9(output, spec)
  complete_now <- all(collected$present)
  passed_now <- complete_now && all(vapply(
    collected$shards, function(shard) {
      all(shard$payload$hard_checks$passed %in% TRUE)
    }, logical(1)
  ))
  if (passed_now) {
    finalized <- lsg_repair_finalize_v9(output, spec, collected$shards)
    return(invisible(list(
      output = output, complete = TRUE, specification = spec,
      acceptance = finalized$acceptance
    )))
  }
  failures <- Filter(function(value) !identical(value$status, "passed"), results)
  if (length(failures)) print(failures)
  cat("V9 repair diagnostic incomplete or a numerical gate failed.\n")
  cat("Valid shards:", sum(collected$present), "/", nrow(spec$tasks), "\n")
  cat("Output directory:", output, "\n")
  stop("V9 repair diagnostic did not pass all frozen-case gates.",
       call. = FALSE)
}


lsg_repair_read_csv_v9 <- function(path, template) {
  value <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                           na.strings = "NA")
  if (!identical(names(value), names(template))) {
    stop("V9 CSV schema differs from the signed RDS template.",
         call. = FALSE)
  }
  for (name in names(template)) {
    reference <- template[[name]]
    if (is.integer(reference)) value[[name]] <- as.integer(value[[name]])
    else if (is.logical(reference)) value[[name]] <- as.logical(value[[name]])
    else if (is.numeric(reference)) value[[name]] <- as.numeric(value[[name]])
    else value[[name]] <- as.character(value[[name]])
  }
  value
}


lsg_verify_numerical_repair_v9 <- function(
    root,
    version = lsg_repair_defaults_v9()$default_version,
    source_root = NULL,
    require_source = FALSE,
    quiet = FALSE
) {
  root <- normalizePath(root, mustWork = TRUE)
  version <- lsg_repair_safe_version_v9(version)
  output <- lsg_repair_output_v9(root, version)
  spec_path <- file.path(output, "study_specification.rds")
  if (!file.exists(spec_path)) stop("Missing V9 study specification.", call. = FALSE)
  spec <- readRDS(spec_path)
  audit <- if (isTRUE(require_source)) lsg_v8_io_audit_source(
    root, source_root = source_root,
    require_runtime = identical(spec$runtime, lsg_runtime_v7())
  ) else NULL
  lsg_repair_validate_spec_v9(
    spec, root, audit = audit,
    require_runtime = isTRUE(require_source)
  )
  manifest_name <- paste0("output_manifest_", version, ".csv")
  marker_name <- paste0("COMPLETED_", version, ".txt")
  manifest_path <- file.path(output, manifest_name)
  marker_path <- file.path(output, marker_name)
  if (!file.exists(manifest_path) || !file.exists(marker_path)) {
    stop("V9 output manifest or completion marker is absent.", call. = FALSE)
  }
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  lsg_repair_verify_inventory_v9(output, manifest)
  actual <- sort(list.files(
    output, recursive = TRUE, all.files = TRUE, include.dirs = FALSE,
    no.. = TRUE
  ), method = "radix")
  expected <- sort(c(manifest$file, manifest_name, marker_name),
                   method = "radix")
  if (!identical(actual, expected)) {
    stop("V9 output contains missing or unexpected files.", call. = FALSE)
  }
  collected <- lsg_repair_collect_v9(output, spec)
  if (!all(collected$present)) stop("V9 shard set is incomplete.", call. = FALSE)
  tables <- lsg_repair_final_tables_v9(collected$shards)
  integrity <- lsg_repair_integrity_checks_v9(tables, spec)
  if (!all(integrity$passed %in% TRUE)) {
    stop("Recomputed V9 final integrity checks failed.", call. = FALSE)
  }
  final_paths <- lsg_repair_final_names_v9(version)
  record <- readRDS(file.path(output, final_paths[["record"]]))
  acceptance <- readRDS(file.path(output, final_paths[["acceptance"]]))
  table_names <- intersect(base::names(tables), base::names(final_paths))
  archived_tables <- stats::setNames(lapply(table_names, function(name) {
    lsg_repair_read_csv_v9(
      file.path(output, final_paths[[name]]), tables[[name]]
    )
  }), table_names)
  csv_tables_match <- all(vapply(table_names, function(name) {
    lsg_portable_numeric_equal_v9(
      archived_tables[[name]], tables[[name]], tolerance = 1e-12
    )
  }, logical(1)))
  archived_integrity <- lsg_repair_read_csv_v9(
    file.path(output, final_paths[["integrity_checks"]]), integrity
  )
  if (!csv_tables_match || !lsg_portable_numeric_equal_v9(
        archived_integrity, integrity, tolerance = 1e-14
      ) || !lsg_portable_numeric_equal_v9(
        record$tables, tables, tolerance = 1e-12
      ) || !identical(record$scientific_signature,
                      spec$scientific_signature) ||
      !identical(acceptance$scientific_signature,
                 spec$scientific_signature) ||
      !isTRUE(acceptance$accepted) ||
      !isTRUE(acceptance$all_hard_gates_passed) ||
      !isTRUE(acceptance$all_integrity_gates_passed)) {
    stop("V9 final record or acceptance replay failed.", call. = FALSE)
  }
  marker <- readLines(marker_path, warn = FALSE)
  expected_marker <- c(
    paste("version", version),
    paste("scientific_signature", spec$scientific_signature),
    paste("output_manifest_sha256", lsg_file_hash_v7(manifest_path)),
    paste("completed_cases", nrow(spec$tasks)),
    "accepted TRUE"
  )
  if (!identical(marker, expected_marker)) {
    stop("V9 completion marker is not bound to the output manifest.",
         call. = FALSE)
  }
  mode <- if (identical(spec$runtime, lsg_runtime_v7())) {
    "exact_runtime"
  } else "cross_runtime_artifact_audit"
  if (!isTRUE(quiet)) {
    cat("V9 numerical-repair diagnostic accepted: TRUE\n")
    cat("Verification mode:", mode, "\n")
    cat("Cases:", nrow(spec$tasks), "\n")
  }
  invisible(list(
    accepted = TRUE, verification_mode = mode,
    output = output, specification = spec,
    acceptance = acceptance, integrity_checks = integrity,
    tables = tables, manifest = manifest
  ))
}


lsg_repair_io_unit_checks_v9 <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  temporary <- tempfile("lsg_v9_io_")
  dir.create(temporary)
  on.exit(unlink(temporary, recursive = TRUE, force = TRUE), add = TRUE)
  file <- file.path(temporary, "immutable.txt")
  writeLines("immutable", file, useBytes = TRUE)
  manifest <- lsg_repair_inventory_absolute_v9(temporary, "immutable.txt")
  manifest_ok <- isTRUE(lsg_repair_verify_inventory_v9(temporary, manifest))
  changed <- manifest
  changed$sha256 <- paste0("0", substring(changed$sha256, 2L))
  rejects <- function(expr) inherits(try(force(expr), silent = TRUE), "try-error")
  corruption_rejected <- rejects(
    lsg_repair_verify_inventory_v9(temporary, changed)
  )
  cases <- lsg_read_failure_cases_v8(file.path(
    root, "config", "logistic_sglasso_v7_failure_cases_v8.csv"
  ))
  tasks <- lsg_repair_task_grid_v9(cases)
  mock_output <- file.path(temporary, "finalize")
  dir.create(mock_output)
  mock_spec <- list(
    version = "lsg_v9_mock_finalize",
    tasks = tasks,
    cases = cases,
    scientific_signature = paste(rep("9", 64L), collapse = "")
  )
  mock_shards <- lapply(seq_len(nrow(tasks)), function(i) {
    task <- tasks[i, , drop = FALSE]
    hard_checks <- data.frame(
      case_id = task$case_id,
      task_id = task$task_id,
      check = paste0("mock_gate_", i),
      passed = TRUE,
      stringsAsFactors = FALSE
    )
    if (i <= 4L) {
      payload <- list(
        schema_version = "logistic_sglasso_repair_case_v9",
        finite_path = data.frame(
          case_id = task$case_id,
          lambda_index = seq_len(39L),
          kkt = rep(1e-7, 39L),
          stringsAsFactors = FALSE
        ),
        candidates = data.frame(
          case_id = task$case_id,
          selected = TRUE,
          validation_log_loss = 0.5,
          stringsAsFactors = FALSE
        ),
        selected_stability = data.frame(
          case_id = task$case_id,
          accepted = TRUE,
          stringsAsFactors = FALSE
        ),
        summary = data.frame(
          case_id = task$case_id,
          eligible_finite_points = 39L,
          invalid_candidate_competitive = FALSE,
          selected_finite_point_interior = TRUE,
          stability_accepted = TRUE,
          stringsAsFactors = FALSE
        ),
        hard_checks = hard_checks
      )
    } else {
      payload <- list(
        schema_version = "logistic_group_lasso_safe_prefix_case_v9",
        points = data.frame(
          case_id = task$case_id,
          lambda_index = 1:3,
          validation_log_loss = c(0.6, 0.5, 0.7),
          stringsAsFactors = FALSE
        ),
        summary = data.frame(
          case_id = task$case_id,
          group_lasso_safe_prefix = TRUE,
          selected_strictly_above_usable_lower_boundary = TRUE,
          stringsAsFactors = FALSE
        ),
        hard_checks = hard_checks
      )
    }
    list(task = task, payload = payload)
  })
  finalization <- try(
    lsg_repair_finalize_v9(mock_output, mock_spec, mock_shards),
    silent = TRUE
  )
  finalization_round_trip <- !inherits(finalization, "try-error")
  if (finalization_round_trip) {
    final_manifest <- utils::read.csv(
      file.path(
        mock_output,
        paste0("output_manifest_", mock_spec$version, ".csv")
      ),
      stringsAsFactors = FALSE
    )
    final_acceptance <- readRDS(file.path(
      mock_output,
      "final",
      paste0("ACCEPTED_", mock_spec$version, ".rds")
    ))
    finalization_round_trip <-
      isTRUE(lsg_repair_verify_inventory_v9(mock_output, final_manifest)) &&
      isTRUE(final_acceptance$accepted) &&
      identical(final_acceptance$completed_cases, 7L)
  }
  checks <- c(
    frozen_seven_case_inventory = nrow(cases) == 7L &&
      identical(cases$task_id, c(153L, 143L, 158L, 125L, 126L, 134L, 150L)),
    unique_atomic_task_names = nrow(tasks) == 7L &&
      !anyDuplicated(tasks$case_id) && !anyDuplicated(tasks$shard_file),
    v9_lambda_grid_exact = identical(
      lsg_lambda_relative_grid_v9(),
      c(4096, 2048, 1024, 512, 256, 128, 64, 32, 16,
        exp(seq(log(8), log(0.05), length.out = 30L)))
    ),
    sha256_manifest_round_trip = manifest_ok,
    sha256_manifest_corruption_rejected = corruption_rejected,
    seven_shard_finalization_round_trip = finalization_round_trip,
    unsafe_version_rejected = rejects(lsg_repair_safe_version_v9("../bad")),
    portable_csv_last_bit_tolerance =
      lsg_portable_numeric_equal_v9(
        data.frame(value = 1, flag = TRUE, label = "a"),
        data.frame(value = 1 + 5e-17, flag = TRUE, label = "a"),
        1e-14
      ),
    material_numeric_change_rejected =
      !lsg_portable_numeric_equal_v9(1, 1 + 1e-6, 1e-14)
  )
  data.frame(check = names(checks),
             passed = unname(vapply(checks, isTRUE, logical(1))),
             stringsAsFactors = FALSE)
}
