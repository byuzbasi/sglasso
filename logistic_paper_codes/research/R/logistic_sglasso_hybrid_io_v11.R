# Checkpoint, manifest, and verification layer for the three-case V11 hybrid
# diagnostic. All artifacts are versioned; signature-mismatching shards are
# rejected rather than skipped or overwritten.

lsg_hybrid_source_files_v11 <- function() {
  if (!exists("lsg_source_files_v7", mode = "function", inherits = TRUE)) {
    stop("Source the frozen V7 workflow before building a V11 inventory.",
         call. = FALSE)
  }
  unique(c(
    lsg_source_files_v7(),
    "R/logistic_sglasso_grpreg_diagnostic_v8.R",
    "R/logistic_sglasso_v7_failure_diagnostics_v8.R",
    "R/logistic_sglasso_v7_failure_diagnostics_io_v8.R",
    "config/logistic_sglasso_v7_failure_cases_v8.csv",
    "LOGISTIC_SGLASSO_V7_FAILURE_DIAGNOSTIC_PROTOCOL_V8.md",
    "R/logistic_sglasso_numerical_repair_v9.R",
    "LOGISTIC_SGLASSO_NUMERICAL_REPAIR_PROTOCOL_V9.md",
    "R/logistic_sglasso_numerical_repair_v10.R",
    "LOGISTIC_SGLASSO_NUMERICAL_REPAIR_PROTOCOL_V10.md",
    "src/logistic_sglasso_joint_solver_v9.cpp",
    "src/logistic_sglasso_block_kernel_v11.hpp",
    "src/logistic_sglasso_hybrid_solver_v11.cpp",
    "R/logistic_sglasso_hybrid_solver_v11.R",
    "R/logistic_sglasso_hybrid_diagnostic_v11.R",
    "R/logistic_sglasso_hybrid_io_v11.R",
    "R/logistic_sglasso_hybrid_integration_v11.R",
    "LOGISTIC_SGLASSO_HYBRID_SOLVER_PROTOCOL_V11.md",
    "scripts/49_validate_logistic_sglasso_hybrid_v11.R",
    "scripts/50_run_logistic_sglasso_hybrid_diagnostic_v11.R",
    "truba/README_TRUBA_LOGISTIC_SGLASSO_V11.md",
    "truba/run_logistic_sglasso_hybrid_smoke_v11.slurm",
    "truba/run_logistic_sglasso_hybrid_diagnostic_v11.slurm",
    "truba/build_logistic_sglasso_hybrid_bundle_v11.sh"
  ))
}


lsg_hybrid_safe_version_v11 <- function(version) {
  if (length(version) != 1L || is.na(version) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", version)) {
    stop("V11 version must contain only letters, digits, _ or -.",
         call. = FALSE)
  }
  version
}


lsg_hybrid_output_v11 <- function(root, version) {
  file.path(
    normalizePath(root, mustWork = TRUE), "outputs", "study",
    lsg_hybrid_safe_version_v11(version)
  )
}


lsg_hybrid_runtime_v11 <- function() {
  required <- c(
    "Rcpp", "RcppArmadillo", "digest", "sglasso", "adelie", "grpreg",
    "logistf", "mltools"
  )
  missing <- required[!vapply(
    required, requireNamespace, logical(1), quietly = TRUE
  )]
  if (length(missing)) {
    stop("Missing installed package(s): ", paste(missing, collapse = ", "),
         ". No automatic installation is permitted.", call. = FALSE)
  }
  list(
    r_version = R.version.string,
    platform = R.version$platform,
    package_versions = stats::setNames(vapply(required, function(package) {
      as.character(utils::packageVersion(package))
    }, character(1)), required),
    library_paths = .libPaths(),
    rng_kind = c("Mersenne-Twister", "Inversion", "Rejection")
  )
}


lsg_hybrid_inventory_v11 <- function(root, files) {
  files <- sort(unique(files), method = "radix")
  paths <- file.path(root, files)
  if (!all(file.exists(paths)) || any(dir.exists(paths)) ||
      any(nzchar(Sys.readlink(paths)))) {
    stop("A required regular V11 source file is absent or a symlink.",
         call. = FALSE)
  }
  data.frame(
    file = files,
    bytes = unname(as.numeric(file.info(paths)$size)),
    sha256 = unname(vapply(paths, lsg_file_hash_v7, character(1))),
    stringsAsFactors = FALSE
  )
}


lsg_hybrid_task_grid_v11 <- function(cases) {
  data.frame(
    diagnostic_task_id = seq_len(nrow(cases)),
    case_id = cases$case_id,
    source_task_id = as.integer(cases$task_id),
    scenario = cases$scenario,
    replication = as.integer(cases$replication),
    seed = as.integer(cases$seed),
    shard_file = sprintf(
      "hybrid_shard_%02d_%s.rds", seq_len(nrow(cases)), cases$case_id
    ),
    stringsAsFactors = FALSE
  )
}


# Resolve and validate the complete dispatch graph before compilation/forking.
# V7 calls its table `design`; V11 retains `scenarios` only internally.
lsg_hybrid_resolve_tasks_v11 <- function(spec) {
  required <- list(
    cases = c("case_id", "task_id", "scenario", "scenario_index",
              "replication", "seed", "alpha", "d"),
    tasks = c("diagnostic_task_id", "case_id", "source_task_id", "scenario",
              "replication", "seed", "shard_file"),
    scenarios = c("scenario", "scenario_index", "n_train", "n_validation",
                  "n_test", "groups", "group_size", "active_group_count",
                  "rho_within", "rho_between", "prevalence",
                  "linear_predictor_sd", "signal_pattern")
  )
  for (field in names(required)) {
    value <- spec[[field]]
    if (!is.data.frame(value) || nrow(value) == 0L ||
        !all(required[[field]] %in% names(value)) ||
        anyNA(value[required[[field]]])) {
      stop("Invalid V11 dispatch table: ", field,
           ". Expected a nonempty data frame with complete required columns.",
           call. = FALSE)
    }
  }
  if (anyDuplicated(spec$scenarios$scenario) ||
      anyDuplicated(spec$scenarios$scenario_index) ||
      anyDuplicated(spec$cases$case_id) ||
      !identical(spec$tasks, lsg_hybrid_task_grid_v11(spec$cases)) ||
      !is.list(spec$blinded_sources) ||
      anyDuplicated(names(spec$blinded_sources))) {
    stop("V11 dispatch identities are duplicate or inconsistent.", call. = FALSE)
  }
  lapply(seq_len(nrow(spec$tasks)), function(index) {
    task <- spec$tasks[index, , drop = FALSE]
    case <- spec$cases[spec$cases$case_id == task$case_id, , drop = FALSE]
    scenario <- spec$scenarios[
      spec$scenarios$scenario == case$scenario, , drop = FALSE
    ]
    source <- spec$blinded_sources[[as.character(case$task_id)]]
    if (nrow(case) != 1L || nrow(scenario) != 1L ||
        !isTRUE(scenario$scenario_index == case$scenario_index) ||
        !is.list(source) ||
        !identical(names(source), c("data_sha256", "firth_target_original",
                                    "sglasso_tuning")) ||
        !identical(names(source$data_sha256), c("training", "validation")) ||
        !all(vapply(source$data_sha256, function(x) {
          is.character(x) && length(x) == 1L && !is.na(x) &&
            grepl("^[a-f0-9]{64}$", x)
        }, logical(1))) ||
        !is.numeric(source$firth_target_original) ||
        length(source$firth_target_original) !=
          scenario$groups * scenario$group_size ||
        any(!is.finite(source$firth_target_original))) {
      stop("V11 task cannot resolve its case, scenario, or blinded source: ",
           task$case_id, call. = FALSE)
    }
    list(task = task, case = case, scenario = scenario, source = source)
  })
}


lsg_hybrid_prepare_runtime_v11 <- function(root) {
  compile_lsg_core(rebuild = FALSE, quiet = TRUE)
  lsg_compile_hybrid_solver_v11(root, rebuild = FALSE, quiet = TRUE)
  required <- c("lsg_kkt_cpp", "lsg_objective_cpp", "lsg_lambda_start_cpp",
                "lsg_fit_one_hybrid_v11_cpp", "lsg_path_hybrid_v11_cpp",
                "lsg_v8_tuning_data", "lsg_v8_data_fingerprints")
  missing <- required[!vapply(
    required, exists, logical(1), envir = environment(lsg_run_hybrid_case_v11),
    inherits = TRUE, mode = "function"
  )]
  if (length(missing)) stop("V11 startup dependency missing: ",
                            paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}


lsg_hybrid_make_spec_v11 <- function(root, audit, version) {
  defaults <- lsg_hybrid_diagnostic_defaults_v11()
  if (!identical(audit$specification$scientific_signature,
                 defaults$source_v7_signature)) {
    stop("The V7 source signature is not the frozen V11 source.",
         call. = FALSE)
  }
  cases <- audit$cases[
    match(defaults$case_ids, audit$cases$case_id), , drop = FALSE
  ]
  if (!identical(cases$case_id, defaults$case_ids)) {
    stop("The three frozen V11 cases are missing or reordered.",
         call. = FALSE)
  }
  tasks <- lsg_hybrid_task_grid_v11(cases)
  source_manifest <- lsg_hybrid_inventory_v11(
    root, lsg_hybrid_source_files_v11()
  )
  runtime <- lsg_hybrid_runtime_v11()
  spec <- list(
    schema_version = defaults$schema_version,
    version = lsg_hybrid_safe_version_v11(version),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    source_v7_version = defaults$source_v7_version,
    source_v7_signature = defaults$source_v7_signature,
    source_v7_root = audit$source_root,
    configuration = defaults$configuration,
    cases = cases,
    tasks = tasks,
    scenarios = audit$specification$design,
    source_configuration = audit$specification$configuration,
    blinded_sources = audit$blinded_sources,
    source_input_manifest = audit$input_manifest,
    source_manifest = source_manifest,
    runtime = runtime
  )
  invisible(lsg_hybrid_resolve_tasks_v11(spec))
  spec$scientific_signature <- lsg_hash_v7(spec[c(
    "schema_version", "version", "source_v7_version",
    "source_v7_signature", "configuration", "cases", "tasks",
    "scenarios", "source_configuration", "blinded_sources",
    "source_input_manifest", "source_manifest", "runtime"
  )])
  spec
}


lsg_hybrid_validate_spec_v11 <- function(spec, root, audit = NULL) {
  invisible(lsg_hybrid_resolve_tasks_v11(spec))
  defaults <- lsg_hybrid_diagnostic_defaults_v11()
  current_sources <- lsg_hybrid_inventory_v11(
    root, lsg_hybrid_source_files_v11()
  )
  expected_signature <- lsg_hash_v7(spec[c(
    "schema_version", "version", "source_v7_version",
    "source_v7_signature", "configuration", "cases", "tasks",
    "scenarios", "source_configuration", "blinded_sources",
    "source_input_manifest", "source_manifest", "runtime"
  )])
  valid <- is.list(spec) &&
    identical(spec$schema_version, defaults$schema_version) &&
    identical(spec$source_v7_signature, defaults$source_v7_signature) &&
    identical(spec$configuration, defaults$configuration) &&
    identical(spec$cases$case_id, defaults$case_ids) &&
    nrow(spec$tasks) == 3L &&
    identical(spec$source_manifest, current_sources) &&
    identical(spec$runtime, lsg_hybrid_runtime_v11()) &&
    identical(spec$scientific_signature, expected_signature)
  if (!is.null(audit)) {
    valid <- valid &&
      identical(spec$source_v7_root, audit$source_root) &&
      identical(spec$source_input_manifest, audit$input_manifest) &&
      identical(spec$blinded_sources, audit$blinded_sources) &&
      identical(spec$scenarios, audit$specification$design) &&
      identical(spec$source_configuration, audit$specification$configuration)
  }
  if (!isTRUE(valid)) {
    stop("Stored V11 specification, runtime, sources, or signature changed.",
         call. = FALSE)
  }
  invisible(TRUE)
}


lsg_hybrid_initialize_v11 <- function(root, audit, version) {
  output <- lsg_hybrid_output_v11(root, version)
  specification_path <- file.path(output, "study_specification.rds")
  if (dir.exists(output)) {
    if (nzchar(Sys.readlink(output)) || !file.exists(specification_path)) {
      stop("Existing V11 output is a symlink or lacks its specification.",
           call. = FALSE)
    }
    spec <- readRDS(specification_path)
    lsg_hybrid_validate_spec_v11(spec, root, audit)
    return(list(output = output, specification = spec))
  }
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  lock <- paste0(output, ".initialize.lock")
  if (!dir.create(lock, showWarnings = FALSE)) {
    stop("A concurrent V11 initialization lock exists.", call. = FALSE)
  }
  on.exit(unlink(lock, recursive = TRUE, force = FALSE), add = TRUE)
  spec <- lsg_hybrid_make_spec_v11(root, audit, version)
  stage <- tempfile(".v11_initialize_", tmpdir = dirname(output))
  dir.create(stage, showWarnings = FALSE)
  on.exit(if (dir.exists(stage)) unlink(stage, recursive = TRUE), add = TRUE)
  dir.create(file.path(stage, "shards"), recursive = TRUE)
  saveRDS(spec, file.path(stage, "study_specification.rds"), version = 3)
  utils::write.csv(
    spec$tasks, file.path(stage, "task_grid.csv"), row.names = FALSE
  )
  utils::write.csv(
    spec$source_manifest, file.path(stage, "source_manifest.csv"),
    row.names = FALSE
  )
  if (!file.rename(stage, output)) {
    stop("Unable to atomically publish the V11 output directory.",
         call. = FALSE)
  }
  list(output = output, specification = spec)
}


lsg_hybrid_shard_paths_v11 <- function(output, task) {
  shard <- file.path(output, "shards", task$shard_file)
  list(shard = shard, receipt = paste0(shard, ".receipt.rds"))
}


lsg_hybrid_expected_checks_v11 <- function() {
  c("test_fields_forbidden", "frozen_training_validation_fingerprints_exact",
    "exact_frozen_63_point_grid", "all_finite_outputs",
    "all_finite_points_numerically_eligible",
    "objective_reconstruction_within_32_ulp", "kkt_reconstructed",
    "intercept_profiled", "prediction_reconstructed",
    "selected_candidate_numerically_eligible",
    "no_competitive_invalid_candidate", "selected_solution_three_start_stable")
}


lsg_hybrid_payload_columns_v11 <- function() {
  path <- c("point_type", "lambda_index", "lambda_ratio", "lambda", "alpha", "d",
    "validation_log_loss", "objective", "objective_reconstruction_error",
    "objective_reconstruction_ulp_units", "objective_reconstruction_pass",
    "kkt", "kkt_reconstruction_error", "intercept_kkt", "group_kkt",
    "converged", "solver_route", "fallback_used", "block_stalled",
    "block_termination", "termination_reason", "block_sweeps", "apg_iterations",
    "selected_groups", "reconstruction_error", "numerically_eligible")
  list(
    finite_path = path, candidates = c(path, "selected"),
    selected_stability = c("selected_lambda", "selected_lambda_ratio",
      "selected_lambda_index", "forward_kkt_before_polish", "forward_kkt",
      "cold_kkt", "reverse_kkt", "forward_converged", "cold_converged",
      "reverse_converged", "forward_solver_route", "cold_solver_route",
      "reverse_solver_route", "maximum_objective_error",
      "maximum_validation_loss_error", "maximum_validation_probability_error",
      "forward_polish_objective_shift", "forward_polish_validation_loss_shift",
      "forward_polish_probability_shift", "polished_validation_log_loss", "accepted"),
    summary = c("alpha", "d", "finite_points", "eligible_finite_points",
      "fallback_finite_points", "stalled_block_points", "selected_point_type",
      "selected_lambda_index", "selected_lambda_ratio", "selected_validation_log_loss",
      "polished_selected_validation_log_loss", "selected_lower_lambda_boundary",
      "invalid_candidate_competitive", "maximum_finite_kkt",
      "maximum_objective_reconstruction_ulp_units", "maximum_reconstruction_error",
      "stability_accepted"),
    hard_checks = c("check", "passed"),
    scientific_outcomes = c("outcome", "value", "interpretation")
  )
}


lsg_hybrid_validate_shard_v11 <- function(shard, spec, task) {
  identity <- is.list(shard) &&
    identical(shard$schema_version,
              "logistic_sglasso_hybrid_diagnostic_shard_v11") &&
    identical(shard$scientific_signature, spec$scientific_signature) &&
    identical(shard$task, task) &&
    is.list(shard$payload) &&
    identical(shard$payload$schema_version,
              "logistic_sglasso_hybrid_case_v11")
  if (!isTRUE(identity)) return(FALSE)
  payload <- shard$payload
  counts <- c(finite_path = 63L, candidates = 64L, selected_stability = 1L,
              summary = 1L, hard_checks = 12L, scientific_outcomes = 4L)
  columns <- lsg_hybrid_payload_columns_v11()
  tables_valid <- all(vapply(names(counts), function(field) {
    value <- payload[[field]]
    is.data.frame(value) && nrow(value) == counts[[field]] &&
      all(c("case_id", "task_id", columns[[field]]) %in% names(value)) &&
      all(value$case_id %in% task$case_id) &&
      all(value$task_id %in% task$source_task_id)
  }, logical(1)))
  if (!tables_valid || !is.data.frame(payload$case) ||
      nrow(payload$case) != 1L) return(FALSE)
  case_fields <- c("case_id", "scenario", "replication", "seed")
  if (!all(case_fields %in% names(payload$case)) ||
      !identical(lapply(payload$case[case_fields], unname),
                 lapply(task[case_fields], unname))) return(FALSE)
  source <- spec$blinded_sources[[as.character(task$source_task_id)]]
  if (!is.list(source) || !identical(payload$data_sha256, source$data_sha256) ||
      !all(vapply(payload$data_sha256, function(x) {
        is.character(x) && length(x) == 1L && !is.na(x) &&
          grepl("^[a-f0-9]{64}$", x)
      }, logical(1)))) return(FALSE)
  checks <- payload$hard_checks
  is.logical(checks$passed) &&
    identical(checks$check, lsg_hybrid_expected_checks_v11()) &&
    identical(as.integer(payload$finite_path$lambda_index), seq_len(63L)) &&
    identical(as.character(payload$candidates$point_type),
              c("penalty_limit", rep("finite", 63L))) &&
    is.logical(payload$candidates$selected) &&
    !anyNA(payload$candidates$selected) &&
    sum(payload$candidates$selected) == 1L &&
    identical(payload$evaluation_data_used, FALSE) &&
    identical(names(payload$data_sha256), c("training", "validation"))
}


lsg_hybrid_read_shard_v11 <- function(output, spec, task) {
  paths <- lsg_hybrid_shard_paths_v11(output, task)
  if (!file.exists(paths$shard)) return(NULL)
  if (dir.exists(paths$shard) || nzchar(Sys.readlink(paths$shard))) {
    stop("A V11 shard path is not a regular file.", call. = FALSE)
  }
  shard <- readRDS(paths$shard)
  if (!lsg_hybrid_validate_shard_v11(shard, spec, task)) {
    stop("A stored V11 shard has the wrong identity or signature.",
         call. = FALSE)
  }
  expected_hash <- lsg_file_hash_v7(paths$shard)
  if (file.exists(paths$receipt)) {
    receipt <- readRDS(paths$receipt)
    if (!is.list(receipt) ||
        !identical(receipt$schema_version,
                   "logistic_sglasso_hybrid_receipt_v11") ||
        !identical(receipt$scientific_signature, spec$scientific_signature) ||
        !identical(receipt$shard_sha256, expected_hash) ||
        !identical(receipt$task_id, as.integer(task$diagnostic_task_id)) ||
        !identical(receipt$case_id, as.character(task$case_id))) {
      stop("A stored V11 shard receipt is invalid.", call. = FALSE)
    }
  } else {
    receipt <- list(
      schema_version = "logistic_sglasso_hybrid_receipt_v11",
      scientific_signature = spec$scientific_signature,
      task_id = as.integer(task$diagnostic_task_id),
      case_id = as.character(task$case_id),
      shard_sha256 = expected_hash
    )
    temporary <- tempfile(".receipt_", tmpdir = dirname(paths$receipt))
    saveRDS(receipt, temporary, version = 3)
    if (!file.rename(temporary, paths$receipt)) {
      stop("Unable to publish a recovered V11 receipt.", call. = FALSE)
    }
  }
  shard
}


lsg_hybrid_publish_shard_v11 <- function(output, spec, task, payload) {
  paths <- lsg_hybrid_shard_paths_v11(output, task)
  existing <- lsg_hybrid_read_shard_v11(output, spec, task)
  if (!is.null(existing)) return(existing)
  shard <- list(
    schema_version = "logistic_sglasso_hybrid_diagnostic_shard_v11",
    scientific_signature = spec$scientific_signature,
    task = task,
    payload = payload
  )
  if (!isTRUE(lsg_hybrid_validate_shard_v11(shard, spec, task))) {
    stop("Refusing to publish a malformed V11 task payload: ",
         task$case_id, call. = FALSE)
  }
  temporary <- tempfile(".hybrid_shard_", tmpdir = dirname(paths$shard))
  saveRDS(shard, temporary, version = 3)
  if (!file.rename(temporary, paths$shard)) {
    stop("Unable to atomically publish a V11 shard.", call. = FALSE)
  }
  lsg_hybrid_read_shard_v11(output, spec, task)
}


lsg_hybrid_bind_v11 <- function(shards, field) {
  values <- lapply(shards, function(shard) shard$payload[[field]])
  answer <- do.call(rbind, values)
  rownames(answer) <- NULL
  answer
}


lsg_hybrid_final_tables_v11 <- function(shards) {
  list(
    finite_paths = lsg_hybrid_bind_v11(shards, "finite_path"),
    candidates = lsg_hybrid_bind_v11(shards, "candidates"),
    selected_stability = lsg_hybrid_bind_v11(shards, "selected_stability"),
    case_summary = lsg_hybrid_bind_v11(shards, "summary"),
    hard_checks = lsg_hybrid_bind_v11(shards, "hard_checks"),
    scientific_outcomes = lsg_hybrid_bind_v11(
      shards, "scientific_outcomes"
    )
  )
}


lsg_hybrid_finalize_v11 <- function(root, output, spec, shards) {
  final <- file.path(output, "final")
  if (dir.exists(final)) return(invisible(final))
  tables <- lsg_hybrid_final_tables_v11(shards)
  accepted <- nrow(tables$hard_checks) > 0L &&
    all(tables$hard_checks$passed %in% TRUE) &&
    !anyDuplicated(tables$hard_checks[c("case_id", "check")]) &&
    nrow(tables$case_summary) == 3L &&
    nrow(tables$finite_paths) == 3L * 63L
  stage <- tempfile(".v11_final_", tmpdir = output)
  dir.create(stage, showWarnings = FALSE)
  on.exit(if (dir.exists(stage)) unlink(stage, recursive = TRUE), add = TRUE)
  file_names <- paste0(names(tables), ".csv")
  for (index in seq_along(tables)) {
    utils::write.csv(
      tables[[index]], file.path(stage, file_names[index]), row.names = FALSE
    )
  }
  saveRDS(
    list(
      schema_version = "logistic_sglasso_hybrid_acceptance_v11",
      version = spec$version,
      scientific_signature = spec$scientific_signature,
      accepted = accepted,
      completed_tasks = length(shards),
      failed_checks = tables$hard_checks[!tables$hard_checks$passed, ]
    ),
    file.path(stage, paste0("DIAGNOSTIC_ACCEPTANCE_", spec$version, ".rds")),
    version = 3
  )
  inventory_files <- c(
    file_names, paste0("DIAGNOSTIC_ACCEPTANCE_", spec$version, ".rds")
  )
  manifest <- data.frame(
    file = inventory_files,
    bytes = as.numeric(file.info(file.path(stage, inventory_files))$size),
    sha256 = vapply(
      file.path(stage, inventory_files), lsg_file_hash_v7, character(1)
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(manifest, file.path(stage, "output_manifest.csv"),
                   row.names = FALSE)
  completion <- c(
    paste("version", spec$version),
    paste("scientific_signature", spec$scientific_signature),
    paste("completed_tasks", length(shards)),
    paste("accepted", accepted),
    paste("manifest_sha256",
          lsg_file_hash_v7(file.path(stage, "output_manifest.csv")))
  )
  writeLines(
    completion,
    file.path(stage, paste0("DIAGNOSTIC_COMPLETE_", spec$version, ".txt")),
    useBytes = TRUE
  )
  if (accepted) {
    saveRDS(
      list(version = spec$version,
           scientific_signature = spec$scientific_signature,
           accepted = TRUE),
      file.path(stage, paste0("DIAGNOSTIC_ACCEPTED_", spec$version, ".rds")),
      version = 3
    )
  }
  if (!file.rename(stage, final)) {
    stop("Unable to atomically publish V11 final diagnostics.",
         call. = FALSE)
  }
  invisible(final)
}


lsg_run_hybrid_diagnostic_v11 <- function(
    root,
    source_root,
    version = lsg_hybrid_diagnostic_defaults_v11()$default_version,
    cores = 1L,
    max_seconds = 6600
) {
  root <- normalizePath(root, mustWork = TRUE)
  cores <- as.integer(cores)
  if (length(cores) != 1L || is.na(cores) || cores < 1L ||
      length(max_seconds) != 1L || is.na(max_seconds) ||
      !is.finite(max_seconds) || max_seconds <= 0) {
    stop("Invalid V11 run controls.", call. = FALSE)
  }
  audit <- lsg_v8_io_audit_source(
    root, source_root = source_root, require_runtime = TRUE
  )
  initialized <- lsg_hybrid_initialize_v11(root, audit, version)
  lsg_hybrid_execute_v11(root, initialized, cores, max_seconds)
}


# Shared by the real dispatcher and the explicitly synthetic integration test.
# Only the outer real-run entry point accepts/audits a frozen V7 source.
lsg_hybrid_execute_v11 <- function(root, initialized, cores, max_seconds,
                                   task_limit = Inf, data_inputs = NULL) {
  output <- initialized$output
  spec <- initialized$specification
  resolved <- lsg_hybrid_resolve_tasks_v11(spec)
  if (!is.null(data_inputs)) {
    if (!is.list(data_inputs) || !identical(names(data_inputs), spec$cases$case_id)) {
      stop("Local replay inputs do not match the frozen case identities.")
    }
    for (index in seq_along(resolved)) {
      data <- data_inputs[[index]]
      if (!identical(names(data), c("X_train", "y_train", "X_validation", "y_validation", "group")) ||
          !identical(lsg_v8_data_fingerprints(data), resolved[[index]]$source$data_sha256)) {
        stop("Local replay input fields/fingerprints differ: ", resolved[[index]]$case$case_id)
      }
    }
  }
  lock <- file.path(output, ".run_lock")
  if (!dir.create(lock, showWarnings = FALSE)) {
    stop("V11 run lock exists; inspect its owner before resuming: ", lock,
         call. = FALSE)
  }
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  saveRDS(list(pid = Sys.getpid(), host = Sys.info()[["nodename"]],
               slurm_job = Sys.getenv("SLURM_JOB_ID")),
          file.path(lock, "owner.rds"))
  start_time <- proc.time()[["elapsed"]]
  pending <- which(vapply(seq_len(nrow(spec$tasks)), function(index) {
    is.null(lsg_hybrid_read_shard_v11(
      output, spec, spec$tasks[index, , drop = FALSE]
    ))
  }, logical(1)))
  if (length(pending)) lsg_hybrid_prepare_runtime_v11(root)
  pending <- head(pending, min(length(pending), task_limit))
  run_one <- function(index) {
    if (proc.time()[["elapsed"]] - start_time > max_seconds) {
      return(list(index = index, status = "budget_deferred"))
    }
    item <- resolved[[index]]
    tryCatch({
      cat("Running V11 case:", item$case$case_id, "\n")
      set.seed(as.integer(item$case$seed), kind = "Mersenne-Twister",
               normal.kind = "Inversion", sample.kind = "Rejection")
      payload <- lsg_run_hybrid_case_v11(
        item$case, item$scenario, item$source,
        spec$source_configuration, spec$configuration,
        data = if (is.null(data_inputs)) NULL else data_inputs[[index]]
      )
      lsg_hybrid_publish_shard_v11(output, spec, item$task, payload)
      list(index = index, case_id = item$case$case_id, status = "completed")
    }, error = function(e) {
      list(index = index, case_id = item$case$case_id, status = "error",
           message = conditionMessage(e))
    })
  }
  if (length(pending)) {
    workers <- min(cores, length(pending))
    results <- if (.Platform$OS.type == "unix" && workers > 1L) {
      parallel::mclapply(
        pending, run_one, mc.cores = workers, mc.preschedule = FALSE
      )
    } else {
      lapply(pending, run_one)
    }
    failed <- vapply(results, inherits, logical(1), what = "try-error")
    if (any(failed)) {
      messages <- vapply(results[failed], as.character, character(1))
      cat("V11 task error(s):\n", paste(messages, collapse = "\n"), "\n",
          file = stderr())
    }
    for (result in results[!failed]) {
      if (identical(result$status, "error")) {
        cat("V11 case", result$case_id, "failed:", result$message, "\n",
            file = stderr())
      }
    }
  }
  shards <- lapply(seq_len(nrow(spec$tasks)), function(index) {
    lsg_hybrid_read_shard_v11(
      output, spec, spec$tasks[index, , drop = FALSE]
    )
  })
  complete <- all(!vapply(shards, is.null, logical(1)))
  if (complete) lsg_hybrid_finalize_v11(root, output, spec, shards)
  list(
    output = output,
    version = spec$version,
    scientific_signature = spec$scientific_signature,
    completed_tasks = sum(!vapply(shards, is.null, logical(1))),
    total_tasks = nrow(spec$tasks),
    complete = complete
  )
}


lsg_verify_hybrid_diagnostic_v11 <- function(
    root,
    version,
    source_root = NULL,
    require_source = TRUE,
    quiet = FALSE
) {
  root <- normalizePath(root, mustWork = TRUE)
  output <- lsg_hybrid_output_v11(root, version)
  specification_path <- file.path(output, "study_specification.rds")
  if (!file.exists(specification_path)) {
    stop("Missing V11 study specification.", call. = FALSE)
  }
  spec <- readRDS(specification_path)
  audit <- if (isTRUE(require_source)) {
    lsg_v8_io_audit_source(
      root, source_root = source_root, require_runtime = TRUE
    )
  } else NULL
  lsg_hybrid_validate_spec_v11(spec, root, audit)
  shards <- lapply(seq_len(nrow(spec$tasks)), function(index) {
    shard <- lsg_hybrid_read_shard_v11(
      output, spec, spec$tasks[index, , drop = FALSE]
    )
    if (is.null(shard)) stop("A required V11 shard is absent.", call. = FALSE)
    shard
  })
  lsg_hybrid_verify_final_v11(output, spec, shards, quiet)
}


lsg_hybrid_verify_final_v11 <- function(output, spec, shards, quiet = FALSE) {
  final <- file.path(output, "final")
  manifest_path <- file.path(final, "output_manifest.csv")
  acceptance_path <- file.path(
    final, paste0("DIAGNOSTIC_ACCEPTANCE_", spec$version, ".rds")
  )
  completion_path <- file.path(
    final, paste0("DIAGNOSTIC_COMPLETE_", spec$version, ".txt")
  )
  if (!all(file.exists(c(manifest_path, acceptance_path, completion_path)))) {
    stop("V11 final diagnostics are incomplete.", call. = FALSE)
  }
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  tables <- lsg_hybrid_final_tables_v11(shards)
  expected_files <- c(paste0(names(tables), ".csv"), basename(acceptance_path))
  if (!identical(names(manifest), c("file", "bytes", "sha256")) ||
      !identical(manifest$file, expected_files) ||
      any(!file.exists(file.path(final, expected_files))) ||
      any(nzchar(Sys.readlink(file.path(final, expected_files))))) {
    stop("V11 final manifest has an invalid file inventory.", call. = FALSE)
  }
  for (field in names(tables)) {
    observed <- utils::read.csv(file.path(final, paste0(field, ".csv")),
                                stringsAsFactors = FALSE)
    if (!isTRUE(all.equal(tables[[field]], observed, tolerance = 1e-12,
                          check.attributes = FALSE))) {
      stop("V11 final table disagrees with validated shards: ", field,
           call. = FALSE)
    }
  }
  current <- data.frame(
    file = manifest$file,
    bytes = as.numeric(file.info(file.path(final, manifest$file))$size),
    sha256 = vapply(
      file.path(final, manifest$file), lsg_file_hash_v7, character(1)
    ),
    stringsAsFactors = FALSE
  )
  acceptance <- readRDS(acceptance_path)
  checks <- utils::read.csv(
    file.path(final, "hard_checks.csv"), stringsAsFactors = FALSE
  )
  completion <- readLines(completion_path, warn = FALSE)
  expected_completion <- c(
    paste("version", spec$version),
    paste("scientific_signature", spec$scientific_signature),
    "completed_tasks 3",
    paste("accepted", acceptance$accepted),
    paste("manifest_sha256", lsg_file_hash_v7(manifest_path))
  )
  accepted_marker <- file.path(
    final, paste0("DIAGNOSTIC_ACCEPTED_", spec$version, ".rds")
  )
  marker_valid <- if (isTRUE(acceptance$accepted)) {
    file.exists(accepted_marker) &&
      identical(readRDS(accepted_marker), list(
        version = spec$version,
        scientific_signature = spec$scientific_signature,
        accepted = TRUE
      ))
  } else !file.exists(accepted_marker)
  valid <- isTRUE(all.equal(
      manifest, current, tolerance = 0, check.attributes = FALSE
    )) &&
    identical(acceptance$scientific_signature, spec$scientific_signature) &&
    identical(acceptance$completed_tasks, 3L) &&
    identical(acceptance$accepted, all(checks$passed %in% TRUE)) &&
    identical(completion, expected_completion) && marker_valid
  if (!isTRUE(valid)) {
    stop("V11 final manifest, acceptance, or hard checks failed.",
         call. = FALSE)
  }
  if (!isTRUE(quiet)) {
    cat("V11 completed tasks: 3 / 3\n")
    cat("V11 numerical acceptance:", acceptance$accepted, "\n")
    cat("Final directory:", final, "\n")
  }
  list(
    accepted = acceptance$accepted, output = output, final = final,
    specification = spec, shards = shards, checks = checks
  )
}


lsg_hybrid_io_unit_checks_v11 <- function() {
  temporary_root <- tempfile("lsg_v11_io_")
  dir.create(file.path(temporary_root, "shards"), recursive = TRUE)
  on.exit(unlink(temporary_root, recursive = TRUE, force = FALSE), add = TRUE)
  task <- data.frame(
    diagnostic_task_id = 1L, case_id = "synthetic_case",
    source_task_id = 1L, scenario = "synthetic", replication = 1L,
    seed = 11L, shard_file = "hybrid_shard_01_synthetic_case.rds",
    stringsAsFactors = FALSE
  )
  payload <- lsg_hybrid_io_fixture_v11(task)
  spec <- list(scientific_signature = paste(rep("a", 64L), collapse = ""),
               blinded_sources = list("1" = list(data_sha256 = payload$data_sha256)))
  published <- lsg_hybrid_publish_shard_v11(
    temporary_root, spec, task, payload
  )
  resumed <- lsg_hybrid_read_shard_v11(temporary_root, spec, task)
  paths <- lsg_hybrid_shard_paths_v11(temporary_root, task)
  receipt <- readRDS(paths$receipt)
  wrong <- spec
  wrong$scientific_signature <- paste(rep("b", 64L), collapse = "")
  mismatch_rejected <- inherits(try(
    lsg_hybrid_read_shard_v11(temporary_root, wrong, task), silent = TRUE
  ), "try-error")
  unsafe_version_rejected <- inherits(try(
    lsg_hybrid_safe_version_v11("../unsafe"), silent = TRUE
  ), "try-error")
  checks <- c(
    atomic_shard_published = file.exists(paths$shard),
    checksum_receipt_published = file.exists(paths$receipt) &&
      identical(receipt$shard_sha256, lsg_file_hash_v7(paths$shard)),
    matching_shard_resumes = identical(published, resumed),
    signature_mismatch_rejected = mismatch_rejected,
    unsafe_version_rejected = unsafe_version_rejected
  )
  data.frame(
    check = names(checks),
    passed = unname(vapply(checks, isTRUE, logical(1))),
    stringsAsFactors = FALSE
  )
}
