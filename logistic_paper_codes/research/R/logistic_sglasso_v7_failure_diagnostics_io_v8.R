# V8 diagnostic-only orchestration for the frozen V7 failure set.
#
# This module never reads test/truth artifacts into a diagnostic payload.  It
# audits the immutable V7 evidence, publishes seven atomic diagnostic shards,
# and keeps numerical hard gates separate from scientific outcomes.

lsg_v8_io_defaults <- function() {
  failures <- data.frame(
    task_id = c(125L, 126L, 134L, 141L, 143L, 144L, 146L, 147L,
                150L, 151L, 153L, 155L, 157L, 158L),
    status = c("gate_failed", "error", "error", rep("gate_failed", 5L),
               "error", rep("gate_failed", 5L)),
    message = c(
      "sglasso_finite_range_resolved",
      "No numerically eligible validation candidate was available for Logistic Group Lasso (grpreg).",
      "No numerically eligible validation candidate was available for Logistic Group Lasso (grpreg).",
      rep("all_finite_sglasso_paths_represented", 5L),
      "No numerically eligible validation candidate was available for Logistic Group Lasso (grpreg).",
      rep("all_finite_sglasso_paths_represented", 5L)
    ), stringsAsFactors = FALSE
  )
  list(
    schema_version = "logistic_sglasso_v7_failure_diagnostic_v8",
    shard_schema = "logistic_sglasso_v7_failure_diagnostic_shard_v8",
    receipt_schema = "logistic_sglasso_v7_failure_diagnostic_shard_receipt_v8",
    attempt_schema = "logistic_sglasso_v7_failure_diagnostic_attempt_v8",
    acceptance_schema = "logistic_sglasso_v7_failure_diagnostic_acceptance_v8",
    default_version = "logistic_sglasso_v7_failure_diagnostic_v8",
    source_version = "logistic_sglasso_prediction_selection_pilot_r20_v7",
    source_schema = "prediction_selection_v7",
    source_signature =
      "93ce8540e29856c38b2fa0ea144a038c68e1b80ceb4ea5adca3094c34f7850ec",
    source_attempt = "attempt_20260831T232447_741476.rds",
    source_task_count = 160L,
    source_completed_count = 146L,
    source_failed_shard_ids = c(125L, 141L, 143L, 144L, 146L, 147L,
                                151L, 153L, 155L, 157L, 158L),
    source_error_ids = c(126L, 134L, 150L),
    source_failures = failures,
    ranked_zero_path_count = 76L,
    case_ids = c("sg_near", "sg_lower_median", "sg_worst", "tail_125",
                 "grlasso_126", "grlasso_134", "grlasso_150"),
    case_task_ids = c(153L, 143L, 158L, 125L, 126L, 134L, 150L),
    tolerances = list(
      lambda_reference = 1e-10,
      replay_lambda = 1e-10,
      replay_loss = 1e-8,
      replay_kkt = 1e-8,
      replay_intercept = 1e-8,
      reconstruction = 2e-6,
      kkt_decomposition = 1e-10,
      profile_intercept = 1e-10,
      profile_objective_increase = 1e-12,
      tail_anchor_loss = 1e-5,
      tail_endpoint_loss = 1e-10
    )
  )
}

lsg_v8_io_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

lsg_v8_io_equal <- function(x, y, tolerance = 1e-10) {
  isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE))
}

lsg_v8_io_require_functions <- function(functions) {
  envir <- environment(lsg_v8_io_require_functions)
  missing <- functions[!vapply(functions, exists, logical(1), mode = "function",
                               envir = envir, inherits = TRUE)]
  lsg_v8_io_assert(!length(missing), paste(
    "Source the frozen V7 workflow and V8 numerical modules first; missing:",
    paste(missing, collapse = ", ")))
  invisible(TRUE)
}

lsg_v8_io_source_files <- function() {
  c(
    lsg_source_files_v7(),
    "R/logistic_sglasso_grpreg_diagnostic_v8.R",
    "R/logistic_sglasso_v7_failure_diagnostics_v8.R",
    "R/logistic_sglasso_v7_failure_diagnostics_io_v8.R",
    "src/logistic_sglasso_intercept_profile_v8.cpp",
    "config/logistic_sglasso_v7_failure_cases_v8.csv",
    "LOGISTIC_SGLASSO_V7_FAILURE_DIAGNOSTIC_PROTOCOL_V8.md",
    "scripts/43_validate_logistic_sglasso_v7_failure_diagnostics_v8.R",
    "scripts/44_run_logistic_sglasso_v7_failure_diagnostics_v8.R",
    "truba/run_logistic_sglasso_failure_smoke_v8.slurm",
    "truba/run_logistic_sglasso_failure_diagnostic_v8.slurm"
  )
}

lsg_v8_io_safe_version <- function(version) {
  lsg_v8_io_assert(length(version) == 1L && !is.na(version) &&
    grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", version),
    "V8 version must contain only letters, digits, _ or -.")
  version
}

lsg_v8_io_output <- function(root, version) {
  lsg_path_v7(root, file.path("outputs", "study", lsg_v8_io_safe_version(version)))
}

# source_root is the V7 study output directory, not the project/code root.
lsg_v8_io_source_root <- function(root, source_root = NULL) {
  defaults <- lsg_v8_io_defaults()
  path <- if (is.null(source_root) || !nzchar(source_root)) {
    file.path(root, "outputs", "study", defaults$source_version)
  } else source_root
  lsg_v8_io_assert(length(path) == 1L && dir.exists(path),
                   paste("Missing V7 source study directory:", path))
  normalizePath(path, mustWork = TRUE)
}

lsg_v8_io_source_path <- function(source_root, relative) {
  lsg_path_v7(source_root, relative)
}

lsg_v8_io_inventory_absolute <- function(base, files) {
  files <- sort(unique(files), method = "radix")
  paths <- vapply(files, function(file) lsg_path_v7(base, file), character(1))
  lsg_v8_io_assert(all(file.exists(paths)) && !any(dir.exists(paths)),
                   "A required V8 inventory file is absent.")
  data.frame(
    file = files,
    bytes = as.numeric(file.info(paths)$size),
    sha256 = unname(vapply(paths, lsg_file_hash_v7, character(1))),
    stringsAsFactors = FALSE
  )
}

lsg_v8_io_verify_inventory <- function(base, inventory) {
  lsg_v8_io_assert(is.data.frame(inventory) &&
    identical(names(inventory), c("file", "bytes", "sha256")) &&
    !anyNA(inventory) && !anyDuplicated(inventory$file),
    "Invalid V8 inventory schema.")
  current <- lsg_v8_io_inventory_absolute(base, inventory$file)
  lsg_v8_io_assert(lsg_v8_io_equal(current, inventory, 0),
                   "V8 inventory bytes or SHA-256 hashes changed.")
  invisible(TRUE)
}

lsg_v8_io_failure_table <- function(attempt) {
  if (!length(attempt$failures)) {
    return(data.frame(task_id = integer(), key = character(), status = character(),
                      message = character(), stringsAsFactors = FALSE))
  }
  rows <- lapply(attempt$failures, function(item) data.frame(
    task_id = as.integer(item$task_id), key = as.character(item$key),
    status = as.character(item$status), message = as.character(item$message),
    stringsAsFactors = FALSE
  ))
  out <- do.call(rbind, rows)
  out[order(out$task_id), , drop = FALSE]
}

lsg_v8_io_read_source_shard <- function(source_root, spec, task) {
  relative <- file.path("shards", task$shard_file)
  path <- lsg_v8_io_source_path(source_root, relative)
  receipt_path <- lsg_v8_io_source_path(source_root,
                                        paste0(relative, ".receipt.rds"))
  lsg_v8_io_assert(file.exists(path) && file.exists(receipt_path), paste(
    "Missing V7 failed shard/receipt pair for", task$key))
  shard <- readRDS(path)
  receipt <- readRDS(receipt_path)
  lsg_v8_io_assert(
    identical(shard$schema_version, "prediction_selection_shard_v7") &&
      identical(shard$scientific_signature, spec$scientific_signature) &&
      identical(shard$task, task) &&
      identical(receipt, lsg_shard_receipt_v7(shard, path)),
    paste("V7 failed shard identity/receipt mismatch for", task$key)
  )
  lsg_v8_io_assert(is.data.frame(shard$checks) &&
    identical(names(shard$checks), c("check", "passed")) &&
    !anyNA(shard$checks$check) && !anyDuplicated(shard$checks$check),
    paste("Invalid archived V7 check schema for", task$key))
  shard
}

lsg_v8_io_zero_path_ranking <- function(shards) {
  rows <- list()
  for (shard in shards) {
    tuning <- shard$payload$sglasso_tuning
    lsg_v8_io_assert(is.data.frame(tuning) && all(c(
      "point_type", "alpha_index", "d_index", "alpha", "d", "kkt",
      "numerically_eligible"
    ) %in% names(tuning)), "Archived V7 SGLASSO tuning schema changed.")
    tuning <- tuning[tuning$point_type == "finite", , drop = FALSE]
    key <- interaction(tuning$alpha_index, tuning$d_index, drop = TRUE,
                       lex.order = TRUE)
    pieces <- split(tuning, key)
    for (piece in pieces) {
      if (!any(piece$numerically_eligible %in% TRUE)) {
        lsg_v8_io_assert(all(is.finite(piece$kkt)),
                         "A ranked V7 finite path has nonfinite KKT evidence.")
        rows[[length(rows) + 1L]] <- data.frame(
          task_id = as.integer(shard$task$task_id),
          alpha_index = as.integer(piece$alpha_index[[1L]]),
          d_index = as.integer(piece$d_index[[1L]]),
          alpha = as.numeric(piece$alpha[[1L]]),
          d = as.numeric(piece$d[[1L]]),
          finite_points = nrow(piece),
          eligible_points = sum(piece$numerically_eligible %in% TRUE),
          minimum_kkt = min(piece$kkt),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  out <- out[order(out$minimum_kkt, out$task_id, out$alpha_index, out$d_index),
             , drop = FALSE]
  out$representative_rank <- seq_len(nrow(out))
  rownames(out) <- NULL
  out
}

lsg_v8_io_audit_source <- function(root, source_root = NULL,
                                   require_runtime = TRUE) {
  lsg_v8_io_require_functions(c(
    "lsg_source_files_v7", "lsg_inventory_v7", "lsg_spec_signature_v7",
    "lsg_configuration_v7", "lsg_design_for_stage_v7", "lsg_task_grid_v7",
    "lsg_shard_receipt_v7", "lsg_read_failure_cases_v8"
  ))
  defaults <- lsg_v8_io_defaults()
  source_root <- lsg_v8_io_source_root(root, source_root)
  spec <- readRDS(lsg_v8_io_source_path(source_root, "study_specification.rds"))
  lsg_v8_io_assert(
    identical(spec$schema_version, defaults$source_schema) &&
      identical(spec$version, defaults$source_version) &&
      identical(spec$stage, "pilot") &&
      identical(spec$scientific_signature, defaults$source_signature) &&
      identical(spec$scientific_signature, lsg_spec_signature_v7(spec)),
    "The V7 source study identity or scientific signature is not exact."
  )
  lsg_v8_io_assert(lsg_v8_io_equal(
    spec$source_manifest, lsg_inventory_v7(root, lsg_source_files_v7()), 0),
    "Frozen V7 source files no longer match the archived source manifest.")
  lsg_v8_io_assert(
    lsg_v8_io_equal(spec$configuration, lsg_configuration_v7("pilot"), 1e-13) &&
      lsg_v8_io_equal(spec$design, lsg_design_for_stage_v7(root, "pilot"), 1e-13) &&
      identical(spec$tasks, lsg_task_grid_v7(spec$design, spec$configuration)) &&
      nrow(spec$tasks) == defaults$source_task_count,
    "Frozen V7 design/configuration/160-task grid changed."
  )
  archived_tasks <- utils::read.csv(lsg_v8_io_source_path(
    source_root, "task_grid.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  archived_design <- utils::read.csv(lsg_v8_io_source_path(
    source_root, "design.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  archived_sources <- utils::read.csv(lsg_v8_io_source_path(
    source_root, "source_manifest.csv"), stringsAsFactors = FALSE,
    check.names = FALSE)
  lsg_v8_io_assert(
    lsg_v8_io_equal(archived_tasks, spec$tasks, 0) &&
      lsg_v8_io_equal(archived_design, spec$design, 1e-13) &&
      lsg_v8_io_equal(archived_sources, spec$source_manifest, 0),
    "Archived V7 task/design/source CSVs disagree with the signed specification."
  )
  if (isTRUE(require_runtime)) {
    lsg_v8_io_assert(identical(spec$runtime, lsg_runtime_v7()),
                     "Exact V7 diagnostic replay requires the archived runtime.")
  }

  attempt_path <- lsg_v8_io_source_path(
    source_root, file.path("attempts", defaults$source_attempt))
  lsg_v8_io_assert(file.exists(attempt_path),
                   "The frozen original V7 attempt record is absent.")
  attempt <- readRDS(attempt_path)
  failures <- lsg_v8_io_failure_table(attempt)
  expected <- defaults$source_failures
  expected$key <- spec$tasks$key[match(expected$task_id, spec$tasks$task_id)]
  expected <- expected[c("task_id", "key", "status", "message")]
  expected <- expected[order(expected$task_id), , drop = FALSE]
  rownames(expected) <- NULL
  rownames(failures) <- NULL
  lsg_v8_io_assert(
    identical(attempt$scientific_signature, defaults$source_signature) &&
      identical(attempt$runtime, spec$runtime) &&
      length(attempt$unstarted_tasks) == 0L &&
      length(attempt$completed_now) == defaults$source_completed_count &&
      setequal(as.integer(attempt$completed_now),
               setdiff(spec$tasks$task_id, expected$task_id)) &&
      identical(failures, expected),
    "The original V7 146-pass/11-gate/3-error attempt record changed."
  )

  shard_tasks <- spec$tasks[match(defaults$source_failed_shard_ids,
                                 spec$tasks$task_id), , drop = FALSE]
  shards <- setNames(lapply(seq_len(nrow(shard_tasks)), function(index) {
    shard <- lsg_v8_io_read_source_shard(source_root, spec,
                                         shard_tasks[index, , drop = FALSE])
    failed <- shard$checks$check[!shard$checks$passed]
    expected_message <- expected$message[expected$task_id == shard$task$task_id]
    lsg_v8_io_assert(identical(failed, expected_message), paste(
      "The archived failed hard gate changed for task", shard$task$task_id))
    shard
  }), as.character(shard_tasks$task_id))

  for (task_id in defaults$source_error_ids) {
    task <- spec$tasks[spec$tasks$task_id == task_id, , drop = FALSE]
    path <- lsg_v8_io_source_path(source_root,
                                  file.path("shards", task$shard_file))
    lsg_v8_io_assert(!file.exists(path) && !file.exists(paste0(path, ".receipt.rds")),
      paste("A V7 source shard unexpectedly exists for original error task", task_id))
  }

  ranking <- lsg_v8_io_zero_path_ranking(shards)
  lsg_v8_io_assert(nrow(ranking) == defaults$ranked_zero_path_count,
                   "The archived zero-eligible SGLASSO path count is not 76.")
  cases <- lsg_read_failure_cases_v8(file.path(
    root, "config", "logistic_sglasso_v7_failure_cases_v8.csv"))
  lsg_v8_io_assert(identical(cases$case_id, defaults$case_ids) &&
    identical(as.integer(cases$task_id), defaults$case_task_ids),
    "The seven V8 cases or their frozen order changed.")
  for (case_id in c("sg_near", "sg_lower_median", "sg_worst")) {
    case <- cases[cases$case_id == case_id, , drop = FALSE]
    ranked <- ranking[ranking$representative_rank == case$representative_rank,
                      , drop = FALSE]
    lsg_v8_io_assert(nrow(ranked) == 1L &&
      ranked$task_id == case$task_id &&
      abs(ranked$alpha - case$alpha) <= 1e-12 &&
      abs(ranked$d - case$d) <= 1e-12 &&
      lsg_v8_io_equal(ranked$minimum_kkt, case$archived_anchor_kkt, 1e-14),
      paste("The archived ranked SGLASSO path changed for", case_id))
  }
  for (index in seq_len(nrow(cases))) {
    case <- cases[index, , drop = FALSE]
    task <- spec$tasks[spec$tasks$task_id == case$task_id, , drop = FALSE]
    failure <- expected[expected$task_id == case$task_id, , drop = FALSE]
    lsg_v8_io_assert(nrow(task) == 1L && task$scenario == case$scenario &&
      task$replication == case$replication && task$seed == case$seed,
      paste("Frozen case/task identity mismatch for", case$case_id))
    lsg_v8_io_assert(nrow(failure) == 1L &&
      failure$status == case$expected_v7_status &&
      failure$message == case$expected_v7_message,
      paste("Frozen attempt outcome mismatch for", case$case_id))
    if (case$diagnostic_type != "grpreg_group_lasso_raw") {
      lsg_v8_io_assert(identical(as.character(case$source_shard_file),
                                 as.character(task$shard_file)), paste(
        "Frozen source shard filename mismatch for", case$case_id))
      tuning <- shards[[as.character(case$task_id)]]$payload$sglasso_tuning
      anchor <- tuning[
        tuning$point_type == "finite" &
          abs(tuning$alpha - case$alpha) <= 1e-12 &
          abs(tuning$d - case$d) <= 1e-12 &
          tuning$lambda_index == case$archived_anchor_lambda_index,
        , drop = FALSE
      ]
      lsg_v8_io_assert(nrow(anchor) == 1L &&
        abs(anchor$lambda_relative_to_reference -
              case$archived_anchor_lambda_relative) <= 1e-12 &&
        abs(anchor$lambda_reference - case$archived_lambda_reference) <= 1e-12 &&
        abs(anchor$kkt - case$archived_anchor_kkt) <= 1e-14 &&
        abs(anchor$intercept_score -
              case$archived_anchor_intercept_score) <= 1e-14 &&
        abs(anchor$validation_log_loss -
              case$archived_anchor_validation_log_loss) <= 1e-13,
        paste("Frozen archived tuning anchor mismatch for", case$case_id))
    } else {
      lsg_v8_io_assert(is.na(case$source_shard_file), paste(
        "An original-error Group Lasso case must not claim a V7 shard:",
        case$case_id))
    }
  }

  input_files <- c(
    "study_specification.rds", "task_grid.csv", "design.csv", "source_manifest.csv",
    file.path("attempts", defaults$source_attempt),
    as.vector(rbind(
      file.path("shards", shard_tasks$shard_file),
      paste0(file.path("shards", shard_tasks$shard_file), ".receipt.rds")
    ))
  )
  blinded_sources <- lapply(shards, lsg_v8_io_blind_source)
  names(blinded_sources) <- names(shards)
  rm(shards)
  list(
    source_root = source_root, specification = spec, attempt = attempt,
    failure_table = failures, blinded_sources = blinded_sources,
    zero_path_ranking = ranking, cases = cases,
    input_manifest = lsg_v8_io_inventory_absolute(source_root, input_files)
  )
}

lsg_v8_io_blind_source <- function(shard) {
  lsg_v8_io_assert(is.list(shard) && is.list(shard$payload) &&
    is.list(shard$payload$artifacts) &&
    is.list(shard$payload$artifacts$data_sha256),
    "A verified V7 source shard is required for blinding.")
  hashes <- shard$payload$artifacts$data_sha256
  allowed_tuning <- c(
    "candidate_id", "alpha_index", "d_index", "alpha", "d",
    "lambda_index", "lambda", "lambda_reference", "lambda_reference_type",
    "lambda_relative_to_reference", "point_type", "validation_log_loss",
    "kkt", "intercept_score", "passes", "converged",
    "numerically_eligible", "selected_free_d"
  )
  lsg_v8_io_assert(is.data.frame(shard$payload$sglasso_tuning) &&
    all(allowed_tuning %in% names(shard$payload$sglasso_tuning)),
    "The archived V7 tuning frame lacks a required blinded column.")
  blinded <- list(
    data_sha256 = list(training = hashes$training,
                       validation = hashes$validation),
    firth_target_original = shard$payload$artifacts$firth_target_original,
    sglasso_tuning = shard$payload$sglasso_tuning[allowed_tuning]
  )
  lsg_v8_io_assert(identical(names(blinded), c(
    "data_sha256", "firth_target_original", "sglasso_tuning")) &&
    identical(names(blinded$data_sha256), c("training", "validation")) &&
    all(vapply(blinded$data_sha256, function(hash) {
      is.character(hash) && length(hash) == 1L && grepl("^[a-f0-9]{64}$", hash)
    }, logical(1))), "Blinded source schema is invalid.")
  blinded
}

lsg_v8_io_task_grid <- function(cases) {
  slugs <- gsub("[^A-Za-z0-9_-]", "_", cases$case_id)
  data.frame(
    diagnostic_task_id = seq_len(nrow(cases)),
    case_id = cases$case_id,
    diagnostic_type = cases$diagnostic_type,
    source_task_id = as.integer(cases$task_id),
    source_key = paste(cases$scenario, cases$replication, sep = "::"),
    shard_file = sprintf("diagnostic_shard_%02d_%s.rds",
                         seq_len(nrow(cases)), slugs),
    stringsAsFactors = FALSE
  )
}

lsg_v8_io_spec_signature <- function(spec) {
  lsg_hash_v7(spec[c(
    "schema_version", "version", "source_version", "source_schema",
    "source_scientific_signature", "source_runtime", "configuration",
    "design", "cases", "tasks", "source_input_manifest", "code_manifest",
    "runtime", "tolerances"
  )])
}

lsg_v8_io_make_spec <- function(root, audit, version) {
  defaults <- lsg_v8_io_defaults()
  spec <- list(
    schema_version = defaults$schema_version,
    version = lsg_v8_io_safe_version(version),
    source_version = defaults$source_version,
    source_schema = defaults$source_schema,
    source_scientific_signature = defaults$source_signature,
    source_runtime = audit$specification$runtime,
    configuration = audit$specification$configuration,
    design = audit$specification$design,
    cases = audit$cases,
    tasks = lsg_v8_io_task_grid(audit$cases),
    source_input_manifest = audit$input_manifest,
    code_manifest = lsg_inventory_v7(root, lsg_v8_io_source_files()),
    runtime = lsg_runtime_v7(),
    tolerances = defaults$tolerances,
    source_location_at_creation = audit$source_root,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  spec$scientific_signature <- lsg_v8_io_spec_signature(spec)
  spec
}

lsg_v8_io_validate_spec <- function(spec, root, audit,
                                    require_runtime = TRUE) {
  defaults <- lsg_v8_io_defaults()
  lsg_v8_io_assert(
    identical(spec$schema_version, defaults$schema_version) &&
      identical(spec$source_version, defaults$source_version) &&
      identical(spec$source_schema, defaults$source_schema) &&
      identical(spec$source_scientific_signature, defaults$source_signature) &&
      identical(spec$scientific_signature, lsg_v8_io_spec_signature(spec)),
    "V8 diagnostic specification identity/signature mismatch."
  )
  lsg_v8_io_assert(
    identical(spec$source_runtime, audit$specification$runtime) &&
      lsg_v8_io_equal(spec$configuration,
                      audit$specification$configuration, 1e-13) &&
      lsg_v8_io_equal(spec$design, audit$specification$design, 1e-13) &&
      identical(spec$cases, audit$cases) &&
      identical(spec$tasks, lsg_v8_io_task_grid(audit$cases)) &&
      identical(spec$source_input_manifest, audit$input_manifest) &&
      identical(spec$code_manifest,
                lsg_inventory_v7(root, lsg_v8_io_source_files())),
    "V8 frozen inputs, cases, code, or source evidence changed."
  )
  if (isTRUE(require_runtime)) {
    lsg_v8_io_assert(identical(spec$runtime, lsg_runtime_v7()),
                     "Run/resume requires the exact V8 runtime.")
  }
  invisible(TRUE)
}

lsg_v8_io_check_rows <- function(case, checks) {
  data.frame(
    case_id = as.character(case$case_id),
    task_id = as.integer(case$task_id),
    check = names(checks),
    passed = unname(vapply(checks, isTRUE, logical(1))),
    stringsAsFactors = FALSE
  )
}

lsg_v8_io_value <- function(frame, name, default = NA) {
  if (is.data.frame(frame) && nrow(frame) == 1L && name %in% names(frame)) {
    frame[[name]][[1L]]
  } else default
}

lsg_v8_io_hard_checks_intercept <- function(payload, case, configuration,
                                             tolerances) {
  comparison <- payload$comparison
  summary <- payload$summary
  expected_ratios <- configuration$lambda_relative_grid
  checks <- list(
    payload_schema = identical(payload$schema_version,
                               "sglasso_intercept_failure_diagnostic_v8"),
    case_identity = is.data.frame(payload$case) && nrow(payload$case) == 1L &&
      identical(as.character(payload$case$case_id), as.character(case$case_id)),
    blinded_no_test = identical(payload$test_used, FALSE) &&
      identical(names(payload$data_sha256), c("training", "validation")),
    unique_comparison_columns = is.data.frame(comparison) &&
      !anyDuplicated(names(comparison)),
    exact_36_point_path = is.data.frame(comparison) &&
      nrow(comparison) == 36L &&
      lsg_v8_io_equal(comparison$lambda_ratio, expected_ratios, 1e-12),
    finite_path_outputs = is.data.frame(comparison) && all(is.finite(c(
      comparison$lambda, comparison$baseline_intercept,
      comparison$baseline_training_log_loss,
      comparison$baseline_validation_log_loss, comparison$baseline_objective,
      comparison$baseline_kkt, comparison$baseline_intercept_kkt,
      comparison$baseline_group_kkt, comparison$profiled_intercept,
      comparison$profiled_training_log_loss,
      comparison$profiled_validation_log_loss, comparison$profiled_objective,
      comparison$profiled_kkt, comparison$profiled_intercept_kkt,
      comparison$profiled_group_kkt
    ))),
    finite_coefficients_and_probabilities = is.data.frame(comparison) &&
      all(comparison$baseline_beta_finite) &&
      all(comparison$baseline_probability_finite) &&
      all(comparison$profiled_probability_finite),
    profiled_beta_unchanged = is.data.frame(comparison) &&
      all(comparison$profile_beta_unchanged),
    lambda_reference_replay = is.finite(lsg_v8_io_value(
      summary, "lambda_reference_error")) &&
      lsg_v8_io_value(summary, "lambda_reference_error") <=
        tolerances$lambda_reference,
    archived_lambda_replay = is.finite(lsg_v8_io_value(
      summary, "baseline_replay_max_lambda_error")) &&
      lsg_v8_io_value(summary, "baseline_replay_max_lambda_error") <=
        tolerances$replay_lambda,
    archived_validation_loss_replay = is.finite(lsg_v8_io_value(
      summary, "baseline_replay_max_validation_loss_error")) &&
      lsg_v8_io_value(summary, "baseline_replay_max_validation_loss_error") <=
        tolerances$replay_loss,
    archived_kkt_replay = is.finite(lsg_v8_io_value(
      summary, "baseline_replay_max_kkt_error")) &&
      lsg_v8_io_value(summary, "baseline_replay_max_kkt_error") <=
        tolerances$replay_kkt,
    archived_intercept_replay = is.finite(lsg_v8_io_value(
      summary, "baseline_replay_max_intercept_kkt_error")) &&
      lsg_v8_io_value(summary, "baseline_replay_max_intercept_kkt_error") <=
        tolerances$replay_intercept,
    archived_passes_replay = isTRUE(lsg_v8_io_value(
      summary, "baseline_replay_passes_exact", FALSE)),
    archived_convergence_replay = isTRUE(lsg_v8_io_value(
      summary, "baseline_replay_convergence_exact", FALSE)),
    baseline_reconstruction = is.data.frame(comparison) &&
      all(is.finite(comparison$baseline_reconstruction_error)) &&
      max(comparison$baseline_reconstruction_error) <= tolerances$reconstruction,
    profiled_reconstruction = is.data.frame(comparison) &&
      all(is.finite(comparison$profiled_reconstruction_error)) &&
      max(comparison$profiled_reconstruction_error) <= tolerances$reconstruction,
    kkt_decomposition = is.data.frame(comparison) && all(is.finite(c(
      comparison$baseline_kkt_decomposition_error,
      comparison$profiled_kkt_decomposition_error
    ))) && max(c(
      comparison$baseline_kkt_decomposition_error,
      comparison$profiled_kkt_decomposition_error
    )) <= tolerances$kkt_decomposition,
    profile_roots_converged = isTRUE(lsg_v8_io_value(
      summary, "profile_root_all_converged", FALSE)),
    profiled_intercept_stationary = is.finite(lsg_v8_io_value(
      summary, "profile_max_intercept_kkt")) &&
      lsg_v8_io_value(summary, "profile_max_intercept_kkt") <=
        tolerances$profile_intercept,
    profiled_objective_nonincreasing = is.finite(lsg_v8_io_value(
      summary, "profile_max_objective_increase")) &&
      lsg_v8_io_value(summary, "profile_max_objective_increase") <=
        tolerances$profile_objective_increase
  )
  lsg_v8_io_check_rows(case, checks)
}

lsg_v8_io_hard_checks_tail <- function(payload, case, configuration,
                                        tolerances) {
  curve <- payload$curve
  replay <- payload$v7_replay
  selection <- payload$selection
  finite <- if (is.data.frame(curve) && "point_type" %in% names(curve)) {
    curve[curve$point_type == "finite", , drop = FALSE]
  } else data.frame()
  endpoint <- if (is.data.frame(curve) && "point_type" %in% names(curve)) {
    curve[curve$point_type == "penalty_limit", , drop = FALSE]
  } else data.frame()
  expected_ratios <- c(4096, 2048, 1024, configuration$lambda_relative_grid)
  checks <- list(
    payload_schema = identical(payload$schema_version,
                               "sglasso_lambda_tail_failure_diagnostic_v8"),
    case_identity = is.data.frame(payload$case) && nrow(payload$case) == 1L &&
      identical(as.character(payload$case$case_id), as.character(case$case_id)),
    blinded_no_test = identical(payload$test_used, FALSE) &&
      identical(names(payload$data_sha256), c("training", "validation")),
    unique_curve_columns = is.data.frame(curve) && !anyDuplicated(names(curve)),
    exact_upward_tail = nrow(finite) == 39L &&
      lsg_v8_io_equal(finite$lambda_ratio, expected_ratios, 1e-12),
    exact_direct_v7_replay = is.data.frame(replay) &&
      nrow(replay) == 36L && !anyDuplicated(names(replay)) &&
      lsg_v8_io_equal(replay$lambda_ratio,
                      configuration$lambda_relative_grid, 1e-12),
    one_analytic_endpoint = nrow(endpoint) == 1L &&
      is.infinite(endpoint$lambda_ratio),
    analytic_endpoint_outputs_finite = nrow(endpoint) == 1L &&
      all(is.finite(c(
        endpoint$baseline_intercept,
        endpoint$baseline_training_log_loss,
        endpoint$baseline_validation_log_loss,
        endpoint$baseline_intercept_kkt,
        endpoint$baseline_group_kkt,
        endpoint$profiled_intercept,
        endpoint$profiled_training_log_loss,
        endpoint$profiled_validation_log_loss,
        endpoint$profiled_kkt,
        endpoint$profiled_intercept_kkt,
        endpoint$profiled_group_kkt,
        endpoint$baseline_reconstruction_error,
        endpoint$profiled_reconstruction_error
      ))),
    finite_tail_outputs = nrow(finite) == 39L && all(is.finite(c(
      finite$lambda, finite$baseline_intercept,
      finite$baseline_training_log_loss, finite$baseline_validation_log_loss,
      finite$baseline_objective, finite$baseline_kkt,
      finite$baseline_intercept_kkt, finite$baseline_group_kkt,
      finite$profiled_intercept, finite$profiled_training_log_loss,
      finite$profiled_validation_log_loss, finite$profiled_objective,
      finite$profiled_kkt, finite$profiled_intercept_kkt,
      finite$profiled_group_kkt
    ))),
    finite_coefficients_and_probabilities = is.data.frame(curve) &&
      all(curve$baseline_beta_finite) &&
      all(curve$baseline_probability_finite) &&
      all(curve$profiled_probability_finite),
    profiled_beta_unchanged = is.data.frame(curve) &&
      all(curve$profile_beta_unchanged),
    profiled_tail_roots_converged = nrow(finite) == 39L &&
      all(finite$profile_root_converged),
    profiled_tail_intercepts_stationary = nrow(finite) == 39L &&
      max(finite$profiled_intercept_kkt) <= tolerances$profile_intercept,
    profiled_tail_objective_nonincreasing = nrow(finite) == 39L &&
      max(finite$profiled_objective_minus_baseline) <=
        tolerances$profile_objective_increase,
    old_v7_grid_embedded = isTRUE(lsg_v8_io_value(
      selection, "old_v7_grid_embedded_exactly", FALSE)),
    lambda_reference_replay = is.finite(lsg_v8_io_value(
      selection, "lambda_reference_error")) &&
      lsg_v8_io_value(selection, "lambda_reference_error") <=
        tolerances$lambda_reference,
    old_v7_validation_loss_replay = is.finite(lsg_v8_io_value(
      selection, "old_v7_grid_max_validation_loss_error")) &&
      lsg_v8_io_value(selection, "old_v7_grid_max_validation_loss_error") <=
        tolerances$replay_loss,
    old_v7_kkt_replay = is.finite(lsg_v8_io_value(
      selection, "old_v7_grid_max_kkt_error")) &&
      lsg_v8_io_value(selection, "old_v7_grid_max_kkt_error") <=
        tolerances$replay_kkt,
    old_v7_intercept_replay = is.finite(lsg_v8_io_value(
      selection, "old_v7_grid_max_intercept_kkt_error")) &&
      lsg_v8_io_value(selection, "old_v7_grid_max_intercept_kkt_error") <=
        tolerances$replay_intercept,
    old_v7_passes_replay = isTRUE(lsg_v8_io_value(
      selection, "old_v7_grid_passes_exact", FALSE)),
    old_v7_convergence_replay = isTRUE(lsg_v8_io_value(
      selection, "old_v7_grid_convergence_exact", FALSE)),
    direct_replay_reconstruction = is.data.frame(replay) && all(is.finite(c(
      replay$baseline_reconstruction_error,
      replay$profiled_reconstruction_error
    ))) && max(c(
      replay$baseline_reconstruction_error,
      replay$profiled_reconstruction_error
    )) <= tolerances$reconstruction,
    ratio_512_anchor_replay = is.finite(lsg_v8_io_value(
      selection, "ratio_512_anchor_loss_error")) &&
      lsg_v8_io_value(selection, "ratio_512_anchor_loss_error") <=
        tolerances$tail_anchor_loss,
    analytic_endpoint_replay = is.finite(lsg_v8_io_value(
      selection, "analytic_endpoint_loss_error")) &&
      lsg_v8_io_value(selection, "analytic_endpoint_loss_error") <=
        tolerances$tail_endpoint_loss,
    reconstruction = is.data.frame(curve) && all(is.finite(c(
      curve$baseline_reconstruction_error,
      curve$profiled_reconstruction_error
    ))) && max(c(
      curve$baseline_reconstruction_error,
      curve$profiled_reconstruction_error
    )) <= tolerances$reconstruction,
    finite_kkt_decomposition = nrow(finite) == 39L && all(is.finite(c(
      finite$baseline_kkt_decomposition_error,
      finite$profiled_kkt_decomposition_error
    ))) && max(c(
      finite$baseline_kkt_decomposition_error,
      finite$profiled_kkt_decomposition_error
    )) <= tolerances$kkt_decomposition
  )
  lsg_v8_io_check_rows(case, checks)
}

lsg_v8_io_hard_checks_grpreg <- function(payload, case) {
  summary <- payload$summary
  points <- payload$points
  raw <- payload$raw
  expected_error <- as.character(case$expected_v7_message)
  reason <- as.character(lsg_v8_io_value(summary, "primary_reason", ""))
  raw_lambda <- if (is.list(raw)) as.numeric(raw$lambda) else numeric()
  raw_iterations <- if (is.list(raw)) as.numeric(raw$iterations) else numeric()
  raw_beta <- tryCatch(as.matrix(raw$fit$beta), error = function(condition) NULL)
  raw_probability <- tryCatch(
    as.matrix(raw$probability), error = function(condition) NULL
  )
  classification <- if (is.list(raw)) raw$path_classification else NULL
  path_length <- if (is.data.frame(points)) nrow(points) else 0L
  raw_coefficient_finite <- if (!is.null(raw_beta) &&
      ncol(raw_beta) == path_length) {
    apply(raw_beta, 2L, function(value) all(is.finite(value)))
  } else logical()
  raw_probability_finite <- if (!is.null(raw_probability) &&
      ncol(raw_probability) == path_length) {
    apply(raw_probability, 2L, function(value) {
      all(is.finite(value)) && all(value >= 0 & value <= 1)
    })
  } else logical()
  recomputed_flags <- tryCatch(
    lsg_grpreg_reason_flags_v8(payload), error = function(condition) NULL
  )
  recomputed_reason <- tryCatch(
    lsg_grpreg_primary_reason_v8(payload), error = function(condition) ""
  )
  controls <- payload$controls
  data_summary <- payload$data_summary
  prefix_replay <- tryCatch(
    lsg_grpreg_usable_prefix_v8(
      iterations = raw_iterations,
      max_iterations = controls$max_iterations,
      finite_coefficient = raw_coefficient_finite,
      finite_probability = raw_probability_finite,
      finite_validation_loss = if (is.data.frame(points) &&
          "validation_log_loss" %in% names(points)) {
        is.finite(points$validation_log_loss)
      } else logical(),
      warnings = unique(c(raw$fit_warnings, raw$prediction_warnings))
    ),
    error = function(condition) NULL
  )
  replay_prefix_rows <- if (is.list(prefix_replay)) {
    which(prefix_replay$member)
  } else integer()
  replay_best <- if (length(replay_prefix_rows) && is.data.frame(points) &&
      "validation_log_loss" %in% names(points)) {
    replay_prefix_rows[which.min(points$validation_log_loss[replay_prefix_rows])]
  } else NA_integer_
  raw_finite_count <- if (is.data.frame(points) && all(c(
      "finite_coefficient", "finite_probability", "finite_validation_loss",
      "legacy_numerically_eligible"
    ) %in% names(points))) {
    sum(points$finite_coefficient %in% TRUE &
          points$finite_probability %in% TRUE &
          points$finite_validation_loss %in% TRUE)
  } else NA_integer_
  checks <- list(
    payload_schema = identical(payload$schema_version,
                               "raw_grpreg_group_lasso_diagnostic_v8"),
    case_identity = is.data.frame(summary) && nrow(summary) == 1L &&
      identical(as.integer(summary$task_id), as.integer(case$task_id)) &&
      identical(as.character(summary$case_id), as.character(case$case_id)),
    blinded_no_test = identical(payload$test_used, FALSE) &&
      identical(names(payload$data_sha256), c("training", "validation")),
    raw_fit_evidence_saved_before_selection =
      isTRUE(payload$raw_evidence_constructed_before_selection) &&
      is.list(payload$raw) &&
      all(c("fit", "lambda", "iterations", "fit_warnings",
            "prediction_warnings", "path_classification") %in%
          names(payload$raw)),
    raw_schema_consistent = isTRUE(lsg_v8_io_value(
      summary, "raw_schema_compatible", FALSE)),
    point_count_matches_raw_path = is.data.frame(points) &&
      nrow(points) == as.integer(lsg_v8_io_value(
        summary, "returned_path_length", -1L)),
    unique_ordered_path_keys = is.data.frame(points) && path_length > 0L &&
      all(c("case_id", "lambda_index") %in% names(points)) &&
      !anyNA(points[c("case_id", "lambda_index")]) &&
      !anyDuplicated(points[c("case_id", "lambda_index")]) &&
      identical(as.integer(points$lambda_index), seq_len(path_length)),
    raw_lambda_and_iterations_match_points = is.data.frame(points) &&
      all(c("lambda", "iterations") %in% names(points)) &&
      lsg_v8_io_equal(as.numeric(points$lambda), raw_lambda, 0) &&
      lsg_v8_io_equal(as.numeric(points$iterations), raw_iterations, 0),
    raw_matrix_shapes_match_path = !is.null(raw_beta) &&
      !is.null(raw_probability) && ncol(raw_beta) == path_length &&
      ncol(raw_probability) == path_length,
    frozen_grpreg_controls_replayed = is.list(controls) &&
      identical(controls$penalty, "grLasso") &&
      identical(controls$family, "binomial") &&
      lsg_v8_io_equal(controls$alpha, 1, 0) &&
      lsg_v8_io_equal(controls$gamma_argument, 3, 0) &&
      identical(as.integer(controls$nlambda), 30L) &&
      lsg_v8_io_equal(controls$lambda_min_ratio, 0.05, 1e-15) &&
      identical(controls$log_lambda, TRUE) &&
      lsg_v8_io_equal(controls$tolerance, 1e-7, 1e-15) &&
      identical(as.integer(controls$max_iterations), 1000000L) &&
      identical(as.integer(controls$dfmax), 600L) &&
      identical(as.integer(controls$gmax), 200L) &&
      identical(controls$returnX, FALSE),
    raw_rows_match_blinded_data_summary = is.data.frame(data_summary) &&
      nrow(data_summary) == 1L && !is.null(raw_beta) &&
      !is.null(raw_probability) &&
      nrow(raw_beta) == as.integer(data_summary$predictor_count) + 1L &&
      nrow(raw_probability) == as.integer(data_summary$validation_sample_size) &&
      identical(as.integer(controls$dfmax),
                as.integer(data_summary$predictor_count)) &&
      identical(as.integer(controls$gmax),
                as.integer(data_summary$group_count)) &&
      length(raw$group_multiplier) == as.integer(data_summary$group_count) &&
      all(is.finite(raw$group_multiplier)) &&
      max(abs(raw$group_multiplier - sqrt(3))) <= 1e-14,
    raw_shapes_match_summary = !is.null(raw_beta) &&
      !is.null(raw_probability) &&
      identical(as.integer(lsg_v8_io_value(
        summary, "coefficient_path_columns", -1L)), ncol(raw_beta)) &&
      identical(as.integer(lsg_v8_io_value(
        summary, "prediction_path_columns", -1L)), ncol(raw_probability)) &&
      identical(as.integer(lsg_v8_io_value(
        summary, "iteration_vector_length", -1L)), length(raw_iterations)),
    finite_flags_recomputed_from_raw = is.data.frame(points) &&
      all(c("finite_coefficient", "finite_probability") %in% names(points)) &&
      identical(as.logical(points$finite_coefficient),
                as.logical(raw_coefficient_finite)) &&
      identical(as.logical(points$finite_probability),
                as.logical(raw_probability_finite)),
    validation_finiteness_self_consistent = is.data.frame(points) &&
      all(c("validation_log_loss", "finite_validation_loss") %in%
          names(points)) &&
      identical(as.logical(points$finite_validation_loss),
                is.finite(points$validation_log_loss)),
    raw_classification_matches_points = is.list(classification) &&
      all(c("legacy_point_converged",
            "legacy_returned_lower_boundary_excluded") %in% names(points)) &&
      identical(as.logical(points$legacy_point_converged),
        as.logical(lsg_align_logical_v8(
          classification$point_converged, path_length))) &&
      identical(as.logical(points$legacy_returned_lower_boundary_excluded),
        as.logical(lsg_align_logical_v8(
          classification$returned_lower_boundary_excluded, path_length))),
    raw_classification_matches_summary = is.list(classification) &&
      identical(isTRUE(lsg_v8_io_value(
        summary, "path_complete", FALSE)),
        path_length == as.integer(lsg_v8_io_value(
          summary, "requested_path_length", -1L))) &&
      identical(isTRUE(lsg_v8_io_value(
        summary, "total_iteration_limit_reached", FALSE)),
        isTRUE(classification$total_iteration_limit_reached)) &&
      identical(isTRUE(lsg_v8_io_value(
        summary, "legacy_path_termination_acceptable", FALSE)),
        isTRUE(classification$path_termination_acceptable)) &&
      identical(as.character(lsg_v8_io_value(
        summary, "legacy_unexpected_warnings", "")),
        paste(classification$unexpected_warnings, collapse = " | ")),
    raw_counts_match_points = is.finite(raw_finite_count) &&
      identical(as.integer(lsg_v8_io_value(
        summary, "raw_finite_point_count", -1L)),
        as.integer(raw_finite_count)) &&
      identical(as.integer(lsg_v8_io_value(
        summary, "legacy_numerically_eligible_count", -1L)),
        sum(points$legacy_numerically_eligible %in% TRUE)),
    raw_iteration_summaries_recomputed = length(raw_iterations) == path_length &&
      all(is.finite(raw_iterations)) &&
      lsg_v8_io_equal(lsg_v8_io_value(summary, "iteration_sum"),
                      sum(raw_iterations), 0) &&
      lsg_v8_io_equal(lsg_v8_io_value(summary, "iteration_max"),
                      max(raw_iterations), 0),
    usable_prefix_recomputed_from_raw = is.list(prefix_replay) &&
      all(c("cumulative_iterations",
            "completed_before_total_iteration_budget",
            "usable_prefix_candidate", "usable_prefix_member") %in%
          names(points)) &&
      lsg_v8_io_equal(points$cumulative_iterations,
                      prefix_replay$cumulative_iterations, 0) &&
      identical(as.logical(points$completed_before_total_iteration_budget),
                as.logical(prefix_replay$completed_before_total_iteration_budget)) &&
      identical(as.logical(points$usable_prefix_candidate),
                as.logical(prefix_replay$candidate)) &&
      identical(as.logical(points$usable_prefix_member),
                as.logical(prefix_replay$member)) &&
      identical(as.integer(lsg_v8_io_value(
        summary, "usable_prefix_length", -1L)),
        as.integer(prefix_replay$length)) &&
      identical(as.integer(lsg_v8_io_value(
        summary, "usable_prefix_best_lambda_index", -1L)),
        as.integer(replay_best)) &&
      (is.na(replay_best) || lsg_v8_io_equal(lsg_v8_io_value(
        summary, "usable_prefix_best_validation_log_loss"),
        points$validation_log_loss[replay_best], 0)) &&
      identical(payload$usable_prefix_rule, prefix_replay$rule),
    original_zero_eligible_reproduced = identical(as.integer(lsg_v8_io_value(
      summary, "legacy_numerically_eligible_count", -1L)), 0L),
    original_selection_error_reproduced = identical(
      as.character(payload$selection_error), expected_error) &&
      identical(as.character(lsg_v8_io_value(summary, "selection_error", "")),
                expected_error),
    primary_reason_classified = length(reason) == 1L && nzchar(reason) &&
      !reason %in% c("unknown", "original_failure_not_reproduced"),
    primary_reason_recomputed_from_raw_summary =
      identical(reason, recomputed_reason),
    nonexclusive_reason_flags_saved = is.data.frame(payload$reason_flags) &&
      nrow(payload$reason_flags) == 1L && ncol(payload$reason_flags) >= 1L &&
      all(vapply(payload$reason_flags, is.logical, logical(1))) &&
      !anyNA(payload$reason_flags) && any(unlist(payload$reason_flags)),
    nonexclusive_reason_flags_recomputed = is.data.frame(recomputed_flags) &&
      identical(payload$reason_flags, recomputed_flags)
  )
  lsg_v8_io_check_rows(case, checks)
}

lsg_v8_io_hard_checks <- function(payload, case, configuration, tolerances) {
  type <- as.character(case$diagnostic_type)
  checks <- switch(
    type,
    sglasso_intercept = lsg_v8_io_hard_checks_intercept(
      payload, case, configuration, tolerances),
    sglasso_lambda_tail = lsg_v8_io_hard_checks_tail(
      payload, case, configuration, tolerances),
    grpreg_group_lasso_raw = lsg_v8_io_hard_checks_grpreg(payload, case),
    stop("Unknown V8 diagnostic type.", call. = FALSE)
  )
  lsg_v8_io_assert(is.data.frame(checks) &&
    identical(names(checks), c("case_id", "task_id", "check", "passed")) &&
    !anyNA(checks$check) && !anyDuplicated(checks$check),
    "Invalid V8 hard-check schema.")
  checks
}

lsg_v8_io_outcome_row <- function(case, outcome, value) {
  type <- if (is.logical(value)) "logical" else if (is.numeric(value)) {
    "numeric"
  } else "character"
  data.frame(
    case_id = as.character(case$case_id), task_id = as.integer(case$task_id),
    outcome = outcome, value_type = type,
    value_logical = if (type == "logical") as.logical(value) else NA,
    value_numeric = if (type == "numeric") as.numeric(value) else NA_real_,
    value_character = if (type == "character") as.character(value) else NA_character_,
    acceptance_gate = FALSE, stringsAsFactors = FALSE
  )
}

lsg_v8_io_scientific_outcomes <- function(payload, case) {
  type <- as.character(case$diagnostic_type)
  values <- switch(type,
    sglasso_intercept = {
      comparison <- payload$comparison
      valid <- is.data.frame(comparison) && nrow(comparison) > 0L && all(c(
        "baseline_kkt", "profiled_kkt", "profiled_objective_minus_baseline",
        "baseline_validation_log_loss", "profiled_validation_log_loss",
        "baseline_group_kkt", "profiled_group_kkt", "profile_beta_unchanged"
      ) %in% names(comparison))
      n <- if (valid) nrow(comparison) else NA_real_
      brought_internal <- if (valid) sum(
        comparison$baseline_kkt > 2e-6 & comparison$profiled_kkt <= 2e-6
      ) else NA_real_
      brought_study <- if (valid) sum(
        comparison$baseline_kkt > 2.05e-6 &
          comparison$profiled_kkt <= 2.05e-6
      ) else NA_real_
      list(
        finite_path_points = n,
        internal_kkt_brought_below_count = brought_internal,
        internal_kkt_brought_below_proportion = brought_internal / n,
        study_kkt_brought_below_count = brought_study,
        study_kkt_brought_below_proportion = brought_study / n,
        profiled_internal_kkt_pass_count = if (valid) {
          sum(comparison$profiled_kkt <= 2e-6)
        } else NA_real_,
        profiled_internal_kkt_pass_proportion = if (valid) {
          mean(comparison$profiled_kkt <= 2e-6)
        } else NA_real_,
        profiled_study_kkt_pass_count = if (valid) {
          sum(comparison$profiled_kkt <= 2.05e-6)
        } else NA_real_,
        profiled_study_kkt_pass_proportion = if (valid) {
          mean(comparison$profiled_kkt <= 2.05e-6)
        } else NA_real_,
        mean_objective_change_after_intercept_profile = if (valid) {
          mean(comparison$profiled_objective_minus_baseline)
        } else NA_real_,
        maximum_objective_change_after_intercept_profile = if (valid) {
          max(comparison$profiled_objective_minus_baseline)
        } else NA_real_,
        mean_validation_log_loss_change_after_intercept_profile = if (valid) {
          mean(comparison$profiled_validation_log_loss -
                 comparison$baseline_validation_log_loss)
        } else NA_real_,
        mean_group_kkt_change_after_intercept_profile = if (valid) {
          mean(comparison$profiled_group_kkt - comparison$baseline_group_kkt)
        } else NA_real_,
        maximum_selected_group_count_change = if (valid &&
            all(comparison$profile_beta_unchanged)) 0 else NA_real_,
        profile_path_has_study_kkt_candidate = lsg_v8_io_value(
          payload$summary, "profile_path_has_study_kkt_candidate", FALSE),
        source_path_min_kkt = lsg_v8_io_value(
          payload$summary, "source_path_min_kkt", NA_real_)
      )
    },
    sglasso_lambda_tail = list(
      diagnostic_range_resolved = lsg_v8_io_value(
        payload$selection, "diagnostic_range_resolved", FALSE),
      diagnostic_selected_is_interior = lsg_v8_io_value(
        payload$selection, "diagnostic_selected_is_interior", FALSE),
      selected_lambda_ratio = lsg_v8_io_value(
        payload$selection, "diagnostic_selected_lambda_ratio", NA_real_),
      selected_validation_log_loss = lsg_v8_io_value(
        payload$selection, "diagnostic_selected_validation_log_loss", NA_real_),
      loss_gap_to_next_candidate = lsg_v8_io_value(
        payload$selection, "diagnostic_loss_gap_to_next_candidate", NA_real_)
    ),
    grpreg_group_lasso_raw = c(list(
        primary_reason = lsg_v8_io_value(
          payload$summary, "primary_reason", "unknown"),
        returned_path_length = lsg_v8_io_value(
          payload$summary, "returned_path_length", NA_real_),
        usable_prefix_length = lsg_v8_io_value(
          payload$summary, "usable_prefix_length", NA_real_),
        iteration_budget_reached = lsg_v8_io_value(
          payload$summary, "total_iteration_limit_reached", FALSE)
      ), if (is.data.frame(payload$reason_flags) &&
          nrow(payload$reason_flags) == 1L) {
        stats::setNames(
          as.list(payload$reason_flags[1L, , drop = TRUE]),
          paste0("reason_", names(payload$reason_flags))
        )
      } else list())
  )
  rows <- Map(function(name, value) lsg_v8_io_outcome_row(
    case, name, value), names(values), values)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

lsg_v8_io_forbidden_payload_paths <- function(value, prefix = "") {
  if (!is.list(value)) return(character())
  output <- character()
  for (name in names(value)) {
    path <- if (nzchar(prefix)) paste(prefix, name, sep = "/") else name
    if (!identical(name, "test_used") && grepl(
      "^(X_test|y_test|n_test|truth|true_|active_group|signal_pattern|test_probability|test_eta|observation_id)",
      name, ignore.case = TRUE
    )) output <- c(output, path)
    output <- c(output, lsg_v8_io_forbidden_payload_paths(value[[name]], path))
  }
  unique(output)
}

lsg_v8_io_shard_receipt <- function(shard, path) {
  defaults <- lsg_v8_io_defaults()
  list(
    schema_version = defaults$receipt_schema,
    case_id = shard$task$case_id,
    scientific_signature = shard$scientific_signature,
    file = basename(path), bytes = as.numeric(file.info(path)$size),
    sha256 = lsg_file_hash_v7(path)
  )
}

lsg_v8_io_publish_shard <- function(output, shard) {
  relative <- file.path("shards", shard$task$shard_file)
  path <- lsg_atomic_v7(shard, output, relative)
  lsg_atomic_v7(lsg_v8_io_shard_receipt(shard, path), output,
                paste0(relative, ".receipt.rds"))
  invisible(path)
}

lsg_v8_io_read_shard <- function(output, spec, task, deep = TRUE) {
  defaults <- lsg_v8_io_defaults()
  relative <- file.path("shards", task$shard_file)
  path <- lsg_path_v7(output, relative)
  receipt_path <- lsg_path_v7(output, paste0(relative, ".receipt.rds"))
  lsg_v8_io_assert(file.exists(path) && file.exists(receipt_path), paste(
    "Incomplete V8 diagnostic shard publication for", task$case_id))
  shard <- readRDS(path)
  receipt <- readRDS(receipt_path)
  lsg_v8_io_assert(
    identical(shard$schema_version, defaults$shard_schema) &&
      identical(shard$scientific_signature, spec$scientific_signature) &&
      identical(shard$task, task) &&
      identical(receipt, lsg_v8_io_shard_receipt(shard, path)),
    paste("V8 diagnostic shard/receipt mismatch for", task$case_id)
  )
  if (isTRUE(deep)) {
    case <- spec$cases[spec$cases$case_id == task$case_id, , drop = FALSE]
    checks <- lsg_v8_io_hard_checks(
      shard$payload, case, spec$configuration, spec$tolerances)
    forbidden <- lsg_v8_io_forbidden_payload_paths(shard$payload)
    checks <- rbind(checks, data.frame(
      case_id = case$case_id, task_id = case$task_id,
      check = "no_test_or_truth_payload_fields",
      passed = length(forbidden) == 0L, stringsAsFactors = FALSE
    ))
    outcomes <- lsg_v8_io_scientific_outcomes(shard$payload, case)
    lsg_v8_io_assert(identical(shard$hard_checks, checks) &&
      identical(shard$scientific_outcomes, outcomes), paste(
        "Stored V8 checks/outcomes disagree with recomputation for", task$case_id))
  }
  shard
}

lsg_v8_io_collect <- function(output, spec, deep = TRUE) {
  tasks <- spec$tasks
  present <- logical(nrow(tasks))
  shards <- vector("list", nrow(tasks))
  for (index in seq_len(nrow(tasks))) {
    relative <- file.path("shards", tasks$shard_file[index])
    path <- lsg_path_v7(output, relative)
    receipt <- lsg_path_v7(output, paste0(relative, ".receipt.rds"))
    if (xor(file.exists(path), file.exists(receipt))) stop(paste(
      "Incomplete V8 atomic shard/receipt pair; preserve and inspect:",
      tasks$case_id[index]), call. = FALSE)
    present[index] <- file.exists(path) && file.exists(receipt)
    if (present[index]) shards[[index]] <- lsg_v8_io_read_shard(
      output, spec, tasks[index, , drop = FALSE], deep = deep)
  }
  names(shards) <- tasks$case_id
  list(present = present, shards = shards)
}

lsg_v8_io_dispatch <- function(root, audit, spec, task) {
  case <- spec$cases[spec$cases$case_id == task$case_id, , drop = FALSE]
  lsg_v8_io_assert(nrow(case) == 1L, "V8 diagnostic task has no frozen case.")
  scenario <- spec$design[spec$design$scenario_index == case$scenario_index,
                          , drop = FALSE]
  lsg_v8_io_assert(nrow(scenario) == 1L,
                   paste("Missing frozen scenario for", case$case_id))
  payload <- switch(
    as.character(case$diagnostic_type),
    sglasso_intercept = {
      source <- audit$blinded_sources[[as.character(case$task_id)]]
      lsg_run_intercept_case_v8(case, scenario, source, spec$configuration)
    },
    sglasso_lambda_tail = {
      source <- audit$blinded_sources[[as.character(case$task_id)]]
      lsg_run_tail_case_v8(case, scenario, source, spec$configuration)
    },
    grpreg_group_lasso_raw = lsg_run_grpreg_case_v8(root, case),
    stop("Unknown V8 diagnostic dispatch type.", call. = FALSE)
  )
  forbidden <- lsg_v8_io_forbidden_payload_paths(payload)
  hard_checks <- lsg_v8_io_hard_checks(
    payload, case, spec$configuration, spec$tolerances)
  hard_checks <- rbind(hard_checks, data.frame(
    case_id = case$case_id, task_id = case$task_id,
    check = "no_test_or_truth_payload_fields",
    passed = length(forbidden) == 0L, stringsAsFactors = FALSE
  ))
  scientific_outcomes <- lsg_v8_io_scientific_outcomes(payload, case)
  list(payload = payload, hard_checks = hard_checks,
       scientific_outcomes = scientific_outcomes)
}

lsg_v8_io_bind <- function(items) {
  items <- Filter(function(item) is.data.frame(item) && nrow(item), items)
  if (!length(items)) return(data.frame())
  columns <- unique(unlist(lapply(items, names), use.names = FALSE))
  items <- lapply(items, function(item) {
    absent <- setdiff(columns, names(item))
    for (name in absent) item[[name]] <- NA
    item[columns]
  })
  output <- do.call(rbind, items)
  rownames(output) <- NULL
  output
}

lsg_v8_io_unique_key <- function(frame, columns) {
  is.data.frame(frame) && nrow(frame) > 0L &&
    all(columns %in% names(frame)) && !anyNA(frame[columns]) &&
    !any(duplicated(frame[columns]))
}

lsg_v8_io_final_integrity_checks <- function(tables) {
  required <- c(
    "intercept_profile_comparison", "intercept_profile_summary",
    "lambda_tail_curve", "lambda_tail_v7_replay",
    "lambda_tail_selection", "grpreg_raw_points", "grpreg_raw_summary",
    "grpreg_reason_flags", "hard_checks", "scientific_outcomes",
    "source_zero_eligible_path_ranking"
  )
  complete_schema <- is.list(tables) && all(required %in% names(tables))
  checks <- c(
    required_final_tables_present = complete_schema,
    intercept_profile_comparison_has_108_rows = complete_schema &&
      nrow(tables$intercept_profile_comparison) == 108L,
    intercept_profile_summary_has_3_rows = complete_schema &&
      nrow(tables$intercept_profile_summary) == 3L,
    lambda_tail_curve_has_40_rows = complete_schema &&
      nrow(tables$lambda_tail_curve) == 40L,
    lambda_tail_v7_replay_has_36_rows = complete_schema &&
      nrow(tables$lambda_tail_v7_replay) == 36L,
    lambda_tail_selection_has_1_row = complete_schema &&
      nrow(tables$lambda_tail_selection) == 1L,
    grpreg_raw_points_nonempty = complete_schema &&
      nrow(tables$grpreg_raw_points) > 0L,
    grpreg_raw_summary_has_3_rows = complete_schema &&
      nrow(tables$grpreg_raw_summary) == 3L,
    grpreg_reason_flags_has_3_rows = complete_schema &&
      nrow(tables$grpreg_reason_flags) == 3L,
    source_zero_path_ranking_has_76_rows = complete_schema &&
      nrow(tables$source_zero_eligible_path_ranking) == 76L,
    no_duplicate_columns = complete_schema && all(vapply(
      tables[required], function(frame) {
        is.data.frame(frame) && !anyDuplicated(names(frame))
      }, logical(1))),
    unique_intercept_profile_comparison_keys = complete_schema &&
      lsg_v8_io_unique_key(
        tables$intercept_profile_comparison, c("case_id", "lambda_index")),
    unique_intercept_profile_summary_keys = complete_schema &&
      lsg_v8_io_unique_key(tables$intercept_profile_summary, "case_id"),
    unique_lambda_tail_curve_keys = complete_schema &&
      lsg_v8_io_unique_key(
        tables$lambda_tail_curve, c("case_id", "point_type", "lambda_ratio")),
    unique_lambda_tail_v7_replay_keys = complete_schema &&
      lsg_v8_io_unique_key(
        tables$lambda_tail_v7_replay, c("case_id", "lambda_index")),
    unique_lambda_tail_selection_keys = complete_schema &&
      lsg_v8_io_unique_key(tables$lambda_tail_selection, "case_id"),
    unique_grpreg_raw_point_keys = complete_schema &&
      lsg_v8_io_unique_key(
        tables$grpreg_raw_points, c("case_id", "lambda_index")),
    unique_grpreg_summary_keys = complete_schema &&
      lsg_v8_io_unique_key(tables$grpreg_raw_summary, "case_id"),
    unique_grpreg_reason_flag_keys = complete_schema &&
      lsg_v8_io_unique_key(tables$grpreg_reason_flags, "case_id"),
    unique_hard_check_keys = complete_schema &&
      lsg_v8_io_unique_key(tables$hard_checks, c("case_id", "check")),
    unique_scientific_outcome_keys = complete_schema &&
      lsg_v8_io_unique_key(
        tables$scientific_outcomes, c("case_id", "outcome")),
    unique_source_zero_path_ranking_keys = complete_schema &&
      lsg_v8_io_unique_key(
        tables$source_zero_eligible_path_ranking,
        c("task_id", "alpha_index", "d_index"))
  )
  data.frame(check = names(checks), passed = unname(checks),
             stringsAsFactors = FALSE)
}

lsg_v8_io_final_tables <- function(collected, audit) {
  shards <- collected$shards
  types <- vapply(shards, function(shard) {
    as.character(shard$task$diagnostic_type)
  }, character(1))
  intercept <- shards[types == "sglasso_intercept"]
  tail <- shards[types == "sglasso_lambda_tail"]
  grpreg <- shards[types == "grpreg_group_lasso_raw"]
  tables <- list(
    intercept_profile_comparison = lsg_v8_io_bind(lapply(
      intercept, function(shard) shard$payload$comparison)),
    intercept_profile_summary = lsg_v8_io_bind(lapply(
      intercept, function(shard) shard$payload$summary)),
    lambda_tail_curve = lsg_v8_io_bind(lapply(
      tail, function(shard) shard$payload$curve)),
    lambda_tail_v7_replay = lsg_v8_io_bind(lapply(
      tail, function(shard) shard$payload$v7_replay)),
    lambda_tail_selection = lsg_v8_io_bind(lapply(
      tail, function(shard) shard$payload$selection)),
    grpreg_raw_points = lsg_v8_io_bind(lapply(
      grpreg, function(shard) shard$payload$points)),
    grpreg_raw_summary = lsg_v8_io_bind(lapply(
      grpreg, function(shard) shard$payload$summary)),
    grpreg_reason_flags = lsg_v8_io_bind(lapply(
      grpreg, function(shard) {
        flags <- shard$payload$reason_flags
        flags$case_id <- shard$task$case_id
        flags$task_id <- shard$task$source_task_id
        flags[c("case_id", "task_id", setdiff(names(flags),
          c("case_id", "task_id")))]
      })),
    hard_checks = lsg_v8_io_bind(lapply(
      shards, function(shard) shard$hard_checks)),
    scientific_outcomes = lsg_v8_io_bind(lapply(
      shards, function(shard) shard$scientific_outcomes)),
    source_zero_eligible_path_ranking = audit$zero_path_ranking
  )
  tables$final_integrity_checks <- lsg_v8_io_final_integrity_checks(tables)
  tables
}

lsg_v8_io_shard_manifest <- function(output, spec) {
  rows <- lapply(seq_len(nrow(spec$tasks)), function(index) {
    task <- spec$tasks[index, , drop = FALSE]
    relative <- file.path("shards", task$shard_file)
    path <- lsg_path_v7(output, relative)
    receipt_path <- lsg_path_v7(output, paste0(relative, ".receipt.rds"))
    receipt <- readRDS(receipt_path)
    data.frame(
      diagnostic_task_id = task$diagnostic_task_id,
      case_id = task$case_id, source_task_id = task$source_task_id,
      shard_file = relative, shard_bytes = receipt$bytes,
      shard_sha256 = receipt$sha256,
      receipt_file = paste0(relative, ".receipt.rds"),
      receipt_bytes = as.numeric(file.info(receipt_path)$size),
      receipt_sha256 = lsg_file_hash_v7(receipt_path),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

lsg_v8_io_final_names <- function(version) {
  c(
    intercept_profile_comparison = "final/intercept_profile_comparison.csv",
    intercept_profile_summary = "final/intercept_profile_summary.csv",
    lambda_tail_curve = "final/lambda_tail_curve.csv",
    lambda_tail_v7_replay = "final/lambda_tail_v7_replay.csv",
    lambda_tail_selection = "final/lambda_tail_selection.csv",
    grpreg_raw_points = "final/grpreg_raw_points.csv",
    grpreg_raw_summary = "final/grpreg_raw_summary.csv",
    grpreg_reason_flags = "final/grpreg_reason_flags.csv",
    hard_checks = "final/hard_checks.csv",
    scientific_outcomes = "final/scientific_outcomes.csv",
    final_integrity_checks = "final/final_integrity_checks.csv",
    source_zero_eligible_path_ranking =
      "final/source_zero_eligible_path_ranking.csv",
    diagnostic_record = paste0("final/diagnostic_record_", version, ".rds"),
    shard_manifest = paste0("final/shard_manifest_", version, ".csv"),
    acceptance = paste0("final/ACCEPTED_", version, ".rds")
  )
}

lsg_v8_io_validate_attempt_files <- function(output, spec) {
  attempt_directory <- lsg_path_v7(output, "attempts")
  if (!dir.exists(attempt_directory)) return(character())
  attempts <- list.files(
    attempt_directory, recursive = FALSE, all.files = TRUE,
    full.names = FALSE, include.dirs = FALSE, no.. = TRUE
  )
  if (!length(attempts)) return(character())
  lsg_v8_io_assert(all(grepl(
    "^attempt_[0-9]{8}T[0-9]{6}_[0-9]+\\.rds$", attempts
  )), "An unexpected file exists in the V8 attempts directory.")
  relative <- file.path("attempts", sort(attempts, method = "radix"))
  for (file in relative) {
    attempt <- readRDS(lsg_path_v7(output, file))
    lsg_v8_io_assert(is.list(attempt) &&
      identical(attempt$schema_version, lsg_v8_io_defaults()$attempt_schema) &&
      identical(attempt$scientific_signature, spec$scientific_signature) &&
      identical(attempt$runtime, spec$runtime) &&
      is.list(attempt$owner) && is.list(attempt$failures) &&
      is.numeric(attempt$requested_cores) &&
      length(attempt$requested_cores) == 1L &&
      is.finite(attempt$requested_cores) && attempt$requested_cores >= 1,
      paste("Invalid V8 attempt receipt:", file))
  }
  relative
}

lsg_v8_io_assert_exact_files <- function(output, expected, ignored = character()) {
  actual <- list.files(output, recursive = TRUE, all.files = TRUE,
                       include.dirs = FALSE, no.. = TRUE)
  actual <- sort(setdiff(actual, ignored), method = "radix")
  expected <- sort(unique(expected), method = "radix")
  lsg_v8_io_assert(identical(actual, expected), paste(
    "Unexpected, missing, or foreign file in V8 output.",
    "Found:", paste(setdiff(actual, expected), collapse = ","),
    "Missing:", paste(setdiff(expected, actual), collapse = ",")
  ))
  invisible(TRUE)
}

lsg_v8_io_manifest_candidates <- function(output, spec) {
  version <- spec$version
  manifest_name <- paste0("output_manifest_", version, ".csv")
  marker_name <- paste0("COMPLETED_", version, ".txt")
  base <- c(
    "study_specification.rds", "diagnostic_case_inventory.csv",
    "diagnostic_task_grid.csv", "source_input_manifest.csv",
    "code_manifest.csv", "runtime.rds"
  )
  shards <- as.vector(rbind(
    file.path("shards", spec$tasks$shard_file),
    paste0(file.path("shards", spec$tasks$shard_file), ".receipt.rds")
  ))
  attempts <- lsg_v8_io_validate_attempt_files(output, spec)
  final <- unname(lsg_v8_io_final_names(version))
  expected <- sort(c(base, shards, attempts, final), method = "radix")
  lock_files <- list.files(
    lsg_path_v7(output, ".run_lock"), recursive = TRUE, all.files = TRUE,
    include.dirs = FALSE, no.. = TRUE
  )
  lsg_v8_io_assert(identical(lock_files, "owner.rds"),
                   "The active V8 run lock has unexpected contents.")
  lsg_v8_io_assert_exact_files(
    output, expected,
    ignored = c(".run_lock/owner.rds", manifest_name, marker_name)
  )
  expected
}

lsg_v8_io_completion_text <- function(spec, output_manifest_path) {
  paste(
    paste("version", spec$version),
    paste("schema", spec$schema_version),
    paste("scientific_signature", spec$scientific_signature),
    "completed_tasks 7",
    "all_hard_gates_passed TRUE",
    paste("output_manifest_sha256", lsg_file_hash_v7(output_manifest_path)),
    sep = "\n"
  )
}

lsg_v8_io_finalize <- function(output, spec, collected, audit) {
  lsg_v8_io_assert(all(collected$present),
                   "All seven V8 diagnostic shards are required.")
  tables <- lsg_v8_io_final_tables(collected, audit)
  lsg_v8_io_assert(nrow(tables$hard_checks) > 0L &&
    all(tables$hard_checks$passed),
    "A V8 hard gate failed; finalization is blocked and shards are preserved.")
  lsg_v8_io_assert(nrow(tables$scientific_outcomes) > 0L &&
    all(!tables$scientific_outcomes$acceptance_gate),
    "Scientific outcomes must never become acceptance gates.")
  lsg_v8_io_assert(nrow(tables$final_integrity_checks) > 0L &&
    all(tables$final_integrity_checks$passed),
    "A V8 final-table dimension or unique-key integrity gate failed.")
  names <- lsg_v8_io_final_names(spec$version)
  table_names <- setdiff(names(tables), character())
  for (name in table_names) {
    lsg_atomic_v7(tables[[name]], output, names[[name]], "csv",
                  identical_ok = TRUE)
  }
  record <- list(
    schema_version = spec$schema_version,
    version = spec$version,
    scientific_signature = spec$scientific_signature,
    source_scientific_signature = spec$source_scientific_signature,
    tables = tables
  )
  lsg_atomic_v7(record, output, names[["diagnostic_record"]], "rds",
                identical_ok = TRUE)
  shard_manifest <- lsg_v8_io_shard_manifest(output, spec)
  lsg_atomic_v7(shard_manifest, output, names[["shard_manifest"]], "csv",
                identical_ok = TRUE)
  acceptance <- list(
    schema_version = lsg_v8_io_defaults()$acceptance_schema,
    version = spec$version,
    scientific_signature = spec$scientific_signature,
    source_scientific_signature = spec$source_scientific_signature,
    source_runtime_exact_at_fit = identical(spec$runtime, spec$source_runtime),
    completed_tasks = 7L,
    hard_checks = nrow(tables$hard_checks),
    final_integrity_checks = nrow(tables$final_integrity_checks),
    all_hard_gates_passed = all(tables$hard_checks$passed) &&
      all(tables$final_integrity_checks$passed),
    scientific_outcomes_are_not_gates =
      all(!tables$scientific_outcomes$acceptance_gate)
  )
  lsg_v8_io_assert(isTRUE(acceptance$source_runtime_exact_at_fit),
                   "V8 fitting did not use the exact archived V7 runtime.")
  lsg_atomic_v7(acceptance, output, names[["acceptance"]], "rds",
                identical_ok = TRUE)

  manifest_relative <- paste0("output_manifest_", spec$version, ".csv")
  marker_relative <- paste0("COMPLETED_", spec$version, ".txt")
  manifest <- lsg_v8_io_inventory_absolute(
    output, lsg_v8_io_manifest_candidates(output, spec))
  lsg_atomic_v7(manifest, output, manifest_relative, "csv", identical_ok = TRUE)
  manifest_path <- lsg_path_v7(output, manifest_relative)
  lsg_atomic_v7(lsg_v8_io_completion_text(spec, manifest_path), output,
                marker_relative, "text", identical_ok = TRUE)
  invisible(list(tables = tables, acceptance = acceptance,
                 manifest = manifest))
}

lsg_v8_io_initialize <- function(root, source_root, version) {
  audit <- lsg_v8_io_audit_source(root, source_root, require_runtime = TRUE)
  output <- lsg_v8_io_output(root, version)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  spec_path <- lsg_path_v7(output, "study_specification.rds")
  spec <- if (file.exists(spec_path)) readRDS(spec_path) else {
    value <- lsg_v8_io_make_spec(root, audit, version)
    lsg_atomic_v7(value, output, "study_specification.rds")
    value
  }
  lsg_v8_io_validate_spec(spec, root, audit, require_runtime = TRUE)
  lsg_atomic_v7(spec$cases, output, "diagnostic_case_inventory.csv", "csv", TRUE)
  lsg_atomic_v7(spec$tasks, output, "diagnostic_task_grid.csv", "csv", TRUE)
  lsg_atomic_v7(spec$source_input_manifest, output,
                "source_input_manifest.csv", "csv", TRUE)
  lsg_atomic_v7(spec$code_manifest, output, "code_manifest.csv", "csv", TRUE)
  lsg_atomic_v7(spec$runtime, output, "runtime.rds", "rds", TRUE)
  list(output = output, audit = audit, specification = spec)
}

lsg_run_failure_diagnostic_v8 <- function(
    root,
    source_root = NULL,
    version = lsg_v8_io_defaults()$default_version,
    cores = 1L,
    max_seconds = Inf,
    max_tasks = Inf
) {
  start <- proc.time()[["elapsed"]]
  root <- normalizePath(root, mustWork = TRUE)
  cores <- as.integer(cores)
  max_seconds <- as.numeric(max_seconds)
  max_tasks <- as.numeric(max_tasks)
  lsg_v8_io_assert(length(cores) == 1L && !is.na(cores) && cores >= 1L,
                   "cores must be one positive integer.")
  lsg_v8_io_assert(length(max_seconds) == 1L && !is.na(max_seconds) &&
    max_seconds > 0, "max_seconds must be positive.")
  lsg_v8_io_assert(length(max_tasks) == 1L && !is.na(max_tasks) &&
    max_tasks > 0, "max_tasks must be positive.")
  if (nzchar(Sys.getenv("SLURM_CPUS_PER_TASK"))) {
    lsg_v8_io_assert(cores <= as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")),
                     "V8 workers exceed allocated SLURM CPUs.")
  }
  lsg_v8_io_assert(.Platform$OS.type == "unix" || cores == 1L,
                   "Parallel V8 diagnostics require Unix fork workers.")
  lsg_v8_io_require_functions(c(
    "lsg_run_intercept_case_v8", "lsg_run_tail_case_v8",
    "lsg_run_grpreg_case_v8", "lsg_compile_intercept_profile_v8",
    "compile_lsg_core"
  ))

  output <- lsg_v8_io_output(root, version)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  complete <- lsg_path_v7(output, paste0("COMPLETED_", version, ".txt"))
  if (file.exists(complete)) {
    verified <- lsg_verify_failure_diagnostic_v8(
      root, source_root = source_root, version = version, quiet = TRUE)
    lsg_v8_io_assert(identical(
      verified$specification$runtime, lsg_runtime_v7()),
      "A completed V8 run has a different runtime; use verify only.")
    cat("Valid completed V8 diagnostic reused; no fits or files changed.\n")
    return(invisible(verified))
  }

  lock <- lsg_path_v7(output, ".run_lock")
  lsg_v8_io_assert(dir.create(lock, showWarnings = FALSE), paste(
    "V8 run lock exists. Inspect the active/stale owner before removal:", lock))
  owner <- list(
    pid = Sys.getpid(), slurm_job = Sys.getenv("SLURM_JOB_ID"),
    host = Sys.info()[["nodename"]],
    started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  on.exit(lsg_release_lock_v7(lock, owner), add = TRUE)
  lsg_atomic_v7(owner, output, ".run_lock/owner.rds")

  initialized <- lsg_v8_io_initialize(root, source_root, version)
  audit <- initialized$audit
  spec <- initialized$specification
  collected <- lsg_v8_io_collect(output, spec, deep = TRUE)
  pending <- which(!collected$present)
  if (length(pending)) {
    compile_lsg_core(quiet = TRUE)
    lsg_compile_intercept_profile_v8(root, quiet = TRUE)
  }
  limit <- min(length(pending), as.integer(min(max_tasks, .Machine$integer.max)))
  queue <- if (limit) head(pending, limit) else integer()
  attempt_id <- if (length(pending)) paste0("attempt_", format(
    Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"), "_", Sys.getpid()) else NA_character_
  completed_now <- integer()
  failures <- list()
  unstarted <- setdiff(pending, queue)
  while (length(queue) && proc.time()[["elapsed"]] - start < max_seconds) {
    wave <- head(queue, min(cores, 7L))
    queue <- queue[-seq_along(wave)]
    worker <- function(index) {
      task <- spec$tasks[index, , drop = FALSE]
      task_start <- proc.time()[["elapsed"]]
      tryCatch({
        cat("Running V8 diagnostic", task$case_id, "(source task",
            task$source_task_id, ")\n")
        result <- lsg_v8_io_dispatch(root, audit, spec, task)
        shard <- list(
          schema_version = lsg_v8_io_defaults()$shard_schema,
          scientific_signature = spec$scientific_signature,
          task = task, payload = result$payload,
          hard_checks = result$hard_checks,
          scientific_outcomes = result$scientific_outcomes,
          attempt = attempt_id,
          runtime_seconds = proc.time()[["elapsed"]] - task_start
        )
        lsg_v8_io_publish_shard(output, shard)
        list(
          diagnostic_task_id = task$diagnostic_task_id,
          case_id = task$case_id, source_task_id = task$source_task_id,
          status = if (all(result$hard_checks$passed)) "passed" else "gate_failed",
          message = paste(result$hard_checks$check[!result$hard_checks$passed],
                          collapse = ";")
        )
      }, error = function(condition) list(
        diagnostic_task_id = task$diagnostic_task_id,
        case_id = task$case_id, source_task_id = task$source_task_id,
        status = "error", message = conditionMessage(condition)
      ))
    }
    outcome <- if (length(wave) > 1L && cores > 1L) {
      parallel::mclapply(
        wave, worker, mc.cores = min(cores, length(wave)),
        mc.preschedule = FALSE, mc.set.seed = FALSE)
    } else lapply(wave, worker)
    for (item in outcome) {
      if (identical(item$status, "passed")) {
        completed_now <- c(completed_now, item$diagnostic_task_id)
      } else failures[[length(failures) + 1L]] <- item
    }
  }
  unstarted <- c(unstarted, queue)
  if (length(pending)) {
    attempt <- list(
      schema_version = lsg_v8_io_defaults()$attempt_schema,
      scientific_signature = spec$scientific_signature,
      owner = owner, requested_cores = cores,
      max_seconds_soft = max_seconds, max_tasks_soft = max_tasks,
      elapsed_seconds = proc.time()[["elapsed"]] - start,
      completed_now = as.integer(completed_now),
      unstarted_diagnostic_tasks = as.integer(unstarted),
      failures = failures, runtime = lsg_runtime_v7(),
      library_paths = .libPaths(),
      session_info = capture.output(utils::sessionInfo()),
      finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
    lsg_atomic_v7(attempt, output,
                  file.path("attempts", paste0(attempt_id, ".rds")))
  }

  collected <- lsg_v8_io_collect(output, spec, deep = TRUE)
  cat("Valid V8 shard/receipt pairs:", sum(collected$present), "/ 7\n")
  if (length(failures)) print(failures)
  if (!all(collected$present)) stop(paste(
    "V8 diagnostic is incomplete. Valid signature-matching shards are resumable",
    "with the identical command and --source-root."), call. = FALSE)
  hard <- lsg_v8_io_bind(lapply(collected$shards, function(shard) {
    shard$hard_checks
  }))
  if (!nrow(hard) || any(!hard$passed)) {
    print(hard[!hard$passed, , drop = FALSE], row.names = FALSE)
    stop(paste(
      "One or more V8 hard gates failed; finalization is blocked.",
      "All diagnostic shards and receipts are preserved."), call. = FALSE)
  }
  lsg_v8_io_finalize(output, spec, collected, audit)
  cat("Output directory:", output, "\nComplete: TRUE\n")
  invisible(list(
    output = output, complete = TRUE, accepted = TRUE,
    specification = spec
  ))
}

lsg_v8_io_read_csv <- function(path, template = NULL) {
  column_classes <- if (is.null(template)) NA else vapply(
    template,
    function(column) {
      if (is.character(column)) "character" else if (is.logical(column)) {
        "logical"
      } else if (is.integer(column)) "integer" else if (is.numeric(column)) {
        "numeric"
      } else stop("Unsupported V8 exported column type.", call. = FALSE)
    },
    character(1)
  )
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                  na.strings = "NA", colClasses = column_classes)
}

lsg_verify_failure_diagnostic_v8 <- function(
    root,
    source_root = NULL,
    version = lsg_v8_io_defaults()$default_version,
    quiet = FALSE
) {
  root <- normalizePath(root, mustWork = TRUE)
  output <- lsg_v8_io_output(root, version)
  lsg_v8_io_assert(dir.exists(output), "Missing V8 diagnostic output directory.")
  audit <- lsg_v8_io_audit_source(root, source_root, require_runtime = FALSE)
  spec <- readRDS(lsg_path_v7(output, "study_specification.rds"))
  lsg_v8_io_validate_spec(spec, root, audit, require_runtime = FALSE)
  lsg_v8_io_assert(identical(spec$runtime, spec$source_runtime),
                   "The accepted V8 fit runtime was not the exact V7 runtime.")

  manifest_relative <- paste0("output_manifest_", version, ".csv")
  marker_relative <- paste0("COMPLETED_", version, ".txt")
  manifest_path <- lsg_path_v7(output, manifest_relative)
  marker_path <- lsg_path_v7(output, marker_relative)
  lsg_v8_io_assert(file.exists(manifest_path) && file.exists(marker_path),
                   "V8 output manifest or completion marker is absent.")
  manifest <- lsg_v8_io_read_csv(manifest_path)
  lsg_v8_io_verify_inventory(output, manifest)
  actual_files <- sort(list.files(
    output, recursive = TRUE, all.files = TRUE,
    include.dirs = FALSE, no.. = TRUE), method = "radix")
  expected_files <- sort(c(manifest$file, manifest_relative, marker_relative),
                         method = "radix")
  lsg_v8_io_assert(identical(actual_files, expected_files),
                   "Unexpected or missing files exist in completed V8 output.")
  expected_marker <- lsg_v8_io_completion_text(spec, manifest_path)
  lsg_v8_io_assert(identical(readLines(marker_path, warn = FALSE),
                             strsplit(expected_marker, "\n", fixed = TRUE)[[1L]]),
                   "The V8 completion marker is not bound to the manifest.")

  collected <- lsg_v8_io_collect(output, spec, deep = TRUE)
  lsg_v8_io_assert(all(collected$present), "Not all seven V8 shards verify.")
  tables <- lsg_v8_io_final_tables(collected, audit)
  lsg_v8_io_assert(all(tables$hard_checks$passed) &&
    all(tables$final_integrity_checks$passed) &&
    all(!tables$scientific_outcomes$acceptance_gate),
    "V8 hard/scientific gate separation failed verification.")
  paths <- lsg_v8_io_final_names(version)
  record <- readRDS(lsg_path_v7(output, paths[["diagnostic_record"]]))
  lsg_v8_io_assert(identical(record$tables, tables) &&
    identical(record$scientific_signature, spec$scientific_signature),
    "V8 final RDS record does not replay exactly from shards.")
  for (name in names(tables)) {
    csv <- lsg_v8_io_read_csv(
      lsg_path_v7(output, paths[[name]]), template = tables[[name]]
    )
    lsg_v8_io_assert(lsg_v8_io_equal(csv, tables[[name]], 1e-10), paste(
      "V8 final CSV does not replay from shards:", name))
  }
  expected_shard_manifest <- lsg_v8_io_shard_manifest(output, spec)
  shard_manifest <- lsg_v8_io_read_csv(
    lsg_path_v7(output, paths[["shard_manifest"]]),
    template = expected_shard_manifest
  )
  lsg_v8_io_assert(lsg_v8_io_equal(
    shard_manifest, expected_shard_manifest, 0),
    "V8 shard manifest disagrees with receipts.")
  acceptance <- readRDS(lsg_path_v7(output, paths[["acceptance"]]))
  lsg_v8_io_assert(
    identical(acceptance$schema_version,
              lsg_v8_io_defaults()$acceptance_schema) &&
      identical(acceptance$scientific_signature, spec$scientific_signature) &&
      identical(acceptance$completed_tasks, 7L) &&
      identical(acceptance$final_integrity_checks,
                nrow(tables$final_integrity_checks)) &&
      isTRUE(acceptance$source_runtime_exact_at_fit) &&
      isTRUE(acceptance$all_hard_gates_passed) &&
      isTRUE(acceptance$scientific_outcomes_are_not_gates),
    "V8 acceptance receipt is invalid."
  )
  if (!quiet) {
    cat("Version:", version,
        "\nScientific signature:", spec$scientific_signature,
        "\nVerified shards: 7 / 7",
        "\nHard checks:", nrow(tables$hard_checks), "passed",
        "\nScientific outcomes:", nrow(tables$scientific_outcomes),
        "(not gates)", "\nOutput directory:", output,
        "\nComplete: TRUE\n")
  }
  invisible(list(
    specification = spec, tables = tables, acceptance = acceptance,
    output_manifest = manifest,
    output = output, complete = TRUE, accepted = TRUE,
    verification_mode = if (identical(spec$runtime, lsg_runtime_v7())) {
      "exact_runtime"
    } else "cross_runtime_artifact_audit"
  ))
}

lsg_v8_io_unit_checks <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  rejects <- function(expression) inherits(try(force(expression), silent = TRUE),
                                           "try-error")
  defaults <- lsg_v8_io_defaults()
  cases <- lsg_read_failure_cases_v8(file.path(
    root, "config", "logistic_sglasso_v7_failure_cases_v8.csv"))
  failure_alignment <- nrow(defaults$source_failures) == 14L &&
    identical(defaults$source_failures$task_id, c(
      125L, 126L, 134L, 141L, 143L, 144L, 146L, 147L,
      150L, 151L, 153L, 155L, 157L, 158L)) &&
    identical(defaults$source_failures$status,
      c("gate_failed", "error", "error", rep("gate_failed", 5L),
        "error", rep("gate_failed", 5L))) &&
    identical(defaults$source_failures$message, c(
      "sglasso_finite_range_resolved",
      rep("No numerically eligible validation candidate was available for Logistic Group Lasso (grpreg).", 2L),
      rep("all_finite_sglasso_paths_represented", 5L),
      "No numerically eligible validation candidate was available for Logistic Group Lasso (grpreg).",
      rep("all_finite_sglasso_paths_represented", 5L))) &&
    sum(defaults$source_failures$status == "gate_failed") == 11L &&
    sum(defaults$source_failures$status == "error") == 3L &&
    identical(defaults$source_failures$message[
      defaults$source_failures$status == "error"], rep(
        "No numerically eligible validation candidate was available for Logistic Group Lasso (grpreg).",
        3L))

  mock_source <- list(payload = list(
    artifacts = list(
      data_sha256 = list(
        training = paste(rep("a", 64L), collapse = ""),
        validation = paste(rep("b", 64L), collapse = ""),
        test = paste(rep("c", 64L), collapse = ""),
        truth = paste(rep("d", 64L), collapse = "")
      ),
      firth_target_original = c(0.1, -0.2),
      y_test = c(0, 1), true_coefficient = c(0, 1)
    ),
    sglasso_tuning = data.frame(
      candidate_id = "mock", alpha_index = 1L, d_index = 1L,
      alpha = 0.1, d = 0.2, lambda_index = 1L, lambda = 1,
      lambda_reference = 1, lambda_reference_type = "d0_null_kkt",
      lambda_relative_to_reference = 1, point_type = "finite",
      validation_log_loss = 0.5, kkt = 1e-6,
      intercept_score = 1e-6, passes = 1L, converged = TRUE,
      numerically_eligible = TRUE, selected_free_d = TRUE,
      n_test = 999L, active_group_count = 3L,
      signal_pattern = "deliberate_truth_leak",
      stringsAsFactors = FALSE
    )
  ))
  blinded <- lsg_v8_io_blind_source(mock_source)
  allowed_tuning <- c(
    "candidate_id", "alpha_index", "d_index", "alpha", "d",
    "lambda_index", "lambda", "lambda_reference", "lambda_reference_type",
    "lambda_relative_to_reference", "point_type", "validation_log_loss",
    "kkt", "intercept_score", "passes", "converged",
    "numerically_eligible", "selected_free_d"
  )
  blinding <- identical(names(blinded), c(
    "data_sha256", "firth_target_original", "sglasso_tuning")) &&
    identical(names(blinded$data_sha256), c("training", "validation")) &&
    identical(names(blinded$sglasso_tuning), allowed_tuning) &&
    length(lsg_v8_io_forbidden_payload_paths(blinded)) == 0L &&
    !any(c("n_test", "active_group_count", "signal_pattern") %in%
      names(blinded$sglasso_tuning))

  group_case <- cases[cases$case_id == "grlasso_126", , drop = FALSE]
  expected_error <- group_case$expected_v7_message
  group_payload <- list(
    schema_version = "raw_grpreg_group_lasso_diagnostic_v8",
    data_sha256 = list(
      training = paste(rep("a", 64L), collapse = ""),
      validation = paste(rep("b", 64L), collapse = "")
    ),
    summary = data.frame(
      task_id = group_case$task_id, case_id = group_case$case_id,
      raw_schema_compatible = TRUE, requested_path_length = 30L,
      returned_path_length = 1L, coefficient_path_columns = 1L,
      iteration_vector_length = 1L, prediction_path_columns = 1L,
      path_complete = FALSE,
      fit_error = "", prediction_error = "",
      legacy_path_termination_acceptable = TRUE,
      legacy_unexpected_warnings = "",
      raw_finite_point_count = 1L,
      legacy_numerically_eligible_count = 0L,
      selection_error = expected_error,
      primary_reason = "incomplete_path",
      usable_prefix_length = 1L,
      usable_prefix_best_lambda_index = 1L,
      usable_prefix_best_validation_log_loss = 0.5,
      total_iteration_limit_reached = FALSE,
      iteration_sum = 1, iteration_max = 1,
      stringsAsFactors = FALSE
    ),
    points = data.frame(
      case_id = group_case$case_id, lambda_index = 1L, lambda = 1,
      iterations = 1, finite_coefficient = TRUE,
      finite_probability = TRUE, validation_log_loss = 0.5,
      finite_validation_loss = TRUE, legacy_point_converged = FALSE,
      legacy_returned_lower_boundary_excluded = FALSE,
      legacy_numerically_eligible = FALSE,
      cumulative_iterations = 1,
      completed_before_total_iteration_budget = TRUE,
      usable_prefix_candidate = TRUE, usable_prefix_member = TRUE,
      stringsAsFactors = FALSE
    ),
    raw = list(fit = list(beta = matrix(0, nrow = 601L, ncol = 1L)),
               probability = matrix(c(0.4, 0.6), ncol = 1L),
               lambda = 1, iterations = 1,
               group_multiplier = rep(sqrt(3), 200L),
               fit_warnings = character(),
               prediction_warnings = character(),
               path_classification = list(
                 point_converged = FALSE,
                 returned_lower_boundary_excluded = FALSE,
                 total_iteration_limit_reached = FALSE,
                 path_termination_acceptable = TRUE,
                 unexpected_warnings = character()
               )),
    controls = list(
      penalty = "grLasso", family = "binomial", alpha = 1,
      gamma_argument = 3, nlambda = 30L, lambda_min_ratio = 0.05,
      log_lambda = TRUE, tolerance = 1e-7,
      max_iterations = 1000000L, dfmax = 600L, gmax = 200L,
      returnX = FALSE
    ),
    data_summary = data.frame(
      predictor_count = 600L, group_count = 200L,
      validation_sample_size = 2L, stringsAsFactors = FALSE
    ),
    usable_prefix_rule = paste(
      "maximal leading sequence with finite coefficients, probabilities and",
      "validation loss, and cumulative grpreg iterations strictly below max.iter;",
      "only the grpreg total-iteration warning is prefix-compatible"
    ),
    raw_evidence_constructed_before_selection = TRUE,
    selection_error = expected_error,
    test_used = FALSE
  )
  group_payload$reason_flags <- lsg_grpreg_reason_flags_v8(group_payload)
  hard <- lsg_v8_io_hard_checks_grpreg(group_payload, group_case)
  science <- lsg_v8_io_scientific_outcomes(group_payload, group_case)
  hard_scientific_separation <- all(hard$passed) &&
    all(!science$acceptance_gate) &&
    !any(science$outcome %in% hard$check)

  temporary <- tempfile("lsg_v8_io_mock_")
  asymmetric <- tempfile("lsg_v8_io_asymmetric_")
  corrupt <- tempfile("lsg_v8_io_corrupt_")
  file_scope <- tempfile("lsg_v8_io_file_scope_")
  dir.create(temporary)
  dir.create(asymmetric)
  dir.create(corrupt)
  dir.create(file_scope)
  on.exit(unlink(c(temporary, asymmetric, corrupt, file_scope), recursive = TRUE,
                 force = TRUE), add = TRUE)
  task <- data.frame(
    diagnostic_task_id = 1L, case_id = group_case$case_id,
    diagnostic_type = group_case$diagnostic_type,
    source_task_id = group_case$task_id,
    source_key = paste(group_case$scenario, group_case$replication, sep = "::"),
    shard_file = "diagnostic_shard_01_mock.rds",
    stringsAsFactors = FALSE
  )
  mock_spec <- list(
    scientific_signature = paste(rep("e", 64L), collapse = ""),
    tasks = task, cases = group_case,
    configuration = list(), tolerances = defaults$tolerances
  )
  shard <- list(
    schema_version = defaults$shard_schema,
    scientific_signature = mock_spec$scientific_signature,
    task = task, payload = group_payload,
    hard_checks = rbind(hard, data.frame(
      case_id = group_case$case_id, task_id = group_case$task_id,
      check = "no_test_or_truth_payload_fields", passed = TRUE,
      stringsAsFactors = FALSE)),
    scientific_outcomes = science,
    attempt = "mock", runtime_seconds = 0
  )
  lsg_v8_io_publish_shard(temporary, shard)
  atomic_pair <- identical(
    lsg_v8_io_read_shard(temporary, mock_spec, task, deep = TRUE), shard)
  changed <- shard
  changed$attempt <- "changed"
  no_clobber <- rejects(lsg_v8_io_publish_shard(temporary, changed))

  relative <- file.path("shards", task$shard_file)
  lsg_atomic_v7(shard, asymmetric, relative)
  asymmetric_rejected <- rejects(lsg_v8_io_collect(
    asymmetric, mock_spec, deep = FALSE))

  lsg_v8_io_publish_shard(corrupt, shard)
  corrupt_receipt_path <- lsg_path_v7(
    corrupt, paste0(relative, ".receipt.rds"))
  corrupt_receipt <- readRDS(corrupt_receipt_path)
  corrupt_receipt$sha256 <- paste(rep("0", 64L), collapse = "")
  saveRDS(corrupt_receipt, corrupt_receipt_path, version = 3)
  corrupt_receipt_rejected <- rejects(lsg_v8_io_read_shard(
    corrupt, mock_spec, task, deep = FALSE))

  lsg_atomic_v7("immutable", temporary, "manifest_fixture.txt", "text")
  inventory <- lsg_v8_io_inventory_absolute(
    temporary, "manifest_fixture.txt")
  manifest_verifies <- isTRUE(lsg_v8_io_verify_inventory(temporary, inventory))
  writeLines("changed", lsg_path_v7(temporary, "manifest_fixture.txt"),
             useBytes = TRUE)
  manifest_corruption_rejected <- rejects(lsg_v8_io_verify_inventory(
    temporary, inventory))

  typed_template <- data.frame(
    integer_column = c(NA_integer_, 2L),
    numeric_column = c(NA_real_, 0.25),
    logical_column = c(NA, TRUE),
    character_column = c(NA_character_, "value"),
    all_na_numeric = c(NA_real_, NA_real_),
    all_na_character = c(NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )
  typed_path <- lsg_path_v7(temporary, "typed_csv_fixture.csv")
  utils::write.csv(typed_template, typed_path, row.names = FALSE,
                   na = "NA", quote = TRUE)
  typed_round_trip <- identical(
    lsg_v8_io_read_csv(typed_path, template = typed_template),
    typed_template
  )

  lsg_atomic_v7("allowed", file_scope, "allowed.txt", "text")
  exact_file_scope_accepts_declared <- !rejects(
    lsg_v8_io_assert_exact_files(file_scope, "allowed.txt")
  )
  lsg_atomic_v7("foreign", file_scope, "foreign.tmp", "text")
  exact_file_scope_rejects_foreign <- rejects(
    lsg_v8_io_assert_exact_files(file_scope, "allowed.txt")
  )
  key_fixture <- data.frame(
    case_id = c("a", "a", "b"), lambda_index = c(1L, 2L, 1L),
    stringsAsFactors = FALSE
  )
  unique_key_accepts_valid <- lsg_v8_io_unique_key(
    key_fixture, c("case_id", "lambda_index")
  )
  duplicate_key_rejected <- !lsg_v8_io_unique_key(
    rbind(key_fixture, key_fixture[1L, , drop = FALSE]),
    c("case_id", "lambda_index")
  )

  checks <- c(
    exact_14_failure_status_message_alignment = failure_alignment,
    exact_seven_case_schema_and_order = identical(cases$case_id,
      defaults$case_ids) && identical(cases$task_id, defaults$case_task_ids),
    blinded_helper_drops_test_and_truth = blinding,
    hard_and_scientific_outcomes_separated = hard_scientific_separation,
    atomic_shard_receipt_round_trip = atomic_pair,
    no_clobber_rejects_changed_shard = no_clobber,
    asymmetric_shard_receipt_rejected = asymmetric_rejected,
    corrupted_receipt_rejected = corrupt_receipt_rejected,
    sha256_inventory_round_trip = manifest_verifies,
    sha256_inventory_corruption_rejected = manifest_corruption_rejected,
    typed_csv_round_trip_preserves_all_na_columns = typed_round_trip,
    exact_output_file_scope_accepts_declared = exact_file_scope_accepts_declared,
    exact_output_file_scope_rejects_foreign = exact_file_scope_rejects_foreign,
    final_table_unique_key_accepts_valid = unique_key_accepts_valid,
    final_table_duplicate_key_rejected = duplicate_key_rejected
  )
  data.frame(check = names(checks), passed = unname(checks),
             stringsAsFactors = FALSE)
}
