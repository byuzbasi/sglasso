# V11.1 executable integration tests. All fitted data here are tiny, synthetic,
# and separately labelled. No frozen V7 task or scientific output is fitted.

lsg_hybrid_io_fixture_v11 <- function(task) {
  frame <- function(n) data.frame(
    case_id = rep(task$case_id, n), task_id = rep(task$source_task_id, n),
    stringsAsFactors = FALSE
  )
  finite <- frame(63L)
  finite$lambda_index <- seq_len(63L)
  candidates <- frame(64L)
  candidates$point_type <- c("penalty_limit", rep("finite", 63L))
  candidates$selected <- c(TRUE, rep(FALSE, 63L))
  checks <- frame(12L)
  checks$check <- lsg_hybrid_expected_checks_v11()
  checks$passed <- TRUE
  fixture <- list(
    schema_version = "logistic_sglasso_hybrid_case_v11",
    case = task[c("case_id", "scenario", "replication", "seed")],
    finite_path = finite, candidates = candidates,
    selected_stability = frame(1L), summary = frame(1L), hard_checks = checks,
    scientific_outcomes = frame(4L), evaluation_data_used = FALSE,
    data_sha256 = list(training = strrep("a", 64L),
                       validation = strrep("b", 64L))
  )
  # Structural I/O double, not numerical evidence; the integration below uses
  # actual fitted payloads and verifies their numerical hard gates separately.
  for (field in names(lsg_hybrid_payload_columns_v11())) {
    for (column in setdiff(lsg_hybrid_payload_columns_v11()[[field]],
                           names(fixture[[field]]))) fixture[[field]][[column]] <- 0
  }
  fixture
}


lsg_hybrid_integration_checks_v11 <- function(root, require_clean = TRUE) {
  rejects <- function(expression) {
    inherits(try(force(expression), silent = TRUE), "try-error")
  }
  runtime_scope <- environment(lsg_run_hybrid_case_v11)
  initial_exports <- vapply(
    c("lsg_kkt_cpp", "lsg_fit_one_hybrid_v11_cpp"), exists, logical(1),
    envir = runtime_scope, inherits = TRUE, mode = "function"
  )
  if (isTRUE(require_clean) && any(initial_exports)) {
    stop("Integration gate must precede numerical unit tests in a fresh Rscript.",
         call. = FALSE)
  }
  scratch <- tempfile("lsg_v11_1_integration_")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)

  # This fixture mimics the *V7 schema*, not a V7 study or source audit. The
  # constructor is tested with synthetic identities; the public run entry point
  # still requires the untouched real V7 audit/signature/runtime.
  scenario <- lsg_design_for_stage_v7(root, "smoke")
  scenario$scenario <- "synthetic_integration_v11_1"
  scenario$n_train <- 60L; scenario$n_validation <- 40L
  scenario$n_test <- 20L; scenario$groups <- 4L; scenario$group_size <- 2L
  scenario$active_group_count <- 2L; scenario$prevalence <- 0.5
  defaults <- lsg_hybrid_diagnostic_defaults_v11()
  cases <- data.frame(
    case_id = defaults$case_ids, task_id = 1:3,
    scenario_index = 1L, scenario = scenario$scenario, replication = 1:3,
    seed = 1110100L + 1:3, alpha = c(0.2, 0.2, 0.1), d = c(0.6, 1, 1),
    stringsAsFactors = FALSE
  )
  sources <- lapply(seq_len(nrow(cases)), function(i) {
    data <- lsg_v8_tuning_data(scenario, cases$seed[i])
    firth <- estimate_groupwise_logistic_target(
      data$X_train, data$y_train, data$group, method = "firth"
    )
    if (!isTRUE(firth$success) || firth$failed_groups != 0L) {
      stop("Synthetic integration Firth target failed; no fallback is allowed.")
    }
    list(data_sha256 = lsg_v8_data_fingerprints(data),
         firth_target_original = firth$target_original,
         sglasso_tuning = data.frame())
  })
  names(sources) <- as.character(cases$task_id)
  audit <- list(
    specification = list(scientific_signature = defaults$source_v7_signature,
                         design = scenario,
                         configuration = lsg_configuration_v7("pilot")),
    source_root = file.path(scratch, "synthetic_source_not_V7"),
    cases = cases, blinded_sources = sources, input_manifest = data.frame()
  )
  spec <- lsg_hybrid_make_spec_v11(root, audit, "synthetic_integration_v11_1")
  resolved <- lsg_hybrid_resolve_tasks_v11(spec)
  missing_design <- audit
  missing_design$specification$design <- NULL
  duplicate_design <- spec
  duplicate_design$scenarios <- rbind(scenario, scenario)
  missing_source <- spec
  missing_source$blinded_sources[[1L]] <- NULL
  missing_row <- spec
  missing_row$scenarios$scenario <- "not_the_requested_scenario"
  wrong_seed <- spec
  wrong_seed$tasks$seed[1L] <- wrong_seed$tasks$seed[1L] + 1L
  forbidden_source <- spec
  forbidden_source$blinded_sources[[1L]]$y_test <- 1
  checks <- c(
    clean_process_before_runtime_initialization = !any(initial_exports),
    V7_schema_design_field_maps_to_three_tasks =
      identical(spec$scenarios, audit$specification$design) &&
      length(resolved) == 3L,
    missing_design_rejected_before_dispatch = rejects(
      lsg_hybrid_make_spec_v11(root, missing_design, "synthetic_bad_design")),
    duplicate_scenario_rejected = rejects(lsg_hybrid_resolve_tasks_v11(duplicate_design)),
    absent_scenario_rejected = rejects(lsg_hybrid_resolve_tasks_v11(missing_row)),
    missing_blinded_source_rejected = rejects(lsg_hybrid_resolve_tasks_v11(missing_source)),
    changed_task_seed_rejected = rejects(lsg_hybrid_resolve_tasks_v11(wrong_seed)),
    forbidden_test_source_rejected = rejects(lsg_hybrid_resolve_tasks_v11(forbidden_source))
  )
  output <- file.path(scratch, "synthetic_run")
  dir.create(file.path(output, "shards"), recursive = TRUE)
  initialized <- list(output = output, specification = spec)
  partial <- lsg_hybrid_execute_v11(root, initialized, 1L, 120, task_limit = 1L)
  first <- lsg_hybrid_shard_paths_v11(output, spec$tasks[1L, , drop = FALSE])$shard
  fingerprint_before <- lsg_file_hash_v7(first)
  mtime_before <- file.info(first)$mtime
  checks <- c(checks,
    actual_dispatch_initializes_base_and_hybrid_cpp = all(vapply(
      c("lsg_kkt_cpp", "lsg_objective_cpp", "lsg_lambda_start_cpp",
        "lsg_fit_one_hybrid_v11_cpp"), exists, logical(1), envir = runtime_scope,
      inherits = TRUE, mode = "function")),
    partial_run_is_not_complete = !partial$complete && partial$completed_tasks == 1L,
    partial_run_has_no_final_acceptance = !dir.exists(file.path(output, "final"))
  )
  # Two tiny workers exercise the same forked dispatcher as TRUBA. The first
  # completed shard must not be fitted or rewritten by the resume.
  resumed <- lsg_hybrid_execute_v11(root, initialized, 2L, 120)
  shards <- lapply(seq_len(3L), function(i) {
    lsg_hybrid_read_shard_v11(output, spec, spec$tasks[i, , drop = FALSE])
  })
  verified <- lsg_hybrid_verify_final_v11(output, spec, shards, quiet = TRUE)
  again <- lsg_hybrid_execute_v11(root, initialized, 1L, 120)
  malformed <- shards[[1L]]
  malformed$payload$finite_path <- NULL
  wrong_task <- shards[[1L]]
  wrong_task$task$seed <- wrong_task$task$seed + 1L
  missing_check <- shards[[1L]]
  missing_check$payload$hard_checks <- head(missing_check$payload$hard_checks, -1L)
  missing_metric <- shards[[1L]]
  missing_metric$payload$finite_path$kkt <- NULL
  wrong_fingerprint <- shards[[1L]]
  wrong_fingerprint$payload$data_sha256$training <- "not-a-hash"
  wrong_signature <- spec
  wrong_signature$scientific_signature <- strrep("f", 64L)
  checks <- c(checks,
    serial_then_parallel_resume_completes_three_tasks =
      resumed$complete && resumed$completed_tasks == 3L,
    saved_shard_is_not_overwritten = identical(fingerprint_before, lsg_file_hash_v7(first)) &&
      identical(mtime_before, file.info(first)$mtime),
    actual_fit_to_final_manifest_verification = isTRUE(verified$accepted),
    completed_resume_performs_no_refit = again$complete &&
      identical(mtime_before, file.info(first)$mtime),
    malformed_payload_rejected = !lsg_hybrid_validate_shard_v11(
      malformed, spec, spec$tasks[1L, , drop = FALSE]),
    missing_hard_check_rejected = !lsg_hybrid_validate_shard_v11(
      missing_check, spec, spec$tasks[1L, , drop = FALSE]),
    missing_numerical_column_rejected = !lsg_hybrid_validate_shard_v11(
      missing_metric, spec, spec$tasks[1L, , drop = FALSE]),
    changed_data_fingerprint_rejected = !lsg_hybrid_validate_shard_v11(
      wrong_fingerprint, spec, spec$tasks[1L, , drop = FALSE]),
    mismatched_shard_task_rejected = !lsg_hybrid_validate_shard_v11(
      wrong_task, spec, spec$tasks[1L, , drop = FALSE]),
    stale_shard_signature_rejected = rejects(lsg_hybrid_read_shard_v11(
      output, wrong_signature, spec$tasks[1L, , drop = FALSE]))
  )
  rejected_output <- file.path(scratch, "synthetic_numerical_rejection")
  dir.create(rejected_output)
  rejected_shards <- shards
  rejected_shards[[1L]]$payload$hard_checks$passed[1L] <- FALSE
  lsg_hybrid_finalize_v11(root, rejected_output, spec, rejected_shards)
  rejection <- lsg_hybrid_verify_final_v11(rejected_output, spec,
                                          rejected_shards, quiet = TRUE)
  checks <- c(checks,
    failed_gate_preserves_results_without_acceptance = !rejection$accepted &&
      file.exists(file.path(rejected_output, "final", "hard_checks.csv")) &&
      !file.exists(file.path(rejected_output, "final", paste0(
        "DIAGNOSTIC_ACCEPTED_", spec$version, ".rds")))
  )
  # Damage only a disposable synthetic output, never a real result.
  writeLines("corrupt synthetic table", file.path(output, "final", "case_summary.csv"))
  checks <- c(checks, corrupted_final_output_rejected = rejects(
    lsg_hybrid_verify_final_v11(output, spec, shards, quiet = TRUE)))
  data.frame(check = paste0("integration__", names(checks)),
             passed = unname(vapply(checks, isTRUE, logical(1))),
             stringsAsFactors = FALSE)
}
