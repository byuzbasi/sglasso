# V8 changes acceptance/reporting only; all V7 fitting functions remain frozen.
allb_configuration_v8 <- function(e, stage) {
  x <- allb_configuration_v7(e, stage)
  x$schema_version <- "all_binary_configuration_v8"
  x$fit_configuration$integration_version <- "all_binary_finite_grid_v8"
  x$fit_configuration$fair_lower_boundary_policy <- "report_not_reject"
  x$fair_tuning$boundary_rule <- "finite_grid_endpoint_is_valid_with_range_warning"
  x$fair_tuning$solver_boundary_rule <- "verified_truncation_is_reported_separately"
  x
}

allb_source_files_v8 <- function() unique(c(allb_source_files_v7(),
  "R/all_binary_workflow_v8.R", "scripts/103_run_all_binary_v8.R",
  "tests/test_all_binary_v8.R", "ALL_BINARY_FINITE_GRID_PROTOCOL_V8.md"))

allb_range_report_v8 <- function(payload, task, configuration) {
  a <- payload$fair_tuning_audit
  data.frame(task_id = task$task_id, replication = task$replication,
    method = a$method,
    requested_common_lower = tail(configuration$fit_configuration$
      fair_common_lambda_relative_grid, 1L),
    selected_lambda_ratio = a$selected_lambda_relative_to_reference,
    selected_common_lower = a$selected_lower_boundary,
    path_truncated_any = a$path_truncated_any,
    selected_returned_endpoint = a$selected_on_solver_lower_boundary,
    selected_returned_endpoint_folds = a$solver_lower_boundary_fold_count,
    warning = ifelse(a$selected_lower_boundary,
      "common_lower_endpoint_selected_beyond_grid_unassessed", ""),
    stringsAsFactors = FALSE)
}

allb_validate_payload_v8 <- function(payload, task, configuration) {
  old <- allb_validate_payload_v4(payload, task, configuration)
  removed <- old$check == "true_common_lambda_lower_boundary_resolved"
  allb_assert(sum(removed) == 1L, "Expected exactly one V7 range gate.")
  a <- payload$fair_tuning_audit
  selected <- payload$inner_candidates
  selected <- selected[selected$selected %in% TRUE, , drop = FALSE]
  lower <- tail(configuration$fit_configuration$fair_common_lambda_relative_grid, 1L)
  expected_boundary <- selected$point_type == "finite" &
    is.finite(selected$lambda_relative_to_reference) &
    abs(selected$lambda_relative_to_reference - lower) <= 1e-10
  expected_boundary[is.na(expected_boundary)] <- FALSE
  matches <- match(a$method, selected$method)
  best <- vapply(allb_methods(), function(method) {
    q <- payload$inner_candidates
    q <- q[q$method == method & q$eligible %in% TRUE &
      is.finite(q$inner_log_loss), , drop = FALSE]
    if (!nrow(q)) return(FALSE)
    q <- q[allb_order_aggregated_candidates_v2(q), , drop = FALSE]
    isTRUE(q$selected[[1L]])
  }, logical(1))
  extra <- c(
    finite_grid_policy_explicit = identical(configuration$fit_configuration$
      fair_lower_boundary_policy, "report_not_reject"),
    lower_boundary_flags_match_selection = !anyNA(matches) &&
      identical(as.logical(a$selected_lower_boundary), expected_boundary[matches]),
    selected_cv_candidates_eligible_and_minimal = all(best) &&
      nrow(selected) == 6L && all(selected$eligible %in% TRUE),
    refit_candidates_numerically_eligible = nrow(payload$selected_rows) == 6L &&
      all(payload$selected_rows$numerically_eligible %in% TRUE),
    range_report_exact = identical(payload$lambda_range_warnings,
      allb_range_report_v8(payload, task, configuration))
  )
  rbind(old[!removed, , drop = FALSE],
    data.frame(check = names(extra), passed = unname(extra)))
}

allb_run_task_v8 <- function(e, data, task, configuration) {
  p <- allb_run_task_v4(e, data, task, configuration)
  p$lambda_range_warnings <- allb_range_report_v8(p, task, configuration)
  p
}

# Rebind only the orchestration functions in an isolated environment. Frozen
# functions and namespaces are not modified; V6 storage/progress schemas persist.
allb_context_v8 <- function() {
  b <- new.env(parent = environment(allb_context_v8))
  originals <- c("allb_valid_shard_v6", "allb_worker_v6",
    "allb_progress_snapshot_v6", "allb_run_pending_v6",
    "allb_finalize_v6", "allb_verify_outputs_v6", "allb_run_v6")
  for (name in originals) {
    f <- get(name, envir = parent.env(b))
    environment(f) <- b
    assign(name, f, envir = b)
  }
  b$allb_load_environment_v6 <- allb_load_environment_v7
  b$allb_make_spec_v6 <- allb_make_spec_v8
  b$allb_validate_spec_v6 <- allb_validate_spec_v8
  b$allb_validate_release_v6 <- allb_validate_release_v8
  b$allb_validate_payload_v6 <- allb_validate_payload_v8
  b$allb_run_task_v6 <- allb_run_task_v8
  b$allb_final_names_v3 <- function() c(allb_final_names_v3(), "lambda_range_warnings")
  b
}

allb_make_spec_v8 <- function(e, root, stage, version) {
  f <- allb_clone_with_bindings_v2(allb_make_spec_v7, list(
    allb_configuration_v7 = allb_configuration_v8,
    allb_source_files_v7 = allb_source_files_v8))
  x <- f(e, root, stage, version)
  x$schema_version <- "all_binary_study_v8"
  x$scientific_signature <- NULL
  x$scientific_signature <- allb_hash(e,
    x[!names(x) %in% c("runtime", "created_utc")])
  x
}

allb_validate_spec_v8 <- function(e, spec, root, stage, version) {
  current <- allb_make_spec_v8(e, root, stage, version)
  keep <- setdiff(names(current), "created_utc")
  allb_assert(identical(spec[keep], current[keep]),
    "V8 identity, sources, data, runtime or configuration changed.")
  invisible(TRUE)
}

allb_validate_release_v8 <- function(e, root, version) {
  receipt <- readRDS(file.path(root, "release/all_binary_v8",
    paste0("LOCAL_VALIDATED_", version, ".rds")))
  allb_assert(isTRUE(receipt$accepted) && all(receipt$checks$passed),
    "V8 local validation is not accepted.")
  allb_validate_spec_v8(e, receipt$identity, root, "production", version)
  invisible(receipt)
}

allb_import_v7_v8 <- function(e, root, spec, parent_version) {
  parent_out <- allb_output_directory(root, parent_version)
  parent_file <- file.path(parent_out, "study_specification.rds")
  parent <- readRDS(parent_file)
  allb_validate_spec_v7(e, parent, root, spec$stage, parent_version)
  expected <- allb_configuration_v8(e, spec$stage)
  allb_assert(identical(spec$configuration, expected) &&
    identical(parent$configuration, allb_configuration_v7(e, spec$stage)) &&
    identical(parent$tasks, spec$tasks) &&
    identical(parent$data_identity, spec$data_identity) &&
    identical(parent$runtime, spec$runtime), "V7 import fitting identity mismatch.")
  out <- allb_output_directory(root, spec$version)
  dir.create(file.path(out, "shards"), recursive = TRUE, showWarnings = FALSE)
  sf <- file.path(out, "study_specification.rds")
  if (file.exists(sf)) {
    spec <- readRDS(sf)
    allb_validate_spec_v8(e, spec, root, spec$stage, spec$version)
  } else allb_atomic_rds(spec, sf)
  imported <- 0L
  for (i in seq_len(nrow(spec$tasks))) {
    task <- spec$tasks[i, , drop = FALSE]
    path <- file.path(parent_out, "shards", task$shard_file)
    if (!file.exists(path)) next
    allb_assert(allb_valid_shard_v6(path, task, parent$scientific_signature,
      parent$configuration), paste("Invalid parent shard:", path))
    s <- readRDS(path)
    s$provenance <- list(operation = "V7_payload_revalidated_under_V8_range_policy",
      parent_version = parent_version, parent_signature = parent$scientific_signature,
      parent_spec_sha256 = digest::digest(file = parent_file, algo = "sha256"),
      parent_shard = normalizePath(path),
      parent_shard_sha256 = digest::digest(file = path, algo = "sha256"),
      model_refits = 0L)
    s$payload$lambda_range_warnings <- allb_range_report_v8(s$payload, task, spec$configuration)
    s$checks <- allb_validate_payload_v8(s$payload, task, spec$configuration)
    allb_assert(all(s$checks$passed), "Imported payload failed V8 checks.")
    s$scientific_signature <- spec$scientific_signature
    allb_atomic_rds(s, file.path(out, "shards", task$shard_file), identical_ok = TRUE)
    imported <- imported + 1L
  }
  b <- allb_context_v8()
  valid <- vapply(seq_len(nrow(spec$tasks)), function(i) {
    b$allb_valid_shard_v6(file.path(out, "shards", spec$tasks$shard_file[i]),
      spec$tasks[i, , drop = FALSE], spec$scientific_signature, spec$configuration)
  }, logical(1))
  status <- ifelse(valid, "passed", "pending")
  snapshot <- b$allb_progress_snapshot_v6(spec, out, status,
    proc.time()[["elapsed"]], "checkpoints_imported",
    paste0("import_", format(Sys.time(), "%Y%m%dT%H%M%S")), 0L)
  allb_write_progress_v3(snapshot, out)
  cat("Revalidated V7 shards:", imported, "\nModel fits: 0\n")
  invisible(imported)
}

allb_report_ranges_v8 <- function(out) {
  x <- readRDS(file.path(out, "final/lambda_range_warnings.rds"))
  counts <- aggregate(cbind(selected_common_lower, path_truncated_any,
    selected_returned_endpoint) ~ method, x, sum)
  counts$replications <- as.integer(table(x$method)[counts$method])
  counts$common_lower_rate <- counts$selected_common_lower / counts$replications
  cat("Finite-grid range warnings and verified truncation counts:\n")
  print(counts, row.names = FALSE)
  invisible(counts)
}
