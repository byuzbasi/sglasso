# ALL BCR/ABL versus NEG V3 workflow.
#
# V3 changes only the interpretation and audit of grpreg's returned path and
# adds fail-fast execution in the I/O layer.  The data, splits, estimators,
# lambda grids, tuning loss, seeds, and six comparison methods remain frozen.

allb_source_files_v3 <- function() unique(c(
  allb_source_files_v2(),
  "R/all_binary_fair_tuning_v3.R",
  "R/all_binary_workflow_v3.R",
  "R/all_binary_io_v3.R",
  "ALL_BINARY_RETURNED_PATH_PROTOCOL_V3.md",
  "scripts/82_run_all_binary_v3.R",
  "scripts/83_validate_all_binary_v3.R",
  "scripts/84_summarize_all_binary_v3.R",
  "scripts/85_preflight_all_binary_v3.R",
  "scripts/86_validate_all_binary_grpreg_returned_path_v3.R",
  "tests/test_all_binary_v3.R",
  "truba/run_all_binary_production_v3.slurm",
  "truba/build_all_binary_bundle_v3.R",
  "truba/README_TRUBA_ALL_BINARY_V3.md"
))

allb_load_environment_v3 <- function(root) {
  e <- allb_load_environment_v2(root)
  source(file.path(root, "R/all_binary_fair_tuning_v3.R"), local = e)
  e$allb_install_fair_tuning_v3(e)
  e
}

allb_fit_configuration_v3 <- function(e, stage) {
  cfg <- allb_fit_configuration_v2(e, stage)
  cfg$fair_grpreg_path_policy <- paste(
    "requested_grid_exact_returned_prefix",
    "all_finite_returned_points_below_max_iter_eligible",
    "solver_lower_boundary_selection_reported_not_rejected",
    sep = "_"
  )
  cfg$integration_version <- "all_binary_returned_path_v3"
  cfg
}

allb_configuration_v3 <- function(e, stage) {
  configuration <- allb_configuration_v2(e, stage)
  configuration$schema_version <- "all_binary_configuration_v3"
  configuration$fit_configuration <- allb_fit_configuration_v3(e, stage)
  configuration$fair_tuning$grpreg_path_rule <- paste(
    "returned_exact_requested_prefix",
    "finite_converged_points_eligible",
    "solver_boundary_selection_audited",
    sep = "_"
  )
  configuration$fair_tuning$solver_boundary_rule <- paste(
    "truncated_returned_endpoint_is_descriptive",
    "only_true_common_endpoint_0.003125_is_rejected",
    sep = "_"
  )
  configuration$execution_policy <-
    "stop_launching_and_terminate_running_workers_after_first_failure"
  configuration
}

allb_inner_select_v3 <- function(e, X, y, group, inner_fold,
                                  configuration) {
  records <- list()
  diagnostics <- list()
  folds <- sort(unique(inner_fold))
  for (fold in folds) {
    train <- inner_fold != fold
    validation <- !train
    bundle <- gab_fit_bundle(
      e, X[train, , drop = FALSE], y[train],
      X[validation, , drop = FALSE], y[validation], group, configuration
    )
    for (method in gab_methods()) {
      candidate <- gab_candidate_frame(bundle, method)
      candidate$fold <- fold
      candidate$fold_size <- sum(validation)
      records[[length(records) + 1L]] <- candidate
    }
    diagnostics[[length(diagnostics) + 1L]] <- data.frame(
      fold = fold, train_n = sum(train), validation_n = sum(validation),
      train_events = sum(y[train]), validation_events = sum(y[validation]),
      firth_failed_groups = bundle$firth$failed_groups,
      firth_maximum_adjusted_score = max(
        bundle$firth$diagnostics$maximum_adjusted_score
      ),
      firth_maximum_fisher_step = max(
        bundle$firth$diagnostics$maximum_fisher_step
      ),
      firth_replay_error = bundle$firth$replay_error,
      stringsAsFactors = FALSE
    )
    rm(bundle)
    invisible(gc(FALSE))
  }

  tuning <- do.call(e$lsg_bind_tuning_rows_v7, records)
  selections <- list()
  summaries <- list()
  for (method in gab_methods()) {
    z <- tuning[tuning$method == method, , drop = FALSE]
    keys <- unique(z$tuning_key)
    candidate <- lapply(keys, function(key) {
      q <- z[z$tuning_key == key, , drop = FALSE]
      represented <- nrow(q) == length(folds) && !anyDuplicated(q$fold)
      eligible <- represented && all(q$numerically_eligible %in% TRUE) &&
        all(is.finite(q$validation_log_loss))
      template <- q[1L, , drop = FALSE]
      is_grpreg <- all(q$engine %in% "grpreg")
      truncated <- if (is_grpreg) {
        q$path_complete %in% FALSE |
          q$saturated_path_truncation %in% TRUE |
          q$iteration_budget_truncation %in% TRUE
      } else rep(FALSE, nrow(q))
      solver_boundary <- if (is_grpreg) {
        truncated & q$lambda_index == q$returned_path_length
      } else rep(FALSE, nrow(q))
      returned_path_evidence <- if (is_grpreg) {
        all(q$requested_grid_prefix_aligned %in% TRUE) &&
          all(q$path_termination_acceptable %in% TRUE) &&
          all(!is.na(q$unexpected_solver_warning) &
            !nzchar(q$unexpected_solver_warning)) &&
          !any(q$returned_lower_boundary_excluded %in% TRUE) &&
          all(q$returned_path_length >= q$lambda_index)
      } else TRUE
      data.frame(
        method = method, tuning_key = key,
        alpha = template$alpha[[1L]],
        d = if ("d" %in% names(template)) template$d[[1L]] else NA_real_,
        gamma = template$gamma[[1L]],
        point_type = template$point_type[[1L]],
        lambda_index = template$lambda_index[[1L]],
        lambda_relative_to_reference =
          template$lambda_relative_to_reference[[1L]],
        folds_represented = length(unique(q$fold)),
        eligible = eligible,
        grpreg_returned_path_evidence = returned_path_evidence,
        path_truncated_any = any(truncated),
        solver_lower_boundary_fold_count = sum(solver_boundary),
        selected_on_solver_lower_boundary = any(solver_boundary),
        invalid_returned_candidate_competitive = if (is_grpreg) {
          any(q$invalid_validation_contender %in% TRUE)
        } else FALSE,
        minimum_returned_path_length = if (is_grpreg) {
          min(q$returned_path_length)
        } else NA_integer_,
        median_returned_path_length = if (is_grpreg) {
          stats::median(q$returned_path_length)
        } else NA_real_,
        maximum_returned_path_length = if (is_grpreg) {
          max(q$returned_path_length)
        } else NA_integer_,
        inner_log_loss = if (eligible) weighted.mean(
          q$validation_log_loss, q$fold_size
        ) else NA_real_,
        stringsAsFactors = FALSE
      )
    })
    candidate <- do.call(rbind, candidate)
    valid <- candidate[candidate$eligible, , drop = FALSE]
    allb_assert(nrow(valid) > 0L,
      paste("No common eligible inner-CV candidate for", method))
    valid <- valid[allb_order_aggregated_candidates_v2(valid), , drop = FALSE]
    selections[[method]] <- valid[1L, , drop = FALSE]
    candidate$selected <- candidate$tuning_key == valid$tuning_key[[1L]]
    summaries[[method]] <- candidate
  }
  list(
    selections = selections,
    candidate_summary = do.call(rbind, summaries),
    fold_diagnostics = do.call(rbind, diagnostics)
  )
}

allb_fair_tuning_audit_v3 <- function(e, inner, configuration) {
  common <- configuration$fit_configuration$fair_common_lambda_relative_grid
  lower <- tail(common, 1L)
  do.call(rbind, lapply(allb_methods(), function(method) {
    candidates <- inner$candidate_summary[
      inner$candidate_summary$method == method, , drop = FALSE
    ]
    selected <- candidates[candidates$selected, , drop = FALSE]
    allb_assert(nrow(selected) == 1L,
      paste("V3 audit expected one selected candidate for", method))
    finite <- candidates[candidates$point_type == "finite" &
      is.finite(candidates$lambda_relative_to_reference), , drop = FALSE]
    grpreg_method <- grepl("(grpreg)", method, fixed = TRUE)
    eligible_finite <- finite[finite$eligible %in% TRUE, , drop = FALSE]
    coverage_frame <- if (grpreg_method) eligible_finite else finite
    path_count <- if (nrow(finite)) {
      length(unique(allb_candidate_path_key_v2(finite)))
    } else 0L
    selected_ratio <- selected$lambda_relative_to_reference[[1L]]
    selected_lower <- identical(selected$point_type[[1L]], "finite") &&
      is.finite(selected_ratio) && abs(selected_ratio - lower) <= 1e-10
    returned <- candidates$minimum_returned_path_length[
      is.finite(candidates$minimum_returned_path_length)
    ]
    data.frame(
      method = method,
      common_grid_points = length(common),
      finite_candidate_paths = path_count,
      full_common_grid_coverage =
        allb_common_path_coverage_v2(candidates, common),
      common_grid_or_returned_prefix_coverage =
        allb_common_path_or_safe_prefix_v2(
          coverage_frame, common, allow_safe_prefix = grpreg_method
        ),
      recognized_returned_path_evidence =
        all(candidates$grpreg_returned_path_evidence %in% TRUE),
      path_truncated_any = any(candidates$path_truncated_any %in% TRUE),
      invalid_returned_candidate_competitive =
        any(candidates$invalid_returned_candidate_competitive %in% TRUE),
      minimum_returned_path_length = if (length(returned)) {
        min(returned)
      } else NA_integer_,
      median_returned_path_length = if (length(returned)) {
        stats::median(returned)
      } else NA_real_,
      maximum_returned_path_length = if (grpreg_method) {
        max(candidates$maximum_returned_path_length, na.rm = TRUE)
      } else NA_integer_,
      selected_point_type = selected$point_type[[1L]],
      selected_lambda_relative_to_reference = selected_ratio,
      selected_lower_boundary = selected_lower,
      selected_on_solver_lower_boundary =
        selected$selected_on_solver_lower_boundary[[1L]],
      solver_lower_boundary_fold_count =
        selected$solver_lower_boundary_fold_count[[1L]],
      common_grid_sha256 = e$lsg_hash_v7(common),
      stringsAsFactors = FALSE
    )
  }))
}

allb_run_task_v3 <- function(e, data, task, configuration) {
  runner <- allb_run_task
  child <- new.env(parent = environment(runner))
  child$gab_inner_select <- allb_inner_select_v3
  environment(runner) <- child
  payload <- runner(e, data, task, configuration)
  inner <- list(candidate_summary = payload$inner_candidates)
  payload$fair_tuning_audit <- allb_fair_tuning_audit_v3(
    e, inner, configuration
  )
  payload
}

allb_validate_payload_v3 <- function(payload, task, configuration) {
  inherited <- allb_validate_payload(payload, task, configuration)
  audit <- payload$fair_tuning_audit
  common <- configuration$fit_configuration$fair_common_lambda_relative_grid
  checks <- c(
    fair_tuning_audit_schema = is.data.frame(audit) &&
      nrow(audit) == length(allb_methods()) &&
      identical(as.character(audit$method), allb_methods()),
    common_grid_request_and_returned_prefix_policy = is.data.frame(audit) &&
      all(audit$common_grid_or_returned_prefix_coverage %in% TRUE) &&
      all(audit$recognized_returned_path_evidence %in% TRUE) &&
      all(audit$common_grid_points == length(common)),
    common_grid_hash_exact = is.data.frame(audit) &&
      length(unique(audit$common_grid_sha256)) == 1L,
    true_common_lambda_lower_boundary_resolved = is.data.frame(audit) &&
      !any(audit$selected_lower_boundary %in% TRUE),
    solver_lower_boundary_is_diagnostic_only = is.data.frame(audit) &&
      all(!audit$selected_on_solver_lower_boundary |
        audit$path_truncated_any),
    test_fields_absent_from_inner_candidates =
      !any(grepl("(^test_|_test$|test_y|y_test)",
        names(payload$inner_candidates), ignore.case = TRUE))
  )
  rbind(inherited, data.frame(
    check = names(checks), passed = unname(checks),
    stringsAsFactors = FALSE
  ))
}
