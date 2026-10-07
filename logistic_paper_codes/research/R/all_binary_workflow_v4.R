# V4 is an overlay on the frozen V3 study: only cumulative grpreg iteration
# budget evidence changes. Data, methods, splits, grids, and seeds are intact.

allb_source_files_v4 <- function() unique(c(
  allb_source_files_v3(),
  "R/all_binary_fair_tuning_v4.R",
  "R/all_binary_workflow_v4.R",
  "R/all_binary_io_v4.R",
  "ALL_BINARY_BUDGET_PREFIX_PROTOCOL_V4.md",
  "scripts/87_run_all_binary_v4.R",
  "scripts/88_validate_all_binary_v4.R",
  "scripts/89_summarize_all_binary_v4.R",
  "scripts/90_preflight_all_binary_v4.R",
  "scripts/91_validate_all_binary_budget_prefix_v4.R",
  "tests/test_all_binary_v4.R",
  "truba/run_all_binary_production_v4.slurm",
  "truba/build_all_binary_bundle_v4.R",
  "truba/README_TRUBA_ALL_BINARY_V4.md"
))

allb_load_environment_v4 <- function(root) {
  e <- allb_load_environment_v3(root)
  source(file.path(root, "R/all_binary_fair_tuning_v4.R"), local = e)
  e$allb_install_fair_tuning_v4(e)
  e
}

allb_fit_configuration_v4 <- function(e, stage) {
  cfg <- allb_fit_configuration_v3(e, stage)
  cfg$fair_grpreg_path_policy <- paste(
    "saturation_all_finite_returned_points_eligible",
    "cumulative_budget_terminal_excluded_if_noncompetitive",
    sep = "_"
  )
  cfg$integration_version <- "all_binary_budget_prefix_v4"
  cfg
}

allb_configuration_v4 <- function(e, stage) {
  cfg <- allb_configuration_v3(e, stage)
  cfg$schema_version <- "all_binary_configuration_v4"
  cfg$fit_configuration <- allb_fit_configuration_v4(e, stage)
  cfg$fair_tuning$grpreg_path_rule <- paste(
    "exact_requested_grid_prefix",
    "saturation_all_finite_returned_points",
    "cumulative_budget_last_point_excluded_if_noncompetitive",
    sep = "_"
  )
  cfg
}

allb_inner_select_v4 <- function(e, X, y, group, inner_fold,
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
    grpreg <- bundle$external$tuning[
      bundle$external$tuning$engine == "grpreg", , drop = FALSE
    ]
    allb_assert(!any(grpreg$invalid_validation_contender %in% TRUE),
      paste("A budget-excluded grpreg terminal is competitive in fold", fold))
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
        truncated & q$lambda_index == q$returned_path_length &
          q$numerically_eligible %in% TRUE
      } else rep(FALSE, nrow(q))
      returned_path_evidence <- if (is_grpreg) {
        budget_exclusion <- q$iteration_budget_truncation %in% TRUE &
          q$lambda_index == q$returned_path_length
        all(q$requested_grid_prefix_aligned %in% TRUE) &&
          all(q$path_termination_acceptable %in% TRUE) &&
          all(!is.na(q$unexpected_solver_warning) &
            !nzchar(q$unexpected_solver_warning)) &&
          identical(as.logical(q$returned_lower_boundary_excluded),
            as.logical(budget_exclusion)) &&
          all(q$returned_path_length >= q$lambda_index) &&
          !any(q$invalid_validation_contender %in% TRUE)
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

allb_run_task_v4 <- function(e, data, task, configuration) {
  runner <- allb_clone_with_bindings_v2(allb_run_task, list(
    gab_inner_select = allb_inner_select_v4
  ))
  payload <- runner(e, data, task, configuration)
  payload$fair_tuning_audit <- allb_fair_tuning_audit_v3(
    e, list(candidate_summary = payload$inner_candidates), configuration
  )
  payload
}

allb_validate_payload_v4 <- function(payload, task, configuration) {
  inherited <- allb_validate_payload_v3(payload, task, configuration)
  audit <- payload$fair_tuning_audit
  rbind(inherited, data.frame(
    check = "budget_excluded_terminal_not_competitive",
    passed = is.data.frame(audit) &&
      !any(audit$invalid_returned_candidate_competitive %in% TRUE),
    stringsAsFactors = FALSE
  ))
}
