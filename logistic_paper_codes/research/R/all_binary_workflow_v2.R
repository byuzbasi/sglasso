# ALL BCR/ABL versus NEG V2 workflow with method-comparable lambda tuning.
# V1, V19, the installed package, and the processed input remain immutable.

allb_source_files_v2 <- function() unique(c(
  allb_source_files(),
  "R/all_binary_fair_tuning_v2.R",
  "R/all_binary_workflow_v2.R",
  "R/all_binary_io_v2.R",
  "ALL_BINARY_FAIR_TUNING_PROTOCOL_V2.md",
  "scripts/78_run_all_binary_v2.R",
  "scripts/79_validate_all_binary_v2.R",
  "scripts/80_summarize_all_binary_v2.R",
  "scripts/81_preflight_all_binary_v2.R",
  "tests/test_all_binary_v2.R",
  "truba/run_all_binary_production_v2.slurm",
  "truba/build_all_binary_bundle_v2.R",
  "truba/README_TRUBA_ALL_BINARY_V2.md"
))

allb_load_environment_v2 <- function(root) {
  e <- allb_load_environment(root)
  source(file.path(root, "R/all_binary_fair_tuning_v2.R"), local = e)
  e$allb_validate_lambda_grids_v2()
  e$allb_install_fair_tuning_v2(e)
  e
}

allb_fit_configuration_v2 <- function(e, stage) {
  cfg <- allb_fit_configuration(e, stage)
  common <- e$allb_common_lambda_grid_v2()
  sglasso <- e$allb_sglasso_lambda_grid_v2()
  cfg$nlambda <- length(sglasso)
  cfg$lambda_relative_grid <- sglasso
  cfg$lambda_base_nlambda <- length(common)
  cfg$lambda_extension_multipliers <- sglasso[sglasso > 1]
  cfg$lambda_upper_multiplier <- sglasso[[1L]]
  cfg$lambda_min_ratio <- tail(common, 1L)
  cfg$benchmark_nlambda <- length(common)
  cfg$benchmark_lambda_min_ratio <- tail(common, 1L)
  cfg$fair_common_lambda_relative_grid <- common
  cfg$fair_sglasso_upper_relative_grid <- sglasso[sglasso > 1]
  cfg$fair_lambda_reference_policy <-
    "method_native_training_only_lambda_max_then_shared_relative_grid"
  cfg$fair_lower_boundary_policy <-
    "selected_0.003125_blocks_task_and_finalization"
  cfg$fair_grpreg_path_policy <- paste(
    "same_requested_relative_grid_recognized_saturation_safe_prefix",
    "terminal_excluded_and_must_be_noncompetitive",
    sep = "_"
  )
  cfg$selection_policy <- paste(
    "weighted_inner_cv_log_loss_exact_ties",
    "lambda_relative_descending_then_parameters_ascending",
    sep = "_"
  )
  cfg$integration_version <- "all_binary_fair_tuning_v2"
  cfg
}

allb_configuration_v2 <- function(e, stage) {
  configuration <- allb_configuration(e, stage)
  configuration$schema_version <- "all_binary_configuration_v2"
  configuration$fit_configuration <- allb_fit_configuration_v2(e, stage)
  configuration$fair_tuning <- list(
    common_relative_grid = e$allb_common_lambda_grid_v2(),
    sglasso_relative_grid = e$allb_sglasso_lambda_grid_v2(),
    selection_metric = "weighted_mean_inner_validation_log_loss",
    exact_tie_rule =
      "larger_relative_lambda_then_alpha_d_gamma_ascending",
    boundary_rule = "lower_endpoint_selection_is_unresolved_and_rejected",
    grpreg_path_rule = paste(
      "recognized_saturation_safe_prefix_terminal_excluded",
      "competitive_terminal_rejected",
      sep = "_"
    ),
    test_data_used_for_tuning = FALSE
  )
  configuration
}

allb_candidate_path_key_v2 <- function(frame) {
  alpha <- ifelse(is.na(frame$alpha), "NA", sprintf("%.12g", frame$alpha))
  d <- ifelse(is.na(frame$d), "NA", sprintf("%.12g", frame$d))
  gamma <- ifelse(is.na(frame$gamma), "NA", sprintf("%.12g", frame$gamma))
  paste(alpha, d, gamma, sep = "|")
}

allb_common_path_coverage_v2 <- function(frame, common) {
  finite <- frame[frame$point_type == "finite" &
    is.finite(frame$lambda_relative_to_reference), , drop = FALSE]
  if (!nrow(finite)) return(FALSE)
  keys <- allb_candidate_path_key_v2(finite)
  all(vapply(split(finite$lambda_relative_to_reference, keys), function(x) {
    all(vapply(common, function(value) {
      any(abs(x - value) <= 1e-10)
    }, logical(1)))
  }, logical(1)))
}

allb_common_path_or_safe_prefix_v2 <- function(
    frame, common, allow_safe_prefix = FALSE
) {
  finite <- frame[frame$point_type == "finite" &
    is.finite(frame$lambda_relative_to_reference), , drop = FALSE]
  if (!nrow(finite)) return(FALSE)
  keys <- allb_candidate_path_key_v2(finite)
  all(vapply(split(finite$lambda_relative_to_reference, keys), function(x) {
    x <- sort(unique(as.numeric(x[x <= common[[1L]] + 1e-10])),
      decreasing = TRUE)
    if (!length(x)) return(FALSE)
    if (isTRUE(allow_safe_prefix)) {
      length(x) <= length(common) &&
        all(abs(x - head(common, length(x))) <= 1e-10)
    } else {
      all(vapply(common, function(value) {
        any(abs(x - value) <= 1e-10)
      }, logical(1)))
    }
  }, logical(1)))
}

allb_order_aggregated_candidates_v2 <- function(valid) {
  ratio <- valid$lambda_relative_to_reference
  ratio[is.na(ratio)] <- -Inf
  d_order <- valid$d
  d_order[is.na(d_order)] <- 0
  gamma_order <- valid$gamma
  gamma_order[is.na(gamma_order)] <- 0
  order(
    valid$inner_log_loss, -ratio, valid$alpha, d_order, gamma_order,
    valid$point_type != "penalty_limit", valid$lambda_index
  )
}

allb_inner_select_v2 <- function(e, X, y, group, inner_fold,
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
    grpreg_tuning <- bundle$external$tuning[
      bundle$external$tuning$engine == "grpreg", , drop = FALSE
    ]
    allb_assert(
      !any(grpreg_tuning$invalid_validation_contender %in% TRUE),
      paste(
        "A terminal grpreg candidate excluded by the safe-prefix policy",
        "was competitive in inner fold", fold
      )
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
      safe_prefix_evidence <- if (is_grpreg) {
        all(q$requested_grid_prefix_aligned %in% TRUE) &&
          all(q$path_termination_acceptable %in% TRUE) &&
          all(!is.na(q$unexpected_solver_warning) &
            !nzchar(q$unexpected_solver_warning)) &&
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
        folds_represented = length(unique(q$fold)), eligible = eligible,
        grpreg_safe_prefix_evidence = safe_prefix_evidence,
        path_truncated_any = if (is_grpreg) {
          any(q$saturated_path_truncation %in% TRUE |
            q$iteration_budget_truncation %in% TRUE)
        } else FALSE,
        terminal_exclusion_any = if (is_grpreg) {
          any(q$returned_lower_boundary_excluded %in% TRUE)
        } else FALSE,
        invalid_terminal_competitive = if (is_grpreg) {
          any(q$invalid_validation_contender %in% TRUE)
        } else FALSE,
        minimum_returned_path_length = if (is_grpreg) {
          min(q$returned_path_length)
        } else NA_integer_,
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

allb_fair_tuning_audit_v2 <- function(e, inner, configuration) {
  common <- configuration$fit_configuration$fair_common_lambda_relative_grid
  lower <- tail(common, 1L)
  do.call(rbind, lapply(allb_methods(), function(method) {
    candidates <- inner$candidate_summary[
      inner$candidate_summary$method == method, , drop = FALSE
    ]
    selected <- candidates[candidates$selected, , drop = FALSE]
    allb_assert(nrow(selected) == 1L,
      paste("V2 audit expected one selected candidate for", method))
    finite <- candidates[candidates$point_type == "finite" &
      is.finite(candidates$lambda_relative_to_reference), , drop = FALSE]
    grpreg_method <- grepl("(grpreg)", method, fixed = TRUE)
    path_count <- if (nrow(finite)) {
      length(unique(allb_candidate_path_key_v2(finite)))
    } else 0L
    selected_ratio <- selected$lambda_relative_to_reference[[1L]]
    selected_lower <- identical(selected$point_type[[1L]], "finite") &&
      is.finite(selected_ratio) && abs(selected_ratio - lower) <= 1e-10
    data.frame(
      method = method,
      common_grid_points = length(common),
      finite_candidate_paths = path_count,
      all_common_path_coverage =
        allb_common_path_coverage_v2(candidates, common),
      common_grid_or_safe_prefix_coverage =
        allb_common_path_or_safe_prefix_v2(
          candidates, common, allow_safe_prefix = grpreg_method
        ),
      recognized_safe_prefix_evidence =
        all(candidates$grpreg_safe_prefix_evidence %in% TRUE),
      path_truncated_any = any(candidates$path_truncated_any %in% TRUE),
      terminal_candidate_competitive =
        any(candidates$invalid_terminal_competitive %in% TRUE),
      minimum_returned_path_length = if (grpreg_method) {
        min(candidates$minimum_returned_path_length, na.rm = TRUE)
      } else NA_integer_,
      maximum_returned_path_length = if (grpreg_method) {
        max(candidates$maximum_returned_path_length, na.rm = TRUE)
      } else NA_integer_,
      selected_point_type = selected$point_type[[1L]],
      selected_lambda_relative_to_reference = selected_ratio,
      selected_lower_boundary = selected_lower,
      common_grid_sha256 = e$lsg_hash_v7(common),
      stringsAsFactors = FALSE
    )
  }))
}

allb_run_task_v2 <- function(e, data, task, configuration) {
  runner <- allb_run_task
  child <- new.env(parent = environment(runner))
  child$gab_inner_select <- allb_inner_select_v2
  environment(runner) <- child
  payload <- runner(e, data, task, configuration)
  inner <- list(
    candidate_summary = payload$inner_candidates,
    selections = lapply(allb_methods(), function(method) {
      payload$inner_candidates[
        payload$inner_candidates$method == method &
          payload$inner_candidates$selected, , drop = FALSE
      ]
    })
  )
  names(inner$selections) <- allb_methods()
  payload$fair_tuning_audit <- allb_fair_tuning_audit_v2(
    e, inner, configuration
  )
  payload
}

allb_validate_payload_v2 <- function(payload, task, configuration) {
  inherited <- allb_validate_payload(payload, task, configuration)
  audit <- payload$fair_tuning_audit
  common <- configuration$fit_configuration$fair_common_lambda_relative_grid
  checks <- c(
    fair_tuning_audit_schema = is.data.frame(audit) &&
      nrow(audit) == length(allb_methods()) &&
      identical(as.character(audit$method), allb_methods()),
    common_grid_request_and_safe_prefix_policy = is.data.frame(audit) &&
      all(audit$common_grid_or_safe_prefix_coverage %in% TRUE) &&
      all(audit$recognized_safe_prefix_evidence %in% TRUE) &&
      !any(audit$terminal_candidate_competitive %in% TRUE) &&
      all(audit$common_grid_points == length(common)),
    common_grid_hash_exact = is.data.frame(audit) &&
      length(unique(audit$common_grid_sha256)) == 1L,
    selected_lambda_lower_boundary_resolved = is.data.frame(audit) &&
      !any(audit$selected_lower_boundary %in% TRUE),
    test_fields_absent_from_inner_candidates =
      !any(grepl("(^test_|_test$|test_y|y_test)",
        names(payload$inner_candidates), ignore.case = TRUE))
  )
  rbind(inherited, data.frame(
    check = names(checks), passed = unname(checks),
    stringsAsFactors = FALSE
  ))
}
