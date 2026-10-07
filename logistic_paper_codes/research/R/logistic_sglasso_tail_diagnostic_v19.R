# Focused local validation for the only V18 SGLASSO upper-boundary failure.
# It regenerates training/validation data only and fits 11 short path prefixes;
# it never evaluates a test outcome or launches the production study.

lsg_tail_diagnostic_constants_v19 <- function() {
  list(
    schema = "logistic_sglasso_task234_tail_diagnostic_v19_1",
    version = "logistic_sglasso_task234_tail_diagnostic_v19_a02",
    source_version = "logistic_sglasso_joint_r50_v18_a01",
    source_signature =
      "11974ca4da2b2f45d5a167cc28d40526ff65a34b2b875a60fc3b878094a3c892",
    source_slurm_job = "6325255",
    task_id = 234L,
    scenario = "weak_strong_mixed_rhob_0.3",
    scenario_index = 5L,
    replication = 34L,
    seed = 730310865L,
    shard_file = "shard_0234.rds",
    selected_alpha = 0.2,
    selected_d = 0.4,
    new_ratios = c(8192, 4096),
    replay_ratios = c(8192, 4096, 2048, 1024, 512, 256, 128),
    finite_anchor_ratios = c(2048, 1024, 512, 256, 128),
    anchor_tolerance = 1e-10,
    minimum_winner_gap = 1e-9,
    expected_path_calls = 11L,
    expected_all_rows = 277L,
    expected_new_rows = 222L
  )
}

lsg_tail_fit_prefix_v19 <- function(
    data, preprocess, target_original, configuration, alpha, ratios
) {
  alpha_index <- match(alpha, configuration$alpha_grid)
  d_values <- if (abs(alpha - 1) <= 1e-12) 0 else configuration$d_grid
  reference <- lsg_tail_reference_v6(
    preprocess, data$y_train, target_original, alpha
  )
  fit <- lsg_fit_hybrid_path_v11(
    data$X_train, data$y_train, data$group,
    lambda = reference * ratios, d = d_values, alpha = alpha,
    target_original = target_original, preprocess = preprocess,
    controls = configuration$path_controls, compile = FALSE,
    use_active_set = TRUE, enable_fallback = TRUE,
    warm_start_d = TRUE, lambda_order = "decreasing"
  )
  fit$passes <- fit$block_sweeps + fit$apg_iterations
  probability <- predict_logistic_sglasso(
    fit, data$X_validation, type = "response"
  )
  rows <- vector("list", length(ratios) * length(d_values))
  position <- 0L
  for (di in seq_along(d_values)) {
    for (li in seq_along(ratios)) {
      position <- position + 1L
      p <- probability[, li, di]
      beta <- fit$beta_solver[, li, di]
      coefficient <- fit$coefficients[, li, di]
      rows[[position]] <- data.frame(
        candidate_id = paste("task234", alpha_index, di, ratios[li], sep = ":"),
        alpha_index = as.integer(alpha_index), d_index = as.integer(di),
        alpha = alpha, d = d_values[di], lambda_index = as.integer(li),
        lambda = fit$lambda[li], lambda_reference = reference,
        lambda_reference_type = if (abs(alpha) <= 1e-12) {
          "null_score_ridge_boundary"
        } else {
          "d0_null_kkt"
        },
        lambda_relative_to_reference = ratios[li], point_type = "finite",
        validation_log_loss = binary_log_loss(data$y_validation, p),
        finite_validation_loss = is.finite(binary_log_loss(data$y_validation, p)),
        finite_coefficient = all(is.finite(coefficient)),
        finite_probability = all(is.finite(p)),
        kkt = as.numeric(fit$kkt[li, di]),
        converged = isTRUE(fit$converged[li, di]),
        numerically_eligible = isTRUE(fit$converged[li, di]) &&
          all(is.finite(c(coefficient, p, fit$kkt[li, di]))) &&
          fit$kkt[li, di] <= configuration$path_controls$kkt_tolerance,
        passes = as.numeric(fit$passes[li, di]),
        solver_route = as.character(fit$solver_route[li, di]),
        fallback_used = isTRUE(fit$fallback_used[li, di]),
        stringsAsFactors = FALSE
      )
    }
  }
  list(rows = do.call(rbind, rows), fit = fit, reference = reference)
}

lsg_tail_selection_gap_v19 <- function(pool, selection) {
  eligible <- pool[pool$numerically_eligible %in% TRUE, , drop = FALSE]
  losses <- sort(eligible$validation_log_loss)
  if (length(losses) < 2L) return(NA_real_)
  losses[2L] - selection$validation_log_loss
}

lsg_run_tail_diagnostic_v19 <- function(root) {
  constants <- lsg_tail_diagnostic_constants_v19()
  configuration <- lsg_configuration_v7("production")
  design <- lsg_design_for_stage_v7(root, "production")
  tasks <- lsg_task_grid_v7(design, configuration)
  task <- tasks[tasks$task_id == constants$task_id, , drop = FALSE]
  expected_task <- data.frame(
    task_id = constants$task_id,
    scenario = constants$scenario,
    scenario_index = constants$scenario_index,
    replication = constants$replication,
    seed = constants$seed,
    key = paste0(constants$scenario, "::", constants$replication),
    shard_file = constants$shard_file,
    stringsAsFactors = FALSE
  )
  rownames(task) <- NULL
  task_ok <- nrow(task) == 1L &&
    identical(task[, names(expected_task), drop = FALSE], expected_task)
  lsg_assert_v7(task_ok, "The frozen task-234 identity changed.")

  anchors <- utils::read.csv(
    file.path(root, "config/logistic_sglasso_task234_v18_anchors_v19.csv"),
    stringsAsFactors = FALSE
  )
  anchor_identity <- nrow(anchors) == 7L &&
    all(anchors$source_version == constants$source_version) &&
    all(anchors$source_scientific_signature == constants$source_signature) &&
    all(as.character(anchors$slurm_job_id) == constants$source_slurm_job) &&
    all(anchors$task_id == constants$task_id) &&
    all(anchors$scenario == constants$scenario) &&
    all(anchors$replication == constants$replication) &&
    all(anchors$seed == constants$seed) &&
    all(anchors$provenance == "user_supplied_truba_console_transcript")
  lsg_assert_v7(anchor_identity, "The frozen task-234 transcript anchors changed.")

  scenario <- design[
    design$scenario_index == constants$scenario_index, , drop = FALSE
  ]
  generated <- lsg_with_rng_v7(
    constants$seed,
    simulate_logistic_sglasso_two_design_v1(scenario, constants$seed)
  )
  allowed <- c("X_train", "y_train", "X_validation", "y_validation", "group")
  data <- generated[allowed]
  rm(generated)
  lsg_assert_v7(identical(names(data), allowed),
    "The task-234 diagnostic used a field outside training/validation.")
  lsg_tail_validate_data_v6(data)
  fingerprints <- list(
    training = lsg_data_fingerprint_v7(list(
      X = data$X_train, y = data$y_train, group = data$group
    )),
    validation = lsg_data_fingerprint_v7(list(
      X = data$X_validation, y = data$y_validation
    ))
  )
  firth <- estimate_groupwise_logistic_target(
    data$X_train, data$y_train, data$group, method = "firth",
    max_iterations = configuration$target_max_iterations,
    tolerance = configuration$target_tolerance
  )
  lsg_assert_v7(isTRUE(firth$success) && firth$failed_groups == 0L &&
    all(is.finite(firth$target_original)),
    "The task-234 training-only Firth target failed.")
  preprocess <- prepare_lsg_design(data$X_train, data$group)

  compile_lsg_core(rebuild = FALSE, quiet = TRUE)
  lsg_compile_hybrid_solver_v11(root, rebuild = FALSE, quiet = TRUE)

  # Fail fast: replay the selected alpha path and validate its five V18
  # finite anchors before fitting the ten remaining alpha prefixes.
  selected_prefix <- lsg_tail_fit_prefix_v19(
    data, preprocess, firth$target_original, configuration,
    constants$selected_alpha, constants$replay_ratios
  )
  selected_rows <- selected_prefix$rows[
    abs(selected_prefix$rows$d - constants$selected_d) <= 1e-12,
    , drop = FALSE
  ]
  finite_anchors <- anchors[
    is.finite(anchors$lambda_relative_to_reference) &
      anchors$evidence_role != "nested_d0_old_global_incumbent",
    c("lambda_relative_to_reference", "validation_log_loss"),
    drop = FALSE
  ]
  names(finite_anchors)[2L] <- "validation_log_loss_v18"
  replay <- selected_rows[
    selected_rows$lambda_relative_to_reference %in%
      constants$finite_anchor_ratios,
    c("lambda_relative_to_reference", "validation_log_loss"),
    drop = FALSE
  ]
  names(replay)[2L] <- "validation_log_loss_replay"
  anchor_comparison <- merge(
    finite_anchors, replay, by = "lambda_relative_to_reference", all = TRUE
  )
  anchor_comparison$absolute_error <- abs(
    anchor_comparison$validation_log_loss_replay -
      anchor_comparison$validation_log_loss_v18
  )
  anchor_comparison$within_tolerance <-
    anchor_comparison$absolute_error <= constants$anchor_tolerance
  lsg_assert_v7(
    nrow(anchor_comparison) == length(constants$finite_anchor_ratios) &&
      setequal(anchor_comparison$lambda_relative_to_reference,
               constants$finite_anchor_ratios) &&
      all(anchor_comparison$within_tolerance),
    "The local task-234 replay did not reproduce the V18 finite anchors."
  )
  endpoint <- lsg_tail_limit_fit_v6(
    data, preprocess, firth$target_original,
    constants$selected_alpha, constants$selected_d, configuration
  )
  endpoint_anchor <- anchors[
    anchors$evidence_role == "free_path_analytic_limit", , drop = FALSE
  ]
  endpoint_error <- abs(
    endpoint$point$validation_log_loss - endpoint_anchor$validation_log_loss
  )
  lsg_assert_v7(
    nrow(endpoint_anchor) == 1L && endpoint_error <= constants$anchor_tolerance,
    "The local task-234 analytic endpoint did not reproduce V18."
  )

  paths <- list(selected_prefix$rows)
  remaining_alpha <- configuration$alpha_grid[
    abs(configuration$alpha_grid - constants$selected_alpha) > 1e-12
  ]
  for (alpha in remaining_alpha) {
    paths[[length(paths) + 1L]] <- lsg_tail_fit_prefix_v19(
      data, preprocess, firth$target_original, configuration,
      alpha, constants$new_ratios
    )$rows
  }
  candidates <- do.call(rbind, paths)
  candidates <- candidates[order(
    candidates$alpha_index, candidates$d_index, candidates$lambda_index
  ), , drop = FALSE]
  rownames(candidates) <- NULL
  new_tail <- candidates[
    candidates$lambda_relative_to_reference %in% constants$new_ratios,
    , drop = FALSE
  ]
  candidate_key <- function(alpha_index, d_index, lambda_index, ratio) {
    sprintf("%d:%d:%d:%.17g", alpha_index, d_index, lambda_index, ratio)
  }
  expected_new_keys <- unlist(lapply(
    seq_along(configuration$alpha_grid),
    function(ai) {
      d_count <- if (abs(configuration$alpha_grid[ai] - 1) <= 1e-12) {
        1L
      } else {
        length(configuration$d_grid)
      }
      as.vector(outer(
        seq_len(d_count), constants$new_ratios,
        function(di, ratio) candidate_key(
          ai, di, match(ratio, configuration$lambda_relative_grid), ratio
        )
      ))
    }
  ), use.names = FALSE)
  observed_new_keys <- candidate_key(
    new_tail$alpha_index, new_tail$d_index, new_tail$lambda_index,
    new_tail$lambda_relative_to_reference
  )
  selected_prefix_keys <- candidate_key(
    selected_prefix$rows$alpha_index, selected_prefix$rows$d_index,
    selected_prefix$rows$lambda_index,
    selected_prefix$rows$lambda_relative_to_reference
  )
  expected_selected_prefix_keys <- as.vector(outer(
    seq_along(configuration$d_grid), constants$replay_ratios,
    function(di, ratio) candidate_key(
      match(constants$selected_alpha, configuration$alpha_grid), di,
      match(ratio, constants$replay_ratios), ratio
    )
  ))

  old_free <- selected_rows[
    selected_rows$lambda_relative_to_reference == 2048, , drop = FALSE
  ]
  free_pool <- rbind(new_tail, old_free)
  free <- lsg_select_sglasso_candidates_v7(
    free_pool, configuration$d_grid
  )
  free_boundary <- lsg_selection_boundary_v7(free$row, configuration)
  free_gap <- lsg_tail_selection_gap_v19(free_pool, free)

  nested_anchor <- anchors[
    anchors$evidence_role == "nested_d0_old_global_incumbent", , drop = FALSE
  ]
  alpha_index <- match(nested_anchor$alpha, configuration$alpha_grid)
  alpha_reference <- unique(new_tail$lambda_reference[
    abs(new_tail$alpha - nested_anchor$alpha) <= 1e-12 &
      abs(new_tail$d) <= 1e-12
  ])
  nested_old <- new_tail[1L, , drop = FALSE]
  nested_old$candidate_id <- "task234:nested_d0_v18_incumbent"
  nested_old$alpha_index <- as.integer(alpha_index)
  nested_old$d_index <- 1L
  nested_old$alpha <- nested_anchor$alpha
  nested_old$d <- 0
  nested_index <- which.min(abs(
    configuration$lambda_relative_grid -
      nested_anchor$lambda_relative_to_reference
  ))
  lsg_assert_v7(
    length(nested_index) == 1L &&
      abs(configuration$lambda_relative_grid[nested_index] -
            nested_anchor$lambda_relative_to_reference) <= 1e-13,
    "The V18 nested-d0 incumbent is absent from the V19 finite grid."
  )
  nested_old$lambda_index <- as.integer(nested_index)
  nested_old$lambda_reference <- alpha_reference
  nested_old$lambda_reference_type <- "d0_null_kkt"
  nested_old$lambda_relative_to_reference <-
    nested_anchor$lambda_relative_to_reference
  nested_old$lambda <- alpha_reference * nested_old$lambda_relative_to_reference
  nested_old$validation_log_loss <- nested_anchor$validation_log_loss
  nested_old$finite_validation_loss <- TRUE
  nested_old$finite_coefficient <- TRUE
  nested_old$finite_probability <- TRUE
  nested_old$kkt <- NA_real_
  nested_old$converged <- TRUE
  nested_old$numerically_eligible <- TRUE
  nested_old$passes <- NA_real_
  nested_old$solver_route <- "v18_reported_global_incumbent"
  nested_old$fallback_used <- FALSE
  nested_pool <- rbind(
    new_tail[abs(new_tail$d) <= 1e-12, , drop = FALSE], nested_old
  )
  nested <- lsg_select_sglasso_candidates_v7(nested_pool, 0)
  nested_boundary <- lsg_selection_boundary_v7(nested$row, configuration)
  nested_gap <- lsg_tail_selection_gap_v19(nested_pool, nested)

  checks <- data.frame(
    check = c(
      "frozen_task234_identity", "transcript_anchor_identity",
      "training_validation_only", "training_firth_target_valid",
      "five_finite_anchors_reproduced", "analytic_endpoint_reproduced",
      "eleven_production_prefix_paths", "diagnostic_candidate_count",
      "all_new_tail_candidates_represented",
      "new_tail_exact_key_set_without_duplicates",
      "selected_alpha_exact_seven_ratio_prefix",
      "all_diagnostic_candidates_finite_converged_precise",
      "free_exact_minimum_is_interior", "free_winner_gap_exceeds_tolerance",
      "nested_d0_exact_minimum_is_interior",
      "nested_d0_winner_gap_exceeds_tolerance"
    ),
    passed = c(
      task_ok, anchor_identity, identical(names(data), allowed),
      isTRUE(firth$success) && firth$failed_groups == 0L,
      nrow(anchor_comparison) == 5L && all(anchor_comparison$within_tolerance),
      endpoint_error <= constants$anchor_tolerance,
      length(paths) == constants$expected_path_calls,
      nrow(candidates) == constants$expected_all_rows,
      nrow(new_tail) == constants$expected_new_rows &&
        identical(as.integer(table(factor(
          new_tail$lambda_relative_to_reference,
          levels = constants$new_ratios
        ))), c(111L, 111L)),
      !anyDuplicated(observed_new_keys) &&
        setequal(observed_new_keys, expected_new_keys) &&
        length(observed_new_keys) == length(expected_new_keys),
      !anyDuplicated(selected_prefix_keys) &&
        setequal(selected_prefix_keys, expected_selected_prefix_keys) &&
        length(selected_prefix_keys) == length(expected_selected_prefix_keys),
      all(candidates$finite_validation_loss %in% TRUE) &&
        all(candidates$finite_coefficient %in% TRUE) &&
        all(candidates$finite_probability %in% TRUE) &&
        all(candidates$converged %in% TRUE) &&
        all(candidates$numerically_eligible %in% TRUE) &&
        all(candidates$kkt <= configuration$path_controls$kkt_tolerance),
      !free_boundary$unresolved &&
        free$row$lambda_relative_to_reference !=
          configuration$lambda_upper_multiplier &&
        free$row$lambda_relative_to_reference != configuration$lambda_min_ratio,
      is.finite(free_gap) && free_gap > constants$minimum_winner_gap,
      !nested_boundary$unresolved &&
        nested$row$lambda_relative_to_reference !=
          configuration$lambda_upper_multiplier &&
        nested$row$lambda_relative_to_reference != configuration$lambda_min_ratio,
      is.finite(nested_gap) && nested_gap > constants$minimum_winner_gap
    ),
    stringsAsFactors = FALSE
  )
  selections <- rbind(
    data.frame(
      role = "free_d", free$row,
      runner_up_gap = free_gap,
      upper_boundary = free_boundary$upper,
      lower_boundary = free_boundary$lower,
      stringsAsFactors = FALSE
    ),
    data.frame(
      role = "nested_d0", nested$row,
      runner_up_gap = nested_gap,
      upper_boundary = nested_boundary$upper,
      lower_boundary = nested_boundary$lower,
      stringsAsFactors = FALSE
    )
  )
  rownames(selections) <- NULL
  list(
    schema = constants$schema,
    accepted = all(checks$passed),
    checks = checks,
    constants = constants,
    candidates = candidates,
    new_tail_candidates = new_tail,
    anchor_comparison = anchor_comparison,
    analytic_endpoint = data.frame(
      validation_log_loss_v18 = endpoint_anchor$validation_log_loss,
      validation_log_loss_replay = endpoint$point$validation_log_loss,
      absolute_error = endpoint_error,
      within_tolerance = endpoint_error <= constants$anchor_tolerance
    ),
    selections = selections,
    fingerprints = fingerprints,
    path_calls = length(paths),
    test_fields_used = FALSE,
    source_version = constants$source_version,
    source_signature = constants$source_signature,
    sources = lsg_inventory_v7(root, lsg_source_files_v7()),
    runtime = lsg_runtime_v7(),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
}
