# Diagnostic-only support for the four failed V15 production tasks.
# This file does not change a production estimator, grid, seed or gate.

lsg_v16_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

lsg_v16_constants <- function() {
  list(
    schema = "log_sglasso_v15_failure_packet_v16_1",
    version = "logistic_sglasso_joint_r50_v15_a01",
    signature = "292b830187e0588a966bb9b4a0e069f9307ac42827dae76e131c063c0415b49e",
    attempt = "20260904T102754.578248_1109015",
    boundary_ids = c(240L, 334L),
    input_ids = c(321L, 350L),
    failures = data.frame(
      task_id = c(240L, 321L, 334L, 350L),
      status = c("gate_failed", "error", "gate_failed", "error"),
      message = c(
        "sglasso_finite_range_resolved",
        "Invalid Group Lasso prefix evidence.",
        "sglasso_finite_range_resolved;adelie_finite_range_resolved",
        "Invalid Group Lasso prefix evidence."
      ),
      stringsAsFactors = FALSE
    )
  )
}

lsg_v16_source_files <- function() {
  c(
    "R/logistic_sglasso_failure_diagnostic_v16.R",
    "scripts/61_export_logistic_sglasso_v15_failures_v16.R",
    "scripts/62_diagnose_logistic_sglasso_v15_failures_v16.R",
    "scripts/63_validate_logistic_sglasso_failure_diagnostic_v16.R",
    "tests/test_failure_diagnostic_v16.R",
    "LOGISTIC_SGLASSO_FAILURE_DIAGNOSTIC_PROTOCOL_V16.md",
    "truba/README_TRUBA_LOGISTIC_SGLASSO_V16_DIAGNOSTIC.md",
    "truba/run_logistic_sglasso_v15_failure_export_v16.slurm",
    "truba/build_logistic_sglasso_failure_diagnostic_v16.R"
  )
}

lsg_v16_load_source <- function(source_root) {
  lsg_v16_assert(length(source_root) == 1L && dir.exists(source_root),
                 "Missing V15 logistic_prework source root.")
  source_root <- normalizePath(source_root, mustWork = TRUE)
  entry <- file.path(source_root, "R/logistic_sglasso_workflow_v15.R")
  lsg_v16_assert(file.exists(entry), "Missing V15 workflow entry point.")
  scope <- new.env(parent = globalenv())
  sys.source(entry, envir = scope)
  lsg_v16_assert(exists("lsg_load_v15", envir = scope, inherits = FALSE),
                 "V15 loader was not defined.")
  e <- scope$lsg_load_v15(source_root)
  required <- c(
    "lsg_validate_spec_v7", "lsg_read_shard_v7", "lsg_with_rng_v7",
    "simulate_logistic_sglasso_two_design_v1", "lsg_tail_validate_data_v6",
    "lsg_data_fingerprint_v7", "lsg_file_hash_v7", "lsg_runtime_v7",
    "validate_penalized_benchmark_grid_v1", "binary_log_loss",
    "lsg_classify_grpreg_path_v13", "lsg_group_lasso_metadata_v13"
  )
  lsg_v16_assert(all(vapply(required, exists, logical(1), envir = e,
                            inherits = FALSE)),
                 "The V15 source environment is incomplete.")
  e
}

lsg_v16_attempt_outcomes <- function(attempt, task_ids) {
  rows <- lapply(task_ids, function(i) {
    item <- attempt$outcomes[[as.character(i)]]
    lsg_v16_assert(is.list(item) && length(item$status) == 1L,
                   paste("Missing recorded outcome for task", i))
    data.frame(
      task_id = as.integer(i),
      status = as.character(item$status),
      message = if (is.null(item$message)) "" else
        paste(as.character(item$message), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

lsg_v16_tuning_data <- function(e, scenario, seed) {
  generated <- e$lsg_with_rng_v7(
    as.integer(seed),
    e$simulate_logistic_sglasso_two_design_v1(scenario, as.integer(seed))
  )
  allowed <- c("X_train", "y_train", "X_validation", "y_validation", "group")
  data <- generated[allowed]
  rm(generated)
  lsg_v16_assert(identical(names(data), allowed),
                 "Tuning-only input schema changed.")
  e$lsg_tail_validate_data_v6(data)
  data
}

lsg_v16_input_fingerprints <- function(e, data) {
  list(
    training = e$lsg_data_fingerprint_v7(list(
      X = data$X_train, y = data$y_train, group = data$group
    )),
    validation = e$lsg_data_fingerprint_v7(list(
      X = data$X_validation, y = data$y_validation
    ))
  )
}

lsg_v16_boundary_evidence <- function(e, output, spec, task) {
  shard <- e$lsg_read_shard_v7(output, spec, task, deep = TRUE)
  failed <- shard$checks$check[!shard$checks$passed]
  selected_methods <- c(
    "Logistic SGLASSO", "Logistic SGLASSO (d=0 boundary)",
    "Logistic Group Elastic Net (adelie)"
  )
  selected_columns <- c(
    "method", "selected_alpha", "selected_d", "selected_lambda",
    "selected_lambda_reference", "selected_lambda_relative_to_reference",
    "selected_lambda_reference_type", "selected_point_type",
    "validation_log_loss", "selected_finite_upper_boundary",
    "selected_finite_lower_boundary"
  )
  selected <- shard$payload$results[
    shard$payload$results$method %in% selected_methods,
    selected_columns, drop = FALSE
  ]
  lsg_v16_assert(nrow(selected) == 3L && !anyDuplicated(selected$method),
                 paste("Invalid selected boundary evidence for task", task$task_id))

  sg_keys <- selected[selected$method %in% selected_methods[1:2],
                      c("selected_alpha", "selected_d"), drop = FALSE]
  sg <- shard$payload$sglasso_tuning
  keep_sg <- rep(FALSE, nrow(sg))
  for (i in seq_len(nrow(sg_keys))) {
    keep_sg <- keep_sg |
      abs(sg$alpha - sg_keys$selected_alpha[i]) <= 1e-12 &
      abs(sg$d - sg_keys$selected_d[i]) <= 1e-12
  }
  sg_columns <- c(
    "candidate_id", "alpha", "d", "lambda_index", "lambda",
    "lambda_reference", "lambda_reference_type",
    "lambda_relative_to_reference", "point_type", "validation_log_loss",
    "numerically_eligible", "selected_free_d", "selected_d0_boundary",
    "converged", "kkt", "finite_coefficient"
  )
  sg_paths <- sg[keep_sg, sg_columns, drop = FALSE]

  et <- shard$payload$external_tuning
  adelie_alpha <- selected$selected_alpha[
    selected$method == "Logistic Group Elastic Net (adelie)"]
  keep_et <- et$method_path == "Logistic Group Elastic Net (adelie)" &
    abs(et$alpha - adelie_alpha) <= 1e-12
  et_columns <- c(
    "candidate_id", "engine", "method_path", "alpha", "lambda_index",
    "lambda", "lambda_reference", "lambda_reference_type",
    "lambda_relative_to_reference", "point_type", "validation_log_loss",
    "numerically_eligible", "selected", "point_converged",
    "finite_coefficient", "finite_probability"
  )
  adelie_path <- et[keep_et, et_columns, drop = FALSE]

  path <- file.path(output, "shards", task$shard_file)
  receipt_path <- paste0(path, ".receipt.rds")
  list(
    task = task,
    failed_checks = failed,
    selected = selected,
    sglasso_paths = sg_paths,
    adelie_path = adelie_path,
    shard_sha256 = e$lsg_file_hash_v7(path),
    receipt_sha256 = e$lsg_file_hash_v7(receipt_path),
    archival_checks = shard$archival_checks,
    replay_checks = shard$checks
  )
}

lsg_v16_build_packet <- function(source_root, version, stage, signature,
                                 attempt_id, expected_outcomes,
                                 boundary_ids, input_ids) {
  lsg_v16_assert(stage %in% c("smoke", "production"), "Invalid source stage.")
  e <- lsg_v16_load_source(source_root)
  output <- file.path(source_root, "outputs", "study", version)
  spec_path <- file.path(output, "study_specification.rds")
  attempt_path <- file.path(output, "attempts", paste0(attempt_id, ".rds"))
  lsg_v16_assert(file.exists(spec_path) && file.exists(attempt_path),
                 "Missing V15 specification or frozen attempt.")
  spec <- readRDS(spec_path)
  e$lsg_validate_spec_v7(spec, source_root, stage, version, resume = TRUE)
  lsg_v16_assert(identical(spec$scientific_signature, signature),
                 "V15 scientific signature mismatch.")
  attempt <- readRDS(attempt_path)
  lsg_v16_assert(identical(attempt$scientific_signature, signature),
                 "Attempt scientific signature mismatch.")
  observed <- lsg_v16_attempt_outcomes(attempt, expected_outcomes$task_id)
  lsg_v16_assert(identical(observed, expected_outcomes),
                 "Recorded V15 failures differ from the frozen diagnostic scope.")
  if (stage == "production") {
    lsg_v16_assert(length(attempt$status) == 400L &&
                     sum(attempt$status == "completed") == 396L &&
                     sum(attempt$status == "failed") == 4L &&
                     !any(attempt$status %in% c("pending", "running")),
                   "The frozen production attempt is not the reported 396/4 state.")
  }
  lsg_v16_assert(all(boundary_ids %in% expected_outcomes$task_id) &&
                   all(input_ids %in% expected_outcomes$task_id),
                 "Diagnostic task IDs are outside the frozen outcomes.")
  tasks <- spec$tasks[match(expected_outcomes$task_id, spec$tasks$task_id), , drop = FALSE]
  lsg_v16_assert(nrow(tasks) == nrow(expected_outcomes) && !anyNA(tasks$task_id),
                 "Frozen diagnostic tasks are absent from the V15 specification.")

  boundaries <- lapply(boundary_ids, function(id) {
    task <- spec$tasks[spec$tasks$task_id == id, , drop = FALSE]
    lsg_v16_boundary_evidence(e, output, spec, task)
  })
  names(boundaries) <- as.character(boundary_ids)

  inputs <- lapply(input_ids, function(id) {
    task <- spec$tasks[spec$tasks$task_id == id, , drop = FALSE]
    scenario <- spec$design[
      spec$design$scenario_index == task$scenario_index, , drop = FALSE]
    first <- lsg_v16_tuning_data(e, scenario, task$seed)
    second <- lsg_v16_tuning_data(e, scenario, task$seed)
    lsg_v16_assert(identical(first, second),
                   paste("Input regeneration was not exact for task", id))
    list(task = task, scenario = scenario, data = first,
         fingerprints = lsg_v16_input_fingerprints(e, first))
  })
  names(inputs) <- as.character(input_ids)

  list(
    schema_version = lsg_v16_constants()$schema,
    source_version = version,
    source_stage = stage,
    source_scientific_signature = signature,
    source_attempt_id = attempt_id,
    source_attempt_sha256 = e$lsg_file_hash_v7(attempt_path),
    source_specification_sha256 = e$lsg_file_hash_v7(spec_path),
    source_runtime = spec$runtime,
    export_runtime = e$lsg_runtime_v7(),
    configuration = spec$configuration,
    task_outcomes = observed,
    boundary_evidence = boundaries,
    group_lasso_inputs = inputs,
    model_fits_performed = 0L,
    evaluation_data_exported = FALSE,
    test_fields_exported = FALSE,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
}

lsg_v16_validate_packet <- function(packet, production = FALSE) {
  lsg_v16_assert(is.list(packet) &&
                   identical(packet$schema_version, lsg_v16_constants()$schema),
                 "Invalid V16 diagnostic packet schema.")
  lsg_v16_assert(identical(packet$model_fits_performed, 0L) &&
                   identical(packet$evaluation_data_exported, FALSE) &&
                   identical(packet$test_fields_exported, FALSE),
                 "The V16 packet must be input/evidence only.")
  lsg_v16_assert(length(packet$boundary_evidence) > 0L &&
                   length(packet$group_lasso_inputs) > 0L,
                 "The V16 packet has no diagnostic cases.")
  allowed <- c("X_train", "y_train", "X_validation", "y_validation", "group")
  for (item in packet$group_lasso_inputs) {
    lsg_v16_assert(identical(names(item$data), allowed),
                   "A forbidden field entered the tuning-data packet.")
    lsg_v16_assert(is.matrix(item$data$X_train) &&
                     is.matrix(item$data$X_validation) &&
                     ncol(item$data$X_train) == length(item$data$group) &&
                     ncol(item$data$X_validation) == length(item$data$group) &&
                     nrow(item$data$X_train) == length(item$data$y_train) &&
                     nrow(item$data$X_validation) == length(item$data$y_validation) &&
                     all(item$data$y_train %in% 0:1) &&
                     all(item$data$y_validation %in% 0:1),
                   "Invalid training/validation dimensions in the V16 packet.")
  }
  if (production) {
    constants <- lsg_v16_constants()
    lsg_v16_assert(identical(packet$source_version, constants$version) &&
                     identical(packet$source_stage, "production") &&
                     identical(packet$source_scientific_signature,
                               constants$signature) &&
                     identical(packet$source_attempt_id, constants$attempt) &&
                     identical(packet$task_outcomes, constants$failures) &&
                     identical(as.integer(names(packet$boundary_evidence)),
                               constants$boundary_ids) &&
                     identical(as.integer(names(packet$group_lasso_inputs)),
                               constants$input_ids),
                   "The packet is not the frozen V15 production failure evidence.")
    for (id in constants$boundary_ids) {
      evidence <- packet$boundary_evidence[[as.character(id)]]
      wanted <- strsplit(constants$failures$message[
        constants$failures$task_id == id], ";", fixed = TRUE)[[1L]]
      lsg_v16_assert(setequal(evidence$failed_checks, wanted),
                     paste("Unexpected failed checks for task", id))
    }
  }
  invisible(TRUE)
}

lsg_v16_publish_packet <- function(packet, output, file_hash) {
  lsg_v16_validate_packet(packet, production = identical(packet$source_stage, "production"))
  lsg_v16_assert(length(output) == 1L && startsWith(output, "/") &&
                   dir.exists(dirname(output)),
                 "Use an absolute new packet path with an existing parent.")
  checksum <- paste0(output, ".sha256")
  lock <- paste0(output, ".lock")
  links <- Sys.readlink(c(output, checksum))
  lsg_v16_assert(!file.exists(output) && !file.exists(checksum) &&
                   all(is.na(links) | !nzchar(links)),
                 "Refusing to overwrite diagnostic packet evidence.")
  lsg_v16_assert(dir.create(lock, showWarnings = FALSE),
                 "Diagnostic packet export lock exists.")
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  temporary <- tempfile(".v16_packet_", tmpdir = dirname(output))
  sidecar <- tempfile(".v16_sha_", tmpdir = dirname(output))
  on.exit(unlink(c(temporary, sidecar)), add = TRUE)
  saveRDS(packet, temporary, compress = "gzip", version = 3)
  restored <- readRDS(temporary)
  lsg_v16_validate_packet(restored,
                          production = identical(packet$source_stage, "production"))
  lsg_v16_assert(identical(packet, restored), "Packet read-back changed its content.")
  hash <- file_hash(temporary)
  writeLines(paste(hash, basename(output), sep = "  "), sidecar, useBytes = TRUE)
  lsg_v16_assert(file.rename(temporary, output) && file.rename(sidecar, checksum),
                 "Atomic packet publication failed.")
  cat("Input/evidence only: TRUE\nModel fits: 0\nTest fields exported: FALSE\n")
  cat("Output:", normalizePath(output), "\nSHA256:", hash, "\n")
  invisible(hash)
}

lsg_v16_read_packet <- function(input, production = TRUE) {
  lsg_v16_assert(file.exists(input) && file.exists(paste0(input, ".sha256")),
                 "Missing packet or SHA-256 sidecar.")
  lsg_v16_assert(requireNamespace("digest", quietly = TRUE),
                 "The installed digest package is required.")
  line <- readLines(paste0(input, ".sha256"), warn = FALSE)
  expected <- strsplit(line, "  ", fixed = TRUE)[[1L]][1L]
  observed <- digest::digest(file = input, algo = "sha256")
  lsg_v16_assert(length(line) == 1L && grepl("^[a-f0-9]{64}  ", line) &&
                   identical(expected, observed), "Packet SHA-256 mismatch.")
  packet <- readRDS(input)
  lsg_v16_validate_packet(packet, production)
  packet
}

lsg_v16_fit_group_lasso_raw <- function(e, data, configuration) {
  lsg_v16_assert(requireNamespace("grpreg", quietly = TRUE),
                 "The installed grpreg package is required.")
  input <- e$validate_penalized_benchmark_grid_v1(
    data$X_train, data$y_train, data$group,
    data$X_validation, data$y_validation,
    alpha_grid = 1, nlambda = configuration$benchmark_nlambda,
    lambda_min_ratio = configuration$benchmark_lambda_min_ratio
  )
  max_iterations <- as.integer(configuration$grpreg_max_iterations)
  tolerance <- configuration$grpreg_tolerance
  group_factor <- factor(input$group, levels = unique(input$group))
  group_multiplier <- sqrt(as.numeric(table(group_factor)))
  captured <- character()
  elapsed <- system.time({
    fit <- withCallingHandlers(
      grpreg::grpreg(
        X = input$X, y = input$y, group = input$group,
        penalty = "grLasso", family = "binomial",
        nlambda = input$nlambda, lambda.min = input$lambda_min_ratio,
        log.lambda = TRUE, alpha = 1, eps = tolerance,
        max.iter = max_iterations, dfmax = ncol(input$X),
        gmax = length(unique(input$group)), gamma = 3,
        group.multiplier = group_multiplier, warn = TRUE, returnX = FALSE
      ),
      warning = function(condition) {
        captured <<- c(captured, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    )
    probability <- as.matrix(stats::predict(
      fit, X = input$X_validation, type = "response"
    ))
  })[["elapsed"]]
  n <- length(fit$lambda)
  lsg_v16_assert(ncol(fit$beta) == n && length(fit$iter) == n &&
                   nrow(probability) == length(input$y_validation) &&
                   ncol(probability) == n,
                 "grpreg returned incompatible diagnostic arrays.")
  status <- e$lsg_classify_grpreg_path_v13(
    "grLasso", fit$iter, input$nlambda, max_iterations, captured
  )
  finite_coefficient <- apply(fit$beta, 2L, function(x) all(is.finite(x)))
  finite_probability <- apply(probability, 2L, function(x) all(is.finite(x)))
  loss <- vapply(seq_len(n), function(i) {
    if (finite_probability[i]) e$binary_log_loss(input$y_validation, probability[, i])
    else NA_real_
  }, numeric(1))
  eligible <- status$path_termination_acceptable & status$point_converged &
    !status$returned_lower_boundary_excluded & finite_coefficient &
    finite_probability & is.finite(loss)
  tuning <- data.frame(
    engine = "grpreg", method_path = "Logistic Group Lasso (grpreg)",
    penalty_family = "group_lasso", fit_index = 1L, alpha_index = 1L,
    alpha = 1, gamma = NA_real_, lambda_index = seq_len(n),
    lambda = fit$lambda, lambda_fraction = fit$lambda / fit$lambda[1L],
    lambda_reference = fit$lambda[1L],
    lambda_reference_type = "native_grpreg_path_start",
    lambda_relative_to_reference = fit$lambda / fit$lambda[1L],
    alpha_zero_ridge_boundary = FALSE, group_selection_capable = TRUE,
    validation_log_loss = loss, passes = as.numeric(fit$iter),
    solver_warning = paste(status$unique_warnings, collapse = " | "),
    unexpected_solver_warning = paste(status$unexpected_warnings, collapse = " | "),
    requested_path_length = input$nlambda, returned_path_length = n,
    path_complete = status$path_complete,
    saturated_path_truncation = status$saturated_path_truncation,
    iteration_budget_truncation = status$iteration_budget_truncation,
    total_iterations = status$total_iterations,
    total_iteration_limit_reached = status$total_iteration_limit_reached,
    returned_lower_boundary_excluded = status$returned_lower_boundary_excluded,
    path_termination_acceptable = status$path_termination_acceptable,
    point_converged = status$point_converged,
    finite_validation_loss = is.finite(loss), numerically_eligible = eligible,
    finite_coefficient = finite_coefficient,
    finite_probability = finite_probability,
    iteration_budget = max_iterations,
    stringsAsFactors = FALSE
  )
  list(tuning = tuning, status = status, warnings = unique(captured),
       elapsed_seconds = elapsed,
       production_metadata_valid = e$lsg_group_lasso_metadata_v13(tuning))
}

lsg_v16_group_lasso_checks <- function(e, fit) {
  z <- fit$tuning[order(fit$tuning$lambda_index), , drop = FALSE]
  status <- fit$status
  n <- nrow(z)
  add <- function(name, passed, detail = "") {
    data.frame(check = name, passed = isTRUE(passed), detail = as.character(detail),
               stringsAsFactors = FALSE)
  }
  fields <- c(
    "path_complete", "saturated_path_truncation", "iteration_budget_truncation",
    "total_iteration_limit_reached", "returned_lower_boundary_excluded",
    "path_termination_acceptable", "point_converged"
  )
  field_match <- all(vapply(fields, function(field) {
    identical(as.logical(z[[field]]), rep_len(as.logical(status[[field]]), n))
  }, logical(1)))
  expected_eligible <- status$point_converged &
    !status$returned_lower_boundary_excluded
  checks <- rbind(
    add("production_metadata_validator", fit$production_metadata_valid),
    add("returned_indices_contiguous",
        identical(as.numeric(z$lambda_index), as.numeric(seq_len(n))),
        paste("returned", n, "requested", z$requested_path_length[1L])),
    add("finite_coefficients_probabilities_losses",
        all(z$finite_coefficient & z$finite_probability &
              is.finite(z$validation_log_loss))),
    add("termination_classified_acceptable", status$path_termination_acceptable,
        paste("complete", status$path_complete,
              "budget_truncated", status$iteration_budget_truncation)),
    add("cumulative_budget_not_reached_before_final_returned_point",
        !isTRUE(status$group_lasso_safe_prefix) || n <= 1L ||
          !any(head(cumsum(z$passes), -1L) >= z$iteration_budget[1L]),
        paste("first_cumulative_budget_index",
          if (any(cumsum(z$passes) >= z$iteration_budget[1L]))
            which(cumsum(z$passes) >= z$iteration_budget[1L])[1L] else "none")),
    add("stored_classification_fields_match", field_match),
    add("total_iterations_equal_sum_of_passes",
        z$total_iterations[1L] == sum(z$passes),
        paste(z$total_iterations[1L], sum(z$passes))),
    add("no_unexpected_solver_warning",
        all(!nzchar(z$unexpected_solver_warning)),
        paste(unique(z$unexpected_solver_warning), collapse = " | ")),
    add("eligibility_matches_prefix_policy",
        identical(as.logical(z$numerically_eligible), expected_eligible))
  )
  rownames(checks) <- NULL
  checks
}

lsg_v16_boundary_tables <- function(packet) {
  selected <- paths <- list()
  k <- 0L
  for (id in names(packet$boundary_evidence)) {
    evidence <- packet$boundary_evidence[[id]]
    selected[[length(selected) + 1L]] <- cbind(
      task_id = as.integer(id), evidence$selected, row.names = NULL
    )
    boundary <- evidence$selected[
      evidence$selected$selected_finite_upper_boundary |
        evidence$selected$selected_finite_lower_boundary, , drop = FALSE]
    for (i in seq_len(nrow(boundary))) {
      row <- boundary[i, ]
      if (row$method %in% c("Logistic SGLASSO",
                            "Logistic SGLASSO (d=0 boundary)")) {
        x <- evidence$sglasso_paths[
          abs(evidence$sglasso_paths$alpha - row$selected_alpha) <= 1e-12 &
            abs(evidence$sglasso_paths$d - row$selected_d) <= 1e-12, , drop = FALSE]
      } else {
        x <- evidence$adelie_path[
          abs(evidence$adelie_path$alpha - row$selected_alpha) <= 1e-12, , drop = FALSE]
      }
      x <- x[order(-x$lambda_relative_to_reference), , drop = FALSE]
      k <- k + 1L
      paths[[k]] <- data.frame(
        task_id = as.integer(id), method = row$method,
        alpha = row$selected_alpha, d = row$selected_d,
        lambda_relative_to_reference = x$lambda_relative_to_reference,
        point_type = x$point_type,
        validation_log_loss = x$validation_log_loss,
        numerically_eligible = x$numerically_eligible,
        stringsAsFactors = FALSE
      )
    }
  }
  list(selected = do.call(rbind, selected), paths = do.call(rbind, paths))
}

lsg_v16_run_diagnostics <- function(source_root, packet) {
  lsg_v16_validate_packet(packet, production = TRUE)
  e <- lsg_v16_load_source(source_root)
  for (id in names(packet$group_lasso_inputs)) {
    item <- packet$group_lasso_inputs[[id]]
    lsg_v16_assert(identical(lsg_v16_input_fingerprints(e, item$data),
                             item$fingerprints),
                   paste("Input fingerprint mismatch for task", id))
  }
  fits <- lapply(names(packet$group_lasso_inputs), function(id) {
    fitted <- lsg_v16_fit_group_lasso_raw(
      e, packet$group_lasso_inputs[[id]]$data, packet$configuration
    )
    list(task_id = as.integer(id), fit = fitted,
         checks = lsg_v16_group_lasso_checks(e, fitted))
  })
  names(fits) <- names(packet$group_lasso_inputs)
  list(boundary = lsg_v16_boundary_tables(packet), group_lasso = fits,
       diagnostic_runtime = e$lsg_runtime_v7())
}
