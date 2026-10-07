# One external-cohort pCR evaluation using the frozen ALL/V19 six-method layer.
# A work unit is one inner-CV fold or the final full-development refit.

gse_assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  invisible(TRUE)
}

gse_atomic_rds <- function(object, path) {
  gse_assert(!file.exists(path), paste("Refusing overwrite:", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp <- tempfile("gse25066_atomic_", tmpdir = dirname(path))
  on.exit(if (file.exists(temp)) unlink(temp), add = TRUE)
  saveRDS(object, temp)
  gse_assert(file.rename(temp, path), paste("Atomic save failed:", path))
  invisible(path)
}

gse_atomic_rds_or_identical <- function(object, path) {
  if (file.exists(path)) {
    gse_assert(identical(readRDS(path), object),
      paste("Existing finalized output differs; refusing overwrite:", path))
  } else gse_atomic_rds(object, path)
  invisible(path)
}

gse_processed_path <- function(root) file.path(root, "data_processed",
  "gse25066_pcr_three_probe_v1.rds")

gse_load_data <- function(root, stage) {
  path <- gse_processed_path(root)
  gse_assert(file.exists(path), "Missing versioned processed GSE25066 data")
  x <- readRDS(path)
  gse_assert(identical(x$schema_version, "gse25066_pcr_three_probe_v1") &&
    identical(dim(x$train$X), c(306L, 3486L)) &&
    identical(dim(x$test$X), c(182L, 3486L)) &&
    identical(colnames(x$train$X), colnames(x$test$X)) &&
    sum(x$train$y) == 57L && sum(x$test$y) == 42L &&
    length(x$group_name) == 1162L &&
    identical(x$group, rep(seq_len(1162L), each = 3L)) &&
    !anyDuplicated(c(x$train$sample_id, x$test$sample_id)) &&
    all(is.finite(x$train$X)) && all(is.finite(x$test$X)),
    "GSE25066 processed identity or values changed")
  if (stage == "smoke") {
    # Explicit execution test only: 20+20 development and 10+10 external.
    keep_train <- unlist(lapply(0:1, function(y)
      head(which(x$train$y == y), 20L)), use.names = FALSE)
    keep_test <- unlist(lapply(0:1, function(y)
      head(which(x$test$y == y), 10L)), use.names = FALSE)
    columns <- seq_len(24L)
    x$train$X <- x$train$X[keep_train, columns, drop = FALSE]
    x$train$y <- x$train$y[keep_train]
    x$train$sample_id <- x$train$sample_id[keep_train]
    x$train$source <- x$train$source[keep_train]
    x$test$X <- x$test$X[keep_test, columns, drop = FALSE]
    x$test$y <- x$test$y[keep_test]
    x$test$sample_id <- x$test$sample_id[keep_test]
    x$test$source <- x$test$source[keep_test]
    x$group <- x$group[columns]
    x$group_name <- x$group_name[seq_len(8L)]
    x$probe_id <- x$probe_id[columns]
  }
  x$processed_sha256 <- digest::digest(file = path, algo = "sha256")
  x
}

gse_fold_assignment <- function(y, folds, seed) {
  set.seed(as.integer(seed))
  assignment <- integer(length(y))
  for (value in 0:1) {
    index <- sample(which(y == value), replace = FALSE)
    assignment[index] <- rep(seq_len(folds), length.out = length(index))
  }
  gse_assert(all(vapply(seq_len(folds), function(k) {
    all(0:1 %in% y[assignment == k]) &&
      all(0:1 %in% y[assignment != k])
  }, logical(1))), "An inner fold lacks an outcome class")
  assignment
}

gse_source_inventory <- function(root) {
  own <- c("R/gse25066_binary_workflow_v1.R",
    "scripts/106_prepare_gse25066_binary_v1.R",
    "scripts/107_run_gse25066_binary_v1.R",
    "tests/test_gse25066_binary_v1.R",
    "GSE25066_BINARY_PROTOCOL_V1.md")
  inherited <- c(
    list.files(file.path(root, "R"), pattern =
      "^(all_binary_|genathum_binary_|logistic_sglasso_).*\\.R$",
      full.names = FALSE),
    "src/logistic_sglasso_block_kernel_local_v3.hpp",
    "src/logistic_sglasso_hybrid_solver_local_apg_v4.cpp")
  inherited <- c(paste0("R/", inherited[grepl("[.]R$", inherited)]),
    inherited[!grepl("[.]R$", inherited)])
  files <- sort(unique(c(own, inherited)), method = "radix")
  paths <- file.path(root, files)
  gse_assert(all(file.exists(paths)), "GSE25066 source inventory has missing files")
  data.frame(file = files, bytes = unname(file.info(paths)$size),
    sha256 = unname(vapply(paths, digest::digest, character(1),
      file = TRUE, algo = "sha256")), stringsAsFactors = FALSE)
}

gse_specification <- function(root, e, stage, version) {
  gse_assert(stage %in% c("smoke", "production"), "Invalid stage")
  data <- gse_load_data(root, stage)
  cfg <- allb_configuration_v8(e, stage)
  cfg$schema_version <- "gse25066_binary_configuration_v1"
  cfg$endpoint_id <- "observed_pcr_vs_rd"
  cfg$repeats <- 1L
  cfg$inner_folds <- if (stage == "smoke") 2L else 5L
  cfg$seed_base <- if (stage == "smoke") 101260101L else 101260102L
  cfg$response_threshold <- 0.5
  fold <- gse_fold_assignment(data$train$y, cfg$inner_folds, cfg$seed_base)
  identity <- list(schema_version = "gse25066_external_study_v1",
    stage = stage, version = version, configuration = cfg,
    fold = fold,
    data_identity = list(processed_sha256 = data$processed_sha256,
      train_sample_id = data$train$sample_id,
      test_sample_id = data$test$sample_id,
      train_response_sha256 = digest::digest(data$train$y, algo = "sha256"),
      test_response_sha256 = digest::digest(data$test$y, algo = "sha256"),
      group_sha256 = digest::digest(data$group, algo = "sha256")),
    source_inventory = gse_source_inventory(root),
    runtime = list(r = R.version.string, platform = R.version$platform,
      package_versions = vapply(c("Rcpp", "RcppArmadillo", "digest",
        "adelie", "grpreg", "logistf", "mltools", "jsonlite"), function(p)
          as.character(utils::packageVersion(p)), character(1))))
  identity$scientific_signature <- digest::digest(identity, algo = "sha256")
  identity
}

gse_output <- function(root, version) file.path(root, "outputs", "study", version)
gse_shard <- function(out, unit) file.path(out, "shards", paste0(unit, ".rds"))

gse_valid_shard <- function(path, spec, unit) {
  if (!file.exists(path)) return(FALSE)
  tryCatch({
    x <- readRDS(path)
    identical(x$schema_version, "gse25066_shard_v1") &&
      identical(x$scientific_signature, spec$scientific_signature) &&
      identical(x$unit, unit) && isTRUE(x$validated) &&
      if (startsWith(unit, "fold_")) {
        p <- x$payload
        identical(as.integer(p$fold), as.integer(sub("fold_", "", unit))) &&
          identical(sort(names(p$candidates)), sort(gab_methods())) &&
          all(vapply(p$candidates, function(y) is.data.frame(y) &&
            nrow(y) > 0L, logical(1))) &&
          p$firth$failed_groups == 0L
      } else {
        identical(unit, "refit") &&
          identical(sort(names(x$payload$models)), sort(gab_methods())) &&
          nrow(x$payload$selected_rows) == 6L
      }
  }, error = function(error) FALSE)
}

gse_fit_fold <- function(e, data, spec, k) {
  train <- spec$fold != k
  validation <- !train
  bundle <- gab_fit_bundle(e, data$train$X[train, , drop = FALSE],
    data$train$y[train], data$train$X[validation, , drop = FALSE],
    data$train$y[validation], data$group,
    spec$configuration$fit_configuration)
  tuning <- lapply(gab_methods(), function(method) {
    frame <- gab_candidate_frame(bundle, method)
    frame$fold <- k
    frame$fold_size <- sum(validation)
    frame
  })
  names(tuning) <- gab_methods()
  gse_assert(!any(bundle$external$tuning$invalid_validation_contender %in% TRUE),
    "A grpreg excluded candidate was validation-competitive")
  gse_assert(bundle$firth$failed_groups == 0L &&
    all(bundle$firth$diagnostics$success), "Firth target failed")
  list(fold = as.integer(k), candidates = tuning,
    firth = list(failed_groups = bundle$firth$failed_groups,
      diagnostics = bundle$firth$diagnostics,
      replay_error = bundle$firth$replay_error),
    train_sample_id = data$train$sample_id[train],
    validation_sample_id = data$train$sample_id[validation],
    train_events = sum(data$train$y[train]),
    validation_events = sum(data$train$y[validation]))
}

# Reuse the existing CV selector, feeding it validated checkpointed fold paths.
# Its candidate scoring and exact tie rule are not reimplemented here.
gse_select_from_folds <- function(e, data, spec, out) {
  stored <- lapply(seq_len(spec$configuration$inner_folds), function(k)
    readRDS(gse_shard(out, paste0("fold_", k)))$payload)
  index <- 0L
  cached <- function(e, X_train, y_train, X_validation, y_validation,
                     group, configuration) {
    index <<- index + 1L
    p <- stored[[index]]
    gse_assert(identical(rownames(X_train), p$train_sample_id) &&
      identical(rownames(X_validation), p$validation_sample_id) &&
      identical(group, data$group), "CV checkpoint/sample mismatch")
    # allb_inner_select_v4 calls gab_candidate_frame on these two tuning tables.
    list(sglasso = list(tuning = p$candidates[[gab_methods()[1L]]]),
      external = list(tuning = do.call(e$lsg_bind_tuning_rows_v7,
        p$candidates[gab_methods()[-c(1L, 2L)]])),
      firth = p$firth)
  }
  selector <- allb_clone_with_bindings_v2(allb_inner_select_v4,
    list(gab_fit_bundle = cached))
  selector(e, data$train$X, data$train$y, data$group, spec$fold,
    spec$configuration$fit_configuration)
}

gse_fit_refit <- function(e, data, spec, selected) {
  X_train <- data$train$X
  y_train <- data$train$y
  X_test <- data$test$X
  bundle <- gab_fit_bundle(e, X_train, y_train, X_train, y_train,
    data$group, spec$configuration$fit_configuration)
  interface <- list(X_train = X_train, y_train = y_train,
    X_validation = X_train, y_validation = y_train,
    X_test = X_test, group = data$group)
  models <- selected_rows <- list()
  for (method in gab_methods()) {
    key <- selected$selections[[method]]$tuning_key[[1L]]
    candidates <- gab_candidate_frame(bundle, method)
    row <- candidates[candidates$tuning_key == key, , drop = FALSE]
    gse_assert(nrow(row) == 1L && row$numerically_eligible[[1L]] %in% TRUE,
      paste("Selected candidate could not be refit for", method))
    if (grepl("SGLASSO", method, fixed = TRUE)) {
      model <- e$lsg_selected_sglasso_model_v7(bundle$sglasso,
        list(row = row), interface)
    } else if (identical(method, "Logistic Group Elastic Net (adelie)")) {
      if (identical(row$point_type[[1L]], "penalty_limit")) {
        # The analytic alpha=0 endpoint has no finite coefficient-path index.
        # Match lsg_selected_external_model_v7's original-scale solution.
        endpoint <- bundle$external$adelie$endpoint
        model <- list(coefficient = as.numeric(endpoint$coefficients),
          probability = as.numeric(stats::plogis(endpoint$intercept_solver +
            drop(e$transform_lsg_newx(bundle$external$adelie$preprocess,
              X_test) %*% endpoint$beta_solver))))
      } else {
        grid <- bundle$external$adelie
        grid$selection$row <- row
        solution <- e$extract_adelie_group_en_solution_v2(grid, X_test)
        model <- list(coefficient = solution$coefficient,
          probability = solution$probability)
      }
    } else {
      grid <- bundle$external$grpreg
      penalty <- grid$penalty_specification$penalty[
        grid$penalty_specification$method_path == method]
      grid$selections[[penalty]]$row <- row
      solution <- e$extract_grpreg_group_solution_v2(grid, penalty, X_test)
      model <- list(coefficient = solution$coefficient,
        probability = solution$probability)
    }
    reconstruction <- max(abs(stats::plogis(model$coefficient[1L] +
      drop(X_test %*% model$coefficient[-1L])) - model$probability))
    gse_assert(length(model$coefficient) == ncol(X_train) + 1L &&
      all(is.finite(model$coefficient)) &&
      length(model$probability) == nrow(X_test) &&
      all(is.finite(model$probability)) &&
      all(model$probability >= 0 & model$probability <= 1) &&
      is.finite(reconstruction) && reconstruction <= 2e-6,
      paste("Invalid refitted model for", method))
    models[[method]] <- model
    row$method <- method
    selected_rows[[method]] <- row
  }
  result <- list(models = models,
    selected_rows = do.call(e$lsg_bind_tuning_rows_v7, selected_rows),
    firth = bundle$firth)
  gse_assert(all(result$firth$diagnostics$success) &&
    nrow(result$selected_rows) == 6L, "Full-development refit invalid")
  result
}

gse_progress <- function(out, spec, status, phase, started, current = NULL) {
  units <- c(paste0("fold_", seq_len(spec$configuration$inner_folds)), "refit")
  valid <- vapply(units, function(unit)
    gse_valid_shard(gse_shard(out, unit), spec, unit), logical(1))
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  completed <- sum(valid)
  completed_this_invocation <- sum(status[units] == "passed" & valid)
  rate <- if (elapsed > 0) completed_this_invocation / elapsed else 0
  running <- if (!is.null(current) && !valid[[match(current, units)]]) 1L else 0L
  fold_units <- head(units, -1L)
  finished_folds <- fold_units[valid[fold_units]]
  durations <- vapply(finished_folds, function(unit) {
    value <- readRDS(gse_shard(out, unit))$runtime_seconds
    if (is.numeric(value) && length(value) == 1L && is.finite(value))
      value else NA_real_
  }, numeric(1))
  remaining_folds <- sum(!valid[fold_units])
  eta <- if (length(durations) >= 2L && all(is.finite(durations)) &&
    remaining_folds > 0L) mean(durations) * remaining_folds else NULL
  now <- Sys.time()
  slurm_end <- suppressWarnings(as.numeric(Sys.getenv("SLURM_JOB_END_TIME", "")))
  snapshot <- list(schema_version = "gse25066_progress_v1",
    run_id = spec$version, slurm_job_id = Sys.getenv("SLURM_JOB_ID", "local"),
    phase = phase, work_unit = "cv_fold_or_full_development_refit",
    total = length(units), completed = completed,
    running = running,
    failed = sum(status[units] == "failed"),
    pending = length(units) - completed -
      running -
      sum(status[units] == "failed"),
    percent_completed = 100 * completed / length(units),
    elapsed_seconds = elapsed,
    throughput_units_per_second = rate,
    eta_seconds = eta,
    estimated_completion_utc = if (!is.null(eta))
      format(now + eta, tz = "UTC", usetz = TRUE) else NULL,
    eta_status = if (!is.null(eta)) "measured_cv_fold_durations" else
      "estimating_or_unknown",
    eta_scope = "remaining_cv_folds_only_excludes_refit_and_final_aggregation",
    slurm_remaining_seconds = if (is.finite(slurm_end))
      max(0, slurm_end - as.numeric(now)) else NULL,
    heartbeat_utc = format(now, tz = "UTC", usetz = TRUE),
    most_recently_completed_work_unit = if (completed)
      tail(units[valid], 1L) else NULL,
    current_work_unit = current)
  temp <- tempfile("gse_progress_", tmpdir = out)
  jsonlite::write_json(snapshot, temp, auto_unbox = TRUE, pretty = TRUE,
    null = "null")
  gse_assert(file.rename(temp, file.path(out, "progress.json")),
    "Atomic progress write failed")
  row <- data.frame(timestamp_utc = snapshot$heartbeat_utc,
    phase = phase, completed = completed, total = length(units),
    running = snapshot$running, failed = snapshot$failed,
    pending = snapshot$pending, percent_completed = snapshot$percent_completed,
    elapsed_seconds = elapsed, throughput_units_per_second = rate,
    eta_seconds = if (is.null(eta)) NA_real_ else eta,
    slurm_remaining_seconds = if (is.null(snapshot$slurm_remaining_seconds))
      NA_real_ else snapshot$slurm_remaining_seconds)
  history <- file.path(out, "progress.tsv")
  utils::write.table(row, history, sep = "\t", row.names = FALSE,
    quote = FALSE, col.names = !file.exists(history), append = file.exists(history))
  cat(sprintf("[%s] %s: %d/%d validated, running=%d failed=%d pending=%d; ETA=%s; SLURM remaining=%s\n",
    snapshot$heartbeat_utc, phase, completed, length(units),
    snapshot$running, snapshot$failed, snapshot$pending,
    if (is.null(eta)) "estimating" else paste0(round(eta), "s"),
    if (is.null(snapshot$slurm_remaining_seconds)) "unavailable" else
      paste0(round(snapshot$slurm_remaining_seconds), "s")))
  flush.console()
  invisible(snapshot)
}

gse_run_unit <- function(unit, out, spec, work, status, started,
                         interval_seconds = 30) {
  path <- gse_shard(out, unit)
  if (gse_valid_shard(path, spec, unit)) return(invisible("resumed"))
  gse_assert(!file.exists(path), paste("Invalid existing shard; no overwrite:", path))
  status[[unit]] <- "running"
  gse_progress(out, spec, status, "fitting", started, unit)
  child <- parallel::mcparallel({
    tryCatch({
      unit_started <- Sys.time()
      payload <- work()
      gse_atomic_rds(list(schema_version = "gse25066_shard_v1",
        scientific_signature = spec$scientific_signature, unit = unit,
        validated = TRUE, payload = payload,
        runtime_seconds = as.numeric(difftime(Sys.time(), unit_started,
          units = "secs")),
        completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)), path)
      list(status = "passed", message = "")
    }, error = function(error) list(status = "failed",
      message = conditionMessage(error)))
  }, silent = TRUE)
  result <- NULL
  repeat {
    collected <- parallel::mccollect(child, wait = FALSE)
    if (!is.null(collected)) {
      result <- collected[[1L]]
      break
    }
    Sys.sleep(interval_seconds)
    gse_progress(out, spec, status, "fitting", started, unit)
  }
  if (!is.list(result) || !identical(result$status, "passed") ||
      !gse_valid_shard(path, spec, unit)) {
    status[[unit]] <- "failed"
    gse_progress(out, spec, status, "failed", started)
    stop("Work unit failed; valid checkpoints preserved: ", unit, "; ",
      if (is.list(result)) result$message else "worker ended unexpectedly",
      call. = FALSE)
  }
  status[[unit]] <- "passed"
  gse_progress(out, spec, status, "unit_validated", started)
  invisible("passed")
}

gse_finalize <- function(e, data, spec, out, selected) {
  gse_assert(!file.exists(file.path(out, "COMPLETE.rds")),
    "Completed output already exists; use verify")
  refit <- readRDS(gse_shard(out, "refit"))$payload
  rows <- predictions <- supports <- list()
  for (method in gab_methods()) {
    model <- refit$models[[method]]
    metric <- gab_response_metrics(e, data$test$y, model$probability, 0.5)
    chosen <- selected$selections[[method]]
    active <- abs(model$coefficient[-1L]) > 1e-8
    selected_groups <- vapply(seq_along(data$group_name), function(g)
      any(active[data$group == g]), logical(1))
    rows[[method]] <- cbind(data.frame(method = method,
      selected_alpha = chosen$alpha[[1L]], selected_d = chosen$d[[1L]],
      selected_gamma = chosen$gamma[[1L]],
      selected_lambda_ratio = chosen$lambda_relative_to_reference[[1L]],
      selected_lambda_index = chosen$lambda_index[[1L]],
      inner_cv_log_loss = chosen$inner_log_loss[[1L]],
      train_n = nrow(data$train$X), test_n = nrow(data$test$X),
      train_events = sum(data$train$y), test_events = sum(data$test$y),
      selected_predictors = sum(active), selected_groups = sum(selected_groups)),
      metric)
    predictions[[method]] <- data.frame(method = method,
      sample_id = data$test$sample_id, y = data$test$y,
      probability = model$probability)
    supports[[method]] <- data.frame(method = method,
      group_id = data$group_name, selected = selected_groups)
  }
  results <- do.call(rbind, rows)
  gse_assert(nrow(results) == 6L && all(is.finite(results$log_loss)) &&
    all(is.finite(results$auc)) &&
    nrow(do.call(rbind, predictions)) == 6L * nrow(data$test$X),
    "External-test output validation failed")
  final <- file.path(out, "final")
  dir.create(final, recursive = TRUE, showWarnings = FALSE)
  outputs <- list(results = results,
    predictions = do.call(rbind, predictions),
    selected_groups = do.call(rbind, supports),
    cv_candidates = selected$candidate_summary,
    cv_fold_diagnostics = selected$fold_diagnostics,
    selected_refit_rows = refit$selected_rows,
    firth_refit_diagnostics = refit$firth$diagnostics)
  for (name in names(outputs))
    gse_atomic_rds_or_identical(outputs[[name]],
      file.path(final, paste0(name, ".rds")))
  made <- file.path(final, paste0(names(outputs), ".rds"))
  manifest <- data.frame(file = basename(made), bytes = file.info(made)$size,
    sha256 = unname(vapply(made, digest::digest, character(1),
      file = TRUE, algo = "sha256")))
  gse_atomic_rds_or_identical(manifest, file.path(final, "manifest.rds"))
  gse_atomic_rds(list(schema_version = "gse25066_complete_v1",
    scientific_signature = spec$scientific_signature,
    validated_units = spec$configuration$inner_folds + 1L,
    method_rows = nrow(results), manifest_sha256 = digest::digest(
      file = file.path(final, "manifest.rds"), algo = "sha256")),
    file.path(out, "COMPLETE.rds"))
  invisible(results)
}

gse_verify <- function(root, e, spec) {
  out <- gse_output(root, spec$version)
  stored <- readRDS(file.path(out, "specification.rds"))
  gse_assert(identical(stored, spec), "GSE25066 scientific identity changed")
  units <- c(paste0("fold_", seq_len(spec$configuration$inner_folds)), "refit")
  gse_assert(all(vapply(units, function(unit)
    gse_valid_shard(gse_shard(out, unit), spec, unit), logical(1))),
    "Missing or invalid fold/refit checkpoint")
  complete <- readRDS(file.path(out, "COMPLETE.rds"))
  final <- file.path(out, "final")
  manifest_path <- file.path(final, "manifest.rds")
  manifest <- readRDS(manifest_path)
  paths <- file.path(final, manifest$file)
  expected_files <- paste0(c("results", "predictions", "selected_groups",
    "cv_candidates", "cv_fold_diagnostics", "selected_refit_rows",
    "firth_refit_diagnostics"), ".rds")
  gse_assert(identical(complete$scientific_signature, spec$scientific_signature) &&
    complete$validated_units == length(units) &&
    complete$method_rows == 6L &&
    identical(complete$manifest_sha256,
      digest::digest(file = manifest_path, algo = "sha256")) &&
    identical(as.character(manifest$file), expected_files) &&
    all(file.exists(paths)) &&
    identical(unname(file.info(paths)$size), unname(manifest$bytes)) &&
    identical(unname(vapply(paths, digest::digest, character(1),
      file = TRUE, algo = "sha256")), as.character(manifest$sha256)),
    "Final manifest/checksum/completeness validation failed")
  results <- readRDS(file.path(final, "results.rds"))
  predictions <- readRDS(file.path(final, "predictions.rds"))
  refit_rows <- readRDS(file.path(final, "selected_refit_rows.rds"))
  data <- gse_load_data(root, spec$stage)
  gse_assert(identical(as.character(results$method), gab_methods()) &&
    all(results$train_n == nrow(data$train$X)) &&
    all(results$test_n == nrow(data$test$X)) &&
    all(results$train_events == sum(data$train$y)) &&
    all(results$test_events == sum(data$test$y)) &&
    nrow(predictions) == 6L * nrow(data$test$X) &&
    all(predictions$probability >= 0 & predictions$probability <= 1) &&
    all(is.finite(predictions$probability)) &&
    all(vapply(split(predictions, predictions$method), function(p)
      identical(as.character(p$sample_id), data$test$sample_id) &&
        identical(as.integer(p$y), data$test$y), logical(1))),
    "Final scientific output dimension or test identity invalid")
  gse_assert(nrow(refit_rows) == 6L &&
    identical(as.character(refit_rows$method), gab_methods()) &&
    all(refit_rows$numerically_eligible %in% TRUE),
    "Selected refit row failed numerical eligibility")
  for (method in gab_methods()) {
    p <- predictions[predictions$method == method, , drop = FALSE]
    check <- gab_response_metrics(e, data$test$y, p$probability, 0.5)
    reported <- results[results$method == method, , drop = FALSE]
    for (metric in names(check))
      gse_assert(isTRUE(all.equal(reported[[metric]], check[[metric]],
        tolerance = 1e-10, check.attributes = FALSE)),
        paste("External-test metric reconstruction failed:", method, metric))
  }
  invisible(TRUE)
}
