# Binary GenAtHum real-data workflow built on the frozen V19 fitting layer.
# This file does not modify package code or any frozen simulation source.

gab_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

gab_methods <- function() c(
  "Logistic SGLASSO",
  "Logistic SGLASSO (d=0 boundary)",
  "Logistic Group Elastic Net (adelie)",
  "Logistic Group Lasso (grpreg)",
  "Logistic Group MCP (grpreg)",
  "Logistic Group SCAD (grpreg)"
)

gab_source_files <- function() c(
  "R/genathum_binary_workflow_v1.R",
  "R/genathum_binary_io_v1.R",
  "R/genathum_binary_firth_v3.R",
  "config/genathum_binary_outcomes_v1.csv",
  "GENATHUM_BINARY_PROTOCOL_V1.md",
  "scripts/70_run_genathum_binary_v1.R",
  "scripts/71_validate_genathum_binary_v1.R",
  "scripts/72_summarize_genathum_binary_v1.R",
  "tests/test_genathum_binary_v1.R",
  "tests/test_genathum_firth_v3.R",
  "truba/run_genathum_binary_production_v1.slurm"
)

gab_load_environment <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  scope <- new.env(parent = globalenv())
  source(file.path(root, "R/logistic_sglasso_workflow_v19.R"), local = scope)
  e <- scope$lsg_load_v19(root)
  e$compile_lsg_core(quiet = TRUE)
  e$lsg_compile_hybrid_solver_v11(root, rebuild = FALSE, quiet = TRUE)
  source(file.path(root, "R/genathum_binary_firth_v3.R"), local = e)
  e
}

gab_read_outcomes <- function(root, stage) {
  x <- utils::read.csv(
    file.path(root, "config/genathum_binary_outcomes_v1.csv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  required <- c("outcome_id", "role", "cutpoint", "operator", "description")
  gab_assert(identical(names(x), required), "Invalid GenAtHum outcome schema.")
  gab_assert(nrow(x) == 2L && !anyDuplicated(x$outcome_id) &&
    identical(x$role, c("primary", "sensitivity")) &&
    all(x$operator == "ge") && all(is.finite(x$cutpoint)),
    "The two frozen GenAtHum endpoints changed.")
  if (identical(stage, "smoke")) x <- x[x$role == "primary", , drop = FALSE]
  rownames(x) <- NULL
  x
}

gab_configuration <- function(e, root, stage) {
  gab_assert(stage %in% c("smoke", "production"), "Invalid GenAtHum stage.")
  cfg <- e$lsg_configuration_v7(if (stage == "smoke") "smoke" else "production")
  if (stage == "smoke") {
    cfg$alpha_grid <- c(0, 1)
    cfg$d_grid <- c(0, 1)
    cfg$benchmark_alpha_grid <- c(0, 1)
  }
  cfg$target_solver <- "training_solver_coordinates_firth_bfgs_v3"
  cfg$target_max_iterations <- 5000L
  cfg$target_tolerance <- 1e-6
  cfg$target_bfgs_reltol <- 1e-15
  cfg$target_bfgs_passes <- 4L
  cfg$target_polish_steps <- 5L
  cfg$target_hessian_step <- 1e-4
  cfg$target_replay_tolerance <- 1e-9
  list(
    schema_version = "genathum_binary_configuration_v1",
    stage = stage,
    repeats = if (stage == "smoke") 1L else 50L,
    outer_train_fraction = 0.70,
    inner_folds = if (stage == "smoke") 2L else 5L,
    seed_base = if (stage == "smoke") 920260921L else 930260921L,
    response_threshold = 0.5,
    selection_threshold = 1e-8,
    predictor_groups_smoke = 79L,
    methods = gab_methods(),
    outcomes = gab_read_outcomes(root, stage),
    fit_configuration = cfg
  )
}

gab_load_data <- function(project_root, configuration) {
  path <- file.path(project_root, "data/GenAtHum.rda")
  gab_assert(file.exists(path), "Missing package data/GenAtHum.rda.")
  box <- new.env(parent = emptyenv())
  loaded <- load(path, envir = box)
  gab_assert(identical(loaded, "GenAtHum"), "Unexpected GenAtHum RDA contents.")
  x <- box$GenAtHum
  required <- c("X", "y", "group", "gene_code", "gene_name", "groups_name")
  gab_assert(is.list(x) && all(required %in% names(x)), "Invalid GenAtHum object.")
  X <- as.matrix(x$X)
  y_continuous <- as.numeric(x$y)
  sample_id <- rownames(X)
  gab_assert(is.numeric(X) && identical(dim(X), c(158L, 2045L)) &&
    length(y_continuous) == nrow(X) && all(is.finite(X)) &&
    all(is.finite(y_continuous)) && length(x$group) == ncol(X) &&
    length(x$groups_name) == 79L && !anyDuplicated(x$groups_name) &&
    length(sample_id) == nrow(X) && !anyDuplicated(sample_id) &&
    all(grepl("^GSM[0-9]+$", sample_id)),
    "GenAtHum dimensions, values, groups, or GSM identifiers changed.")
  gab_assert(nrow(X) == 2L * length(x$groups_name),
    "The expected two samples per tissue/cell-type block are absent.")
  sample_block <- rep(as.character(x$groups_name), each = 2L)
  predictor_group <- as.integer(factor(x$group, levels = unique(x$group)))
  predictor_group_name <- as.character(x$groups_name)[predictor_group]
  if (configuration$stage == "smoke") {
    keep <- predictor_group <= configuration$predictor_groups_smoke
    X <- X[, keep, drop = FALSE]
    predictor_group <- predictor_group[keep]
    predictor_group_name <- predictor_group_name[keep]
  }
  list(
    X = X,
    y_continuous = y_continuous,
    sample_id = sample_id,
    sample_block = sample_block,
    group = predictor_group,
    predictor_group_name = predictor_group_name,
    groups_name = as.character(x$groups_name),
    source_path = normalizePath(path, mustWork = TRUE),
    source_sha256 = digest::digest(file = path, algo = "sha256")
  )
}

gab_make_labels <- function(y_continuous, outcome) {
  gab_assert(identical(outcome$operator[[1L]], "ge"), "Unsupported endpoint operator.")
  as.integer(y_continuous >= outcome$cutpoint[[1L]])
}

gab_stratified_block_partition <- function(block, y, train_fraction, folds,
                                            seed) {
  block <- as.character(block)
  ids <- unique(block)
  event_count <- vapply(ids, function(id) sum(y[block == id]), integer(1))
  size <- vapply(ids, function(id) sum(block == id), integer(1))
  gab_assert(all(size == 2L) && all(event_count %in% 0:2),
    "Each tissue block must contain exactly two binary outcomes.")
  strata <- split(ids, event_count)
  set.seed(as.integer(seed))
  shuffled <- unlist(lapply(strata, sample), use.names = FALSE)
  stratum <- event_count[match(shuffled, ids)]
  if (!is.null(folds)) {
    assignment <- integer(length(shuffled))
    for (s in sort(unique(stratum))) {
      at <- which(stratum == s)
      assignment[at] <- rep(seq_len(folds), length.out = length(at))
    }
    names(assignment) <- shuffled
    return(assignment)
  }
  selected <- character(0)
  for (s in sort(unique(stratum))) {
    candidates <- shuffled[stratum == s]
    number <- round(length(candidates) * train_fraction)
    number <- min(max(number, 1L), length(candidates) - 1L)
    selected <- c(selected, candidates[seq_len(number)])
  }
  selected
}

gab_split_task <- function(data, outcome, replication, seed, configuration) {
  y <- gab_make_labels(data$y_continuous, outcome)
  train_blocks <- gab_stratified_block_partition(
    data$sample_block, y, configuration$outer_train_fraction, NULL, seed
  )
  train <- data$sample_block %in% train_blocks
  test <- !train
  fold_by_block <- gab_stratified_block_partition(
    data$sample_block[train], y[train], NA_real_, configuration$inner_folds,
    seed + 100000L
  )
  inner_fold <- unname(fold_by_block[data$sample_block[train]])
  gab_assert(all(c(0L, 1L) %in% y[train]) && all(c(0L, 1L) %in% y[test]) &&
    length(unique(inner_fold)) == configuration$inner_folds &&
    all(vapply(seq_len(configuration$inner_folds), function(k) {
      all(c(0L, 1L) %in% y[train][inner_fold == k]) &&
        all(c(0L, 1L) %in% y[train][inner_fold != k])
    }, logical(1))), "A frozen split lacks a response class.")
  list(
    outcome_id = outcome$outcome_id[[1L]], replication = as.integer(replication),
    seed = as.integer(seed), y = y, train = train, test = test,
    inner_fold = inner_fold, train_blocks = sort(train_blocks),
    test_blocks = sort(setdiff(unique(data$sample_block), train_blocks))
  )
}

gab_tasks <- function(configuration) {
  rows <- lapply(seq_len(nrow(configuration$outcomes)), function(i) {
    data.frame(
      outcome_id = configuration$outcomes$outcome_id[i],
      outcome_index = i,
      replication = seq_len(configuration$repeats),
      stringsAsFactors = FALSE
    )
  })
  tasks <- do.call(rbind, rows)
  tasks$task_id <- seq_len(nrow(tasks))
  tasks$seed <- as.integer(configuration$seed_base +
    tasks$outcome_index * 10000L + tasks$replication)
  tasks$key <- paste(tasks$outcome_id, tasks$replication, sep = "::")
  tasks$shard_file <- sprintf("shard_%04d.rds", tasks$task_id)
  tasks <- tasks[c("task_id", "outcome_index", "outcome_id", "replication",
    "seed", "key", "shard_file")]
  rownames(tasks) <- NULL
  gab_assert(!anyDuplicated(tasks$key) && !anyDuplicated(tasks$seed),
    "Duplicate GenAtHum task key or seed.")
  tasks
}

gab_tuning_key <- function(frame, method) {
  ratio <- frame$lambda_relative_to_reference
  ratio_text <- ifelse(is.infinite(ratio), "Inf", sprintf("%.12g", ratio))
  if (grepl("SGLASSO", method, fixed = TRUE)) {
    paste(method, frame$alpha_index, frame$d_index, frame$point_type,
      frame$lambda_index, ratio_text, sep = "|")
  } else {
    paste(method, frame$alpha_index, frame$fit_index, frame$point_type,
      frame$lambda_index, ratio_text, sep = "|")
  }
}

gab_candidate_frame <- function(bundle, method) {
  if (identical(method, "Logistic SGLASSO")) {
    x <- bundle$sglasso$tuning
  } else if (identical(method, "Logistic SGLASSO (d=0 boundary)")) {
    x <- bundle$sglasso$tuning[abs(bundle$sglasso$tuning$d) <= 1e-12, , drop = FALSE]
  } else {
    x <- bundle$external$tuning[bundle$external$tuning$method_path == method, , drop = FALSE]
  }
  x$method <- method
  x$tuning_key <- gab_tuning_key(x, method)
  x
}

gab_fit_bundle <- function(e, X_train, y_train, X_validation, y_validation,
                           group, configuration) {
  tuning_data <- list(
    X_train = X_train, y_train = y_train,
    X_validation = X_validation, y_validation = y_validation, group = group
  )
  e$lsg_tail_validate_data_v6(tuning_data)
  firth <- e$gab_estimate_training_firth_target(
    e, X_train, y_train, group, configuration)
  gab_assert(isTRUE(firth$success) && firth$failed_groups == 0L &&
    all(is.finite(firth$target_original)),
    "Training-only groupwise Firth target failed; no fallback is allowed.")
  fitters <- e$lsg_finite_fitters_v7(configuration)
  list(
    sglasso = e$lsg_fit_sglasso_joint_v7(
      tuning_data, firth$target_original, configuration, fitters
    ),
    external = e$lsg_fit_external_joint_v7(tuning_data, configuration, fitters),
    firth = firth
  )
}

gab_inner_select <- function(e, X, y, group, inner_fold, configuration) {
  records <- list()
  diagnostics <- list()
  for (fold in sort(unique(inner_fold))) {
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
        bundle$firth$diagnostics$maximum_adjusted_score),
      firth_maximum_fisher_step = max(
        bundle$firth$diagnostics$maximum_fisher_step),
      firth_replay_error = bundle$firth$replay_error,
      stringsAsFactors = FALSE
    )
    rm(bundle); invisible(gc(FALSE))
  }
  tuning <- do.call(e$lsg_bind_tuning_rows_v7, records)
  selections <- list()
  summaries <- list()
  for (method in gab_methods()) {
    z <- tuning[tuning$method == method, , drop = FALSE]
    keys <- unique(z$tuning_key)
    candidate <- lapply(keys, function(key) {
      q <- z[z$tuning_key == key, , drop = FALSE]
      represented <- nrow(q) == length(unique(inner_fold)) &&
        !anyDuplicated(q$fold)
      eligible <- represented && all(q$numerically_eligible %in% TRUE) &&
        all(is.finite(q$validation_log_loss))
      template <- q[1L, , drop = FALSE]
      data.frame(
        method = method, tuning_key = key,
        alpha = template$alpha[[1L]],
        d = if ("d" %in% names(template)) template$d[[1L]] else NA_real_,
        gamma = template$gamma[[1L]], point_type = template$point_type[[1L]],
        lambda_index = template$lambda_index[[1L]],
        lambda_relative_to_reference =
          template$lambda_relative_to_reference[[1L]],
        folds_represented = length(unique(q$fold)), eligible = eligible,
        inner_log_loss = if (eligible) weighted.mean(
          q$validation_log_loss, q$fold_size) else NA_real_,
        stringsAsFactors = FALSE
      )
    })
    candidate <- do.call(rbind, candidate)
    valid <- candidate[candidate$eligible, , drop = FALSE]
    gab_assert(nrow(valid) > 0L,
      paste("No common eligible inner-CV candidate for", method))
    valid <- valid[order(
      valid$inner_log_loss, valid$alpha,
      ifelse(is.na(valid$d), 0, valid$d),
      valid$point_type != "penalty_limit",
      -valid$lambda_relative_to_reference, valid$lambda_index
    ), , drop = FALSE]
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

gab_refit_selected <- function(e, X_train, y_train, X_test, group,
                               configuration, inner_selection) {
  bundle <- gab_fit_bundle(
    e, X_train, y_train, X_train, y_train, group, configuration
  )
  data_interface <- list(
    X_train = X_train, y_train = y_train,
    X_validation = X_train, y_validation = y_train,
    X_test = X_test, group = group
  )
  models <- list()
  selected_rows <- list()
  for (method in gab_methods()) {
    target_key <- inner_selection$selections[[method]]$tuning_key[[1L]]
    candidates <- gab_candidate_frame(bundle, method)
    row <- candidates[candidates$tuning_key == target_key, , drop = FALSE]
    gab_assert(nrow(row) == 1L && row$numerically_eligible[[1L]] %in% TRUE,
      paste("Selected candidate could not be refit for", method))
    if (grepl("SGLASSO", method, fixed = TRUE)) {
      model <- e$lsg_selected_sglasso_model_v7(
        bundle$sglasso, list(row = row), data_interface
      )
    } else if (identical(method, "Logistic Group Elastic Net (adelie)")) {
      bundle$external$adelie$selection$row <- row
      solution <- e$extract_adelie_group_en_solution_v2(
        bundle$external$adelie, X_test
      )
      model <- list(coefficient = solution$coefficient,
        probability = solution$probability)
    } else {
      spec <- bundle$external$grpreg$penalty_specification
      penalty <- spec$penalty[spec$method_path == method]
      bundle$external$grpreg$selections[[penalty]]$row <- row
      solution <- e$extract_grpreg_group_solution_v2(
        bundle$external$grpreg, penalty, X_test
      )
      model <- list(coefficient = solution$coefficient,
        probability = solution$probability)
    }
    gab_assert(length(model$coefficient) == ncol(X_train) + 1L &&
      all(is.finite(model$coefficient)) && length(model$probability) == nrow(X_test) &&
      all(is.finite(model$probability)) &&
      all(model$probability >= 0 & model$probability <= 1),
      paste("Invalid refitted model for", method))
    models[[method]] <- model
    row$method <- method
    selected_rows[[method]] <- row
  }
  list(models = models,
    selected_rows = do.call(e$lsg_bind_tuning_rows_v7, selected_rows),
    firth = bundle$firth)
}

gab_auc <- function(y, probability) {
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  (sum(rank(probability, ties.method = "average")[y == 1L]) -
    n1 * (n1 + 1) / 2) / (n1 * n0)
}

gab_calibration <- function(y, probability) {
  link <- stats::qlogis(pmin(pmax(probability, 1e-8), 1 - 1e-8))
  intercept <- slope <- NA_real_
  if (length(unique(y)) == 2L) {
    intercept <- suppressWarnings(tryCatch(
      unname(stats::coef(stats::glm(y ~ 1 + offset(link), family = stats::binomial()))[1L]),
      error = function(e) NA_real_))
    if (length(unique(link)) > 1L) {
      slope <- suppressWarnings(tryCatch(
        unname(stats::coef(stats::glm(y ~ link, family = stats::binomial()))[2L]),
        error = function(e) NA_real_))
    }
  }
  c(intercept = intercept, slope = slope)
}

gab_response_metrics <- function(e, y, probability, threshold = 0.5) {
  predicted <- probability >= threshold
  scores <- e$lsg_binary_scores_v7(y == 1L, predicted)
  calibration <- gab_calibration(y, probability)
  data.frame(
    log_loss = e$binary_log_loss(y, probability),
    brier = mean((y - probability)^2), auc = gab_auc(y, probability),
    average_precision = e$lsg_average_precision_v7(y, probability),
    classification_error = mean(predicted != (y == 1L)),
    sensitivity = scores$tpr, specificity = scores$tnr,
    balanced_accuracy = scores$balanced_accuracy, mcc = scores$mcc,
    mcc_degenerate = scores$mcc_degenerate,
    predicted_positive_rate = mean(predicted),
    calibration_intercept = calibration[["intercept"]],
    calibration_slope = calibration[["slope"]],
    auc_defined = is.finite(gab_auc(y, probability)),
    average_precision_defined = sum(y == 1L) > 0L,
    balanced_accuracy_defined = scores$balanced_accuracy_defined,
    calibration_intercept_defined = is.finite(calibration[["intercept"]]),
    calibration_slope_defined = is.finite(calibration[["slope"]]),
    stringsAsFactors = FALSE
  )
}

gab_run_task <- function(e, data, outcome, task, configuration) {
  split <- gab_split_task(
    data, outcome, task$replication[[1L]], task$seed[[1L]], configuration
  )
  X_train <- data$X[split$train, , drop = FALSE]
  y_train <- split$y[split$train]
  X_test <- data$X[split$test, , drop = FALSE]
  y_test <- split$y[split$test]
  inner <- gab_inner_select(
    e, X_train, y_train, data$group, split$inner_fold,
    configuration$fit_configuration
  )
  refit <- gab_refit_selected(
    e, X_train, y_train, X_test, data$group,
    configuration$fit_configuration, inner
  )
  results <- probabilities <- supports <- list()
  for (method in gab_methods()) {
    model <- refit$models[[method]]
    metric <- gab_response_metrics(
      e, y_test, model$probability, configuration$response_threshold
    )
    beta <- model$coefficient[-1L]
    active_predictor <- abs(beta) > configuration$selection_threshold
    group_levels <- unique(data$group)
    active_group <- vapply(group_levels, function(g) {
      any(active_predictor[data$group == g])
    }, logical(1))
    selected <- inner$selections[[method]]
    results[[method]] <- cbind(data.frame(
      task_id = task$task_id[[1L]], outcome_id = task$outcome_id[[1L]],
      outcome_role = outcome$role[[1L]], replication = task$replication[[1L]],
      seed = task$seed[[1L]], method = method,
      selected_alpha = selected$alpha[[1L]], selected_d = selected$d[[1L]],
      selected_gamma = selected$gamma[[1L]],
      selected_point_type = selected$point_type[[1L]],
      selected_lambda_index = selected$lambda_index[[1L]],
      selected_lambda_relative_to_reference =
        selected$lambda_relative_to_reference[[1L]],
      inner_cv_log_loss = selected$inner_log_loss[[1L]],
      train_n = length(y_train), test_n = length(y_test),
      train_events = sum(y_train), test_events = sum(y_test),
      selected_predictors = sum(active_predictor),
      selected_groups = sum(active_group), stringsAsFactors = FALSE
    ), metric)
    probabilities[[method]] <- data.frame(
      task_id = task$task_id[[1L]], outcome_id = task$outcome_id[[1L]],
      replication = task$replication[[1L]], method = method,
      sample_id = data$sample_id[split$test],
      sample_block = data$sample_block[split$test], y = y_test,
      probability = model$probability, stringsAsFactors = FALSE
    )
    supports[[method]] <- data.frame(
      task_id = task$task_id[[1L]], outcome_id = task$outcome_id[[1L]],
      replication = task$replication[[1L]], method = method,
      predictor_group = group_levels,
      predictor_group_name = data$groups_name[group_levels],
      selected = active_group, stringsAsFactors = FALSE
    )
  }
  list(
    results = do.call(rbind, results),
    predictions = do.call(rbind, probabilities),
    selected_groups = do.call(rbind, supports),
    inner_candidates = inner$candidate_summary,
    inner_folds = inner$fold_diagnostics,
    split = split[c("outcome_id", "replication", "seed", "train_blocks", "test_blocks")],
    selected_rows = refit$selected_rows,
    firth_refit_diagnostics = refit$firth$diagnostics,
    firth_refit_replay_error = refit$firth$replay_error,
    data_fingerprint = e$lsg_data_fingerprint_v7(list(
      train_sample_id = data$sample_id[split$train],
      test_sample_id = data$sample_id[split$test],
      y_train = y_train, y_test = y_test
    ))
  )
}

gab_validate_payload <- function(payload, task, configuration) {
  checks <- c(
    six_methods = is.data.frame(payload$results) &&
      identical(as.character(payload$results$method), gab_methods()),
    finite_primary_metrics = all(is.finite(payload$results$log_loss)) &&
      all(is.finite(payload$results$brier)),
    probability_rows = nrow(payload$predictions) ==
      nrow(payload$results) * payload$results$test_n[[1L]],
    test_only_probabilities = all(payload$predictions$y %in% c(0L, 1L)) &&
      all(payload$predictions$probability >= 0 &
        payload$predictions$probability <= 1),
    grouped_split = !length(intersect(
      payload$split$train_blocks, payload$split$test_blocks)),
    all_inner_folds_valid = nrow(payload$inner_folds) ==
      configuration$inner_folds &&
      all(payload$inner_folds$train_events > 0) &&
      all(payload$inner_folds$validation_events > 0),
    all_training_firth_targets_valid =
      all(payload$inner_folds$firth_failed_groups == 0L) &&
      all(payload$inner_folds$firth_maximum_adjusted_score <=
        configuration$fit_configuration$target_tolerance) &&
      all(payload$inner_folds$firth_maximum_fisher_step <=
        configuration$fit_configuration$target_tolerance) &&
      is.data.frame(payload$firth_refit_diagnostics) &&
      all(payload$firth_refit_diagnostics$success),
    one_inner_selection_per_method = all(vapply(gab_methods(), function(m) {
      sum(payload$inner_candidates$method == m &
        payload$inner_candidates$selected) == 1L
    }, logical(1))),
    task_identity = all(payload$results$task_id == task$task_id[[1L]]) &&
      all(payload$results$outcome_id == task$outcome_id[[1L]]) &&
      all(payload$results$replication == task$replication[[1L]])
  )
  data.frame(check = names(checks), passed = unname(checks),
    stringsAsFactors = FALSE)
}
