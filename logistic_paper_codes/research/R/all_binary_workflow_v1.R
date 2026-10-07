# ALL BCR/ABL versus NEG real-data workflow built on the frozen V19 layer.
# The sglasso package, the V19 fitting code, and raw Bioconductor data are not
# modified by this workflow.

allb_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

allb_methods <- function() c(
  "Logistic SGLASSO",
  "Logistic SGLASSO (d=0 boundary)",
  "Logistic Group Elastic Net (adelie)",
  "Logistic Group Lasso (grpreg)",
  "Logistic Group MCP (grpreg)",
  "Logistic Group SCAD (grpreg)"
)

allb_source_files <- function() c(
  "R/genathum_binary_workflow_v1.R",
  "R/genathum_binary_io_v1.R",
  "R/genathum_binary_firth_v3.R",
  "R/all_binary_workflow_v1.R",
  "R/all_binary_io_v1.R",
  "data_processed/all_bcrabl_neg_gene_cytoband_v1.rds",
  "data_processed/all_bcrabl_neg_gene_cytoband_v1.rds.sha256",
  "ALL_BINARY_PROTOCOL_V1.md",
  "scripts/73_prepare_all_binary_v1.R",
  "scripts/74_run_all_binary_v1.R",
  "scripts/75_validate_all_binary_v1.R",
  "scripts/76_summarize_all_binary_v1.R",
  "scripts/77_preflight_all_binary_v1.R",
  "tests/test_all_binary_v1.R",
  "truba/run_all_binary_production_v1.slurm",
  "truba/build_all_binary_bundle_v1.R",
  "truba/README_TRUBA_ALL_BINARY_V1.md"
)

allb_load_environment <- function(root) {
  gab_load_environment(root)
}

allb_fit_configuration <- function(e, stage) {
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
  cfg
}

allb_configuration <- function(e, stage) {
  allb_assert(stage %in% c("smoke", "production"), "Invalid ALL binary stage.")
  list(
    schema_version = "all_binary_configuration_v1",
    stage = stage,
    repeats = if (stage == "smoke") 1L else 50L,
    outer_train_fraction = 0.70,
    inner_folds = if (stage == "smoke") 2L else 5L,
    seed_base = if (stage == "smoke") 940260921L else 950260921L,
    response_threshold = 0.5,
    selection_threshold = 1e-8,
    smoke_largest_groups = 8L,
    methods = allb_methods(),
    endpoint_id = "bcrabl_vs_neg",
    event_label = "BCR/ABL",
    reference_label = "NEG",
    fit_configuration = allb_fit_configuration(e, stage)
  )
}

allb_processed_path <- function(root) {
  file.path(root, "data_processed", "all_bcrabl_neg_gene_cytoband_v1.rds")
}

allb_load_data <- function(root, configuration) {
  path <- allb_processed_path(root)
  allb_assert(file.exists(path), "Missing processed ALL binary data.")
  x <- readRDS(path)
  allb_assert(identical(x$schema_version,
      "all_bcrabl_neg_gene_cytoband_v1") &&
    identical(dim(x$X), c(79L, 8545L)) &&
    length(x$y) == 79L && sum(x$y) == 37L && sum(x$y == 0L) == 42L &&
    identical(rownames(x$X), x$sample_id) &&
    identical(colnames(x$X), x$gene_annotation$entrez_id) &&
    length(x$group) == ncol(x$X) && length(x$group_name) == 762L &&
    all(is.finite(x$X)) && all(x$y %in% c(0L, 1L)) &&
    !anyDuplicated(x$sample_id) && !anyDuplicated(colnames(x$X)),
    "Processed ALL binary schema, values, dimensions, or identities changed.")
  X <- x$X
  group <- x$group
  gene_annotation <- x$gene_annotation
  group_name <- x$group_name
  if (identical(configuration$stage, "smoke")) {
    size <- tabulate(group, nbins = length(group_name))
    chosen <- order(-size, group_name)[seq_len(configuration$smoke_largest_groups)]
    keep <- group %in% chosen
    X <- X[, keep, drop = FALSE]
    gene_annotation <- gene_annotation[keep, , drop = FALSE]
    retained_names <- group_name[sort(chosen)]
    group <- match(gene_annotation$cytoband, retained_names)
    group_name <- retained_names
  }
  allb_assert(identical(colnames(X), gene_annotation$entrez_id) &&
    identical(group_name[group], gene_annotation$cytoband),
    "ALL gene/group mapping failed after stage-specific selection.")
  list(
    X = X, y = as.integer(x$y), response_label = x$response_label,
    sample_id = x$sample_id, group = as.integer(group),
    group_name = group_name, gene_annotation = gene_annotation,
    phenotype_audit = x$phenotype_audit, preprocessing = x$preprocessing,
    audit = x$audit, provenance = x$provenance,
    source_path = normalizePath(path, mustWork = TRUE),
    source_sha256 = digest::digest(file = path, algo = "sha256")
  )
}

allb_stratified_train <- function(y, fraction, seed) {
  set.seed(as.integer(seed))
  selected <- unlist(lapply(c(0L, 1L), function(value) {
    index <- which(y == value)
    number <- round(length(index) * fraction)
    number <- min(max(number, 1L), length(index) - 1L)
    sample(index, number, replace = FALSE)
  }), use.names = FALSE)
  sort(as.integer(selected))
}

allb_stratified_folds <- function(y, folds, seed) {
  assignment <- integer(length(y))
  set.seed(as.integer(seed))
  for (value in c(0L, 1L)) {
    index <- sample(which(y == value), replace = FALSE)
    assignment[index] <- rep(seq_len(folds), length.out = length(index))
  }
  assignment
}

allb_split_task <- function(data, replication, seed, configuration) {
  train_index <- allb_stratified_train(
    data$y, configuration$outer_train_fraction, seed
  )
  train <- seq_along(data$y) %in% train_index
  test <- !train
  inner_fold <- allb_stratified_folds(
    data$y[train], configuration$inner_folds, seed + 100000L
  )
  allb_assert(sum(train) == 55L && sum(test) == 24L &&
    identical(as.integer(table(data$y[train])), c(29L, 26L)) &&
    identical(as.integer(table(data$y[test])), c(13L, 11L)) &&
    length(unique(inner_fold)) == configuration$inner_folds &&
    all(vapply(seq_len(configuration$inner_folds), function(k) {
      all(c(0L, 1L) %in% data$y[train][inner_fold == k]) &&
        all(c(0L, 1L) %in% data$y[train][inner_fold != k])
    }, logical(1))), "A frozen ALL split or inner fold lacks a response class.")
  list(
    endpoint_id = configuration$endpoint_id,
    replication = as.integer(replication), seed = as.integer(seed),
    train = train, test = test, inner_fold = inner_fold,
    train_sample_id = data$sample_id[train],
    test_sample_id = data$sample_id[test]
  )
}

allb_tasks <- function(configuration) {
  tasks <- data.frame(
    task_id = seq_len(configuration$repeats),
    endpoint_id = configuration$endpoint_id,
    replication = seq_len(configuration$repeats),
    seed = as.integer(configuration$seed_base + seq_len(configuration$repeats)),
    stringsAsFactors = FALSE
  )
  tasks$key <- paste(tasks$endpoint_id, tasks$replication, sep = "::")
  tasks$shard_file <- sprintf("shard_%04d.rds", tasks$task_id)
  allb_assert(!anyDuplicated(tasks$key) && !anyDuplicated(tasks$seed),
    "Duplicate ALL task key or seed.")
  tasks
}

allb_run_task <- function(e, data, task, configuration) {
  split <- allb_split_task(
    data, task$replication[[1L]], task$seed[[1L]], configuration
  )
  X_train <- data$X[split$train, , drop = FALSE]
  y_train <- data$y[split$train]
  X_test <- data$X[split$test, , drop = FALSE]
  y_test <- data$y[split$test]
  inner <- gab_inner_select(
    e, X_train, y_train, data$group, split$inner_fold,
    configuration$fit_configuration
  )
  refit <- gab_refit_selected(
    e, X_train, y_train, X_test, data$group,
    configuration$fit_configuration, inner
  )
  results <- predictions <- supports <- list()
  for (method in allb_methods()) {
    model <- refit$models[[method]]
    metric <- gab_response_metrics(
      e, y_test, model$probability, configuration$response_threshold
    )
    beta <- model$coefficient[-1L]
    active_predictor <- abs(beta) > configuration$selection_threshold
    group_levels <- seq_along(data$group_name)
    active_group <- vapply(group_levels, function(g) {
      any(active_predictor[data$group == g])
    }, logical(1))
    selected <- inner$selections[[method]]
    results[[method]] <- cbind(data.frame(
      task_id = task$task_id[[1L]], endpoint_id = configuration$endpoint_id,
      replication = task$replication[[1L]], seed = task$seed[[1L]],
      method = method, selected_alpha = selected$alpha[[1L]],
      selected_d = selected$d[[1L]], selected_gamma = selected$gamma[[1L]],
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
    predictions[[method]] <- data.frame(
      task_id = task$task_id[[1L]], endpoint_id = configuration$endpoint_id,
      replication = task$replication[[1L]], method = method,
      sample_id = data$sample_id[split$test], y = y_test,
      probability = model$probability, stringsAsFactors = FALSE
    )
    supports[[method]] <- data.frame(
      task_id = task$task_id[[1L]], endpoint_id = configuration$endpoint_id,
      replication = task$replication[[1L]], method = method,
      predictor_group = group_levels,
      predictor_group_name = data$group_name,
      selected = active_group, stringsAsFactors = FALSE
    )
  }
  list(
    results = do.call(rbind, results),
    predictions = do.call(rbind, predictions),
    selected_groups = do.call(rbind, supports),
    inner_candidates = inner$candidate_summary,
    inner_folds = inner$fold_diagnostics,
    split = split[c("endpoint_id", "replication", "seed",
      "train_sample_id", "test_sample_id")],
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

allb_validate_payload <- function(payload, task, configuration) {
  checks <- c(
    six_methods = is.data.frame(payload$results) &&
      identical(as.character(payload$results$method), allb_methods()),
    finite_primary_metrics = all(is.finite(payload$results$log_loss)) &&
      all(is.finite(payload$results$brier)) && all(is.finite(payload$results$auc)),
    probability_rows = nrow(payload$predictions) ==
      nrow(payload$results) * payload$results$test_n[[1L]],
    test_only_probabilities = all(payload$predictions$y %in% c(0L, 1L)) &&
      all(payload$predictions$probability >= 0 &
        payload$predictions$probability <= 1),
    independent_split = !length(intersect(
      payload$split$train_sample_id, payload$split$test_sample_id)),
    split_sizes = all(payload$results$train_n == 55L) &&
      all(payload$results$test_n == 24L) &&
      all(payload$results$train_events == 26L) &&
      all(payload$results$test_events == 11L),
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
    one_inner_selection_per_method = all(vapply(allb_methods(), function(m) {
      sum(payload$inner_candidates$method == m &
        payload$inner_candidates$selected) == 1L
    }, logical(1))),
    task_identity = all(payload$results$task_id == task$task_id[[1L]]) &&
      all(payload$results$endpoint_id == task$endpoint_id[[1L]]) &&
      all(payload$results$replication == task$replication[[1L]])
  )
  data.frame(check = names(checks), passed = unname(checks),
    stringsAsFactors = FALSE)
}

allb_firth_feasibility <- function(e, root, cores = 2L) {
  configuration <- allb_configuration(e, "production")
  data <- allb_load_data(root, configuration)
  task <- allb_tasks(configuration)[1L, , drop = FALSE]
  split <- allb_split_task(data, task$replication, task$seed, configuration)
  X <- data$X[split$train, , drop = FALSE]
  y <- data$y[split$train]
  units <- c(seq_len(configuration$inner_folds), 0L)
  one <- function(fold) {
    keep <- if (fold == 0L) rep(TRUE, length(y)) else split$inner_fold != fold
    started <- proc.time()[["elapsed"]]
    fit <- tryCatch(e$gab_estimate_training_firth_target(
      e, X[keep, , drop = FALSE], y[keep], data$group,
      configuration$fit_configuration
    ), error = function(error) error)
    success <- !inherits(fit, "error") && isTRUE(fit$success)
    data.frame(
      fold = fold, training_n = sum(keep), predictors = ncol(X),
      groups = length(data$group_name), success = success,
      failed_groups = if (success) fit$failed_groups else NA_integer_,
      maximum_score = if (success)
        max(fit$diagnostics$maximum_adjusted_score) else NA_real_,
      maximum_step = if (success)
        max(fit$diagnostics$maximum_fisher_step) else NA_real_,
      replay_error = if (success) fit$replay_error else NA_real_,
      elapsed_seconds = proc.time()[["elapsed"]] - started,
      message = if (success) "" else conditionMessage(fit),
      stringsAsFactors = FALSE
    )
  }
  rows <- parallel::mclapply(units, one,
    mc.cores = min(as.integer(cores), length(units)),
    mc.set.seed = FALSE, mc.preschedule = FALSE)
  rows <- do.call(rbind, rows)
  accepted <- nrow(rows) == 6L && all(rows$success) &&
    all(rows$failed_groups == 0L) &&
    all(rows$maximum_score <= configuration$fit_configuration$target_tolerance) &&
    all(rows$maximum_step <= configuration$fit_configuration$target_tolerance)
  list(
    schema_version = "all_binary_firth_feasibility_v1",
    accepted = accepted, task = task, split = split,
    results = rows, processed_data_sha256 = data$source_sha256,
    target_configuration = configuration$fit_configuration,
    scope = "six training-only Firth targets; no penalized path or test metric fitted"
  )
}
