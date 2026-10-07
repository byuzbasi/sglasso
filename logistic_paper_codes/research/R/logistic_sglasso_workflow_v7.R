# Standalone V7 orchestration. Frozen package/core and earlier outputs are read-only.

lsg_required_packages_v7 <- function() {
  c("Rcpp", "RcppArmadillo", "digest", "adelie", "grpreg", "logistf", "mltools")
}

lsg_source_files_v7 <- function() {
  c(file.path("R", c(
    "logistic_sglasso.R", "pilot_utils.R", "study_utils.R",
    "groupwise_logistic_targets.R", "correlated_logistic_design.R",
    "correlated_target_study_utils.R", "sglasso_rhob_design.R",
    "logistic_penalized_benchmarks_v1.R", "sglasso_rhob_benchmark_study_utils_v1.R",
    "logistic_group_benchmarks_v2.R", "sglasso_rhob_group_benchmark_study_utils_v2.R",
    "logistic_sglasso_two_design_v1.R", "logistic_sglasso_two_design_study_utils_v1.R",
    "logistic_sglasso_lambda_extension_v5.R", "logistic_sglasso_lambda_tail_v6.R",
    "logistic_sglasso_lambda_tail_io_v6.R", "logistic_sglasso_design_metrics_v7.R",
    "logistic_sglasso_prediction_selection_v7.R", "logistic_sglasso_workflow_v7.R")),
    "src/logistic_sglasso_core.cpp",
    "config/logistic_sglasso_prediction_selection_v7.csv",
    "LOGISTIC_SGLASSO_PREDICTION_SELECTION_PROTOCOL_V7.md",
    "scripts/42_run_logistic_sglasso_prediction_selection_v7.R")
}

lsg_source_v7 <- function(root, envir = environment(lsg_source_v7)) {
  options(sglasso.logistic_prework.root = normalizePath(root, mustWork = TRUE))
  files <- lsg_source_files_v7()
  files <- files[grepl("^R/", files) & basename(files) != "logistic_sglasso_workflow_v7.R"]
  for (file in files) source(file.path(root, file), local = envir)
  invisible(envir)
}

lsg_assert_v7 <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

lsg_equal_v7 <- function(x, y, tolerance = 1e-10) {
  isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE))
}

# Typed canonical values avoid R serialization writer-version headers in
# scientific signatures. All real bytes have a fixed endian order.
lsg_canonical_v7 <- function(value) {
  characters <- function(x) paste(vapply(x, function(s) {
    if (is.na(s)) return("NA;")
    s <- enc2utf8(s)
    paste0(nchar(s, type = "bytes"), ":", s, ";")
  }, character(1)), collapse = "")
  if (is.null(value)) return("null;")
  if (is.factor(value) || is.environment(value) || is.function(value)) {
    stop("Unsupported object in a V7 scientific signature.", call. = FALSE)
  }
  attributes_text <- if (is.data.frame(value)) {
    paste0("rows:", nrow(value), ";cols:", characters(names(value)))
  } else paste0("names:", characters(names(value)), ";dim:",
                 paste(dim(value), collapse = ","), ";")
  body <- if (is.list(value)) {
    paste(vapply(value, lsg_canonical_v7, character(1)), collapse = "")
  } else if (is.character(value)) characters(value) else if (is.logical(value)) {
    paste(ifelse(is.na(value), "N", ifelse(value, "T", "F")), collapse = "")
  } else if (is.integer(value)) {
    paste(ifelse(is.na(value), "N", as.character(value)), collapse = ",")
  } else if (is.double(value)) {
    paste(format(writeBin(as.numeric(value), raw(), size = 8L, endian = "little")),
          collapse = "")
  } else stop("Unsupported scalar type in a V7 signature.", call. = FALSE)
  paste0(typeof(value), "[", length(value), "]{", attributes_text,
         nchar(body, type = "bytes"), ":", body, "}")
}

lsg_hash_v7 <- function(value) {
  digest::digest(lsg_canonical_v7(value), algo = "sha256", serialize = FALSE)
}

lsg_file_hash_v7 <- function(path) digest::digest(file = path, algo = "sha256")

# Data fingerprints hash numeric bytes directly rather than expanding a large
# design matrix into a hexadecimal string or hashing an R serialization header.
lsg_data_fingerprint_v7 <- function(value) {
  if (is.factor(value) || is.environment(value) || is.function(value))
    stop("Unsupported data fingerprint type.", call. = FALSE)
  if (is.list(value)) return(lsg_hash_v7(list(type = "list", names = names(value),
    children = unname(vapply(value, lsg_data_fingerprint_v7, character(1))))))
  if (is.numeric(value)) {
    body <- if (is.integer(value)) writeBin(as.integer(value), raw(), size = 4L,
      endian = "little") else writeBin(as.numeric(value), raw(), size = 8L, endian = "little")
    return(lsg_hash_v7(list(type = typeof(value), length = length(value),
      dimensions = dim(value), dimnames = dimnames(value), names = names(value),
      data_sha256 = digest::digest(body, serialize = FALSE, algo = "sha256"))))
  }
  lsg_hash_v7(value)
}

lsg_path_v7 <- function(root, relative) {
  lsg_assert_v7(length(root) == 1L && dir.exists(root), "Missing V7 root directory.")
  lsg_assert_v7(is.character(relative) && length(relative) == 1L &&
    nzchar(relative) && !grepl("(^/|\\\\|(^|/)\\.\\.?(/|$))", relative),
    "Unsafe relative output path.")
  root <- normalizePath(root, mustWork = TRUE)
  path <- root
  for (part in strsplit(relative, "/", fixed = TRUE)[[1L]]) {
    path <- file.path(path, part)
    link <- Sys.readlink(path)
    lsg_assert_v7(is.na(link) || !nzchar(link), paste("Symlink path refused:", path))
  }
  path
}

lsg_inventory_v7 <- function(root, files) {
  files <- sort(unique(files), method = "radix")
  paths <- vapply(files, function(f) lsg_path_v7(root, f), character(1))
  lsg_assert_v7(all(file.exists(paths)) && !any(dir.exists(paths)),
                "A required V7 manifest file is absent.")
  data.frame(file = files, bytes = as.numeric(file.info(paths)$size),
    sha256 = unname(vapply(paths, lsg_file_hash_v7, character(1))),
    stringsAsFactors = FALSE)
}

lsg_atomic_v7 <- function(value, directory, relative, kind = "rds",
                          identical_ok = FALSE) {
  path <- lsg_path_v7(directory, relative)
  lsg_tail_atomic_v6(value, path, kind, allow_identical = identical_ok)
}

lsg_runtime_v7 <- function() {
  packages <- lsg_required_packages_v7()
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  lsg_assert_v7(!length(missing), paste("Missing installed packages (no automatic install):",
                                      paste(missing, collapse = ", ")))
  list(r_version = R.version.string, platform = R.version$platform,
       package_versions = stats::setNames(vapply(packages, function(p) {
         as.character(utils::packageVersion(p))
       }, character(1)), packages),
       rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"))
}

lsg_default_version_v7 <- function(stage) {
  switch(stage, smoke = "logistic_sglasso_prediction_selection_smoke_v7",
    pilot = "logistic_sglasso_prediction_selection_pilot_r20_v7",
    stop("Only smoke and pilot stages are implemented.", call. = FALSE))
}

lsg_configuration_v7 <- function(stage) {
  smoke <- identical(stage, "smoke")
  lsg_assert_v7(stage %in% c("smoke", "pilot"), "Invalid V7 stage.")
  alpha <- if (smoke) c(0, 0.5, 1) else seq(0, 1, by = 0.1)
  ratios <- c(512, 256, 128, 64, 32, 16,
              exp(seq(log(8), log(0.05), length.out = 30L)))
  list(stage = stage, solver = "abgd", methods = c(
    "Logistic SGLASSO", "Logistic SGLASSO (d=0 boundary)",
    "Logistic Group Elastic Net (adelie)", "Logistic Group Lasso (grpreg)",
    "Logistic Group MCP (grpreg)", "Logistic Group SCAD (grpreg)"),
    replications_per_scenario = if (smoke) 1L else 20L,
    seed_base = if (smoke) 720260831L else 730260831L, seed_stride = 10000L,
    alpha_grid = alpha, d_grid = alpha, nlambda = length(ratios),
    lambda_relative_grid = ratios, lambda_base_nlambda = 30L,
    lambda_extension_multipliers = c(512, 256, 128, 64, 32, 16),
    lambda_upper_multiplier = 512, lambda_min_ratio = 0.05,
    analytic_penalty_limit = TRUE,
    selection_policy = "validation_log_loss_only_exact_ties_alpha_d_ascending_limit_first_lambda_descending",
    target_mode = "groupwise_firth_no_fallback", target_max_iterations = 50L,
    target_tolerance = 1e-7, max_passes = if (smoke) 500L else 1000L,
    max_inner = 2500L, tolerance = 1e-6, inner_tolerance = 1e-8,
    full_path_kkt_limit = 2.05e-6, benchmark_alpha_grid = alpha,
    benchmark_nlambda = 30L, benchmark_lambda_min_ratio = 0.05,
    adelie_tolerance = 1e-7, adelie_max_iterations = 100000L,
    adelie_irls_tolerance = 1e-7, adelie_irls_max_iterations = 10000L,
    grpreg_tolerance = 1e-7, grpreg_max_iterations = 1000000L,
    grpreg_penalty_specification = default_grpreg_penalty_specification_v2(),
    selection_threshold = 1e-8, response_threshold = 0.5,
    mcc_implementation = "mltools::mcc_zero_denominator_zero_with_degeneracy_flag",
    average_precision_definition = "noninterpolated_tied_score_blocks",
    reconstruction_tolerance = 2e-6, endpoint_tolerance = 1e-10,
    tuning_sample = "independent_validation_only",
    test_sample = "evaluation_only_never_tuning")
}

lsg_design_for_stage_v7 <- function(root, stage) {
  full <- lsg_read_design_v7(file.path(root, "config",
                                      "logistic_sglasso_prediction_selection_v7.csv"))
  if (stage == "pilot") return(full)
  toy <- full[2L, , drop = FALSE]
  toy$design_id <- "smoke"
  toy$scenario <- "smoke_homogeneous_rhob_0.3"
  toy$n_train <- 60L; toy$n_validation <- 60L; toy$n_test <- 120L
  toy$groups <- 8L; toy$active_group_count <- 2L
  toy$scenario_index <- 1L
  toy$description <- "Tiny executable smoke only; not scientific simulation evidence"
  rownames(toy) <- NULL
  toy
}

lsg_task_grid_v7 <- function(design, configuration) {
  tasks <- do.call(rbind, lapply(seq_len(nrow(design)), function(i) {
    rep <- seq_len(configuration$replications_per_scenario)
    data.frame(scenario_index = design$scenario_index[i], scenario = design$scenario[i],
      replication = rep, seed = as.integer(configuration$seed_base +
        design$scenario_index[i] * configuration$seed_stride + rep),
      stringsAsFactors = FALSE)
  }))
  tasks$task_id <- seq_len(nrow(tasks))
  tasks$key <- paste(tasks$scenario, tasks$replication, sep = "::")
  tasks$shard_file <- sprintf("shard_%04d.rds", tasks$task_id)
  tasks <- tasks[c("task_id", "scenario_index", "scenario", "replication", "seed",
                   "key", "shard_file")]
  rownames(tasks) <- NULL
  lsg_assert_v7(!anyDuplicated(tasks$seed) && !anyDuplicated(tasks$key),
                "Duplicate frozen task seed/key.")
  tasks
}

lsg_spec_signature_v7 <- function(spec) {
  lsg_hash_v7(spec[c("schema_version", "version", "stage", "configuration", "design",
                     "tasks", "runtime", "source_manifest")])
}

lsg_make_spec_v7 <- function(root, stage, version) {
  configuration <- lsg_configuration_v7(stage)
  design <- lsg_design_for_stage_v7(root, stage)
  spec <- list(schema_version = "prediction_selection_v7", version = version,
    stage = stage, configuration = configuration, design = design,
    tasks = lsg_task_grid_v7(design, configuration), runtime = lsg_runtime_v7(),
    source_manifest = lsg_inventory_v7(root, lsg_source_files_v7()),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  spec$scientific_signature <- lsg_spec_signature_v7(spec)
  spec
}

lsg_validate_spec_v7 <- function(spec, root, stage, version, resume = FALSE) {
  lsg_assert_v7(identical(spec$schema_version, "prediction_selection_v7") &&
    identical(spec$stage, stage) && identical(spec$version, version) &&
    identical(spec$scientific_signature, lsg_spec_signature_v7(spec)),
    "Metadata identity/scientific signature mismatch.")
  lsg_assert_v7(lsg_equal_v7(spec$source_manifest,
    lsg_inventory_v7(root, lsg_source_files_v7()), 0),
    "Source files changed; do not resume or rebaseline old V7 outputs.")
  lsg_assert_v7(lsg_equal_v7(spec$configuration, lsg_configuration_v7(stage), 1e-13) &&
    lsg_equal_v7(spec$design, lsg_design_for_stage_v7(root, stage), 1e-13) &&
    identical(spec$tasks, lsg_task_grid_v7(spec$design, spec$configuration)),
    "Frozen design, numerical controls or task/seed grid changed.")
  if (resume) lsg_assert_v7(identical(spec$runtime, lsg_runtime_v7()),
    "Run/resume requires identical R, platform, RNG and package versions.")
  invisible(TRUE)
}

lsg_output_v7 <- function(root, version) {
  lsg_assert_v7(length(version) == 1L && grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", version),
    "Version must start with an alphanumeric and contain only letters, digits, _ or -.")
  lsg_path_v7(root, file.path("outputs", "study", version))
}

lsg_with_rng_v7 <- function(seed, expression) {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  force(expression)
}

lsg_payload_checks_v7 <- function(payload, task, scenario, configuration,
                                  exact_replay = TRUE) {
  lsg_assert_v7(is.list(payload) && is.data.frame(payload$results) &&
    is.list(payload$artifacts), paste("Missing payload artifacts:", task$key))
  result <- payload$results
  artifact <- payload$artifacts
  checks <- lsg_validate_tuning_payload_v7(payload, configuration)
  lsg_assert_v7(is.data.frame(checks) && all(c("check", "passed") %in% names(checks)),
    "Invalid numerical-check schema.")
  checks <- checks[c("check", "passed")]
  add <- function(name, passed) {
    checks[nrow(checks) + 1L, ] <<- list(name, isTRUE(passed))
  }
  methods <- configuration$methods
  add("exactly_six_unique_methods", nrow(result) == length(methods) &&
      !anyDuplicated(result$method) && setequal(result$method, methods) &&
      !anyDuplicated(names(artifact$models)) && setequal(names(artifact$models), methods))
  add("task_identity", all(result$scenario == task$scenario) &&
      all(result$replication == task$replication) && all(result$seed == task$seed))
  lsg_assert_v7(all(tail(checks$passed, 2L)), "Selected-result identity mismatch.")
  data <- lsg_with_rng_v7(task$seed,
    simulate_logistic_sglasso_two_design_v1(scenario, task$seed))
  true_coefficient <- c(data$intercept, data$beta)
  true_eta <- true_coefficient[1L] + drop(data$X_test %*% true_coefficient[-1L])
  add("saved_test_labels_and_order", identical(as.numeric(artifact$y_test),
      as.numeric(data$y_test)) && identical(as.integer(artifact$observation_id),
        seq_len(scenario$n_test)))
  add("saved_validation_labels_and_order", identical(as.numeric(artifact$y_validation),
      as.numeric(data$y_validation)) &&
      identical(as.integer(artifact$validation_observation_id), seq_len(scenario$n_validation)))
  add("saved_training_labels", identical(as.numeric(artifact$y_train), as.numeric(data$y_train)))
  add("saved_training_sample_size", identical(as.integer(artifact$training_sample_size),
                                             as.integer(scenario$n_train)))
  for (sample in c("train", "validation", "test")) {
    y <- data[[paste0("y_", sample)]]
    add(paste0(sample, "_event_counts_recomputed"),
      all(result[[paste0(sample, "_event_count")]] == sum(y)) &&
      all(abs(result[[paste0(sample, "_observed_prevalence")]] - mean(y)) <= 1e-12))
  }
  add("saved_truth_group_mapping", lsg_equal_v7(artifact$true_coefficient, true_coefficient) &&
      identical(as.character(artifact$group), as.character(data$group)) &&
      setequal(as.character(artifact$active_groups), as.character(data$active_groups)))
  add("true_support_agrees_with_active_groups", setequal(unique(as.character(
      data$group[data$beta != 0])), as.character(data$active_groups)))
  add("true_test_predictions_reconstruct", lsg_equal_v7(artifact$true_test_eta, true_eta) &&
      lsg_equal_v7(artifact$true_test_probability, stats::plogis(true_eta)))
  regenerated_fingerprints <- list(
    training = lsg_data_fingerprint_v7(list(X = data$X_train, y = data$y_train, group = data$group)),
    validation = lsg_data_fingerprint_v7(list(X = data$X_validation, y = data$y_validation)),
    test = lsg_data_fingerprint_v7(list(X = data$X_test, y = data$y_test)),
    truth = lsg_data_fingerprint_v7(list(coefficient = true_coefficient, group = data$group,
                                        active_groups = data$active_groups)))
  fingerprint_schema <- is.list(artifact$data_sha256) &&
    identical(names(artifact$data_sha256), names(regenerated_fingerprints)) &&
    all(vapply(artifact$data_sha256, function(x) is.character(x) && length(x) == 1L &&
      !is.na(x) && grepl("^[a-f0-9]{64}$", x), logical(1)))
  add("portable_fingerprint_schema_and_saved_truth", fingerprint_schema &&
    identical(artifact$data_sha256$truth, lsg_data_fingerprint_v7(list(
      coefficient = artifact$true_coefficient, group = artifact$group,
      active_groups = artifact$active_groups))))
  if (exact_replay) add("portable_data_fingerprints_recomputed",
                        identical(artifact$data_sha256, regenerated_fingerprints))
  for (method in methods) {
    row <- result[result$method == method, , drop = FALSE]
    model <- artifact$models[[method]]
    coefficient <- as.numeric(model$coefficient)
    prefix <- paste0(match(method, methods), "_")
    add(paste0(prefix, "artifact_selected_candidate_identity"),
      identical(model$candidate_id, row$selected_candidate_id) &&
      identical(model$point_type, row$selected_point_type) &&
      lsg_equal_v7(model$alpha, row$selected_alpha) && lsg_equal_v7(model$d, row$selected_d) &&
      lsg_equal_v7(model$lambda, row$selected_lambda) &&
      lsg_equal_v7(model$lambda_relative, row$selected_lambda_relative_to_reference))
    dimension_ok <- length(coefficient) == length(true_coefficient) &&
      all(is.finite(coefficient)) && length(model$probability) == scenario$n_test &&
      length(model$validation_probability) == scenario$n_validation
    add(paste0(prefix, "artifact_dimensions"), dimension_ok)
    if (!dimension_ok) next
    eta <- coefficient[1L] + drop(data$X_test %*% coefficient[-1L])
    validation_eta <- coefficient[1L] + drop(data$X_validation %*% coefficient[-1L])
    reconstruction_error <- max(abs(c(model$probability - stats::plogis(eta),
      model$validation_probability - stats::plogis(validation_eta))))
    add(paste0(prefix, "coefficient_probability_reconstruction"),
        is.finite(reconstruction_error) && reconstruction_error <= 2e-6 &&
        lsg_equal_v7(model$eta, eta, 2e-6) && lsg_equal_v7(model$validation_eta, validation_eta, 2e-6))
    loss <- binary_log_loss(data$y_validation, model$validation_probability)
    add(paste0(prefix, "validation_loss_recomputed"),
        lsg_equal_v7(loss, row$validation_log_loss, 1e-9))
    # Recompute reported metrics from saved evidence. Regeneration is a
    # separate discrete-identity and numerical-reconstruction audit.
    metrics <- lsg_metrics_v7(artifact$y_test, model$probability, coefficient,
      artifact$true_coefficient, artifact$group, model$eta, artifact$true_test_eta,
      threshold = configuration$selection_threshold)
    add(paste0(prefix, "every_metric_recomputed"),
        all(names(metrics) %in% names(row)) &&
        lsg_equal_v7(row[names(metrics)], metrics, 1e-8))
  }
  rownames(checks) <- NULL
  attr(checks, "fingerprint_replay") <- data.frame(
    component = names(regenerated_fingerprints),
    exact_match = vapply(names(regenerated_fingerprints), function(name) {
      identical(artifact$data_sha256[[name]], regenerated_fingerprints[[name]])
    }, logical(1)), stringsAsFactors = FALSE)
  checks
}

lsg_shard_receipt_v7 <- function(shard, path) {
  list(schema_version = "prediction_selection_shard_receipt_v7",
    key = shard$task$key, scientific_signature = shard$scientific_signature,
    file = basename(path), bytes = as.numeric(file.info(path)$size),
    sha256 = lsg_file_hash_v7(path))
}

lsg_read_shard_v7 <- function(output, spec, task, deep = TRUE) {
  file <- file.path("shards", task$shard_file)
  path <- lsg_path_v7(output, file)
  receipt_path <- lsg_path_v7(output, paste0(file, ".receipt.rds"))
  lsg_assert_v7(file.exists(path) && file.exists(receipt_path),
                paste("Incomplete shard publication; preserve and inspect:", task$key))
  shard <- readRDS(path)
  lsg_assert_v7(identical(readRDS(receipt_path), lsg_shard_receipt_v7(shard, path)) &&
    identical(shard$schema_version, "prediction_selection_shard_v7") &&
    identical(shard$scientific_signature, spec$scientific_signature) &&
    identical(shard$task, task), paste("Shard hash/signature/key mismatch:", task$key))
  shard$archival_checks <- shard$checks
  shard$verification_mode <- if (identical(spec$runtime, lsg_runtime_v7()))
    "exact_replay" else "cross_runtime_artifact_audit"
  if (deep) {
    scenario <- spec$design[spec$design$scenario_index == task$scenario_index, , drop = FALSE]
    strict <- identical(shard$verification_mode, "exact_replay")
    checks <- lsg_payload_checks_v7(shard$payload, task, scenario, spec$configuration,
                                    exact_replay = strict)
    archived <- shard$checks
    replay_gate <- archived$check == "portable_data_fingerprints_recomputed"
    lsg_assert_v7(is.data.frame(archived) && !anyNA(archived$check) &&
      !anyDuplicated(archived$check) && sum(replay_gate) == 1L,
      "Archival gate schema is missing or duplicates its exact replay check.")
    if (!strict) {
      lsg_assert_v7(isTRUE(archived$passed[replay_gate]),
                    "The original execution never passed exact data replay.")
      archived <- archived[!replay_gate, ]
    }
    lsg_assert_v7(identical(as.character(checks$check), as.character(archived$check)) &&
      identical(as.logical(checks$passed), as.logical(archived$passed)),
      paste("Stored gates disagree with artifact recomputation:", task$key))
    shard$checks <- checks
  }
  shard
}

lsg_collect_v7 <- function(output, spec, deep = TRUE, require_all = TRUE) {
  tasks <- spec$tasks
  expected <- as.vector(rbind(tasks$shard_file, paste0(tasks$shard_file, ".receipt.rds")))
  listed <- list.files(file.path(output, "shards"), all.files = TRUE,
                       no.. = TRUE, recursive = TRUE)
  # Interrupted atomic publications can leave only their own temporary file.
  unexpected <- setdiff(listed, expected)
  lsg_assert_v7(!length(unexpected), paste("Unexpected/incomplete shard files; inspect:",
                                         paste(unexpected, collapse = ", ")))
  present <- file.exists(file.path(output, "shards", tasks$shard_file))
  receipts_present <- file.exists(file.path(output, "shards", paste0(tasks$shard_file, ".receipt.rds")))
  lsg_assert_v7(identical(present, receipts_present),
                "Asymmetric shard/receipt pair; inspect the interrupted publication before resume.")
  if (require_all) lsg_assert_v7(all(present),
    sprintf("Study incomplete: %d/%d task shards are present.", sum(present), nrow(tasks)))
  rows <- gates <- archival_gates <- replay <- vector("list", sum(present))
  j <- 0L
  for (i in which(present)) {
    j <- j + 1L
    shard <- lsg_read_shard_v7(output, spec, tasks[i, , drop = FALSE], deep)
    rows[[j]] <- shard$payload$results
    gates[[j]] <- cbind(tasks[i, c("task_id", "scenario", "replication")], shard$checks,
                       row.names = NULL)
    archival_gates[[j]] <- cbind(tasks[i, c("task_id", "scenario", "replication")],
                                shard$archival_checks, row.names = NULL)
    replay[[j]] <- cbind(tasks[i, c("task_id", "scenario", "replication")],
                         attr(shard$checks, "fingerprint_replay"), row.names = NULL)
  }
  list(results = if (length(rows)) do.call(rbind, rows) else data.frame(),
       checks = if (length(gates)) do.call(rbind, gates) else data.frame(),
       archival_checks = if (length(archival_gates)) do.call(rbind, archival_gates) else data.frame(),
       fingerprint_replay = if (length(replay)) do.call(rbind, replay) else data.frame(),
       present = present)
}

lsg_summary_frames_v7 <- function(results) {
  registry <- lsg_metric_registry_v7()
  lsg_assert_v7(!anyDuplicated(registry$metric) &&
    all(registry$metric %in% names(results)), "Summary metric schema mismatch.")
  moment <- function(value) {
    finite <- is.finite(value)
    value_finite <- as.numeric(value[finite])
    n <- length(value_finite)
    se <- if (n > 1L) stats::sd(value_finite) / sqrt(n) else NA_real_
    center <- if (n) mean(value_finite) else NA_real_
    width <- if (n > 1L) stats::qt(0.975, n - 1L) * se else NA_real_
    data.frame(n_total = length(value), n_defined = n,
      n_undefined = sum(!finite), mean = center, mcse = se,
      interval_low = center - width, interval_high = center + width)
  }
  summaries <- reasons <- paired <- list()
  k <- u <- j <- 0L
  for (scenario in unique(results$scenario)) {
    cell <- results[results$scenario == scenario, , drop = FALSE]
    for (method in unique(cell$method)) {
      method_rows <- cell[cell$method == method, , drop = FALSE]
      for (m in seq_len(nrow(registry))) {
        metric <- registry$metric[m]
        value <- method_rows[[metric]]
        k <- k + 1L
        summaries[[k]] <- cbind(data.frame(scenario = scenario, method = method),
          registry[m, , drop = FALSE], moment(value))
        reason_field <- paste0(metric, "_undefined_reason")
        why <- if (reason_field %in% names(method_rows)) method_rows[[reason_field]] else
          ifelse(is.finite(value), "", "nonfinite_without_reason")
        if (any(!is.finite(value))) {
          counts <- table(why[!is.finite(value)], useNA = "ifany")
          u <- u + 1L
          reasons[[u]] <- data.frame(scenario = scenario, method = method,
            metric = metric, reason = names(counts), count = as.integer(counts))
        }
      }
    }
    reference <- cell[cell$method == "Logistic SGLASSO", , drop = FALSE]
    for (other in setdiff(unique(cell$method), "Logistic SGLASSO")) {
      comparator <- cell[cell$method == other, , drop = FALSE]
      index <- match(reference$replication, comparator$replication)
      lsg_assert_v7(!anyNA(index) && !anyDuplicated(comparator$replication) &&
        identical(reference$seed, comparator$seed[index]), "Unpaired comparison detected.")
      for (m in seq_len(nrow(registry))) {
        metric <- registry$metric[m]
        difference <- reference[[metric]] - comparator[[metric]][index]
        finite <- is.finite(difference)
        direction <- registry$direction[m]
        better <- if (direction == "lower") difference < 0 else
          if (direction == "higher") difference > 0 else rep(NA, length(difference))
        j <- j + 1L
        paired[[j]] <- cbind(data.frame(scenario = scenario,
          contrast = paste("Logistic SGLASSO minus", other),
          metric = metric, direction = direction), moment(difference),
          data.frame(better_rate = if (any(finite) && direction %in% c("lower", "higher"))
            mean(better[finite]) else NA_real_,
            tie_rate = if (any(finite)) mean(difference[finite] == 0) else NA_real_))
      }
    }
  }
  keys <- paste(results$scenario, results$method, sep = "::")
  tuning <- do.call(rbind, lapply(unique(keys), function(key) {
      x <- results[keys == key, , drop = FALSE]
      data.frame(scenario = x$scenario[1L], method = x$method[1L], n = nrow(x),
        alpha_zero = sum(x$selected_alpha == 0, na.rm = TRUE),
        d_positive = sum(x$selected_d > 0, na.rm = TRUE),
        penalty_limit = sum(x$selected_point_type == "penalty_limit"),
        all_groups_selected = sum(x$group_selected_all),
        no_groups_selected = sum(x$group_selected_none),
        response_mcc_degenerate = sum(x$y_mcc_degenerate),
        group_mcc_degenerate = sum(x$group_mcc_degenerate))
    }))
  rows <- results[!duplicated(paste(results$scenario, results$replication)), , drop = FALSE]
  runtime <- do.call(rbind, lapply(unique(rows$scenario), function(scenario) {
    x <- rows[rows$scenario == scenario, , drop = FALSE]
    data.frame(scenario = x$scenario[1L], tasks = nrow(x),
      task_seconds_median = stats::median(x$task_runtime_seconds),
      task_seconds_max = max(x$task_runtime_seconds))
  }))
  list(metric_summary = do.call(rbind, summaries), paired_differences = do.call(rbind, paired),
    undefined_metrics = if (length(reasons)) do.call(rbind, reasons) else
      data.frame(scenario = character(), method = character(), metric = character(),
                  reason = character(), count = integer()),
    selection_diagnostics = tuning, runtime_summary = runtime,
    metric_registry = registry)
}

lsg_manifest_files_v7 <- function(output, version) {
  files <- list.files(output, all.files = TRUE, no.. = TRUE, recursive = TRUE)
  files <- files[!grepl("^\\.run_lock(/|$)", files) &
    files != paste0("output_manifest_", version, ".csv")]
  sort(files, method = "radix")
}

lsg_verify_manifest_v7 <- function(output, version) {
  path <- lsg_path_v7(output, paste0("output_manifest_", version, ".csv"))
  lsg_assert_v7(file.exists(path), "Missing complete output checksum manifest.")
  manifest <- utils::read.csv(path, stringsAsFactors = FALSE)
  lsg_assert_v7(identical(names(manifest), c("file", "bytes", "sha256")) &&
    !anyNA(manifest) && !anyDuplicated(manifest$file) &&
    identical(manifest$file, lsg_manifest_files_v7(output, version)),
    "Output manifest inventory differs from actual files.")
  actual <- lsg_inventory_v7(output, manifest$file)
  lsg_assert_v7(lsg_equal_v7(actual, manifest, 0), "Output SHA-256/byte manifest failed.")
  invisible(manifest)
}

lsg_finalize_v7 <- function(output, spec, collected) {
  lsg_assert_v7(all(collected$present) && nrow(collected$checks) > 0L &&
    all(collected$checks$passed), "Scientific/numerical gates failed; no accepted final output.")
  summaries <- lsg_summary_frames_v7(collected$results)
  exports <- c(list(selected_results = collected$results, validation_checks = collected$checks),
               summaries)
  for (name in names(exports)) lsg_atomic_v7(exports[[name]], output,
    file.path("final", paste0(name, "_", spec$version, ".csv")), "csv", TRUE)
  receipt <- list(schema_version = "prediction_selection_acceptance_v7",
    version = spec$version, scientific_signature = spec$scientific_signature,
    source_signature = lsg_hash_v7(spec$source_manifest), runtime = spec$runtime,
    task_count = nrow(spec$tasks), selected_rows = nrow(collected$results),
    all_gates_passed = TRUE, completed_utc = spec$created_utc)
  # Time above is the immutable specification timestamp, not a claimed finish
  # timestamp; the completed job times are retained in each attempt record.
  names(receipt)[names(receipt) == "completed_utc"] <- "specification_created_utc"
  if (spec$stage == "smoke") {
    unit_file <- lsg_path_v7(output, "unit_checks.rds")
    lsg_assert_v7(file.exists(unit_file) && all(readRDS(unit_file)$passed),
                  "Smoke unit checks are missing or failed.")
    lsg_atomic_v7(receipt, output, paste0("SMOKE_ACCEPTED_", spec$version, ".rds"),
                  identical_ok = TRUE)
  }
  lsg_atomic_v7(receipt, output, file.path("final", paste0("ACCEPTED_", spec$version, ".rds")),
                identical_ok = TRUE)
  marker <- c(paste("version", spec$version), paste("stage", spec$stage),
    paste("scientific_signature", spec$scientific_signature),
    paste("completed_tasks", nrow(spec$tasks)), "all_gates_passed TRUE")
  lsg_atomic_v7(marker, output, paste0("COMPLETED_", spec$version, ".txt"), "text", TRUE)
  manifest <- lsg_inventory_v7(output, lsg_manifest_files_v7(output, spec$version))
  lsg_atomic_v7(manifest, output, paste0("output_manifest_", spec$version, ".csv"), "csv", TRUE)
  invisible(receipt)
}

lsg_verify_v7 <- function(root, stage, version, quiet = FALSE) {
  output <- lsg_output_v7(root, version)
  lsg_assert_v7(dir.exists(output), paste("No V7 output directory:", output))
  lsg_verify_manifest_v7(output, version)
  spec <- readRDS(lsg_path_v7(output, "study_specification.rds"))
  lsg_validate_spec_v7(spec, root, stage, version)
  mode <- if (identical(spec$runtime, lsg_runtime_v7())) "exact_replay" else
    "cross_runtime_artifact_audit"
  collected <- lsg_collect_v7(output, spec, deep = TRUE)
  lsg_assert_v7(nrow(collected$results) == nrow(spec$tasks) * 6L &&
    nrow(collected$checks) > 0L && all(collected$checks$passed) &&
    nrow(collected$archival_checks) > 0L && all(collected$archival_checks$passed),
    "Not all tasks, selected results and scientific/numerical gates passed.")
  accepted <- readRDS(lsg_path_v7(output, file.path("final", paste0("ACCEPTED_", version, ".rds"))))
  lsg_assert_v7(identical(accepted$scientific_signature, spec$scientific_signature) &&
    identical(accepted$source_signature, lsg_hash_v7(spec$source_manifest)) &&
    identical(accepted$runtime, spec$runtime) && isTRUE(accepted$all_gates_passed) &&
    identical(accepted$task_count, nrow(spec$tasks)) &&
    identical(accepted$selected_rows, nrow(collected$results)), "Invalid acceptance receipt.")
  if (stage == "smoke") {
    lsg_assert_v7(identical(accepted, readRDS(lsg_path_v7(output,
      paste0("SMOKE_ACCEPTED_", version, ".rds")))) &&
      all(readRDS(lsg_path_v7(output, "unit_checks.rds"))$passed),
      "Smoke receipt or unit-check results failed.")
  }
  expected_marker <- c(paste("version", version), paste("stage", stage),
    paste("scientific_signature", spec$scientific_signature),
    paste("completed_tasks", nrow(spec$tasks)), "all_gates_passed TRUE")
  lsg_assert_v7(identical(readLines(lsg_path_v7(output, paste0("COMPLETED_", version, ".txt"))),
    expected_marker), "Completion marker content mismatch.")
  expected <- c(list(selected_results = collected$results, validation_checks = collected$archival_checks),
                 lsg_summary_frames_v7(collected$results))
  for (name in names(expected)) {
    target <- expected[[name]]
    column_classes <- vapply(target, function(column) {
      if (is.character(column)) "character" else if (is.logical(column)) "logical" else
        if (is.integer(column)) "integer" else if (is.numeric(column)) "numeric" else
          stop("Unsupported exported column type.", call. = FALSE)
    }, character(1))
    actual <- utils::read.csv(lsg_path_v7(output, file.path("final",
      paste0(name, "_", version, ".csv"))), stringsAsFactors = FALSE,
      check.names = FALSE, na.strings = "NA", colClasses = column_classes)
    # Explicit schema preserves all-blank strings and all-NA numeric columns.
    lsg_assert_v7(identical(names(actual), names(target)) && nrow(actual) == nrow(target),
                  paste("Final table shape mismatch:", name))
    for (column in names(target)) {
      lsg_assert_v7(lsg_equal_v7(actual[[column]], target[[column]], 1e-8),
                    paste("Final table cannot be recomputed:", name, column))
    }
  }
  if (!quiet) {
    cat("Verification mode:", mode, "\n")
    if (mode == "exact_replay") cat("Artifact, exact data replay, metric, raw tuning and SHA-256 verification: TRUE\n")
    else {
      cat("Original execution integrity and cross-runtime artifact/metric audit: TRUE\n")
      cat("Byte-identical data replay: NOT ASSERTED (different execution runtime)\n")
      cat("Observed regenerated fingerprint matches:", sum(collected$fingerprint_replay$exact_match),
          "/", nrow(collected$fingerprint_replay), "(diagnostic only in this mode)\n")
    }
    cat("Completed tasks:", nrow(spec$tasks), "/", nrow(spec$tasks), "\n")
    cat("Selected method rows:", nrow(collected$results), "\n")
    cat("Original R:", spec$runtime$r_version, "\nVerifier R:", R.version.string, "\n")
    cat("Output directory:", output, "\nComplete: TRUE\n")
  }
  invisible(list(specification = spec, results = collected$results, checks = collected$checks,
                 verification_mode = mode, fingerprint_replay = collected$fingerprint_replay))
}

lsg_release_lock_v7 <- function(lock, owner) {
  if (!dir.exists(lock)) return(invisible(TRUE))
  link <- Sys.readlink(lock)
  contents <- list.files(lock, all.files = TRUE, no.. = TRUE)
  safe <- (is.na(link) || !nzchar(link)) && basename(lock) == ".run_lock" &&
    (length(contents) == 0L || identical(contents, "owner.rds") &&
      isTRUE(tryCatch(identical(readRDS(file.path(lock, "owner.rds")), owner),
                      error = function(e) FALSE)))
  if (safe) {
    # The exact directory was exclusively created by this invocation; no
    # unknown files or another owner's receipt may be removed.
    unlink(lock, recursive = TRUE)
  } else warning("Lock ownership/content changed; leaving it untouched: ", lock)
  invisible(safe && !dir.exists(lock))
}

lsg_io_unit_checks_v7 <- function() {
  temporary <- tempfile("lsg_v7_io_unit_", tmpdir = tempdir())
  dir.create(temporary)
  # Only this newly owned fixture directory is removed.
  on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
  rejects <- function(code) inherits(try(force(code), silent = TRUE), "try-error")
  sample <- list(integer = 1:3, value = c(0, 0.1, Inf, NA_real_),
                 labels = c("a", "b"), nested = list(flag = TRUE))
  hash <- lsg_hash_v7(sample)
  lsg_atomic_v7(sample, temporary, "sample.rds")
  same <- lsg_file_hash_v7(file.path(temporary, "sample.rds"))
  no_clobber <- rejects(lsg_atomic_v7(list(changed = TRUE), temporary, "sample.rds"))
  manifest <- lsg_inventory_v7(temporary, "sample.rds")
  lsg_atomic_v7("intentionally corrupted fixture", temporary, "different.txt", "text")
  original <- readRDS(file.path(temporary, "sample.rds"))
  lsg_atomic_v7(sample, temporary, "sample.rds", identical_ok = TRUE)
  lock <- file.path(temporary, ".run_lock")
  dir.create(lock)
  owner <- list(pid = Sys.getpid(), fixture = TRUE)
  lsg_atomic_v7(owner, temporary, ".run_lock/owner.rds")
  lock_released <- lsg_release_lock_v7(lock, owner)
  csv_fixture <- data.frame(reason = c("", ""), missing = rep(NA_real_, 2L))
  lsg_atomic_v7(csv_fixture, temporary, "schema.csv", "csv")
  csv_restored <- utils::read.csv(file.path(temporary, "schema.csv"),
    colClasses = c("character", "numeric"), na.strings = "NA")
  checks <- c(
    portable_signature_roundtrip = identical(hash, lsg_hash_v7(original)),
    typed_signature_distinguishes_integer_double =
      !identical(lsg_hash_v7(1L), lsg_hash_v7(1)),
    signature_detects_value_change = !identical(hash, lsg_hash_v7(c(sample, extra = 1))),
    portable_data_fingerprint_roundtrip = identical(lsg_data_fingerprint_v7(sample),
      lsg_data_fingerprint_v7(original)),
    data_fingerprint_preserves_dimensions = !identical(
      lsg_data_fingerprint_v7(matrix(1:6, 2L)), lsg_data_fingerprint_v7(matrix(1:6, 3L))),
    atomic_no_clobber = no_clobber &&
      identical(same, lsg_file_hash_v7(file.path(temporary, "sample.rds"))),
    byte_manifest_matches = identical(manifest$sha256, same) &&
      manifest$bytes == file.info(file.path(temporary, "sample.rds"))$size,
    changed_file_hash_detected = !identical(same,
      lsg_file_hash_v7(file.path(temporary, "different.txt"))),
    owned_lock_removed_for_resume = lock_released && !dir.exists(lock),
    csv_empty_strings_and_numeric_na_preserved = identical(csv_fixture, csv_restored),
    relative_parent_path_rejected = rejects(lsg_path_v7(temporary, "../outside.rds")),
    absolute_path_rejected = rejects(lsg_path_v7(temporary, "/outside.rds")))
  link <- file.path(temporary, "linked")
  if (.Platform$OS.type == "unix" && file.symlink(tempdir(), link)) {
    checks["symlink_path_rejected"] <- rejects(lsg_path_v7(temporary, "linked/a.rds"))
    unlink(link)
  }
  data.frame(check = names(checks), passed = unname(checks))
}

lsg_unit_checks_v7 <- function(root) {
  design <- lsg_read_design_v7(file.path(root, "config",
    "logistic_sglasso_prediction_selection_v7.csv"))
  checks <- list(lsg_design_checks_v7(design), lsg_metrics_unit_checks_v7(),
                 lsg_io_unit_checks_v7(), lsg_tuning_unit_checks_v7())
  checks <- do.call(rbind, lapply(checks, function(x) x[c("check", "passed")]))
  pilot <- lsg_configuration_v7("pilot")
  tasks <- lsg_task_grid_v7(design, pilot)
  smoke_tasks <- lsg_task_grid_v7(lsg_design_for_stage_v7(root, "smoke"),
                                 lsg_configuration_v7("smoke"))
  additions <- c(eight_by_twenty_tasks = nrow(tasks) == 160L &&
      all(table(tasks$scenario) == 20L),
    independent_new_seed_ranges = !any(tasks$seed %in% smoke_tasks$seed) &&
      min(tasks$seed) > 120760876L,
    exact_36_finite_lambda_grid = length(pilot$lambda_relative_grid) == 36L &&
      all(diff(pilot$lambda_relative_grid) < 0),
    v5_33_point_tail_preserved = identical(pilot$lambda_relative_grid[-(1:3)],
      c(64, 32, 16, exp(seq(log(8), log(0.05), length.out = 30L)))))
  rbind(checks, data.frame(check = names(additions), passed = unname(additions)))
}

lsg_run_v7 <- function(root, stage, version, cores, max_seconds) {
  start <- proc.time()[["elapsed"]]
  lsg_assert_v7(.Platform$OS.type == "unix" || cores == 1L,
                "Multicore V7 tasks require a Unix fork implementation.")
  if (nzchar(Sys.getenv("SLURM_CPUS_PER_TASK"))) lsg_assert_v7(
    cores <= as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")), "Workers exceed allocated CPUs.")
  if (stage == "pilot") {
    smoke <- lsg_verify_v7(root, "smoke", lsg_default_version_v7("smoke"), quiet = TRUE)
    lsg_assert_v7(identical(smoke$specification$runtime, lsg_runtime_v7()),
                  "Pilot requires a passed smoke with identical runtime/packages.")
  }
  output <- lsg_output_v7(root, version)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  complete_manifest <- file.path(output, paste0("output_manifest_", version, ".csv"))
  if (file.exists(complete_manifest)) {
    existing <- readRDS(lsg_path_v7(output, "study_specification.rds"))
    lsg_validate_spec_v7(existing, root, stage, version, resume = TRUE)
    result <- lsg_verify_v7(root, stage, version)
    lsg_assert_v7(identical(result$specification$runtime, lsg_runtime_v7()),
                  "An already complete run has a different runtime; use action=verify only.")
    cat("Valid completed output reused; no models refitted or files changed.\n")
    return(invisible(result))
  }
  lock <- lsg_path_v7(output, ".run_lock")
  lsg_assert_v7(dir.create(lock, showWarnings = FALSE), paste(
    "Run lock exists. Check its owner and active jobs before removing a stale lock:", lock))
  owner <- list(pid = Sys.getpid(), slurm_job = Sys.getenv("SLURM_JOB_ID"),
                host = Sys.info()[["nodename"]], started_utc = format(Sys.time(), tz = "UTC"))
  on.exit(lsg_release_lock_v7(lock, owner), add = TRUE)
  lsg_atomic_v7(owner, output, ".run_lock/owner.rds")
  spec_path <- lsg_path_v7(output, "study_specification.rds")
  spec <- if (file.exists(spec_path)) readRDS(spec_path) else {
    spec <- lsg_make_spec_v7(root, stage, version)
    lsg_atomic_v7(spec, output, "study_specification.rds")
    spec
  }
  lsg_validate_spec_v7(spec, root, stage, version, resume = TRUE)
  lsg_atomic_v7(spec$tasks, output, "task_grid.csv", "csv", TRUE)
  lsg_atomic_v7(spec$source_manifest, output, "source_manifest.csv", "csv", TRUE)
  lsg_atomic_v7(spec$design, output, "design.csv", "csv", TRUE)
  if (stage == "smoke") {
    checks <- lsg_unit_checks_v7(root)
    print(checks, row.names = FALSE)
    lsg_assert_v7(all(checks$passed), "V7 unit checks failed; no study fitting started.")
    lsg_atomic_v7(checks, output, "unit_checks.rds", identical_ok = TRUE)
  }
  collected <- lsg_collect_v7(output, spec, deep = TRUE, require_all = FALSE)
  pending <- which(!collected$present)
  if (nrow(collected$checks) && any(!collected$checks$passed)) {
    print(collected$checks[!collected$checks$passed, ], row.names = FALSE)
    stop("Existing shards failed scientific gates. Preserve them; do not silently retune or drop tasks.",
         call. = FALSE)
  }
  cat("Stage:", stage, "\nVersion:", version,
      "\nScientific signature:", spec$scientific_signature,
      "\nPending tasks:", length(pending), "/", nrow(spec$tasks),
      "\nRequested workers:", cores, "\n")
  if (length(pending)) compile_lsg_core(quiet = TRUE)
  attempt <- paste0("attempt_", format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"),
                    "_", Sys.getpid())
  failures <- list()
  completed_now <- integer()
  while (length(pending) && proc.time()[["elapsed"]] - start < max_seconds) {
    wave <- head(pending, cores)
    pending <- pending[-seq_along(wave)]
    worker <- function(i) {
      task <- spec$tasks[i, , drop = FALSE]
      scenario <- spec$design[spec$design$scenario_index == task$scenario_index, , drop = FALSE]
      task_start <- proc.time()[["elapsed"]]
      tryCatch({
        cat("Running", task$key, "\n")
        payload <- lsg_with_rng_v7(task$seed,
          lsg_run_task_v7(scenario, task$replication, task$seed, spec$configuration))
        payload$results$task_runtime_seconds <- proc.time()[["elapsed"]] - task_start
        checks <- lsg_payload_checks_v7(payload, task, scenario, spec$configuration)
        shard <- list(schema_version = "prediction_selection_shard_v7",
          scientific_signature = spec$scientific_signature, task = task,
          payload = payload, checks = checks, attempt = attempt)
        relative <- file.path("shards", task$shard_file)
        path <- lsg_atomic_v7(shard, output, relative)
        lsg_atomic_v7(lsg_shard_receipt_v7(shard, path), output,
                      paste0(relative, ".receipt.rds"))
        list(task_id = i, key = task$key, status = if (all(checks$passed)) "passed" else
          "gate_failed", message = paste(checks$check[!checks$passed], collapse = ";"))
      }, error = function(e) list(task_id = i, key = task$key,
        status = "error", message = conditionMessage(e)))
    }
    outcome <- if (length(wave) > 1L && cores > 1L) {
      parallel::mclapply(wave, worker, mc.cores = min(cores, length(wave)),
        mc.preschedule = FALSE, mc.set.seed = FALSE)
    } else lapply(wave, worker)
    for (k in seq_along(outcome)) {
      item <- outcome[[k]]
      if (!is.list(item) || is.null(item$status)) item <- list(task_id = wave[k],
        key = spec$tasks$key[wave[k]], status = "worker_error", message = as.character(item))
      if (identical(item$status, "passed")) completed_now <- c(completed_now, item$task_id)
      else failures[[length(failures) + 1L]] <- item
    }
    cat("New accepted shards:", length(completed_now), "; queued for later:", length(pending), "\n")
  }
  record <- list(owner = owner, scientific_signature = spec$scientific_signature,
    requested_workers = cores, max_seconds_soft = max_seconds,
    elapsed_seconds = proc.time()[["elapsed"]] - start, completed_now = completed_now,
    unstarted_tasks = pending, failures = failures, runtime = lsg_runtime_v7(),
    library_paths = .libPaths(), session_info = capture.output(utils::sessionInfo()),
    finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
  lsg_atomic_v7(record, output, file.path("attempts", paste0(attempt, ".rds")))
  collected <- lsg_collect_v7(output, spec, deep = TRUE, require_all = FALSE)
  cat("Task shards:", sum(collected$present), "/", nrow(spec$tasks), "\n")
  if (length(failures)) print(failures)
  if (!all(collected$present) || !nrow(collected$checks) || any(!collected$checks$passed)) {
    if (nrow(collected$checks)) print(collected$checks[!collected$checks$passed, ], row.names = FALSE)
    stop(paste("Run incomplete or scientific/numerical gate failed. Valid shards are preserved.",
      "Resume with the identical command; never change a frozen grid under this version."),
      call. = FALSE)
  }
  lsg_finalize_v7(output, spec, collected)
  cat("Output directory:", output, "\nComplete: TRUE\n")
  invisible(spec)
}

lsg_cli_v7 <- function(args, root) {
  allowed <- c("stage", "version", "cores", "max-seconds", "action")
  values <- list()
  for (argument in args) {
    lsg_assert_v7(grepl("^--[^=]+=.+$", argument), "Use --name=value arguments.")
    name <- sub("^--([^=]+)=.*$", "\\1", argument)
    value <- sub("^--[^=]+=", "", argument)
    lsg_assert_v7(name %in% allowed && is.null(values[[name]]),
                  paste("Unknown or repeated argument:", name))
    values[[name]] <- value
  }
  stage <- if (is.null(values$stage)) "smoke" else values$stage
  version <- if (is.null(values$version)) lsg_default_version_v7(stage) else values$version
  action <- if (is.null(values$action)) "run" else values$action
  lsg_assert_v7(stage %in% c("smoke", "pilot") &&
    action %in% c("run", "verify", "summarize", "design", "unit"), "Unsupported stage/action.")
  integer_value <- function(value, default) {
    if (is.null(value)) return(default)
    lsg_assert_v7(grepl("^[1-9][0-9]*$", value) && is.finite(as.numeric(value)) &&
                   as.numeric(value) < .Machine$integer.max, "Invalid positive integer option.")
    as.integer(value)
  }
  cores <- integer_value(values$cores, 1L)
  seconds <- integer_value(values[["max-seconds"]], if (stage == "smoke") 300L else 3000L)
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
    BLIS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
  if (action == "design") {
    configuration <- lsg_configuration_v7(stage)
    design <- lsg_design_for_stage_v7(root, stage)
    print(design, row.names = FALSE)
    tasks <- lsg_task_grid_v7(design, configuration)
    cat("Frozen tasks:", nrow(tasks), "\nSeed range:", range(tasks$seed),
        "\nFinite lambda points:", configuration$nlambda,
        "+ analytic limit per canonical alpha,d\nNo models fitted; no outputs written.\n")
    return(invisible(tasks))
  }
  runtime <- lsg_runtime_v7()
  cat("R:", runtime$r_version, "\nPlatform:", runtime$platform, "\n")
  print(runtime$package_versions)
  cat("Library paths:\n", paste(.libPaths(), collapse = "\n"), "\n")
  if (action == "unit") {
    checks <- lsg_unit_checks_v7(root)
    print(checks, row.names = FALSE)
    lsg_assert_v7(all(checks$passed), "V7 unit checks failed.")
    return(invisible(checks))
  }
  if (action %in% c("verify", "summarize")) {
    verified <- lsg_verify_v7(root, stage, version)
    if (action == "summarize") cat("All existing summaries independently recomputed and verified; no refits or overwrites.\n")
    return(invisible(verified))
  }
  lsg_run_v7(root, stage, version, cores, seconds)
}
