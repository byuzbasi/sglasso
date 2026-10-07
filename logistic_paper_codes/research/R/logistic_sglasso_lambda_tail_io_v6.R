# File/signature infrastructure only. No fitting or V5-output mutation occurs here.

lsg_tail_assert_v6 <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

lsg_tail_constants_v6 <- function() {
  list(
    source_version = "logistic_sglasso_two_design_r50_v5",
    source_scientific_signature =
      "92c663c8cf66cfdf7febaef9e5a18595472cc28b5f8bbb18cb4e4e0e07800717",
    source_tuning_signature =
      "92a3bad483b4d90033336c19d3c8395d3952e8bc0960c30b53d31c9dc4b1eed5",
    diagnostic_version = "logistic_sglasso_lambda_tail_diagnostic_v6",
    smoke_version = "logistic_sglasso_lambda_tail_smoke_v6",
    packages = c("Rcpp", "RcppArmadillo", "digest", "adelie", "grpreg", "logistf"),
    rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"),
    methods = c("Logistic SGLASSO", "Logistic SGLASSO (d=0 boundary)",
                "Logistic Group Elastic Net (adelie)"),
    case_columns = c("case_id", "task_id", "scenario", "replication", "seed",
                     "method", "alpha", "d", "lambda_reference",
                     "validation_log_loss_v5", "source_shard_file")
  )
}

lsg_tail_relative_v6 <- function(paths) {
  out <- sub("^.*?/logistic_prework/", "logistic_prework/", paths)
  valid <- grepl("^logistic_prework/", out) &
    !grepl("(^|/)\\.\\.?(/|$)", out) & !grepl("\\\\", out)
  lsg_tail_assert_v6(all(valid), "Unsafe or unmappable project-relative path.")
  out
}

lsg_tail_manifest_v6 <- function(root, relative_paths) {
  relative_paths <- sort(unique(lsg_tail_relative_v6(relative_paths)))
  paths <- file.path(dirname(normalizePath(root, mustWork = TRUE)), relative_paths)
  lsg_tail_assert_v6(all(file.exists(paths)) && !any(dir.exists(paths)),
                     "A required manifest input file is missing.")
  data.frame(file = relative_paths, bytes = as.numeric(file.info(paths)$size),
             sha256 = unname(vapply(paths, digest::digest, character(1),
                                    file = TRUE, algo = "sha256")),
             stringsAsFactors = FALSE)
}

lsg_tail_equal_v6 <- function(x, y, tolerance = 0) {
  isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE))
}

lsg_tail_atomic_v6 <- function(value, path, kind, allow_identical = FALSE) {
  directory <- dirname(path)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(".lsg_tail_tmp_", tmpdir = directory)
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  switch(kind,
         rds = saveRDS(value, temporary, version = 3L),
         csv = utils::write.csv(value, temporary, row.names = FALSE, na = "NA"),
         text = writeLines(value, temporary, useBytes = TRUE),
         stop("Unknown atomic output kind.", call. = FALSE))
  if (file.exists(path)) {
    same <- identical(digest::digest(file = path, algo = "sha256"),
                      digest::digest(file = temporary, algo = "sha256"))
    lsg_tail_assert_v6(allow_identical && same,
                       paste("Refusing to replace existing output:", path))
  } else {
    # A same-filesystem hard link atomically publishes without replacing a target.
    lsg_tail_assert_v6(file.link(temporary, path),
                       paste("Atomic no-clobber publication failed:", path))
  }
  invisible(path)
}

lsg_tail_atomic_rds_v6 <- function(value, path, allow_identical = FALSE) {
  lsg_tail_atomic_v6(value, path, "rds", allow_identical)
}

lsg_tail_runtime_v6 <- function() {
  packages <- lsg_tail_constants_v6()$packages
  list(r_version = R.version.string, platform = R.version$platform,
       package_versions = stats::setNames(vapply(packages, function(package) {
         if (requireNamespace(package, quietly = TRUE)) {
           as.character(utils::packageVersion(package))
         } else NA_character_
       }, character(1)), packages), rng_kind = RNGkind())
}

lsg_tail_runtime_matches_v6 <- function(runtime, source_configuration) {
  identical(runtime$r_version, source_configuration$r_version) &&
    identical(runtime$platform, source_configuration$platform) &&
    identical(runtime$package_versions, source_configuration$package_versions) &&
    identical(runtime$rng_kind, lsg_tail_constants_v6()$rng_kind)
}

lsg_tail_v5_signatures_v6 <- function(metadata) {
  configuration <- metadata$configuration
  source_sha256 <- metadata$source_sha256
  tuning_names <- c(
    "logistic_sglasso_core.cpp", "logistic_sglasso.R",
    "groupwise_logistic_targets.R", "correlated_target_study_utils.R",
    "logistic_penalized_benchmarks_v1.R", "sglasso_rhob_benchmark_study_utils_v1.R",
    "logistic_group_benchmarks_v2.R", "sglasso_rhob_group_benchmark_study_utils_v2.R",
    "logistic_sglasso_two_design_study_utils_v1.R",
    "logistic_sglasso_lambda_extension_v5.R")
  # Preserve original TRUBA names/order here: these names entered the V5 hash.
  tuning_signature <- digest::digest(list(
    alpha_grid = configuration$alpha_grid, d_grid = configuration$d_grid,
    nlambda = configuration$nlambda,
    lambda_relative_grid = configuration$lambda_relative_grid,
    lambda_base_nlambda = configuration$lambda_base_nlambda,
    lambda_extension_multipliers = configuration$lambda_extension_multipliers,
    lambda_upper_multiplier = configuration$lambda_upper_multiplier,
    lambda_min_ratio = configuration$lambda_min_ratio,
    max_passes = configuration$max_passes, max_inner = configuration$max_inner,
    tolerance = configuration$tolerance, inner_tolerance = configuration$inner_tolerance,
    full_path_kkt_limit = configuration$full_path_kkt_limit,
    benchmark_alpha_grid = configuration$benchmark_alpha_grid,
    benchmark_nlambda = configuration$benchmark_nlambda,
    benchmark_lambda_min_ratio = configuration$benchmark_lambda_min_ratio,
    numerical_controls = configuration[c(
      "adelie_tolerance", "adelie_max_iterations", "adelie_irls_tolerance",
      "adelie_irls_max_iterations", "grpreg_tolerance", "grpreg_max_iterations",
      "target_max_iterations", "target_tolerance")],
    package_versions = configuration$package_versions,
    r_version = configuration$r_version, platform = configuration$platform,
    source_sha256 = source_sha256[basename(metadata$source_files) %in% tuning_names]
  ), algo = "sha256", serialize = TRUE)
  scientific_signature <- digest::digest(list(
    configuration = configuration, design = metadata$design,
    source_sha256 = source_sha256, tuning_signature = tuning_signature
  ), algo = "sha256", serialize = TRUE)
  list(tuning = tuning_signature, scientific = scientific_signature)
}

lsg_tail_selection_v6 <- function(shard, task, method, which_selected) {
  frame <- if (method == "Logistic Group Elastic Net (adelie)") {
    frame <- shard$payload$external_tuning
    frame[frame$selected & frame$method_path == method, , drop = FALSE]
  } else {
    frame <- shard$payload$sglasso_tuning
    frame[frame[[which_selected]], , drop = FALSE]
  }
  lsg_tail_assert_v6(nrow(frame) == 1L && isTRUE(frame$numerically_eligible) &&
                     is.finite(frame$validation_log_loss),
                     paste("Invalid V5 selected candidate:", task$key, method))
  upper <- frame$lambda_index == 1L ||
    abs(frame$lambda_relative_to_reference - 64) <= 1e-10
  data.frame(task_id = task$task_id, scenario = task$scenario,
             replication = task$replication, seed = task$seed, method = method,
             alpha = frame$alpha,
             d = if (method == "Logistic Group Elastic Net (adelie)") NA_real_ else frame$d,
             lambda_reference = frame$lambda_reference,
             validation_log_loss_v5 = frame$validation_log_loss,
             source_shard_file = task$shard_file,
             lambda_index = frame$lambda_index,
             lambda_relative_to_reference = frame$lambda_relative_to_reference,
             upper_boundary = upper, stringsAsFactors = FALSE)
}

lsg_tail_read_source_v6 <- function(root) {
  constants <- lsg_tail_constants_v6()
  version <- constants$source_version
  output <- file.path(root, "outputs", "study", version)
  metadata_path <- file.path(output, paste0("study_specification_", version, ".rds"))
  audit_path <- file.path(output, "final", paste0("LAMBDA_BOUNDARY_AUDIT_", version, ".rds"))
  lsg_tail_assert_v6(all(file.exists(c(metadata_path, audit_path))),
                     "The frozen V5 metadata/boundary audit is missing.")
  metadata <- readRDS(metadata_path)
  lsg_tail_assert_v6(identical(metadata$schema_version, "logistic_sglasso_truba_shards_v1") &&
                     identical(metadata$version, version) && identical(metadata$stage, "production"),
                     "Unexpected source schema/stage/version.")
  recomputed <- lsg_tail_v5_signatures_v6(metadata)
  lsg_tail_assert_v6(identical(recomputed$scientific, constants$source_scientific_signature) &&
                     identical(recomputed$tuning, constants$source_tuning_signature) &&
                     identical(metadata$scientific_signature, recomputed$scientific) &&
                     identical(metadata$tuning_signature, recomputed$tuning),
                     "The frozen V5 scientific/tuning signatures failed reconstruction.")
  source_relative <- lsg_tail_relative_v6(metadata$source_files)
  lsg_tail_assert_v6(length(source_relative) == 30L && !anyDuplicated(source_relative) &&
                     identical(names(metadata$source_sha256), metadata$source_files),
                     "Unexpected V5 source hash inventory.")
  source_manifest <- lsg_tail_manifest_v6(root, source_relative)
  original_sha <- unname(metadata$source_sha256[match(source_manifest$file, source_relative)])
  lsg_tail_assert_v6(identical(source_manifest$sha256, original_sha),
                     "A frozen V5 source file differs from the executed TRUBA code.")
  tasks <- metadata$task_grid
  lsg_tail_assert_v6(nrow(tasks) == 200L && nrow(metadata$design) == 4L &&
                     identical(tasks$task_id, seq_len(200L)) && !anyDuplicated(tasks$key) &&
                     !anyDuplicated(tasks$shard_file) &&
                     all(basename(tasks$shard_file) == tasks$shard_file) &&
                     all(tasks$scenario == metadata$design$scenario[tasks$scenario_index]) &&
                     identical(tasks$key, paste(tasks$scenario, tasks$replication, sep = "::")) &&
                     all(tasks$seed == metadata$configuration$seed_base +
                           metadata$design$scenario_index[tasks$scenario_index] * 100000L + tasks$replication) &&
                     all(table(tasks$scenario) == 50L) &&
                     all(vapply(split(tasks$replication, tasks$scenario), function(reps) {
                       !anyDuplicated(reps) && setequal(reps, seq_len(50L))
                     }, logical(1))) &&
                     identical(tasks$shard_file, sprintf("shard_%04d_%s_r%03d.rds", tasks$task_id,
                                gsub("[^A-Za-z0-9_.-]", "_", tasks$scenario), tasks$replication)),
                     "The frozen 200-task seed grid is invalid.")
  shard_paths <- file.path(output, "shards", tasks$shard_file)
  listed <- list.files(file.path(output, "shards"), pattern = "\\.rds$", full.names = FALSE)
  lsg_tail_assert_v6(setequal(listed, tasks$shard_file) && all(file.exists(shard_paths)),
                     "V5 requires exactly 200 expected source shards.")
  input_manifest <- lsg_tail_manifest_v6(root, c(metadata_path, audit_path, shard_paths))
  payload_names <- c("results", "sglasso_tuning", "external_tuning", "targets",
                     "diagnostics", "external_diagnostics", "group_diagnostics")
  selections <- vector("list", 200L)
  for (index in seq_len(200L)) {
    shard <- readRDS(shard_paths[index])
    task <- tasks[index, , drop = FALSE]
    valid <- identical(shard$schema_version, "logistic_sglasso_truba_shard_v1") &&
      identical(shard$scientific_signature, metadata$scientific_signature) &&
      identical(shard$tuning_signature, metadata$tuning_signature) &&
      identical(shard$key, task$key[[1L]]) && identical(shard$task, task) &&
      is.list(shard$payload) && all(vapply(payload_names, function(name) {
        is.data.frame(shard$payload[[name]]) && nrow(shard$payload[[name]]) > 0L
      }, logical(1)))
    lsg_tail_assert_v6(valid, paste("Invalid source shard:", task$key))
    result <- shard$payload$results
    lsg_tail_assert_v6(nrow(result) == 6L && !anyDuplicated(result$method) &&
                       setequal(result$method, metadata$configuration$methods) &&
                       all(result$scenario == task$scenario) &&
                       all(result$replication == task$replication) && all(result$seed == task$seed),
                       paste("Source task/method mismatch:", task$key))
    for (name in c("sglasso_tuning", "external_tuning", "targets", "diagnostics")) {
      frame <- shard$payload[[name]]
      lsg_tail_assert_v6(all(frame$scenario == task$scenario) &&
                         all(frame$replication == task$replication) && all(frame$seed == task$seed),
                         paste("Source payload task/seed mismatch:", task$key, name))
    }
    selections[[index]] <- rbind(
      lsg_tail_selection_v6(shard, task, constants$methods[1L], "selected_free_d"),
      lsg_tail_selection_v6(shard, task, constants$methods[2L], "selected_d0_boundary"),
      lsg_tail_selection_v6(shard, task, constants$methods[3L], "selected"))
  }
  selections <- do.call(rbind, selections)
  rownames(selections) <- NULL
  audit <- readRDS(audit_path)
  boundary_check_names <- c(
    "all_tasks_present", "frozen_augmented_grid_exact", "one_free_sglasso_selection_per_task",
    "one_d0_sglasso_selection_per_task", "one_adelie_selection_per_task",
    "selected_candidates_numerically_eligible", "free_sglasso_not_upper_lambda_boundary",
    "d0_sglasso_not_upper_lambda_boundary", "adelie_not_upper_lambda_boundary",
    "sglasso_not_lower_lambda_boundary")
  lsg_tail_assert_v6(identical(audit$schema_version, "logistic_sglasso_lambda_boundary_gate_v5") &&
                     identical(audit$version, version) && identical(audit$stage, "production") &&
                     identical(audit$scientific_signature, metadata$scientific_signature) &&
                     identical(audit$tuning_signature, metadata$tuning_signature) &&
                     identical(audit$accepted, FALSE) && is.data.frame(audit$selections) &&
                     nrow(audit$selections) == 600L,
                     "V5 boundary audit is not the frozen rejected production audit.")
  lsg_tail_assert_v6(is.data.frame(audit$checks) &&
                     !anyDuplicated(audit$checks$check) &&
                     setequal(audit$checks$check, boundary_check_names) &&
                     identical(audit$checks$passed[match(boundary_check_names, audit$checks$check)],
                               c(rep(TRUE, 6L), rep(FALSE, 3L), TRUE)),
                     "V5 rejection is not restricted to the three known upper-boundary checks.")
  keys <- function(frame) paste(frame$scenario, frame$replication, frame$method, sep = "::")
  lsg_tail_assert_v6(!anyDuplicated(keys(audit$selections)) &&
                     setequal(keys(selections), keys(audit$selections)),
                     "Boundary audit selection keys differ from source shards.")
  audit_selection <- audit$selections[match(keys(selections), keys(audit$selections)), ]
  lsg_tail_assert_v6(lsg_tail_equal_v6(selections$alpha, audit_selection$selected_alpha) &&
                     lsg_tail_equal_v6(selections$d, audit_selection$selected_d) &&
                     lsg_tail_equal_v6(selections$validation_log_loss_v5,
                                       audit_selection$validation_log_loss) &&
                     lsg_tail_equal_v6(selections$lambda_index, audit_selection$lambda_index) &&
                     lsg_tail_equal_v6(selections$lambda_relative_to_reference,
                                       audit_selection$lambda_relative_to_reference) &&
                     identical(selections$upper_boundary, audit_selection$upper_boundary),
                     "Boundary audit values differ from the underlying V5 candidates.")
  boundary <- selections[selections$upper_boundary, , drop = FALSE]
  boundary$case_id <- seq_len(nrow(boundary))
  boundary <- boundary[, constants$case_columns]
  rownames(boundary) <- NULL
  lsg_tail_assert_v6(nrow(boundary) == 16L && length(unique(boundary$task_id)) == 11L &&
                     identical(as.integer(table(factor(boundary$method, constants$methods))), c(8L, 4L, 4L)),
                     "The frozen diagnostic must contain 16 selections in 11 tasks.")
  case_path <- file.path(root, "config", "logistic_sglasso_lambda_tail_cases_v6.csv")
  lsg_tail_assert_v6(file.exists(case_path), "The frozen V6 case inventory is missing.")
  cases <- utils::read.csv(case_path, stringsAsFactors = FALSE, check.names = FALSE)
  lsg_tail_assert_v6(identical(names(cases), constants$case_columns) &&
                     lsg_tail_equal_v6(cases, boundary, tolerance = 1e-12),
                     "The frozen case CSV does not equal all V5 upper-boundary selections.")
  lsg_tail_verify_manifest_v6(root, input_manifest)
  list(metadata = metadata, cases = cases,
       task_grid = tasks[tasks$task_id %in% cases$task_id, , drop = FALSE],
       input_manifest = input_manifest, source_manifest = source_manifest,
       source_output_dir = output, selection_count = nrow(selections))
}

lsg_tail_required_source_paths_v6 <- function(root) {
  version <- lsg_tail_constants_v6()$source_version
  metadata <- readRDS(file.path(root, "outputs", "study", version,
                                paste0("study_specification_", version, ".rds")))
  fresh <- c(
    "R/logistic_sglasso_lambda_tail_v6.R",
    "R/logistic_sglasso_lambda_tail_io_v6.R",
    "scripts/36_validate_logistic_sglasso_lambda_tail_v6.R",
    "scripts/37_run_logistic_sglasso_lambda_tail_v6.R",
    "config/logistic_sglasso_lambda_tail_cases_v6.csv",
    "LOGISTIC_SGLASSO_LAMBDA_TAIL_PROTOCOL_V6.md",
    "truba/run_logistic_sglasso_tail_smoke_v6.slurm",
    "truba/run_logistic_sglasso_tail_diagnostic_v6.slurm",
    "truba/README_TRUBA_LOGISTIC_SGLASSO_V6.md")
  sort(c(lsg_tail_relative_v6(metadata$source_files),
         file.path("logistic_prework", fresh)))
}

lsg_tail_verify_manifest_v6 <- function(root, manifest, expected_relative = NULL) {
  lsg_tail_assert_v6(is.data.frame(manifest) &&
                     identical(names(manifest), c("file", "bytes", "sha256")) &&
                     nrow(manifest) > 0L && !anyDuplicated(manifest$file) &&
                     all(is.finite(manifest$bytes)) && all(manifest$bytes >= 0) &&
                     all(grepl("^[0-9a-f]{64}$", manifest$sha256)),
                     "Invalid SHA-256/byte manifest schema.")
  if (!is.null(expected_relative)) {
    lsg_tail_assert_v6(setequal(manifest$file, lsg_tail_relative_v6(expected_relative)),
                       "Manifest file inventory differs from the expected artifacts.")
  }
  observed <- lsg_tail_manifest_v6(root, manifest$file)
  manifest <- manifest[order(manifest$file), , drop = FALSE]
  rownames(manifest) <- rownames(observed) <- NULL
  lsg_tail_assert_v6(lsg_tail_equal_v6(manifest, observed),
                     "Manifest SHA-256 or byte-size validation failed.")
  invisible(TRUE)
}

lsg_tail_smoke_checks_v6 <- function() {
  c("closed_form_penalty_limit", "boundary_alpha_d_zero_target",
    "unequal_group_weights_cancel", "training_intercept_profiled",
    "coefficient_prediction_reconstruction", "large_finite_lambda_approaches_limit",
    "tail_paths_and_anchor_schema", "firth_training_only_no_fallback",
    "test_fields_rejected", "no_clobber_and_manifest_checks",
    "frozen_v5_sources_unchanged")
}

lsg_tail_payload_checks_v6 <- function() {
  c("fixed_cases_preserved", "five_points_per_case", "finite_points_eligible",
    "limit_points_eligible", "lambda_reference_reproduced", "v5_anchor_reproduced",
    "target_reproduced", "training_validation_only", "coefficient_reconstruction")
}

lsg_tail_check_frame_v6 <- function(checks, expected) {
  is.data.frame(checks) && identical(names(checks), c("check", "passed")) &&
    !anyDuplicated(checks$check) && setequal(checks$check, expected) &&
    is.logical(checks$passed) && isTRUE(all(checks$passed))
}

lsg_tail_verify_smoke_v6 <- function(root, code_manifest, source_configuration,
                                    required = FALSE, strict_runtime = FALSE) {
  version <- lsg_tail_constants_v6()$smoke_version
  directory <- file.path(root, "outputs", "study", version)
  paths <- file.path(directory, c(paste0("SMOKE_ACCEPTED_", version, ".rds"),
                                  paste0("smoke_checks_", version, ".csv"),
                                  paste0("smoke_output_manifest_", version, ".csv")))
  if (!all(file.exists(paths))) {
    lsg_tail_assert_v6(!any(file.exists(paths)) && !required,
                       "The required V6 smoke receipt/manifest is missing or incomplete.")
    return(list(present = FALSE, runtime_matches_source = FALSE,
                runtime_matches_current = FALSE))
  }
  manifest <- utils::read.csv(paths[3L], stringsAsFactors = FALSE)
  lsg_tail_verify_manifest_v6(root, manifest, paths[1:2])
  receipt <- readRDS(paths[1L])
  checks <- utils::read.csv(paths[2L], stringsAsFactors = FALSE)
  expected_sha <- stats::setNames(code_manifest$sha256, code_manifest$file)
  recorded_sha <- receipt$source_sha256
  lsg_tail_assert_v6(identical(receipt$schema_version, "logistic_sglasso_lambda_tail_smoke_v6") &&
                     identical(receipt$version, version) && isTRUE(receipt$passed) &&
                     lsg_tail_check_frame_v6(receipt$checks, lsg_tail_smoke_checks_v6()) &&
                     lsg_tail_equal_v6(receipt$checks, checks) &&
                     !anyDuplicated(names(recorded_sha)) &&
                     setequal(names(recorded_sha), names(expected_sha)) &&
                     identical(unname(recorded_sha[names(expected_sha)]), unname(expected_sha)),
                     "V6 smoke receipt checks or frozen code hashes are invalid.")
  runtime <- receipt[c("r_version", "platform", "package_versions", "rng_kind")]
  lsg_tail_assert_v6(length(runtime$r_version) == 1L && nzchar(runtime$r_version) &&
                     length(runtime$platform) == 1L && nzchar(runtime$platform) &&
                     identical(names(runtime$package_versions), lsg_tail_constants_v6()$packages) &&
                     !anyNA(runtime$package_versions) && all(nzchar(runtime$package_versions)) &&
                     identical(runtime$rng_kind, lsg_tail_constants_v6()$rng_kind),
                     "V6 smoke receipt runtime inventory is incomplete.")
  matches_source <- lsg_tail_runtime_matches_v6(runtime, source_configuration)
  matches_current <- identical(runtime, lsg_tail_runtime_v6())
  lsg_tail_assert_v6(!strict_runtime || (matches_source && matches_current),
                     "Diagnostic requires a V6 smoke receipt from the exact V5/current runtime.")
  list(present = TRUE, runtime_matches_source = matches_source,
       runtime_matches_current = matches_current,
       manifest = lsg_tail_manifest_v6(root, paths), receipt = receipt)
}

lsg_tail_validate_payload_v6 <- function(payload, cases, configuration) {
  lsg_tail_assert_v6(is.list(payload) &&
                     setequal(names(payload), c("curves", "anchors", "target_diagnostics",
                                                "training_validation_fingerprint", "checks", "test_used")) &&
                     identical(payload$test_used, FALSE) &&
                     is.character(payload$training_validation_fingerprint) &&
                     length(payload$training_validation_fingerprint) == 1L &&
                     grepl("^[0-9a-f]{64}$", payload$training_validation_fingerprint) &&
                     lsg_tail_check_frame_v6(payload$checks, lsg_tail_payload_checks_v6()),
                     "Diagnostic payload schema, test exclusion or named checks failed.")
  curves <- payload$curves
  anchors <- payload$anchors
  targets <- payload$target_diagnostics
  curve_columns <- c("case_id", "scenario", "replication", "method", "alpha", "d",
                     "lambda_relative", "lambda", "point_type", "validation_log_loss",
                     "training_log_loss", "selected_groups", "finite_kkt", "penalty_kkt",
                     "intercept_score", "reconstruction_error", "numerically_eligible")
  anchor_columns <- c("case_id", "lambda_reference_stored", "lambda_reference_recomputed",
                      "lambda_reference_error", "validation_log_loss_v5",
                      "validation_log_loss_replay", "anchor_loss_error", "anchor_within_tolerance")
  lsg_tail_assert_v6(is.data.frame(curves) && identical(names(curves), curve_columns) &&
                     is.data.frame(anchors) && identical(names(anchors), anchor_columns) &&
                     nrow(curves) == 5L * nrow(cases) && nrow(anchors) == nrow(cases) &&
                     !anyDuplicated(anchors$case_id) && setequal(anchors$case_id, cases$case_id) &&
                     setequal(unique(curves$case_id), cases$case_id),
                     "Diagnostic curve/anchor dimensions or fields differ from the frozen cases.")
  expected_groups <- configuration$expected_group_count
  lsg_tail_assert_v6(length(expected_groups) == 1L && is.finite(expected_groups) && expected_groups > 0,
                     "Expected training Firth group count is missing.")
  for (i in seq_len(nrow(cases))) {
    case <- cases[i, , drop = FALSE]
    points <- curves[curves$case_id == case$case_id, , drop = FALSE]
    anchor <- anchors[anchors$case_id == case$case_id, , drop = FALSE]
    lsg_tail_assert_v6(nrow(points) == 5L && !anyDuplicated(points$lambda_relative) &&
                       setequal(points$lambda_relative, c(configuration$ratios, Inf)) &&
                       all(points$scenario == case$scenario) &&
                       all(points$replication == case$replication) && all(points$method == case$method) &&
                       all(points$alpha == case$alpha) &&
                       (if (is.na(case$d)) all(is.na(points$d)) else all(points$d == case$d)),
                       "A diagnostic case changed alpha/d, identity or its five lambda points.")
    finite <- is.finite(points$lambda_relative)
    limit <- !finite
    expected_type <- ifelse(limit, "penalty_limit",
                            ifelse(points$lambda_relative == 64, "v5_anchor_replay", "finite_extension"))
    lsg_tail_assert_v6(identical(points$point_type, expected_type) &&
                       all(is.finite(points$lambda[finite])) && all(points$lambda[finite] > 0) &&
                       identical(points$lambda[limit], Inf) &&
                       max(abs(points$lambda[finite] -
                                 points$lambda_relative[finite] * anchor$lambda_reference_recomputed)) <=
                         configuration$lambda_reference_tolerance * max(1, abs(points$lambda[finite])),
                       "Diagnostic lambda/reference or limit labels are inconsistent.")
    lsg_tail_assert_v6(all(is.finite(points$validation_log_loss)) &&
                       all(is.finite(points$training_log_loss)) &&
                       all(points$validation_log_loss >= 0) && all(points$training_log_loss >= 0) &&
                       is.logical(points$numerically_eligible) && all(points$numerically_eligible) &&
                       all(is.finite(points$selected_groups)) &&
                       all(points$selected_groups == as.integer(points$selected_groups)) &&
                       all(points$selected_groups >= 0 & points$selected_groups <= expected_groups) &&
                       all(is.finite(points$intercept_score)) && all(points$intercept_score >= 0) &&
                       all(is.finite(points$reconstruction_error)) && all(points$reconstruction_error >= 0) &&
                       all(points$reconstruction_error <= configuration$reconstruction_tolerance) &&
                       all(is.finite(points$finite_kkt[finite])) && all(points$finite_kkt[finite] >= 0) &&
                       all(is.na(points$finite_kkt[limit])) && all(is.na(points$penalty_kkt[finite])) &&
                       is.finite(points$penalty_kkt[limit]) && points$penalty_kkt[limit] >= 0 &&
                       points$penalty_kkt[limit] <= configuration$endpoint_tolerance &&
                       points$intercept_score[limit] <= configuration$endpoint_tolerance,
                       "Diagnostic losses, reconstruction, eligibility or analytic-limit KKT failed.")
    if (case$method != "Logistic Group Elastic Net (adelie)") {
      lsg_tail_assert_v6(all(points$finite_kkt[finite] <= configuration$kkt_limit),
                         "SGLASSO finite-path KKT exceeds the frozen tolerance.")
    }
    reference_error <- abs(anchor$lambda_reference_recomputed - case$lambda_reference) /
      max(1, abs(case$lambda_reference))
    replay <- points$validation_log_loss[points$lambda_relative == 64]
    loss_error <- abs(replay - case$validation_log_loss_v5)
    numeric_anchor_fields <- c("lambda_reference_stored", "lambda_reference_recomputed",
                               "lambda_reference_error", "validation_log_loss_v5",
                               "validation_log_loss_replay", "anchor_loss_error")
    lsg_tail_assert_v6(all(is.finite(unlist(anchor[1L, numeric_anchor_fields], use.names = FALSE))) &&
                       anchor$lambda_reference_recomputed > 0 &&
                       lsg_tail_equal_v6(anchor$lambda_reference_stored, case$lambda_reference, 1e-12) &&
                       lsg_tail_equal_v6(anchor$validation_log_loss_v5, case$validation_log_loss_v5, 1e-12) &&
                       lsg_tail_equal_v6(anchor$validation_log_loss_replay, replay, 1e-12) &&
                       lsg_tail_equal_v6(anchor$lambda_reference_error, reference_error, 1e-12) &&
                       lsg_tail_equal_v6(anchor$anchor_loss_error, loss_error, 1e-12) &&
                       reference_error <= configuration$lambda_reference_tolerance &&
                       loss_error <= configuration$anchor_loss_tolerance &&
                       isTRUE(anchor$anchor_within_tolerance),
                       "V5 lambda reference or 64-anchor replay did not reproduce within tolerance.")
  }
  target_columns <- c("group", "success", "converged", "finite_coefficients",
                      "target_norm_stored", "target_norm_recomputed", "target_norm_relative_error")
  lsg_tail_assert_v6(is.data.frame(targets) && all(target_columns %in% names(targets)) &&
                     nrow(targets) == expected_groups && !anyDuplicated(targets$group) &&
                     setequal(as.character(targets$group), as.character(seq_len(expected_groups))) &&
                     all(targets$success) && all(targets$converged) && all(targets$finite_coefficients) &&
                     all(is.finite(targets$target_norm_stored)) &&
                     all(is.finite(targets$target_norm_recomputed)) &&
                     all(is.finite(targets$target_norm_relative_error)),
                     "The training-only Firth target diagnostics are incomplete or failed.")
  expected_error <- abs(targets$target_norm_stored - targets$target_norm_recomputed) /
    pmax(1, abs(targets$target_norm_stored))
  lsg_tail_assert_v6(lsg_tail_equal_v6(targets$target_norm_relative_error, expected_error, 1e-12) &&
                     all(expected_error <= configuration$target_norm_tolerance),
                     "The recomputed Firth target norm does not reproduce V5.")
  invisible(TRUE)
}

lsg_tail_object_hash_v6 <- function(value) {
  digest::digest(value, algo = "sha256", serialize = TRUE, serializeVersion = 2L)
}

lsg_tail_specification_v6 <- function(source, code_manifest, smoke, runtime, version) {
  configuration <- utils::modifyList(source$metadata$configuration, lsg_tail_defaults_v6())
  configuration$version <- version
  configuration$stage <- "lambda_tail_diagnostic"
  configuration$source_configuration <- source$metadata$configuration
  configuration$methods <- lsg_tail_constants_v6()$methods
  configuration$external_methods <- "Logistic Group Elastic Net (adelie)"
  configuration$primary_endpoint <- "conditional_validation_log_loss_tail"
  configuration$evaluation_sample <- "training_and_validation_only_no_test"
  configuration$tuning_sample <- "none_alpha_d_frozen_from_V5"
  configuration$nlambda <- length(configuration$ratios)
  configuration$lambda_relative_grid <- configuration$ratios
  configuration$lambda_upper_multiplier <- max(configuration$ratios)
  configuration$include_analytic_limit <- TRUE
  configuration$alpha_grid <- configuration$d_grid <- NULL
  configuration$replications_per_scenario <- NULL
  configuration$benchmark_alpha_grid <- NULL
  configuration$lambda_base_nlambda <- configuration$lambda_extension_multipliers <- NULL
  configuration$lambda_original_upper_multiplier <- configuration$lambda_min_ratio <- NULL
  configuration$lambda_grid_policy <- "fixed_512_256_128_64_plus_analytic_penalty_limit"
  configuration$sglasso_selection_policy <- "none_fixed_cases_conditional_curve_only"
  configuration$nested_d0_policy <- "fixed_V5_selected_d0_case_only"
  configuration$upper_lambda_selection_policy <- "record_all_tail_points_no_interior_or_improvement_gate"
  configuration$expected_group_count <- unique(source$metadata$design$groups)
  configuration$interpretation <- "conditional_fixed_alpha_d_validation_tail_diagnostic_only"
  configuration$test_used <- FALSE
  configuration$source_rng_assumption <- "V5_Rscript_vanilla_default_not_recorded_in_V5_metadata"
  tasks <- source$task_grid
  tasks$diagnostic_shard_file <- sprintf("tail_task_%04d.rds", tasks$task_id)
  rownames(tasks) <- NULL
  value <- list(
    schema_version = "logistic_sglasso_lambda_tail_diagnostic_v6", version = version,
    configuration = configuration, cases = source$cases, task_grid = tasks,
    source_scientific_signature = source$metadata$scientific_signature,
    source_tuning_signature = source$metadata$tuning_signature,
    input_manifest = source$input_manifest, code_manifest = code_manifest,
    smoke_manifest = smoke$manifest, runtime = runtime)
  value$scientific_signature <- lsg_tail_object_hash_v6(value)
  value
}

lsg_tail_specification_valid_v6 <- function(value) {
  is.list(value) && identical(value$schema_version, "logistic_sglasso_lambda_tail_diagnostic_v6") &&
    identical(value$scientific_signature,
              lsg_tail_object_hash_v6(value[setdiff(names(value), c("scientific_signature", "initial_execution"))]))
}

lsg_tail_output_paths_v6 <- function(root, version) {
  directory <- file.path(root, "outputs", "study", version)
  list(directory = directory, shard_directory = file.path(directory, "shards"),
       metadata = file.path(directory, paste0("study_specification_", version, ".rds")),
       completion = file.path(directory, paste0("DIAGNOSTIC_COMPLETED_", version, ".txt")),
       output_manifest = file.path(directory, paste0("output_manifest_", version, ".csv")))
}

lsg_tail_initialize_v6 <- function(root, specification, requested_cores) {
  paths <- lsg_tail_output_paths_v6(root, specification$version)
  if (file.exists(paths$metadata)) {
    previous <- readRDS(paths$metadata)
    lsg_tail_assert_v6(lsg_tail_specification_valid_v6(previous) &&
                       identical(previous$scientific_signature, specification$scientific_signature),
                       "Existing V6 metadata has a different source/code/runtime/case signature; use a new version.")
    specification <- previous
  } else {
    lsg_tail_assert_v6(!dir.exists(paths$directory) ||
                       !length(list.files(paths$directory, all.files = TRUE, no.. = TRUE)),
                       "A non-empty V6 output directory lacks compatible immutable metadata.")
    dir.create(paths$shard_directory, recursive = TRUE, showWarnings = FALSE)
    specification$initial_execution <- list(
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      requested_cores = requested_cores, slurm_job = Sys.getenv("SLURM_JOB_ID", unset = "not-set"),
      host = Sys.info()[["nodename"]], library_paths = .libPaths())
    lsg_tail_atomic_rds_v6(specification, paths$metadata)
  }
  dir.create(paths$shard_directory, recursive = TRUE, showWarnings = FALSE)
  tables <- list(input_manifest = specification$input_manifest,
                 source_manifest = specification$code_manifest,
                 smoke_manifest = specification$smoke_manifest,
                 case_inventory = specification$cases)
  for (name in names(tables)) {
    lsg_tail_atomic_v6(tables[[name]], file.path(paths$directory,
                       paste0(name, "_", specification$version, ".csv")), "csv", TRUE)
  }
  specification
}

lsg_tail_shard_v6 <- function(task, payload, specification, runtime_seconds) {
  lsg_tail_validate_payload_v6(payload, specification$cases[
    specification$cases$task_id == task$task_id, , drop = FALSE], specification$configuration)
  value <- list(schema_version = "logistic_sglasso_lambda_tail_diagnostic_shard_v6",
                version = specification$version, scientific_signature = specification$scientific_signature,
                key = task$key[[1L]], task = task, runtime_seconds = runtime_seconds,
                created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), payload = payload)
  value$content_sha256 <- lsg_tail_object_hash_v6(value)
  value
}

lsg_tail_read_shard_v6 <- function(path, task, specification) {
  if (!file.exists(path)) return(NULL)
  value <- tryCatch(readRDS(path), error = function(condition) {
    stop("Unreadable existing V6 shard; not overwriting: ", path, call. = FALSE)
  })
  lsg_tail_assert_v6(is.list(value) &&
                     identical(value$schema_version, "logistic_sglasso_lambda_tail_diagnostic_shard_v6") &&
                     identical(value$version, specification$version) &&
                     identical(value$scientific_signature, specification$scientific_signature) &&
                     identical(value$key, task$key[[1L]]) && identical(value$task, task) &&
                     is.finite(value$runtime_seconds) && value$runtime_seconds >= 0 &&
                     identical(value$content_sha256,
                               lsg_tail_object_hash_v6(value[setdiff(names(value), "content_sha256")])),
                     paste("Conflicting or invalid V6 shard; not overwriting:", path))
  lsg_tail_validate_payload_v6(value$payload, specification$cases[
    specification$cases$task_id == task$task_id, , drop = FALSE], specification$configuration)
  value
}

lsg_tail_load_shards_v6 <- function(root, specification) {
  paths <- lsg_tail_output_paths_v6(root, specification$version)
  tasks <- specification$task_grid
  listed <- list.files(paths$shard_directory, pattern = "\\.rds$", full.names = FALSE)
  lsg_tail_assert_v6(all(listed %in% tasks$diagnostic_shard_file),
                     "Unexpected diagnostic shards exist; refusing to silently ignore them.")
  lapply(seq_len(nrow(tasks)), function(i) {
    lsg_tail_read_shard_v6(file.path(paths$shard_directory, tasks$diagnostic_shard_file[i]),
                          tasks[i, , drop = FALSE], specification)
  })
}

lsg_tail_blind_source_shard_v6 <- function(source_shard) {
  keep <- c("scenario", "replication", "seed", "target_mode", "target_available",
            "target_l2_norm", "failed_groups")
  target <- source_shard$payload$targets
  lsg_tail_assert_v6(all(keep %in% names(target)), "Missing source target diagnostics.")
  out <- source_shard[c("schema_version", "scientific_signature", "tuning_signature", "key", "task")]
  out$payload <- list(targets = target[, keep, drop = FALSE])
  out
}

lsg_tail_summary_v6 <- function(curves, cases) {
  rows <- lapply(seq_len(nrow(cases)), function(i) {
    case <- cases[i, , drop = FALSE]
    points <- curves[curves$case_id == case$case_id, , drop = FALSE]
    loss <- function(ratio) points$validation_log_loss[points$lambda_relative == ratio]
    out <- case[, c("case_id", "task_id", "scenario", "replication", "method", "alpha", "d")]
    out$interpretation <- "conditional_fixed_alpha_d_no_test_or_method_ranking"
    out$validation_loss_v5_64 <- case$validation_log_loss_v5
    out$validation_loss_replayed_64 <- loss(64)
    for (ratio in c(128, 256, 512, Inf)) {
      label <- if (is.infinite(ratio)) "limit" else as.character(ratio)
      out[[paste0("validation_loss_", label)]] <- loss(ratio)
      out[[paste0("delta_", label, "_minus_replayed_64")]] <- loss(ratio) - loss(64)
      out[[paste0("delta_", label, "_minus_v5_64")]] <- loss(ratio) - case$validation_log_loss_v5
    }
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

lsg_tail_final_tables_v6 <- function(specification, shards) {
  lsg_tail_assert_v6(length(shards) == 11L && all(vapply(shards, is.list, logical(1))),
                     "Diagnostic finalization requires all 11 valid task shards.")
  bind <- function(name) {
    out <- do.call(rbind, lapply(shards, function(shard) shard$payload[[name]]))
    rownames(out) <- NULL
    out
  }
  curves <- bind("curves")
  anchors <- bind("anchors")
  targets <- do.call(rbind, lapply(shards, function(shard) {
    frame <- shard$payload$target_diagnostics
    cbind(shard$task[rep(1L, nrow(frame)), c("task_id", "scenario", "replication", "seed")], frame)
  }))
  task_checks <- do.call(rbind, lapply(shards, function(shard) {
    cbind(task_id = shard$task$task_id, shard$payload$checks)
  }))
  rownames(targets) <- rownames(task_checks) <- NULL
  checks <- data.frame(
    check = c("all_11_signature_matching_task_shards", "all_16_fixed_cases", "all_80_curve_points",
              "all_16_anchor_replays", "all_2200_training_firth_group_audits",
              "all_99_named_task_checks", "test_not_used", "not_a_production_or_calibration_approval"),
    passed = c(length(shards) == 11L, setequal(unique(curves$case_id), specification$cases$case_id) &&
                 length(unique(curves$case_id)) == 16L,
               nrow(curves) == 80L && all(table(curves$case_id) == 5L),
               nrow(anchors) == 16L && !anyDuplicated(anchors$case_id),
               nrow(targets) == 2200L && all(table(targets$task_id) == 200L),
               nrow(task_checks) == 99L && all(task_checks$passed),
               all(vapply(shards, function(shard) identical(shard$payload$test_used, FALSE), logical(1))),
               identical(specification$configuration$interpretation,
                         "conditional_fixed_alpha_d_validation_tail_diagnostic_only")),
    stringsAsFactors = FALSE)
  lsg_tail_assert_v6(all(checks$passed), "A diagnostic completeness/dimension check failed.")
  list(tail_curves = curves, tail_anchors = anchors, tail_target_diagnostics = targets,
       tail_task_checks = task_checks, tail_summary = lsg_tail_summary_v6(curves, specification$cases),
       diagnostic_checks = checks)
}

lsg_tail_expected_artifacts_v6 <- function(root, specification) {
  version <- specification$version
  paths <- lsg_tail_output_paths_v6(root, version)
  tables <- c("input_manifest", "source_manifest", "smoke_manifest", "case_inventory",
              "tail_curves", "tail_anchors", "tail_target_diagnostics", "tail_task_checks",
              "tail_summary", "diagnostic_checks", "shard_manifest")
  c(paths$metadata, file.path(paths$directory, paste0(tables, "_", version, ".csv")),
    file.path(paths$shard_directory, specification$task_grid$diagnostic_shard_file),
    list.files(paths$directory, pattern = "^run_attempts_.*\\.csv$", full.names = TRUE))
}

lsg_tail_read_table_v6 <- function(path, expected) {
  # Automatic CSV conversion turns character group IDs into integers and an
  # all-empty character message column into logical NA. Preserve the schema
  # of the independently reconstructed shard table when verifying exports.
  classes <- vapply(expected, function(column) {
    if (is.character(column)) return("character")
    if (is.logical(column)) return("logical")
    if (is.integer(column)) return("integer")
    if (is.numeric(column)) return("numeric")
    stop("Unsupported diagnostic CSV column class.", call. = FALSE)
  }, character(1))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                  colClasses = classes)
}

lsg_tail_finalize_v6 <- function(root, specification, shards) {
  paths <- lsg_tail_output_paths_v6(root, specification$version)
  lsg_tail_verify_manifest_v6(root, specification$input_manifest)
  lsg_tail_verify_manifest_v6(root, specification$code_manifest)
  lsg_tail_verify_manifest_v6(root, specification$smoke_manifest)
  tables <- lsg_tail_final_tables_v6(specification, shards)
  for (name in names(tables)) {
    lsg_tail_atomic_v6(tables[[name]], file.path(paths$directory,
                       paste0(name, "_", specification$version, ".csv")), "csv", TRUE)
  }
  shard_paths <- file.path(paths$shard_directory, specification$task_grid$diagnostic_shard_file)
  shard_manifest <- lsg_tail_manifest_v6(root, shard_paths)
  lsg_tail_atomic_v6(shard_manifest, file.path(paths$directory,
                     paste0("shard_manifest_", specification$version, ".csv")), "csv", TRUE)
  output_manifest <- lsg_tail_manifest_v6(root, lsg_tail_expected_artifacts_v6(root, specification))
  lsg_tail_atomic_v6(output_manifest, paths$output_manifest, "csv", TRUE)
  lsg_tail_verify_manifest_v6(root, output_manifest)
  lsg_tail_atomic_v6(c(
    paste("version", specification$version), "stage conditional_lambda_tail_diagnostic",
    paste("scientific_signature", specification$scientific_signature),
    paste("source_scientific_signature", specification$source_scientific_signature),
    paste("source_tuning_signature", specification$source_tuning_signature),
    "completed_tasks 11", "fixed_cases 16", "curve_points 80", "test_used FALSE",
    "calibration_or_production_approved FALSE",
    paste("output_manifest_sha256", digest::digest(file = paths$output_manifest, algo = "sha256"))
  ), paths$completion, "text", TRUE)
  invisible(tables)
}

lsg_tail_verify_output_v6 <- function(root, specification) {
  paths <- lsg_tail_output_paths_v6(root, specification$version)
  lsg_tail_assert_v6(file.exists(paths$completion) && file.exists(paths$output_manifest),
                     "V6 diagnostic is incomplete: its final marker/manifest is missing.")
  stored <- readRDS(paths$metadata)
  lsg_tail_assert_v6(lsg_tail_specification_valid_v6(stored) &&
                     identical(stored$scientific_signature, specification$scientific_signature),
                     "V6 diagnostic metadata does not match its immutable inputs/code/runtime.")
  for (name in c("input_manifest", "code_manifest", "smoke_manifest")) {
    lsg_tail_verify_manifest_v6(root, stored[[name]])
  }
  manifest <- utils::read.csv(paths$output_manifest, stringsAsFactors = FALSE)
  expected <- lsg_tail_expected_artifacts_v6(root, stored)
  lsg_tail_verify_manifest_v6(root, manifest, expected)
  actual_files <- list.files(paths$directory, recursive = TRUE, full.names = TRUE)
  lsg_tail_assert_v6(setequal(lsg_tail_relative_v6(actual_files),
                             lsg_tail_relative_v6(c(expected, paths$output_manifest, paths$completion))),
                     "The completed diagnostic directory has an unexpected visible file inventory.")
  marker <- readLines(paths$completion, warn = FALSE)
  required <- c(paste("version", stored$version), "stage conditional_lambda_tail_diagnostic",
                paste("scientific_signature", stored$scientific_signature),
                paste("source_scientific_signature", stored$source_scientific_signature),
                paste("source_tuning_signature", stored$source_tuning_signature),
                "completed_tasks 11", "fixed_cases 16", "curve_points 80", "test_used FALSE",
                "calibration_or_production_approved FALSE",
                paste("output_manifest_sha256", digest::digest(file = paths$output_manifest, algo = "sha256")))
  lsg_tail_assert_v6(identical(marker, required), "V6 completion marker does not bind the verified manifest.")
  shards <- lsg_tail_load_shards_v6(root, stored)
  tables <- lsg_tail_final_tables_v6(stored, shards)
  for (name in names(tables)) {
    path <- file.path(paths$directory, paste0(name, "_", stored$version, ".csv"))
    saved <- lsg_tail_read_table_v6(path, tables[[name]])
    lsg_tail_assert_v6(lsg_tail_equal_v6(saved, tables[[name]], tolerance = 1e-12),
                       paste("Saved diagnostic table does not reproduce its task shards:", name))
  }
  list(complete = TRUE, tasks = length(shards), cases = nrow(stored$cases), points = nrow(tables$tail_curves),
       checks = tables$diagnostic_checks, summary = tables$tail_summary)
}
