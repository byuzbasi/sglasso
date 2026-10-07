# GenAtHum median-only replay with the locally validated IRLS/APG core.
# Frozen V1/V3 scientific functions and package code are not modified.

gab_source_files_v4 <- function() unique(c(
  gab_source_files(),
  "src/logistic_sglasso_block_kernel_local_v3.hpp",
  "src/logistic_sglasso_hybrid_solver_local_apg_v4.cpp",
  "R/genathum_binary_irls_v4.R",
  "scripts/98_run_genathum_binary_irls_v4.R",
  "GENATHUM_BINARY_IRLS_PROTOCOL_V4.md"
))

gab_configuration_v4 <- function(e, root, stage) {
  configuration <- gab_configuration(e, root, stage)
  configuration$schema_version <- "genathum_binary_configuration_v4"
  configuration$outcomes <- configuration$outcomes[
    configuration$outcomes$role == "primary", , drop = FALSE
  ]
  rownames(configuration$outcomes) <- NULL
  gab_assert(nrow(configuration$outcomes) == 1L &&
    identical(configuration$outcomes$outcome_id, "morf4l1_median"),
    "Frozen GenAtHum median endpoint is missing.")
  # Exact IEEE-754 value from the accepted V3 specification. R 4.6's CSV
  # parser otherwise reads the printed decimal one ULP higher than R 4.3.
  bytes <- as.raw(c(0xde, 0xb3, 0xb8, 0xa1, 0xb9, 0x4d, 0x1f, 0x40))
  configuration$outcomes$cutpoint <- readBin(bytes, what = "double",
    n = 1L, size = 8L, endian = "little")
  gab_assert(is.finite(configuration$outcomes$cutpoint) &&
    abs(configuration$outcomes$cutpoint - 7.8259034413322) < 1e-12,
    "Frozen median threshold representation is invalid.")
  configuration$fit_configuration$solver_implementation <-
    "local_cached_block_irls_apg_v4"
  configuration
}

gab_clone_v4 <- function(fun, bindings) {
  scope <- list2env(bindings, parent = environment(fun))
  environment(fun) <- scope
  fun
}

gab_make_spec_v4 <- gab_clone_v4(gab_make_spec, list(
  gab_configuration = gab_configuration_v4,
  gab_source_files = gab_source_files_v4
))

gab_validate_spec_v4 <- gab_clone_v4(gab_validate_spec, list(
  gab_make_spec = gab_make_spec_v4
))

gab_load_environment_v4 <- function(root) {
  e <- gab_load_environment(root)
  gab_assert(requireNamespace("Rcpp", quietly = TRUE) &&
    requireNamespace("RcppArmadillo", quietly = TRUE),
    "Rcpp and RcppArmadillo must already be installed.")
  root <- normalizePath(root, mustWork = TRUE)
  paths <- file.path(root, c(
    "src/logistic_sglasso_hybrid_solver_local_apg_v4.cpp",
    "src/logistic_sglasso_block_kernel_local_v3.hpp"
  ))
  gab_assert(all(file.exists(paths)) && !any(dir.exists(paths)) &&
    !any(nzchar(Sys.readlink(paths))), "Missing regular IRLS core sources.")
  old_cppflags <- Sys.getenv("PKG_CPPFLAGS", unset = NA_character_)
  old_makevars <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  include <- paste0("-I", shQuote(file.path(root, "src")))
  Sys.setenv(PKG_CPPFLAGS = if (is.na(old_cppflags) || !nzchar(old_cppflags))
    include else paste(old_cppflags, include))
  local_makevars <- file.path(root, "config/Makevars.local")
  local_fortran <- paste0("/usr/local/gfortran/lib/gcc/",
    "aarch64-apple-darwin23/14.1.0/libemutls_w.a")
  if (file.exists(local_makevars) && file.exists(local_fortran))
    Sys.setenv(R_MAKEVARS_USER = local_makevars)
  on.exit({
    if (is.na(old_cppflags)) Sys.unsetenv("PKG_CPPFLAGS")
    else Sys.setenv(PKG_CPPFLAGS = old_cppflags)
    if (is.na(old_makevars)) Sys.unsetenv("R_MAKEVARS_USER")
    else Sys.setenv(R_MAKEVARS_USER = old_makevars)
  }, add = TRUE)
  Rcpp::sourceCpp(paths[1L], rebuild = FALSE, showOutput = FALSE,
    verbose = FALSE, env = e)
  required <- c("lsg_fit_one_hybrid_v11_cpp", "lsg_path_hybrid_v11_cpp")
  gab_assert(all(vapply(required, exists, logical(1), envir = e,
    inherits = FALSE, mode = "function")), "IRLS core exports are missing.")
  one <- e$lsg_fit_one_hybrid_v11_cpp
  path <- e$lsg_path_hybrid_v11_cpp
  e$lsg_fit_one_hybrid_v11_cpp <- function(...) one(...,
    reuse_apg_offsets = TRUE, use_irls = TRUE)
  e$lsg_path_hybrid_v11_cpp <- function(...) path(...,
    reuse_apg_offsets = TRUE, use_irls = TRUE)
  e
}

gab_validate_release_v4 <- function(e, project_root, root, version) {
  path <- file.path(root, "release", "genathum_binary_irls_v4",
    "LOCAL_VALIDATED_genathum_binary_irls_v4.rds")
  gab_assert(file.exists(path), "Missing accepted local V4 release receipt.")
  receipt <- readRDS(path)
  current <- gab_make_spec_v4(e, project_root, root, "production", version)
  gab_assert(identical(receipt$schema_version, "genathum_binary_irls_release_v4") &&
    isTRUE(receipt$accepted) && all(receipt$checks$passed) &&
    identical(receipt$production_version, version) &&
    identical(receipt$production_signature, current$scientific_signature) &&
    identical(receipt$source_manifest, current$source_manifest) &&
    identical(receipt$data_identity, current$data_identity) &&
    identical(receipt$configuration, current$configuration) &&
    identical(receipt$tasks, current$tasks) &&
    identical(receipt$runtime, current$runtime),
    "Local V4 release identity changed; production is blocked.")
  invisible(receipt)
}

gab_run_pending_v4 <- function(spec, data, e, out, cores, interval) {
  valid <- vapply(seq_len(nrow(spec$tasks)), function(i)
    gab_valid_shard(file.path(out, "shards", spec$tasks$shard_file[i]),
      spec$tasks[i, , drop = FALSE], spec$scientific_signature,
      spec$configuration), logical(1))
  pending <- which(!valid)
  status <- rep("pending", nrow(spec$tasks)); status[valid] <- "passed"
  started <- proc.time()[["elapsed"]]
  invocation <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "_", Sys.getpid())
  outcomes <- list(); active <- list(); next_task <- 1L
  completed_now <- 0L; last_update <- -Inf; failed <- FALSE
  launch <- function(i) {
    status[i] <<- "running"
    job <- parallel::mcparallel(gab_worker(i, spec, data, e, out), silent = TRUE)
    active[[as.character(job$pid)]] <<- list(job = job, index = i)
  }
  while (next_task <= length(pending) && length(active) < cores) {
    launch(pending[next_task]); next_task <- next_task + 1L
  }
  repeat {
    collected <- if (length(active)) parallel::mccollect(
      lapply(active, `[[`, "job"), wait = FALSE) else NULL
    if (length(collected)) for (pid in names(collected)) {
      info <- active[[pid]]; result <- collected[[pid]]
      outcomes[[as.character(info$index)]] <- result
      if (inherits(result, "try-error") || is.null(result$status) ||
          !identical(result$status, "passed")) {
        status[info$index] <- "failed"; failed <- TRUE
        cat("Task ", info$index, " failed: ",
          if (is.list(result)) result$message else as.character(result), "\n",
          sep = "")
      } else {
        status[info$index] <- "passed"
        completed_now <- completed_now + 1L
      }
      active[[pid]] <- NULL
    }
    while (!failed && next_task <= length(pending) && length(active) < cores) {
      launch(pending[next_task]); next_task <- next_task + 1L
    }
    elapsed <- proc.time()[["elapsed"]] - started
    if (elapsed - last_update >= interval || !length(active)) {
      phase <- if (failed && !length(active)) "failed" else if (!length(active))
        "fitting_complete" else "fitting"
      gab_write_progress(gab_progress_snapshot(spec, out, status, started,
        phase, invocation, completed_now), out)
      last_update <- elapsed
    }
    if (!length(active)) break
    Sys.sleep(1)
  }
  gab_atomic_rds(list(invocation_id = invocation, status = status,
    outcomes = outcomes, scientific_signature = spec$scientific_signature),
    file.path(out, "attempts", paste0(invocation, ".rds")))
  status
}

gab_run_v4 <- function(project_root, root, stage, version, cores = 1L,
                       interval = 30) {
  gab_assert(stage %in% c("smoke", "production") &&
    length(cores) == 1L && is.finite(cores) && cores >= 1L &&
    cores == as.integer(cores) && is.finite(interval) && interval >= 1,
    "Invalid stage, cores, or progress interval.")
  e <- gab_load_environment_v4(root)
  if (stage == "production")
    gab_validate_release_v4(e, project_root, root, version)
  proposed <- gab_make_spec_v4(e, project_root, root, stage, version)
  out <- gab_output_directory(root, version)
  old_shards <- file.path(out, "shards", proposed$tasks$shard_file)
  present <- which(file.exists(old_shards))
  gab_assert(all(vapply(present, function(i) gab_valid_shard(old_shards[i],
    proposed$tasks[i, , drop = FALSE], proposed$scientific_signature,
    proposed$configuration), logical(1))),
    "Existing invalid or signature-mismatched shard blocks resume; preserved for diagnosis.")
  dir.create(file.path(out, "shards"), recursive = TRUE, showWarnings = FALSE)
  spec_path <- file.path(out, "study_specification.rds")
  if (file.exists(spec_path)) {
    spec <- readRDS(spec_path)
    gab_validate_spec_v4(e, spec, project_root, root, stage, version)
  } else {
    spec <- proposed
    gab_atomic_rds(spec, spec_path)
  }
  cat("Stage:", stage, "\nVersion:", version,
    "\nScientific signature:", spec$scientific_signature,
    "\nWork units:", nrow(spec$tasks), "\nRequested workers:", cores, "\n")
  print(spec$runtime)
  data <- gab_load_data(project_root, spec$configuration)
  status <- gab_run_pending_v4(spec, data, e, out, as.integer(cores), interval)
  gab_assert(!any(status == "failed"),
    "One or more tasks failed; valid signature-matching shards are resumable.")
  marker_path <- file.path(out, paste0("COMPLETED_", version, ".rds"))
  if (!file.exists(marker_path)) gab_finalize(e, spec, out)
  gab_verify_outputs(spec, out)
  progress <- jsonlite::read_json(file.path(out, "progress.json"),
    simplifyVector = TRUE)
  progress$phase <- "complete"
  progress$eta_seconds <- 0
  progress$eta_status <- "complete_verified"
  progress$heartbeat_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  progress$estimated_completion_utc <- progress$heartbeat_utc
  gab_write_progress(progress, out)
  cat("Complete: TRUE\nOutput directory:", normalizePath(out),
    "\nCompleted tasks:", nrow(spec$tasks), "\n")
  invisible(list(spec = spec, output = out))
}
