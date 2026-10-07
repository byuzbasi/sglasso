# Same dispatcher for tiny local smoke and future frozen workloads. No pilot
# launch is exposed until the separate release/range decision is approved.
lsg_dispatch_v12 <- function(p, worker, accept, max_seconds, task_limit = Inf) {
  jobs <- list(); launched <- 0L; outcomes <- list()
  on.exit({
    # Only subprocesses created by this invocation are terminated on an R error.
    for (job in jobs) try(tools::pskill(job$pid), silent = TRUE)
    if (length(jobs)) try(parallel::mccollect(jobs, wait = FALSE), silent = TRUE)
  }, add = TRUE)
  repeat {
    now <- as.numeric(Sys.time())
    ready <- which(p$status == "pending")
    while (length(ready) && length(jobs) < p$cores && launched < task_limit &&
           now - p$compute_started < max_seconds) {
      i <- ready[1L]; ready <- ready[-1L]
      p$status[i] <- "running"; p$task_started[i] <- now; launched <- launched + 1L
      job <- parallel::mcparallel(local({ index <- i; worker(index) }), mc.set.seed = FALSE, silent = TRUE)
      attr(job, "task_index") <- i; jobs[[as.character(job$pid)]] <- job
      lsg_progress_emit_v12(p, "fitting", TRUE)
    }
    if (!length(jobs)) break
    done <- parallel::mccollect(jobs, wait = FALSE)
    if (length(done)) for (pid in names(done)) {
      i <- attr(jobs[[pid]], "task_index"); jobs[[pid]] <- NULL
      item <- done[[pid]]
      if (!is.list(item) || is.null(item$status)) item <- list(status = "error", message = as.character(item))
      item <- tryCatch(accept(i, item), error = function(e) list(status = "error", message = conditionMessage(e)))
      outcomes[[as.character(i)]] <- item
      if (identical(item$status, "passed")) {
        p$status[i] <- "completed"; p$new_completed <- c(p$new_completed, i)
        p$durations[i] <- as.numeric(Sys.time()) - p$task_started[i]
        p$latest <- p$spec$tasks$key[i]
      } else p$status[i] <- "failed"
      lsg_progress_emit_v12(p, "fitting", TRUE)
    }
    lsg_progress_emit_v12(p, "fitting")
    if (length(jobs)) Sys.sleep(min(0.2, p$interval))
  }
  outcomes
}

lsg_run_v12 <- function(root, version, cores = 1L, max_seconds = 120,
                       task_limit = Inf, interval = 30) {
  stage <- "smoke"
  lsg_assert_v7(.Platform$OS.type == "unix", "The validated dispatcher currently requires Unix forks.")
  lsg_assert_v7(cores %in% 1:2 && is.finite(max_seconds) && max_seconds > 0 && max_seconds <= 120,
    "Local integration is restricted to at most two workers and a 120-second soft budget.")
  allocated <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA_character_)))
  if (!is.na(allocated)) lsg_assert_v7(cores <= allocated, "Workers exceed allocated CPUs.")
  runtime <- lsg_runtime_v7() # dependency preflight: never installs packages
  cat("R:", runtime$r_version, "\nPlatform:", runtime$platform, "\n")
  print(runtime$package_versions)
  cat("Library paths:\n", paste(.libPaths(), collapse = "\n"), "\n")
  for (variable in c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")) {
    value <- Sys.getenv(variable)
    lsg_assert_v7(identical(value, "1"), paste("Preflight requires", variable, "=1 before R starts."))
  }
  output <- lsg_output_v7(root, version)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(file.path(output, paste0("output_manifest_", version, ".csv")))) {
    spec <- readRDS(file.path(output, "study_specification.rds"))
    lsg_validate_spec_v7(spec, root, stage, version, resume = TRUE)
    return(invisible(lsg_verify_v7(root, stage, version)))
  }
  lock <- lsg_path_v7(output, ".run_lock")
  lsg_assert_v7(dir.create(lock, showWarnings = FALSE), "Run lock exists; inspect its owner before retrying.")
  owner <- list(pid = Sys.getpid(), host = Sys.info()[["nodename"]], job = Sys.getenv("SLURM_JOB_ID"))
  on.exit(lsg_release_lock_v7(lock, owner), add = TRUE)
  lsg_atomic_v7(owner, output, ".run_lock/owner.rds")
  spec_path <- file.path(output, "study_specification.rds")
  spec <- if (file.exists(spec_path)) readRDS(spec_path) else {
    x <- lsg_make_spec_v7(root, stage, version)
    lsg_atomic_v7(x, output, "study_specification.rds"); x
  }
  lsg_validate_spec_v7(spec, root, stage, version, resume = TRUE)
  for (name in c("tasks", "source_manifest", "design"))
    lsg_atomic_v7(spec[[name]], output, paste0(if (name == "tasks") "task_grid" else name, ".csv"), "csv", TRUE)
  units <- rbind(lsg_design_checks_v7(lsg_design_for_stage_v7(root, "pilot")), lsg_metrics_unit_checks_v7())
  lsg_assert_v7(all(units$passed), "Design/metric unit checks failed before fitting.")
  lsg_atomic_v7(units, output, "unit_checks.rds", identical_ok = TRUE)
  old <- lsg_collect_v7(output, spec, deep = TRUE, require_all = FALSE)
  lsg_assert_v7(!nrow(old$checks) || all(old$checks$passed), "Existing gate-failed evidence is preserved; do not silently retry it.")
  p <- lsg_progress_new_v12(output, spec, which(old$present), cores, interval)
  closed <- FALSE
  on.exit(if (!closed) {
    p$status[p$status == "running"] <- "pending"
    try(lsg_progress_emit_v12(p, "interrupted", TRUE), silent = TRUE)
  }, add = TRUE)
  lsg_progress_emit_v12(p, "preflight", TRUE)
  if (any(!old$present) && task_limit > 0) {
    compile_lsg_core(quiet = TRUE); lsg_compile_hybrid_solver_v11(root)
  }
  p$compute_started <- as.numeric(Sys.time())
  worker <- function(i) {
    task <- spec$tasks[i, , drop = FALSE]
    scenario <- spec$design[spec$design$scenario_index == task$scenario_index, , drop = FALSE]
    started <- proc.time()[["elapsed"]]
    tryCatch({
      payload <- lsg_with_rng_v7(task$seed,
        lsg_run_task_v7(scenario, task$replication, task$seed, spec$configuration))
      payload$results$task_runtime_seconds <- proc.time()[["elapsed"]] - started
      checks <- lsg_payload_checks_v7(payload, task, scenario, spec$configuration)
      list(status = if (all(checks$passed)) "passed" else "gate_failed", payload = payload, checks = checks,
           message = paste(checks$check[!checks$passed], collapse = ";"))
    }, error = function(e) list(status = "error", message = conditionMessage(e)))
  }
  accept <- function(i, item) {
    task <- spec$tasks[i, , drop = FALSE]
    if (is.null(item$payload) && !identical(item$status, "error"))
      stop("A successful worker must supply a validated publishable payload.", call. = FALSE)
    if (!is.null(item$payload)) {
      shard <- list(schema_version = "prediction_selection_shard_v7",
        scientific_signature = spec$scientific_signature, task = task,
        payload = item$payload, checks = item$checks, attempt = p$id)
      relative <- file.path("shards", task$shard_file)
      path <- lsg_atomic_v7(shard, output, relative)
      lsg_atomic_v7(lsg_shard_receipt_v7(shard, path), output, paste0(relative, ".receipt.rds"))
      readback <- lsg_read_shard_v7(output, spec, task, deep = TRUE)
      item$status <- if (all(readback$checks$passed)) "passed" else "gate_failed"
    }
    item[c("status", "message")]
  }
  outcomes <- lsg_dispatch_v12(p, worker, accept, max_seconds, task_limit)
  record <- list(invocation_id = p$id, scientific_signature = spec$scientific_signature,
    runtime = lsg_runtime_v7(), library_paths = .libPaths(), session = capture.output(sessionInfo()),
    outcomes = outcomes, status = p$status, new_completed = p$new_completed,
    started_epoch = p$started, finished_epoch = as.numeric(Sys.time()),
    cores = cores, max_seconds_soft = max_seconds)
  lsg_atomic_v7(record, output, file.path("attempts", paste0(p$id, ".rds")))
  if (!all(p$status == "completed")) {
    lsg_progress_emit_v12(p, if (any(p$status == "failed")) "failed" else "paused", TRUE)
    closed <- TRUE
    return(invisible(list(complete = FALSE, outcomes = outcomes, output = output)))
  }
  lsg_progress_emit_v12(p, "aggregation_and_validation", TRUE)
  collected <- lsg_collect_v7(output, spec, deep = TRUE)
  lsg_finalize_v7(output, spec, collected)
  result <- lsg_verify_v7(root, stage, version)
  lsg_progress_emit_v12(p, "completed", TRUE); closed <- TRUE
  invisible(c(list(complete = TRUE), result))
}
