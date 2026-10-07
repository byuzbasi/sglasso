# Only the dispatcher writes these mutable telemetry files. Scientific outputs
# use immutable V7 publication; telemetry does not certify acceptance.
lsg_progress_new_v12 <- function(output, spec, completed, cores, interval = 30, now = as.numeric(Sys.time())) {
  lsg_assert_v7(is.finite(interval) && interval > 0 && cores >= 1L, "Invalid progress controls.")
  p <- new.env(parent = emptyenv())
  p$output <- output; p$spec <- spec; p$cores <- cores; p$interval <- interval
  p$started <- now; p$compute_started <- now; p$last_emit <- -Inf
  p$id <- paste0(format(as.POSIXct(now, origin = "1970-01-01", tz = "UTC"), "%Y%m%dT%H%M%OS6"), "_", Sys.getpid())
  p$status <- rep("pending", nrow(spec$tasks)); p$status[completed] <- "completed"
  p$task_started <- rep(NA_real_, length(p$status)); p$durations <- p$task_started
  p$new_completed <- integer(); p$latest <- NA_character_
  previous <- lsg_path_v7(output, "progress.json")
  if (file.exists(previous)) {
    old <- jsonlite::fromJSON(previous)
    lsg_assert_v7(identical(old$run_id, spec$version) &&
      identical(old$scientific_signature, spec$scientific_signature), "Foreign progress file refused.")
    if (!is.null(old$last_completed_work_unit) && old$last_completed_work_unit %in% spec$tasks$key[completed])
      p$latest <- old$last_completed_work_unit
  }
  p
}

lsg_progress_eta_v12 <- function(p, now) {
  if (any(p$status == "failed")) return(NA_real_)
  todo <- which(p$status %in% c("running", "pending"))
  if (!length(todo)) return(0)
  if (length(p$new_completed) < 2L) return(NA_real_)
  scenario <- p$spec$tasks$scenario
  means <- tapply(p$durations[p$new_completed], scenario[p$new_completed], mean)
  cost <- as.numeric(means[scenario[todo]])
  if (any(!is.finite(cost))) return(NA_real_)
  running <- p$status[todo] == "running"
  running_left <- pmax(0, cost[running] - (now - p$task_started[todo[running]]))
  load <- c(running_left, rep(0, max(0L, p$cores - length(running_left))))
  for (seconds in cost[!running]) {
    k <- which.min(load); load[k] <- load[k] + seconds
  }
  max(load)
}

lsg_progress_emit_v12 <- function(p, phase, force = FALSE, now = as.numeric(Sys.time())) {
  if (!force && now - p$last_emit < p$interval) return(invisible(NULL))
  utc <- function(x) if (!is.finite(x)) NA_character_ else
    format(as.POSIXct(x, origin = "1970-01-01", tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")
  counts <- table(factor(p$status, levels = c("completed", "running", "failed", "pending")))
  total <- nrow(p$spec$tasks)
  lsg_assert_v7(sum(counts) == total, "Progress counts do not partition the work units.")
  elapsed <- max(0, now - p$started)
  processing_elapsed <- max(0, now - p$compute_started)
  eta <- if (phase == "fitting") lsg_progress_eta_v12(p, now) else if (phase == "completed") 0 else NA_real_
  deadline <- suppressWarnings(as.numeric(Sys.getenv("SLURM_JOB_END_TIME", NA_character_)))
  remaining <- if (length(deadline) == 1L && is.finite(deadline)) max(0, deadline - now) else NA_real_
  snapshot <- list(schema_version = "lsg_progress_v12", run_id = p$spec$version,
    scientific_signature = p$spec$scientific_signature,
    slurm_job_id = Sys.getenv("SLURM_JOB_ID", NA_character_), invocation_id = p$id,
    phase = phase, work_unit = "scenario_replication_all_methods_validated",
    total = total, completed = unname(counts[1L]), running = unname(counts[2L]),
    failed = unname(counts[3L]), pending = unname(counts[4L]),
    percent_completed = 100 * unname(counts[1L]) / total,
    elapsed_seconds = elapsed, processing_elapsed_seconds = processing_elapsed,
    completed_this_invocation = length(p$new_completed),
    throughput_tasks_per_second = if (processing_elapsed > 0 && length(p$new_completed)) length(p$new_completed) / processing_elapsed else NA_real_,
    eta_seconds = eta, estimated_completion_utc = utc(now + eta),
    eta_status = if (is.finite(eta)) "estimated" else "estimating_or_unknown",
    eta_scope = "fitting_and_task_validation_only_excludes_aggregation_and_final_verification",
    slurm_remaining_seconds = remaining, heartbeat_utc = utc(now),
    last_completed_work_unit = p$latest)
  path <- lsg_path_v7(p$output, "progress.json")
  temporary <- tempfile(".progress_", tmpdir = p$output)
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  writeLines(jsonlite::toJSON(snapshot, auto_unbox = TRUE, na = "null", digits = 16), temporary, useBytes = TRUE)
  lsg_assert_v7(file.rename(temporary, path), "Atomic progress snapshot publication failed.")
  history <- lsg_path_v7(p$output, "progress.tsv")
  utils::write.table(as.data.frame(snapshot, stringsAsFactors = FALSE), history,
    sep = "\t", row.names = FALSE, col.names = !file.exists(history), append = file.exists(history),
    quote = TRUE, na = "NA")
  cat(sprintf("[%s] run=%s job=%s phase=%s complete=%d/%d (%.1f%%) running=%d failed=%d pending=%d elapsed=%.1fs rate=%s tasks/s fitting_ETA=%s wall_remaining=%s last=%s\n",
    snapshot$heartbeat_utc, snapshot$run_id, snapshot$slurm_job_id, phase,
    snapshot$completed, total, snapshot$percent_completed, snapshot$running, snapshot$failed,
    snapshot$pending, elapsed, if (is.finite(snapshot$throughput_tasks_per_second)) sprintf("%.4f", snapshot$throughput_tasks_per_second) else "unknown",
    if (is.finite(eta)) sprintf("%.1fs", eta) else "estimating/unknown",
    if (is.finite(remaining)) sprintf("%.1fs", remaining) else "unknown", p$latest))
  flush.console(); p$last_emit <- now
  invisible(snapshot)
}
