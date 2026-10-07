# Telemetry only. No objective, seed, candidate, acceptance rule or worker changes.
lsg_progress_new_v15 <- function(output, spec, completed, cores, interval = 30,
                                now = as.numeric(Sys.time())) {
  p <- v15_base$lsg_progress_new_v12(output, spec, completed, cores, interval, now)
  p$phase <- "preflight"
  v15_progress <<- p
  p
}

lsg_progress_emit_v15 <- function(p, phase, force = FALSE, now = as.numeric(Sys.time())) {
  p$phase <- phase
  v15_base$lsg_progress_emit_v12(p, phase, force, now)
}

lsg_progress_eta_v15 <- function(p, now) {
  if (any(p$status == "failed")) return(NA_real_)
  todo <- which(p$status %in% c("running", "pending"))
  if (!length(todo)) return(0)
  if (length(p$new_completed) < 2L) return(NA_real_)
  scenario <- p$spec$tasks$scenario
  means <- tapply(p$durations[p$new_completed], scenario[p$new_completed], mean)
  cost <- as.numeric(means[scenario[todo]])
  if (any(!is.finite(cost)) || any(cost <= 0)) return(NA_real_)
  running <- p$status[todo] == "running"
  remaining <- cost[running] - (now - p$task_started[todo[running]])
  # An overrun invalidates the mean-duration ETA; never claim zero time left
  # while unfinished work remains, nor invent a positive time floor.
  if (any(!is.finite(remaining)) || any(remaining <= 0)) return(NA_real_)
  load <- c(remaining, rep(0, max(0L, p$cores - length(remaining))))
  for (seconds in cost[!running]) {
    k <- which.min(load); load[k] <- load[k] + seconds
  }
  if (!length(load) || max(load) <= 0) NA_real_ else max(load)
}

lsg_read_shard_progress_v15 <- function(output, spec, task, deep = TRUE) {
  shard <- v15_base$lsg_read_shard_v7(output, spec, task, deep)
  p <- v15_progress
  if (is.environment(p) && identical(p$output, output) &&
      identical(p$spec$scientific_signature, spec$scientific_signature)) {
    if (identical(p$phase, "checkpoint_validation")) {
      i <- match(task$key, p$spec$tasks$key)
      if (!is.na(i)) {
        p$status[i] <- if (all(shard$checks$passed)) "completed" else "failed"
        if (p$status[i] == "completed") p$latest <- task$key
      }
    }
    # Collection and final verification are per-shard: emit only at the
    # configured interval, from the same parent that already owns telemetry.
    if (p$phase %in% c("checkpoint_validation", "aggregation_and_validation"))
      lsg_progress_emit_v12(p, p$phase)
  }
  shard
}

lsg_scan_existing_v15 <- function(output, spec, cores, interval) {
  p <- lsg_progress_new_v12(output, spec, integer(), cores, interval)
  lsg_progress_emit_v12(p, "checkpoint_validation", TRUE)
  answer <- v15_base$lsg_collect_v7(output, spec, deep = TRUE, require_all = FALSE)
  lsg_progress_emit_v12(p, "checkpoint_validation", TRUE)
  answer
}
