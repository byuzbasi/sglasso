# Local, diagnostic-only ALL V5 prototype. V1--V4 functions remain unchanged.

allb_timed_environment_v5 <- function(e, elapsed) {
  stopifnot(is.environment(e), is.environment(elapsed))
  child <- new.env(parent = parent.env(e))
  for (name in ls(e, all.names = TRUE)) {
    child[[name]] <- get(name, envir = e, inherits = FALSE)
  }
  original_fitters <- e$lsg_finite_fitters_v7
  child$lsg_finite_fitters_v7 <- function(configuration) {
    fitters <- original_fitters(configuration)
    adelie <- fitters$adelie
    fitters$adelie <- function(...) {
      started <- proc.time()[["elapsed"]]
      result <- adelie(...)
      elapsed$adelie <- elapsed$adelie +
        proc.time()[["elapsed"]] - started
      result
    }
    fitters
  }
  external <- e$lsg_fit_external_joint_v7
  original_grpreg <- get("fit_grpreg_validation_paths_v2",
    envir = environment(external), inherits = TRUE)
  external_scope <- new.env(parent = environment(external))
  external_scope$fit_grpreg_validation_paths_v2 <- function(...) {
    started <- proc.time()[["elapsed"]]
    result <- original_grpreg(...)
    elapsed$grpreg <- elapsed$grpreg +
      proc.time()[["elapsed"]] - started
    result
  }
  environment(external) <- external_scope
  child$lsg_fit_external_joint_v7 <- function(...) {
    started <- proc.time()[["elapsed"]]
    result <- external(...)
    elapsed$external <- elapsed$external +
      proc.time()[["elapsed"]] - started
    result
  }
  for (entry in list(
    c("gab_estimate_training_firth_target", "firth"),
    c("lsg_fit_sglasso_joint_v7", "sglasso")
  )) {
    name <- entry[[1L]]
    stage <- entry[[2L]]
    original <- get(name, envir = e, inherits = TRUE)
    child[[name]] <- local({
      f <- original
      key <- stage
      function(...) {
        started <- proc.time()[["elapsed"]]
        result <- f(...)
        elapsed[[key]] <- elapsed[[key]] +
          proc.time()[["elapsed"]] - started
        result
      }
    })
  }
  child
}

allb_reduced_refit_configuration_v5 <- function(configuration, selections) {
  stopifnot(is.list(configuration), is.list(selections))
  get_one <- function(method) {
    row <- selections[[method]]
    stopifnot(is.data.frame(row), nrow(row) == 1L,
      row$eligible %in% TRUE)
    row
  }
  free <- get_one("Logistic SGLASSO")
  d0 <- get_one("Logistic SGLASSO (d=0 boundary)")
  adelie <- get_one("Logistic Group Elastic Net (adelie)")
  stopifnot(is.finite(free$alpha), is.finite(free$d),
    is.finite(d0$alpha), abs(d0$d) <= 1e-12,
    is.finite(adelie$alpha))
  reduced <- configuration
  reduced$alpha_grid <- sort(unique(c(free$alpha, d0$alpha)))
  reduced$d_grid <- sort(unique(c(0, free$d)))
  reduced$benchmark_alpha_grid <- sort(unique(c(0, adelie$alpha)))
  reduced
}

allb_refit_row_match_v5 <- function(candidates, selected, method) {
  stopifnot(is.data.frame(candidates), is.data.frame(selected),
    nrow(selected) == 1L)
  close <- function(x, y) is.finite(x) & is.finite(y) & abs(x - y) <= 1e-12
  match <- close(candidates$alpha, selected$alpha[[1L]]) &
    candidates$point_type == selected$point_type[[1L]] &
    candidates$lambda_index == selected$lambda_index[[1L]]
  if (grepl("SGLASSO", method, fixed = TRUE)) {
    match <- match & close(candidates$d, selected$d[[1L]])
  }
  if (identical(selected$point_type[[1L]], "finite")) {
    match <- match & close(candidates$lambda_relative_to_reference,
      selected$lambda_relative_to_reference[[1L]])
  } else {
    match <- match & is.infinite(candidates$lambda_relative_to_reference)
  }
  stopifnot(sum(match, na.rm = TRUE) == 1L)
  match[is.na(match)] <- FALSE
  match
}

allb_refit_selected_v5 <- function(e, X_train, y_train, X_test, group,
                                    configuration, inner_selection) {
  selections <- inner_selection$selections
  stopifnot(setequal(names(selections), gab_methods()))
  reduced <- allb_reduced_refit_configuration_v5(
    configuration, selections)
  original_candidate_frame <- gab_candidate_frame
  matched_candidate_frame <- function(bundle, method) {
    candidate <- original_candidate_frame(bundle, method)
    matched <- allb_refit_row_match_v5(
      candidate, selections[[method]], method)
    candidate$tuning_key <- paste0("v5_refit|", candidate$tuning_key)
    candidate$tuning_key[matched] <-
      selections[[method]]$tuning_key[[1L]]
    candidate
  }
  original_bundle <- gab_fit_bundle
  reduced_bundle <- function(e, X_train, y_train, X_validation,
                             y_validation, group, configuration) {
    original_bundle(e, X_train, y_train, X_validation,
      y_validation, group, reduced)
  }
  runner <- allb_clone_with_bindings_v2(gab_refit_selected, list(
    gab_fit_bundle = reduced_bundle,
    gab_candidate_frame = matched_candidate_frame
  ))
  result <- runner(e, X_train, y_train, X_test, group,
    configuration, inner_selection)
  for (method in gab_methods()) {
    row <- result$selected_rows$method == method
    stopifnot(sum(row) == 1L)
    for (name in c("alpha_index", "d_index", "fit_index")) {
      if (name %in% names(result$selected_rows) &&
          name %in% names(selections[[method]])) {
        result$selected_rows[[name]][row] <-
          selections[[method]][[name]][[1L]]
      }
    }
  }
  result$refit_grid_sizes <- c(
    sglasso_alpha = length(reduced$alpha_grid),
    sglasso_d = length(reduced$d_grid),
    adelie_alpha = length(reduced$benchmark_alpha_grid))
  result
}
