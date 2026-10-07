# Private, checked adapters around frozen code. No source file is rewritten.
# An exact one-node AST match is required: upstream drift fails at load time.
lsg_replace_expression_v13 <- function(fun, from, to) {
  count <- 0L
  walk <- function(x) {
    if (identical(x, from)) { count <<- count + 1L; return(to) }
    if (is.call(x)) for (i in seq_along(x)) {
      # [[<- NULL would delete an argument (e.g. rownames(x) <- NULL).
      if (!is.null(x[[i]]) && !identical(x[[i]], quote(expr = ))) x[[i]] <- walk(x[[i]])
    }
    x
  }
  updated <- walk(body(fun))
  if (count != 1L) stop("V13 adapter expected exactly one frozen expression; found ", count)
  body(fun) <- updated
  fun
}

lsg_classify_grpreg_path_v13 <- function(penalty, iterations, requested_path_length,
                                         max_iterations, warnings = character()) {
  if (identical(penalty, "grLasso")) {
    lsg_classify_grpreg_path_v9(penalty, iterations, requested_path_length, max_iterations, warnings)
  } else {
    v13_original$classify_grpreg_path_v2(penalty, iterations, requested_path_length, max_iterations, warnings)
  }
}

# Recompute the Group Lasso metadata from raw iteration/finite-point evidence.
# Unlike nonconvex saturation, this policy allows NO nonfinite returned point.
lsg_group_lasso_metadata_v13 <- function(tuning, require_selection_state = FALSE) {
  rows <- tuning$engine %in% "grpreg" & tuning$penalty_family %in% "group_lasso"
  x <- tuning[rows, , drop = FALSE]
  if (!nrow(x)) return(TRUE)
  required <- c("method_path", "fit_index", "lambda_index", "passes", "iteration_budget",
    "finite_coefficient", "finite_probability", "validation_log_loss", "solver_warning",
    "unexpected_solver_warning", "requested_path_length", "returned_path_length",
    "path_complete", "saturated_path_truncation", "iteration_budget_truncation",
    "total_iterations", "total_iteration_limit_reached", "returned_lower_boundary_excluded",
    "path_termination_acceptable", "point_converged", "finite_validation_loss", "numerically_eligible")
  if (!all(required %in% names(x))) return(FALSE)
  key <- do.call(paste, c(unname(x[c(intersect(c("scenario", "replication"), names(x)),
    "method_path", "fit_index")]), sep = "::"))
  all(vapply(split(x, key), function(z) {
    tryCatch({
      z <- z[order(z$lambda_index), , drop = FALSE]
      n <- nrow(z)
      scalar <- c("requested_path_length", "returned_path_length", "iteration_budget",
        "solver_warning", "unexpected_solver_warning", "total_iterations")
      if (!all(vapply(z[scalar], function(v) length(unique(v)) == 1L && !anyNA(v), logical(1))) ||
          !identical(as.numeric(z$lambda_index), as.numeric(seq_len(n))) ||
          z$returned_path_length[1L] != n || any(!is.finite(z$passes)) ||
          z$requested_path_length[1L] != floor(z$requested_path_length[1L]) ||
          z$iteration_budget[1L] != floor(z$iteration_budget[1L]) ||
          any(z$passes != floor(z$passes)) ||
          !all(z$finite_coefficient %in% TRUE & z$finite_probability %in% TRUE) ||
          !all(is.finite(z$validation_log_loss))) return(FALSE)
      warnings <- if (nzchar(z$solver_warning[1L])) strsplit(z$solver_warning[1L], " | ", fixed = TRUE)[[1L]] else character()
      status <- lsg_classify_grpreg_path_v9("grLasso", z$passes,
        z$requested_path_length[1L], z$iteration_budget[1L], warnings)
      if (!isTRUE(status$path_termination_acceptable)) return(FALSE)
      if (isTRUE(status$group_lasso_safe_prefix) &&
          any(head(cumsum(z$passes), -1L) >= z$iteration_budget[1L])) return(FALSE)
      fields <- c("path_complete", "saturated_path_truncation", "iteration_budget_truncation",
        "total_iteration_limit_reached", "returned_lower_boundary_excluded",
        "path_termination_acceptable", "point_converged")
      if (!all(vapply(fields, function(field) {
        identical(as.logical(z[[field]]), rep_len(as.logical(status[[field]]), n))
      }, logical(1))) || z$total_iterations[1L] != sum(z$passes) ||
          anyNA(z$unexpected_solver_warning) || any(nzchar(z$unexpected_solver_warning)) ||
          !all(z$finite_validation_loss %in% TRUE)) return(FALSE)
      eligible <- status$point_converged & !status$returned_lower_boundary_excluded
      if (!identical(as.logical(z$numerically_eligible), eligible)) return(FALSE)
      if (isTRUE(require_selection_state)) {
        if (!all(c("selected", "invalid_validation_contender") %in% names(z)) ||
            anyNA(z$selected) || anyNA(z$invalid_validation_contender) ||
            any(z$selected & !eligible) || any(z$invalid_validation_contender)) return(FALSE)
      }
      TRUE
    }, error = function(error) FALSE)
  }, logical(1)))
}

lsg_audit_external_v13 <- function(tuning, require_selection_state = FALSE) {
  result <- v13_extended_audit(tuning, require_selection_state)
  strict <- lsg_group_lasso_metadata_v13(tuning, require_selection_state)
  result$metadata_valid <- result$metadata_valid && strict
  result$audit_passes <- result$audit_passes && strict
  result
}

lsg_select_external_v13 <- function(tuning, method_path) {
  candidates <- tuning[tuning$method_path == method_path, , drop = FALSE]
  if (!lsg_group_lasso_metadata_v13(candidates)) stop("Invalid Group Lasso prefix evidence.")
  selected <- v13_original$select_valid_external_grid_v2(tuning, method_path)
  row <- selected$row
  if (row$engine == "grpreg" && row$penalty_family == "group_lasso" &&
      isTRUE(row$iteration_budget_truncation)) {
    selected$truncated_path_selection_interior <-
      row$lambda_index < row$returned_path_length - 1L
  }
  selected
}

lsg_external_path_rules_v13 <- function(et) {
  all(et$path_complete[et$engine == "adelie"] %in% TRUE) &&
    lsg_group_lasso_metadata_v13(et, require_selection_state = TRUE)
}

lsg_install_prefix_adapters_v13 <- function(e) {
  # Keep the native grpreg call, prediction, losses, order and tie policy exact.
  # Only retain evidence already computed inside its frame for downstream audit.
  # Locate the unique assignment rather than duplicating its large frozen body.
  matches <- list()
  visit <- function(x) {
    if (is.call(x)) {
      if (identical(x[[1L]], as.name("<-")) && identical(x[[2L]], quote(tuning[[pi]])) &&
          is.call(x[[3L]]) && identical(x[[3L]][[1L]], as.name("data.frame"))) {
        matches[[length(matches) + 1L]] <<- x
      }
      for (i in seq_along(x)) if (!identical(x[[i]], quote(expr = ))) visit(x[[i]])
    }
  }
  visit(body(e$v13_original$fit_grpreg_validation_paths_v2))
  if (length(matches) != 1L) stop("Native grpreg tuning assignment changed.")
  from <- matches[[1L]]; to <- from
  to[[3L]]$finite_coefficient <- quote(finite_coefficient)
  to[[3L]]$finite_probability <- quote(finite_probability)
  to[[3L]]$iteration_budget <- quote(max_iterations)
  e$fit_grpreg_validation_paths_v2 <- lsg_replace_expression_v13(
    e$v13_original$fit_grpreg_validation_paths_v2, from, to)
  # The extended structural audit is always paired with the stricter raw GL
  # evidence audit above. MCP/SCAD rows see the same expressions as before.
  e$v13_extended_audit <- lsg_replace_expression_v13(
    e$v13_original$audit_external_nonfinite_validation_v3,
    quote(tuning$penalty_family %in% c("group_mcp", "group_scad")),
    quote(tuning$penalty_family %in% c("group_mcp", "group_scad", "group_lasso")))
  e$v7_original$lsg_validate_tuning_payload_v7 <- lsg_replace_expression_v13(
    e$v7_original$lsg_validate_tuning_payload_v7,
    quote(all(et$path_complete[et$engine == "adelie" | et$penalty_family == "group_lasso"] %in% TRUE)),
    quote(lsg_external_path_rules_v13(et)))
  e$classify_grpreg_path_v2 <- e$lsg_classify_grpreg_path_v13
  e$audit_external_nonfinite_validation_v3 <- e$lsg_audit_external_v13
  e$select_valid_external_grid_v2 <- e$lsg_select_external_v13
  invisible(e)
}
