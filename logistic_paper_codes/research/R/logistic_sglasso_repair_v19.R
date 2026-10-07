# V19 production repair overlay. Frozen V15 and V18 files are never rewritten.

lsg_lambda_relative_grid_v19 <- function() {
  c(8192, 4096, v19_base$lsg_lambda_relative_grid_v7())
}

lsg_configuration_v19 <- function(stage) {
  cfg <- v19_base$lsg_configuration_v7(stage)
  grid <- lsg_lambda_relative_grid_v19()
  cfg$nlambda <- length(grid)
  cfg$lambda_relative_grid <- grid
  cfg$lambda_extension_multipliers <- grid[seq_len(10L)]
  cfg$lambda_upper_multiplier <- grid[1L]
  cfg$integration_version <- "v19_fixed_8192_tail_complete_grlasso_audit"
  cfg$group_lasso_prefix_policy <-
    "complete_path_final_total_budget_point_excluded_audited_v19"
  cfg
}

lsg_finite_fitters_v19 <- function(configuration) {
  expected <- lsg_lambda_relative_grid_v19()
  if (!isTRUE(all.equal(as.numeric(configuration$lambda_relative_grid),
                        expected, tolerance = 1e-13)) ||
      !identical(as.integer(configuration$nlambda), length(expected)) ||
      !isTRUE(all.equal(as.numeric(configuration$lambda_upper_multiplier),
                        8192)) ||
      !isTRUE(all.equal(as.numeric(configuration$lambda_min_ratio), 0.05))) {
    stop("V19 requires the frozen 40-point finite lambda grid.",
         call. = FALSE)
  }
  clone <- function(original) {
    child <- new.env(parent = environment(original))
    child$logistic_sglasso_lambda_relative_grid_v5 <- local({
      frozen_grid <- expected
      function() frozen_grid
    })
    result <- original
    environment(result) <- child
    result
  }
  fitters <- list(
    sglasso = clone(fit_extended_sglasso_validation_grid_v5),
    adelie = clone(fit_adelie_validation_grid_v5)
  )

  # Preserve the V18 production dispatch: the frozen V5 grid builder calls
  # the V11 Rcpp/RcppArmadillo hybrid solver through this isolated binding.
  hybrid <- new.env(parent = environment(fitters$sglasso))
  hybrid$fit_logistic_sglasso <- function(X, y, group, lambda, d, alpha,
      max_outer, max_inner, tolerance, inner_tolerance, target_original,
      compile, preprocess, use_active_set, warm_start_d, solver) {
    lsg_assert_v7(
      identical(solver, "hybrid_v11_2"),
      "Unexpected V19 integration solver."
    )
    fit <- lsg_fit_hybrid_path_v11(
      X, y, group, lambda, d, alpha, target_original,
      preprocess = preprocess, controls = configuration$path_controls,
      compile = FALSE, use_active_set = use_active_set,
      warm_start_d = warm_start_d
    )
    fit$passes <- fit$block_sweeps + fit$apg_iterations
    fit
  }
  environment(fitters$sglasso) <- hybrid
  fitters
}

lsg_selection_boundary_v19 <- function(selected, configuration) {
  finite <- identical(selected$point_type[[1L]], "finite")
  reference_type <- selected$lambda_reference_type[[1L]]
  extended <- reference_type %in% c("null_score_ridge_boundary", "d0_null_kkt")
  ratio <- selected$lambda_relative_to_reference[[1L]]
  upper <- finite && extended && is.finite(ratio) &&
    abs(ratio - configuration$lambda_upper_multiplier) <= 1e-10
  lower <- finite && extended && is.finite(ratio) &&
    abs(ratio - configuration$lambda_min_ratio) <= 1e-10
  list(
    upper = upper, lower = lower,
    native_upper = finite && !extended && is.finite(ratio) &&
      abs(ratio - 1) <= 1e-10,
    native_lower = finite && !extended && is.finite(ratio) &&
      abs(ratio - configuration$lambda_min_ratio) <= 1e-10,
    unresolved = upper || lower
  )
}

# Return TRUE only for rows belonging to the exact complete-path Group Lasso
# case introduced in V18: the total iteration budget is reached for the first
# time at the final requested point, and that final point alone is excluded.
lsg_group_lasso_complete_final_budget_rows_v19 <- function(
    tuning, require_selection_state = FALSE
) {
  result <- rep(FALSE, nrow(tuning))
  rows <- tuning$engine %in% "grpreg" &
    tuning$penalty_family %in% "group_lasso"
  if (!any(rows)) return(result)
  key_columns <- c(
    intersect(c("scenario", "replication"), names(tuning)),
    "method_path", "fit_index"
  )
  keys <- do.call(paste, c(unname(tuning[key_columns]), list(sep = "::")))
  for (index in split(which(rows), keys[rows])) {
    z <- tuning[index, , drop = FALSE]
    accepted <- tryCatch({
      z <- z[order(z$lambda_index), , drop = FALSE]
      n <- nrow(z)
      budget <- unique(z$iteration_budget)
      requested <- unique(z$requested_path_length)
      returned <- unique(z$returned_path_length)
      isTRUE(lsg_group_lasso_metadata_v13(
        z, require_selection_state = require_selection_state
      )) &&
        length(budget) == 1L && length(requested) == 1L &&
        length(returned) == 1L && is.finite(budget) && budget > 0 &&
        identical(as.numeric(requested), as.numeric(n)) &&
        identical(as.numeric(returned), as.numeric(n)) &&
        all(z$path_complete %in% TRUE) &&
        all(z$iteration_budget_truncation %in% TRUE) &&
        all(z$saturated_path_truncation %in% FALSE) &&
        all(z$total_iteration_limit_reached %in% TRUE) &&
        all(is.finite(z$passes)) && all(z$passes >= 0) &&
        all(z$passes == floor(z$passes)) &&
        identical(which(z$returned_lower_boundary_excluded %in% TRUE), n) &&
        identical(as.numeric(z$total_iterations), rep(sum(z$passes), n)) &&
        identical(as.numeric(sum(z$passes)), as.numeric(budget)) &&
        all(head(cumsum(z$passes), -1L) < budget)
    }, error = function(error) FALSE)
    if (isTRUE(accepted)) result[index] <- TRUE
  }
  result
}

# Patch one named assignment in the inherited structural audit. This keeps its
# return schema, row order and all non-Group-Lasso logic byte-for-byte intact.
lsg_build_extended_audit_v19 <- function() {
  original <- v19_base$v13_extended_audit
  matches <- list()
  visit <- function(x) {
    if (is.call(x)) {
      if (identical(x[[1L]], as.name("<-")) &&
          identical(x[[2L]], as.name("truncation_metadata_valid"))) {
        matches[[length(matches) + 1L]] <<- x
      }
      for (i in seq_along(x)) {
        if (!is.null(x[[i]]) && !identical(x[[i]], quote(expr = ))) {
          visit(x[[i]])
        }
      }
    }
  }
  visit(body(original))
  if (length(matches) != 1L) {
    stop("V19 expected one inherited truncation audit assignment.",
         call. = FALSE)
  }
  from <- matches[[1L]]
  to <- substitute(
    truncation_metadata_valid <- OLD |
      lsg_group_lasso_complete_final_budget_rows_v19(
        tuning, require_selection_state
      ),
    list(OLD = from[[3L]])
  )
  patched <- lsg_replace_expression_v13(original, from, to)
  child <- new.env(parent = environment(patched))
  child$lsg_group_lasso_complete_final_budget_rows_v19 <-
    lsg_group_lasso_complete_final_budget_rows_v19
  environment(patched) <- child
  patched
}

lsg_build_external_audit_v19 <- function() {
  extended <- lsg_build_extended_audit_v19()
  audit <- v19_base$audit_external_nonfinite_validation_v3
  child <- new.env(parent = environment(audit))
  child$v13_extended_audit <- extended
  environment(audit) <- child
  audit
}

lsg_tuning_unit_checks_v19 <- function() {
  check <- function(name, passed) data.frame(
    check = name, passed = isTRUE(passed), stringsAsFactors = FALSE
  )
  inherited <- v19_base$lsg_tuning_unit_checks_v7()
  grid <- lsg_lambda_relative_grid_v19()
  configuration <- list(
    lambda_relative_grid = grid, nlambda = length(grid),
    lambda_upper_multiplier = 8192, lambda_min_ratio = 0.05
  )
  old_sg <- fit_extended_sglasso_validation_grid_v5
  old_ad <- fit_adelie_validation_grid_v5
  old_options <- options()
  fitters <- lsg_finite_fitters_v19(configuration)
  checks <- list(
    check("v19_inherited_v18_tuning_checks", all(inherited$passed)),
    check("v19_cloned_finite_fitters_preserve_bodies",
      identical(body(fitters$sglasso), body(old_sg)) &&
        identical(body(fitters$adelie), body(old_ad))),
    check("v19_no_global_frozen_binding_or_option_changes",
      identical(old_options, options()) &&
        identical(old_sg, fit_extended_sglasso_validation_grid_v5) &&
        identical(old_ad, fit_adelie_validation_grid_v5)),
    check("v19_isolated_40_point_grid",
      identical(get("logistic_sglasso_lambda_relative_grid_v5",
                    envir = environment(fitters$sglasso),
                    inherits = TRUE)(), grid) &&
        identical(get("logistic_sglasso_lambda_relative_grid_v5",
                      envir = environment(fitters$adelie),
                      inherits = TRUE)(), grid) &&
        identical(grid[3:40], v19_base$lsg_lambda_relative_grid_v7())),
    check("v19_hybrid_dispatch_layer_present",
      exists("fit_logistic_sglasso", envir = environment(fitters$sglasso),
             inherits = FALSE) &&
        grepl("lsg_fit_hybrid_path_v11", paste(deparse(body(get(
          "fit_logistic_sglasso", envir = environment(fitters$sglasso),
          inherits = FALSE))), collapse = "\n"), fixed = TRUE))
  )
  toy <- data.frame(
    point_type = c("finite", "finite", "penalty_limit"),
    lambda_reference_type = "d0_null_kkt",
    lambda_relative_to_reference = c(8192, 4096, Inf),
    stringsAsFactors = FALSE
  )
  checks[[length(checks) + 1L]] <- check(
    "v19_only_8192_is_finite_upper_boundary",
    lsg_selection_boundary_v19(toy[1L, ], configuration)$unresolved &&
      !lsg_selection_boundary_v19(toy[2L, ], configuration)$unresolved &&
      !lsg_selection_boundary_v19(toy[3L, ], configuration)$unresolved
  )
  bad <- configuration
  bad$lambda_relative_grid[1L] <- 16384
  checks[[length(checks) + 1L]] <- check(
    "v19_unapproved_grid_rejected",
    inherits(tryCatch(lsg_finite_fitters_v19(bad), error = identity), "error")
  )
  do.call(rbind, c(list(inherited), checks))
}
