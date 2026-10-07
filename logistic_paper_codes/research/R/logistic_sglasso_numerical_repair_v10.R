# V10 diagnostics extend the V9 finite lambda path and audit selected
# solutions more tightly.  The mathematical objective and V9 RcppArmadillo
# joint solver are unchanged.

lsg_lambda_relative_grid_v10 <- function() {
  frozen_v9 <- lsg_lambda_relative_grid_v9()
  if (length(frozen_v9) != 39L ||
      !isTRUE(all.equal(tail(frozen_v9, 1L), 0.05,
                        tolerance = 1e-14))) {
    stop("The frozen V9 lambda grid is unavailable or changed.",
         call. = FALSE)
  }
  log_step <- frozen_v9[length(frozen_v9)] /
    frozen_v9[length(frozen_v9) - 1L]
  lower_tail <- tail(frozen_v9, 1L) * log_step^seq_len(24L)
  c(frozen_v9, lower_tail)
}


lsg_objective_reconstruction_audit_v10 <- function(
    observed,
    reference,
    ulp_factor = 32
) {
  observed <- as.numeric(observed)
  reference <- as.numeric(reference)
  if (length(observed) != length(reference) || !length(observed) ||
      length(ulp_factor) != 1L || is.na(ulp_factor) ||
      !is.finite(ulp_factor) || ulp_factor <= 0) {
    stop("Invalid V10 objective-reconstruction audit input.",
         call. = FALSE)
  }
  scale <- 1 + abs(reference)
  unit <- .Machine$double.eps * scale
  error <- abs(observed - reference)
  tolerance <- ulp_factor * unit
  data.frame(
    objective_reconstruction_error = error,
    objective_reconstruction_tolerance = tolerance,
    objective_reconstruction_ulp_units = error / unit,
    objective_reconstruction_pass = is.finite(observed) &
      is.finite(reference) & is.finite(error) & error <= tolerance,
    stringsAsFactors = FALSE
  )
}


# V10 deliberately reuses the already validated V9 RcppArmadillo core.  These
# wrappers give the new diagnostic an explicit version boundary without
# copying or changing the estimator implementation.
lsg_compile_joint_solver_v10 <- function(...) {
  lsg_compile_joint_solver_v9(...)
}


lsg_fit_joint_path_v10 <- function(...) {
  lsg_fit_joint_path_v9(...)
}


lsg_fit_one_joint_v10_cpp <- function(...) {
  lsg_fit_one_joint_v9_cpp(...)
}


lsg_classify_grpreg_path_v10 <- function(...) {
  lsg_classify_grpreg_path_v9(...)
}


lsg_select_grpreg_prefix_v10 <- function(...) {
  lsg_select_grpreg_prefix_v9(...)
}


lsg_portable_numeric_equal_v10 <- function(...) {
  lsg_portable_numeric_equal_v9(...)
}


lsg_joint_solver_unit_checks_v10 <- function(root = lsg_prework_root()) {
  inherited <- lsg_joint_solver_unit_checks_v9(root)
  grid <- lsg_lambda_relative_grid_v10()
  frozen_v9 <- lsg_lambda_relative_grid_v9()
  reference <- c(1, 48000, 5e8)
  near <- reference + 7 * .Machine$double.eps * (1 + abs(reference))
  material <- reference + 128 * .Machine$double.eps *
    (1 + abs(reference))
  near_audit <- lsg_objective_reconstruction_audit_v10(near, reference)
  material_audit <- lsg_objective_reconstruction_audit_v10(
    material, reference
  )
  defaults <- lsg_repair_defaults_v10()$configuration
  checks <- c(
    v9_solver_unit_checks_preserved = all(inherited$passed %in% TRUE),
    frozen_v9_grid_is_exact_prefix = length(grid) == 63L &&
      identical(grid[seq_along(frozen_v9)], frozen_v9),
    lower_tail_has_24_fixed_points = length(grid) - length(frozen_v9) == 24L,
    extended_grid_strictly_decreasing = !is.unsorted(
      -grid, strictly = TRUE
    ),
    extended_grid_positive_and_finite = all(is.finite(grid)) && all(grid > 0),
    extended_grid_endpoint_frozen = isTRUE(all.equal(
      tail(grid, 1L), 0.0007496709952762384, tolerance = 2e-14
    )),
    seven_ulp_roundoff_is_accepted = all(
      near_audit$objective_reconstruction_pass
    ),
    material_ulp_change_is_rejected = all(
      !material_audit$objective_reconstruction_pass
    ),
    objective_gate_is_scale_aware = all(
      near_audit$objective_reconstruction_tolerance > 0
    ) && max(near_audit$objective_reconstruction_ulp_units) < 16,
    polishing_kkt_is_stricter =
      defaults$polish_kkt_tolerance < defaults$joint_kkt_tolerance &&
      defaults$polish_kkt_limit < defaults$study_kkt_limit,
    polishing_budget_is_larger =
      defaults$polish_max_sweeps > defaults$joint_max_sweeps
  )
  out <- rbind(
    transform(inherited, check = paste0("inherited_v9__", check)),
    data.frame(
      check = names(checks),
      passed = unname(vapply(checks, isTRUE, logical(1))),
      stringsAsFactors = FALSE
    )
  )
  rownames(out) <- NULL
  out
}
