# V7 extends only the shared lower relative-lambda range. The V6 solver,
# observations, preprocessing, methods, folds, seeds, and gates are unchanged.

allb_source_files_v7 <- function() unique(c(
  allb_source_files_v6(),
  "R/all_binary_workflow_v7.R", "R/all_binary_io_v7.R",
  "scripts/99_run_all_binary_v7.R",
  "scripts/100_validate_all_binary_v7.R",
  "tests/test_all_binary_v7.R",
  "ALL_BINARY_LOWER_LAMBDA_PROTOCOL_V7.md"
))

allb_load_environment_v7 <- function(root) {
  e <- allb_load_environment_v6(root)
  original <- e$allb_common_lambda_grid_v2()
  allb_assert(length(original) == 34L &&
    isTRUE(all.equal(tail(original, 1L), 0.003125)),
    "The frozen V6 common lambda grid changed unexpectedly.")
  extended <- c(original, 0.0015625, 0.00078125)
  e$allb_common_lambda_grid_v2 <- function() extended
  e$allb_validate_lambda_grids_v2 <- function() {
    common <- e$allb_common_lambda_grid_v2()
    sglasso <- e$allb_sglasso_lambda_grid_v2()
    allb_assert(identical(common, extended) && length(sglasso) == 49L &&
      all(is.finite(common)) && all(common > 0) && all(diff(common) < 0) &&
      all(is.finite(sglasso)) && all(sglasso > 0) &&
      all(diff(sglasso) < 0) &&
      all(vapply(common, function(value) {
        any(abs(sglasso - value) <= 1e-13)
      }, logical(1))), "Invalid V7 extended lambda grids.")
    invisible(TRUE)
  }
  e$allb_validate_lambda_grids_v2()
  allb_assert(identical(e$allb_sglasso_lambda_grid_v2(),
    c(head(e$allb_sglasso_lambda_grid_v2(), 13L), extended)),
    "The SGLASSO path did not inherit the extended common grid.")
  e
}

allb_configuration_v7 <- function(e, stage) {
  configuration <- allb_configuration_v6(e, stage)
  common <- e$allb_common_lambda_grid_v2()
  sglasso <- e$allb_sglasso_lambda_grid_v2()
  allb_assert(length(common) == 36L &&
    identical(tail(common, 3L), c(0.003125, 0.0015625, 0.00078125)) &&
    length(sglasso) == 49L &&
    e$allb_grid_equal_v2(configuration$fit_configuration$lambda_relative_grid,
      sglasso) &&
    e$allb_grid_equal_v2(
      configuration$fit_configuration$fair_common_lambda_relative_grid,
      common), "V7 configuration did not propagate the extended grid.")
  configuration$schema_version <- "all_binary_configuration_v7"
  configuration$fit_configuration$fair_lower_boundary_policy <-
    "selected_0.00078125_blocks_task_and_finalization"
  configuration$fit_configuration$integration_version <-
    "all_binary_lower_lambda_v7"
  configuration$fair_tuning$common_relative_grid <- common
  configuration$fair_tuning$sglasso_relative_grid <- sglasso
  configuration$fair_tuning$solver_boundary_rule <- paste(
    "truncated_returned_endpoint_is_descriptive",
    "only_true_common_endpoint_0.00078125_is_rejected", sep = "_"
  )
  configuration
}
