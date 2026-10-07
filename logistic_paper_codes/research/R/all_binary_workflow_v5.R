# ALL V5 is a local validation candidate for the optimized V11-equivalent core.
# It inherits the frozen V4 data, design, estimators, tuning, and diagnostics.

allb_source_files_v5 <- function() unique(c(
  allb_source_files_v4(),
  "src/logistic_sglasso_block_kernel_allb_v5.hpp",
  "src/logistic_sglasso_hybrid_solver_allb_v5.cpp",
  "R/all_binary_workflow_v5.R",
  "R/all_binary_io_v5.R",
  "scripts/94_run_all_binary_v5.R",
  "scripts/95_validate_all_binary_v5.R",
  "tests/test_all_binary_v5.R",
  "ALL_BINARY_OPTIMIZED_CORE_PROTOCOL_V5.md"
))

allb_compile_optimized_core_v5 <- function(e, root) {
  allb_assert(is.environment(e), "ALL V5 requires a workflow environment.")
  allb_assert(requireNamespace("Rcpp", quietly = TRUE) &&
    requireNamespace("RcppArmadillo", quietly = TRUE),
    "Rcpp and RcppArmadillo are required; no installation is attempted.")
  root <- normalizePath(root, mustWork = TRUE)
  source_file <- file.path(root, "src",
    "logistic_sglasso_hybrid_solver_allb_v5.cpp")
  header_file <- file.path(root, "src",
    "logistic_sglasso_block_kernel_allb_v5.hpp")
  paths <- c(source_file, header_file)
  allb_assert(all(file.exists(paths)) && !any(dir.exists(paths)) &&
    !any(nzchar(Sys.readlink(paths))),
    "Missing regular ALL V5 optimized core source or header.")

  previous_cppflags <- Sys.getenv("PKG_CPPFLAGS", unset = NA_character_)
  include_flag <- paste0("-I", shQuote(file.path(root, "src")))
  Sys.setenv(PKG_CPPFLAGS = if (is.na(previous_cppflags) ||
    !nzchar(previous_cppflags)) include_flag else
      paste(previous_cppflags, include_flag))
  previous_makevars <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  local_makevars <- file.path(root, "config", "Makevars.local")
  local_gfortran_runtime <- paste0(
    "/usr/local/gfortran/lib/gcc/",
    "aarch64-apple-darwin23/14.1.0/libemutls_w.a")
  if (file.exists(local_makevars) && file.exists(local_gfortran_runtime)) {
    Sys.setenv(R_MAKEVARS_USER = local_makevars)
  }
  on.exit({
    if (is.na(previous_cppflags)) Sys.unsetenv("PKG_CPPFLAGS")
    else Sys.setenv(PKG_CPPFLAGS = previous_cppflags)
    if (is.na(previous_makevars)) Sys.unsetenv("R_MAKEVARS_USER")
    else Sys.setenv(R_MAKEVARS_USER = previous_makevars)
  }, add = TRUE)

  # The V4 loader has already populated e with the original C++ entry points.
  # sourceCpp replaces only those entry points in e; the R fitting API stays
  # unchanged, including V11 numerical controls and the selected-fit audit.
  Rcpp::sourceCpp(source_file, rebuild = FALSE, showOutput = FALSE,
    verbose = FALSE, env = e)
  required <- c("lsg_shifted_group_prox_v11_cpp",
    "lsg_fit_one_hybrid_v11_cpp", "lsg_path_hybrid_v11_cpp")
  allb_assert(all(vapply(required, exists, logical(1), mode = "function",
    envir = e, inherits = FALSE)),
    "ALL V5 optimized C++ entry points were not loaded.")
  e$allb_v5_core_source_sha256 <- digest::digest(file = source_file,
    algo = "sha256")
  e$allb_v5_core_header_sha256 <- digest::digest(file = header_file,
    algo = "sha256")
  invisible(TRUE)
}

allb_load_environment_v5 <- function(root) {
  e <- allb_load_environment_v4(root)
  allb_compile_optimized_core_v5(e, root)
  e
}

# Preserve V4 statistical configuration exactly; implementation identity is
# carried by the V5 source manifest and study scientific signature.
allb_configuration_v5 <- function(e, stage) {
  allb_configuration_v4(e, stage)
}

allb_run_task_v5 <- function(e, data, task, configuration) {
  allb_run_task_v4(e, data, task, configuration)
}

allb_validate_payload_v5 <- function(payload, task, configuration) {
  allb_validate_payload_v4(payload, task, configuration)
}
