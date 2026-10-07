# ALL V6 changes only the Logistic SGLASSO numerical implementation.
# The V4 endpoint, data, splits, methods, grids, seeds and gates are frozen.

allb_source_files_v6 <- function() unique(c(
  allb_source_files_v4(),
  "src/logistic_sglasso_block_kernel_local_v3.hpp",
  "src/logistic_sglasso_hybrid_solver_local_apg_v4.cpp",
  "R/all_binary_workflow_v6.R",
  "R/all_binary_io_v6.R",
  "scripts/96_run_all_binary_v6.R",
  "scripts/97_validate_all_binary_v6.R",
  "tests/test_all_binary_v6.R",
  "ALL_BINARY_IRLS_PROTOCOL_V6.md"
))

allb_load_environment_v6 <- function(root) {
  e <- allb_load_environment_v4(root)
  allb_assert(requireNamespace("Rcpp", quietly = TRUE) &&
    requireNamespace("RcppArmadillo", quietly = TRUE),
    "Rcpp and RcppArmadillo are required; no package is installed.")
  root <- normalizePath(root, mustWork = TRUE)
  source_file <- file.path(root, "src",
    "logistic_sglasso_hybrid_solver_local_apg_v4.cpp")
  header_file <- file.path(root, "src",
    "logistic_sglasso_block_kernel_local_v3.hpp")
  paths <- c(source_file, header_file)
  allb_assert(all(file.exists(paths)) && !any(dir.exists(paths)) &&
    !any(nzchar(Sys.readlink(paths))),
    "Missing regular ALL V6 IRLS source or header.")

  previous_cppflags <- Sys.getenv("PKG_CPPFLAGS", unset = NA_character_)
  previous_makevars <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
  include_flag <- paste0("-I", shQuote(file.path(root, "src")))
  Sys.setenv(PKG_CPPFLAGS = if (is.na(previous_cppflags) ||
    !nzchar(previous_cppflags)) include_flag else
      paste(previous_cppflags, include_flag))
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

  # Both R entry points resolve their C++ counterparts in e. Wrapping the
  # compiled functions therefore covers the finite path and all three
  # selected-solution starts without altering the frozen R fitting API.
  Rcpp::sourceCpp(source_file, rebuild = FALSE, showOutput = FALSE,
    verbose = FALSE, env = e)
  required <- c("lsg_fit_one_hybrid_v11_cpp", "lsg_path_hybrid_v11_cpp")
  allb_assert(all(vapply(required, exists, logical(1), mode = "function",
    envir = e, inherits = FALSE)), "ALL V6 IRLS entry points were not loaded.")
  path_cpp <- e$lsg_path_hybrid_v11_cpp
  single_cpp <- e$lsg_fit_one_hybrid_v11_cpp
  e$lsg_path_hybrid_v11_cpp <- function(...) path_cpp(...,
    reuse_apg_offsets = TRUE, use_irls = TRUE)
  e$lsg_fit_one_hybrid_v11_cpp <- function(...) single_cpp(...,
    reuse_apg_offsets = TRUE, use_irls = TRUE)
  e$allb_v6_core_sha256 <- vapply(paths, digest::digest, character(1),
    file = TRUE, algo = "sha256")
  e
}

allb_configuration_v6 <- function(e, stage) {
  allb_configuration_v4(e, stage)
}

allb_run_task_v6 <- function(e, data, task, configuration) {
  allb_run_task_v4(e, data, task, configuration)
}

allb_validate_payload_v6 <- function(payload, task, configuration) {
  allb_validate_payload_v4(payload, task, configuration)
}
