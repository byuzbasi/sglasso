# V7 reuses the validated V6 checkpoint format. Scientific signatures and
# release receipts differ; V6 shards cannot be resumed under V7.

allb_make_spec_v7 <- function(e, root, stage, version) {
  inherited <- allb_clone_with_bindings_v2(allb_make_spec_v3, list(
    allb_configuration_v3 = allb_configuration_v7,
    allb_source_files_v3 = allb_source_files_v7,
    allb_inventory = function(e, root, files) {
      e$lsg_inventory_v7(root, sort(unique(files), method = "radix"))
    }
  ))
  spec <- inherited(e, root, stage, version)
  spec$schema_version <- "all_binary_study_v7"
  spec$scientific_signature <- NULL
  identity <- spec[!names(spec) %in% c("runtime", "created_utc")]
  spec$scientific_signature <- allb_hash(e, identity)
  spec
}

allb_validate_spec_v7 <- function(e, spec, root, stage, version) {
  current <- allb_make_spec_v7(e, root, stage, version)
  fields <- c("schema_version", "version", "stage", "configuration",
    "tasks", "data_identity", "source_manifest", "scientific_signature",
    "runtime")
  allb_assert(identical(spec[fields], current[fields]),
    "ALL V7 identity, sources, data, runtime, grids, or seeds changed.")
  invisible(TRUE)
}

allb_validate_release_v7 <- function(e, root, version) {
  path <- file.path(root, "release", "all_binary_v7",
    paste0("LOCAL_VALIDATED_", version, ".rds"))
  allb_assert(file.exists(path),
    "ALL V7 local smoke receipt is missing; production is blocked.")
  receipt <- readRDS(path)
  current <- allb_make_spec_v7(e, root, "production", version)
  fields <- c("schema_version", "version", "stage", "configuration",
    "tasks", "data_identity", "source_manifest", "scientific_signature",
    "runtime")
  allb_assert(identical(receipt$schema_version,
    "all_binary_local_release_v7") && isTRUE(receipt$accepted) &&
    all(receipt$checks$passed) &&
    identical(receipt$identity[fields], current[fields]),
    "ALL V7 local release does not match the production identity.")
  invisible(receipt)
}

allb_run_v7 <- allb_clone_with_bindings_v2(allb_run_v6, list(
  allb_load_environment_v6 = allb_load_environment_v7,
  allb_make_spec_v6 = allb_make_spec_v7,
  allb_validate_spec_v6 = allb_validate_spec_v7,
  allb_validate_release_v6 = allb_validate_release_v7
))
