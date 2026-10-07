# Read-only compatibility audit for the immutable V6 diagnostic archive.
# This file neither fits models nor writes, replaces, or re-signs any artifact.

lsg_portable_assert_v1 <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
  invisible(TRUE)
}

lsg_portable_manifest_schema_v1 <- function(manifest) {
  lsg_portable_assert_v1(
    is.data.frame(manifest) && identical(names(manifest), c("file", "bytes", "sha256")) &&
      nrow(manifest) > 0L && is.character(manifest$file) && !anyNA(manifest$file) &&
      !anyDuplicated(manifest$file) && is.numeric(manifest$bytes) &&
      all(is.finite(manifest$bytes)) && all(manifest$bytes >= 0) &&
      all(manifest$bytes == floor(manifest$bytes)) && is.character(manifest$sha256) &&
      !anyNA(manifest$sha256) && all(grepl("^[0-9a-f]{64}$", manifest$sha256)),
    "Invalid portable manifest filename/byte/SHA-256 schema.")
  relative <- lsg_tail_relative_v6(manifest$file)
  lsg_portable_assert_v1(identical(relative, manifest$file),
                          "Manifest paths must be project-relative, not host-specific.")
  invisible(TRUE)
}

lsg_portable_manifest_v1 <- function(root, paths) {
  relative <- lsg_tail_relative_v6(paths)
  lsg_portable_assert_v1(length(relative) > 0L && !anyNA(relative) && !anyDuplicated(relative),
                          "Manifest inputs must be nonempty, unique file paths.")
  relative <- sort(relative, method = "radix")
  actual <- file.path(dirname(normalizePath(root, mustWork = TRUE)), relative)
  lsg_portable_assert_v1(all(file.exists(actual)) && !any(dir.exists(actual)),
                          "A required portable-audit file is missing.")
  result <- data.frame(
    file = relative, bytes = as.numeric(file.info(actual)$size),
    sha256 = unname(vapply(actual, digest::digest, character(1), file = TRUE, algo = "sha256")),
    stringsAsFactors = FALSE)
  lsg_portable_manifest_schema_v1(result)
  result
}

lsg_portable_same_manifest_v1 <- function(left, right) {
  lsg_portable_manifest_schema_v1(left)
  lsg_portable_manifest_schema_v1(right)
  lsg_portable_assert_v1(setequal(left$file, right$file), "Manifest file inventories differ.")
  index <- match(left$file, right$file)
  lsg_portable_assert_v1(all(left$bytes == right$bytes[index]) &&
                          identical(left$sha256, right$sha256[index]),
                          "Manifest bytes or SHA-256 values differ for a matched filename.")
  invisible(TRUE)
}

lsg_portable_verify_manifest_v1 <- function(root, manifest, expected_relative = NULL) {
  lsg_portable_manifest_schema_v1(manifest)
  if (!is.null(expected_relative)) {
    expected <- lsg_tail_relative_v6(expected_relative)
    lsg_portable_assert_v1(!anyDuplicated(expected) && setequal(manifest$file, expected),
                            "Manifest inventory differs from the required artifacts.")
  }
  observed <- lsg_portable_manifest_v1(root, manifest$file)
  lsg_portable_same_manifest_v1(manifest, observed)
  invisible(TRUE)
}

lsg_portable_read_raw_v1 <- function(path) {
  connection <- gzfile(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  value <- readBin(connection, what = "raw", n = 10000001L)
  lsg_portable_assert_v1(length(value) <= 10000000L, "Serialized shard exceeds the audit size bound.")
  value
}

lsg_portable_parse_serialized_v1 <- function(value) {
  # XDR only; parse structural boundaries before treating any bytes as flags.
  lsg_portable_assert_v1(is.raw(value) && length(value) >= 14L &&
                          identical(value[1:2], charToRaw("X\n")),
                          "Unsupported serialization: an XDR stream is required.")
  position <- 3L
  headers <- vector("list", 20000L)
  header_count <- 0L
  reference_count <- 0L
  consume <- function(count) {
    lsg_portable_assert_v1(length(count) == 1L && is.finite(count) && count >= 0 &&
                            count == floor(count) && position + count - 1 <= length(value),
                            "Truncated or invalid serialized item length.")
    start <- position
    position <<- position + as.integer(count)
    if (!count) raw() else value[seq.int(start, length.out = count)]
  }
  unsigned <- function() {
    sum(as.double(as.integer(consume(4L))) * c(16777216, 65536, 256, 1))
  }
  bounded_length <- function(allow_na = FALSE) {
    count <- unsigned()
    if (allow_na && count == 4294967295) return(-1L)
    lsg_portable_assert_v1(count <= 1000000L, "Unsupported long or oversized serialized vector.")
    as.integer(count)
  }
  format_version <- unsigned()
  writer_positions <- seq.int(position, length.out = 4L)
  writer_version <- unsigned()
  minimum_reader_version <- unsigned()
  lsg_portable_assert_v1(format_version %in% c(2, 3), "Only serialization versions 2 and 3 are supported.")
  encoding <- NULL
  if (format_version == 3) {
    count <- bounded_length()
    lsg_portable_assert_v1(count > 0L && count <= 100L, "Unsupported native-encoding header.")
    encoding <- rawToChar(consume(count))
    lsg_portable_assert_v1(identical(encoding, "UTF-8"),
                            "Only the recorded UTF-8 version-3 encoding is supported.")
  }
  body_start <- position
  node <- function(path, depth = 0L) {
    lsg_portable_assert_v1(depth <= 100L && header_count < length(headers),
                            "Serialized object depth/item count exceeds the audit bound.")
    offset <- position
    flags <- unsigned()
    lsg_portable_assert_v1(flags <= .Machine$integer.max, "Unsupported high-bit serialization flags.")
    type <- bitwAnd(as.integer(flags), 255L)
    attribute <- bitwAnd(as.integer(flags), 512L) != 0L
    tag <- bitwAnd(as.integer(flags), 1024L) != 0L
    header_count <<- header_count + 1L
    current <- header_count
    headers[[current]] <<- list(offset = offset, type = type, flags = flags,
                                length = NA_integer_, path = path)
    set_length <- function(count) headers[[current]]$length <<- count
    if (type == 255L) {
      reference <- bitwShiftR(as.integer(flags), 8L)
      if (reference == 0L) reference <- unsigned()
      lsg_portable_assert_v1(reference > 0L && reference <= reference_count,
                              "Invalid serialization symbol reference.")
      return(invisible(NULL))
    }
    if (type == 254L) return(invisible(NULL))
    if (type == 1L) {
      lsg_portable_assert_v1(!attribute && !tag, "Attributed or tagged symbols are unsupported.")
      node(paste0(path, "/symbol"), depth + 1L)
      reference_count <<- reference_count + 1L
      return(invisible(NULL))
    }
    if (type == 2L) {
      if (attribute) node(paste0(path, "/@attributes"), depth + 1L)
      if (tag) node(paste0(path, "/@tag"), depth + 1L)
      node(paste0(path, "/car"), depth + 1L)
      node(paste0(path, "/cdr"), depth + 1L)
      return(invisible(NULL))
    }
    supported_vectors <- c(9L, 10L, 13L, 14L, 16L, 19L)
    lsg_portable_assert_v1(type %in% supported_vectors && !tag,
                            paste("Unsupported serialized SEXP type/flags:", type))
    count <- bounded_length(allow_na = type == 9L)
    set_length(count)
    if (type == 9L) {
      if (count >= 0L) consume(count)
    } else if (type %in% c(10L, 13L)) {
      consume(4L * count)
    } else if (type == 14L) {
      consume(8L * count)
    } else {
      for (i in seq_len(count)) node(paste0(path, "/", i), depth + 1L)
    }
    if (attribute) node(paste0(path, "/@attributes"), depth + 1L)
    invisible(NULL)
  }
  node("$")
  lsg_portable_assert_v1(position == length(value) + 1L, "Trailing or unparsed serialization bytes remain.")
  frame <- do.call(rbind, lapply(headers[seq_len(header_count)], as.data.frame,
                                 stringsAsFactors = FALSE))
  rownames(frame) <- NULL
  lsg_portable_assert_v1(!anyDuplicated(frame$path), "Serialized structural paths are ambiguous.")
  list(format_version = format_version, writer_version = writer_version,
       minimum_reader_version = minimum_reader_version, encoding = encoding,
       writer_positions = writer_positions, body_start = body_start, headers = frame)
}

lsg_portable_content_hash_v1 <- function(original_raw, object, expected_hash) {
  lsg_portable_assert_v1(is.list(object) &&
                          identical(names(object), c("schema_version", "version", "scientific_signature",
                                                     "key", "task", "runtime_seconds", "created_utc",
                                                     "payload", "content_sha256")) &&
                          identical(object$content_sha256, expected_hash) &&
                          length(expected_hash) == 1L && grepl("^[0-9a-f]{64}$", expected_hash),
                          "Unexpected immutable shard/hash schema.")
  archived <- lsg_portable_parse_serialized_v1(original_raw)
  local_raw <- serialize(object, NULL, version = 3L)
  native <- lsg_portable_parse_serialized_v1(local_raw)
  lsg_portable_assert_v1(archived$format_version == 3L &&
                          identical(archived$encoding, native$encoding) &&
                          identical(archived$minimum_reader_version, native$minimum_reader_version) &&
                          length(original_raw) == length(local_raw) &&
                          identical(archived$headers[, c("offset", "type", "length", "path")],
                                    native$headers[, c("offset", "type", "length", "path")]),
                          "Serialized shard structure/format differs beyond the supported representation.")
  old <- archived$headers
  now <- native$headers
  differing <- which(old$flags != now$flags)
  expected_target_names <- c("method", "group", "group_size", "success", "converged",
                             "finite_coefficients", "boundary", "separation_suspected",
                             "iterations", "maximum_absolute_slope", "message")
  lsg_portable_assert_v1(is.data.frame(object$payload$target_diagnostics) &&
                          identical(names(object$payload$target_diagnostics)[seq_len(11L)], expected_target_names) &&
                          nrow(object$payload$target_diagnostics) == 200L,
                          "The archived Firth diagnostic vector schema is unsupported.")
  allowed_paths <- paste0("$/8/3/", seq_len(11L))
  if (length(differing)) {
    lsg_portable_assert_v1(
      all(old$path[differing] %in% allowed_paths) &&
        all(old$type[differing] %in% c(10L, 13L, 14L, 16L)) &&
        all(old$length[differing] == 200L) &&
        all(bitwXor(as.integer(old$flags[differing]), as.integer(now$flags[differing])) == 131072L) &&
        all(bitwAnd(as.integer(old$flags[differing]), 131072L) == 131072L) &&
        all(bitwAnd(as.integer(now$flags[differing]), 131072L) == 0L),
      "An unsupported internal flag change occurred; no general flag normalization is allowed.")
  }
  restored <- local_raw
  restored[native$writer_positions] <- original_raw[archived$writer_positions]
  for (i in differing) {
    positions <- seq.int(old$offset[i], length.out = 4L)
    restored[positions] <- original_raw[positions]
  }
  lsg_portable_assert_v1(identical(restored, original_raw),
                          "Serialized bytes differ outside the permitted writer/header flags.")
  content <- object[setdiff(names(object), "content_sha256")]
  content_raw <- serialize(content, NULL, version = 2L)
  content_structure <- lsg_portable_parse_serialized_v1(content_raw)
  body <- seq.int(content_structure$body_start, length(content_raw))
  native_hash <- digest::digest(content_raw[body], algo = "sha256", serialize = FALSE)
  lsg_portable_assert_v1(identical(native_hash, lsg_tail_object_hash_v6(content)),
                          "Native digest serialization/header semantics are unsupported.")
  for (i in differing) {
    target <- match(old$path[i], content_structure$headers$path)
    header <- content_structure$headers[target, , drop = FALSE]
    lsg_portable_assert_v1(!is.na(target) && header$type == now$type[i] &&
                            header$length == now$length[i] && header$flags == now$flags[i],
                            "Cannot map a proven flag loss to the version-2 signed content structure.")
    content_raw[seq.int(header$offset, length.out = 4L)] <-
      original_raw[seq.int(old$offset[i], length.out = 4L)]
  }
  recovered <- digest::digest(content_raw[body], algo = "sha256", serialize = FALSE)
  lsg_portable_assert_v1(identical(recovered, expected_hash),
                          "Original archived content SHA-256 did not reproduce exactly.")
  list(passed = TRUE, native_content_hash_matches = identical(native_hash, expected_hash),
       restored_header_count = length(differing), original_content_hash_recovered = TRUE,
       archived_writer_version = archived$writer_version, native_writer_version = native$writer_version,
       reconstructed_content_sha256 = recovered)
}

lsg_portable_match_specification_v1 <- function(fresh, stored) {
  lsg_portable_assert_v1(lsg_tail_specification_valid_v6(stored),
                          "The original V6 metadata self-signature is invalid.")
  excluded <- c("scientific_signature", "initial_execution")
  candidate <- fresh[setdiff(names(fresh), excluded)]
  original <- stored[setdiff(names(stored), excluded)]
  lsg_portable_assert_v1(identical(names(candidate), names(original)),
                          "Fresh and archived metadata field inventories differ.")
  manifests <- c("input_manifest", "code_manifest", "smoke_manifest")
  for (name in manifests) {
    # Row order and R attribute order are representation, not file identity.
    # Reuse the exact archived representation only AFTER filename/size/hash proof.
    lsg_portable_same_manifest_v1(candidate[[name]], original[[name]])
    candidate[[name]] <- original[[name]]
  }
  for (name in setdiff(names(candidate), manifests)) {
    lsg_portable_assert_v1(
      identical(candidate[[name]], original[[name]], num.eq = FALSE,
                single.NA = FALSE, attrib.as.set = FALSE) &&
        identical(lsg_tail_object_hash_v6(candidate[[name]]), lsg_tail_object_hash_v6(original[[name]])),
      paste("An immutable scientific/runtime metadata field changed:", name))
  }
  lsg_portable_assert_v1(identical(lsg_tail_object_hash_v6(candidate), stored$scientific_signature),
                          "The original scientific signature did not reconstruct exactly.")
  invisible(TRUE)
}

lsg_portable_r_version_code_v1 <- function(version) {
  parts <- regmatches(version, regexec("^R version ([0-9]+)\\.([0-9]+)\\.([0-9]+)", version))[[1L]]
  lsg_portable_assert_v1(length(parts) == 4L, "Unsupported recorded R version string.")
  components <- as.integer(parts[-1L])
  lsg_portable_assert_v1(all(!is.na(components)) && all(components >= 0 & components <= 255),
                          "Recorded R version components are out of bounds.")
  sum(components * c(65536, 256, 1))
}

lsg_portable_verify_shard_v1 <- function(path, task, specification) {
  lsg_portable_assert_v1(is.data.frame(task) && nrow(task) == 1L &&
                          identical(basename(path), task$diagnostic_shard_file[[1L]]) && file.exists(path),
                          "Diagnostic shard path does not match its expected task filename.")
  value <- readRDS(path)
  lsg_portable_assert_v1(
    is.list(value) && identical(value$schema_version, "logistic_sglasso_lambda_tail_diagnostic_shard_v6") &&
      identical(value$version, specification$version) &&
      identical(value$scientific_signature, specification$scientific_signature) &&
      identical(value$key, task$key[[1L]]) && identical(value$task, task) &&
      length(value$runtime_seconds) == 1L && is.finite(value$runtime_seconds) && value$runtime_seconds >= 0,
    "Diagnostic shard schema/version/signature/task/runtime invariant failed.")
  proof <- lsg_portable_content_hash_v1(lsg_portable_read_raw_v1(path), value, value$content_sha256)
  lsg_portable_assert_v1(
    proof$archived_writer_version == lsg_portable_r_version_code_v1(specification$runtime$r_version),
    "RDS writer version differs from the recorded diagnostic runtime.")
  cases <- specification$cases[specification$cases$task_id == task$task_id, , drop = FALSE]
  lsg_tail_validate_payload_v6(value$payload, cases, specification$configuration)
  checks <- data.frame(
    task_id = task$task_id, key = task$key, shard_file = basename(path),
    signature_and_task_verified = TRUE, archived_writer_verified = TRUE,
    archived_writer_version = proof$archived_writer_version,
    native_writer_version = proof$native_writer_version,
    native_content_hash_matches = proof$native_content_hash_matches,
    restored_header_count = proof$restored_header_count,
    original_content_hash_recovered = proof$original_content_hash_recovered,
    original_content_sha256 = value$content_sha256,
    reconstructed_content_sha256 = proof$reconstructed_content_sha256,
    numerical_checks_passed = nrow(value$payload$checks),
    curve_points = nrow(value$payload$curves), fixed_cases = nrow(cases),
    firth_group_rows = nrow(value$payload$target_diagnostics), test_used = FALSE,
    stringsAsFactors = FALSE)
  list(passed = TRUE, shard = value, checks = checks)
}

lsg_portable_read_manifest_v1 <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                  colClasses = c(file = "character", bytes = "numeric", sha256 = "character"))
}

lsg_portable_audit_v1 <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  constants <- lsg_tail_constants_v6()
  version <- constants$diagnostic_version
  paths <- lsg_tail_output_paths_v6(root, version)
  required <- c(paths$metadata, paths$completion, paths$output_manifest)
  lsg_portable_assert_v1(all(file.exists(required)), "The complete archived V6 evidence is missing.")
  stored <- readRDS(paths$metadata)
  lsg_portable_assert_v1(identical(stored$version, version) && lsg_tail_specification_valid_v6(stored),
                          "The archived diagnostic version or original metadata self-hash is invalid.")

  # All physical files are checked independently of object reserialization.
  expected_artifacts <- lsg_tail_expected_artifacts_v6(root, stored)
  output_manifest <- lsg_portable_read_manifest_v1(paths$output_manifest)
  lsg_portable_verify_manifest_v1(root, output_manifest, expected_artifacts)
  visible <- list.files(paths$directory, recursive = TRUE, full.names = TRUE)
  lsg_portable_assert_v1(setequal(lsg_tail_relative_v6(visible),
                                  lsg_tail_relative_v6(c(expected_artifacts, paths$output_manifest, paths$completion))),
                          "Unexpected or missing files in the archived diagnostic directory.")
  expected_marker <- c(
    paste("version", stored$version), "stage conditional_lambda_tail_diagnostic",
    paste("scientific_signature", stored$scientific_signature),
    paste("source_scientific_signature", stored$source_scientific_signature),
    paste("source_tuning_signature", stored$source_tuning_signature),
    "completed_tasks 11", "fixed_cases 16", "curve_points 80", "test_used FALSE",
    "calibration_or_production_approved FALSE",
    paste("output_manifest_sha256", digest::digest(file = paths$output_manifest, algo = "sha256")))
  lsg_portable_assert_v1(identical(readLines(paths$completion, warn = FALSE), expected_marker),
                          "The original completion marker does not bind the verified output manifest.")
  for (name in c("input_manifest", "code_manifest", "smoke_manifest")) {
    lsg_portable_verify_manifest_v1(root, stored[[name]])
  }

  # Reconstruct the frozen V5 study and exact boundary-case inventory using the
  # unmodified, hash-verified legacy functions; no solver is sourced or invoked.
  source_study <- lsg_tail_read_source_v6(root)
  code_manifest <- lsg_portable_manifest_v1(root, lsg_tail_required_source_paths_v6(root))
  lsg_portable_assert_v1(nrow(code_manifest) == 39L &&
                          lsg_tail_runtime_matches_v6(stored$runtime, source_study$metadata$configuration),
                          "The original 39-file inventory or recorded fitting runtime differs from V5.")
  smoke <- lsg_tail_verify_smoke_v6(root, code_manifest, source_study$metadata$configuration,
                                    required = TRUE, strict_runtime = FALSE)
  lsg_portable_assert_v1(isTRUE(smoke$present) && isTRUE(smoke$runtime_matches_source),
                          "The archived V6 smoke receipt is not from the frozen V5 runtime.")
  fresh <- lsg_tail_specification_v6(source_study, code_manifest, smoke, stored$runtime, version)
  lsg_portable_match_specification_v1(fresh, stored)
  manifest_representation_changes <- vapply(c("input_manifest", "code_manifest", "smoke_manifest"), function(name) {
    !identical(fresh[[name]], stored[[name]], num.eq = FALSE, single.NA = FALSE, attrib.as.set = FALSE) ||
      !identical(lsg_tail_object_hash_v6(fresh[[name]]), lsg_tail_object_hash_v6(stored[[name]]))
  }, logical(1))

  tasks <- stored$task_grid
  shard_paths <- file.path(paths$shard_directory, tasks$diagnostic_shard_file)
  verified <- lapply(seq_len(nrow(tasks)), function(i) {
    lsg_portable_verify_shard_v1(shard_paths[i], tasks[i, , drop = FALSE], stored)
  })
  shards <- lapply(verified, `[[`, "shard")
  shard_checks <- do.call(rbind, lapply(verified, `[[`, "checks"))
  rownames(shard_checks) <- NULL
  tables <- lsg_tail_final_tables_v6(stored, shards)
  for (name in names(tables)) {
    path <- file.path(paths$directory, paste0(name, "_", version, ".csv"))
    exported <- lsg_tail_read_table_v6(path, tables[[name]])
    lsg_portable_assert_v1(lsg_tail_equal_v6(exported, tables[[name]], tolerance = 1e-12),
                            paste("Exported CSV does not reproduce the original shard values:", name))
  }
  for (pair in list(c("input_manifest", "input_manifest"), c("source_manifest", "code_manifest"),
                    c("smoke_manifest", "smoke_manifest"))) {
    exported <- lsg_portable_read_manifest_v1(file.path(paths$directory, paste0(pair[1L], "_", version, ".csv")))
    lsg_portable_same_manifest_v1(exported, stored[[pair[2L]]])
  }
  case_export <- lsg_tail_read_table_v6(file.path(paths$directory, paste0("case_inventory_", version, ".csv")), stored$cases)
  lsg_portable_assert_v1(lsg_tail_equal_v6(case_export, stored$cases, tolerance = 1e-12),
                          "Exported case inventory differs from the original frozen cases.")
  shard_manifest <- lsg_portable_read_manifest_v1(file.path(paths$directory, paste0("shard_manifest_", version, ".csv")))
  lsg_portable_verify_manifest_v1(root, shard_manifest, shard_paths)

  # Recheck provenance after reading/validating every payload (fail on mid-audit
  # modifications); include the marker and manifest themselves in the report.
  for (name in c("input_manifest", "code_manifest", "smoke_manifest")) {
    lsg_portable_verify_manifest_v1(root, stored[[name]])
  }
  lsg_portable_verify_manifest_v1(root, output_manifest, expected_artifacts)
  marker_evidence <- lsg_portable_manifest_v1(root, c(paths$completion, paths$output_manifest))
  union <- do.call(rbind, list(stored$input_manifest, stored$code_manifest, stored$smoke_manifest,
                                output_manifest, marker_evidence))
  for (file in unique(union$file[duplicated(union$file)])) {
    rows <- union[union$file == file, , drop = FALSE]
    lsg_portable_assert_v1(length(unique(rows$bytes)) == 1L && length(unique(rows$sha256)) == 1L,
                            "Repeated evidence filename has conflicting byte/hash records.")
  }
  union <- union[!duplicated(union$file), , drop = FALSE]
  union <- union[order(union$file, method = "radix"), , drop = FALSE]
  rownames(union) <- NULL
  lsg_portable_verify_manifest_v1(root, union)

  checks <- rbind(data.frame(
    check = c("original_v6_metadata_self_hash", "original_completion_manifest_binding",
              "all_202_v5_input_file_bytes_and_hashes", "all_39_frozen_source_file_bytes_and_hashes",
              "all_3_original_smoke_evidence_files", "exact_v5_scientific_and_tuning_signatures",
              "exact_200_v5_tasks_and_16_boundary_selections", "original_v6_scientific_signature_reproduced",
              "archived_r_and_package_runtime_preserved", "all_11_original_shard_content_hashes_recovered",
              "all_99_original_payload_checks", "all_2200_firth_group_rows", "all_6_exported_tables_reproduced",
              "all_original_output_file_bytes_and_hashes", "all_evidence_rechecked_after_audit"),
    passed = c(TRUE, TRUE, nrow(stored$input_manifest) == 202L, nrow(code_manifest) == 39L,
               nrow(stored$smoke_manifest) == 3L, TRUE, nrow(source_study$metadata$task_grid) == 200L &&
                 nrow(source_study$cases) == 16L, TRUE, TRUE,
               nrow(shard_checks) == 11L && all(shard_checks$original_content_hash_recovered),
               sum(shard_checks$numerical_checks_passed) == 99L,
               sum(shard_checks$firth_group_rows) == 2200L, length(tables) == 6L, TRUE, TRUE),
    stringsAsFactors = FALSE), tables$diagnostic_checks)
  lsg_portable_assert_v1(!anyDuplicated(checks$check) && all(checks$passed),
                          "Portable audit has an incomplete or failed required check.")
  list(checks = checks, shard_checks = shard_checks,
       counts = list(tasks = length(shards), cases = nrow(stored$cases), points = nrow(tables$tail_curves)),
       source_manifest = stored$code_manifest, input_manifest = stored$input_manifest,
       smoke_manifest = stored$smoke_manifest, output_manifest = output_manifest, verified_manifest = union,
       scientific_signature = stored$scientific_signature, summary = tables$tail_summary,
       verification_metadata = list(
         schema_version = "logistic_sglasso_tail_portable_audit_v1",
         verified_version = version, verifier_r_version = R.version.string,
         verifier_platform = R.version$platform, verifier_locale = Sys.getlocale(),
         verifier_rng_kind = RNGkind(), verifier_package_versions = lsg_tail_runtime_v6()$package_versions,
         recorded_fitting_runtime = stored$runtime,
         manifest_representation_changes = manifest_representation_changes,
         representation_policy = "exact_original_hash_recovery_with_proven_manifest_and_vector_header_differences_only",
         fitting_performed = FALSE, outputs_modified = FALSE, scientific_policy_changed = FALSE,
         calibration_or_production_approved = FALSE))
}
