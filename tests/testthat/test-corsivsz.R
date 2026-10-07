corsiv_fixture <- function() {
  group <- rep(seq_len(409L), c(12L, 12L, rep(3L, 269L), rep(2L, 138L)))
  probe <- paste0("cg", seq_len(1107L))
  part <- function(n, cases, prefix) {
    id <- paste0(prefix, seq_len(n))
    list(X = matrix(0.4, n, 1107L, dimnames = list(id, probe)),
      y = c(rep(1L, cases), rep(0L, n - cases)), sample_id = id,
      geo_sample_id = paste0("GSM", id))
  }
  structure(list(schema_version = "CoRSIVSZ_v1",
    group = as.integer(group), group_name = paste0("g", seq_len(409L)),
    probe_id = probe, development = part(847L, 414L, "D"),
    external = part(675L, 353L, "E")), class = c("CoRSIVSZ", "list"))
}

corsiv_test_file <- function(x, file) {
  saveRDS(x, file)
  list(Bytes = as.character(file.info(file)$size),
    SHA256 = digest::digest(file = file, algo = "sha256"))
}

test_that("CoRSIVSZ metadata is small and pins the external release", {
  m <- sglasso:::.corsivsz_manifest()
  expect_identical(m$Schema, "CoRSIVSZ_v1")
  expect_identical(m$File, "CoRSIVSZ_v1.rds")
  expect_match(m$SHA256, "^[a-f0-9]{64}$")
  expect_gt(as.numeric(m$Bytes), 5e6)
  expect_identical(system.file("data", m$File, package = "sglasso"), "")
})

test_that("CoRSIVSZ structure and identities are validated", {
  x <- corsiv_fixture()
  expect_silent(sglasso:::.corsivsz_validate(x))
  bad <- x; bad$external$y[1] <- 0L
  expect_error(sglasso:::.corsivsz_validate(bad), "structure")
  bad <- x; colnames(bad$development$X)[1] <- "wrong"
  expect_error(sglasso:::.corsivsz_validate(bad), "identities")
  bad <- x; bad$group[1] <- 410L
  expect_error(sglasso:::.corsivsz_validate(bad), "structure")
  bad <- x; bad$development$X[1] <- NA_real_
  expect_error(sglasso:::.corsivsz_validate(bad), "structure")
  bad <- x; bad$external$geo_sample_id[1] <- x$development$geo_sample_id[1]
  expect_error(sglasso:::.corsivsz_validate(bad), "identities")
})

test_that("checksum precedes deserialization and values survive round trip", {
  f <- tempfile(fileext = ".rds"); on.exit(unlink(f))
  x <- corsiv_fixture(); m <- corsiv_test_file(x, f)
  expect_identical(sglasso:::.corsivsz_read_verified(f, m), x)
  m$SHA256 <- paste(rep("0", 64), collapse = "")
  expect_error(sglasso:::.corsivsz_read_verified(f, m), "SHA-256 mismatch")
  writeLines("not an RDS file", f)
  expect_error(sglasso:::.corsivsz_read_verified(f, m), "SHA-256 mismatch")
  expect_error(load_CoRSIVSZ(f), "SHA-256 mismatch")
  expect_error(load_CoRSIVSZ(character()), "nonempty")
  expect_error(load_CoRSIVSZ(paste0(f, "missing")), "does not exist")
})

test_that("downloads reject unsafe inputs without network access", {
  d <- tempfile(); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  dest <- file.path(d, "data.rds")
  testthat::local_mocked_bindings(.corsivsz_download = function(...) stop("NETWORK"),
    .package = "sglasso")
  expect_error(download_CoRSIVSZ("http://example.org/file.rds", dest), "HTTPS")
  expect_error(download_CoRSIVSZ("https://user:pass@example.org/f", dest), "credentials")
  expect_error(download_CoRSIVSZ("https://example.org/f", dest, NA), "quiet")
  expect_error(download_CoRSIVSZ("https://example.org/f", file.path(d, "missing/f")),
    "directory")
  writeLines("preserve", dest)
  expect_error(download_CoRSIVSZ("https://example.org/f", dest), "refusing overwrite")
  expect_identical(readLines(dest), "preserve")
})

test_that("verified downloads publish once and preserve an existing destination", {
  d <- tempfile(); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  input <- file.path(d, "input.rds"); dest <- file.path(d, "release.rds")
  x <- corsiv_fixture(); m <- corsiv_test_file(x, input)
  testthat::local_mocked_bindings(.corsivsz_manifest = function() m,
    .corsivsz_download = function(url, destfile, quiet) {
      stopifnot(file.copy(input, destfile, overwrite = FALSE)); 0L
    }, .package = "sglasso")
  expect_identical(download_CoRSIVSZ("https://example.org/f", dest), dest)
  expect_identical(load_CoRSIVSZ(dest), x)
  expect_error(download_CoRSIVSZ("https://example.org/f", dest), "refusing overwrite")
  expect_false(any(startsWith(list.files(d), "CoRSIVSZ_download_")))
})

test_that("failed or corrupted downloads do not publish partial data", {
  d <- tempfile(); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  dest <- file.path(d, "release.rds")
  testthat::local_mocked_bindings(.corsivsz_download = function(url, destfile, quiet) {
    writeLines("incomplete", destfile); 0L
  }, .package = "sglasso")
  expect_error(download_CoRSIVSZ("https://example.org/f", dest), "SHA-256 mismatch")
  expect_false(file.exists(dest))
  expect_length(list.files(d), 0L)
})
