# Prepare/verify the separate data release. No download, fitting or RNG use.
# From the repository root: Rscript --vanilla data-raw/CoRSIVSZ.R --prepare
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L, args %in% c("--prepare", "--verify"))
stopifnot(file.exists("DESCRIPTION"))
source("R/corsivsz.R")
input <- "logistic_prework/data_processed/corsiv_sz_complete_groups_v1"
out <- "output/data/CoRSIVSZ_v1"
metadata_file <- "inst/extdata/CoRSIVSZ_v1.dcf"
hash <- function(p) digest::digest(file = p, algo = "sha256")
m <- readRDS(file.path(input, "MANIFEST.rds"))
stopifnot(!anyDuplicated(m$file), !any(grepl("(^/|[.][.])", m$file)))
for (i in seq_len(nrow(m))) {
  f <- file.path(input, m$file[i])
  stopifnot(file.exists(f), file.info(f)$size == m$bytes[i], hash(f) == m$sha256[i])
}
original <- readRDS(file.path(input, "data.rds"))
cohort <- function(z) z[c("X", "y", "sample_id", "geo_sample_id")]
CoRSIVSZ <- structure(list(schema_version = "CoRSIVSZ_v1",
  development = cohort(original$train), external = cohort(original$test),
  group = original$group, group_name = original$group_name,
  probe_id = original$probe_id, preprocessing = original$preprocessing,
  provenance = list(development = "GSE84727", external = "GSE80417",
    publication = "https://doi.org/10.1038/s41398-021-01496-3",
    annotation = "https://github.com/waterlandlab/CoRSIV-Methylation-based-SZ-Risk-Score",
    annotation_sha256 = hash(paste0("logistic_prework/data_external/corsiv_sz_v1/",
      "CoRSIV_ESS_SIV_CG_sites_clusters_hg38.csv")),
    processed_source_sha256 = hash(file.path(input, "data.rds")))),
  class = c("CoRSIVSZ", "list"))
.corsivsz_validate(CoRSIVSZ)
stopifnot(identical(CoRSIVSZ$development$X, original$train$X),
  identical(CoRSIVSZ$external$X, original$test$X),
  identical(CoRSIVSZ$group, original$group))
chars <- function(x) c(if (is.character(x)) x else if (is.list(x))
  unlist(lapply(x, chars), use.names = FALSE),
  if (length(attributes(x))) unlist(lapply(attributes(x), chars), use.names = FALSE))
stopifnot(!any(grepl("/Users/|/arf/|codex|bahadir|yuzbasi|byuzbasi|inonu.edu.tr",
  chars(CoRSIVSZ), ignore.case = TRUE)))
filename <- "CoRSIVSZ_v1.rds"
if (args == "--prepare") {
  if (dir.exists(out) || file.exists(metadata_file)) stop("Release exists; use --verify")
  dir.create(out, recursive = TRUE)
  path <- file.path(out, filename)
  saveRDS(CoRSIVSZ, path, version = 2, compress = "xz")
  info <- c(Dataset = "CoRSIVSZ", Schema = "CoRSIVSZ_v1", File = filename,
    Bytes = as.character(file.info(path)$size), SHA256 = hash(path))
  dir.create(dirname(metadata_file), recursive = TRUE, showWarnings = FALSE)
  write.dcf(as.data.frame(as.list(info)), metadata_file, width = 100L)
  stopifnot(file.copy(metadata_file, file.path(out, basename(metadata_file))),
    file.copy("inst/CoRSIVSZ-NOTICE.txt", file.path(out, "NOTICE.txt")))
  writeLines(paste(info[["SHA256"]], filename, sep = "  "),
    file.path(out, paste0(filename, ".sha256")))
  writeLines(c("CoRSIVSZ version 1 - separate data release",
    "",
    "This folder is a prepared local release, not evidence of a public deposit.",
    "An HTTPS release-asset URL will be needed for downloading after publication.",
    "With the package version providing the loader, use:",
    "  CoRSIVSZ <- sglasso::load_CoRSIVSZ(\"CoRSIVSZ_v1.rds\")",
    "",
    "The RDS contains development and external cohorts; their original split is",
    "preserved. y=0 control; y=1 schizophrenia. See help(\"CoRSIVSZ\").",
    "No extra observations, predictors, imputation, scaling or recoding are",
    "introduced by packaging. Model fitting and CV are not part of loading.",
    "Public sample IDs are retained; researcher-local paths are not included.",
    "The exact size and SHA-256 are in CoRSIVSZ_v1.dcf; attribution is in NOTICE.txt."
  ), file.path(out, "README.txt"))
  files <- sort(list.files(out))
  write.csv(data.frame(file = files, bytes = file.info(file.path(out, files))$size,
    sha256 = vapply(file.path(out, files), hash, "")),
    file.path(out, "MANIFEST.csv"), row.names = FALSE)
}
manifest <- as.list(read.dcf(metadata_file)[1L, ])
restored <- .corsivsz_read_verified(file.path(out, filename), manifest)
stopifnot(identical(restored, CoRSIVSZ))
released <- read.csv(file.path(out, "MANIFEST.csv"))
for (i in seq_len(nrow(released))) {
  f <- file.path(out, released$file[i])
  stopifnot(file.info(f)$size == released$bytes[i], hash(f) == released$sha256[i])
}
cat("VERIFIED: CoRSIVSZ, exact values/labels/groups, privacy and checksums.\n",
  "Output: ", out, "\nNo model fits; no public upload.\n", sep = "")
