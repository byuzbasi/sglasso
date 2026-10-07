#' Load the CoRSIVSZ methylation dataset
#'
#' Read the separately distributed, versioned CoRSIVSZ RDS file. The file size
#' and SHA-256 digest are checked against package metadata before deserialization.
#' No download, model fitting, recoding or preprocessing is performed.
#'
#' @param file Path to the unmodified \code{CoRSIVSZ_v1.rds} file.
#' @return A list of class \code{CoRSIVSZ}, described in \code{\link{CoRSIVSZ}}.
#' @seealso \code{\link{download_CoRSIVSZ}}, \code{\link{CoRSIVSZ}}
#' @examples
#' # Does not access the network or fit a model.
#' if (file.exists("CoRSIVSZ_v1.rds")) {
#'   CoRSIVSZ <- load_CoRSIVSZ("CoRSIVSZ_v1.rds")
#'   dim(CoRSIVSZ$development$X)
#'   table(CoRSIVSZ$development$y)
#' }
#' @export
load_CoRSIVSZ <- function(file) {
  .corsivsz_read_verified(file, .corsivsz_manifest())
}

#' Download the separately distributed CoRSIVSZ dataset
#'
#' Download only on an explicit call, using an HTTPS URL supplied by the user.
#' The expected version, size and SHA-256 digest are pinned in the package.
#' A temporary file is verified before publishing it at \code{destfile}.
#' Existing files are never overwritten. No model is fitted.
#'
#' @param url An HTTPS URL for the exact \code{CoRSIVSZ_v1.rds} release asset.
#'   Obtain this from the dataset distribution notice. No default hosting
#'   endpoint is assumed, and credentials in URLs are rejected.
#' @param destfile Destination filename. Its parent directory must exist.
#' @param quiet Logical; passed to \code{utils::download.file}.
#' @return Invisibly, the destination path. Load it with
#'   \code{\link{load_CoRSIVSZ}}.
#' @seealso \code{\link{CoRSIVSZ}}, \code{\link{load_CoRSIVSZ}}
#' @examples
#' # Explicit opt-in only; package checks do not download data.
#' \dontrun{
#' # dataset_url must be an HTTPS release-asset URL, not a repository page.
#' download_CoRSIVSZ(dataset_url, "CoRSIVSZ_v1.rds")
#' CoRSIVSZ <- load_CoRSIVSZ("CoRSIVSZ_v1.rds")
#' }
#' @export
download_CoRSIVSZ <- function(url, destfile, quiet = TRUE) {
  .corsivsz_string(url, "url")
  .corsivsz_string(destfile, "destfile")
  if (!grepl("^https://[^/[:space:]]+/.+", url) ||
      grepl("^https://[^/]*@", url) || grepl("[[:space:]]", url))
    stop("url must be an HTTPS release-asset URL without credentials.", call. = FALSE)
  if (!is.logical(quiet) || length(quiet) != 1L || is.na(quiet))
    stop("quiet must be TRUE or FALSE.", call. = FALSE)
  link_target <- Sys.readlink(destfile)
  if (file.exists(destfile) || (!is.na(link_target) && nzchar(link_target)))
    stop("Destination exists; refusing overwrite.", call. = FALSE)
  if (!dir.exists(dirname(destfile)))
    stop("The destination directory must already exist.", call. = FALSE)
  manifest <- .corsivsz_manifest()
  temporary <- tempfile("CoRSIVSZ_download_", tmpdir = dirname(destfile))
  on.exit(unlink(temporary), add = TRUE)
  status <- .corsivsz_download(url, temporary, quiet)
  if (!identical(status, 0L)) stop("Download failed.", call. = FALSE)
  .corsivsz_read_verified(temporary, manifest)
  # A hard link publishes without replacing a path created concurrently.
  if (!file.link(temporary, destfile))
    stop("Cannot publish verified file without overwrite; check filesystem support.",
      call. = FALSE)
  invisible(destfile)
}

.corsivsz_download <- function(url, destfile, quiet) {
  utils::download.file(url, destfile, mode = "wb", quiet = quiet)
}

.corsivsz_string <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x))
    stop(name, " must be one nonempty character string.", call. = FALSE)
}

.corsivsz_manifest <- function() {
  file <- system.file("extdata", "CoRSIVSZ_v1.dcf", package = "sglasso")
  if (!nzchar(file)) stop("CoRSIVSZ release metadata is missing.", call. = FALSE)
  x <- read.dcf(file)
  required <- c("Dataset", "Schema", "File", "Bytes", "SHA256")
  if (nrow(x) != 1L || !all(required %in% colnames(x)) ||
      x[1, "Dataset"] != "CoRSIVSZ" || x[1, "Schema"] != "CoRSIVSZ_v1" ||
      !grepl("^[0-9a-f]{64}$", x[1, "SHA256"]) ||
      !is.finite(as.numeric(x[1, "Bytes"])))
    stop("Invalid CoRSIVSZ release metadata.", call. = FALSE)
  as.list(x[1, required])
}

.corsivsz_read_verified <- function(file, manifest) {
  .corsivsz_string(file, "file")
  if (!file.exists(file) || dir.exists(file))
    stop("CoRSIVSZ file does not exist or is not a regular file.", call. = FALSE)
  if (file.info(file)$size != as.numeric(manifest$Bytes) ||
      !identical(digest::digest(file = file, algo = "sha256"), manifest$SHA256))
    stop("CoRSIVSZ size/SHA-256 mismatch; file is incomplete or a different version.",
      call. = FALSE)
  x <- readRDS(file)
  .corsivsz_validate(x)
  x
}

.corsivsz_validate <- function(x) {
  fail <- function(ok) {
    if (!isTRUE(ok)) stop("Invalid CoRSIVSZ dataset structure or identities.", call. = FALSE)
  }
  fail(is.list(x) && inherits(x, "CoRSIVSZ") &&
    identical(x$schema_version, "CoRSIVSZ_v1"))
  fail(is.integer(x$group) && length(x$group) == 1107L &&
    !anyNA(x$group) && identical(sort(unique(x$group)), seq_len(409L)))
  fail(is.character(x$group_name) && length(x$group_name) == 409L &&
    !anyNA(x$group_name) && !anyDuplicated(x$group_name))
  fail(is.character(x$probe_id) && length(x$probe_id) == 1107L &&
    !anyNA(x$probe_id) && !anyDuplicated(x$probe_id))
  fail(min(tabulate(x$group)) == 2L && max(tabulate(x$group)) == 12L)
  expected <- list(development = c(847L, 414L), external = c(675L, 353L))
  for (name in names(expected)) {
    z <- x[[name]]; n <- expected[[name]][1L]
    fail(is.list(z) && is.matrix(z$X) && is.double(z$X) &&
      identical(dim(z$X), c(n, 1107L)))
    fail(all(is.finite(z$X)) && all(z$X >= 0 & z$X <= 1))
    fail(identical(colnames(z$X), x$probe_id))
    fail(is.integer(z$y) && length(z$y) == n && !anyNA(z$y) &&
      all(z$y %in% 0:1) && sum(z$y) == expected[[name]][2L])
    for (id in c("sample_id", "geo_sample_id"))
      fail(is.character(z[[id]]) && length(z[[id]]) == n &&
        !anyNA(z[[id]]) && !anyDuplicated(z[[id]]))
    fail(identical(rownames(z$X), z$sample_id))
  }
  fail(!anyDuplicated(c(x$development$geo_sample_id, x$external$geo_sample_id)))
  fail(!anyDuplicated(c(x$development$sample_id, x$external$sample_id)))
  invisible(TRUE)
}
