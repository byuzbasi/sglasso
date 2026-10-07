# Download the separately distributed CoRSIVSZ dataset

Download only on an explicit call, using an HTTPS URL supplied by the
user. The expected version, size and SHA-256 digest are pinned in the
package. A temporary file is verified before publishing it at
`destfile`. Existing files are never overwritten. No model is fitted.

## Usage

``` r
download_CoRSIVSZ(url, destfile, quiet = TRUE)
```

## Arguments

- url:

  An HTTPS URL for the exact `CoRSIVSZ_v1.rds` release asset. Obtain
  this from the dataset distribution notice. No default hosting endpoint
  is assumed, and credentials in URLs are rejected.

- destfile:

  Destination filename. Its parent directory must exist.

- quiet:

  Logical; passed to
  [`utils::download.file`](https://rdrr.io/r/utils/download.file.html).

## Value

Invisibly, the destination path. Load it with
[`load_CoRSIVSZ`](https://byuzbasi.github.io/sglasso/reference/load_CoRSIVSZ.md).

## See also

[`CoRSIVSZ`](https://byuzbasi.github.io/sglasso/reference/CoRSIVSZ.md),
[`load_CoRSIVSZ`](https://byuzbasi.github.io/sglasso/reference/load_CoRSIVSZ.md)

## Examples

``` r
# Explicit opt-in only; package checks do not download data.
if (FALSE) { # \dontrun{
# dataset_url must be an HTTPS release-asset URL, not a repository page.
download_CoRSIVSZ(dataset_url, "CoRSIVSZ_v1.rds")
CoRSIVSZ <- load_CoRSIVSZ("CoRSIVSZ_v1.rds")
} # }
```
