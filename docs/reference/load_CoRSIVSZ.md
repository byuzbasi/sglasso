# Load the CoRSIVSZ methylation dataset

Read the separately distributed, versioned CoRSIVSZ RDS file. The file
size and SHA-256 digest are checked against package metadata before
deserialization. No download, model fitting, recoding or preprocessing
is performed.

## Usage

``` r
load_CoRSIVSZ(file)
```

## Arguments

- file:

  Path to the unmodified `CoRSIVSZ_v1.rds` file.

## Value

A list of class `CoRSIVSZ`, described in
[`CoRSIVSZ`](https://byuzbasi.github.io/sglasso/reference/CoRSIVSZ.md).

## See also

[`download_CoRSIVSZ`](https://byuzbasi.github.io/sglasso/reference/download_CoRSIVSZ.md),
[`CoRSIVSZ`](https://byuzbasi.github.io/sglasso/reference/CoRSIVSZ.md)

## Examples

``` r
# Does not access the network or fit a model.
if (file.exists("CoRSIVSZ_v1.rds")) {
  CoRSIVSZ <- load_CoRSIVSZ("CoRSIVSZ_v1.rds")
  dim(CoRSIVSZ$development$X)
  table(CoRSIVSZ$development$y)
}
```
