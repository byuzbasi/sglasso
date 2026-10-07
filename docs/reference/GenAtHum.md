# Gene Atlas Human (GenAtHum) Data

The microarray data were obtained from the GDS DataSet of the Gene
Expression Omnibus (GEO) repository in the NCBI archives.

## Format

A list of length 6 containing:

- X:

  Design matrix of dimension 158 x 2045.

- y:

  Response vector of length 158.

- group:

  Group membership code for each predictor.

- gene_code:

  Vector of gene codes.

- gene_name:

  Vector of gene names.

- groups_name:

  Vector of group names.

## Source

<https://www.ncbi.nlm.nih.gov/geo/>

## Details

This data set contains 158 samples, 2045 predictors, and 79 groups, as
described in Yuzbasi and Cao (2025).

## Examples

``` r
data(GenAtHum)
X <- GenAtHum$X
y <- GenAtHum$y
group <- GenAtHum$group
```
