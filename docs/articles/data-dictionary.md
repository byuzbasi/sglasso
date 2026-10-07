# Data dictionary and source terms

## GenAtHum

The package’s Gene Atlas Human example has 158 observations, 2,045
predictors and 79 groups. It is used in the Gaussian examples and is not
the CoRSIVSZ schizophrenia dataset. Load it without a network request:

``` r

library(sglasso)
```

    ## Loading required package: Matrix

``` r

data(GenAtHum, package="sglasso")
c(observations=nrow(GenAtHum$X), predictors=ncol(GenAtHum$X),
  groups=length(unique(GenAtHum$group)))
```

    ## observations   predictors       groups 
    ##          158         2045           79

`X` is the numeric design matrix; `y` is the response; `group` maps
columns to groups. `gene_code`, `gene_name` and `groups_name` give
annotation. See
[`help("GenAtHum")`](https://byuzbasi.github.io/sglasso/reference/GenAtHum.md)
for the package’s source description.

## CoRSIVSZ

CoRSIVSZ is a processed blood-methylation panel derived from public
cohorts, not a newly collected cohort or a random train/test split.
Original measurements are in GEO
[GSE84727](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84727)
(development) and
[GSE80417](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80417)
(independent external test).

| Component | Meaning |
|:---|:---|
| `development$X` | 847 by 1,107 deposited methylation values |
| `development$y` | 414 cases and 433 controls |
| `external$X` | 675 by 1,107 measurements in the same column order |
| `external$y` | 353 cases and 322 controls |
| `group` | One integer group label per column; 409 groups |
| `group_name` | Original annotation cluster labels |
| `probe_id` | CpG identifiers matching matrix columns |
| `sample_id`, `geo_sample_id` | Public cohort/array and GEO sample identifiers |
| `preprocessing`, `provenance` | Construction rules, accessions and source records |

For `y`, zero denotes control and one denotes schizophrenia, derived
from the deposited diagnosis codes. It is not a dichotomized
gene-expression outcome. Case-control proportions do not estimate
population prevalence.

Groups are annotation-defined, not discovered from response performance.
Each retained group has at least two original probes and every member is
measured in both cohorts. Incomplete groups are excluded whole. There is
no extra correlation cutoff, imputation, probe averaging, group merging
or response-association filter. Deposited normalized measurements remain
unchanged. Training-specific transformations belong inside model
fitting/CV.

## Load the separate file

The exact `CoRSIVSZ_v1.rds` asset is prepared locally; a public hosting
endpoint has not yet been claimed. It is not inside the package tarball
or the logistic run-code directory. `data(CoRSIVSZ)` is not supported.

``` r

# Obtain the exact versioned file from its separate distribution notice.
CoRSIVSZ <- load_CoRSIVSZ("CoRSIVSZ_v1.rds")
dim(CoRSIVSZ$development$X)
dim(CoRSIVSZ$external$X)
```

[`load_CoRSIVSZ()`](https://byuzbasi.github.io/sglasso/reference/load_CoRSIVSZ.md)
checks the expected size/SHA-256 before reading the file.
`download_CoRSIVSZ(url, destfile)` requires an explicit HTTPS asset URL
and refuses overwrites. Neither helper fits models or changes
measurements. No data download is performed by this guide.

## Original attribution

Gunasekara CJ, Hannon E, MacKay H, et al. (2021). A machine learning
case-control classifier for schizophrenia based on DNA methylation in
blood. *Translational Psychiatry* 11, 412.
[doi:10.1038/s41398-021-01496-3](https://doi.org/10.1038/s41398-021-01496-3).

The [original annotation
repository](https://github.com/waterlandlab/CoRSIV-Methylation-based-SZ-Risk-Score)
supplies cluster definitions. See the installed `CoRSIVSZ-NOTICE.txt`
for source and license notices. The package’s 409-group analysis panel
is not an exact reproduction of the original paper’s 1,982-region panel.
Public access alone does not establish unrestricted third-party
licensing or institutional ethics exemption. Retain original attribution
when using these resources.
