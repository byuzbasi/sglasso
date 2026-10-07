# CoRSIVSZ: grouped methylation measurements for schizophrenia classification

A processed analysis panel derived from public blood DNA methylation
cohorts GSE84727 and GSE80417 and the published CoRSIV/ESS/SIV cluster
annotation. `CoRSIVSZ` names this derived panel, not a newly collected
dataset or an exact reproduction of the original publication's
1,982-region panel. The full RDS file is distributed separately, not
under package `data/`. Use
[`load_CoRSIVSZ`](https://byuzbasi.github.io/sglasso/reference/load_CoRSIVSZ.md);
`data(CoRSIVSZ)` is not supported.

## Format

A list of class `CoRSIVSZ` with:

- schema_version:

  The string `CoRSIVSZ_v1`.

- development:

  A list with `X` (847 by 1,107 matrix), `y` (integer response),
  `sample_id` (public array identifiers) and `geo_sample_id` (GSM
  accessions). There are 414 cases and 433 controls.

- external:

  The same fields for 675 independent-test participants: 353 cases and
  322 controls.

- group:

  Integer group membership for the 1,107 columns, from 1 to 409.

- group_name:

  The 409 original annotation cluster labels.

- probe_id:

  The 1,107 CpG identifiers, in matrix-column order.

- preprocessing:

  An explicit description of panel construction and the unchanged
  deposited measurements.

- provenance:

  Public accessions, source links and checksums.

## Source

<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84727>

<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80417>

Annotation:
<https://github.com/waterlandlab/CoRSIV-Methylation-based-SZ-Risk-Score>

## Details

The outcome is the deposited diagnosis: `0` denotes a control and `1` a
schizophrenia case, corresponding to original GEO codes 1 and 2. It is
not a dichotomized gene-expression outcome. The independent cohorts are
not a random split, and their case-control proportions are not
population prevalence estimates.

Of 2,409 annotation-defined groups, retain only groups having at least
two original probes and every member measured in both cohorts.
Incomplete groups are excluded whole. This produces 409 groups with
sizes 2–12. No participants are removed. Columns are ordered by
radix-sorted group names, then probe IDs. No response association,
performance or correlation cutoff selects groups. There is no
imputation, probe averaging, merging or extra smoking-probe filter.

Deposited cohort-level pfilter/dasen-normalized methylation values are
unchanged. No additional scaling or covariate adjustment is applied
here. Model-specific transformations and targets must be estimated only
within training partitions. Public IDs establish common measurement
availability, not external-outcome-driven feature selection.

Package metadata pins the exact release file size and SHA-256 digest.
The complete file is approximately 11 MiB after compression and is
excluded from the package tarball. Package examples never download data
automatically. See the installed `CoRSIVSZ-NOTICE.txt` for source
attribution and terms.

## References

Gunasekara CJ, Hannon E, MacKay H, et al. (2021). A machine learning
case-control classifier for schizophrenia based on DNA methylation in
blood. Translational Psychiatry, 11, 412.
[doi:10.1038/s41398-021-01496-3](https://doi.org/10.1038/s41398-021-01496-3)
.

## See also

[`load_CoRSIVSZ`](https://byuzbasi.github.io/sglasso/reference/load_CoRSIVSZ.md),
[`download_CoRSIVSZ`](https://byuzbasi.github.io/sglasso/reference/download_CoRSIVSZ.md)
