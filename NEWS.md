# sglasso 1.2.0 — development source

This describes GitHub development source, not a tagged or CRAN release.

## Binomial interface

- `sglasso(..., family = "binomial")` accepts numeric/logical 0/1 responses.
- Training-only groupwise Firth targets and the accelerated RcppArmadillo
  hybrid solver support target-directed quadratic shrinkage.
- Binomial CV recomputes preprocessing and targets inside each training fold,
  selects by log-loss and keeps alpha fixed within a call.
- Response predictions are probabilities of the outcome coded as one.
- Finite-grid endpoint warnings, required-target checks and numerical diagnostics
  are documented. No universal positive-target null-model lambda is claimed.
- Gaussian numerical routines are unchanged by the binomial integration.

## Data helpers

- `CoRSIVSZ` documentation describes independent development and external
  methylation cohorts and annotation-defined predictor groups.
- Explicit loaders/download helpers verify the separately distributed file
  against its recorded size and SHA-256. Nothing is downloaded automatically.
- Public hosting of the release asset remains pending.

## Documentation and research code

- SVG identity, revised README and a `pkgdown` documentation website.
- Small executable guides, a data dictionary, gallery and reproduction guide
  distinguish package usage from article-specific runs.
- `CITATION.cff` supplies software metadata and the published Gaussian article
  as the preferred citation; no DOI is assigned to the unpublished logistic study.
- The published article citation is completed consistently across documentation:
  Yüzbaşı, B. and Cao, J. (2026). **Collinear Groupwise Selection via Scaled Group
  Lasso.** *The American Statistician*, 1–23. Advance online publication,
  17 September 2026.
  [doi:10.1080/00031305.2026.2709494](https://doi.org/10.1080/00031305.2026.2709494).
- Logistic run code is separate from `paper_codes/`, which belongs to the
  article in *The American Statistician*, and contains code/settings only.
- Issue templates request minimal, non-sensitive reproductions rather than data.
- No new production fits, CV study or bootstrap study were run for this
  documentation; only small documentation examples are executed.
