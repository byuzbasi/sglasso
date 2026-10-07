# Select tuning parameters for sglasso

Select tuning parameters for sglasso

## Usage

``` r
# S3 method for class 'sglasso'
select(
  obj,
  criterion = c("AIC", "BIC", "EBIC", "GCV", "AICc"),
  ebic_gamma = 0.5,
  ebic_level = c("group", "feature"),
  tol = 1e-08,
  ...
)
```

## Arguments

- obj:

  A fitted `sglasso` object.

- criterion:

  Selection criterion.

- ebic_gamma:

  EBIC gamma parameter. Default is 0.5.

- ebic_level:

  EBIC penalty level. Use `"group"` for group-level EBIC or `"feature"`
  for feature-level EBIC.

- tol:

  Tolerance used to determine nonzero coefficients.

- ...:

  Additional arguments.
