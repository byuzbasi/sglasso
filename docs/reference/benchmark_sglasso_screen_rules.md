# Benchmark SGLASSO Screening Rules

Runs a focused benchmark comparing `"none"`, `"SSR"`, and `"SSR_fast"`
screening rules on block-correlated synthetic designs.

## Usage

``` r
benchmark_sglasso_screen_rules(
  scenarios = data.frame(n = c(200L, 300L), p = c(20000L, 50000L), J = c(2000L, 5000L)),
  nlambda = 100,
  active_groups = 20,
  rho_within = 0.7,
  rho_between = 0.3,
  alpha = 0.5,
  d = 0.5,
  screens = c("none", "SSR", "SSR_fast"),
  reps = 1,
  seed = 2026,
  standardize = TRUE,
  eps = 1e-04,
  verbose = TRUE
)
```

## Arguments

- scenarios:

  A data frame with columns `n`, `p`, and `J`.

- nlambda:

  Number of lambda values.

- active_groups:

  Number of truly active groups.

- rho_within:

  Within-group correlation.

- rho_between:

  Between-group correlation.

- alpha:

  Elastic net mixing parameter passed to
  [`sglasso`](https://byuzbasi.github.io/sglasso/reference/sglasso.md).

- d:

  Scale parameter passed to
  [`sglasso`](https://byuzbasi.github.io/sglasso/reference/sglasso.md).

- screens:

  Screening rules to compare.

- reps:

  Number of repetitions per scenario.

- seed:

  Random seed.

- standardize:

  Whether
  [`sglasso`](https://byuzbasi.github.io/sglasso/reference/sglasso.md)
  should standardize internally.

- eps:

  Convergence tolerance.

- verbose:

  Print progress messages.

## Value

A list containing raw run summaries, per-screen timing summaries, and
pairwise accuracy comparisons against `screen = "none"`.
