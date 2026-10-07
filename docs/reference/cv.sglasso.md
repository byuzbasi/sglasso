# Cross-validation for scaled group lasso

Performs K-fold cross-validation for the scaled group lasso over a grid
of lambda and d values.

## Usage

``` r
cv.sglasso(
  X,
  Y,
  group = 1:ncol(X),
  lambda,
  nlambda = 100,
  d,
  nd = 11,
  alpha = 0.5,
  nfolds = 10,
  fold,
  beta_start = NULL,
  family = "gaussian",
  bilevel = FALSE,
  max_iter = 1e+08,
  eps = 1e-04,
  standardize = TRUE,
  screen = c("SSR", "none", "SSR_fast"),
  ...
)
```

## Arguments

- X:

  Design matrix.

- Y:

  Response vector.

- group:

  Group membership vector.

- lambda:

  Optional lambda sequence.

- nlambda:

  Number of lambda values.

- d:

  Optional d sequence.

- nd:

  Number of d values.

- alpha:

  Mixing parameter.

- nfolds:

  Number of cross-validation folds.

- fold:

  Optional fold assignment vector.

- beta_start:

  Optional initial coefficient vector.

- family:

  Model family: `"gaussian"` or `"binomial"`.

- bilevel:

  Logical; whether bilevel selection is used.

- max_iter:

  Maximum number of iterations.

- eps:

  Convergence tolerance.

- standardize:

  Logical flag indicating whether `X` should be centered and scaled
  internally before fitting.

- screen:

  Screening rule passed to
  [`sglasso`](https://byuzbasi.github.io/sglasso/reference/sglasso.md).

- ...:

  For binomial fits: `lambda.min.ratio` and `binomial.control`, passed
  to the numerical adapter.

## Value

An object of class `cv.sglasso`.

## Details

Binomial CV minimizes mean out-of-fold log-loss, using raw training rows
to recompute standardization, group orthonormalization and Firth targets
separately in every fold. Automatic paths align by a common relative
grid with fold-specific training-only reference scales; explicit
`lambda` aligns by absolute values. `alpha` is a fixed scalar, not tuned
by this function. Use the same supplied folds for comparisons across
alpha values. The test set must not be supplied to this function. See
[`sglasso-binomial`](https://byuzbasi.github.io/sglasso/reference/sglasso-binomial.md)
for finite-range and numerical limitations.

## See also

[`sglasso`](https://byuzbasi.github.io/sglasso/reference/sglasso.md)

## Examples

``` r
data(Birthwt, package = "grpreg")
X <- Birthwt$X
group <- Birthwt$group
Y <- Birthwt$bwt
CVsglasso_fit <- cv.sglasso(X, Y, group)
```
