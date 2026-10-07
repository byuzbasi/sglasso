# Fit a scaled group lasso regression path

Computes scaled group lasso paths for Gaussian or binary responses.

## Usage

``` r
sglasso(
  X,
  Y,
  group = 1:ncol(X),
  lambda,
  nlambda = 50,
  d,
  nd = 11,
  alpha = 0.5,
  beta_start = NULL,
  family = "gaussian",
  bilevel = FALSE,
  max_iter = 1e+08,
  eps = 1e-04,
  standardize = TRUE,
  screen = c("SSR", "none", "SSR_fast"),
  diagnostics = FALSE,
  profile = FALSE,
  transform = c("eager", "lazy"),
  lambda.min.ratio = 0.005,
  dfmax = p,
  gmax = length(unique(group)),
  binomial.control = NULL
)
```

## Arguments

- X:

  The design matrix, without an intercept. `sglasso` standardizes the
  data and includes an intercept by default.

- Y:

  The response vector.

- group:

  A vector describing the grouping of the coefficients.

- lambda:

  A user supplied sequence of `lambda` values. Typically, this is left
  unspecified, and the function automatically computes a grid of lambda
  values that ranges uniformly on the log scale over the relevant range
  of lambda values.

- nlambda:

  The number of `lambda` values. Default is 50

- d:

  The scale parameter between 0 and 1.

- nd:

  The number of `d` values. Default is 11.

- alpha:

  Elastic Net tuning constant: the value must be between 0 and 1.
  Default is 0.5.

- beta_start:

  Optional initial coefficient vector.

- family:

  Either `"gaussian"` (default) or `"binomial"`. Binomial responses must
  be numeric/logical 0/1 with both classes present.

- bilevel:

  bi-level selection is not supported at this moment.

- max_iter:

  Maximum number of iterations Default is 1e+08.

- eps:

  Convergence threshhold. The algorithm iterates until the BCD for the
  change in linear predictors for each coefficient is less than `eps`.
  Default is `1e-4`.

- standardize:

  Logical flag indicating whether `X` should be centered and scaled
  internally before fitting. The default is `TRUE`, preserving the usual
  package behavior. Set to `FALSE` only when `X` has already been
  centered/scaled using the intended training-data transformation.

- screen:

  Screening rule used before the final KKT checks. One of `"SSR"`,
  `"none"`, or `"SSR_fast"`. The `"SSR_fast"` mode checks the rest-set
  KKT conditions at the first, first-quartile, median, third-quartile,
  and final lambda values. The default is `"SSR"`.

- diagnostics:

  Logical flag indicating whether detailed per-lambda timing and final
  inactive-set KKT residual diagnostics should be computed. The default
  is `FALSE` to avoid diagnostic overhead in ordinary fits.

- profile:

  Logical flag indicating whether coarse R-level runtime components
  should be recorded. The default is `FALSE`.

- transform:

  Coefficient transformation mode. The default `"eager"` returns the
  fitted coefficients on the original data scale during `sglasso()`. The
  experimental `"lazy"` mode stores the solver-scale path and performs
  the back-transformation when
  [`coef()`](https://rdrr.io/r/stats/coef.html) is called, which can
  reduce fit-time overhead in runtime comparisons.

- lambda.min.ratio:

  The smallest value for `lambda`, as a fraction of the starting
  reference. Default is .005. For binomial models this is a finite
  search range, not a guarantee of a null model; see
  [`sglasso-binomial`](https://byuzbasi.github.io/sglasso/reference/sglasso-binomial.md).

- dfmax:

  Limit on the number of parameters allowed to be nonzero. If this limit
  is exceeded, the algorithm will exit early from the regularization
  path.

- gmax:

  Limit on the number of groups allowed to have nonzero elements. If
  this limit is exceeded, the algorithm will exit early from the
  regularization path.

- binomial.control:

  Named list of numerical controls for the binomial solver; see
  [`sglasso-binomial`](https://byuzbasi.github.io/sglasso/reference/sglasso-binomial.md).
  Ignored for Gaussian fits.

## Value

An object with S3 class `"sglasso"` containing:

- beta:

  The fitted matrix of coefficients. The number of rows is equal to the
  number of coefficients, and the number of columns is equal to
  `nlambda`.

- family:

  Same as above.

- group:

  Same as above.

- lambda:

  The sequence of `lambda` values in the path.

- alpha:

  Same as above.

- deviance:

  A vector containing the deviance of the fitted model at each value of
  \`lambda\`.

- n:

  Number of observations.

- penalty:

  Same as above.

- df:

  A vector of length \`nlambda\` containing estimates of effective
  number of model parameters all the points along the regularization
  path. For details on how this is calculated, see Breheny and Huang
  (2009).

- iter:

  A vector of length \`nlambda\` containing the number of iterations
  until convergence at each value of \`lambda\`.

- group.multiplier:

  A named vector containing the multiplicative constant applied to each
  group's penalty.

## See also

[`cv.sglasso`](https://byuzbasi.github.io/sglasso/reference/cv.sglasso.md)

## Author

Bahadir Yuzbasi and Jiguo Cao

## Examples

``` r
data(GenAtHum,package = "sglasso")
X <- GenAtHum$X
y <- GenAtHum$y
group <- GenAtHum$group
n =  nrow(X)
p =  ncol(X)
set.seed(123)
fit <- sglasso(X = X, Y = y, group = group, nlambda = 20, nd = 3)
select(fit,"EBIC")
```
