# Predict method for sglasso objects

Predict method for sglasso objects

## Usage

``` r
# S3 method for class 'sglasso'
predict(
  object,
  newx = NULL,
  type = c("response", "coefficients", "vars", "groups", "link"),
  lambda,
  d,
  s = c("ALL", "opt"),
  opt_beta = NULL,
  which = 1:length(object$lambda),
  drop = TRUE,
  ...
)
```

## Arguments

- object:

  Fitted sglasso object.

- newx:

  New design matrix.

- type:

  Prediction type. For binomial objects, `"response"` returns
  probabilities of Y=1 and `"link"` returns log odds. For Gaussian
  objects both return the linear predictor. No classification threshold
  is fitted.

- lambda:

  Lambda values.

- d:

  Scaling parameter values.

- s:

  Selection mode.

- opt_beta:

  Optional coefficient vector.

- which:

  Indices of lambda values.

- drop:

  Logical; should dimensions be dropped?

- ...:

  Additional arguments.
