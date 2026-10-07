# Binary responses, targets and cross-validation

## What the target changes

Let the training sample contain $`n`$ observations. Observation $`i`$
has binary response $`y_i\in\{0,1\}`$. There are $`G`$ prespecified
groups. After training-only transformations, $`z_{ig}`$ is the retained
predictor vector for group $`g`$, and $`b_g`$ is its coefficient vector
in the same coordinates. The unpenalized intercept is $`b_0`$, and
$`\eta_i=b_0+\sum_{g=1}^{G}z_{ig}^{\mathsf T}b_g`$ is the linear
predictor.

Conditionally on the training-derived groupwise Firth targets $`t_g`$,
Equation (1) defines the binomial estimation objective: mean negative
Bernoulli log-likelihood plus a mixed group-sparsity/target-shrinkage
penalty.

``` math
 \frac{1}{n}\sum_{i=1}^n\{\log(1+e^{\eta_i})-y_i\eta_i\}
 +\lambda\sum_{g=1}^G w_g
 \left\{\alpha\lVert b_g\rVert_2+
 \frac{1-\alpha}{2}\lVert b_g-dt_g\rVert_2^2\right\}. \tag{1}
```

In Equation (1), $`\lambda>0`$ is finite, $`\alpha,d\in[0,1]`$,
$`\lVert\cdot\rVert_2`$ is the Euclidean norm, and $`w_g`$ is the square
root of the group’s retained rank. Targets and transformations are
estimated from the current training partition, not validation or test
responses. The targets are training-data dependent; this conditional
formulation does not assert an unconditional risk guarantee.

At $`d=0`$, Equation (1) reduces to zero-target group elastic net. At
$`\alpha=1`$, the quadratic term vanishes and $`d`$ has no effect.
Positive $`d`$ shifts the shrinkage destination; it does not guarantee
better prediction or exact group-support recovery. The intercept is not
part of the penalty. See the package’s `sglasso-binomial` reference and
the [logistic research
code](https://github.com/byuzbasi/sglasso/tree/main/logistic_paper_codes)
for the study-specific implementation.

## A small fit

``` r

set.seed(19)
x <- matrix(rnorm(240), 60, 4)
y <- rbinom(60, 1, plogis(x[, 1]))
group <- c(1, 1, 2, 2)
fit <- sglasso(x, y, group, family="binomial",
               lambda=c(0.3, 0.1), d=c(0, 0.5))
predict(fit, newx=x[1:3, ], lambda=0.1, d=0.5, type="response")
#> [1] 0.2800003 0.5824570 0.5030327
```

These are probabilities for `y = 1`, not automatically fitted class
labels. For classification, declare a threshold using
training/validation information or an application-specific decision
policy. Never optimize it on the test set. The first three training rows
above illustrate syntax only.

## Training-only targets in CV

Binomial CV recomputes standardization, group orthonormalization and
Firth targets inside each training fold. It minimizes
observation-weighted out-of-fold log-loss. Alpha remains fixed within a
call.

``` r

set.seed(29)
cv <- cv.sglasso(x, y, group, family="binomial",
                 lambda=c(0.3, 0.1), d=c(0, 0.5),
                 alpha=0.5, nfolds=3)
#> Warning: CV selected a finite lambda-grid endpoint; performance beyond the grid
#> is unassessed.
probability <- predict(cv, newx=x[1:3, ], s="opt", type="response")
probability
#> [1] 0.2800003 0.5824570 0.5030327
```

The tiny grid is for documentation testing, not a recommendation for a
study. Endpoint warnings are informative: a better value beyond the
requested grid has not been ruled out. A supplied lambda vector is
absolute; automatic paths are training-fold-specific relative paths
aligned by grid position. No universal positive-target null-model lambda
maximum is claimed.

## Numerical scope

The compiled solver uses IRLS/proximal-Newton updates, profiled block
updates and monotone accelerated proximal-gradient fallback. Required
Firth targets must pass their adjusted-score/Fisher-step checks; they
are not replaced by ordinary MLE or ridge on failure. Requested finite
path points must pass their numerical gates. No guarantee covers every
possible input.

Missing or nonnumeric predictors are rejected, not imputed. Constant
columns and rank-deficient groups are handled by the documented
training-only transformations. Gaussian screening and stopping options
are not simply transplanted into the binomial interface. Read the
reference before changing `eps`, `max_iter` or `binomial.control`.

The across-fold `cvse` is not an independent-test confidence interval.
ROC-AUC measures ranking; log-loss and Brier score assess probabilities;
MCC assesses classification under a declared threshold. Good performance
on one is not automatically good performance on another.
