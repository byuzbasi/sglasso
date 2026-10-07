# Start with sglasso

## From grouped predictors to predictions

`sglasso` fits penalized regression paths for user-defined,
nonoverlapping predictor groups. The grouping is an input, not an
estimated clustering. Gaussian fits predict continuous responses;
binomial fits predict the probability of the outcome coded as one.

The example below is deliberately small and synthetic. It demonstrates
the interface, not comparative performance on a research dataset.

``` r

set.seed(19)
x <- matrix(rnorm(80 * 6), 80, 6)
group <- rep(1:3, each=2)
y <- x[, 1] - x[, 3] + rnorm(80)
fit <- sglasso(x, y, group, family="gaussian",
               lambda=c(0.3, 0.1, 0.03), d=c(0, 0.5),
               alpha=0.5, screen="none")
coef(fit, lambda=0.1, d=0.5)
#>    intercept           X1           X2           X3           X4           X5 
#> -0.005216301  0.830825784 -0.070756932 -0.946693824  0.010329021  0.000000000 
#>           X6 
#>  0.000000000
predict(fit, newx=x[1:3, ], lambda=0.1, d=0.5)
#> , , 1
#> 
#>            [,1]
#> [1,] -2.4620222
#> [2,]  0.4347928
#> [3,] -1.3061643
```

The first three training rows illustrate the prediction call. These
numbers are not out-of-sample performance estimates. For a real
application, keep test observations outside model fitting and tuning.

## Cross-validation

Use the documented `nfolds` argument and fix a seed. Alpha is fixed
within each call; use identical folds when comparing different alpha
values.

``` r

set.seed(29)
cv <- cv.sglasso(x, y, group, family="gaussian",
                 lambda=c(0.3, 0.1, 0.03), d=c(0, 0.5),
                 alpha=0.5, nfolds=3, screen="none")
predict(cv, newx=x[1:3, ], s="opt")
#> [1] -2.5172399  0.4950107 -1.4247301
```

This three-fold example is a documentation smoke test, not the protocol
of either article. See the study-specific run code for the original
grids, seeds, validation checks and test designs.

## Read next

- [Binary responses and
  CV](https://byuzbasi.github.io/sglasso/articles/binomial-and-cv.md):
  Firth targets, fold-local preprocessing, probability prediction and
  numerical limits.
- [Data
  dictionary](https://byuzbasi.github.io/sglasso/articles/data-dictionary.md):
  GenAtHum and CoRSIVSZ.
- [Example
  gallery](https://byuzbasi.github.io/sglasso/articles/gallery.md):
  illustrations and completed external results.
- [Reproducing the
  studies](https://byuzbasi.github.io/sglasso/articles/reproducibility.md):
  package API versus run code.

This website documents GitHub development source 1.2.0, not a tagged or
CRAN release. See the repository for the current source and the
study-specific code.
