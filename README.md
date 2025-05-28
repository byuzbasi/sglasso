
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sglasso

<!-- badges: start -->
<!-- badges: end -->

This package introduces a group-wise variable selection method, called
Scaled Group Lasso (SGLASSO), to consider correlations among groups.
Incorporating group-wise structure into the regression model improves
interpretability, reduces overfitting, and clarifies the relationships
between different variable groups and the outcome variable using group
selection methods.

## Installation

You can install the development version of sglasso like so:

``` r
# install from GitHub
install.packages("devtools") # if you have not installed "devtools" package
devtools::install_github("byuzbasi/sglasso")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(sglasso)
#> Loading required package: Matrix
data(GenAtHum,package = "sglasso")
X <- GenAtHum$X
y <- GenAtHum$y
group <- GenAtHum$group
n =  nrow(X)
p =  ncol(X)
set.seed(123)
model_CV <- cv.sglasso(X,y,group, nlambda=20, nd=3, nfolds =5)
plot(model_CV)
```

<img src="man/figures/README-example-1.png" width="100%" />
