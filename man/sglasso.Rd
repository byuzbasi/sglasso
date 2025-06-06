% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sglasso.R
\name{sglasso}
\alias{sglasso}
\title{Fit a scaled group lasso regression path}
\usage{
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
  bilevel = F,
  max_iter = 1e+08,
  eps = 1e-04,
  lambda.min = 0.005,
  dfmax = p,
  gmax = length(unique(group))
)
}
\arguments{
\item{X}{The design matrix, without an intercept.  \code{sglasso}
standardizes the data and includes an intercept by default.}

\item{group}{A vector describing the grouping of the coefficients.}

\item{lambda}{A user supplied sequence of \code{lambda} values.  Typically,
this is left unspecified, and the function automatically computes a grid of
lambda values that ranges uniformly on the log scale over the relevant range
of lambda values.}

\item{nlambda}{The number of \code{lambda} values.  Default is 50}

\item{d}{The scale parameter between 0 and 1.}

\item{nd}{The number of \code{d} values.  Default is 11.}

\item{alpha}{Elastic Net tuning constant: the value must be between 0 and 1. Default is 0.5.}

\item{family}{is "gaussian", not other option at this moment, depending on the response.}

\item{bilevel}{bi-level selection is not supported at this moment.}

\item{eps}{Convergence threshhold.  The algorithm iterates until the BCD
for the change in linear predictors for each coefficient is less than
\code{eps}.  Default is \code{1e-4}.}

\item{lambda.min}{The smallest value for \code{lambda}, as a fraction of
\code{lambda.max}.  Default is .0005.}

\item{dfmax}{Limit on the number of parameters allowed to be nonzero.  If
this limit is exceeded, the algorithm will exit early from the
regularization path.}

\item{gmax}{Limit on the number of groups allowed to have nonzero elements.
If this limit is exceeded, the algorithm will exit early from the
regularization path.}

\item{y}{The response vector}

\item{max.iter}{Maximum number of iterations (total across entire path).
Default is 1e+08.}
}
\value{
An object with S3 class \code{"sglasso"} containing:
\describe{
\item{beta}{The fitted matrix of coefficients.  The number of rows is equal
to the number of coefficients, and the number of columns is equal to
\code{nlambda}.}
\item{family}{Same as above.}
\item{group}{Same as above.}
\item{lambda}{The sequence of \code{lambda} values in the path.}
\item{alpha}{Same as above.}
\item{deviance}{A vector containing the deviance of the fitted model at each value of `lambda`.}
\item{n}{Number of observations.}
\item{penalty}{Same as above.}
\item{df}{A vector of length `nlambda` containing estimates of effective number of model parameters all the points along the regularization path.  For details on how this is calculated, see Breheny and Huang (2009).}
\item{iter}{A vector of length `nlambda` containing the number of iterations until convergence at each value of `lambda`.}
\item{group.multiplier}{A named vector containing the multiplicative constant applied to each group's penalty.}
}
}
\description{
Computes a scaled group lasso regularized linear models.
}
\examples{
data(GenAtHum,package = "sglasso")
X <- GenAtHum$X
y <- GenAtHum$y
group <- GenAtHum$group
n =  nrow(X)
p =  ncol(X)
set.seed(123)
fit <- sglasso(X,y,group, nlambda = 20, nd = 3)
select(fit,"EBIC")
}
\seealso{
\code{\link{cv.sglasso}}
}
\author{
Bahadir Yuzbasi and Jiguo Cao
}
