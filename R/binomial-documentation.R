#' Logistic SGLASSO for binary responses
#'
#' @name sglasso-binomial
#' @description The binomial API uses the accelerated IRLS/proximal-Newton
#' solver with profiled block updates and monotone accelerated proximal-gradient
#' fallback. The compiled core is a snapshot of the validated research solver,
#' not the older Gaussian coordinate-descent implementation.
#' @details On training-standardized, groupwise orthonormalized predictors,
#' the objective is mean negative Bernoulli log-likelihood plus
#' \deqn{\lambda\sum_g w_g\{\alpha\|b_g\|_2 +
#' (1-\alpha)\|b_g-d t_g\|_2^2/2\}.}
#' Here \eqn{b_g} is the solver-coordinate coefficient vector for group
#' \eqn{g}, \eqn{w_g} is the square root of its retained rank, \eqn{t_g}
#' is its training-only Firth slope estimate, \eqn{d\in[0,1]} scales the
#' target, and \eqn{\alpha\in[0,1]} mixes group sparsity and quadratic shrinkage.
#' The intercept is unpenalized. At \eqn{d=0} this is group elastic net; at
#' \eqn{\alpha=1}, the target has no effect. Firth fitting is skipped whenever
#' the target cannot affect the objective. There is no fallback to MLE or ridge
#' if a required Firth fit fails.
#'
#' Constant columns (RMS centered scale at most \eqn{10^{-8}}) are removed
#' internally and returned with zero coefficients. Each standardized group is
#' reduced by SVD using relative tolerance \eqn{10^{-10}} times its largest
#' dimension. The preprocessing record reports retained ranks, transformations,
#' and excluded columns. Coefficients and predictions are on the original scale.
#' Group labels need not be consecutive or numerically ordered. Missing data
#' and nonnumeric predictors are rejected; no imputation is performed.
#'
#' Automatic finite lambda grids start at the maximum weighted group null-score
#' norm divided by alpha when alpha is positive. For alpha zero, the undivided
#' score is only a reference scale. Neither reference is asserted to produce a
#' null model for a positive target. The default path ends at 0.005 times its
#' start. Specify a broader positive finite lambda grid when appropriate.
#' Zero and infinite lambda are not part of this API. CV reports and warns about
#' endpoint selection; it does not claim optimality outside the requested grid.
#'
#' Binomial fits require \code{standardize=TRUE}, \code{bilevel=FALSE},
#' \code{beta_start=NULL}, and eager coefficient transformation. Gaussian
#' SSR/SSR_fast rules and dfmax/gmax stopping are not transplanted to logistic
#' regression. The default screen setting selects the hybrid active set, with
#' full final KKT checks. \code{SSR_fast} is rejected. No path truncation by
#' training deviance is enabled. All requested points must pass numerical checks.
#'
#' \code{binomial.control} accepts \code{block_max_sweeps} (4000),
#' \code{block_chunk_sweeps} (100), \code{block_stall_window} (200),
#' \code{block_stall_relative_improvement} (0.01),
#' \code{apg_max_iterations} (10000), \code{apg_kkt_check_interval} (5),
#' \code{kkt_tolerance} (\code{eps}), \code{update_tolerance} (1e-10),
#' \code{intercept_tolerance} (1e-12), \code{max_intercept_iterations} (100),
#' \code{irls_max_outer} (50), and \code{irls_max_inner} (2000).
#' When omitted, binomial \code{eps} is 1e-10 (Gaussian default unchanged).
#' \code{max_iter} caps each main iteration budget; it is not a wall-time limit.
#' Training Firth fits use BFGS with at most four 5000-iteration passes and five
#' score-polishing steps, requiring adjusted-score and Fisher-step norms <=1e-6.
#'
#' \code{cv.sglasso} returns fold-local diagnostics and an observation-weighted
#' mean log-loss. Its \code{cvse} describes across-fold variation and is not an
#' independent-test confidence interval. Automatic folds are stratified and use
#' the caller's RNG; set a seed or supply \code{fold} for reproducibility.
#' Supplied lambda values are a fixed absolute grid; automatic CV paths are
#' aligned by relative position and recomputed within each training fold.
#' Selection uses minimum log-loss, breaking exact ties by path order (larger
#' lambda, then smaller d). A final full-training fit supplies predictions.
#' Effective degrees of freedom are not implemented: \code{logLik} returns
#' unpenalized log-likelihood with NA degrees of freedom and \code{select}
#' rejects binomial information-criterion selection. Use log-loss CV instead.
#'
#' @examples
#' set.seed(19)
#' x <- matrix(rnorm(240), 60, 4)
#' y <- rbinom(60, 1, plogis(x[, 1]))
#' fit <- sglasso(x, y, c(1, 1, 2, 2), family = "binomial",
#'                lambda = c(0.3, 0.1), d = c(0, 0.5))
#' predict(fit, x[1:3, , drop = FALSE], type = "response")
#' @seealso \code{\link{sglasso}}, \code{\link{cv.sglasso}}
NULL
