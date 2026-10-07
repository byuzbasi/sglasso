binary_fixture <- function() {
  set.seed(107)
  x <- matrix(rnorm(480), 80, 6)
  x[, 2] <- 0.8 * x[, 1] + 0.2 * x[, 2]
  y <- rbinom(80, 1, plogis(x[, 1] - 0.5 * x[, 3]))
  list(x = x, y = y, g = rep(1:3, each = 2))
}

test_that("binary path reconstructs objective and probabilities", {
  z <- binary_fixture()
  f <- sglasso(z$x, z$y, z$g, family = "binomial", lambda = c(0.3, 0.1), d = c(0, 0.5))
  expect_identical(f$solver, "local_irls_proximal_newton")
  expect_true(all(f$converged))
  expect_lte(max(f$kkt), 1e-10)
  b <- coef(f, drop = FALSE)
  eta <- predict(f, z$x, type = "link", drop = FALSE)
  expect_equal(predict(f, z$x, drop = FALSE), plogis(eta))
  expect_equal(dim(eta), c(80L, 2L, 2L))
  expect_equal(eta[, 1, 1], drop(cbind(1, z$x) %*% b[, 1, 1]))
  for (j in 1:2) for (i in 1:2) {
    v <- eta[, i, j]
    loss <- mean(pmax(v, 0) - z$y * v + log1p(exp(-abs(v))))
    beta <- f$solver_diagnostics$beta[, i, j]
    penalty <- sum(vapply(seq_along(f$preprocess$blocks), function(g) {
      ix <- f$preprocess$blocks[[g]]$solver_index
      f$lambda[i] * f$preprocess$group_weight[g] * (
        f$alpha * sqrt(sum(beta[ix]^2)) + (1 - f$alpha) / 2 *
          sum((beta[ix] - f$d[j] * f$target[ix])^2))
    }, numeric(1)))
    expect_equal(f$objective[i, j], loss + penalty, tolerance = 1e-9)
    expect_equal(as.numeric(logLik(f))[i + (j - 1) * 2], -80 * loss)
  }
  expect_error(select(f), "degrees of freedom")
})

test_that("boundary alpha/d, separation and rank-deficient groups are handled", {
  z <- binary_fixture()
  for (a in c(0, 1)) {
    f <- sglasso(z$x, z$y, z$g, family = "binomial", alpha = a,
      lambda = 0.2, d = c(0, 0.5))
    expect_true(all(f$converged))
    if (a == 1) expect_equal(f$betas[, , 1], f$betas[, , 2], tolerance = 1e-8)
  }
  x <- cbind(z$x[, 1:2], z$x[, 1], constant = 1)
  f <- sglasso(x, as.integer(x[, 1] > 0), c(1, 1, 1, 2),
    family = "binomial", lambda = 0.2, d = 0.5)
  expect_true(all(f$converged))
  expect_equal(f$dropped_constant_columns, 4L)
  expect_true(all(vapply(f$target_diagnostics, function(v) v$success, logical(1))))
  expect_true(all(is.finite(predict(f, x))))
})

test_that("CV recomputes transformations and target within training folds", {
  z <- binary_fixture()
  fold <- rep(1:2, each = 40)
  cv <- suppressWarnings(cv.sglasso(z$x, z$y, z$g, family = "binomial",
    lambda = c(0.3, 0.1), d = c(0, 0.5), fold = fold))
  losses <- array(NA_real_, c(80, 2, 2))
  for (k in 1:2) {
    f <- sglasso(z$x[fold != k, ], z$y[fold != k], z$g,
      family = "binomial", lambda = c(0.3, 0.1), d = c(0, 0.5))
    expect_equal(cv$fold_fits[[k]]$center, colMeans(z$x[fold != k, ]))
    expect_equal(cv$fold_fits[[k]]$target_original, f$target_original)
    eta <- predict(f, z$x[fold == k, ], type = "link", drop = FALSE)
    losses[fold == k, , ] <- pmax(eta, 0) - z$y[fold == k] * eta + log1p(exp(-abs(eta)))
  }
  expect_equal(unname(cv$cve), apply(losses, c(2, 3), mean))
  expect_equal(predict(cv, z$x, s = "opt"),
    plogis(drop(cbind(1, z$x) %*% cv$beta_opt)))
  expect_equal(as.numeric(coef(cv)), as.numeric(cv$beta_opt))
  auto <- suppressWarnings(cv.sglasso(z$x, z$y, z$g, family = "binomial",
    nlambda = 2, d = 0, fold = fold))
  expect_identical(auto$lambda_alignment, "training_fold_relative_scale")
  expect_equal(auto$fold_fits[[1]]$lambda / auto$fold_fits[[1]]$lambda[1], c(1, 0.005))
})

test_that("unsupported inputs and nonconvergence fail explicitly", {
  z <- binary_fixture()
  expect_error(sglasso(z$x, factor(z$y), z$g, family = "binomial"), "0/1")
  expect_error(sglasso(z$x, z$y, z$g, family = "binomial", lambda = 0), "positive finite")
  expect_error(sglasso(z$x, z$y, z$g, family = "binomial", standardize = FALSE), "standardize")
  expect_error(sglasso(z$x, z$y, z$g, family = "binomial", screen = "SSR_fast"), "Gaussian-only")
  expect_error(sglasso(z$x, z$y, z$g, family = "binomial", lambda = 0.001,
    d = 0.5, max_iter = 1, eps = 1e-14), "convergence/KKT")
})

test_that("HD fitting and disabled active sets preserve predictions", {
  set.seed(23)
  x <- matrix(rnorm(30 * 45), 30, 45)
  y <- rep(0:1, 15)
  group <- rep(paste0("group", 1:15), each = 3)
  a <- sglasso(x, y, group, family = "binomial", lambda = 0.3, d = 0.5)
  b <- sglasso(x, y, group, family = "binomial", lambda = 0.3, d = 0.5,
    screen = "none")
  expect_true(all(a$converged))
  expect_equal(predict(a, x), predict(b, x), tolerance = 1e-8)
  expect_false(b$solver_diagnostics$use_active_set)
  order <- rev(seq_len(ncol(x)))
  c <- sglasso(x[, order], y, group[order], family = "binomial",
    lambda = 0.3, d = 0.5)
  expect_equal(predict(a, x), predict(c, x[, order]), tolerance = 1e-8)
})

test_that("unequal folds use observation-weighted log-loss", {
  z <- binary_fixture()
  cv <- suppressWarnings(cv.sglasso(z$x, z$y, z$g, family = "binomial",
    lambda = c(0.3, 0.1), d = 0, fold = c(rep("a", 25), rep("b", 55))))
  expect_equal(as.numeric(cv$cve), as.numeric((cv$fold_errors[, , 1] * 25 +
    cv$fold_errors[, , 2] * 55) / 80))
})
