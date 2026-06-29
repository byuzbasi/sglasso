#' Benchmark SGLASSO Screening Rules
#'
#' Runs a focused benchmark comparing \code{"none"}, \code{"SSR"}, and
#' \code{"SSR_fast"} screening rules on block-correlated synthetic designs.
#'
#' @param scenarios A data frame with columns \code{n}, \code{p}, and \code{J}.
#' @param nlambda Number of lambda values.
#' @param active_groups Number of truly active groups.
#' @param rho_within Within-group correlation.
#' @param rho_between Between-group correlation.
#' @param alpha Elastic net mixing parameter passed to \code{\link{sglasso}}.
#' @param d Scale parameter passed to \code{\link{sglasso}}.
#' @param screens Screening rules to compare.
#' @param reps Number of repetitions per scenario.
#' @param seed Random seed.
#' @param standardize Whether \code{\link{sglasso}} should standardize internally.
#' @param eps Convergence tolerance.
#' @param verbose Print progress messages.
#'
#' @return A list containing raw run summaries, per-screen timing summaries,
#' and pairwise accuracy comparisons against \code{screen = "none"}.
#' @export
benchmark_sglasso_screen_rules <- function(
    scenarios = data.frame(n = c(200L, 300L),
                           p = c(20000L, 50000L),
                           J = c(2000L, 5000L)),
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
    eps = 1e-4,
    verbose = TRUE) {
  if (!is.data.frame(scenarios) || !all(c("n", "p", "J") %in% names(scenarios))) {
    stop("scenarios must be a data frame with columns n, p, and J.", call. = FALSE)
  }
  screens <- match.arg(screens, choices = c("none", "SSR", "SSR_fast"),
                       several.ok = TRUE)
  if (rho_between < 0 || rho_within < rho_between || rho_within >= 1) {
    stop("Require 0 <= rho_between <= rho_within < 1.", call. = FALSE)
  }

  results <- list()
  row_id <- 1L
  set.seed(seed)

  for (s in seq_len(nrow(scenarios))) {
    n <- as.integer(scenarios$n[s])
    p <- as.integer(scenarios$p[s])
    J <- as.integer(scenarios$J[s])
    if (p %% J != 0L) {
      stop("Each scenario must satisfy p %% J == 0.", call. = FALSE)
    }

    for (rep_id in seq_len(reps)) {
      dat <- block_screen_benchmark_data(n = n, p = p, J = J,
                                         active_groups = active_groups,
                                         rho_within = rho_within,
                                         rho_between = rho_between)
      fit_list <- list()
      beta_list <- list()

      for (screen in screens) {
        if (isTRUE(verbose)) {
          message(sprintf("screen benchmark: n=%d p=%d J=%d rep=%d screen=%s",
                          n, p, J, rep_id, screen))
        }
        elapsed <- system.time({
          fit <- sglasso(X = dat$X, Y = dat$y, group = dat$group,
                         alpha = alpha, d = d, nlambda = nlambda,
                         screen = screen,
                         standardize = standardize, eps = eps,
                         diagnostics = TRUE)
        })[["elapsed"]]

        beta <- coef(fit, drop = FALSE)
        fit_list[[screen]] <- fit
        beta_list[[screen]] <- beta

        stats <- fit$screen_stats
        results[[row_id]] <- data.frame(
          scenario = s,
          rep = rep_id,
          n = n,
          p = p,
          J = J,
          group_size = p / J,
          screen = screen,
          elapsed_time = elapsed,
          bcd_time = sum(stats$bcd_time),
          screening_time = sum(stats$screening_time),
          kkt_checking_time = sum(stats$kkt_checking_time),
          KKT_violations = sum(stats$n_kkt_violations),
          screened_groups_mean = mean(stats$n_groups_discarded_by_screen),
          kept_after_screening_mean = mean(stats$n_groups_kept_after_screening),
          screening_rate_mean = mean(stats$screening_rate),
          kept_fraction_mean = mean(stats$kept_fraction),
          violation_rate = mean(stats$has_kkt_violation),
          active_groups_mean = mean(stats$n_groups_active),
          working_set_groups_mean = mean(stats$n_groups_working_set),
          refits = sum(stats$n_refits),
          strong_kkt_checked = sum(stats$n_strong_kkt_checked),
          strong_kkt_violations = sum(stats$n_strong_kkt_violations),
          rest_kkt_checked = sum(stats$n_rest_kkt_checked),
          rest_kkt_violations = sum(stats$n_rest_kkt_violations),
          bcd_iterations = sum(stats$n_bcd_iterations),
          final_inactive_kkt_residual = max(stats$final_inactive_kkt_residual),
          max_abs_diff_vs_none = NA_real_,
          L2_diff_vs_none = NA_real_,
          selected_groups_diff_vs_none = NA_integer_,
          first_lambda_diff_vs_none = NA_real_,
          stringsAsFactors = FALSE
        )
        row_id <- row_id + 1L
      }

      if ("none" %in% names(beta_list)) {
        b_none <- beta_list[["none"]]
        none_groups <- selected_groups_from_beta(b_none, dat$group)
        for (screen in setdiff(names(beta_list), "none")) {
          b_screen <- beta_list[[screen]]
          screen_groups <- selected_groups_from_beta(b_screen, dat$group)
          beta_diff <- drop(b_screen - b_none)
          stats <- fit_list[[screen]]$screen_stats
          results[[row_id]] <- data.frame(
            scenario = s,
            rep = rep_id,
            n = n,
            p = p,
            J = J,
            group_size = p / J,
            screen = paste0(screen, "_vs_none"),
            elapsed_time = NA_real_,
            bcd_time = NA_real_,
            screening_time = NA_real_,
            kkt_checking_time = NA_real_,
            KKT_violations = sum(stats$n_kkt_violations),
            screened_groups_mean = mean(stats$n_groups_discarded_by_screen),
            kept_after_screening_mean = mean(stats$n_groups_kept_after_screening),
            screening_rate_mean = mean(stats$screening_rate),
            kept_fraction_mean = mean(stats$kept_fraction),
            violation_rate = mean(stats$has_kkt_violation),
            active_groups_mean = mean(stats$n_groups_active),
            working_set_groups_mean = mean(stats$n_groups_working_set),
            refits = sum(stats$n_refits),
            strong_kkt_checked = sum(stats$n_strong_kkt_checked),
            strong_kkt_violations = sum(stats$n_strong_kkt_violations),
            rest_kkt_checked = sum(stats$n_rest_kkt_checked),
            rest_kkt_violations = sum(stats$n_rest_kkt_violations),
            bcd_iterations = sum(stats$n_bcd_iterations),
            final_inactive_kkt_residual = max(stats$final_inactive_kkt_residual),
            max_abs_diff_vs_none = max(abs(beta_diff), na.rm = TRUE),
            L2_diff_vs_none = sqrt(sum(beta_diff^2, na.rm = TRUE)),
            selected_groups_diff_vs_none = length(setdiff(screen_groups, none_groups)) +
              length(setdiff(none_groups, screen_groups)),
            first_lambda_diff_vs_none = first_lambda_difference(
              b_screen, b_none, fit_list[[screen]]$lambda),
            stringsAsFactors = FALSE
          )
          row_id <- row_id + 1L
        }
      }
    }
  }

  raw <- do.call(rbind, results)
  timing_rows <- raw[raw$screen %in% screens, , drop = FALSE]
  accuracy_rows <- raw[!raw$screen %in% screens, , drop = FALSE]
  summary <- summarize_screen_benchmark(timing_rows, accuracy_rows)

  list(raw = raw, timing = timing_rows, accuracy = accuracy_rows, summary = summary)
}

block_screen_benchmark_data <- function(n, p, J, active_groups,
                                        rho_within, rho_between) {
  group_size <- p / J
  group <- rep(seq_len(J), each = group_size)
  common_factor <- matrix(rnorm(n), n, 1)
  group_factor <- matrix(rnorm(n * J), n, J)
  noise <- matrix(rnorm(n * p), n, p)
  X <- sqrt(rho_between) * common_factor[, 1] +
    sqrt(rho_within - rho_between) * group_factor[, group] +
    sqrt(1 - rho_within) * noise
  X <- scale(X)

  active <- seq_len(min(active_groups, J))
  beta <- rep(0, p)
  for (g in active) {
    ind <- which(group == g)
    beta[ind] <- 1 / sqrt(length(ind))
  }
  y <- drop(X %*% beta + rnorm(n))
  list(X = X, y = y, group = group, beta = beta, active_groups = active)
}

selected_groups_from_beta <- function(beta_array, group, tol = 1e-8) {
  beta_mat <- drop(beta_array[-1, , , drop = FALSE])
  if (is.null(dim(beta_mat))) beta_mat <- matrix(beta_mat, ncol = 1)
  selected_features <- which(rowSums(abs(beta_mat) > tol) > 0)
  sort(unique(group[selected_features]))
}

first_lambda_difference <- function(beta_a, beta_b, lambda, tol = 1e-8) {
  diff_arr <- abs(beta_a - beta_b)
  dims <- dim(diff_arr)
  if (length(dims) < 2L) return(NA_real_)
  lambda_diff <- apply(diff_arr, 2L, max, na.rm = TRUE)
  hit <- which(lambda_diff > tol)
  if (length(hit) == 0L) NA_real_ else lambda[hit[1L]]
}

summarize_screen_benchmark <- function(timing_rows, accuracy_rows) {
  timing <- aggregate(
    cbind(elapsed_time, bcd_time, screening_time, kkt_checking_time,
          KKT_violations, screened_groups_mean, kept_after_screening_mean,
          screening_rate_mean, kept_fraction_mean, violation_rate,
          working_set_groups_mean, refits, strong_kkt_checked,
          strong_kkt_violations, rest_kkt_checked,
          rest_kkt_violations, bcd_iterations) ~ screen,
    data = timing_rows,
    FUN = median,
    na.rm = TRUE
  )
  none_time <- timing$elapsed_time[timing$screen == "none"]
  timing$speedup_vs_none <- if (length(none_time)) none_time / timing$elapsed_time else NA_real_
  none_kkt_time <- timing$kkt_checking_time[timing$screen == "none"]
  timing$kkt_time_saved_by_ssr <- if (length(none_kkt_time)) none_kkt_time - timing$kkt_checking_time else NA_real_
  timing$screening_overhead <- timing$screening_time
  timing$net_ssr_gain <- timing$kkt_time_saved_by_ssr - timing$screening_overhead
  none_bcd_time <- timing$bcd_time[timing$screen == "none"]
  timing$bcd_time_difference <- if (length(none_bcd_time)) none_bcd_time - timing$bcd_time else NA_real_
  none_bcd_iter <- timing$bcd_iterations[timing$screen == "none"]
  timing$bcd_iteration_difference <- if (length(none_bcd_iter)) none_bcd_iter - timing$bcd_iterations else NA_real_
  total_time <- timing$elapsed_time
  timing$time_share_screening <- timing$screening_time / total_time
  timing$time_share_bcd <- timing$bcd_time / total_time
  timing$time_share_kkt <- timing$kkt_checking_time / total_time
  names(timing)[names(timing) == "elapsed_time"] <- "median_time"

  acc <- accuracy_rows[, c("screen", "max_abs_diff_vs_none", "L2_diff_vs_none",
                           "selected_groups_diff_vs_none", "first_lambda_diff_vs_none",
                           "final_inactive_kkt_residual"), drop = FALSE]
  names(acc)[1L] <- "comparison"
  max_beta_diff <- setNames(rep(0, nrow(timing)), timing$screen)
  for (i in seq_len(nrow(acc))) {
    screen <- sub("_vs_none$", "", acc$comparison[i])
    max_beta_diff[screen] <- acc$max_abs_diff_vs_none[i]
  }
  table <- data.frame(
    screen = timing$screen,
    median_time = timing$median_time,
    speedup_vs_none = timing$speedup_vs_none,
    max_beta_diff_vs_none = unname(max_beta_diff[timing$screen]),
    KKT_violations = timing$KKT_violations,
    screened_groups_mean = timing$screened_groups_mean,
    kept_after_screening_mean = timing$kept_after_screening_mean,
    screening_rate_mean = timing$screening_rate_mean,
    kept_fraction_mean = timing$kept_fraction_mean,
    violation_rate = timing$violation_rate,
    working_set_groups_mean = timing$working_set_groups_mean,
    refits = timing$refits,
    kkt_time_saved_by_ssr = timing$kkt_time_saved_by_ssr,
    screening_overhead = timing$screening_overhead,
    net_ssr_gain = timing$net_ssr_gain,
    bcd_time_difference = timing$bcd_time_difference,
    bcd_iteration_difference = timing$bcd_iteration_difference,
    stringsAsFactors = FALSE
  )
  list(timing = timing, accuracy = acc, table = table)
}
