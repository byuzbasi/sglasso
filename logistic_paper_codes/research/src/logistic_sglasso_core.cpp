#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

constexpr double kMajorizationCurvature = 0.25;

inline double stable_log1pexp(const double x) {
  if (x > 0.0) return x + std::log1p(std::exp(-x));
  return std::log1p(std::exp(x));
}

arma::vec stable_expit(const arma::vec& eta) {
  arma::vec out(eta.n_elem);
  for (arma::uword i = 0; i < eta.n_elem; ++i) {
    const double x = eta[i];
    if (x >= 0.0) {
      const double z = std::exp(-x);
      out[i] = 1.0 / (1.0 + z);
    } else {
      const double z = std::exp(x);
      out[i] = z / (1.0 + z);
    }
  }
  return out;
}

double logistic_loss_impl(const arma::vec& eta, const arma::vec& y) {
  double out = 0.0;
  for (arma::uword i = 0; i < eta.n_elem; ++i) {
    out += stable_log1pexp(eta[i]) - y[i] * eta[i];
  }
  return out / static_cast<double>(eta.n_elem);
}

double penalty_impl(const arma::vec& beta,
                    const arma::uvec& group_start,
                    const arma::uvec& group_end,
                    const arma::vec& group_weight,
                    const arma::vec& target,
                    const double lambda,
                    const double alpha,
                    const double d) {
  double out = 0.0;
  for (arma::uword g = 0; g < group_start.n_elem; ++g) {
    const arma::uword first = group_start[g];
    const arma::uword last = group_end[g];
    const arma::vec bg = beta.subvec(first, last);
    const arma::vec shifted = bg - d * target.subvec(first, last);
    const double lambda1 = lambda * alpha * group_weight[g];
    const double lambda2 = lambda * (1.0 - alpha) * group_weight[g];
    out += lambda1 * arma::norm(bg, 2);
    out += 0.5 * lambda2 * arma::dot(shifted, shifted);
  }
  return out;
}

double objective_from_eta_impl(const arma::vec& eta,
                               const arma::vec& y,
                               const arma::vec& beta,
                               const arma::uvec& group_start,
                               const arma::uvec& group_end,
                               const arma::vec& group_weight,
                               const arma::vec& target,
                               const double lambda,
                               const double alpha,
                               const double d) {
  return logistic_loss_impl(eta, y) + penalty_impl(
    beta, group_start, group_end, group_weight, target, lambda, alpha, d
  );
}

double objective_impl(const arma::mat& X,
                      const arma::vec& y,
                      const arma::vec& beta,
                      const double intercept,
                      const arma::uvec& group_start,
                      const arma::uvec& group_end,
                      const arma::vec& group_weight,
                      const arma::vec& target,
                      const double lambda,
                      const double alpha,
                      const double d) {
  return objective_from_eta_impl(
    intercept + X * beta, y, beta, group_start, group_end, group_weight,
    target, lambda, alpha, d
  );
}

double kkt_residual_impl(const arma::mat& X,
                         const arma::vec& y,
                         const arma::vec& beta,
                         const double intercept,
                         const arma::uvec& group_start,
                         const arma::uvec& group_end,
                         const arma::vec& group_weight,
                         const arma::vec& target,
                         const double lambda,
                         const double alpha,
                         const double d,
                         arma::vec* group_residuals = nullptr,
                         double* intercept_residual = nullptr) {
  const double n = static_cast<double>(X.n_rows);
  const arma::vec probability = stable_expit(intercept + X * beta);
  const arma::vec gradient = X.t() * (probability - y) / n;
  const double intercept_kkt = std::abs(arma::mean(probability - y));
  double max_residual = intercept_kkt;

  if (group_residuals != nullptr) {
    group_residuals->set_size(group_start.n_elem);
    group_residuals->zeros();
  }

  for (arma::uword g = 0; g < group_start.n_elem; ++g) {
    const arma::uword first = group_start[g];
    const arma::uword last = group_end[g];
    const arma::vec bg = beta.subvec(first, last);
    const arma::vec tg = target.subvec(first, last);
    const arma::vec gg = gradient.subvec(first, last);
    const double beta_norm = arma::norm(bg, 2);
    const double lambda1 = lambda * alpha * group_weight[g];
    const double lambda2 = lambda * (1.0 - alpha) * group_weight[g];
    double residual = 0.0;

    if (beta_norm > 1e-10) {
      const arma::vec stationarity =
        gg + lambda2 * (bg - d * tg) + lambda1 * bg / beta_norm;
      residual = arma::norm(stationarity, 2);
    } else {
      const arma::vec smooth_at_zero = gg - lambda2 * d * tg;
      residual = std::max(0.0, arma::norm(smooth_at_zero, 2) - lambda1);
    }

    if (group_residuals != nullptr) (*group_residuals)[g] = residual;
    if (residual > max_residual) max_residual = residual;
  }

  if (intercept_residual != nullptr) *intercept_residual = intercept_kkt;
  return max_residual;
}

struct OneFit {
  arma::vec beta;
  double intercept;
  double objective;
  double kkt;
  double intercept_kkt;
  int outer_iterations;
  int inner_iterations;
  bool converged;
  double max_raw_objective_increase;
  int active_set_scans;
  double group_updates;
  std::vector<double> objective_trace;
};

OneFit fit_one(const arma::mat& X,
               const arma::vec& y,
               const arma::uvec& group_start,
               const arma::uvec& group_end,
               const arma::vec& group_weight,
               const arma::vec& target,
               const double lambda,
               const double alpha,
               const double d,
               arma::vec beta,
               double intercept,
               const int max_outer,
               const int max_inner,
               const double tolerance,
               const double inner_tolerance,
               const bool use_active_set,
               const bool keep_trace) {
  const double n = static_cast<double>(X.n_rows);
  double objective = objective_impl(
    X, y, beta, intercept, group_start, group_end, group_weight,
    target, lambda, alpha, d
  );

  OneFit result;
  result.outer_iterations = 0;
  result.inner_iterations = 0;
  result.converged = false;
  result.max_raw_objective_increase = 0.0;
  result.active_set_scans = 0;
  result.group_updates = 0.0;
  if (keep_trace) result.objective_trace.push_back(objective);

  for (int outer = 0; outer < max_outer; ++outer) {
    const arma::vec beta_old = beta;
    const double intercept_old = intercept;
    const double objective_old = objective;

    const arma::vec eta = intercept + X * beta;
    const arma::vec probability = stable_expit(eta);
    const arma::vec working_response =
      eta + (y - probability) / kMajorizationCurvature;
    arma::vec residual = working_response - intercept - X * beta;

    const arma::uword groups = group_start.n_elem;
    arma::uvec active(groups, arma::fill::ones);
    if (use_active_set) {
      active.zeros();
      for (arma::uword g = 0; g < groups; ++g) {
        if (arma::norm(beta.subvec(group_start[g], group_end[g]), 2) > 1e-12) {
          active[g] = 1u;
        }
      }
    }

    int inner_used = 0;
    bool surrogate_finished = false;
    while (inner_used < max_inner && !surrogate_finished) {
      bool active_converged = false;
      while (inner_used < max_inner) {
        double max_delta = 0.0;

        const double intercept_delta = arma::mean(residual);
        if (intercept_delta != 0.0) {
          intercept += intercept_delta;
          residual -= intercept_delta;
          max_delta = std::abs(intercept_delta);
        }

        for (arma::uword g = 0; g < groups; ++g) {
          if (active[g] == 0u) continue;
          const arma::uword first = group_start[g];
          const arma::uword last = group_end[g];
          const auto Xg = X.cols(first, last);
          const arma::vec beta_current = beta.subvec(first, last);
          const arma::vec partial_score =
            beta_current + Xg.t() * residual / n;
          const double lambda1 = lambda * alpha * group_weight[g];
          const double lambda2 = lambda * (1.0 - alpha) * group_weight[g];
          const arma::vec shifted_score =
            kMajorizationCurvature * partial_score +
            lambda2 * d * target.subvec(first, last);
          const double shifted_norm = arma::norm(shifted_score, 2);

          arma::vec beta_new(beta_current.n_elem, arma::fill::zeros);
          if (shifted_norm > lambda1 && shifted_norm > 0.0) {
            const double multiplier =
              (1.0 - lambda1 / shifted_norm) /
              (kMajorizationCurvature + lambda2);
            beta_new = multiplier * shifted_score;
          }

          const arma::vec delta = beta_new - beta_current;
          const double group_delta = arma::abs(delta).max();
          if (group_delta > 0.0) {
            beta.subvec(first, last) = beta_new;
            residual -= Xg * delta;
            if (group_delta > max_delta) max_delta = group_delta;
          }
          result.group_updates += 1.0;
        }

        ++inner_used;
        if (max_delta <= inner_tolerance *
            (1.0 + std::max(std::abs(intercept), arma::abs(beta).max()))) {
          active_converged = true;
          break;
        }
      }

      if (!use_active_set || !active_converged) {
        surrogate_finished = true;
        continue;
      }

      bool violation = false;
      ++result.active_set_scans;
      for (arma::uword g = 0; g < groups; ++g) {
        if (active[g] != 0u) continue;
        const arma::uword first = group_start[g];
        const arma::uword last = group_end[g];
        const auto Xg = X.cols(first, last);
        const arma::vec partial_score = Xg.t() * residual / n;
        const double lambda1 = lambda * alpha * group_weight[g];
        const double lambda2 = lambda * (1.0 - alpha) * group_weight[g];
        const arma::vec shifted_score =
          kMajorizationCurvature * partial_score +
          lambda2 * d * target.subvec(first, last);
        const double scan_tolerance =
          1e-12 * (1.0 + lambda1 + arma::norm(shifted_score, 2));
        if (arma::norm(shifted_score, 2) > lambda1 + scan_tolerance) {
          active[g] = 1u;
          violation = true;
        }
      }
      surrogate_finished = !violation;
    }
    result.inner_iterations += inner_used;

    double candidate_objective = objective_impl(
      X, y, beta, intercept, group_start, group_end, group_weight,
      target, lambda, alpha, d
    );
    const double raw_increase = candidate_objective - objective_old;
    if (raw_increase > result.max_raw_objective_increase) {
      result.max_raw_objective_increase = raw_increase;
    }

    // The fixed 1/4 bound and a converged inner solve should already provide
    // descent. This safeguard only handles finite precision or an early inner
    // termination without changing the target objective.
    if (!std::isfinite(candidate_objective) ||
        candidate_objective > objective_old + 1e-12) {
      double step = 0.5;
      bool accepted = false;
      for (int line = 0; line < 40; ++line) {
        const arma::vec beta_candidate = beta_old + step * (beta - beta_old);
        const double intercept_candidate =
          intercept_old + step * (intercept - intercept_old);
        const double line_objective = objective_impl(
          X, y, beta_candidate, intercept_candidate, group_start, group_end,
          group_weight, target, lambda, alpha, d
        );
        if (std::isfinite(line_objective) &&
            line_objective <= objective_old + 1e-12) {
          beta = beta_candidate;
          intercept = intercept_candidate;
          candidate_objective = line_objective;
          accepted = true;
          break;
        }
        step *= 0.5;
      }
      if (!accepted) {
        beta = beta_old;
        intercept = intercept_old;
        candidate_objective = objective_old;
      }
    }

    objective = candidate_objective;
    result.outer_iterations = outer + 1;
    if (keep_trace) result.objective_trace.push_back(objective);

    double intercept_kkt = 0.0;
    const double kkt = kkt_residual_impl(
      X, y, beta, intercept, group_start, group_end, group_weight,
      target, lambda, alpha, d, nullptr, &intercept_kkt
    );
    const double coefficient_change = std::max(
      std::abs(intercept - intercept_old), arma::abs(beta - beta_old).max()
    );
    const double relative_objective_change =
      std::abs(objective_old - objective) / (1.0 + std::abs(objective_old));

    if (kkt <= std::max(10.0 * tolerance, 1e-7) &&
        (coefficient_change <= std::sqrt(tolerance) *
           (1.0 + std::max(std::abs(intercept), arma::abs(beta).max())) ||
         relative_objective_change <= tolerance)) {
      result.converged = true;
      break;
    }
  }

  result.beta = beta;
  result.intercept = intercept;
  result.objective = objective;
  result.kkt = kkt_residual_impl(
    X, y, beta, intercept, group_start, group_end, group_weight,
    target, lambda, alpha, d, nullptr, &result.intercept_kkt
  );
  if (result.kkt <= std::max(25.0 * tolerance, 2.5e-7)) {
    result.converged = true;
  }
  return result;
}

struct ABGDFit {
  arma::vec beta;
  double intercept;
  double objective;
  double kkt;
  double intercept_kkt;
  int passes;
  bool converged;
  double max_raw_objective_increase;
  int kkt_scans;
  double group_updates;
  arma::uvec active;
  std::vector<double> objective_trace;
};

arma::vec group_curvature_impl(const arma::mat& X,
                               const arma::uvec& group_start,
                               const arma::uvec& group_end) {
  const double n = static_cast<double>(X.n_rows);
  arma::vec curvature(group_start.n_elem, arma::fill::zeros);
  for (arma::uword g = 0; g < group_start.n_elem; ++g) {
    const arma::mat Xg = X.cols(group_start[g], group_end[g]);
    const arma::mat gram = Xg.t() * Xg / n;
    arma::vec eigenvalues;
    const bool eigensolved = arma::eig_sym(eigenvalues, gram);
    double maximum_eigenvalue = 0.0;
    if (eigensolved && eigenvalues.n_elem > 0u && eigenvalues.is_finite()) {
      maximum_eigenvalue = eigenvalues.max();
    } else {
      maximum_eigenvalue = arma::norm(gram, 2);
    }
    if (!std::isfinite(maximum_eigenvalue) || maximum_eigenvalue < 0.0) {
      Rcpp::stop("Unable to compute a finite group curvature bound.");
    }
    curvature[g] = std::max(
      kMajorizationCurvature * maximum_eigenvalue,
      10.0 * std::numeric_limits<double>::epsilon()
    );
  }
  return curvature;
}

ABGDFit fit_one_abgd(const arma::mat& X,
                     const arma::vec& y,
                     const arma::uvec& group_start,
                     const arma::uvec& group_end,
                     const arma::vec& group_weight,
                     const arma::vec& target,
                     const arma::vec& group_curvature,
                     const double lambda,
                     const double alpha,
                     const double d,
                     arma::vec beta,
                     double intercept,
                     arma::uvec active,
                     const int max_passes,
                     const double tolerance,
                     const double update_tolerance,
                     const bool use_active_set,
                     const bool keep_trace) {
  const double n = static_cast<double>(X.n_rows);
  const arma::uword groups = group_start.n_elem;
  const double kkt_tolerance = std::max(2.0 * tolerance, 2e-8);
  const int scan_frequency = 10;

  if (active.n_elem != groups) active.zeros(groups);
  if (!use_active_set) active.ones(groups);
  for (arma::uword g = 0; g < groups; ++g) {
    if (arma::norm(beta.subvec(group_start[g], group_end[g]), 2) > 1e-12) {
      active[g] = 1u;
    }
  }

  arma::vec eta = intercept + X * beta;
  arma::vec probability = stable_expit(eta);
  arma::vec residual = y - probability;
  double objective = objective_from_eta_impl(
    eta, y, beta, group_start, group_end, group_weight,
    target, lambda, alpha, d
  );

  ABGDFit result;
  result.passes = 0;
  result.converged = false;
  result.max_raw_objective_increase = 0.0;
  result.kkt_scans = 0;
  result.group_updates = 0.0;
  if (keep_trace) result.objective_trace.push_back(objective);

  arma::vec group_kkt;
  double intercept_kkt = 0.0;
  double kkt = kkt_residual_impl(
    X, y, beta, intercept, group_start, group_end, group_weight,
    target, lambda, alpha, d, &group_kkt, &intercept_kkt
  );
  ++result.kkt_scans;
  if (use_active_set) {
    for (arma::uword g = 0; g < groups; ++g) {
      if (group_kkt[g] > kkt_tolerance) active[g] = 1u;
    }
  }
  if (kkt <= kkt_tolerance) result.converged = true;

  for (int pass = 0; pass < max_passes && !result.converged; ++pass) {
    if (pass % 10 == 0) Rcpp::checkUserInterrupt();
    const arma::vec beta_old = beta;
    const double intercept_old = intercept;
    const arma::vec eta_old = eta;
    const double objective_old = objective;

    const double intercept_delta = arma::mean(residual) /
      kMajorizationCurvature;
    if (intercept_delta != 0.0) {
      intercept += intercept_delta;
      eta += intercept_delta;
      probability = stable_expit(eta);
      residual = y - probability;
    }

    for (arma::uword g = 0; g < groups; ++g) {
      if (active[g] == 0u) continue;
      const arma::uword first = group_start[g];
      const arma::uword last = group_end[g];
      const auto Xg = X.cols(first, last);
      const arma::vec beta_current = beta.subvec(first, last);
      const double lambda1 = lambda * alpha * group_weight[g];
      const double lambda2 = lambda * (1.0 - alpha) * group_weight[g];
      const arma::vec score = Xg.t() * residual / n;
      const arma::vec shifted_score =
        group_curvature[g] * beta_current + score +
        lambda2 * d * target.subvec(first, last);
      const double shifted_norm = arma::norm(shifted_score, 2);

      arma::vec beta_new(beta_current.n_elem, arma::fill::zeros);
      if (shifted_norm > lambda1 && shifted_norm > 0.0) {
        beta_new = ((1.0 - lambda1 / shifted_norm) /
          (group_curvature[g] + lambda2)) * shifted_score;
      }

      const arma::vec delta = beta_new - beta_current;
      if (arma::any(delta != 0.0)) {
        beta.subvec(first, last) = beta_new;
        eta += Xg * delta;
        probability = stable_expit(eta);
        residual = y - probability;
      }
      result.group_updates += 1.0;
    }

    double candidate_objective = objective_from_eta_impl(
      eta, y, beta, group_start, group_end, group_weight,
      target, lambda, alpha, d
    );
    const double raw_increase = candidate_objective - objective_old;
    if (raw_increase > result.max_raw_objective_increase) {
      result.max_raw_objective_increase = raw_increase;
    }

    const double objective_slack =
      1e-12 * (1.0 + std::abs(objective_old));
    if (!std::isfinite(candidate_objective) ||
        candidate_objective > objective_old + objective_slack) {
      const arma::vec beta_direction = beta - beta_old;
      const double intercept_direction = intercept - intercept_old;
      const arma::vec eta_direction = eta - eta_old;
      double step = 0.5;
      bool accepted = false;
      for (int line = 0; line < 40; ++line) {
        const arma::vec beta_candidate = beta_old + step * beta_direction;
        const arma::vec eta_candidate = eta_old + step * eta_direction;
        const double line_objective = objective_from_eta_impl(
          eta_candidate, y, beta_candidate, group_start, group_end,
          group_weight, target, lambda, alpha, d
        );
        if (std::isfinite(line_objective) &&
            line_objective <= objective_old + objective_slack) {
          beta = beta_candidate;
          intercept = intercept_old + step * intercept_direction;
          eta = eta_candidate;
          candidate_objective = line_objective;
          accepted = true;
          break;
        }
        step *= 0.5;
      }
      if (!accepted) {
        beta = beta_old;
        intercept = intercept_old;
        eta = eta_old;
        candidate_objective = objective_old;
      }
      probability = stable_expit(eta);
      residual = y - probability;
    }

    objective = candidate_objective;
    result.passes = pass + 1;
    if (keep_trace) result.objective_trace.push_back(objective);

    const double coefficient_change = std::max(
      std::abs(intercept - intercept_old), arma::abs(beta - beta_old).max()
    );
    const double coefficient_scale = 1.0 + std::max(
      std::abs(intercept), arma::abs(beta).max()
    );
    const double relative_objective_change =
      std::abs(objective_old - objective) / (1.0 + std::abs(objective_old));
    const bool small_update =
      coefficient_change <= update_tolerance * coefficient_scale;
    const bool scan_now = small_update ||
      relative_objective_change <= tolerance ||
      ((pass + 1) % scan_frequency == 0) ||
      (pass + 1 == max_passes);

    if (scan_now) {
      kkt = kkt_residual_impl(
        X, y, beta, intercept, group_start, group_end, group_weight,
        target, lambda, alpha, d, &group_kkt, &intercept_kkt
      );
      ++result.kkt_scans;
      bool activated = false;
      if (use_active_set) {
        for (arma::uword g = 0; g < groups; ++g) {
          if (active[g] == 0u && group_kkt[g] > kkt_tolerance) {
            active[g] = 1u;
            activated = true;
          }
        }
      }
      if (!activated && kkt <= kkt_tolerance) {
        result.converged = true;
      }
    }
  }

  kkt = kkt_residual_impl(
    X, y, beta, intercept, group_start, group_end, group_weight,
    target, lambda, alpha, d, nullptr, &intercept_kkt
  );
  ++result.kkt_scans;
  if (kkt <= kkt_tolerance) result.converged = true;

  result.beta = beta;
  result.intercept = intercept;
  result.objective = objective;
  result.kkt = kkt;
  result.intercept_kkt = intercept_kkt;
  result.active = active;
  return result;
}

struct FeasibleInterval {
  bool feasible;
  double lower;
  double upper;
};

FeasibleInterval nonnegative_quadratic_interval(const double A,
                                                const double B,
                                                const double C) {
  const double scale = 1.0 + std::abs(A) + std::abs(B) + std::abs(C);
  const double tol = 1e-12 * scale;
  FeasibleInterval out{false, NA_REAL, NA_REAL};

  if (std::abs(A) <= tol) {
    if (std::abs(B) <= tol) {
      if (C <= tol) return FeasibleInterval{true, 0.0, R_PosInf};
      return out;
    }
    const double root = -C / B;
    if (B < 0.0) {
      return FeasibleInterval{true, std::max(0.0, root), R_PosInf};
    }
    if (root >= 0.0 && C <= tol) {
      return FeasibleInterval{true, 0.0, root};
    }
    return out;
  }

  const double discriminant = B * B - 4.0 * A * C;
  if (discriminant < -tol) {
    if (A < 0.0) return FeasibleInterval{true, 0.0, R_PosInf};
    return out;
  }

  const double sqrt_disc = std::sqrt(std::max(0.0, discriminant));
  double root1 = (-B - sqrt_disc) / (2.0 * A);
  double root2 = (-B + sqrt_disc) / (2.0 * A);
  if (root1 > root2) std::swap(root1, root2);

  if (A > 0.0) {
    const double lower = std::max(0.0, root1);
    const double upper = root2;
    if (upper + tol >= lower && upper >= 0.0) {
      return FeasibleInterval{true, lower, std::max(lower, upper)};
    }
    return out;
  }

  // For A < 0, f(lambda) <= 0 outside [root1, root2]. In the intended
  // SGLASSO setting C = ||score||^2 >= 0, so the relevant nonnegative branch
  // begins at the positive root.
  if (C <= tol && root1 > 0.0) {
    return FeasibleInterval{true, 0.0, root1};
  }
  return FeasibleInterval{true, std::max(0.0, root2), R_PosInf};
}

void validate_groups(const arma::mat& X,
                     const arma::uvec& group_start,
                     const arma::uvec& group_end,
                     const arma::vec& group_weight,
                     const arma::vec& target) {
  if (group_start.n_elem == 0 || group_start.n_elem != group_end.n_elem ||
      group_start.n_elem != group_weight.n_elem) {
    Rcpp::stop("Group boundaries and weights are incompatible.");
  }
  if (target.n_elem != X.n_cols) {
    Rcpp::stop("target must have one value per transformed predictor.");
  }
  if (group_start[0] != 0 || group_end[group_end.n_elem - 1] + 1 != X.n_cols) {
    Rcpp::stop("Groups must cover all transformed predictors contiguously.");
  }
  for (arma::uword g = 0; g < group_start.n_elem; ++g) {
    if (group_start[g] > group_end[g] || group_end[g] >= X.n_cols) {
      Rcpp::stop("Invalid group boundary.");
    }
    if (g > 0 && group_start[g] != group_end[g - 1] + 1) {
      Rcpp::stop("Group boundaries must be contiguous.");
    }
    if (!(group_weight[g] > 0.0) || !std::isfinite(group_weight[g])) {
      Rcpp::stop("All group weights must be finite and positive.");
    }
  }
}

}  // namespace


// [[Rcpp::export]]
Rcpp::List lsg_fisher_target_cpp(const arma::mat& X,
                                 const arma::vec& y,
                                 const arma::uvec& group_start,
                                 const arma::uvec& group_end) {
  if (X.n_rows != y.n_elem) Rcpp::stop("X and y have incompatible dimensions.");
  const double prevalence = arma::mean(y);
  if (!(prevalence > 0.0 && prevalence < 1.0)) {
    Rcpp::stop("Both outcome classes are required.");
  }
  const double n = static_cast<double>(X.n_rows);
  const double null_weight = prevalence * (1.0 - prevalence);
  const double null_intercept = std::log(prevalence / (1.0 - prevalence));
  const arma::vec score = X.t() * (y - prevalence) / n;
  arma::vec target(X.n_cols, arma::fill::zeros);

  for (arma::uword g = 0; g < group_start.n_elem; ++g) {
    const arma::uword first = group_start[g];
    const arma::uword last = group_end[g];
    const arma::mat Xg = X.cols(first, last);
    const arma::mat information = null_weight * (Xg.t() * Xg / n);
    arma::vec tg;
    const bool solved = arma::solve(
      tg, information, score.subvec(first, last),
      arma::solve_opts::likely_sympd + arma::solve_opts::no_approx
    );
    if (!solved || !tg.is_finite()) {
      tg = arma::pinv(information) * score.subvec(first, last);
    }
    target.subvec(first, last) = tg;
  }

  return Rcpp::List::create(
    Rcpp::_["target"] = target,
    Rcpp::_["score"] = score,
    Rcpp::_["prevalence"] = prevalence,
    Rcpp::_["null_weight"] = null_weight,
    Rcpp::_["null_intercept"] = null_intercept
  );
}


// [[Rcpp::export]]
Rcpp::List lsg_lambda_start_cpp(const arma::mat& X,
                                const arma::vec& y,
                                const arma::uvec& group_start,
                                const arma::uvec& group_end,
                                const arma::vec& group_weight,
                                const arma::vec& target,
                                const double alpha,
                                const double d) {
  if (!(alpha >= 0.0 && alpha <= 1.0)) {
    Rcpp::stop("alpha must lie in [0, 1].");
  }
  validate_groups(X, group_start, group_end, group_weight, target);
  const double prevalence = arma::mean(y);
  if (!(prevalence > 0.0 && prevalence < 1.0)) {
    Rcpp::stop("Both outcome classes are required.");
  }
  const double n = static_cast<double>(X.n_rows);
  const arma::vec score = X.t() * (y - prevalence) / n;
  const arma::uword groups = group_start.n_elem;
  arma::vec lower(groups, arma::fill::zeros);
  arma::vec upper(groups);
  upper.fill(R_PosInf);
  Rcpp::LogicalVector feasible(groups, true);
  double common_lower = 0.0;
  double common_upper = R_PosInf;

  for (arma::uword g = 0; g < groups; ++g) {
    const arma::uword first = group_start[g];
    const arma::uword last = group_end[g];
    const arma::vec sg = score.subvec(first, last);
    const arma::vec cg =
      (1.0 - alpha) * group_weight[g] * d * target.subvec(first, last);
    const double threshold_slope = alpha * group_weight[g];
    const double A = arma::dot(cg, cg) - threshold_slope * threshold_slope;
    const double B = 2.0 * arma::dot(sg, cg);
    const double C = arma::dot(sg, sg);
    const FeasibleInterval interval = nonnegative_quadratic_interval(A, B, C);
    feasible[g] = interval.feasible;
    lower[g] = interval.lower;
    upper[g] = interval.upper;
    if (interval.feasible) {
      common_lower = std::max(common_lower, interval.lower);
      common_upper = std::min(common_upper, interval.upper);
    }
  }

  bool common_feasible = true;
  for (arma::uword g = 0; g < groups; ++g) {
    if (!feasible[g]) common_feasible = false;
  }
  if (common_lower > common_upper + 1e-10 * (1.0 + std::abs(common_upper))) {
    common_feasible = false;
  }

  return Rcpp::List::create(
    Rcpp::_["zero_model_feasible"] = common_feasible,
    Rcpp::_["lambda_start"] = common_feasible ? common_lower : NA_REAL,
    Rcpp::_["common_upper"] = common_feasible ? common_upper : NA_REAL,
    Rcpp::_["group_lower"] = lower,
    Rcpp::_["group_upper"] = upper,
    Rcpp::_["group_feasible"] = feasible,
    Rcpp::_["score"] = score,
    Rcpp::_["alpha"] = alpha,
    Rcpp::_["d"] = d
  );
}


// [[Rcpp::export]]
Rcpp::List lsg_kkt_cpp(const arma::mat& X,
                       const arma::vec& y,
                       const arma::vec& beta,
                       const double intercept,
                       const arma::uvec& group_start,
                       const arma::uvec& group_end,
                       const arma::vec& group_weight,
                       const arma::vec& target,
                       const double lambda,
                       const double alpha,
                       const double d) {
  validate_groups(X, group_start, group_end, group_weight, target);
  arma::vec by_group;
  double intercept_residual = 0.0;
  const double maximum = kkt_residual_impl(
    X, y, beta, intercept, group_start, group_end, group_weight,
    target, lambda, alpha, d, &by_group, &intercept_residual
  );
  return Rcpp::List::create(
    Rcpp::_["maximum"] = maximum,
    Rcpp::_["intercept"] = intercept_residual,
    Rcpp::_["by_group"] = by_group
  );
}


// [[Rcpp::export]]
double lsg_objective_cpp(const arma::mat& X,
                         const arma::vec& y,
                         const arma::vec& beta,
                         const double intercept,
                         const arma::uvec& group_start,
                         const arma::uvec& group_end,
                         const arma::vec& group_weight,
                         const arma::vec& target,
                         const double lambda,
                         const double alpha,
                         const double d) {
  validate_groups(X, group_start, group_end, group_weight, target);
  return objective_impl(
    X, y, beta, intercept, group_start, group_end, group_weight,
    target, lambda, alpha, d
  );
}


// [[Rcpp::export]]
Rcpp::List lsg_path_cpp(const arma::mat& X,
                        const arma::vec& y,
                        const arma::uvec& group_start,
                        const arma::uvec& group_end,
                        const arma::vec& group_weight,
                        const arma::vec& target,
                        const arma::vec& lambda,
                        const arma::vec& d,
                        const double alpha,
                        const int max_outer = 250,
                        const int max_inner = 2000,
                        const double tolerance = 1e-7,
                        const double inner_tolerance = 1e-9,
                        const bool use_active_set = true,
                        const bool warm_start_d = true,
                        const bool keep_traces = false) {
  if (X.n_rows != y.n_elem) Rcpp::stop("X and y have incompatible dimensions.");
  if (!(alpha >= 0.0 && alpha <= 1.0)) Rcpp::stop("alpha must lie in [0, 1].");
  if (lambda.n_elem == 0 || d.n_elem == 0) Rcpp::stop("lambda and d cannot be empty.");
  if (arma::any(lambda <= 0.0)) Rcpp::stop("All lambda values must be positive.");
  if (arma::any(d < 0.0) || arma::any(d > 1.0)) Rcpp::stop("d must lie in [0, 1].");
  validate_groups(X, group_start, group_end, group_weight, target);

  const double prevalence = arma::mean(y);
  if (!(prevalence > 0.0 && prevalence < 1.0)) {
    Rcpp::stop("Both outcome classes are required.");
  }
  const double null_intercept = std::log(prevalence / (1.0 - prevalence));
  const arma::uword p = X.n_cols;
  const arma::uword L = lambda.n_elem;
  const arma::uword D = d.n_elem;

  arma::cube beta_path(p, L, D, arma::fill::zeros);
  arma::mat intercept_path(L, D, arma::fill::zeros);
  arma::mat objective_path(L, D, arma::fill::zeros);
  arma::mat kkt_path(L, D, arma::fill::zeros);
  arma::mat intercept_kkt_path(L, D, arma::fill::zeros);
  arma::imat outer_iterations(L, D, arma::fill::zeros);
  arma::imat inner_iterations(L, D, arma::fill::zeros);
  arma::umat converged(L, D, arma::fill::zeros);
  arma::mat max_raw_increase(L, D, arma::fill::zeros);
  arma::mat selected_groups(L, D, arma::fill::zeros);
  arma::imat active_set_scans(L, D, arma::fill::zeros);
  arma::mat group_updates(L, D, arma::fill::zeros);
  Rcpp::List traces(L * D);

  for (arma::uword di = 0; di < D; ++di) {
    arma::vec beta(p, arma::fill::zeros);
    double intercept = null_intercept;
    if (warm_start_d && di > 0u) {
      beta = beta_path.slice(di - 1u).col(0u);
      intercept = intercept_path(0u, di - 1u);
    }
    for (arma::uword li = 0; li < L; ++li) {
      OneFit fit = fit_one(
        X, y, group_start, group_end, group_weight, target,
        lambda[li], alpha, d[di], beta, intercept,
        max_outer, max_inner, tolerance, inner_tolerance, use_active_set,
        keep_traces
      );
      beta = fit.beta;
      intercept = fit.intercept;
      beta_path.slice(di).col(li) = beta;
      intercept_path(li, di) = intercept;
      objective_path(li, di) = fit.objective;
      kkt_path(li, di) = fit.kkt;
      intercept_kkt_path(li, di) = fit.intercept_kkt;
      outer_iterations(li, di) = fit.outer_iterations;
      inner_iterations(li, di) = fit.inner_iterations;
      converged(li, di) = fit.converged ? 1u : 0u;
      max_raw_increase(li, di) = fit.max_raw_objective_increase;
      active_set_scans(li, di) = fit.active_set_scans;
      group_updates(li, di) = fit.group_updates;
      int active = 0;
      for (arma::uword g = 0; g < group_start.n_elem; ++g) {
        if (arma::norm(beta.subvec(group_start[g], group_end[g]), 2) > 1e-8) {
          ++active;
        }
      }
      selected_groups(li, di) = active;
      if (keep_traces) {
        traces[di * L + li] = Rcpp::wrap(fit.objective_trace);
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::_["beta"] = beta_path,
    Rcpp::_["intercept"] = intercept_path,
    Rcpp::_["objective"] = objective_path,
    Rcpp::_["kkt"] = kkt_path,
    Rcpp::_["intercept_kkt"] = intercept_kkt_path,
    Rcpp::_["outer_iterations"] = outer_iterations,
    Rcpp::_["inner_iterations"] = inner_iterations,
    Rcpp::_["converged"] = converged,
    Rcpp::_["max_raw_objective_increase"] = max_raw_increase,
    Rcpp::_["selected_groups"] = selected_groups,
    Rcpp::_["active_set_scans"] = active_set_scans,
    Rcpp::_["group_updates"] = group_updates,
    Rcpp::_["use_active_set"] = use_active_set,
    Rcpp::_["warm_start_d"] = warm_start_d,
    Rcpp::_["objective_traces"] = traces,
    Rcpp::_["majorization_curvature"] = kMajorizationCurvature
  );
}


// [[Rcpp::export]]
Rcpp::List lsg_path_abgd_cpp(const arma::mat& X,
                             const arma::vec& y,
                             const arma::uvec& group_start,
                             const arma::uvec& group_end,
                             const arma::vec& group_weight,
                             const arma::vec& target,
                             const arma::vec& lambda,
                             const arma::vec& d,
                             const double alpha,
                             const int max_passes = 250,
                             const double tolerance = 1e-7,
                             const double update_tolerance = 1e-9,
                             const bool use_active_set = true,
                             const bool warm_start_d = true,
                             const bool keep_traces = false) {
  if (X.n_rows != y.n_elem) Rcpp::stop("X and y have incompatible dimensions.");
  if (!(alpha >= 0.0 && alpha <= 1.0)) Rcpp::stop("alpha must lie in [0, 1].");
  if (lambda.n_elem == 0 || d.n_elem == 0) Rcpp::stop("lambda and d cannot be empty.");
  if (arma::any(lambda <= 0.0)) Rcpp::stop("All lambda values must be positive.");
  if (arma::any(d < 0.0) || arma::any(d > 1.0)) Rcpp::stop("d must lie in [0, 1].");
  if (max_passes < 1) Rcpp::stop("max_passes must be positive.");
  if (!(tolerance > 0.0) || !(update_tolerance > 0.0)) {
    Rcpp::stop("Solver tolerances must be positive.");
  }
  validate_groups(X, group_start, group_end, group_weight, target);

  const double prevalence = arma::mean(y);
  if (!(prevalence > 0.0 && prevalence < 1.0)) {
    Rcpp::stop("Both outcome classes are required.");
  }
  const double null_intercept = std::log(prevalence / (1.0 - prevalence));
  const arma::uword p = X.n_cols;
  const arma::uword L = lambda.n_elem;
  const arma::uword D = d.n_elem;
  const arma::vec group_curvature = group_curvature_impl(
    X, group_start, group_end
  );

  arma::cube beta_path(p, L, D, arma::fill::zeros);
  arma::mat intercept_path(L, D, arma::fill::zeros);
  arma::mat objective_path(L, D, arma::fill::zeros);
  arma::mat kkt_path(L, D, arma::fill::zeros);
  arma::mat intercept_kkt_path(L, D, arma::fill::zeros);
  arma::imat passes(L, D, arma::fill::zeros);
  arma::imat inner_iterations(L, D, arma::fill::zeros);
  arma::umat converged(L, D, arma::fill::zeros);
  arma::mat max_raw_increase(L, D, arma::fill::zeros);
  arma::mat selected_groups(L, D, arma::fill::zeros);
  arma::imat kkt_scans(L, D, arma::fill::zeros);
  arma::mat group_updates(L, D, arma::fill::zeros);
  Rcpp::List traces(L * D);

  for (arma::uword di = 0; di < D; ++di) {
    arma::vec beta(p, arma::fill::zeros);
    double intercept = null_intercept;
    arma::uvec active(group_start.n_elem, arma::fill::zeros);
    if (warm_start_d && di > 0u) {
      beta = beta_path.slice(di - 1u).col(0u);
      intercept = intercept_path(0u, di - 1u);
      for (arma::uword g = 0; g < group_start.n_elem; ++g) {
        if (arma::norm(beta.subvec(group_start[g], group_end[g]), 2) > 1e-12) {
          active[g] = 1u;
        }
      }
    }

    for (arma::uword li = 0; li < L; ++li) {
      ABGDFit fit = fit_one_abgd(
        X, y, group_start, group_end, group_weight, target, group_curvature,
        lambda[li], alpha, d[di], beta, intercept, active,
        max_passes, tolerance, update_tolerance, use_active_set, keep_traces
      );
      beta = fit.beta;
      intercept = fit.intercept;
      active = fit.active;
      beta_path.slice(di).col(li) = beta;
      intercept_path(li, di) = intercept;
      objective_path(li, di) = fit.objective;
      kkt_path(li, di) = fit.kkt;
      intercept_kkt_path(li, di) = fit.intercept_kkt;
      passes(li, di) = fit.passes;
      converged(li, di) = fit.converged ? 1u : 0u;
      max_raw_increase(li, di) = fit.max_raw_objective_increase;
      kkt_scans(li, di) = fit.kkt_scans;
      group_updates(li, di) = fit.group_updates;
      int selected = 0;
      for (arma::uword g = 0; g < group_start.n_elem; ++g) {
        if (arma::norm(beta.subvec(group_start[g], group_end[g]), 2) > 1e-8) {
          ++selected;
        }
      }
      selected_groups(li, di) = selected;
      if (keep_traces) {
        traces[di * L + li] = Rcpp::wrap(fit.objective_trace);
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::_["beta"] = beta_path,
    Rcpp::_["intercept"] = intercept_path,
    Rcpp::_["objective"] = objective_path,
    Rcpp::_["kkt"] = kkt_path,
    Rcpp::_["intercept_kkt"] = intercept_kkt_path,
    Rcpp::_["outer_iterations"] = passes,
    Rcpp::_["inner_iterations"] = inner_iterations,
    Rcpp::_["converged"] = converged,
    Rcpp::_["max_raw_objective_increase"] = max_raw_increase,
    Rcpp::_["selected_groups"] = selected_groups,
    Rcpp::_["active_set_scans"] = kkt_scans,
    Rcpp::_["group_updates"] = group_updates,
    Rcpp::_["use_active_set"] = use_active_set,
    Rcpp::_["warm_start_d"] = warm_start_d,
    Rcpp::_["objective_traces"] = traces,
    Rcpp::_["majorization_curvature"] = kMajorizationCurvature,
    Rcpp::_["passes"] = passes,
    Rcpp::_["kkt_scans"] = kkt_scans,
    Rcpp::_["group_curvature"] = group_curvature,
    Rcpp::_["solver"] = "abgd"
  );
}
