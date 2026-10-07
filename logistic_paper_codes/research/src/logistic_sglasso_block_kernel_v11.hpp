#ifndef LOGISTIC_SGLASSO_BLOCK_KERNEL_V11_HPP
#define LOGISTIC_SGLASSO_BLOCK_KERNEL_V11_HPP

// Immutable numerical snapshot of the validated V9 profiled block kernel.
// V11 keeps the public V9 files unchanged and validates equivalence through
// objective/KKT reconstruction and frozen-source manifests.

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]

// V9 numerical-repair solver.  The statistical objective is identical to the
// frozen Logistic SGLASSO core.  This translation unit is deliberately
// separate so V7/V8 source identities remain immutable.

namespace {

inline double stable_expit_scalar(const double value) {
  if (value >= 0.0) {
    const double z = std::exp(-value);
    return 1.0 / (1.0 + z);
  }
  const double z = std::exp(value);
  return z / (1.0 + z);
}

inline double stable_log1pexp_scalar(const double value) {
  if (value > 0.0) return value + std::log1p(std::exp(-value));
  return std::log1p(std::exp(value));
}

arma::vec stable_expit(const arma::vec& eta) {
  arma::vec probability(eta.n_elem);
  for (arma::uword i = 0; i < eta.n_elem; ++i) {
    probability[i] = stable_expit_scalar(eta[i]);
  }
  return probability;
}

double logistic_loss(const arma::vec& eta, const arma::vec& y) {
  double value = 0.0;
  for (arma::uword i = 0; i < eta.n_elem; ++i) {
    value += stable_log1pexp_scalar(eta[i]) - y[i] * eta[i];
  }
  return value / static_cast<double>(eta.n_elem);
}

double group_penalty(const arma::vec& beta_group,
                     const arma::vec& target_group,
                     const double group_weight,
                     const double lambda,
                     const double alpha,
                     const double d) {
  const double lambda1 = lambda * alpha * group_weight;
  const double lambda2 = lambda * (1.0 - alpha) * group_weight;
  const arma::vec shifted = beta_group - d * target_group;
  return lambda1 * arma::norm(beta_group, 2) +
    0.5 * lambda2 * arma::dot(shifted, shifted);
}

double penalty(const arma::vec& beta,
               const arma::uvec& group_start,
               const arma::uvec& group_end,
               const arma::vec& group_weight,
               const arma::vec& target,
               const double lambda,
               const double alpha,
               const double d) {
  double value = 0.0;
  for (arma::uword g = 0; g < group_start.n_elem; ++g) {
    value += group_penalty(
      beta.subvec(group_start[g], group_end[g]),
      target.subvec(group_start[g], group_end[g]),
      group_weight[g], lambda, alpha, d
    );
  }
  return value;
}

void validate_groups(const arma::mat& X,
                     const arma::uvec& group_start,
                     const arma::uvec& group_end,
                     const arma::vec& group_weight,
                     const arma::vec& target) {
  if (X.n_cols == 0 || X.n_rows < 2 || target.n_elem != X.n_cols ||
      group_start.n_elem == 0 || group_start.n_elem != group_end.n_elem ||
      group_start.n_elem != group_weight.n_elem || !X.is_finite() ||
      !target.is_finite() || !group_weight.is_finite() ||
      arma::any(group_weight <= 0.0)) {
    Rcpp::stop("Invalid V9 design, target, or group metadata.");
  }
  arma::uword expected = 0u;
  for (arma::uword g = 0; g < group_start.n_elem; ++g) {
    if (group_start[g] != expected || group_end[g] < group_start[g] ||
        group_end[g] >= X.n_cols) {
      Rcpp::stop("V9 groups must be ordered, contiguous, and exhaustive.");
    }
    expected = group_end[g] + 1u;
  }
  if (expected != X.n_cols) {
    Rcpp::stop("V9 groups do not cover every solver-coordinate predictor.");
  }
}

void validate_response(const arma::vec& y, const arma::uword n) {
  if (y.n_elem != n || !y.is_finite()) {
    Rcpp::stop("X and y have incompatible dimensions or nonfinite values.");
  }
  bool zero = false;
  bool one = false;
  for (arma::uword i = 0; i < y.n_elem; ++i) {
    if (y[i] == 0.0) zero = true;
    else if (y[i] == 1.0) one = true;
    else Rcpp::stop("y must contain only numeric zero and one.");
  }
  if (!zero || !one) Rcpp::stop("Both response classes are required.");
}

struct InterceptResult {
  double intercept;
  double score;
  int iterations;
  bool converged;
};

struct InterceptEvaluation {
  double score;
  double curvature;
};

InterceptEvaluation evaluate_intercept(const arma::vec& offset,
                                       const arma::vec& y,
                                       const double intercept) {
  double score = 0.0;
  double curvature = 0.0;
  for (arma::uword i = 0; i < offset.n_elem; ++i) {
    const double eta = intercept + offset[i];
    if (!std::isfinite(eta)) {
      return InterceptEvaluation{
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()
      };
    }
    const double probability = stable_expit_scalar(eta);
    score += probability - y[i];
    curvature += probability * (1.0 - probability);
  }
  const double n = static_cast<double>(offset.n_elem);
  return InterceptEvaluation{score / n, curvature / n};
}

InterceptResult profile_intercept(const arma::vec& offset,
                                  const arma::vec& y,
                                  const double initial_intercept,
                                  const double score_tolerance,
                                  const int max_iterations) {
  const double prevalence = arma::mean(y);
  const double null_logit = std::log(prevalence) - std::log1p(-prevalence);
  double lower = null_logit - offset.max();
  double upper = null_logit - offset.min();
  if (!std::isfinite(lower) || !std::isfinite(upper) || lower > upper) {
    return InterceptResult{initial_intercept,
      std::numeric_limits<double>::infinity(), 0, false};
  }

  InterceptEvaluation lower_eval = evaluate_intercept(offset, y, lower);
  InterceptEvaluation upper_eval = evaluate_intercept(offset, y, upper);
  double expansion = std::max(1.0, upper - lower);
  for (int expansion_index = 0;
       expansion_index < 64 &&
         (lower_eval.score > 0.0 || upper_eval.score < 0.0);
       ++expansion_index) {
    if (lower_eval.score > 0.0) lower -= expansion;
    if (upper_eval.score < 0.0) upper += expansion;
    if (!std::isfinite(lower) || !std::isfinite(upper)) break;
    lower_eval = evaluate_intercept(offset, y, lower);
    upper_eval = evaluate_intercept(offset, y, upper);
    expansion *= 2.0;
  }
  if (!std::isfinite(lower_eval.score) || !std::isfinite(upper_eval.score) ||
      lower_eval.score > 0.0 || upper_eval.score < 0.0) {
    return InterceptResult{initial_intercept,
      std::numeric_limits<double>::infinity(), 0, false};
  }

  double intercept = std::min(upper, std::max(lower, initial_intercept));
  if (!std::isfinite(intercept)) intercept = lower + 0.5 * (upper - lower);
  InterceptEvaluation evaluation = evaluate_intercept(offset, y, intercept);
  int iterations = 0;
  bool converged = false;

  for (int iteration = 0; iteration < max_iterations; ++iteration) {
    iterations = iteration + 1;
    evaluation = evaluate_intercept(offset, y, intercept);
    if (!std::isfinite(evaluation.score) ||
        !std::isfinite(evaluation.curvature)) break;
    if (std::abs(evaluation.score) <= score_tolerance) {
      converged = true;
      break;
    }
    if (evaluation.score < 0.0) lower = intercept;
    else upper = intercept;
    const double width = upper - lower;
    if (!(width >= 0.0) || !std::isfinite(width)) break;
    if (width <= score_tolerance * (1.0 + std::abs(intercept))) {
      intercept = lower + 0.5 * width;
      evaluation = evaluate_intercept(offset, y, intercept);
      converged = std::isfinite(evaluation.score) &&
        std::abs(evaluation.score) <= 10.0 * score_tolerance;
      break;
    }
    double candidate = std::numeric_limits<double>::quiet_NaN();
    if (evaluation.curvature > std::numeric_limits<double>::min()) {
      candidate = intercept - evaluation.score / evaluation.curvature;
    }
    intercept = std::isfinite(candidate) && candidate > lower && candidate < upper ?
      candidate : lower + 0.5 * width;
  }
  evaluation = evaluate_intercept(offset, y, intercept);
  if (std::isfinite(evaluation.score) &&
      std::abs(evaluation.score) <= score_tolerance) converged = true;
  return InterceptResult{intercept, std::abs(evaluation.score),
                         iterations, converged};
}

struct KKTResult {
  double maximum;
  double intercept;
  double group_maximum;
  arma::vec group;
};

KKTResult kkt_residual(const arma::mat& X,
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
  const double n = static_cast<double>(X.n_rows);
  const arma::vec probability = stable_expit(intercept + X * beta);
  const arma::vec gradient = X.t() * (probability - y) / n;
  const double intercept_residual = std::abs(arma::mean(probability - y));
  arma::vec group_residual(group_start.n_elem, arma::fill::zeros);
  double group_maximum = 0.0;
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
      residual = arma::norm(
        gg + lambda2 * (bg - d * tg) + lambda1 * bg / beta_norm, 2
      );
    } else {
      residual = std::max(
        0.0, arma::norm(gg - lambda2 * d * tg, 2) - lambda1
      );
    }
    group_residual[g] = residual;
    group_maximum = std::max(group_maximum, residual);
  }
  return KKTResult{
    std::max(intercept_residual, group_maximum),
    intercept_residual, group_maximum, group_residual
  };
}

arma::vec conservative_group_curvature(const arma::mat& X,
                                       const arma::uvec& group_start,
                                       const arma::uvec& group_end) {
  const double n = static_cast<double>(X.n_rows);
  arma::vec curvature(group_start.n_elem, arma::fill::zeros);
  for (arma::uword g = 0; g < group_start.n_elem; ++g) {
    const arma::mat Xg = X.cols(group_start[g], group_end[g]);
    arma::vec eigenvalues;
    const arma::mat gram = Xg.t() * Xg / n;
    double largest = 0.0;
    if (arma::eig_sym(eigenvalues, gram) && eigenvalues.is_finite() &&
        eigenvalues.n_elem > 0u) largest = eigenvalues.max();
    else largest = arma::norm(gram, 2);
    if (!std::isfinite(largest) || largest < 0.0) {
      Rcpp::stop("Unable to compute a finite V9 group curvature.");
    }
    curvature[g] = std::max(
      0.25 * largest, 100.0 * std::numeric_limits<double>::epsilon()
    );
  }
  return curvature;
}

struct JointFit {
  arma::vec beta;
  double intercept;
  double objective;
  double kkt;
  double intercept_kkt;
  double group_kkt;
  int sweeps;
  bool converged;
  std::string termination;
  int intercept_iterations;
  int backtracking_steps;
  int kkt_scans;
  double group_updates;
  double max_raw_objective_increase;
  double max_accepted_objective_increase;
  int active_groups;
  std::vector<double> objective_trace;
};

JointFit fit_one_joint(const arma::mat& X,
                       const arma::vec& y,
                       const arma::uvec& group_start,
                       const arma::uvec& group_end,
                       const arma::vec& group_weight,
                       const arma::vec& target,
                       const arma::vec& conservative_curvature,
                       const double lambda,
                       const double alpha,
                       const double d,
                       arma::vec beta,
                       double intercept,
                       const int max_sweeps,
                       const double kkt_tolerance,
                       const double update_tolerance,
                       const double intercept_tolerance,
                       const int max_intercept_iterations,
                       const bool use_active_set,
                       const bool keep_trace) {
  const double n = static_cast<double>(X.n_rows);
  const arma::uword groups = group_start.n_elem;
  arma::vec adaptive_curvature = conservative_curvature;
  arma::uvec active(groups, arma::fill::ones);
  if (use_active_set) active.zeros();

  InterceptResult profile = profile_intercept(
    X * beta, y, intercept, intercept_tolerance, max_intercept_iterations
  );
  JointFit result;
  result.beta = beta;
  result.intercept = profile.intercept;
  result.objective = std::numeric_limits<double>::infinity();
  result.kkt = std::numeric_limits<double>::infinity();
  result.intercept_kkt = profile.score;
  result.group_kkt = std::numeric_limits<double>::infinity();
  result.sweeps = 0;
  result.converged = false;
  result.termination = profile.converged ? "maximum_sweeps" :
    "initial_intercept_failure";
  result.intercept_iterations = profile.iterations;
  result.backtracking_steps = 0;
  result.kkt_scans = 0;
  result.group_updates = 0.0;
  result.max_raw_objective_increase = 0.0;
  result.max_accepted_objective_increase = 0.0;
  result.active_groups = 0;

  if (!profile.converged) return result;
  intercept = profile.intercept;
  arma::vec eta = intercept + X * beta;
  arma::vec probability = stable_expit(eta);
  arma::vec residual = y - probability;
  double log_loss_value = logistic_loss(eta, y);
  double penalty_value = penalty(
    beta, group_start, group_end, group_weight, target,
    lambda, alpha, d
  );
  double objective_value = log_loss_value + penalty_value;
  if (keep_trace) result.objective_trace.push_back(objective_value);

  KKTResult kkt = kkt_residual(
    X, y, beta, intercept, group_start, group_end, group_weight,
    target, lambda, alpha, d
  );
  ++result.kkt_scans;
  for (arma::uword g = 0; g < groups; ++g) {
    if (!use_active_set || kkt.group[g] > kkt_tolerance ||
        arma::norm(beta.subvec(group_start[g], group_end[g]), 2) > 1e-12) {
      active[g] = 1u;
    }
  }
  if (kkt.maximum <= kkt_tolerance) {
    result.converged = true;
    result.termination = "kkt_tolerance_initial";
  }

  for (int sweep = 0; sweep < max_sweeps && !result.converged; ++sweep) {
    if (sweep % 10 == 0) Rcpp::checkUserInterrupt();
    const arma::vec beta_old = beta;
    const double intercept_old = intercept;
    const double objective_old = objective_value;
    bool block_failure = false;

    for (arma::uword g = 0; g < groups; ++g) {
      if (active[g] == 0u) continue;
      const arma::uword first = group_start[g];
      const arma::uword last = group_end[g];
      const arma::mat Xg = X.cols(first, last);
      const arma::vec beta_current = beta.subvec(first, last);
      const arma::vec target_group = target.subvec(first, last);
      const double old_group_penalty = group_penalty(
        beta_current, target_group, group_weight[g], lambda, alpha, d
      );
      const arma::vec score = Xg.t() * residual / n;
      const double lambda1 = lambda * alpha * group_weight[g];
      const double lambda2 = lambda * (1.0 - alpha) * group_weight[g];
      double trial_curvature = std::max(
        0.5 * adaptive_curvature[g],
        100.0 * std::numeric_limits<double>::epsilon()
      );
      bool accepted = false;

      for (int backtrack = 0; backtrack < 60; ++backtrack) {
        const arma::vec shifted_score =
          trial_curvature * beta_current + score +
          lambda2 * d * target_group;
        const double shifted_norm = arma::norm(shifted_score, 2);
        arma::vec beta_candidate(beta_current.n_elem, arma::fill::zeros);
        if (shifted_norm > lambda1 && shifted_norm > 0.0) {
          beta_candidate = ((1.0 - lambda1 / shifted_norm) /
            (trial_curvature + lambda2)) * shifted_score;
        }
        const arma::vec delta = beta_candidate - beta_current;
        const arma::vec eta_candidate = eta + Xg * delta;
        const double candidate_log_loss = logistic_loss(eta_candidate, y);
        const double candidate_group_penalty = group_penalty(
          beta_candidate, target_group, group_weight[g], lambda, alpha, d
        );
        const double candidate_objective = objective_value +
          (candidate_log_loss - log_loss_value) +
          (candidate_group_penalty - old_group_penalty);
        const double raw_increase = candidate_objective - objective_value;
        result.max_raw_objective_increase = std::max(
          result.max_raw_objective_increase, raw_increase
        );
        const double slack = 1e-12 * (1.0 + std::abs(objective_value));
        if (std::isfinite(candidate_objective) &&
            candidate_objective <= objective_value + slack) {
          result.max_accepted_objective_increase = std::max(
            result.max_accepted_objective_increase,
            candidate_objective - objective_value
          );
          beta.subvec(first, last) = beta_candidate;
          eta = eta_candidate;
          probability = stable_expit(eta);
          residual = y - probability;
          log_loss_value = candidate_log_loss;
          penalty_value += candidate_group_penalty - old_group_penalty;
          objective_value = candidate_objective;
          adaptive_curvature[g] = trial_curvature;
          result.group_updates += 1.0;
          accepted = true;
          break;
        }
        trial_curvature *= 2.0;
        ++result.backtracking_steps;
      }
      if (!accepted) {
        block_failure = true;
        break;
      }
    }

    if (block_failure) {
      beta = beta_old;
      intercept = intercept_old;
      result.termination = "group_backtracking_failure";
      break;
    }

    profile = profile_intercept(
      X * beta, y, intercept, intercept_tolerance,
      max_intercept_iterations
    );
    result.intercept_iterations += profile.iterations;
    if (!profile.converged) {
      beta = beta_old;
      intercept = intercept_old;
      result.termination = "intercept_failure";
      break;
    }
    intercept = profile.intercept;
    eta = intercept + X * beta;
    probability = stable_expit(eta);
    residual = y - probability;
    log_loss_value = logistic_loss(eta, y);
    penalty_value = penalty(
      beta, group_start, group_end, group_weight, target,
      lambda, alpha, d
    );
    const double profiled_objective = log_loss_value + penalty_value;
    // Each accepted block permits only round-off-sized slack.  The full-sweep
    // audit therefore scales that same allowance by the number of groups;
    // this is a numerical comparison tolerance, not an optimization step.
    const double pass_slack =
      (2.0 + static_cast<double>(groups)) * 1e-12 *
      (1.0 + std::abs(objective_old));
    if (!std::isfinite(profiled_objective) ||
        profiled_objective > objective_old + pass_slack) {
      beta = beta_old;
      intercept = intercept_old;
      result.termination = "profiled_objective_increase";
      break;
    }
    objective_value = profiled_objective;
    result.sweeps = sweep + 1;
    if (keep_trace) result.objective_trace.push_back(objective_value);

    kkt = kkt_residual(
      X, y, beta, intercept, group_start, group_end, group_weight,
      target, lambda, alpha, d
    );
    ++result.kkt_scans;
    bool activated = false;
    if (use_active_set) {
      for (arma::uword g = 0; g < groups; ++g) {
        if (active[g] == 0u && kkt.group[g] > kkt_tolerance) {
          active[g] = 1u;
          activated = true;
        }
      }
    }
    if (!activated && kkt.maximum <= kkt_tolerance) {
      result.converged = true;
      result.termination = "kkt_tolerance";
      break;
    }

    const double coefficient_change = std::max(
      std::abs(intercept - intercept_old),
      arma::abs(beta - beta_old).max()
    );
    const double coefficient_scale = 1.0 + std::max(
      std::abs(intercept), arma::abs(beta).max()
    );
    if (coefficient_change <= update_tolerance * coefficient_scale &&
        kkt.maximum > kkt_tolerance) {
      // Do not label a small update as convergence.  Reset conservative
      // curvatures once to escape a numerically timid block step.
      adaptive_curvature = conservative_curvature;
    }
  }

  eta = intercept + X * beta;
  objective_value = logistic_loss(eta, y) + penalty(
    beta, group_start, group_end, group_weight, target,
    lambda, alpha, d
  );
  kkt = kkt_residual(
    X, y, beta, intercept, group_start, group_end, group_weight,
    target, lambda, alpha, d
  );
  ++result.kkt_scans;
  if (kkt.maximum <= kkt_tolerance &&
      kkt.intercept <= std::max(intercept_tolerance, kkt_tolerance)) {
    result.converged = true;
    if (result.termination == "maximum_sweeps") {
      result.termination = "kkt_tolerance_final";
    }
  }

  result.beta = beta;
  result.intercept = intercept;
  result.objective = objective_value;
  result.kkt = kkt.maximum;
  result.intercept_kkt = kkt.intercept;
  result.group_kkt = kkt.group_maximum;
  result.active_groups = static_cast<int>(arma::accu(active));
  return result;
}

}  // namespace

#endif  // LOGISTIC_SGLASSO_BLOCK_KERNEL_V11_HPP
