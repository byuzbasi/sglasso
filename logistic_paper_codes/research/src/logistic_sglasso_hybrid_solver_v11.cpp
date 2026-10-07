#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]

// V11 uses a versioned snapshot of the immutable V9 block kernel. The V11 R
// loader adds logistic_prework/src to the compiler include path; manifests
// retain both the V9 source and this snapshot for provenance checks.
#include "logistic_sglasso_block_kernel_v11.hpp"

// V11 changes only the numerical algorithm for the already frozen convex
// Logistic SGLASSO objective. It first runs the V9 profiled block solver in
// bounded chunks while auditing KKT progress. A monotone accelerated proximal
// gradient fallback is used only when the block stage stalls, exhausts its
// budget, or encounters a numerical termination. The group proximal map below
// includes the shifted quadratic target exactly; it is not an approximation.

namespace {

struct BlockAudit {
  JointFit fit;
  bool stalled;
  int chunks;
  double initial_kkt;
  std::string stage_termination;
};

struct HybridFit {
  JointFit fit;
  int block_sweeps;
  int apg_iterations;
  int block_chunks;
  int block_backtracking_steps;
  int apg_backtracking_steps;
  int apg_restarts;
  bool fallback_used;
  bool block_stalled;
  double block_initial_kkt;
  double fallback_initial_kkt;
  double global_logistic_curvature;
  std::string block_termination;
  std::string solver_route;
};

arma::vec shifted_group_prox(
    const arma::vec& point,
    const arma::vec& gradient,
    const arma::uvec& group_start,
    const arma::uvec& group_end,
    const arma::vec& group_weight,
    const arma::vec& target,
    const double lambda,
    const double alpha,
    const double d,
    const double curvature) {
  arma::vec candidate(point.n_elem, arma::fill::zeros);
  for (arma::uword g = 0; g < group_start.n_elem; ++g) {
    const arma::uword first = group_start[g];
    const arma::uword last = group_end[g];
    const double lambda1 = lambda * alpha * group_weight[g];
    const double lambda2 = lambda * (1.0 - alpha) * group_weight[g];
    const arma::vec shifted_score =
      curvature * point.subvec(first, last) -
      gradient.subvec(first, last) +
      lambda2 * d * target.subvec(first, last);
    const double shifted_norm = arma::norm(shifted_score, 2);
    if (shifted_norm > lambda1 && shifted_norm > 0.0) {
      candidate.subvec(first, last) =
        ((1.0 - lambda1 / shifted_norm) /
         (curvature + lambda2)) * shifted_score;
    }
  }
  return candidate;
}

double global_logistic_curvature(const arma::mat& X) {
  const double n = static_cast<double>(X.n_rows);
  arma::mat gram;
  if (X.n_rows <= X.n_cols) gram = X * X.t() / n;
  else gram = X.t() * X / n;
  arma::vec eigenvalues;
  double largest = 0.0;
  if (arma::eig_sym(eigenvalues, gram) && eigenvalues.is_finite() &&
      eigenvalues.n_elem > 0u) {
    largest = std::max(0.0, eigenvalues.max());
  } else {
    largest = arma::norm(gram, 2);
  }
  if (!std::isfinite(largest) || largest < 0.0) {
    Rcpp::stop("Unable to compute a finite V11 global curvature.");
  }
  return std::max(
    0.25 * largest,
    100.0 * std::numeric_limits<double>::epsilon()
  );
}

void append_trace(std::vector<double>& destination,
                  const std::vector<double>& source) {
  if (source.empty()) return;
  std::size_t first = 0u;
  if (!destination.empty()) {
    const double scale = 1.0 + std::abs(destination.back());
    if (std::abs(destination.back() - source.front()) <= 5e-12 * scale) {
      first = 1u;
    }
  }
  destination.insert(destination.end(), source.begin() + first, source.end());
}

BlockAudit run_block_stage(
    const arma::mat& X,
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
    const int chunk_sweeps,
    const int stall_window,
    const double stall_relative_improvement,
    const double kkt_tolerance,
    const double update_tolerance,
    const double intercept_tolerance,
    const int max_intercept_iterations,
    const bool use_active_set,
    const bool keep_trace) {
  KKTResult initial = kkt_residual(
    X, y, beta, intercept, group_start, group_end, group_weight, target,
    lambda, alpha, d
  );
  const double initial_kkt = initial.maximum;
  double progress_reference = initial_kkt;
  int sweeps_since_progress = 0;
  int total_sweeps = 0;
  int total_intercept_iterations = 0;
  int total_backtracking = 0;
  int total_kkt_scans = 0;
  double total_group_updates = 0.0;
  double max_raw_increase = 0.0;
  double max_accepted_increase = 0.0;
  int chunks = 0;
  bool stalled = false;
  std::vector<double> trace;

  // max_sweeps is validated as positive, so the loop initializes last before
  // any read of its fields.
  JointFit last;

  while (total_sweeps < max_sweeps) {
    Rcpp::checkUserInterrupt();
    const int budget = std::min(chunk_sweeps, max_sweeps - total_sweeps);
    last = fit_one_joint(
      X, y, group_start, group_end, group_weight, target,
      conservative_curvature, lambda, alpha, d, beta, intercept, budget,
      kkt_tolerance, update_tolerance, intercept_tolerance,
      max_intercept_iterations, use_active_set, keep_trace
    );
    ++chunks;
    total_sweeps += last.sweeps;
    total_intercept_iterations += last.intercept_iterations;
    total_backtracking += last.backtracking_steps;
    total_kkt_scans += last.kkt_scans;
    total_group_updates += last.group_updates;
    max_raw_increase = std::max(
      max_raw_increase, last.max_raw_objective_increase
    );
    max_accepted_increase = std::max(
      max_accepted_increase, last.max_accepted_objective_increase
    );
    if (keep_trace) append_trace(trace, last.objective_trace);
    beta = last.beta;
    intercept = last.intercept;

    if (last.converged) break;
    const bool ordinary_budget =
      last.termination == "maximum_sweeps" ||
      last.termination == "kkt_tolerance_final";
    if (!ordinary_budget || last.sweeps < budget) break;

    const double material_target = progress_reference *
      (1.0 - stall_relative_improvement);
    const bool material_progress =
      last.kkt <= material_target || last.kkt <= kkt_tolerance;
    if (material_progress) {
      progress_reference = last.kkt;
      sweeps_since_progress = 0;
    } else {
      sweeps_since_progress += last.sweeps;
    }
    if (total_sweeps >= chunk_sweeps &&
        sweeps_since_progress >= stall_window) {
      stalled = true;
      break;
    }
  }

  last.sweeps = total_sweeps;
  last.intercept_iterations = total_intercept_iterations;
  last.backtracking_steps = total_backtracking;
  last.kkt_scans = total_kkt_scans;
  last.group_updates = total_group_updates;
  last.max_raw_objective_increase = max_raw_increase;
  last.max_accepted_objective_increase = max_accepted_increase;
  if (keep_trace) last.objective_trace = trace;
  std::string termination = last.termination;
  if (!last.converged && stalled) termination = "kkt_progress_stall";
  else if (!last.converged && total_sweeps >= max_sweeps &&
           last.termination == "maximum_sweeps") {
    termination = "block_sweep_budget";
  }
  last.termination = termination;
  return BlockAudit{last, stalled, chunks, initial_kkt, termination};
}

struct ProximalTrial {
  arma::vec beta;
  double intercept;
  double logistic_loss;
  double objective;
  double curvature;
  int intercept_iterations;
  int backtracking_steps;
  bool accepted;
};

ProximalTrial proximal_trial(
    const arma::mat& X,
    const arma::vec& y,
    const arma::uvec& group_start,
    const arma::uvec& group_end,
    const arma::vec& group_weight,
    const arma::vec& target,
    const double lambda,
    const double alpha,
    const double d,
    const arma::vec& point,
    const double point_intercept,
    double curvature,
    const double intercept_tolerance,
    const int max_intercept_iterations,
    const bool require_monotone,
    const double monotone_reference) {
  InterceptResult point_profile = profile_intercept(
    X * point, y, point_intercept, intercept_tolerance,
    max_intercept_iterations
  );
  int intercept_iterations = point_profile.iterations;
  if (!point_profile.converged) {
    return ProximalTrial{
      point, point_intercept, std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(), curvature,
      intercept_iterations, 0, false
    };
  }
  const arma::vec point_eta = point_profile.intercept + X * point;
  const double point_loss = logistic_loss(point_eta, y);
  const arma::vec point_probability = stable_expit(point_eta);
  const arma::vec gradient =
    X.t() * (point_probability - y) / static_cast<double>(X.n_rows);

  for (int backtrack = 0; backtrack < 80; ++backtrack) {
    const arma::vec candidate = shifted_group_prox(
      point, gradient, group_start, group_end, group_weight, target,
      lambda, alpha, d, curvature
    );
    const arma::vec delta = candidate - point;
    InterceptResult candidate_profile = profile_intercept(
      X * candidate, y, point_profile.intercept, intercept_tolerance,
      max_intercept_iterations
    );
    intercept_iterations += candidate_profile.iterations;
    if (candidate_profile.converged) {
      const double candidate_loss = logistic_loss(
        candidate_profile.intercept + X * candidate, y
      );
      const double majorizer = point_loss + arma::dot(gradient, delta) +
        0.5 * curvature * arma::dot(delta, delta);
      const double candidate_objective = candidate_loss + penalty(
        candidate, group_start, group_end, group_weight, target,
        lambda, alpha, d
      );
      const double smooth_slack = 5e-13 * (1.0 + std::abs(point_loss));
      const double objective_slack = 5e-13 *
        (1.0 + std::abs(monotone_reference));
      const bool majorized = std::isfinite(candidate_loss) &&
        candidate_loss <= majorizer + smooth_slack;
      const bool monotone = !require_monotone ||
        (std::isfinite(candidate_objective) &&
         candidate_objective <= monotone_reference + objective_slack);
      if (majorized && monotone) {
        return ProximalTrial{
          candidate, candidate_profile.intercept, candidate_loss,
          candidate_objective, curvature, intercept_iterations,
          backtrack, true
        };
      }
    }
    curvature *= 2.0;
    if (!std::isfinite(curvature)) break;
  }
  return ProximalTrial{
    point, point_profile.intercept, point_loss,
    std::numeric_limits<double>::infinity(), curvature,
    intercept_iterations, 80, false
  };
}

JointFit run_apg_stage(
    const arma::mat& X,
    const arma::vec& y,
    const arma::uvec& group_start,
    const arma::uvec& group_end,
    const arma::vec& group_weight,
    const arma::vec& target,
    const double lambda,
    const double alpha,
    const double d,
    const JointFit& start,
    const int max_iterations,
    const double kkt_tolerance,
    const double update_tolerance,
    const double intercept_tolerance,
    const int max_intercept_iterations,
    const int kkt_check_interval,
    const double global_curvature,
    const bool keep_trace,
    int& apg_backtracking_steps,
    int& apg_restarts) {
  arma::vec current = start.beta;
  arma::vec previous = current;
  double current_intercept = start.intercept;
  InterceptResult initial_profile = profile_intercept(
    X * current, y, current_intercept, intercept_tolerance,
    max_intercept_iterations
  );
  JointFit result = start;
  result.intercept_iterations = initial_profile.iterations;
  result.backtracking_steps = 0;
  result.kkt_scans = 0;
  result.group_updates = 0.0;
  result.sweeps = 0;
  result.max_raw_objective_increase = 0.0;
  result.max_accepted_objective_increase = 0.0;
  result.objective_trace.clear();
  result.converged = false;
  result.termination = "apg_max_iterations";
  apg_backtracking_steps = 0;
  apg_restarts = 0;
  if (!initial_profile.converged) {
    result.termination = "apg_initial_intercept_failure";
    return result;
  }
  current_intercept = initial_profile.intercept;
  double current_objective = logistic_loss(
    current_intercept + X * current, y
  ) + penalty(
    current, group_start, group_end, group_weight, target,
    lambda, alpha, d
  );
  if (keep_trace) result.objective_trace.push_back(current_objective);
  KKTResult kkt = kkt_residual(
    X, y, current, current_intercept, group_start, group_end,
    group_weight, target, lambda, alpha, d
  );
  ++result.kkt_scans;
  if (kkt.maximum <= kkt_tolerance) {
    result.converged = true;
    result.termination = "apg_kkt_tolerance_initial";
  }

  double acceleration = 1.0;
  double curvature = global_curvature;
  // At high accuracy the smooth decrease is below the absolute comparison
  // slack in proximal_trial(). Reducing L on that evidence can accept an
  // oversized step and destroy an already small KKT residual. The global
  // logistic bound remains valid after profiling the intercept. Use it as a
  // floor in this regime; the objective and proximal map are unchanged.
  const bool high_precision = kkt_tolerance < 1e-8;
  const double minimum_curvature = high_precision ? global_curvature :
    100.0 * std::numeric_limits<double>::epsilon();
  for (int iteration = 0;
       iteration < max_iterations && !result.converged;
       ++iteration) {
    if (iteration % 10 == 0) Rcpp::checkUserInterrupt();
    const double next_acceleration =
      0.5 * (1.0 + std::sqrt(1.0 + 4.0 * acceleration * acceleration));
    const double momentum = (acceleration - 1.0) / next_acceleration;
    arma::vec extrapolated = current + momentum * (current - previous);
    double extrapolated_intercept = current_intercept;
    ProximalTrial trial = proximal_trial(
      X, y, group_start, group_end, group_weight, target,
      lambda, alpha, d, extrapolated, extrapolated_intercept,
      std::max(0.8 * curvature, minimum_curvature),
      intercept_tolerance, max_intercept_iterations, false,
      current_objective
    );
    result.intercept_iterations += trial.intercept_iterations;
    apg_backtracking_steps += trial.backtracking_steps;
    if (!trial.accepted) {
      result.termination = "apg_backtracking_failure";
      break;
    }

    const double monotone_slack = 5e-13 *
      (1.0 + std::abs(current_objective));
    bool restarted = trial.objective > current_objective + monotone_slack;
    if (restarted) {
      ++apg_restarts;
      trial = proximal_trial(
        X, y, group_start, group_end, group_weight, target,
        lambda, alpha, d, current, current_intercept,
        std::max(trial.curvature, minimum_curvature),
        intercept_tolerance, max_intercept_iterations, true,
        current_objective
      );
      result.intercept_iterations += trial.intercept_iterations;
      apg_backtracking_steps += trial.backtracking_steps;
      if (!trial.accepted) {
        result.termination = "apg_monotone_restart_failure";
        break;
      }
    }

    const arma::vec old_current = current;
    const double old_objective = current_objective;
    previous = current;
    current = trial.beta;
    current_intercept = trial.intercept;
    current_objective = trial.objective;
    curvature = trial.curvature;
    acceleration = restarted ?
      0.5 * (1.0 + std::sqrt(5.0)) : next_acceleration;
    // Directional restart does not subtract nearly equal objective values.
    if (high_precision &&
        arma::dot(extrapolated - current, current - old_current) > 0.0) {
      previous = current;
      acceleration = 1.0;
      ++apg_restarts;
    }
    result.sweeps = iteration + 1;
    result.group_updates += static_cast<double>(group_start.n_elem);
    result.max_accepted_objective_increase = std::max(
      result.max_accepted_objective_increase,
      current_objective - old_objective
    );
    if (keep_trace) result.objective_trace.push_back(current_objective);

    const bool check_kkt = iteration == 0 ||
      ((iteration + 1) % kkt_check_interval == 0) ||
      iteration + 1 == max_iterations;
    if (check_kkt) {
      kkt = kkt_residual(
        X, y, current, current_intercept, group_start, group_end,
        group_weight, target, lambda, alpha, d
      );
      ++result.kkt_scans;
      if (kkt.maximum <= kkt_tolerance) {
        result.converged = true;
        result.termination = "apg_kkt_tolerance";
        break;
      }
    }

    const double coefficient_change = arma::abs(current - old_current).max();
    const double coefficient_scale = 1.0 + arma::abs(current).max();
    if (!high_precision &&
        coefficient_change <= update_tolerance * coefficient_scale &&
        !result.converged) {
      // A small step is not convergence. Restart acceleration and let local
      // backtracking lower the curvature on the next iteration.
      // High-precision mode uses the directional restart above: resetting on
      // every tiny step would suppress acceleration before reaching its KKT
      // target in weakly curved directions.
      previous = current;
      acceleration = 1.0;
      ++apg_restarts;
    }
  }

  InterceptResult final_profile = profile_intercept(
    X * current, y, current_intercept, intercept_tolerance,
    max_intercept_iterations
  );
  result.intercept_iterations += final_profile.iterations;
  if (final_profile.converged) current_intercept = final_profile.intercept;
  const double final_objective = logistic_loss(
    current_intercept + X * current, y
  ) + penalty(
    current, group_start, group_end, group_weight, target,
    lambda, alpha, d
  );
  kkt = kkt_residual(
    X, y, current, current_intercept, group_start, group_end,
    group_weight, target, lambda, alpha, d
  );
  ++result.kkt_scans;
  if (kkt.maximum <= kkt_tolerance && final_profile.converged) {
    result.converged = true;
    if (result.termination == "apg_max_iterations") {
      result.termination = "apg_kkt_tolerance_final";
    }
  }
  result.beta = current;
  result.intercept = current_intercept;
  result.objective = final_objective;
  result.kkt = kkt.maximum;
  result.intercept_kkt = kkt.intercept;
  result.group_kkt = kkt.group_maximum;
  int active_groups = 0;
  for (arma::uword g = 0; g < group_start.n_elem; ++g) {
    if (arma::norm(current.subvec(group_start[g], group_end[g]), 2) > 1e-12) {
      ++active_groups;
    }
  }
  result.active_groups = active_groups;
  result.backtracking_steps = apg_backtracking_steps;
  return result;
}

HybridFit fit_one_hybrid(
    const arma::mat& X,
    const arma::vec& y,
    const arma::uvec& group_start,
    const arma::uvec& group_end,
    const arma::vec& group_weight,
    const arma::vec& target,
    const arma::vec& conservative_curvature,
    const double global_curvature,
    const double lambda,
    const double alpha,
    const double d,
    const arma::vec& beta_initial,
    const double intercept_initial,
    const int block_max_sweeps,
    const int block_chunk_sweeps,
    const int block_stall_window,
    const double block_stall_relative_improvement,
    const int apg_max_iterations,
    const int apg_kkt_check_interval,
    const double kkt_tolerance,
    const double update_tolerance,
    const double intercept_tolerance,
    const int max_intercept_iterations,
    const bool use_active_set,
    const bool enable_fallback,
    const bool keep_trace) {
  // In high-precision mode the block stage is a warm-up, not the final fit.
  // Avoid spending thousands of sweeps near its round-off-sensitive floor;
  // the APG stage still has to satisfy the caller's original tight tolerance.
  const double block_tolerance = (enable_fallback && kkt_tolerance < 1e-8) ?
    std::max(kkt_tolerance, 2e-6) : kkt_tolerance;
  BlockAudit block = run_block_stage(
    X, y, group_start, group_end, group_weight, target,
    conservative_curvature, lambda, alpha, d, beta_initial,
    intercept_initial, block_max_sweeps, block_chunk_sweeps,
    block_stall_window, block_stall_relative_improvement,
    block_tolerance, update_tolerance, intercept_tolerance,
    max_intercept_iterations, use_active_set, keep_trace
  );
  HybridFit output;
  output.fit = block.fit;
  output.block_sweeps = block.fit.sweeps;
  output.apg_iterations = 0;
  output.block_chunks = block.chunks;
  output.block_backtracking_steps = block.fit.backtracking_steps;
  output.apg_backtracking_steps = 0;
  output.apg_restarts = 0;
  output.fallback_used = false;
  output.block_stalled = block.stalled;
  output.block_initial_kkt = block.initial_kkt;
  output.fallback_initial_kkt = block.fit.kkt;
  output.global_logistic_curvature = global_curvature;
  output.block_termination = block.stage_termination;
  const bool final_converged = block.fit.converged &&
    block.fit.kkt <= kkt_tolerance;
  output.fit.converged = final_converged;
  output.solver_route = final_converged ? "profiled_block" : "block_only";
  if (final_converged || !enable_fallback) return output;

  int apg_backtracking = 0;
  int apg_restarts = 0;
  JointFit fallback = run_apg_stage(
    X, y, group_start, group_end, group_weight, target,
    lambda, alpha, d, block.fit, apg_max_iterations,
    kkt_tolerance, update_tolerance, intercept_tolerance,
    max_intercept_iterations, apg_kkt_check_interval, global_curvature,
    keep_trace, apg_backtracking, apg_restarts
  );
  if (keep_trace) {
    std::vector<double> combined = block.fit.objective_trace;
    append_trace(combined, fallback.objective_trace);
    fallback.objective_trace = combined;
  }
  fallback.intercept_iterations += block.fit.intercept_iterations;
  fallback.kkt_scans += block.fit.kkt_scans;
  fallback.group_updates += block.fit.group_updates;
  fallback.max_raw_objective_increase = std::max(
    fallback.max_raw_objective_increase,
    block.fit.max_raw_objective_increase
  );
  fallback.max_accepted_objective_increase = std::max(
    fallback.max_accepted_objective_increase,
    block.fit.max_accepted_objective_increase
  );
  output.fit = fallback;
  output.apg_iterations = fallback.sweeps;
  output.apg_backtracking_steps = apg_backtracking;
  output.apg_restarts = apg_restarts;
  output.fallback_used = true;
  output.solver_route = "profiled_block_then_monotone_apg";
  return output;
}

void validate_hybrid_controls(
    const int block_max_sweeps,
    const int block_chunk_sweeps,
    const int block_stall_window,
    const double block_stall_relative_improvement,
    const int apg_max_iterations,
    const int apg_kkt_check_interval,
    const double kkt_tolerance,
    const double update_tolerance,
    const double intercept_tolerance,
    const int max_intercept_iterations) {
  if (block_max_sweeps < 1 || block_chunk_sweeps < 1 ||
      block_stall_window < 1 ||
      !(block_stall_relative_improvement > 0.0 &&
        block_stall_relative_improvement < 1.0) ||
      apg_max_iterations < 1 || apg_kkt_check_interval < 1 ||
      !(kkt_tolerance > 0.0) || !(update_tolerance > 0.0) ||
      !(intercept_tolerance > 0.0) || max_intercept_iterations < 1) {
    Rcpp::stop("Invalid V11 hybrid numerical controls.");
  }
}

Rcpp::List wrap_hybrid_fit(const HybridFit& value) {
  const JointFit& fit = value.fit;
  return Rcpp::List::create(
    Rcpp::_ ["schema_version"] = "logistic_sglasso_hybrid_fit_v11",
    Rcpp::_ ["beta"] = fit.beta,
    Rcpp::_ ["intercept"] = fit.intercept,
    Rcpp::_ ["objective"] = fit.objective,
    Rcpp::_ ["kkt"] = fit.kkt,
    Rcpp::_ ["intercept_kkt"] = fit.intercept_kkt,
    Rcpp::_ ["group_kkt"] = fit.group_kkt,
    Rcpp::_ ["converged"] = fit.converged,
    Rcpp::_ ["termination_reason"] = fit.termination,
    Rcpp::_ ["solver_route"] = value.solver_route,
    Rcpp::_ ["fallback_used"] = value.fallback_used,
    Rcpp::_ ["block_stalled"] = value.block_stalled,
    Rcpp::_ ["block_termination"] = value.block_termination,
    Rcpp::_ ["block_initial_kkt"] = value.block_initial_kkt,
    Rcpp::_ ["fallback_initial_kkt"] = value.fallback_initial_kkt,
    Rcpp::_ ["block_sweeps"] = value.block_sweeps,
    Rcpp::_ ["block_chunks"] = value.block_chunks,
    Rcpp::_ ["apg_iterations"] = value.apg_iterations,
    Rcpp::_ ["intercept_iterations"] = fit.intercept_iterations,
    Rcpp::_ ["block_backtracking_steps"] =
      value.block_backtracking_steps,
    Rcpp::_ ["apg_backtracking_steps"] = value.apg_backtracking_steps,
    Rcpp::_ ["apg_restarts"] = value.apg_restarts,
    Rcpp::_ ["kkt_scans"] = fit.kkt_scans,
    Rcpp::_ ["group_updates"] = fit.group_updates,
    Rcpp::_ ["max_raw_objective_increase"] =
      fit.max_raw_objective_increase,
    Rcpp::_ ["max_accepted_objective_increase"] =
      fit.max_accepted_objective_increase,
    Rcpp::_ ["active_groups"] = fit.active_groups,
    Rcpp::_ ["global_logistic_curvature"] =
      value.global_logistic_curvature,
    Rcpp::_ ["objective_trace"] = fit.objective_trace
  );
}

}  // namespace


// Exact V11 proximal map, exported for mathematical unit tests.
// [[Rcpp::export]]
arma::vec lsg_shifted_group_prox_v11_cpp(
    const arma::vec& point,
    const arma::vec& gradient,
    const arma::uvec& group_start,
    const arma::uvec& group_end,
    const arma::vec& group_weight,
    const arma::vec& target,
    const double lambda,
    const double alpha,
    const double d,
    const double curvature) {
  if (point.n_elem == 0 || point.n_elem != gradient.n_elem ||
      point.n_elem != target.n_elem || !point.is_finite() ||
      !gradient.is_finite() || !(lambda > 0.0) ||
      !(alpha >= 0.0 && alpha <= 1.0) ||
      !(d >= 0.0 && d <= 1.0) || !(curvature > 0.0)) {
    Rcpp::stop("Invalid V11 proximal-map input.");
  }
  arma::mat dummy(2, point.n_elem, arma::fill::zeros);
  validate_groups(dummy, group_start, group_end, group_weight, target);
  return shifted_group_prox(
    point, gradient, group_start, group_end, group_weight, target,
    lambda, alpha, d, curvature
  );
}


// Fit one Logistic SGLASSO point from an explicit starting value.
// [[Rcpp::export]]
Rcpp::List lsg_fit_one_hybrid_v11_cpp(
    const arma::mat& X,
    const arma::vec& y,
    const arma::uvec& group_start,
    const arma::uvec& group_end,
    const arma::vec& group_weight,
    const arma::vec& target,
    const double lambda,
    const double alpha,
    const double d,
    const arma::vec& beta_initial,
    const double intercept_initial,
    const int block_max_sweeps = 4000,
    const int block_chunk_sweeps = 100,
    const int block_stall_window = 200,
    const double block_stall_relative_improvement = 0.01,
    const int apg_max_iterations = 10000,
    const int apg_kkt_check_interval = 5,
    const double kkt_tolerance = 2e-6,
    const double update_tolerance = 1e-10,
    const double intercept_tolerance = 1e-12,
    const int max_intercept_iterations = 100,
    const bool use_active_set = true,
    const bool enable_fallback = true,
    const bool keep_trace = false) {
  validate_groups(X, group_start, group_end, group_weight, target);
  validate_response(y, X.n_rows);
  validate_hybrid_controls(
    block_max_sweeps, block_chunk_sweeps, block_stall_window,
    block_stall_relative_improvement, apg_max_iterations,
    apg_kkt_check_interval, kkt_tolerance, update_tolerance,
    intercept_tolerance, max_intercept_iterations
  );
  if (beta_initial.n_elem != X.n_cols || !beta_initial.is_finite() ||
      !std::isfinite(intercept_initial) || !(lambda > 0.0) ||
      !std::isfinite(lambda) || !(alpha >= 0.0 && alpha <= 1.0) ||
      !(d >= 0.0 && d <= 1.0)) {
    Rcpp::stop("Invalid V11 single-fit input or starting value.");
  }
  const arma::vec group_curvature = conservative_group_curvature(
    X, group_start, group_end
  );
  const double global_curvature = global_logistic_curvature(X);
  return wrap_hybrid_fit(fit_one_hybrid(
    X, y, group_start, group_end, group_weight, target,
    group_curvature, global_curvature, lambda, alpha, d,
    beta_initial, intercept_initial, block_max_sweeps,
    block_chunk_sweeps, block_stall_window,
    block_stall_relative_improvement, apg_max_iterations,
    apg_kkt_check_interval, kkt_tolerance, update_tolerance,
    intercept_tolerance, max_intercept_iterations, use_active_set,
    enable_fallback, keep_trace
  ));
}


// Fit a finite lambda-by-d path in the supplied order.
// [[Rcpp::export]]
Rcpp::List lsg_path_hybrid_v11_cpp(
    const arma::mat& X,
    const arma::vec& y,
    const arma::uvec& group_start,
    const arma::uvec& group_end,
    const arma::vec& group_weight,
    const arma::vec& target,
    const arma::vec& lambda,
    const arma::vec& d,
    const double alpha,
    const int block_max_sweeps = 4000,
    const int block_chunk_sweeps = 100,
    const int block_stall_window = 200,
    const double block_stall_relative_improvement = 0.01,
    const int apg_max_iterations = 10000,
    const int apg_kkt_check_interval = 5,
    const double kkt_tolerance = 2e-6,
    const double update_tolerance = 1e-10,
    const double intercept_tolerance = 1e-12,
    const int max_intercept_iterations = 100,
    const bool use_active_set = true,
    const bool enable_fallback = true,
    const bool warm_start_d = true,
    const bool keep_traces = false) {
  validate_groups(X, group_start, group_end, group_weight, target);
  validate_response(y, X.n_rows);
  validate_hybrid_controls(
    block_max_sweeps, block_chunk_sweeps, block_stall_window,
    block_stall_relative_improvement, apg_max_iterations,
    apg_kkt_check_interval, kkt_tolerance, update_tolerance,
    intercept_tolerance, max_intercept_iterations
  );
  if (lambda.n_elem == 0 || d.n_elem == 0 || !lambda.is_finite() ||
      !d.is_finite() || arma::any(lambda <= 0.0) || arma::any(d < 0.0) ||
      arma::any(d > 1.0) || !(alpha >= 0.0 && alpha <= 1.0)) {
    Rcpp::stop("Invalid V11 path input.");
  }

  const arma::uword p = X.n_cols;
  const arma::uword L = lambda.n_elem;
  const arma::uword D = d.n_elem;
  const double prevalence = arma::mean(y);
  const double null_intercept = std::log(prevalence) -
    std::log1p(-prevalence);
  const arma::vec group_curvature = conservative_group_curvature(
    X, group_start, group_end
  );
  const double global_curvature = global_logistic_curvature(X);

  arma::cube beta_path(p, L, D, arma::fill::zeros);
  arma::mat intercept_path(L, D, arma::fill::zeros);
  arma::mat objective_path(L, D, arma::fill::zeros);
  arma::mat kkt_path(L, D, arma::fill::zeros);
  arma::mat intercept_kkt_path(L, D, arma::fill::zeros);
  arma::mat group_kkt_path(L, D, arma::fill::zeros);
  arma::umat converged_path(L, D, arma::fill::zeros);
  arma::umat fallback_path(L, D, arma::fill::zeros);
  arma::umat stalled_path(L, D, arma::fill::zeros);
  arma::imat block_sweeps_path(L, D, arma::fill::zeros);
  arma::imat block_chunks_path(L, D, arma::fill::zeros);
  arma::imat apg_iterations_path(L, D, arma::fill::zeros);
  arma::imat intercept_iterations_path(L, D, arma::fill::zeros);
  arma::imat block_backtracking_path(L, D, arma::fill::zeros);
  arma::imat apg_backtracking_path(L, D, arma::fill::zeros);
  arma::imat apg_restarts_path(L, D, arma::fill::zeros);
  arma::imat kkt_scans_path(L, D, arma::fill::zeros);
  arma::mat group_updates_path(L, D, arma::fill::zeros);
  arma::mat block_initial_kkt_path(L, D, arma::fill::zeros);
  arma::mat fallback_initial_kkt_path(L, D, arma::fill::zeros);
  arma::mat raw_increase_path(L, D, arma::fill::zeros);
  arma::mat accepted_increase_path(L, D, arma::fill::zeros);
  arma::imat active_groups_path(L, D, arma::fill::zeros);
  Rcpp::CharacterMatrix termination_path(L, D);
  Rcpp::CharacterMatrix block_termination_path(L, D);
  Rcpp::CharacterMatrix route_path(L, D);
  Rcpp::List traces(L * D);

  for (arma::uword di = 0; di < D; ++di) {
    arma::vec beta(p, arma::fill::zeros);
    double intercept = null_intercept;
    if (warm_start_d && di > 0u) {
      beta = beta_path.slice(di - 1u).col(0u);
      intercept = intercept_path(0u, di - 1u);
    }
    for (arma::uword li = 0; li < L; ++li) {
      HybridFit value = fit_one_hybrid(
        X, y, group_start, group_end, group_weight, target,
        group_curvature, global_curvature, lambda[li], alpha, d[di],
        beta, intercept, block_max_sweeps, block_chunk_sweeps,
        block_stall_window, block_stall_relative_improvement,
        apg_max_iterations, apg_kkt_check_interval, kkt_tolerance,
        update_tolerance, intercept_tolerance, max_intercept_iterations,
        use_active_set, enable_fallback, keep_traces
      );
      const JointFit& fit = value.fit;
      beta = fit.beta;
      intercept = fit.intercept;
      beta_path.slice(di).col(li) = fit.beta;
      intercept_path(li, di) = fit.intercept;
      objective_path(li, di) = fit.objective;
      kkt_path(li, di) = fit.kkt;
      intercept_kkt_path(li, di) = fit.intercept_kkt;
      group_kkt_path(li, di) = fit.group_kkt;
      converged_path(li, di) = fit.converged ? 1u : 0u;
      fallback_path(li, di) = value.fallback_used ? 1u : 0u;
      stalled_path(li, di) = value.block_stalled ? 1u : 0u;
      block_sweeps_path(li, di) = value.block_sweeps;
      block_chunks_path(li, di) = value.block_chunks;
      apg_iterations_path(li, di) = value.apg_iterations;
      intercept_iterations_path(li, di) = fit.intercept_iterations;
      block_backtracking_path(li, di) =
        value.block_backtracking_steps;
      apg_backtracking_path(li, di) = value.apg_backtracking_steps;
      apg_restarts_path(li, di) = value.apg_restarts;
      kkt_scans_path(li, di) = fit.kkt_scans;
      group_updates_path(li, di) = fit.group_updates;
      block_initial_kkt_path(li, di) = value.block_initial_kkt;
      fallback_initial_kkt_path(li, di) = value.fallback_initial_kkt;
      raw_increase_path(li, di) = fit.max_raw_objective_increase;
      accepted_increase_path(li, di) =
        fit.max_accepted_objective_increase;
      active_groups_path(li, di) = fit.active_groups;
      termination_path(li, di) = fit.termination;
      block_termination_path(li, di) = value.block_termination;
      route_path(li, di) = value.solver_route;
      if (keep_traces) traces[di * L + li] = fit.objective_trace;
    }
  }

  return Rcpp::List::create(
    Rcpp::_ ["schema_version"] = "logistic_sglasso_hybrid_path_v11",
    Rcpp::_ ["beta"] = beta_path,
    Rcpp::_ ["intercept"] = intercept_path,
    Rcpp::_ ["objective"] = objective_path,
    Rcpp::_ ["kkt"] = kkt_path,
    Rcpp::_ ["intercept_kkt"] = intercept_kkt_path,
    Rcpp::_ ["group_kkt"] = group_kkt_path,
    Rcpp::_ ["converged"] = converged_path,
    Rcpp::_ ["termination_reason"] = termination_path,
    Rcpp::_ ["block_termination"] = block_termination_path,
    Rcpp::_ ["solver_route"] = route_path,
    Rcpp::_ ["fallback_used"] = fallback_path,
    Rcpp::_ ["block_stalled"] = stalled_path,
    Rcpp::_ ["block_initial_kkt"] = block_initial_kkt_path,
    Rcpp::_ ["fallback_initial_kkt"] = fallback_initial_kkt_path,
    Rcpp::_ ["block_sweeps"] = block_sweeps_path,
    Rcpp::_ ["block_chunks"] = block_chunks_path,
    Rcpp::_ ["apg_iterations"] = apg_iterations_path,
    Rcpp::_ ["intercept_iterations"] = intercept_iterations_path,
    Rcpp::_ ["block_backtracking_steps"] = block_backtracking_path,
    Rcpp::_ ["apg_backtracking_steps"] = apg_backtracking_path,
    Rcpp::_ ["apg_restarts"] = apg_restarts_path,
    Rcpp::_ ["kkt_scans"] = kkt_scans_path,
    Rcpp::_ ["group_updates"] = group_updates_path,
    Rcpp::_ ["max_raw_objective_increase"] = raw_increase_path,
    Rcpp::_ ["max_accepted_objective_increase"] = accepted_increase_path,
    Rcpp::_ ["active_groups"] = active_groups_path,
    Rcpp::_ ["group_curvature_upper_bound"] = group_curvature,
    Rcpp::_ ["global_logistic_curvature"] = global_curvature,
    Rcpp::_ ["objective_traces"] = traces,
    Rcpp::_ ["solver"] = "profiled_block_monotone_apg_v11",
    Rcpp::_ ["warm_start_d"] = warm_start_d,
    Rcpp::_ ["use_active_set"] = use_active_set,
    Rcpp::_ ["fallback_enabled"] = enable_fallback,
    Rcpp::_ ["kkt_tolerance"] = kkt_tolerance,
    Rcpp::_ ["intercept_tolerance"] = intercept_tolerance
  );
}
