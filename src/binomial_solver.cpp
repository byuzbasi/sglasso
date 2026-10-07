#include <RcppArmadillo.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]

// V11 uses a versioned snapshot of the immutable V9 block kernel. The V11 R
// loader adds logistic_prework/src to the compiler include path; manifests
// retain both the V9 source and this snapshot for provenance checks.
#include "binomial_block_kernel.h"

// Local-only APG/opt-in IRLS prototype; frozen V11 and ALL V5 are unchanged.
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

struct ApgTiming {
  double matvec_seconds;
  double intercept_seconds;
  double gradient_seconds;
  double prox_seconds;
  double kkt_seconds;
};

struct HybridFit {
  double diagnostic_irls_seconds = 0.0;
  int irls_iterations = 0;
  int irls_inner_sweeps = 0;
  std::string irls_termination = "not_used";
  JointFit fit;
  double diagnostic_block_seconds;
  double diagnostic_apg_seconds;
  ApgTiming diagnostic_apg_parts;
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

// Local opt-in proximal Newton implementation. The Hessian is used only to
// construct a search direction; acceptance and convergence always use the
// ORIGINAL logistic objective and KKT residual, including the Firth shift.
struct QuadraticGroup {
  arma::vec eigenvalue;
  arma::mat basis;
  bool ready = false;
};

bool prepare_quadratic_group(const arma::mat& A, QuadraticGroup& cache) {
  // A = sqrt(W) * centered X_g. Keep the null-space target component even
  // when group size exceeds n. No p-by-p full-design Hessian is constructed.
  if (A.n_cols <= A.n_rows) {
    if (!arma::eig_sym(cache.eigenvalue, cache.basis, A.t() * A)) return false;
    cache.eigenvalue = arma::clamp(cache.eigenvalue, 0.0,
      std::numeric_limits<double>::infinity());
  } else {
    arma::mat U;
    arma::vec s;
    if (!arma::svd_econ(U, s, cache.basis, A, "both", "std")) return false;
    cache.eigenvalue = arma::square(s);
  }
  cache.ready = cache.eigenvalue.is_finite() && cache.basis.is_finite();
  return cache.ready;
}

bool solve_quadratic_group(const QuadraticGroup& cache,
                          const arma::vec& score, const double lambda1,
                          const double lambda2, const double initial_norm,
                          arma::vec& beta) {
  const double score_norm = arma::norm(score, 2);
  beta.zeros(score.n_elem);
  if (score_norm <= lambda1) return true;
  const arma::vec projected = cache.basis.t() * score;
  // A full orthonormal basis has no complementary component. Avoid creating
  // a spurious null component from roundoff in that case.
  arma::vec complement(score.n_elem, arma::fill::zeros);
  if (cache.basis.n_cols < score.n_elem) {
    complement = score - cache.basis * projected;
  }
  const double complement_sq = arma::dot(complement, complement);
  const arma::vec diagonal = cache.eigenvalue + lambda2;
  if (lambda1 == 0.0) {
    if (!(lambda2 > 0.0)) return false;
    beta = cache.basis * (projected / diagonal) + complement / lambda2;
    return beta.is_finite();
  }
  // If r = ||beta_g||, the block stationarity equation gives coefficients
  // r*s_j/(r*(h_j+lambda2)+lambda1). Their norm must equal r. This scalar
  // decreasing equation is solved by safeguarded Newton, not by a scalar
  // curvature approximation to the group Hessian.
  const auto equation = [&](const double r, double* derivative) {
    const arma::vec denom = r * diagonal + lambda1;
    const arma::vec ratio = projected / denom;
    const double null_denom = r * lambda2 + lambda1;
    if (derivative) *derivative = -2.0 * (
      arma::accu(diagonal % arma::square(ratio) / denom) +
      lambda2 * complement_sq / (null_denom * null_denom * null_denom));
    return arma::dot(ratio, ratio) +
      complement_sq / (null_denom * null_denom) - 1.0;
  };
  double lower = 0.0;
  double upper = std::max(initial_norm,
    score_norm / std::max(1e-12, diagonal.max()));
  upper = std::max(upper, 1e-12);
  int expansions = 0;
  while (equation(upper, NULL) > 0.0 && expansions++ < 100) upper *= 2.0;
  if (!std::isfinite(upper) || equation(upper, NULL) > 0.0) return false;
  double r = std::min(upper, std::max(0.0, initial_norm));
  for (int it = 0; it < 100; ++it) {
    double derivative;
    const double value = equation(r, &derivative);
    if (std::abs(value) <= 2e-14) break;
    if (value > 0.0) lower = r; else upper = r;
    if (upper - lower <= 2e-14 * std::max(1e-12, upper)) {
      r = 0.5 * (lower + upper);
      break;
    }
    const double next = r - value / derivative;
    r = std::isfinite(next) && next > lower && next < upper ? next :
      0.5 * (lower + upper);
  }
  beta = cache.basis * (r * projected / (r * diagonal + lambda1)) +
    (r / (r * lambda2 + lambda1)) * complement;
  return beta.is_finite();
}

JointFit run_irls_stage(
    const arma::mat& X, const arma::vec& y,
    const arma::uvec& first, const arma::uvec& last,
    const arma::vec& weight, const arma::vec& target,
    const double lambda, const double alpha, const double d,
    arma::vec beta, double intercept, const double kkt_tolerance,
    const double intercept_tolerance, const int max_intercept_iterations,
    const bool use_active_set, const bool keep_trace,
    const int max_outer, const int max_inner,
    int& outer_iterations, int& inner_sweeps) {
  JointFit fit{};
  fit.termination = "irls_outer_budget";
  const double n = static_cast<double>(X.n_rows);
  const arma::uword G = first.n_elem;
  arma::vec offset = X * beta;
  InterceptResult prof = profile_intercept(offset, y, intercept,
    intercept_tolerance, max_intercept_iterations);
  fit.intercept_iterations += prof.iterations;
  if (prof.converged) intercept = prof.intercept;
  double objective = logistic_loss(offset + intercept, y) +
    penalty(beta, first, last, weight, target, lambda, alpha, d);
  if (keep_trace) fit.objective_trace.push_back(objective);
  KKTResult kkt = kkt_residual(X, y, beta, intercept, first, last,
    weight, target, lambda, alpha, d);
  ++fit.kkt_scans;
  outer_iterations = 0;
  inner_sweeps = 0;
  if (!prof.converged) fit.termination = "irls_initial_intercept_failure";
  for (int outer = 0; prof.converged && outer < max_outer &&
       kkt.maximum > kkt_tolerance; ++outer) {
    Rcpp::checkUserInterrupt();
    ++outer_iterations;
    const arma::vec probability = stable_expit(offset + intercept);
    // Positive curvature floor ONLY for the local quadratic approximation;
    // never clip probabilities in the likelihood or add a model penalty.
    const arma::vec w = arma::clamp(probability % (1.0 - probability),
      1e-12, 0.25) / n;
    const double sumw = arma::accu(w);
    const arma::rowvec means = (w.t() * X) / sumw;
    arma::mat centered = X;
    centered.each_row() -= means;
    const arma::vec score0 = (y - probability) / n;
    const arma::vec residual0 = score0 - w * (arma::accu(score0) / sumw);
    arma::vec residual = residual0;
    arma::vec candidate = beta;
    const double inner_tolerance = std::max(0.1 * kkt_tolerance,
      std::min(1e-3, 0.05 * kkt.maximum));
    std::vector<QuadraticGroup> cache(G);
    arma::uvec active(G, arma::fill::zeros);
    for (arma::uword g = 0; g < G; ++g) {
      active[g] = !use_active_set || kkt.group[g] > inner_tolerance ||
        arma::norm(beta.subvec(first[g], last[g]), 2) > 1e-12;
    }
    bool inner_ok = false;
    bool failed = false;
    if (alpha == 0.0) {
      // Exact weighted-ridge subproblem using an n-by-n Woodbury solve.
      // This also retains shifted targets in the design's null space.
      arma::vec inverse_ridge(X.n_cols);
      for (arma::uword g = 0; g < G; ++g)
        inverse_ridge.subvec(first[g], last[g]).fill(1.0 / (lambda * weight[g]));
      arma::mat A = centered;
      A.each_col() %= arma::sqrt(w);
      arma::mat scaled = A;
      scaled.each_row() %= inverse_ridge.t();
      const arma::vec rhs = residual0 / arma::sqrt(w) + A * (beta - d * target);
      arma::vec dual;
      inner_ok = arma::solve(dual, arma::eye<arma::mat>(X.n_rows, X.n_rows) +
        scaled * A.t(), rhs, arma::solve_opts::likely_sympd);
      if (inner_ok) candidate = d * target + scaled.t() * dual;
      ++inner_sweeps;
      fit.group_updates += static_cast<double>(G);
    } else {
      for (int sweep = 0; sweep < max_inner; ++sweep) {
        if (sweep % 10 == 0) Rcpp::checkUserInterrupt();
        ++inner_sweeps;
        for (arma::uword g = 0; g < G; ++g) {
          if (!active[g]) continue;
          const auto Xg = centered.cols(first[g], last[g]);
          if (!cache[g].ready) {
            arma::mat A = Xg;
            A.each_col() %= arma::sqrt(w);
            if (!prepare_quadratic_group(A, cache[g])) { failed = true; break; }
          }
          const arma::vec bg = candidate.subvec(first[g], last[g]);
          const double l1 = lambda * alpha * weight[g];
          const double l2 = lambda * (1.0 - alpha) * weight[g];
          const arma::vec score = Xg.t() * (residual + w % (Xg * bg)) +
            l2 * d * target.subvec(first[g], last[g]);
          arma::vec next;
          if (!solve_quadratic_group(cache[g], score, l1, l2,
              arma::norm(bg, 2), next)) { failed = true; break; }
          residual -= w % (Xg * (next - bg));
          candidate.subvec(first[g], last[g]) = next;
          ++fit.group_updates;
        }
        if (failed) break;
        if (sweep % 5 != 0 && sweep + 1 < max_inner) continue;
        // Reconstruct residual and scan ALL groups, including inactive ones.
        // No strong-rule exclusion is ever accepted without this check.
        residual = residual0 - w % (centered * (candidate - beta));
        const arma::vec gradient = -centered.t() * residual;
        double maximum = 0.0;
        for (arma::uword g = 0; g < G; ++g) {
          const arma::vec bg = candidate.subvec(first[g], last[g]);
          const double norm = arma::norm(bg, 2);
          const double l1 = lambda * alpha * weight[g];
          const double l2 = lambda * (1.0 - alpha) * weight[g];
          const arma::vec gg = gradient.subvec(first[g], last[g]) +
            l2 * (bg - d * target.subvec(first[g], last[g]));
          const double error = norm > 1e-10 ? arma::norm(gg + l1 * bg / norm, 2) :
            std::max(0.0, arma::norm(gg, 2) - l1);
          maximum = std::max(maximum, error);
          if (error > inner_tolerance) active[g] = 1u;
        }
        if (maximum <= inner_tolerance) { inner_ok = true; break; }
      }
    }
    if (failed || !inner_ok || !candidate.is_finite()) {
      fit.termination = failed ? "irls_quadratic_failure" : "irls_inner_budget";
      break;
    }
    // Profile the intercept in the TRUE likelihood during line search.
    const arma::vec direction = candidate - beta;
    const arma::vec direction_offset = X * direction;
    const arma::vec gradient = X.t() * (probability - y) / n;
    const double delta = arma::dot(gradient, direction) +
      penalty(candidate, first, last, weight, target, lambda, alpha, d) -
      penalty(beta, first, last, weight, target, lambda, alpha, d);
    const double slack = 32.0 * std::numeric_limits<double>::epsilon() *
      (1.0 + std::abs(objective));
    bool accepted = false;
    double step = 1.0;
    for (int ls = 0; ls < 50; ++ls) {
      const arma::vec trial = beta + step * direction;
      const arma::vec trial_offset = offset + step * direction_offset;
      InterceptResult trial_prof = profile_intercept(trial_offset, y,
        intercept, intercept_tolerance, max_intercept_iterations);
      fit.intercept_iterations += trial_prof.iterations;
      const double trial_objective = logistic_loss(
        trial_offset + trial_prof.intercept, y) +
        penalty(trial, first, last, weight, target, lambda, alpha, d);
      fit.max_raw_objective_increase = std::max(fit.max_raw_objective_increase,
        trial_objective - objective);
      if (trial_prof.converged && std::isfinite(trial_objective) &&
          trial_objective <= objective + 1e-4 * step * std::min(0.0, delta) + slack) {
        fit.max_accepted_objective_increase = std::max(
          fit.max_accepted_objective_increase, trial_objective - objective);
        beta = trial;
        offset = X * beta; // eliminate recursively accumulated prediction error
        intercept = trial_prof.intercept;
        objective = logistic_loss(offset + intercept, y) +
          penalty(beta, first, last, weight, target, lambda, alpha, d);
        accepted = true;
        if (keep_trace) fit.objective_trace.push_back(objective);
        break;
      }
      step *= 0.5;
      ++fit.backtracking_steps;
    }
    if (!accepted) { fit.termination = "irls_line_search_failure"; break; }
    kkt = kkt_residual(X, y, beta, intercept, first, last, weight,
      target, lambda, alpha, d);
    ++fit.kkt_scans;
  }
  // Independent original-objective KKT check, including the intercept.
  prof = profile_intercept(X * beta, y, intercept, intercept_tolerance,
    max_intercept_iterations);
  fit.intercept_iterations += prof.iterations;
  if (prof.converged) intercept = prof.intercept;
  kkt = kkt_residual(X, y, beta, intercept, first, last, weight,
    target, lambda, alpha, d);
  ++fit.kkt_scans;
  fit.beta = beta;
  fit.intercept = intercept;
  fit.objective = logistic_loss(X * beta + intercept, y) +
    penalty(beta, first, last, weight, target, lambda, alpha, d);
  fit.kkt = kkt.maximum;
  fit.group_kkt = kkt.group_maximum;
  fit.intercept_kkt = kkt.intercept;
  fit.converged = prof.converged && kkt.maximum <= kkt_tolerance;
  if (fit.converged) fit.termination = "irls_kkt_tolerance";
  fit.sweeps = inner_sweeps;
  for (arma::uword g = 0; g < G; ++g)
    fit.active_groups += arma::norm(beta.subvec(first[g], last[g]), 2) > 1e-12;
  return fit;
}

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
  arma::vec offset;
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
    const arma::mat& Xt,
    const arma::vec& y,
    const arma::uvec& group_start,
    const arma::uvec& group_end,
    const arma::vec& group_weight,
    const arma::vec& target,
    const double lambda,
    const double alpha,
    const double d,
    const arma::vec& point,
    const arma::vec& point_offset,
    const double point_intercept,
    double curvature,
    const double intercept_tolerance,
    const int max_intercept_iterations,
    const bool require_monotone,
    const double monotone_reference,
    ApgTiming& timing) {
  auto timed_at = std::chrono::steady_clock::now();
  InterceptResult point_profile = profile_intercept(
    point_offset, y, point_intercept, intercept_tolerance,
    max_intercept_iterations
  );
  timing.intercept_seconds += std::chrono::duration<double>(
    std::chrono::steady_clock::now() - timed_at).count();
  int intercept_iterations = point_profile.iterations;
  if (!point_profile.converged) {
    return ProximalTrial{
      point, point_offset, point_intercept,
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(), curvature,
      intercept_iterations, 0, false
    };
  }
  const arma::vec point_eta = point_profile.intercept + point_offset;
  const double point_loss = logistic_loss(point_eta, y);
  const arma::vec point_probability = stable_expit(point_eta);
  timed_at = std::chrono::steady_clock::now();
  const arma::vec gradient =
    Xt * (point_probability - y) / static_cast<double>(X.n_rows);
  timing.gradient_seconds += std::chrono::duration<double>(
    std::chrono::steady_clock::now() - timed_at).count();

  for (int backtrack = 0; backtrack < 80; ++backtrack) {
    timed_at = std::chrono::steady_clock::now();
    const arma::vec candidate = shifted_group_prox(
      point, gradient, group_start, group_end, group_weight, target,
      lambda, alpha, d, curvature
    );
    timing.prox_seconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - timed_at).count();
    const arma::vec delta = candidate - point;
    timed_at = std::chrono::steady_clock::now();
    const arma::vec candidate_offset = X * candidate;
    timing.matvec_seconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - timed_at).count();
    timed_at = std::chrono::steady_clock::now();
    InterceptResult candidate_profile = profile_intercept(
      candidate_offset, y, point_profile.intercept, intercept_tolerance,
      max_intercept_iterations
    );
    timing.intercept_seconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - timed_at).count();
    intercept_iterations += candidate_profile.iterations;
    if (candidate_profile.converged) {
      const double candidate_loss = logistic_loss(
        candidate_profile.intercept + candidate_offset, y
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
          candidate, candidate_offset, candidate_profile.intercept,
          candidate_loss,
          candidate_objective, curvature, intercept_iterations,
          backtrack, true
        };
      }
    }
    curvature *= 2.0;
    if (!std::isfinite(curvature)) break;
  }
  return ProximalTrial{
    point, point_offset, point_profile.intercept, point_loss,
    std::numeric_limits<double>::infinity(), curvature,
    intercept_iterations, 80, false
  };
}

JointFit run_apg_stage(
    const arma::mat& X,
    const arma::mat& Xt,
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
    int& apg_restarts,
    ApgTiming& timing,
    const bool reuse_apg_offsets) {
  arma::vec current = start.beta;
  arma::vec previous = current;
  double current_intercept = start.intercept;
  auto timed_at = std::chrono::steady_clock::now();
  arma::vec current_offset = X * current;
  arma::vec previous_offset;
  if (reuse_apg_offsets) previous_offset = current_offset;
  timing.matvec_seconds += std::chrono::duration<double>(
    std::chrono::steady_clock::now() - timed_at).count();
  timed_at = std::chrono::steady_clock::now();
  InterceptResult initial_profile = profile_intercept(
    current_offset, y, current_intercept, intercept_tolerance,
    max_intercept_iterations
  );
  timing.intercept_seconds += std::chrono::duration<double>(
    std::chrono::steady_clock::now() - timed_at).count();
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
    current_intercept + current_offset, y
  ) + penalty(
    current, group_start, group_end, group_weight, target,
    lambda, alpha, d
  );
  if (keep_trace) result.objective_trace.push_back(current_objective);
  timed_at = std::chrono::steady_clock::now();
  KKTResult kkt = kkt_residual(
    X, y, current, current_intercept, group_start, group_end,
    group_weight, target, lambda, alpha, d
  );
  timing.kkt_seconds += std::chrono::duration<double>(
    std::chrono::steady_clock::now() - timed_at).count();
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
    arma::vec extrapolated_offset = current_offset;
    if (momentum != 0.0) {
      timed_at = std::chrono::steady_clock::now();
      // Each accepted offset is still computed directly as X * beta.
      // Linearity saves the extra extrapolation matvec without accumulating
      // recursive updates of the accepted predictions.
      if (reuse_apg_offsets) {
        extrapolated_offset = current_offset +
          momentum * (current_offset - previous_offset);
      } else {
        extrapolated_offset = X * extrapolated;
      }
      timing.matvec_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - timed_at).count();
    }
    double extrapolated_intercept = current_intercept;
    ProximalTrial trial = proximal_trial(
      X, Xt, y, group_start, group_end, group_weight, target,
      lambda, alpha, d, extrapolated, extrapolated_offset,
      extrapolated_intercept,
      std::max(0.8 * curvature, minimum_curvature),
      intercept_tolerance, max_intercept_iterations, false,
      current_objective, timing
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
        X, Xt, y, group_start, group_end, group_weight, target,
        lambda, alpha, d, current, current_offset, current_intercept,
        std::max(trial.curvature, minimum_curvature),
        intercept_tolerance, max_intercept_iterations, true,
        current_objective, timing
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
    if (reuse_apg_offsets) previous_offset = current_offset;
    current_offset = trial.offset;
    current_intercept = trial.intercept;
    current_objective = trial.objective;
    curvature = trial.curvature;
    acceleration = restarted ?
      0.5 * (1.0 + std::sqrt(5.0)) : next_acceleration;
    // Directional restart does not subtract nearly equal objective values.
    if (high_precision &&
        arma::dot(extrapolated - current, current - old_current) > 0.0) {
      previous = current;
      if (reuse_apg_offsets) previous_offset = current_offset;
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
      timed_at = std::chrono::steady_clock::now();
      kkt = kkt_residual(
        X, y, current, current_intercept, group_start, group_end,
        group_weight, target, lambda, alpha, d
      );
      timing.kkt_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - timed_at).count();
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
      if (reuse_apg_offsets) previous_offset = current_offset;
      acceleration = 1.0;
      ++apg_restarts;
    }
  }

  timed_at = std::chrono::steady_clock::now();
  InterceptResult final_profile = profile_intercept(
    current_offset, y, current_intercept, intercept_tolerance,
    max_intercept_iterations
  );
  timing.intercept_seconds += std::chrono::duration<double>(
    std::chrono::steady_clock::now() - timed_at).count();
  result.intercept_iterations += final_profile.iterations;
  if (final_profile.converged) current_intercept = final_profile.intercept;
  const double final_objective = logistic_loss(
    current_intercept + current_offset, y
  ) + penalty(
    current, group_start, group_end, group_weight, target,
    lambda, alpha, d
  );
  timed_at = std::chrono::steady_clock::now();
  kkt = kkt_residual(
    X, y, current, current_intercept, group_start, group_end,
    group_weight, target, lambda, alpha, d
  );
  timing.kkt_seconds += std::chrono::duration<double>(
    std::chrono::steady_clock::now() - timed_at).count();
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
    const arma::mat& Xt,
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
    const bool keep_trace,
    const bool reuse_apg_offsets,
    const bool use_irls,
    const int irls_max_outer,
    const int irls_max_inner) {
  if (use_irls) {
    const auto irls_started = std::chrono::steady_clock::now();
    int outer = 0, inner = 0;
    JointFit fit = run_irls_stage(X, y, group_start, group_end, group_weight,
      target, lambda, alpha, d, beta_initial, intercept_initial, kkt_tolerance,
      intercept_tolerance, max_intercept_iterations, use_active_set, keep_trace,
      irls_max_outer, irls_max_inner, outer, inner);
    const double seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - irls_started).count();
    HybridFit output{};
    if (!fit.converged && enable_fallback) {
      output = fit_one_hybrid(X, Xt, y, group_start, group_end, group_weight,
        target, conservative_curvature, global_curvature, lambda, alpha, d,
        fit.beta, fit.intercept, block_max_sweeps, block_chunk_sweeps,
        block_stall_window, block_stall_relative_improvement, apg_max_iterations,
        apg_kkt_check_interval, kkt_tolerance, update_tolerance, intercept_tolerance,
        max_intercept_iterations, use_active_set, enable_fallback, keep_trace,
        reuse_apg_offsets, false, irls_max_outer, irls_max_inner);
      output.solver_route = "irls_then_" + output.solver_route;
      output.fallback_used = true;
      output.fit.intercept_iterations += fit.intercept_iterations;
      output.fit.kkt_scans += fit.kkt_scans;
      output.fit.group_updates += fit.group_updates;
      output.fit.max_raw_objective_increase = std::max(
        output.fit.max_raw_objective_increase, fit.max_raw_objective_increase);
      output.fit.max_accepted_objective_increase = std::max(
        output.fit.max_accepted_objective_increase, fit.max_accepted_objective_increase);
      if (keep_trace) {
        append_trace(fit.objective_trace, output.fit.objective_trace);
        output.fit.objective_trace = fit.objective_trace;
      }
    } else {
      output.fit = fit;
      output.solver_route = "irls_proximal_newton";
      output.global_logistic_curvature = global_curvature;
      output.block_termination = "not_used";
    }
    output.diagnostic_irls_seconds = seconds;
    output.irls_iterations = outer;
    output.irls_inner_sweeps = inner;
    output.irls_termination = fit.termination;
    return output;
  }
  // In high-precision mode the block stage is a warm-up, not the final fit.
  // Avoid spending thousands of sweeps near its round-off-sensitive floor;
  // the APG stage still has to satisfy the caller's original tight tolerance.
  const double block_tolerance = (enable_fallback && kkt_tolerance < 1e-8) ?
    std::max(kkt_tolerance, 2e-6) : kkt_tolerance;
  const auto block_started = std::chrono::steady_clock::now();
  BlockAudit block = run_block_stage(
    X, y, group_start, group_end, group_weight, target,
    conservative_curvature, lambda, alpha, d, beta_initial,
    intercept_initial, block_max_sweeps, block_chunk_sweeps,
    block_stall_window, block_stall_relative_improvement,
    block_tolerance, update_tolerance, intercept_tolerance,
    max_intercept_iterations, use_active_set, keep_trace
  );
  HybridFit output;
  output.diagnostic_block_seconds = std::chrono::duration<double>(
    std::chrono::steady_clock::now() - block_started).count();
  output.diagnostic_apg_seconds = 0.0;
  output.diagnostic_apg_parts = ApgTiming{0.0, 0.0, 0.0, 0.0, 0.0};
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
  const auto apg_started = std::chrono::steady_clock::now();
  JointFit fallback = run_apg_stage(
    X, Xt, y, group_start, group_end, group_weight, target,
    lambda, alpha, d, block.fit, apg_max_iterations,
    kkt_tolerance, update_tolerance, intercept_tolerance,
    max_intercept_iterations, apg_kkt_check_interval, global_curvature,
    keep_trace, apg_backtracking, apg_restarts,
    output.diagnostic_apg_parts, reuse_apg_offsets
  );
  output.diagnostic_apg_seconds = std::chrono::duration<double>(
    std::chrono::steady_clock::now() - apg_started).count();
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
    Rcpp::_ ["diagnostic_irls_seconds"] = value.diagnostic_irls_seconds,
    Rcpp::_ ["irls_iterations"] = value.irls_iterations,
    Rcpp::_ ["irls_inner_sweeps"] = value.irls_inner_sweeps,
    Rcpp::_ ["irls_termination"] = value.irls_termination,
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
    Rcpp::_ ["diagnostic_apg_matvec_seconds"] =
      value.diagnostic_apg_parts.matvec_seconds,
    Rcpp::_ ["diagnostic_apg_intercept_seconds"] =
      value.diagnostic_apg_parts.intercept_seconds,
    Rcpp::_ ["diagnostic_apg_gradient_seconds"] =
      value.diagnostic_apg_parts.gradient_seconds,
    Rcpp::_ ["diagnostic_apg_prox_seconds"] =
      value.diagnostic_apg_parts.prox_seconds,
    Rcpp::_ ["diagnostic_apg_kkt_seconds"] =
      value.diagnostic_apg_parts.kkt_seconds,
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
    const bool keep_trace = false,
    const bool reuse_apg_offsets = false,
    const bool use_irls = false,
    const int irls_max_outer = 50,
    const int irls_max_inner = 2000) {
  validate_groups(X, group_start, group_end, group_weight, target);
  validate_response(y, X.n_rows);
  if (irls_max_outer < 1 || irls_max_inner < 1)
    Rcpp::stop("Local IRLS budgets must be positive.");
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
  const arma::mat Xt = X.t();
  return wrap_hybrid_fit(fit_one_hybrid(
    X, Xt, y, group_start, group_end, group_weight, target,
    group_curvature, global_curvature, lambda, alpha, d,
    beta_initial, intercept_initial, block_max_sweeps,
    block_chunk_sweeps, block_stall_window,
    block_stall_relative_improvement, apg_max_iterations,
    apg_kkt_check_interval, kkt_tolerance, update_tolerance,
    intercept_tolerance, max_intercept_iterations, use_active_set,
    enable_fallback, keep_trace, reuse_apg_offsets, use_irls,
    irls_max_outer, irls_max_inner
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
    const bool keep_traces = false,
    const double early_exit_deviance_fraction = -1.0,
    const bool reuse_apg_offsets = false,
    const bool use_irls = false,
    const int irls_max_outer = 50,
    const int irls_max_inner = 2000) {
  validate_groups(X, group_start, group_end, group_weight, target);
  validate_response(y, X.n_rows);
  if (irls_max_outer < 1 || irls_max_inner < 1)
    Rcpp::stop("Local IRLS budgets must be positive.");
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
  const bool early_exit = early_exit_deviance_fraction > 0.0;
  if (!std::isfinite(early_exit_deviance_fraction) ||
      (early_exit_deviance_fraction != -1.0 &&
       !(early_exit_deviance_fraction > 0.0 &&
         early_exit_deviance_fraction < 1.0))) {
    Rcpp::stop("Local early-exit fraction must be -1 or strictly between 0 and 1.");
  }

  const arma::uword p = X.n_cols;
  const arma::uword L = lambda.n_elem;
  const arma::uword D = d.n_elem;
  const double prevalence = arma::mean(y);
  const double null_intercept = std::log(prevalence) -
    std::log1p(-prevalence);
  const double null_training_loss = early_exit ? logistic_loss(
    arma::vec(y.n_elem, arma::fill::value(null_intercept)), y) : 0.0;
  if (early_exit && !(null_training_loss > 0.0)) {
    Rcpp::stop("Local early exit requires positive null training loss.");
  }
  const arma::vec group_curvature = conservative_group_curvature(
    X, group_start, group_end
  );
  const double global_curvature = global_logistic_curvature(X);
  const arma::mat Xt = X.t();

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
  arma::mat diagnostic_fit_seconds_path(L, D, arma::fill::zeros);
  arma::mat diagnostic_irls_seconds_path(L, D, arma::fill::zeros);
  arma::imat irls_iterations_path(L, D, arma::fill::zeros);
  arma::imat irls_inner_sweeps_path(L, D, arma::fill::zeros);
  arma::mat diagnostic_block_seconds_path(L, D, arma::fill::zeros);
  arma::mat diagnostic_apg_seconds_path(L, D, arma::fill::zeros);
  arma::mat diagnostic_apg_matvec_seconds_path(L, D, arma::fill::zeros);
  arma::mat diagnostic_apg_intercept_seconds_path(L, D, arma::fill::zeros);
  arma::mat diagnostic_apg_gradient_seconds_path(L, D, arma::fill::zeros);
  arma::mat diagnostic_apg_prox_seconds_path(L, D, arma::fill::zeros);
  arma::mat diagnostic_apg_kkt_seconds_path(L, D, arma::fill::zeros);
  arma::imat active_groups_path(L, D, arma::fill::zeros);
  Rcpp::CharacterMatrix termination_path(L, D);
  Rcpp::CharacterMatrix block_termination_path(L, D);
  Rcpp::CharacterMatrix route_path(L, D);
  Rcpp::CharacterMatrix irls_termination_path(L, D);
  Rcpp::List traces(L * D);
  arma::mat training_deviance_path(L, D, arma::fill::zeros);
  Rcpp::IntegerVector computed_count(D);
  Rcpp::LogicalVector threshold_reached(D);

  for (arma::uword di = 0; di < D; ++di) {
    arma::vec beta(p, arma::fill::zeros);
    double intercept = null_intercept;
    if (warm_start_d && di > 0u) {
      beta = beta_path.slice(di - 1u).col(0u);
      intercept = intercept_path(0u, di - 1u);
    }
    for (arma::uword li = 0; li < L; ++li) {
      const auto fit_started = std::chrono::steady_clock::now();
      HybridFit value = fit_one_hybrid(
        X, Xt, y, group_start, group_end, group_weight, target,
        group_curvature, global_curvature, lambda[li], alpha, d[di],
        beta, intercept, block_max_sweeps, block_chunk_sweeps,
        block_stall_window, block_stall_relative_improvement,
        apg_max_iterations, apg_kkt_check_interval, kkt_tolerance,
        update_tolerance, intercept_tolerance, max_intercept_iterations,
        use_active_set, enable_fallback, keep_traces, reuse_apg_offsets,
        use_irls, irls_max_outer, irls_max_inner
      );
      diagnostic_fit_seconds_path(li, di) = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - fit_started).count();
      diagnostic_irls_seconds_path(li, di) = value.diagnostic_irls_seconds;
      irls_iterations_path(li, di) = value.irls_iterations;
      irls_inner_sweeps_path(li, di) = value.irls_inner_sweeps;
      irls_termination_path(li, di) = value.irls_termination;
      diagnostic_block_seconds_path(li, di) =
        value.diagnostic_block_seconds;
      diagnostic_apg_seconds_path(li, di) =
        value.diagnostic_apg_seconds;
      diagnostic_apg_matvec_seconds_path(li, di) =
        value.diagnostic_apg_parts.matvec_seconds;
      diagnostic_apg_intercept_seconds_path(li, di) =
        value.diagnostic_apg_parts.intercept_seconds;
      diagnostic_apg_gradient_seconds_path(li, di) =
        value.diagnostic_apg_parts.gradient_seconds;
      diagnostic_apg_prox_seconds_path(li, di) =
        value.diagnostic_apg_parts.prox_seconds;
      diagnostic_apg_kkt_seconds_path(li, di) =
        value.diagnostic_apg_parts.kkt_seconds;
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
      if (early_exit) {
        const double training_loss = logistic_loss(
          X * fit.beta + fit.intercept, y);
        training_deviance_path(li, di) = 1.0 -
          training_loss / null_training_loss;
        computed_count[di] = static_cast<int>(li + 1u);
        if (fit.converged && std::isfinite(fit.kkt) &&
            fit.kkt <= kkt_tolerance &&
            training_deviance_path(li, di) >=
              early_exit_deviance_fraction) {
          threshold_reached[di] = true;
          break;
        }
      }
    }
  }

  if (early_exit) {
    // A distinct diagnostic schema returns only fitted prefixes. The
    // production wrapper must never mistake omitted lambdas for fitted zeros.
    Rcpp::List paths(D);
    for (arma::uword di = 0; di < D; ++di) {
      const arma::uword n = static_cast<arma::uword>(computed_count[di]);
      if (n == 0u) Rcpp::stop("Local early-exit path is empty.");
      paths[di] = Rcpp::List::create(
        Rcpp::_ ["d"] = d[di],
        Rcpp::_ ["lambda"] = lambda.subvec(0u, n - 1u),
        Rcpp::_ ["beta"] = beta_path.slice(di).cols(0u, n - 1u),
        Rcpp::_ ["intercept"] = intercept_path.col(di).subvec(0u, n - 1u),
        Rcpp::_ ["objective"] = objective_path.col(di).subvec(0u, n - 1u),
        Rcpp::_ ["kkt"] = kkt_path.col(di).subvec(0u, n - 1u),
        Rcpp::_ ["converged"] = converged_path.col(di).subvec(0u, n - 1u),
        Rcpp::_ ["training_deviance_explained"] =
          training_deviance_path.col(di).subvec(0u, n - 1u),
        Rcpp::_ ["diagnostic_fit_seconds"] =
          diagnostic_fit_seconds_path.col(di).subvec(0u, n - 1u)
      );
    }
    return Rcpp::List::create(
      Rcpp::_ ["schema_version"] = "logistic_sglasso_local_early_exit_v1",
      Rcpp::_ ["paths"] = paths,
      Rcpp::_ ["computed_lambda_count"] = computed_count,
      Rcpp::_ ["threshold_reached"] = threshold_reached,
      Rcpp::_ ["requested_lambda_count"] = static_cast<int>(L),
      Rcpp::_ ["threshold"] = early_exit_deviance_fraction,
      Rcpp::_ ["null_training_log_loss"] = null_training_loss
    );
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
    Rcpp::_ ["diagnostic_fit_seconds"] = diagnostic_fit_seconds_path,
    Rcpp::_ ["diagnostic_irls_seconds"] = diagnostic_irls_seconds_path,
    Rcpp::_ ["irls_iterations"] = irls_iterations_path,
    Rcpp::_ ["irls_inner_sweeps"] = irls_inner_sweeps_path,
    Rcpp::_ ["irls_termination"] = irls_termination_path,
    Rcpp::_ ["diagnostic_block_seconds"] = diagnostic_block_seconds_path,
    Rcpp::_ ["diagnostic_apg_seconds"] = diagnostic_apg_seconds_path,
    Rcpp::_ ["diagnostic_apg_matvec_seconds"] =
      diagnostic_apg_matvec_seconds_path,
    Rcpp::_ ["diagnostic_apg_intercept_seconds"] =
      diagnostic_apg_intercept_seconds_path,
    Rcpp::_ ["diagnostic_apg_gradient_seconds"] =
      diagnostic_apg_gradient_seconds_path,
    Rcpp::_ ["diagnostic_apg_prox_seconds"] =
      diagnostic_apg_prox_seconds_path,
    Rcpp::_ ["diagnostic_apg_kkt_seconds"] =
      diagnostic_apg_kkt_seconds_path,
    Rcpp::_ ["objective_traces"] = traces,
    Rcpp::_ ["solver"] = use_irls ? "local_irls_proximal_newton" :
      "profiled_block_monotone_apg_v11",
    Rcpp::_ ["warm_start_d"] = warm_start_d,
    Rcpp::_ ["use_active_set"] = use_active_set,
    Rcpp::_ ["fallback_enabled"] = enable_fallback,
    Rcpp::_ ["kkt_tolerance"] = kkt_tolerance,
    Rcpp::_ ["intercept_tolerance"] = intercept_tolerance
  );
}
