#include <RcppArmadillo.h>
#include <chrono>
#include <cmath>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;




// Check for convergence of any beta vector
bool checkConvergence(const arma::colvec& beta_new, 
                      const arma::colvec& beta_old, 
                      const double& eps,
                      const arma::uword& p) {
  for (arma::uword j=0; j < p; j++) {
    if (beta_new[j]!=0 && beta_old[j]!=0) {
      if (abs(beta_new[j]-beta_old[j]) > eps) {
        return(false);
      }
    } else if (beta_new[j]==0 && beta_old[j]!=0) {
      return(false);
    } else if (beta_new[j]!=0 && beta_old[j]==0) {
      return(false);
    }
  }
  return(true);
}

bool checkConvergenceActive(const arma::colvec& beta_new,
                            const arma::colvec& beta_old,
                            const double& eps,
                            const arma::colvec& K1,
                            const std::vector<arma::uword>& active_groups) {
  for (arma::uword active_indx = 0; active_indx < active_groups.size(); active_indx++) {
    arma::uword g = active_groups[active_indx];
    for (arma::uword j = K1[g]; j < K1[g + 1]; j++) {
      if (beta_new[j] != 0 && beta_old[j] != 0) {
        if (abs(beta_new[j] - beta_old[j]) > eps) {
          return(false);
        }
      } else if (beta_new[j] == 0 && beta_old[j] != 0) {
        return(false);
      } else if (beta_new[j] != 0 && beta_old[j] == 0) {
        return(false);
      }
    }
  }
  return(true);
}

inline void activate_group(arma::uword g,
                           arma::colvec& ever_active,
                           std::vector<arma::uword>& active_groups) {
  if (ever_active[g] == 0) {
    ever_active[g] = 1;
    active_groups.push_back(g);
  }
}

inline double dot_col(const arma::mat& X, const arma::colvec& y, arma::uword j) {
  const double* xj = X.colptr(j);
  const double* yy = y.memptr();
  double out = 0.0;
  for (arma::uword i = 0; i < X.n_rows; i++) {
    out += xj[i] * yy[i];
  }
  return out;
}

inline void add_scaled_col(const arma::mat& X,
                           arma::uword j,
                           double scale,
                           arma::colvec& y) {
  const double* xj = X.colptr(j);
  double* yy = y.memptr();
  for (arma::uword i = 0; i < X.n_rows; i++) {
    yy[i] += scale * xj[i];
  }
}


// Gaussian loss function
double gLoss(const arma::colvec& r, int n) {
  return arma::dot(r, r);
}

// standardize
// [[Rcpp::export]]
List std_c(const arma::mat & x){
  // standardize a matrix
  List output;
  arma::uword n = x.n_rows;
  arma::uword p = x.n_cols;
  double er = 0;
  arma::mat xx  = zeros(n, p);
  arma::vec c = zeros(p, 1);
  arma::vec s = zeros(p, 1);
  for (arma::uword j = 0; j < p; j++){
    //Center
    c[j] = sum(x.col(j))/n;
    xx.col(j) = x.col(j) - c[j];
    er = dot(xx.col(j),xx.col(j))/n;
    s[j] = sqrt(er);
    if (s[j] > 1e-6) {
      xx.col(j) /= s[j];
    } else {
      xx.col(j).zeros();
    }
  }
  output["xx"] = xx;
  output["c"] = c;
  output["s"] = s;
  return(output);
}

// Orthogonalize each group using the same SVD transform as the R helper, but
// keep the loop in C++ to avoid thousands of small R-level SVD calls.
// [[Rcpp::export]]
List orthogonalize_c(const arma::mat& X, const arma::uvec& group) {
  arma::uword n = X.n_rows;
  arma::uword p = X.n_cols;
  arma::uword J = group.max();
  arma::mat XX(n, p, arma::fill::zeros);
  Rcpp::List T(J);

  std::vector<std::vector<arma::uword> > group_index(J + 1);
  for (arma::uword j = 0; j < p; j++) {
    arma::uword g = group[j];
    if (g <= J) group_index[g].push_back(j);
  }

  for (arma::uword idx = 0; idx < group_index[0].size(); idx++) {
    arma::uword col = group_index[0][idx];
    XX.col(col) = X.col(col);
  }

  const double sqrt_n = std::sqrt(static_cast<double>(n));
  for (arma::uword g = 1; g <= J; g++) {
    std::vector<arma::uword>& idx_vec = group_index[g];
    if (idx_vec.empty()) continue;

    arma::uvec ind(idx_vec.size());
    for (arma::uword k = 0; k < idx_vec.size(); k++) ind[k] = idx_vec[k];

    arma::mat Xg = X.cols(ind);
    arma::mat U;
    arma::vec s;
    arma::mat V;
    bool ok = arma::svd_econ(U, s, V, Xg, "both", "std");
    if (!ok) {
      arma::svd(U, s, V, Xg, "std");
    }

    arma::uvec r = arma::find(s > 1e-10);
    if (r.n_elem == 0) continue;

    arma::mat Tg = V.cols(r);
    arma::rowvec scale = arma::trans(sqrt_n / s.elem(r));
    Tg.each_row() %= scale;
    T[g - 1] = Tg;
    XX.cols(ind.elem(r)) = Xg * Tg;
  }

  arma::rowvec col_abs = arma::sum(arma::abs(XX), 0);
  arma::uvec nz = arma::find(arma::trans(col_abs) > 0);
  arma::mat XX_nz = XX.cols(nz);
  arma::uvec group_nz = group.elem(nz);
  return Rcpp::List::create(
    Rcpp::_["X"] = XX_nz,
    Rcpp::_["T"] = T,
    Rcpp::_["group"] = group_nz
  );
}


// Soft Threshold function
double S_c(double a, double lambda){
  if (a > lambda){
    return(a-lambda);
  }else if (a < -lambda){
    return(a + lambda);
  }else{
    return(0);
  }
}



// Loss function of sglasso regression
double loss_c(const arma::mat& Xtilde, 
               const arma::colvec& Ytilde, 
               const arma::colvec& beta,
               const arma::colvec& beta_lse,
               const arma::colvec& gm,
               double lambda,
               double d,
               double n,
               double J,
               double alpha,
               const arma::colvec& K1){
  // Calculate and return objective value
  double pen=0;
  for(int g = 0; g < J; g++){
    for(int j = K1[g]; j < K1[g+1]; j++){
   pen =  pen + lambda*gm[g]*alpha*abs(beta[j]) + .5*lambda*gm[g]*(1-alpha)*pow((beta[j]-d*beta_lse[j]),2);
    }
  }
  return accu(square(Ytilde - Xtilde * beta)/(2 * n)) + pen;
}


// lambda max value
// [[Rcpp::export]]
double lambda_max_c(const arma::mat& Xtilde, 
                       const arma::colvec& Ytilde, 
                       const arma::colvec& K,
                       const arma::colvec& K1){
  
  // Number of data points
  int J = K1.n_rows-1;
  // calculate z
  double zmax = 0;
  for(int g = 0; g < J; g++){
      for(int j = K1[g]; j < K1[g+1]; j++){
        double z = std::abs(dot_col(Xtilde, Ytilde, j) / sqrt(K[g]));
        if(z > zmax) zmax = z;
        }
  }

  // Return value
  return zmax;
  
}  



// sequential strong rule
void ssr_sglasso(const arma::colvec& K,
                const arma::colvec& K1,
                const arma::mat& Xtilde,
                const arma::colvec& r,
                const arma::colvec& gm,
                const arma::colvec& lambda,
                const double lam_max, 
                const arma::uword l_indx,
                const double alpha,
                const double d_val,
                const double n,
                const double J,
                const arma::colvec& xty,
                const arma::colvec& beta_new,
                const arma::colvec& ever_active,
                const int screen_rule,
                arma::colvec& xTr,
                arma::colvec& ever_strong) {
  if (screen_rule == 0) {
    ever_strong.ones();
    for (arma::uword g = 0; g < J; g++) xTr[g] = 0;
    return;
  }

  for(int g = 0; g < J; g++){
    if (ever_active[g] != 0) {
      ever_strong[g] = 1;
      continue;
    }
    double z_norm_sq = 0;
    double lambda_ref = (l_indx != 0) ? lambda[l_indx - 1] : lam_max;
    if (lambda_ref <= 0) lambda_ref = lambda[l_indx];
    double lambda2_ref = lambda_ref * gm[g] * (1 - alpha);
    for(arma::uword j = K1[g]; j < K1[g+1]; j++){
      double z = dot_col(Xtilde, r, j) / n;
      double LZ = z + lambda2_ref * (d_val * xty[j] - beta_new[j]);
      z_norm_sq += LZ * LZ;
    }
    xTr[g] = std::sqrt(z_norm_sq);
    double cutoff;
    double TOLERANCE = 1e-8;
    if (l_indx != 0) {
      cutoff = alpha * gm[g] * (2 * lambda[l_indx] - lambda[l_indx-1]);
    } else {
      if (lam_max > 0) cutoff = alpha * gm[g] * (2 * lambda[l_indx] - lam_max);
      else cutoff = 0;
    }
    if (xTr[g] +  TOLERANCE > cutoff) {
      ever_strong[g] = 1; // not reject, in strong set
    } else {
      ever_strong[g] = 0; //reject
    }
  }
}




// Scan for violations in rest set
bool Check_Rest_Set(const arma::mat& Xtilde, 
                    const arma::colvec& Ytilde,
                    const arma::colvec& r, 
                    const arma::colvec& K1,
                    const arma::colvec& K,
                    const double& lambda_val,
                    const double& alpha,
                    const double& d_val,
                    const double& n, 
                    const double& J, 
                    const arma::colvec& gm,
                    const arma::colvec& xty,
                    const arma::colvec& beta_new,
                    arma::colvec& xTr,
                    arma::colvec& ever_strong, 
                    arma::colvec& ever_active,
                    std::vector<arma::uword>& active_groups,
                    int& kkt_checked,
                    int& kkt_violations,
                    const arma::colvec& bedpp_safe){
  // Write some explanitions
  bool violations = false;
  double TOLERANCE = 1e-8;
  for(arma::uword g = 0; g < J; g++){
    double lambda2 = lambda_val * gm[g] * (1-alpha);
    if(ever_strong[g] == 0 && bedpp_safe[g] == 0){
      kkt_checked += 1;
      double LZ_norm_sq = 0;
      for(arma::uword j = K1[g]; j < K1[g+1]; j++){
        double z = dot_col(Xtilde, r, j) / n;
        double LZ = z + lambda2 * (d_val * xty[j] - beta_new[j]);
        LZ_norm_sq += LZ * LZ;
      }
      xTr[g] = std::sqrt(LZ_norm_sq);
      if (xTr[g]  +  TOLERANCE > alpha * lambda_val * gm[g]) {
        ever_strong[g] = 1;
        activate_group(g, ever_active, active_groups);
        kkt_violations += 1;
        violations = true;
        }
      }
  }
  return(violations);
}

double sglasso_primal_objective(const arma::mat& Xtilde,
                                const arma::colvec& Ytilde,
                                const arma::colvec& r,
                                const arma::colvec& beta,
                                const arma::colvec& K1,
                                const double& lambda_val,
                                const double& alpha,
                                const double& d_val,
                                const double& n,
                                const double& J,
                                const arma::colvec& gm,
                                const arma::colvec& xty) {
  double out = arma::dot(r, r) / (2.0 * n);
  for (arma::uword g = 0; g < J; g++) {
    double beta_norm_sq = 0;
    double shift_norm_sq = 0;
    for (arma::uword j = K1[g]; j < K1[g + 1]; j++) {
      beta_norm_sq += beta[j] * beta[j];
      double shift = beta[j] - d_val * xty[j];
      shift_norm_sq += shift * shift;
    }
    out += lambda_val * alpha * gm[g] * std::sqrt(beta_norm_sq);
    out += 0.5 * lambda_val * (1 - alpha) * gm[g] * shift_norm_sq;
  }
  return out;
}

double sglasso_dual_scale(const arma::mat& Xtilde,
                          const arma::colvec& r,
                          const arma::colvec& beta,
                          const arma::colvec& K1,
                          const double& lambda_val,
                          const double& alpha,
                          const double& d_val,
                          const double& n,
                          const double& J,
                          const arma::colvec& gm,
                          const arma::colvec& xty) {
  double rho = 0;
  if (lambda_val <= 0 || alpha <= 0) return arma::datum::inf;
  for (arma::uword g = 0; g < J; g++) {
    double lambda2 = lambda_val * gm[g] * (1 - alpha);
    double score_norm_sq = 0;
    for (arma::uword j = K1[g]; j < K1[g + 1]; j++) {
      double z = dot_col(Xtilde, r, j) / n;
      double score = z + lambda2 * (d_val * xty[j] - beta[j]);
      score_norm_sq += score * score;
    }
    double tau = lambda_val * alpha * gm[g];
    if (tau > 0) {
      double ratio = std::sqrt(score_norm_sq) / tau;
      if (ratio > rho) rho = ratio;
    }
  }
  return std::max(1.0, rho);
}

double sglasso_dual_objective(const arma::mat& Xtilde,
                              const arma::colvec& Ytilde,
                              const arma::colvec& r,
                              const arma::colvec& beta,
                              const arma::colvec& K1,
                              const double& lambda_val,
                              const double& alpha,
                              const double& d_val,
                              const double& n,
                              const double& J,
                              const arma::colvec& gm,
                              const arma::colvec& xty,
                              const double& dual_scale) {
  double btilde_norm_sq = arma::dot(Ytilde, Ytilde) / n;
  double diff_norm_sq = 0;
  for (arma::uword i = 0; i < Ytilde.n_elem; i++) {
    double diff = Ytilde[i] / std::sqrt(n) - r[i] / (std::sqrt(n) * dual_scale);
    diff_norm_sq += diff * diff;
  }
  for (arma::uword g = 0; g < J; g++) {
    double lambda2 = lambda_val * gm[g] * (1 - alpha);
    for (arma::uword j = K1[g]; j < K1[g + 1]; j++) {
      double center = d_val * xty[j];
      btilde_norm_sq += lambda2 * center * center;
      double diff = std::sqrt(lambda2) * (center - (center - beta[j]) / dual_scale);
      diff_norm_sq += diff * diff;
    }
  }
  return 0.5 * btilde_norm_sq - 0.5 * diff_norm_sq;
}

double exact_group_bedpp_bound(const arma::mat& Xtilde,
                               const arma::colvec& K1,
                               const arma::uword& g,
                               const double& lambda2,
                               const double& n) {
  arma::uword start = static_cast<arma::uword>(K1[g]);
  arma::uword end = static_cast<arma::uword>(K1[g + 1] - 1);
  arma::mat block = Xtilde.cols(start, end).t() * Xtilde / n;
  for (arma::uword j = start; j <= end; j++) {
    block(j - start, j) += lambda2;
  }
  return arma::norm(block, 2);
}

void bedpp_sglasso_safe(const arma::mat& Xtilde,
                        const arma::colvec& Ytilde,
                        const arma::colvec& r,
                        const arma::colvec& beta,
                        const arma::colvec& K,
                        const arma::colvec& K1,
                        const double& lambda_val,
                        const double& alpha,
                        const double& d_val,
                        const double& n,
                        const double& J,
                        const arma::colvec& gm,
                        const arma::colvec& xty,
                        const arma::colvec& ever_active,
                        const arma::colvec& ever_strong,
                        const double& x_fro_norm,
                        const bool exact_bound,
                        double& bedpp_gap,
                        double& bedpp_radius,
                        double& bedpp_mu,
                        double& bedpp_min_margin,
                        arma::colvec& bedpp_safe) {
  bedpp_safe.zeros();
  bedpp_gap = arma::datum::nan;
  bedpp_radius = arma::datum::nan;
  bedpp_mu = arma::datum::nan;
  bedpp_min_margin = arma::datum::nan;
  if (lambda_val <= 0 || alpha <= 0 || alpha >= 1) return;

  double min_gm = gm.min();
  double mu = lambda_val * (1 - alpha) * min_gm;
  bedpp_mu = mu;
  if (mu <= 0) return;

  double scale = sglasso_dual_scale(Xtilde, r, beta, K1, lambda_val, alpha,
                                    d_val, n, J, gm, xty);
  if (!std::isfinite(scale)) return;

  double primal = sglasso_primal_objective(Xtilde, Ytilde, r, beta, K1,
                                           lambda_val, alpha, d_val, n, J,
                                           gm, xty);
  double dual = sglasso_dual_objective(Xtilde, Ytilde, r, beta, K1,
                                       lambda_val, alpha, d_val, n, J,
                                       gm, xty, scale);
  double gap = primal - dual;
  if (gap < 0 && gap > -1e-8) gap = 0;
  if (gap < 0 || !std::isfinite(gap)) return;
  bedpp_gap = gap;

  double radius = std::sqrt(2.0 * gap / mu);
  bedpp_radius = radius;
  double tolerance = 1e-8;
  for (arma::uword g = 0; g < J; g++) {
    if (ever_active[g] != 0 || ever_strong[g] != 0) continue;
    double lambda2 = lambda_val * gm[g] * (1 - alpha);
    double score_norm_sq = 0;
    for (arma::uword j = K1[g]; j < K1[g + 1]; j++) {
      double z = dot_col(Xtilde, r, j) / n;
      double score = z + lambda2 * (d_val * xty[j] - beta[j]);
      score_norm_sq += score * score;
    }
    double bj_bound;
    if (exact_bound) {
      bj_bound = exact_group_bedpp_bound(Xtilde, K1, g, lambda2, n);
    } else {
      double group_x_fro_bound = std::sqrt(K[g]) * x_fro_norm / n;
      bj_bound = group_x_fro_bound + lambda2;
    }
    double safe_bound = std::sqrt(score_norm_sq) + radius * bj_bound;
    double margin = lambda_val * alpha * gm[g] - safe_bound;
    if (!std::isfinite(bedpp_min_margin) || margin < bedpp_min_margin) {
      bedpp_min_margin = margin;
    }
    if (safe_bound + tolerance < lambda_val * alpha * gm[g]) {
      bedpp_safe[g] = 1;
    }
  }
}

double group_design_spectral_norm_sq(const arma::mat& Xtilde,
                                     const arma::colvec& K1,
                                     const arma::uword& g,
                                     const double& n) {
  arma::uword start = static_cast<arma::uword>(K1[g]);
  arma::uword end = static_cast<arma::uword>(K1[g + 1] - 1);
  arma::mat gram = Xtilde.cols(start, end).t() * Xtilde.cols(start, end) / n;
  return arma::eig_sym(gram).max();
}

void pathwise_bedpp_sglasso(const arma::mat& Xtilde,
                            const arma::colvec& r,
                            const arma::colvec& beta,
                            const arma::colvec& K,
                            const arma::colvec& K1,
                            const double& lambda_val,
                            const double& lambda_ref,
                            const double& alpha,
                            const double& d_val,
                            const double& n,
                            const double& J,
                            const arma::colvec& gm,
                            const arma::colvec& xty,
                            const arma::colvec& ever_active,
                            const arma::colvec& ever_strong,
                            double& bedpp_gap,
                            double& bedpp_radius,
                            double& bedpp_mu,
                            double& bedpp_min_margin,
                            arma::colvec& bedpp_safe) {
  bedpp_safe.zeros();
  bedpp_gap = arma::datum::nan;
  bedpp_radius = arma::datum::nan;
  bedpp_mu = arma::datum::nan;
  bedpp_min_margin = arma::datum::nan;
  if (lambda_val <= 0 || lambda_ref <= 0 || alpha <= 0 || alpha >= 1) return;

  double kappa = 1 - alpha;
  bedpp_mu = lambda_val * kappa * gm.min();

  double radius_sq = 0;
  for (arma::uword g = 0; g < J; g++) {
    double root_new = std::sqrt(lambda_val * kappa * gm[g]);
    double root_ref = std::sqrt(lambda_ref * kappa * gm[g]);
    double root_diff = root_new - root_ref;
    for (arma::uword j = K1[g]; j < K1[g + 1]; j++) {
      double eta = d_val * xty[j] - beta[j];
      radius_sq += root_diff * root_diff * eta * eta;
    }
  }
  double radius = std::sqrt(radius_sq);
  bedpp_radius = radius;
  bedpp_gap = radius_sq;

  double tolerance = 1e-8;
  for (arma::uword g = 0; g < J; g++) {
    if (ever_active[g] != 0 || ever_strong[g] != 0) continue;
    double lambda2 = lambda_val * gm[g] * kappa;
    double score_norm_sq = 0;
    for (arma::uword j = K1[g]; j < K1[g + 1]; j++) {
      double z = dot_col(Xtilde, r, j) / n;
      double score = z + lambda2 * (d_val * xty[j] - beta[j]);
      score_norm_sq += score * score;
    }
    double aj_norm = std::sqrt(group_design_spectral_norm_sq(Xtilde, K1, g, n) + lambda2);
    double safe_bound = std::sqrt(score_norm_sq) + radius * aj_norm;
    double margin = lambda_val * alpha * gm[g] - safe_bound;
    if (!std::isfinite(bedpp_min_margin) || margin < bedpp_min_margin) {
      bedpp_min_margin = margin;
    }
    if (margin > tolerance) {
      bedpp_safe[g] = 1;
    }
  }
}


// Scan for violations in strong set
bool Check_Strong_Set(const arma::mat& Xtilde, 
                      const arma::colvec& Ytilde,
                      const arma::colvec& r, 
                      const arma::colvec& K1,
                      const arma::colvec& K,
                      const double& lambda_val,
                      const double& alpha,
                      const double& d_val,
                      const double& n,
                      const double& J, 
                      const arma::colvec& gm,
                      const arma::colvec& xty,
                      const arma::colvec& beta_new,
                      arma::colvec& xTr,
                      arma::colvec& ever_strong, 
                      arma::colvec& ever_active,
                      std::vector<arma::uword>& active_groups,
                      int& kkt_checked,
                      int& kkt_violations){
  // write some..
  bool violations = false;
  for(arma::uword g = 0; g < J; g++){
    double lambda2 = lambda_val * gm[g] * (1-alpha);
    if(ever_active[g] == 0 && ever_strong[g] == 1){
      kkt_checked += 1;
      double LZ_norm_sq = 0;
      for(arma::uword j = K1[g]; j < K1[g+1]; j++){
        double z = dot_col(Xtilde, r, j) / n;
        double LZ = z + lambda2 * (d_val * xty[j] - beta_new[j]);
        LZ_norm_sq += LZ * LZ;
      }
      xTr[g] = std::sqrt(LZ_norm_sq);
      if (xTr[g]  > lambda_val * gm[g] * alpha) {
        activate_group(g, ever_active, active_groups);
        kkt_violations += 1;
        violations = true;
        }
      }
  }
  return(violations);
}

double max_inactive_kkt_residual(const arma::mat& Xtilde,
                                 const arma::colvec& r,
                                 const arma::colvec& K1,
                                 const double& lambda_val,
                                 const double& alpha,
                                 const double& d_val,
                                 const double& n,
                                 const double& J,
                                 const arma::colvec& gm,
                                 const arma::colvec& xty,
                                 const arma::colvec& beta_new,
                                 const arma::colvec& ever_active) {
  double max_residual = 0;
  for (arma::uword g = 0; g < J; g++) {
    if (ever_active[g] != 0) continue;
    double lambda2 = lambda_val * gm[g] * (1 - alpha);
    double LZ_norm_sq = 0;
    for (arma::uword j = K1[g]; j < K1[g + 1]; j++) {
      double z = dot_col(Xtilde, r, j) / n;
      double LZ = z + lambda2 * (d_val * xty[j] - beta_new[j]);
      LZ_norm_sq += LZ * LZ;
    }
    double residual = std::sqrt(LZ_norm_sq) - lambda_val * gm[g] * alpha;
    if (residual > max_residual) max_residual = residual;
  }
  return max_residual;
}

arma::uword count_nonzero_groups(const arma::colvec& beta,
                                 const arma::colvec& K1,
                                 const double& J) {
  arma::uword out = 0;
  for (arma::uword g = 0; g < J; g++) {
    bool any_nonzero = false;
    for (arma::uword j = K1[g]; j < K1[g + 1]; j++) {
      if (beta[j] != 0) {
        any_nonzero = true;
        break;
      }
    }
    if (any_nonzero) out += 1;
  }
  return out;
}

void activate_nonzero_groups(const arma::colvec& beta,
                             const arma::colvec& K1,
                             const double& J,
                             arma::colvec& ever_active,
                             std::vector<arma::uword>& active_groups) {
  for (arma::uword g = 0; g < J; g++) {
    for (arma::uword j = K1[g]; j < K1[g + 1]; j++) {
      if (beta[j] != 0) {
        activate_group(g, ever_active, active_groups);
        break;
      }
    }
  }
}



// Group descent update for sglasso
void gd_sglasso(const arma::mat& Xtilde, 
                const arma::uword g,
                const arma::uword d_indx,
                const double lambda1,
                const double lambda2,
                const arma::colvec& d,
                const arma::colvec& K,
                const arma::colvec& K1,
                const double n,
                const arma::colvec& xty,
                arma::colvec& beta,
                arma::colvec& r,
                arma::colvec& LZ,
                double& max_delta){
  // write some
  arma::uword start = static_cast<arma::uword>(K1[g]);
  arma::uword end = static_cast<arma::uword>(K1[g+1] - 1);
  arma::uword group_size = end - start + 1;
  double d_lambda2 = d[d_indx] * lambda2;
  double LZ_norm_sq = 0;
  const double* rr = r.memptr();
  bool has_nonzero_old = false;

  for (arma::uword offset = 0; offset < group_size; offset++) {
    arma::uword j = start + offset;
    const double* xj = Xtilde.colptr(j);
    double xr = 0;
    for (arma::uword i = 0; i < Xtilde.n_rows; i++) {
      xr += xj[i] * rr[i];
    }
    double z = beta[j] + xr / n;
    LZ[offset] = z + d_lambda2 * xty[j];
    LZ_norm_sq += LZ[offset] * LZ[offset];
    if (beta[j] != 0) {
      has_nonzero_old = true;
    }
  }

  double LZ_norm = std::sqrt(LZ_norm_sq);
  double len = S_c(LZ_norm,lambda1)/(1+lambda2);

  if(len != 0 || has_nonzero_old){
    for (arma::uword offset = 0; offset < group_size; offset++) {
      arma::uword j = start + offset;
      double beta_new_j = 0;
      if (LZ_norm > 0) {
        beta_new_j = len * LZ[offset] / LZ_norm;
      }
      double delta = beta_new_j - beta[j];
      if (delta != 0) {
        double abs_delta = std::abs(delta);
        if (abs_delta > max_delta) max_delta = abs_delta;
        beta[j] = beta_new_j;
        add_scaled_col(Xtilde, j, -delta, r);
      }
    }
  }
  }

double df_sglasso_group(const arma::mat& Xtilde,
                        const arma::uword g,
                        const arma::uword d_indx,
                        const double lambda1,
                        const double lambda2,
                        const arma::colvec& d,
                        const arma::colvec& K,
                        const arma::colvec& K1,
                        const double n,
                        const arma::colvec& xty,
                        const arma::colvec& beta,
                        const arma::colvec& r,
                        arma::colvec& LZ) {
  arma::uword start = static_cast<arma::uword>(K1[g]);
  arma::uword end = static_cast<arma::uword>(K1[g + 1] - 1);
  arma::uword group_size = end - start + 1;
  double d_lambda2 = d[d_indx] * lambda2;
  double LZ_norm_sq = 0;
  const double* rr = r.memptr();

  for (arma::uword offset = 0; offset < group_size; offset++) {
    arma::uword j = start + offset;
    const double* xj = Xtilde.colptr(j);
    double xr = 0;
    for (arma::uword i = 0; i < Xtilde.n_rows; i++) {
      xr += xj[i] * rr[i];
    }
    double z = beta[j] + xr / n;
    LZ[offset] = z + d_lambda2 * xty[j];
    LZ_norm_sq += LZ[offset] * LZ[offset];
  }

  double LZ_norm = std::sqrt(LZ_norm_sq);
  if (LZ_norm <= 0) return 0;
  double len = S_c(LZ_norm, lambda1) / (1 + lambda2);
  if (len <= 0) return 0;
  return K[g] * len / LZ_norm;
}



// [[Rcpp::export]]
List gd_sglasso_ssr(const arma::mat& Xtilde, 
                    const arma::colvec& Ytilde, 
                    const arma::colvec lambda,
                    const double lambda_max,
                    const arma::colvec d, 
                    const double alpha, 
                    const arma::colvec& K,
                    const arma::colvec& K1,
                    const arma::uword& K0,
                    const arma::colvec& gm,
                    const arma::colvec& beta_start, 
                    const double max_iter, 
                    const double eps,
                    const arma::uword& dfmax,
                    const arma::uword& gmax,
                    const int screen_rule,
                    const bool diagnostics){
  // Initiliaze variables
  arma::uword n = Xtilde.n_rows;
  arma::uword p = Xtilde.n_cols;
  arma::uword J = K1.n_rows-1;
  arma::uword D = d.n_elem;
  arma::uword L = lambda.n_elem;
  arma::colvec xty = Xtilde.t() * Ytilde / n;
  double x_fro_norm = arma::norm(Xtilde, "fro");
  // Outcome
  // arma::cube eta_mat  = cube(n, L,D,fill::zeros);
  arma::cube beta_mat = cube(p, L,D,fill::zeros);
  arma::mat df   = zeros(L, D);
  arma::mat loss = zeros(L, D);
  arma::mat iter = zeros(L, D);
  arma::mat rejections = zeros(L, D);
  arma::mat screen_stats = zeros(L * D, 33);
  
  



  
  // Path for d
  for(arma::uword d_indx = 0; d_indx < D; d_indx++){
  
  // Intermediate quantities
    arma::colvec r = Ytilde - Xtilde * beta_start;
    arma::colvec beta = beta_start;
    arma::colvec LZ_scratch = colvec(static_cast<arma::uword>(K.max()), fill::zeros);

    
    // Path for lambda
    for(arma::uword l_indx = 0; l_indx < L; l_indx++){
//      R_CheckUserInterrupt();
      std::chrono::steady_clock::time_point lambda_start;
      if (diagnostics) lambda_start = std::chrono::steady_clock::now();
      arma::uword iteration = 0;
      int kkt_checked = 0;
      int kkt_violations = 0;
      int strong_kkt_checked = 0;
      int strong_kkt_violations = 0;
      int rest_kkt_checked = 0;
      int rest_kkt_violations = 0;
      int n_refits = 0;
      int bedpp_safe_count = 0;
      double screening_time = 0;
      double bedpp_time = 0;
      double bedpp_gap = arma::datum::nan;
      double bedpp_radius = arma::datum::nan;
      double bedpp_mu = arma::datum::nan;
      double bedpp_min_margin = arma::datum::nan;
      double bcd_time = 0;
      double kkt_time = 0;
      arma::colvec xTr = colvec(J, fill::zeros);
      arma::colvec ever_strong = colvec(J, fill::zeros);
      arma::colvec ever_active = colvec(J, fill::zeros);
      arma::colvec bedpp_safe = colvec(J, fill::zeros);
      std::vector<arma::uword> active_groups;
      activate_nonzero_groups(beta, K1, J, ever_active, active_groups);


      // ssr screening
      std::chrono::steady_clock::time_point screen_start;
      if (diagnostics) screen_start = std::chrono::steady_clock::now();
      ssr_sglasso(K,K1,Xtilde, r,gm,lambda,lambda_max, l_indx,alpha, d[d_indx], n, J, xty, beta, ever_active, screen_rule, xTr, ever_strong);
      if (diagnostics) {
        auto screen_end = std::chrono::steady_clock::now();
        screening_time += std::chrono::duration<double>(screen_end - screen_start).count();
      }
      rejections(l_indx,d_indx) = J - sum(ever_strong);

      if (screen_rule == 3 || screen_rule == 4) {
        std::chrono::steady_clock::time_point bedpp_start;
        if (diagnostics) bedpp_start = std::chrono::steady_clock::now();
        bedpp_sglasso_safe(Xtilde, Ytilde, r, beta, K, K1,
                           lambda[l_indx], alpha, d[d_indx], n, J, gm,
                           xty, ever_active, ever_strong, x_fro_norm,
                           screen_rule == 4,
                           bedpp_gap, bedpp_radius, bedpp_mu,
                           bedpp_min_margin,
                           bedpp_safe);
        bedpp_safe_count = static_cast<int>(sum(bedpp_safe));
        if (diagnostics) {
          auto bedpp_end = std::chrono::steady_clock::now();
          bedpp_time += std::chrono::duration<double>(bedpp_end - bedpp_start).count();
        }
      }
      if (screen_rule == 5) {
        std::chrono::steady_clock::time_point bedpp_start;
        if (diagnostics) bedpp_start = std::chrono::steady_clock::now();
        double lambda_ref = (l_indx != 0) ? lambda[l_indx - 1] : lambda_max;
        if (lambda_ref <= 0) lambda_ref = lambda[l_indx];
        pathwise_bedpp_sglasso(Xtilde, r, beta, K, K1,
                               lambda[l_indx], lambda_ref, alpha, d[d_indx],
                               n, J, gm, xty, ever_active, ever_strong,
                               bedpp_gap, bedpp_radius, bedpp_mu,
                               bedpp_min_margin, bedpp_safe);
        bedpp_safe_count = static_cast<int>(sum(bedpp_safe));
        if (diagnostics) {
          auto bedpp_end = std::chrono::steady_clock::now();
          bedpp_time += std::chrono::duration<double>(bedpp_end - bedpp_start).count();
        }
      }
      
      while (iteration < max_iter) {
        while (iteration < max_iter) {
          while (iteration < max_iter) {
            bool converged = false;
            double max_delta = 0;
            iteration += 1;
            iter(l_indx,d_indx) +=1;
            // Update only active groups. Keep the list incrementally to avoid
            // scanning all groups on every BCD sweep.
            std::chrono::steady_clock::time_point bcd_start;
            if (diagnostics) bcd_start = std::chrono::steady_clock::now();
            for (arma::uword active_indx = 0; active_indx < active_groups.size(); active_indx++) {
              arma::uword g = active_groups[active_indx];
              double lambda1 = lambda[l_indx] * gm[g] * alpha;
              double lambda2 = lambda[l_indx] * gm[g] * (1-alpha);
              gd_sglasso(Xtilde, g,d_indx,lambda1,lambda2,d,K,K1,n,xty,beta,
                         r,LZ_scratch,max_delta);
            }
            if (diagnostics) {
              auto bcd_end = std::chrono::steady_clock::now();
              bcd_time += std::chrono::duration<double>(bcd_end - bcd_start).count();
            }
            
            // Check convergence
            if(max_delta <= eps){
              converged = true;
              // update loss
              loss(l_indx,d_indx) = gLoss(r, n);
              if(converged) break;
            }
          } // end for first while loop
          
          
          // Scan for violations in strong set
          std::chrono::steady_clock::time_point kkt_start;
          if (diagnostics) kkt_start = std::chrono::steady_clock::now();
          int checked_before = kkt_checked;
          int violations_before = kkt_violations;
          bool violations = Check_Strong_Set(Xtilde, Ytilde, r, K1, K, lambda[l_indx], alpha, d[d_indx],n, J, gm, xty, beta, xTr, ever_strong, ever_active, active_groups, kkt_checked, kkt_violations);
          strong_kkt_checked += kkt_checked - checked_before;
          strong_kkt_violations += kkt_violations - violations_before;
          if (diagnostics) {
            auto kkt_end = std::chrono::steady_clock::now();
            kkt_time += std::chrono::duration<double>(kkt_end - kkt_start).count();
          }
          if (violations) n_refits += 1;
          if (violations == false) break;
        } // end for second while loop
        
        // Exact modes check the full rest set at every lambda. The fast mode
        // checks only strategic path checkpoints: first, Q1, median, Q3, last.
        arma::uword q1_indx = (L > 1) ? static_cast<arma::uword>(std::floor((L - 1) * 0.25)) : 0;
        arma::uword q2_indx = (L > 1) ? static_cast<arma::uword>(std::floor((L - 1) * 0.50)) : 0;
        arma::uword q3_indx = (L > 1) ? static_cast<arma::uword>(std::floor((L - 1) * 0.75)) : 0;
        bool fast_checkpoint = (l_indx == 0) ||
          (l_indx == q1_indx) ||
          (l_indx == q2_indx) ||
          (l_indx == q3_indx) ||
          (l_indx + 1 == L);
        bool do_rest_kkt = (screen_rule != 2) || fast_checkpoint;
        if (do_rest_kkt) {
          std::chrono::steady_clock::time_point kkt_start;
          if (diagnostics) kkt_start = std::chrono::steady_clock::now();
          int checked_before = kkt_checked;
          int violations_before = kkt_violations;
          bool violations = Check_Rest_Set(Xtilde, Ytilde, r, K1, K, lambda[l_indx], alpha, d[d_indx],n, J, gm, xty, beta, xTr, ever_strong, ever_active, active_groups, kkt_checked, kkt_violations, bedpp_safe);
          rest_kkt_checked += kkt_checked - checked_before;
          rest_kkt_violations += kkt_violations - violations_before;
          if (diagnostics) {
            auto kkt_end = std::chrono::steady_clock::now();
            kkt_time += std::chrono::duration<double>(kkt_end - kkt_start).count();
          }
          if (violations) n_refits += 1;
          if (violations == false) break;
        } else {
          break;
        }
      } // end for third while loop

      // Store the final coefficients for this lambda/d pair once all BCD and
      // KKT passes have finished.
      df(l_indx,d_indx) = 0;
      for (arma::uword active_indx = 0; active_indx < active_groups.size(); active_indx++) {
        arma::uword g = active_groups[active_indx];
        double lambda1 = lambda[l_indx] * gm[g] * alpha;
        double lambda2 = lambda[l_indx] * gm[g] * (1-alpha);
        df(l_indx,d_indx) += df_sglasso_group(Xtilde, g, d_indx, lambda1,
                                              lambda2, d, K, K1, n, xty,
                                              beta, r, LZ_scratch);
      }
      beta_mat.slice(d_indx).col(l_indx) = beta;
      double elapsed = 0;
      if (diagnostics) {
        auto lambda_end = std::chrono::steady_clock::now();
        elapsed = std::chrono::duration<double>(lambda_end - lambda_start).count();
      }
      arma::uword stats_row = d_indx * L + l_indx;
      screen_stats(stats_row, 0) = lambda[l_indx];
      screen_stats(stats_row, 1) = d[d_indx];
      screen_stats(stats_row, 2) = J;
      screen_stats(stats_row, 3) = count_nonzero_groups(beta, K1, J);
      screen_stats(stats_row, 4) = active_groups.size();
      screen_stats(stats_row, 5) = rejections(l_indx, d_indx);
      screen_stats(stats_row, 6) = kkt_checked;
      screen_stats(stats_row, 7) = kkt_violations;
      screen_stats(stats_row, 8) = elapsed;
      screen_stats(stats_row, 9) = bcd_time;
      screen_stats(stats_row, 10) = screening_time;
      screen_stats(stats_row, 11) = kkt_time;
      screen_stats(stats_row, 12) = iter(l_indx, d_indx);
      screen_stats(stats_row, 13) = screen_rule;
      screen_stats(stats_row, 14) = d_indx + 1;
      screen_stats(stats_row, 15) = diagnostics ?
        max_inactive_kkt_residual(Xtilde, r, K1, lambda[l_indx],
                                  alpha, d[d_indx], n, J, gm,
                                  xty, beta, ever_active) :
        arma::datum::nan;
      double kept_after_screening = J - rejections(l_indx, d_indx);
      screen_stats(stats_row, 16) = kept_after_screening;
      screen_stats(stats_row, 17) = J > 0 ? rejections(l_indx, d_indx) / J : 0;
      screen_stats(stats_row, 18) = J > 0 ? kept_after_screening / J : 0;
      screen_stats(stats_row, 19) = kkt_violations > 0 ? 1 : 0;
      screen_stats(stats_row, 20) = n_refits;
      screen_stats(stats_row, 21) = strong_kkt_checked;
      screen_stats(stats_row, 22) = strong_kkt_violations;
      screen_stats(stats_row, 23) = rest_kkt_checked;
      screen_stats(stats_row, 24) = rest_kkt_violations;
      screen_stats(stats_row, 25) = bedpp_safe_count;
      screen_stats(stats_row, 26) = J > 0 ? bedpp_safe_count / static_cast<double>(J) : 0;
      screen_stats(stats_row, 27) = bedpp_time;
      screen_stats(stats_row, 28) = rest_kkt_checked + bedpp_safe_count;
      screen_stats(stats_row, 29) = bedpp_gap;
      screen_stats(stats_row, 30) = bedpp_radius;
      screen_stats(stats_row, 31) = bedpp_mu;
      screen_stats(stats_row, 32) = bedpp_min_margin;
        
   
      } // ends for lambda index 
    } // ends for d index 
  // Return elements of the list
  return List::create(_["beta_mat"]=beta_mat,
                      _["df"]=df,
	                      _["loss"]=loss,
	                      _["iter"]=iter,
	                      _["rejections"]=rejections,
	                      _["screen_stats"]=screen_stats);
  
  } // List gd_sglasso_ssr
