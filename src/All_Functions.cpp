#include <RcppArmadillo.h>
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
                            const arma::uvec& active_groups) {
  for (arma::uword active_indx = 0; active_indx < active_groups.n_elem; active_indx++) {
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
  // Initialize parameters
  arma::mat X = Xtilde;
  arma::colvec Y = Ytilde;
  arma::colvec param = beta;
  // Calculate and return objective value
  double pen=0;
  for(int g = 0; g < J; g++){
    for(int j = K1[g]; j < K1[g+1]; j++){
   pen =  pen + lambda*gm[g]*alpha*abs(beta[j]) + .5*lambda*gm[g]*(1-alpha)*pow((beta[j]-d*beta_lse[j]),2);
    }
  }
  return accu(square(Y - X * param)/(2 * n)) + pen;
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
        double z = arma::as_scalar(abs((Xtilde.col(j).t() * Ytilde) / sqrt(K[g])));
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
                const double l_indx,
                const double alpha,
                const double n,
                const double J,
                arma::colvec& xTr,
                arma::colvec& ever_strong) {
  for(int g = 0; g < J; g++){
    double z_norm_sq = 0;
    for(arma::uword j = K1[g]; j < K1[g+1]; j++){
      double z = arma::dot(Xtilde.col(j), r) / n;
      z_norm_sq += z * z;
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
                    arma::colvec& ever_active){
  // Write some explanitions
  bool violations = false;
  double TOLERANCE = 1e-8;
  for(arma::uword g = 0; g < J; g++){
    double lambda2 = lambda_val * gm[g] * (1-alpha);
    if(ever_strong[g] == 0){
      double LZ_norm_sq = 0;
      for(arma::uword j = K1[g]; j < K1[g+1]; j++){
        double z = arma::dot(Xtilde.col(j), r) / n;
        double LZ = z + lambda2 * (d_val * xty[j] - beta_new[j]);
        LZ_norm_sq += LZ * LZ;
      }
      xTr[g] = std::sqrt(LZ_norm_sq);
      if (xTr[g]  +  TOLERANCE > alpha * lambda_val * gm[g]) {
        ever_active[g] = ever_strong[g] = 1;
        violations = true;
        }
      }
  }
  return(violations);
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
                      arma::colvec& ever_active){
  // write some..
  bool violations = false;
  for(arma::uword g = 0; g < J; g++){
    double lambda2 = lambda_val * gm[g] * (1-alpha);
    if(ever_active[g] == 0 && ever_strong[g] == 1){
      double LZ_norm_sq = 0;
      for(arma::uword j = K1[g]; j < K1[g+1]; j++){
        double z = arma::dot(Xtilde.col(j), r) / n;
        double LZ = z + lambda2 * (d_val * xty[j] - beta_new[j]);
        LZ_norm_sq += LZ * LZ;
      }
      xTr[g] = std::sqrt(LZ_norm_sq);
      if (xTr[g]  > lambda_val * gm[g] * alpha) {
        ever_active[g] = 1;
        violations = true;
        }
      }
  }
  return(violations);
}



// Group descent update for sglasso
void gd_sglasso(const arma::mat& Xtilde, 
                const arma::uword g,
                const arma::uword l_indx,
                const double lambda1,
                const double lambda2,
                const arma::colvec& d,
                const arma::uword d_indx,
                const arma::colvec& K,
                const arma::colvec& K1,
                const double n,
                const arma::colvec& xty,
                arma::colvec& beta_last,
                arma::colvec& beta_new,
                arma::colvec& r,
                arma::mat& df,
                arma::colvec& LZ){
  // write some
  arma::uword start = static_cast<arma::uword>(K1[g]);
  arma::uword end = static_cast<arma::uword>(K1[g+1] - 1);
  arma::uword group_size = end - start + 1;
  double d_lambda2 = d[d_indx] * lambda2;
  double LZ_norm_sq = 0;

  for (arma::uword offset = 0; offset < group_size; offset++) {
    arma::uword j = start + offset;
    double z = beta_last[j] + arma::dot(Xtilde.col(j), r) / n;
    LZ[offset] = z + d_lambda2 * xty[j];
    LZ_norm_sq += LZ[offset] * LZ[offset];
  }

  double LZ_norm = std::sqrt(LZ_norm_sq);
  double len = S_c(LZ_norm,lambda1)/(1+lambda2);

  bool has_nonzero_old = false;
  for (arma::uword j = start; j <= end; j++) {
    if (beta_last[j] != 0) {
      has_nonzero_old = true;
      break;
    }
  }

  if(len != 0 || has_nonzero_old){
    for (arma::uword offset = 0; offset < group_size; offset++) {
      arma::uword j = start + offset;
      double beta_new_j = 0;
      if (LZ_norm > 0) {
        beta_new_j = len * LZ[offset] / LZ_norm;
      }
      double delta = beta_new_j - beta_last[j];
      if (delta != 0) {
        beta_new[j] = beta_new_j;
        r -= Xtilde.col(j) * delta;
      }
    }
  }
  // Update df
  if(len > 0 ) df(l_indx,d_indx) +=  K[g] * len / LZ_norm;
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
                    const arma::uword& gmax){
  // Initiliaze variables
  arma::uword n = Xtilde.n_rows;
  arma::uword p = Xtilde.n_cols;
  arma::uword J = K1.n_rows-1;
  arma::uword D = d.n_elem;
  arma::uword L = lambda.n_elem;
  arma::colvec xty = Xtilde.t() * Ytilde / n;
  // Outcome
  // arma::cube eta_mat  = cube(n, L,D,fill::zeros);
  arma::cube beta_mat = cube(p, L,D,fill::zeros);
  arma::mat df   = zeros(L, D);
  arma::mat loss = zeros(L, D);
  arma::mat iter = zeros(L, D);
  arma::mat rejections = zeros(L, D);
  
  



  
  // Path for d
  for(arma::uword d_indx = 0; d_indx < D; d_indx++){
  
  // Intermediate quantities
    arma::colvec r = Ytilde - Xtilde * beta_start;
    arma::colvec beta_last = beta_start;
    arma::colvec beta_new = beta_start;
    arma::colvec xTr = colvec(J, fill::zeros);
    arma::colvec ever_strong = colvec(J, fill::zeros);
    arma::colvec ever_active = colvec(J, fill::zeros);
    arma::colvec LZ_scratch = colvec(static_cast<arma::uword>(K.max()), fill::zeros);

    
    // Path for lambda
    for(arma::uword l_indx = 0; l_indx < L; l_indx++){
//      R_CheckUserInterrupt();
      arma::uword iteration = 0;


      // ssr screening
      ssr_sglasso(K,K1,Xtilde, r,gm,lambda,lambda_max, l_indx,alpha, n, J, xTr, ever_strong);
      rejections(l_indx,d_indx) = J - sum(ever_strong);
      
      while (iteration < max_iter) {
        while (iteration < max_iter) {
          while (iteration < max_iter) {
            bool converged = false;
            iteration += 1;
            iter(l_indx,d_indx) +=1;
            df(l_indx,d_indx) = 0;


            
            // Update only active groups. This avoids checking every group on
            // every BCD sweep when SSR/KKT screening leaves a sparse active set.
            arma::uvec active_groups = arma::find(ever_active != 0);
            for (arma::uword active_indx = 0; active_indx < active_groups.n_elem; active_indx++) {
              arma::uword g = active_groups[active_indx];
              double lambda1 = lambda[l_indx] * gm[g] * alpha;
              double lambda2 = lambda[l_indx] * gm[g] * (1-alpha);
              gd_sglasso(Xtilde, g,l_indx,lambda1,lambda2,d,d_indx,K,K1,n,xty,beta_last,
                         beta_new,r,df,LZ_scratch);
            }
            
            // Check convergence
            if(checkConvergenceActive(beta_new, beta_last, eps, K1, active_groups)){
              converged = true;
              // update loss
              loss(l_indx,d_indx) = gLoss(r, n);
              beta_last = beta_new;
              if(converged) break;
            }
            beta_last = beta_new;
          } // end for first while loop
          
          
          // Scan for violations in strong set
          bool violations = Check_Strong_Set(Xtilde, Ytilde, r, K1, K, lambda[l_indx], alpha, d[d_indx],n, J, gm, xty, beta_last, xTr, ever_strong, ever_active);
          if (violations == false) break;
        } // end for second while loop
        
        // Scan for violations in strong set
        bool violations = Check_Rest_Set(Xtilde, Ytilde, r, K1, K, lambda[l_indx], alpha, d[d_indx],n, J, gm, xty, beta_last, xTr, ever_strong, ever_active);
        if (violations == false) break;
      } // end for third while loop

      // Store the final coefficients for this lambda/d pair once all BCD and
      // KKT passes have finished.
      beta_mat.slice(d_indx).col(l_indx) = beta_last;
        
   
      } // ends for lambda index 
    } // ends for d index 
  // Return elements of the list
  return List::create(_["beta_mat"]=beta_mat,
                      _["df"]=df,
                      _["loss"]=loss,
                      _["iter"]=iter,
                      _["rejections"]=rejections);
  
  } // List gd_sglasso_ssr
