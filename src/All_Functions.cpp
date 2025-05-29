#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;




// Check for convergence of any beta vector
bool checkConvergence(const arma::colvec& beta_new, 
                      const arma::colvec& beta_old, 
                      const double& eps,
                      const double& J) {
  bool converged = true;
  for (arma::uword j=0; j < J; j++) {
    if (beta_new[j]!=0 && beta_old[j]!=0) {
      if (abs(beta_new[j]-beta_old[j]) > eps) {
        return(false);
        break;
      }
    } else if (beta_new[j]==0 && beta_old[j]!=0) {
      return(false);
      break;
    } else if (beta_new[j]!=0 && beta_old[j]==0) {
      return(false);
      break;
    }
  }
  return(converged);
}


// Gaussian loss function
double gLoss(const arma::colvec r, int n) {
  double val = 0;
  for (int i=0;i<n;i++) val += pow(r[i],2);
  return(val);
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
    xx.col(j) /= s[j];
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
   pen =  pen + lambda*gm[g]*alpha+abs(beta[j]) + .5*lambda*gm[g]*(1-alpha)*pow((beta[j]-d*beta_lse[j]),2);
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
                const arma::colvec& beta_new,
                arma::colvec& xTr,
                arma::colvec& ever_strong) {
  for(int g = 0; g < J; g++){
    arma::colvec z = colvec(K[g], fill::zeros);
    for(int j = K1[g]; j < K1[g+1]; j++){
      z[j-K1[g]] =  arma::as_scalar(Xtilde.col(j).t() * r / n);
    }
    xTr[g] = norm(z,2);
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
      arma::colvec z = colvec(K[g], fill::zeros);
      arma::colvec LZ = colvec(K[g], fill::zeros);
      arma::colvec beta_ols_g = colvec(K[g], fill::zeros);
      for(arma::uword j = K1[g]; j < K1[g+1]; j++){
        z[j-K1[g]] =  arma::as_scalar(Xtilde.col(j).t() * r / n);
        beta_ols_g[j-K1[g]] = arma::as_scalar(Xtilde.col(j).t() * Ytilde / n);
        LZ[j-K1[g]] = z[j-K1[g]] + lambda2*(d_val*beta_ols_g[j-K1[g]] - beta_new[j-K1[g]]);
      }
      xTr[g] = norm(LZ,2);
      if (xTr[g]  +  TOLERANCE > alpha * lambda_val * gm[g]) {
        ever_active[g] = ever_strong[g] = 1;
        return(true);
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
                      const arma::colvec& beta_new,
                      arma::colvec& xTr,
                      arma::colvec& ever_strong, 
                      arma::colvec& ever_active){
  // write some..
  bool violations = false;
  for(arma::uword g = 0; g < J; g++){
    double lambda2 = lambda_val * gm[g] * (1-alpha);
    if(ever_active[g] == 0 && ever_strong[g] == 1){
      arma::colvec z = colvec(K[g], fill::zeros);
      arma::colvec LZ = colvec(K[g], fill::zeros);
      arma::colvec beta_ols_g = colvec(K[g], fill::zeros);
      for(arma::uword j = K1[g]; j < K1[g+1]; j++){
        z[j-K1[g]] =  arma::as_scalar(Xtilde.col(j).t() * r / n);
        beta_ols_g[j-K1[g]] = arma::as_scalar(Xtilde.col(j).t() * Ytilde / n);
        LZ[j-K1[g]] = z[j-K1[g]] + lambda2*(d_val*beta_ols_g[j-K1[g]] - beta_new[j-K1[g]]);
      }
      xTr[g] = norm(LZ,2);
      if (xTr[g]  > lambda_val * gm[g] * alpha) {
        ever_active[g] = 1;
        return(true);
        }
      }
  }
  return(violations);
}



// Group descent update for sglasso
void gd_sglasso(const arma::mat& Xtilde, 
                const arma::colvec& Ytilde, 
                const double g,
                const double l_indx,
                const double lambda1,
                const double lambda2,
                const arma::colvec d,
                const double d_indx,
                const arma::colvec& K,
                const arma::colvec& K1,
                const double n,
                const arma::colvec& gm,
                arma::colvec& beta_last,
                arma::colvec& beta_new,
                arma::colvec& r,
                arma::mat& df){
  // write some
  arma::colvec z = colvec(K[g], fill::zeros);
  arma::colvec LZ = colvec(K[g], fill::zeros);
  arma::colvec beta_ols_g = colvec(K[g], fill::zeros);
  for(arma::uword j = K1[g]; j < K1[g+1]; j++){
    z[j-K1[g]] = arma::as_scalar(beta_last[j] + ((Xtilde.col(j).t() * r) / n));
    beta_ols_g[j-K1[g]] = arma::as_scalar(Xtilde.col(j).t() * Ytilde / n);
    LZ[j-K1[g]] = z[j-K1[g]] + d[d_indx]*lambda2*beta_ols_g[j-K1[g]];
    }
  double LZ_norm = norm(LZ,2);
  double len = S_c(LZ_norm,lambda1)/(1+lambda2);
   if(len !=0 || beta_last[K1[g]] !=0){
    for(arma::uword j = K1[g]; j < K1[g+1]; j++){
      // Update beta
      beta_new[j] = len * LZ[j-K1[g]]/ LZ_norm;
      // Update Partial residual
      r = r - Xtilde.col(j) * (beta_new[j] - beta_last[j]);
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
  arma::uword iteration = 0;
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

    
    // Path for lambda
    for(arma::uword l_indx = 0; l_indx < L; l_indx++){
//      R_CheckUserInterrupt();


      // ssr screening
      ssr_sglasso(K,K1,Xtilde, r,gm,lambda,lambda_max, l_indx,alpha, n, J,beta_last, xTr, ever_strong);
      rejections(l_indx,d_indx) = J - sum(ever_strong);
      
      while (iteration < max_iter) {
        while (iteration < max_iter) {
          while (iteration < max_iter) {
            bool converged = false;
            iteration += 1;
            iter(l_indx,d_indx) +=1;
            df(l_indx,d_indx) = 0;


            
            // update for penalized groups
            for (arma::uword g=0; g<J; g++) {
              double lambda1 = lambda[l_indx] * gm[g] * alpha;
              double lambda2 = lambda[l_indx] * gm[g] * (1-alpha);
              if (ever_active[g]) {
                gd_sglasso(Xtilde, Ytilde, g,l_indx,lambda1,lambda2,d,d_indx,K,K1,n,gm,beta_last,
                           beta_new,r,df);
              }
            }
            
            // update beta values
            beta_last = beta_new;
            // store new beta to beta_mat
            beta_mat.slice(d_indx).col(l_indx) = beta_last;
            
            
            // Check convergence
            if(checkConvergence(beta_new, beta_last, eps,J)){
              converged = true;
              // update loss
              loss(l_indx,d_indx) = gLoss(r, n);
              if(converged) break;
            }
          } // end for first while loop
          
          
          // Scan for violations in strong set
          bool violations = Check_Strong_Set(Xtilde, Ytilde, r, K1, K, lambda[l_indx], alpha, d[d_indx],n, J, gm, beta_last, xTr, ever_strong, ever_active);
          if (violations == false) break;
        } // end for second while loop
        
        // Scan for violations in strong set
        bool violations = Check_Rest_Set(Xtilde, Ytilde, r, K1, K, lambda[l_indx], alpha, d[d_indx],n, J, gm, beta_last, xTr, ever_strong, ever_active);
        if (violations == false) break;
      } // end for third while loop
        
   
      } // ends for lambda index 
    } // ends for d index 
  // Return elements of the list
  return List::create(_["beta_mat"]=beta_mat,
                      _["df"]=df,
                      _["loss"]=loss,
                      _["iter"]=iter,
                      _["rejections"]=rejections);
  
  } // List gd_sglasso_ssr