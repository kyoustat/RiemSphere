#include "RcppArmadillo.h"
#include "riemfactory.hpp"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RiemBase)]]

using namespace Rcpp;
using namespace arma;


double engine_mean_eval(arma::mat tgt, arma::cube data, std::string name, arma::vec weights){
  const int N = data.n_slices;
  
  double output = 0.0;
  double tmpout = 0.0;
  for (int i=0;i<N;i++){
    tmpout = riemfunc_dist(tgt, data.slice(i), name);
    output += tmpout*tmpout*weights(i);
  }
  return(output);
}
double engine_stepsize(arma::mat mold, arma::mat dtmp, arma::cube data, std::string name, arma::vec weights){
  const int np = 20;
  arma::vec steps = arma::linspace<vec>(0, 5, np);
  arma::vec fvals(np,fill::zeros);
  arma::mat mnew; mnew.copy_size(mold); 
  fvals(0)   = engine_mean_eval(mold, data, name, weights);
  for (int i=1;i<np;i++){
    mnew     = riemfunc_exp(mold, dtmp, steps(i), name);
    fvals(i) = engine_mean_eval(mnew, data, name, weights);
  }
  double optval = steps(fvals.index_min());
  return(optval);
}


arma::mat engine_extrinsicmean(arma::cube data, std::string name){
  int dnrow = data.n_rows;
  int dncol = data.n_cols;
  int nslice = data.n_slices;
  
  arma::vec testdata = riemfunc_equiv(data.slice(0),dnrow,dncol,name);
  int L = testdata.n_elem; 
  
  arma::mat Xext(L,nslice,fill::zeros);
  for (int i=0;i<nslice;i++){
    Xext.col(i) = riemfunc_equiv(data.slice(i),dnrow,dncol,name);
  }
  
  arma::vec extmean = arma::mean(Xext, 1); // find mean for each row.
  arma::mat inveqmu = riemfunc_invequiv(extmean,dnrow,dncol,name); // projected mean
  return(inveqmu);
}

// [[Rcpp::export]]
Rcpp::List engine_wmean(arma::cube data, std::string name, int maxiter, double eps, arma::vec weights){
  // get parameters
  int N = data.n_slices;
  int iter = 0;
  
  // initialize
  arma::mat mold = engine_extrinsicmean(data, name); // extrinsic mean as an initializer
  arma::mat mnew;   mnew.copy_size(mold);  mnew.fill(0); 
  arma::cube tvecs; tvecs.copy_size(data); tvecs.fill(0);
  
  arma::mat dtmp; dtmp.copy_size(mold);  dtmp.fill(0); // on TpM
  
  // let's iterate !
  double sqnorm = 10000.00;
  double opstep = 0.0;
  for (int it=0;it<maxiter;it++){
    // 0. update iter
    iter += 1;
    // 1. compute log-pulled vectors & 2. adjust weights
    dtmp.fill(0);
    for (int i=0;i<N;i++){
      dtmp += riemfunc_log(mold, data.slice(i), name)*weights(i);
    }
    dtmp /= arma::norm(dtmp,"fro");
    // 3. update using exponential map and compute
    opstep = engine_stepsize(mold, dtmp, data, name, weights); // sort of greed search
    if (opstep < 1e-10){
      break;
    }
    mnew   = riemfunc_exp(mold, dtmp, opstep, name);
    
    // 4. iteration : update sqnorm
    sqnorm = riemfunc_dist(mold,mnew,name);
    if ((sqnorm < eps)&&(it > 10)){
      break;
    }
    mold = mnew;
  }
  
  return(Rcpp::List::create(Rcpp::Named("x")=mold,
                            Rcpp::Named("iteration")=iter));
}