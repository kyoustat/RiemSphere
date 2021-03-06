#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


/* Auxiliary, Elementary Computation
 * (1) cppdist_int_MtoN
 * (2) cppdist_ext_MtoN
 * (3) cppdist_pair_int
 * (4) cppdist_pair_ext
 * (5) cppdist_int_1toN
 * (6) cppdist_ext_1toN
 */

// (1) cppdist_int_MtoN cppdist_pair_int
//[[Rcpp::export]]
arma::mat cppdist_int_MtoN(arma::mat &X, arma::mat &Y){
  // parameters
  int M = X.n_rows;
  int N = Y.n_rows;
  
  // preliminary
  arma::mat output(M,N,fill::zeros);
  arma::rowvec xvec;
  arma::rowvec yvec;
  
  // iteration
  for (int m=0;m<M;m++){
    xvec = X.row(m);
    for (int n=0;n<N;n++){
      yvec = Y.row(n);
      if (arma::norm(xvec-yvec,2) > 1e-10){
        output(m,n) = acos(arma::dot(xvec, yvec));  
      } else {
        output(m,n) = 0;
      }
    }
  }
  
  // return
  return(output);
}

// (2) cppdist_ext_MtoN
//[[Rcpp::export]]
arma::mat cppdist_ext_MtoN(arma::mat &X, arma::mat &Y){
  // parameters
  int M = X.n_rows;
  int N = Y.n_rows;
  
  // preliminary
  arma::mat output(M,N,fill::zeros);
  arma::rowvec xvec;
  arma::rowvec yvec;
  
  // iteration
  for (int m=0;m<M;m++){
    xvec = X.row(m);
    for (int n=0;n<N;n++){
      yvec = Y.row(n);
      output(m,n) = arma::norm(xvec-yvec, 2);
    }
  }
  
  // return
  return(output);
}

// (3) cppdist_pair_int
//[[Rcpp::export]]
arma::mat cppdist_pair_int(arma::mat &X){
  // parameters
  int M = X.n_rows;
  
  // preliminary
  arma::mat output(M,M,fill::zeros);
  arma::rowvec x1;
  arma::rowvec x2;
  double dval = 0.0;
  
  // iteration
  for (int m=0;m<(M-1);m++){
    x1 = X.row(m);
    for (int n=(m+1);n<M;n++){
      x2 = X.row(n);
      if (arma::norm(x1-x2,2) > 1e-10){
        dval = acos(arma::dot(x1, x2));
      } else {
        dval = 0;
      }
      output(m,n) = dval;
      output(n,m) = dval;
    }
  }
  return(output);
}

// (4) cppdist_pair_ext
//[[Rcpp::export]]
arma::mat cppdist_pair_ext(arma::mat &X){
  // parameters
  int M = X.n_rows;
  
  // preliminary
  arma::mat output(M,M,fill::zeros);
  arma::rowvec x1;
  arma::rowvec x2;
  double dval = 0.0;
  
  // iteration
  for (int m=0;m<(M-1);m++){
    x1 = X.row(m);
    for (int n=(m+1);n<M;n++){
      x2 = X.row(n);
      
      dval = arma::norm(x1-x2, 2);
      output(m,n) = dval;
      output(n,m) = dval;
    }
  }
  return(output);
}

// (5) cppdist_int_1toN
//[[Rcpp::export]]
arma::vec cppdist_int_1toN(arma::vec x, arma::mat &Y){
  // parameters
  int N = Y.n_rows;
  
  // preliminary
  arma::vec output(N,fill::zeros);
  arma::rowvec xvec = x.t();
  arma::rowvec yvec;
  
  // iteration
  for (int n=0;n<N;n++){
    yvec = Y.row(n);
    if (arma::norm(xvec-yvec,2) > 1e-10){
      output(n) = acos(arma::dot(xvec, yvec));  
    } else {
      output(n) = 0;
    }
  }

  // return
  return(output);
}

// (6) cppdist_ext_1toN
//[[Rcpp::export]]
arma::vec cppdist_ext_1toN(arma::vec x, arma::mat &Y){
  // parameters
  int N = Y.n_rows;
  
  // preliminary
  arma::vec output(N,fill::zeros);
  arma::rowvec xvec = x.t();
  arma::rowvec yvec;
  
  // iteration
  for (int n=0;n<N;n++){
    yvec      = Y.row(n);
    output(n) = arma::norm(xvec-yvec, 2);
  }
  
  // return
  return(output);
}