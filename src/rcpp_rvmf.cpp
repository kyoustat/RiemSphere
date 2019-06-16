#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
NumericVector rvmf_h(double n,double ca,double d1,double x0,double m,double k,double b){
  double ta,u,z,tmp=0;
  NumericVector w(n);
  for(int i=0;i<n;++i) {
    for(ta=-1000.0,u=1.0;ta - ca < std::log(u);) {
      z = R::rbeta(m,m);
      u = R::runif(0,1);
      tmp = ( 1 - (1 + b) * z ) / ( 1 - (1 - b) * z );
      ta = k * tmp + d1 * static_cast<double>(std::log(1 - x0 * tmp));
    }
    w[i]=tmp;
  }
  return w;
}