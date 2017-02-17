#include <Rcpp.h>

using namespace Rcpp;

// declare
double ppexpC(int q, NumericVector rate, NumericVector t);

// [[Rcpp::export]]
NumericVector ppexpCV(int q, NumericMatrix rate, NumericVector t){
  std::vector<double> result(rate.nrow());
  for(int i = 0; i < rate.nrow(); i++){
    NumericVector r = rate( i, _);
    result[i] = ppexpC(q,r,t);
  }
  return Rcpp::wrap(result);
}
