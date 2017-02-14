#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector cumsum2(NumericVector x){
  // initialize the result vector
  NumericVector res(x.size());
  std::partial_sum(x.begin(), x.end(), res.begin());
  return res;
}

