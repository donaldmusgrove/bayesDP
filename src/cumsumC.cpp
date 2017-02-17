#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<double> cumsumC(std::vector<double> x){
  // initialize the result vector
  std::vector<double> res(x.size());
  std::partial_sum(x.begin(), x.end(), res.begin());
  return res;
}

