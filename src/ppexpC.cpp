#include <Rcpp.h>
using namespace Rcpp;

// declare
double pexpC(double x, double scale, int lower_tail, int log_p);

// [[Rcpp::export]]
NumericVector ppexpC(NumericVector x){
  //out = cumsumC(x)
  return 1;
}

