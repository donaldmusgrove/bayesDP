#include <Rcpp.h>
using namespace Rcpp;

// declare
NumericVector ppexpC(NumericVector x);

// [[Rcpp::export]]
NumericVector ppexpCV(NumericVector x){
  return ppexpC(1);
}

