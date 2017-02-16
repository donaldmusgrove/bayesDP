#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;


// declare
extern "C" double pexpC(double x, double scale, int lower_tail, int log_p);
LogicalVector equalgreaterouterC(int a, NumericVector b);
NumericVector rowSumsC(NumericMatrix x);

// [[Rcpp::export]]
int ppexpC(int q, NumericVector rate, NumericVector t){
  if(q < 0){
    q = 0;
  }

  LogicalVector eg = equalgreaterouterC(q,t);

  int ind = 0;
  for(int i = 0; i < t.size(); i++){
    if(eg[i]==TRUE){
      ind++;
    }
  }

  double ret = pexpC(q - t[ind], rate[ind],1,0);
  //int mi = min(t.size(), algorithm::max(ind));
/*
  if (t.size() > 1) {

    NumericVector ts = t.erase(0)
    NumericVector tmi = t.erase(mi-1);

    for(int i = 0; i < ts.size(); i++){
      NumericVector dt[i] = ts[i] - tmi[i];
    }

    NumericVector pe = pexpC(dt, rate[-mi-1]);
    NumericVector cps = cumprod(1 - pe);
    NumericVector cp = cp0.insert(0, 1);
    ret1
    ret2

    double ret =
  }
*/

  return ret;
}

