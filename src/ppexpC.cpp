#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;


// declare
extern "C" double pexpC(double x, double scale, int lower_tail, int log_p);
LogicalVector equalgreaterouterC(int a, NumericVector b);
NumericVector rowSumsC(NumericMatrix x);

// [[Rcpp::export]]
NumericVector ppexpC(int q, NumericVector rate, NumericVector t){
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
  //int mi = min(t.size(), max(ind));
  int mi = t.size();
  std::vector<double> dt;
  if (t.size() > 1) {
    //std::vector<double> ts = Rcpp::as<std::vector<double> >(t);
    //std::vector<double> t0 = ts;
    //std::vector<double> tmi = ts;
    NumericVector t0, tmi = t;
    t0.erase(0);
    tmi.erase(mi-1);

    for(int i = 0; i < t.size(); i++){
      dt[i] = t0[i] - tmi[mi];
    }

    //std::vector<double> pe = pexpC(dt, rate[-mi-1]);
    //NumericVector cps = cumprod(1 - pe);
    //NumericVector cp = cp0.insert(0, 1);

  }

  //return Rcpp::wrap< std::vector<double> >( dt ) ;
  //NumericVector IntegerVector = as<NumericVector>( Rcpp::wrap(dt) ) ;
  return Rcpp::wrap(dt);

}

