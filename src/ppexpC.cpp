#include <Rcpp.h>
#include <algorithm>
#include <Rcpp/sugar/functions/cumprod.h>
#include <Rcpp/sugar/functions/cumsum.h>

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

  //double ret = pexpC(q - t[ind], rate[ind],1,0);
  //int mi = min(t.size(), max(ind));
  int mi = t.size();
  std::vector<double> dt;
  std::vector<double> pe;
  std::vector<double> cpi;
  std::vector<double> cp;
  std::vector<double> tC = Rcpp::as<std::vector<double> >(t);
  double ret;

  if (tC.size() > 1) {
    std::vector<double> tC = Rcpp::as<std::vector<double> >(t);
    //std::vector<double> t0 = tC;
    //std::vector<double> tmi = tC;

    std::vector<double> t0, tmi = tC;
    //t0.erase(t0.begin());
    //tmi.erase(mi-1);

    //for(int i = 0; i < t.size(); i++){
    //  dt[i] = t0[i] - tmi[mi];
    //}

    std::fill(pe.begin(), pe.end(), 42);

    //std::vector<double> pe = pexpC(dt, rate[-mi-1]);

    //for(int i = 0; i < tC.size(); i++){
    //  cpi[i] = 1 - pe[i];
    //}

    //NumericVector cps = cumprod(1 - pe);

    std::vector<double> cp;
    cp.insert(cp.begin(), 0);
    //ret <- c(0, cumsum(cp[-length(cp)] * pe))[ind] + ret * cp[ind]



    //NumericVector cp = cp0.insert(0, 1);

  }

  //return Rcpp::wrap< std::vector<double> >( dt ) ;
  //NumericVector IntegerVector = as<NumericVector>( Rcpp::wrap(dt) ) ;
  return Rcpp::wrap(tC);

}

