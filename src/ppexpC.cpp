#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <numeric>
#include <functional>


using namespace Rcpp;


// declare
extern "C" double pexpC(double x, double scale, int lower_tail, int log_p);
LogicalVector equalgreaterouterC(int a, NumericVector b);
NumericVector rowSumsC(NumericMatrix x);
std::vector<double> cumsumC(std::vector<double> x);

// [[Rcpp::export]]
NumericVector ppexpC(int q, NumericVector rate, NumericVector t){

  std::vector<double> dt;
  std::vector<double> tC = Rcpp::as<std::vector<double> >(t);
  std::vector<double> pe(tC.size());
  double result;

    //  q[q < 0] <- 0
  if(q < 0){
    q = 0;
  }

  //  ind <- rowSums(outer(q, t, ">="))
  LogicalVector eg = equalgreaterouterC(q,t);
  int ind = 0;
  for(int i = 0; i < t.size(); i++){
    if(eg[i]==TRUE){
      ind++;
    }
  }

  //  ret <- pexp(q - t[ind], rate[ind])
  /*  double ret = pexpC(q - tC[ind], rate[ind],1,0);  */
  double ret = 0;

  int tCs = tC.size();
  int mi = std::min(tCs, ind);

  if (tC.size() > 1) {
    //  dt <- t[-1] - t[-mi]
    std::vector<double> t0(tC.begin() + 1, tC.end());
    std::vector<double> tmi(tC.begin(), tC.end() - mi - 1);

    //  pe <- pexp(dt, rate[-mi])
    /*  std::vector<double> pe = pexpC(dt, rate[-mi-1]);  */
    std::vector<double> pe(49999);
    std::fill(pe.begin(), pe.end(), 42);

    //  cp <- c(1, cumprod(1 - pe))
    std::vector<double> cp(50000);
    for(int i = 0; i < cp.size(); i++){
      cp[i] = 1 - pe[i];
      //cp.push_back(1 - pe[i]);
    }
    cp.insert(cp.begin(), 1);

    //  ret <- c(0, cumsum(cp[-length(cp)] * pe))[ind] + ret * cp[ind]
    cp.pop_back();

    std::vector<double> multeach(cp.size());
    for(int i = 0; i < multeach.size(); i++){
      multeach[i] = cp[i] * pe[i];
    }

    std::vector<double> cps = cumsumC(multeach);
    std::vector<double> cpsi = cps;
    cpsi.insert(cpsi.begin(), 0);
    result = cpsi[ind] + ret * cp[ind];
    int result = 1;
  }

  return Rcpp::wrap(result);

}

