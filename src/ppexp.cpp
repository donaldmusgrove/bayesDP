#include <Rcpp.h>
using namespace Rcpp;
// Start

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]]
double ppexp(double q, const vec& x, const vec& cuts) {

  int nc = cuts.n_elem;
  double ret;
  double mi = nc;
  vec cp1 = ones<vec>(1);
  vec cp0 = zeros<vec>(1);

  // Find cut point where q >= cuts
  double ind;
  vec indcut(nc);

  for(int i = 0; i < nc; i++){
    indcut(i) = q >= cuts(i);
  }
  ind = sum(indcut) -1;

  // Compute cdf at cuts[q]
  ret = R::pexp(q - cuts(ind), 1/x(ind), TRUE, FALSE);

  if( mi > ind){
    mi = ind;
  }

  if (nc > 1) {
    vec dt = cuts(find(cuts>0)) - cuts(find(cuts != cuts(mi)));

    vec xc = x(find(cuts != cuts(mi)));

    int nc1 = nc-1;
    vec pe(nc1);
    for(int i = 0; i < nc1; i++){
      pe(i) = R::pexp(dt(i), 1/xc(i), TRUE, FALSE);;
    }
    vec cp = cumprod(1 - pe);
    cp = join_vert(cp1,cp);

    //ret <- c(0, cumsum(cp[-length(cp)] * pe))[ind] + ret*cp[ind]
    int ncp = cp.n_elem;
    vec ret0 = cumsum(cp.head(ncp-1)%pe);
    ret0 = join_vert(cp0, ret0);
    double ret1 = ret0(ind) + ret*cp(ind);
    return(ret1);
  } else{
    return(ret);
  }
}


// [[Rcpp::export]]
vec ppexpM(double q, const mat& x, const vec& cuts) {

  mat y = x.t();

  int ny = y.n_cols;
  vec ret(ny);
  for(int i = 0; i < ny; i++){
    ret(i) = ppexp(q, y.col(i), cuts);
  }
  return(ret);
}

// End
