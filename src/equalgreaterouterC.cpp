#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
LogicalVector equalgreaterouterC(int a, NumericVector b) {
  int n = b.size();
  LogicalVector out(n);
  for (int i = 0; i < n; i++) {
    if(a >= b[i]) {
      out[i] = TRUE;
    }
    else{
      out[i] = FALSE;
    }
  }
  return out;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
equalgreaterouterC(42,c(1,4,77,65,32))
*/
