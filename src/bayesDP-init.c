
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

extern SEXP bayesDP_ppexpM(SEXP, SEXP, SEXP);
extern SEXP bayesDP_ppexpV(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"bayesDP_ppexpM", (DL_FUNC) &bayesDP_ppexpM, 3},
  {"bayesDP_ppexpV", (DL_FUNC) &bayesDP_ppexpV, 3},
  {NULL, NULL, 0}
};

void R_init_bayesDP(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
