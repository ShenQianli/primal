#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h> /* for NULL */

/* Declarations for the .C routines */
extern void Dantzig_api(int *, int *, double *, double *,
                        int *, double *, int *, double *,
                        double *, double *);
extern void SparseSVM_api(int *, int *, double *, double *,
                          int *, double *, int *, double *,
                          double *, double *, double *);
extern void QuantileRegression_api(int *, int *, double *, double *,
                                   double *, int *, double *, int *,
                                   double *, double *, double *);
extern void CompressedSensing_api(int *, int *, double *, double *,
                                  int *, double *, int *, double *,
                                  double *, double *);

static const R_CMethodDef CEntries[] = {
    {"Dantzig_api",              (DL_FUNC) &Dantzig_api,              10},
    {"SparseSVM_api",            (DL_FUNC) &SparseSVM_api,            11},
    {"QuantileRegression_api",   (DL_FUNC) &QuantileRegression_api,   11},
    {"CompressedSensing_api",    (DL_FUNC) &CompressedSensing_api,    10},
    {NULL, NULL, 0}
};

void R_init_PRIMAL(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
