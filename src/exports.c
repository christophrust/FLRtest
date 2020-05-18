#include "flrtest.h"

static const R_CallMethodDef CallEntries[] = {
        {"tstatseq_smspl", (DL_FUNC) &tstatseq_smspl, 10},
        {"R_dfGivenRho", (DL_FUNC) &R_dfGivenRho, 7},
        {"R_estmodel_smspl", (DL_FUNC) &R_estmodel_smspl, 11},
        {"R_spline_basis",(DL_FUNC) &R_spline_basis, 4},
        {"R_SplitSplineBasis",(DL_FUNC) &R_SplitSplineBasis, 3},
        {"R_estmodel_spl",(DL_FUNC) &R_estmodel_spl, 9},
        {NULL, NULL, 0}
};



void R_init_FLRtest(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

