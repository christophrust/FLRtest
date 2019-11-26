// #include <float.h>
// #include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include "flrtest.h"



static double dfgrho(double x, struct callinfo *info){
  double res;
  res = dfGivenRho(x, info->npXtX, info->X,
		   info->Amat, *info->n, *info->p);
  res -= *info->df;
  return res;
}


SEXP tstatseq(SEXP y, SEXP X,  SEXP Amats, SEXP p, SEXP n,  SEXP df,
	      SEXP npXtX, SEXP tol, SEXP maxit,SEXP intercept){
  
  SEXP res;
  int np = INTEGER(p)[0];
  res = PROTECT(allocMatrix(REALSXP, np-1, 4));

  int intercpt = INTEGER(intercept)[0];
  
  double rho;
  //double *Amat = REAL(Amats);
  double *fullmodel, *nullmodel;
  
  double *Tol = REAL(tol);
  int *Maxit = INTEGER(maxit);
  
  struct callinfo info;
  
  info.npXtX = REAL(npXtX);
  info.X = REAL(X);
  info.Amat = REAL(Amats);
  info.df = REAL(df);
  info.n = INTEGER(n);
  info.p = INTEGER(p);

  rho = R_zeroin2(-200.0, 500.0,
		    dfgrho(-100, &info),
		    dfgrho(100, &info), (double (*)(double, void*)) dfgrho ,
		    (void *) &info,
		    Tol, Maxit);
  Rprintf("Rho1: %f\n", rho);

  
  for (int j=0; j< (np-1); j++){
    
    // obtain rho s.t. df in splitted model equal to df in original model (df passed to funciton)
    Rprintf("effdf: %f\n", *info.df);
    rho = R_zeroin2(-100.0, 100.0,
		    dfgrho(-100, &info),
		    dfgrho(100, &info), (double (*)(double, void*)) dfgrho ,
		    (void *) &info,
		    Tol, Maxit);

    
    Rprintf("Rho: %f\n", rho);
    // estimate full model
    

    info.selector = *INTEGER(p);
    fullmodel = estmodel(&info, rho);
    Rprintf("df-full: %f\n", fullmodel[0]);

    // estimate null model
    if (intercpt == 1){
      info.selector = j+1;
    } else {
      info.selector = j;
    }
    
    
    nullmodel = estmodel(&info, rho);
    Rprintf("df-null: %f\n", nullmodel[0]);
    
    // store results
    
    // increment
    if (j< (np-2)){
      info.Amat += np * np;
    }
  }
  
  //info.Amat = REAL(Amats);
  
  
  
  /*
    edf = dfGivenRho(0, REAL(npXtX), REAL(X),
    Amat, *INTEGER(n), *INTEGER(p));
  */
  Rprintf("rho: %f\n", rho);
  Rprintf("f(a): %f\n", dfgrho(-100, &info));
  Rprintf("f(b): %f\n", dfgrho(100, &info));
  
  Rprintf("f(rho*): %f\n", dfGivenRho(rho, (&info)->npXtX, (&info)->X,
				      (&info)->Amat, *(&info)->n, *(&info)->p));
  
  /*a = PROTECT(allocVector(REALSXP, 1));
    REAL(a)[0] = 1.0;
    SET_VECTOR_ELT(res, 0, a);
  */
  
  UNPROTECT(1);
  return res;
}




static const R_CallMethodDef CallEntries[] = {
    {"tstatseq", (DL_FUNC) &tstatseq, 10},
    {"R_dfGivenRho", (DL_FUNC) &R_dfGivenRho, 6},
    {"R_estmodel", (DL_FUNC) &R_estmodel, 8},
    {NULL, NULL, 0}
};



void R_init_FLRtest(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}



/*
TODO: test dfGivenRho

implement loop over different splits
*/
