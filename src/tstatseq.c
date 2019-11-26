// #include <float.h>
// #include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include "flrtest.h"



static double dfgrho(double x, struct callinfo *info){
  double res;
  res = dfGivenRho(x, info->npXtX, info->X,
		   info->Amat, *info->n, *info->dim, *info->p);
  res -= *info->df;
  return res;
}


SEXP tstatseq(SEXP y, SEXP X,  SEXP Amats, SEXP p, SEXP n,  SEXP df,
	      SEXP npXtX, SEXP tol, SEXP maxit,SEXP intercept){
  
  SEXP res;
  int np = INTEGER(p)[0];
  res = PROTECT(allocMatrix(REALSXP, np-1, 5));

  int intercpt = INTEGER(intercept)[0];
  int dim = np + intercpt;
  
  double logrho;
  //double *Amat = REAL(Amats);
  double *fullmodel, *nullmodel;

  const double tolerance = *REAL(tol);
  const int maxiter = *INTEGER(maxit);
  int Maxiter = maxiter;
  double Tolerance = tolerance;
  
  double *Tol = &Tolerance;
  int *Maxit = &Maxiter;
  
  struct callinfo info;
  
  info.npXtX = REAL(npXtX);
  info.X = REAL(X);
  info.y = REAL(y);
  info.Amat = REAL(Amats);
  info.df = REAL(df);
  info.n = INTEGER(n);
  info.p = INTEGER(p);
  info.dim = &dim;

   
  for (int j=0; j< (np-1); j++){
    
    // obtain rho s.t. df in splitted model equal to df in original model (df passed to funciton)
    //Rprintf("effdf: %f\n", *info.df);
    logrho = R_zeroin2(-200.0, 500.0,
		       dfgrho(-200, &info),
		       dfgrho(500, &info), (double (*)(double, void*)) dfgrho ,
		       (void *) &info,
		       Tol, Maxit);

    // reset Maxiter and Tol after every iteration
    if (*Maxit<0){
      error("Numerical root-finding not successful!");
    } else {
      *Maxit = maxiter;
      *Tol = tolerance;
    }

    // Rprintf("Rho: %f\n", logrho);
    
    // Rprintf("Amat[10,10], [%i]: %f\n", j, info.Amat[np*9 + 9/2]);
    
    // estimate full model
    info.selector = dim;
    fullmodel = estmodel(&info, logrho);

    Rprintf("[%i] df1: %f, df2: %f, lrho: %f\n",j, dfgrho(logrho, &info)+*info.df,fullmodel[0], logrho);
    
    REAL(res)[j] = fullmodel[0];
    REAL(res)[j + (np-1)] = fullmodel[1];
    
    // estimate null model
    if (intercpt == 1){
      info.selector = j+1;
    } else {
      info.selector = j;
    }
    
    
    nullmodel = estmodel(&info, logrho);

    REAL(res)[j +2*(np-1)] = nullmodel[0];
    REAL(res)[j +3*(np-1)] = nullmodel[1];

    REAL(res)[j +4*(np-1)] = logrho;

    //Rprintf("Rho: %f, effDfFull: %f, effDfNull: %f\n", logrho,
    //	    REAL(res)[j ],
    //	    REAL(res)[j +2*(np-1)]);
    //Rprintf("df-null: %f\n", nullmodel[0]);
    
    // store results
    
    // increment
    if (j< (np-2)){
      info.Amat += dim * dim;
    }
  }
  
  //info.Amat = REAL(Amats);
  
  
  
  /*
    edf = dfGivenRho(0, REAL(npXtX), REAL(X),
    Amat, *INTEGER(n), *INTEGER(p));
  Rprintf("rho: %f\n", logrho);
  Rprintf("f(a): %f\n", dfgrho(-100, &info));
  Rprintf("f(b): %f\n", dfgrho(100, &info));
  
  Rprintf("f(rho*): %f\n", dfGivenRho(logrho, (&info)->npXtX, (&info)->X,
  (&info)->Amat, *(&info)->n, *(&info)->dim, *(&info)->p));
  */

  
  /*a = PROTECT(allocVector(REALSXP, 1));
    REAL(a)[0] = 1.0;
    SET_VECTOR_ELT(res, 0, a);
  */
  
  UNPROTECT(1);
  return res;
}




static const R_CallMethodDef CallEntries[] = {
    {"tstatseq", (DL_FUNC) &tstatseq, 10},
    {"R_dfGivenRho", (DL_FUNC) &R_dfGivenRho, 7},
    {"R_estmodel", (DL_FUNC) &R_estmodel, 10},
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
