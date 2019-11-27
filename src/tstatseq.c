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
  res -= (*info->df + 3);
  return res;
}


SEXP tstatseq(SEXP y, SEXP X,  SEXP Amats, SEXP p, SEXP n,  SEXP df,
	      SEXP npXtX, SEXP tol, SEXP maxit,SEXP intercept){
  
  SEXP res;
  int np = INTEGER(p)[0];
  res = PROTECT(allocMatrix(REALSXP, np-1, 6));

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

  /*
  for (int j=0; j< 7; j++){
    
    // obtain rho s.t. df in splitted model equal to df in original model (df passed to funciton)
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
    Rprintf("logrho = %f\n", logrho);
    info.Amat += dim * dim;
  }

  error("stop");
  */
  
  for (int j=0; j< (np-1); j++){
    
    // obtain rho s.t. df in splitted model equal to df in original model (df passed to funciton)
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
    
    // estimate full model
    info.selector = dim;
    fullmodel = estmodel(&info, logrho);

    // Rprintf("[%i] df1: %f, df2: %f, lrho: %f\n",j, dfgrho(logrho, &info)+*info.df,fullmodel[0], logrho);
    
    REAL(res)[j] = fullmodel[0];
    REAL(res)[j + (np-1)] = fullmodel[1];
    
    // estimate null model
    if (intercpt == 1){
      info.selector = j+2;
    } else {
      info.selector = j+1;
    }
    
    
    nullmodel = estmodel(&info, logrho);

    REAL(res)[j +2*(np-1)] = nullmodel[0];
    REAL(res)[j +3*(np-1)] = nullmodel[1];

    REAL(res)[j +4*(np-1)] = logrho;

    // compute test statistic: ( (rss0 - rss1)/df1-df0) / (rss1 / (n - df1) 
    REAL(res)[j +5*(np-1)] =
      ( (REAL(res)[j +3*(np-1)] - REAL(res)[j +1*(np-1)]) /
	(REAL(res)[j +0*(np-1)] - REAL(res)[j +2*(np-1)])) /
      (REAL(res)[j +1*(np-1)] /
       ( (double) *INTEGER(n) - REAL(res)[j +2*(np-1)]));



    if (j< (np-2)){
      info.Amat += dim * dim;
    }
  }
  
  
  
  
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

