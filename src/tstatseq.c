// #include <float.h>
// #include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include "flrtest.h"



static double dfgrho(double x, struct callinfo_smspl *info){
  double res;
  res = dfGivenRho(x, info->npXtX, info->X,
		   info->Amat, *info->n, *info->dim, *info->p);
  res -= (*info->df + 3);
  return res;
}

/*
  tstatseq computes the test sequence described in the manuscript
  "directed local tesing in functional regression".

  Inputs:
   - y: numeric vector of length n containing observations of dependent variable
   - X: numeric matrix n times (p + intercept) containing discretized functions
        (every row is one observation)
   - Amats: a stacked (as.vector) numeric array (p times (p + intercept)
        times (p+ intercept))
        containing the smoother matrices of smooting spline(see manuscript).
   - p: integer, number of discretization points.
   - n: number of observations.
   - df: number of degrees of freedom for both models
   - npXtX: matrix used in the model fitting, is available in the return of EstFLM()
   - tol: tolerance for the root finding procedure
   - maxit: maximum number of iterations of the root finding procedure

  Output:
   - A numeric matrix with p rows and 6 columns:
     1. rss of the full model
     2. edf of the full model
     3. rss of null model
     4. edf of null model
     5. logarithm of smoothing parameter rho
     6. test statistic (see manuscript)
 */
SEXP tstatseq_smspl(SEXP y, SEXP X,  SEXP Amats, SEXP p, SEXP n,  SEXP df,
	      SEXP npXtX, SEXP tol, SEXP maxit,SEXP intercept){

  // obtain input variables and initialize pointers to these objects
  int np = INTEGER(p)[0];
  int intercpt = INTEGER(intercept)[0];
  int dim = np + intercpt;
  const double tolerance = *REAL(tol);
  const int maxiter = *INTEGER(maxit);
  int Maxiter = maxiter;  // Maxiter is changed by R_zeroin2 in every interation
  double Tolerance = tolerance; // same for Tolerance

  double *Tol = &Tolerance; // pointers passed to R_zeroin2
  int *Maxit = &Maxiter;


  // initialize return object
  SEXP res;
  res = PROTECT(allocMatrix(REALSXP, np-1, 6));

  // containers for intermediate results
  double logrho;
  double *fullmodel, *nullmodel;

  // structure info containing the neccessary information for the call to dfgrho
  struct callinfo_smspl info;

  info.npXtX = REAL(npXtX);
  info.X = REAL(X);
  info.y = REAL(y);
  info.Amat = REAL(Amats);
  info.df = REAL(df);
  info.n = INTEGER(n);
  info.p = INTEGER(p);
  info.dim = &dim;

  // main iteration
  for (int j=0; j< (np-1); j++){

    // obtain rho s.t. df in splitted model equal to df in original
    // model (df passed to function dfgrho)
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
    fullmodel = estmodel_smspl(&info, logrho,0);

    // Rprintf("[%i] df1: %f, df2: %f, lrho: %f\n",j, dfgrho(logrho, &info)+*info.df,fullmodel[0], logrho);

    REAL(res)[j] = fullmodel[0];          // 1st col: rss of full model
    REAL(res)[j + (np-1)] = fullmodel[1]; // 2nd col: edf of full model

    // estimate null model
    if (intercpt == 1){
      info.selector = j+2;
    } else {
      info.selector = j+1;
    }


    nullmodel = estmodel_smspl(&info, logrho,0);

    REAL(res)[j +2*(np-1)] = nullmodel[0]; // 3rd col: rss of null model
    REAL(res)[j +3*(np-1)] = nullmodel[1]; // 4th col: edf of null model

    REAL(res)[j +4*(np-1)] = logrho;       // 5th col: log of rho

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
    {"tstatseq_smspl", (DL_FUNC) &tstatseq_smspl, 10},
    {"R_dfGivenRho", (DL_FUNC) &R_dfGivenRho, 7},
    {"R_estmodel_smspl", (DL_FUNC) &R_estmodel_smspl, 11},
    {"R_spline_basis",(DL_FUNC) &R_spline_basis, 4},
    {"R_SplitSplineBasis",(DL_FUNC) &R_SplitSplineBasis, 3},
    {"R_estmodel_spl",(DL_FUNC) &R_estmodel_spl, 6},
    {NULL, NULL, 0}
};



void R_init_FLRtest(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

