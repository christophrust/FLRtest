#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include "flrtest.h"





/* 

Function to compute effective degrees of freedom of the smspl 
fit given a smoothing paramter rho 

*/
double dfGivenRho(double x, double *npXtX, double *X, double *Amat, int Nobs, int dim, int p){
  
  
  
  int i,j;
  int ll, ur;
  
  double alpha = 1/((double) Nobs);
  double beta = 0, res = 0;
  
  int info;
  
  // initialize arrays
  double *npXtXplusA, *npXtXplusA1Xt;
  npXtXplusA = (double *) Calloc(pow(dim,2), double);
  npXtXplusA1Xt = (double *) Calloc( dim*Nobs, double);
  
  
  // sum of npXtX + rho * A
  for (i=0; i< (dim * dim); i++){
    npXtXplusA[i] = npXtX[i] + exp(x) * Amat[i];
  }
  
  // compute cholesky
  F77_CALL(dpotrf)("U",&dim, npXtXplusA, &dim, &info);
  
  // compute inverse from cholesky
  if (info){
    error("Error when using cholesky decomposition, Info = %i!", info);
  } else {
    
    F77_CALL(dpotri)("U",&dim, npXtXplusA, &dim, &info);
  }
  
  // compute product of inv(XtX + rho * A) * X'
  if (info){
    
    error("Error when inverting after cholesky!");
  } else {
    
    // make matrix symmetric
    for (j = 0; j<dim; j++){
      for (i = j+1; i<dim; i++){
	ll = j*dim + i;
	ur = i*dim + j;
	npXtXplusA[ll] = npXtXplusA[ur];
      }
    }
    
    F77_CALL(dgemm)("n","t",&dim, &Nobs, &dim, &alpha, npXtXplusA, &dim,
		    X, &Nobs, &beta, npXtXplusA1Xt, &dim); 
  }
  
  // compute trace of X * npXtXplusA1Xt
  for (int i=0; i<Nobs; i++){
    for (int j=0; j<dim; j++){
      res += X[j*Nobs + i] * npXtXplusA1Xt[i * dim + j];
    }
  }

  Free(npXtXplusA1Xt);
  Free(npXtXplusA);
    
  res = res/ ((double) p);
  
  return res;
}



/* R wrapper around the above function */
SEXP R_dfGivenRho(SEXP rho, SEXP npXtX, SEXP X, SEXP Amat, SEXP Nobs, SEXP dim, SEXP p){

  SEXP res = PROTECT(allocVector(REALSXP, 1));

  //double df = 0;
  
  
  REAL(res)[0] = dfGivenRho(REAL(rho)[0],
			    REAL(npXtX),
			    REAL(X),
			    REAL(Amat),
			    *INTEGER(Nobs),
			    *INTEGER(dim),
			    *INTEGER(p));
  
  //REAL(res)[0] = df;
  UNPROTECT(1);
  return res;
}
