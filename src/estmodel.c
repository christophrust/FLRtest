#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include "flrtest.h"


double * estmodel(struct callinfo *model, double logrho){
  

  static double res[2];
  double * npXtXplusA, *npXtXplusA1Xt;
  int i, j, ll, ur, cnt, add, info;
  int dim = *model->p;
  int selector = model->selector;
  int n= *model->n;
  double exprho = exp(logrho);
  double edf = 0;
  
  double alpha = 1/ ((double) *model->n);
  double beta = 0;
    
  npXtXplusA = (double *) malloc( ( selector * selector) * sizeof(double));
  npXtXplusA1Xt = (double *) malloc( (selector * (*model->n)) * sizeof(double));

  if (selector > dim){

    error("selector can be at max as large as dim!");

  } else if (selector >= dim){

    for (i=0; i< (dim * dim); i++){
      npXtXplusA[i] = model->npXtX[i] + exprho * model->Amat[i];
    }
    
  } else {
    
    cnt = 1;
    add = 0;
    
    for (i=0; i< (selector * selector); i++, cnt++){

      // Rprintf("cnt: %i; selector: %i; small: %i; large: %i\n",cnt, selector, i, i+add);
      npXtXplusA[i] = model->npXtX[i + add] +  exprho * model->Amat[i+add];

      // jump to next column
      if (cnt == selector){
	add += (dim - selector);
	cnt = 0;
      } 
      
    }
  }
  
  // compute cholesky
  F77_CALL(dpotrf)("U",&selector, npXtXplusA, &selector, &info);
  
  if (info){
    error("Error when using cholesky decomposition, Info = %i!", info);
  } else {
    
    F77_CALL(dpotri)("U",&selector, npXtXplusA, &selector, &info);
  }

  if (info){
    
    error("Error when inverting after cholesky!");
  } else {
    
    // make matrix symmetric
    for (j = 0; j<selector; j++){
      for (i = j+1; i<selector; i++){
	ll = j*selector + i;
	ur = i*selector + j;
	npXtXplusA[ll] = npXtXplusA[ur];
      }
    }
      
    
    F77_CALL(dgemm)("n","t",&selector, &n, &selector, &alpha, npXtXplusA, &selector,
		    model->X, &n, &beta, npXtXplusA1Xt, &selector); 

  
  }
  

  
  // compute trace of X * npXtXplusA1Xt
  for (i=0; i<n; i++){
    for (j=0; j<selector; j++){
      edf += model->X[j*n + i] * npXtXplusA1Xt[i * selector + j];
    }
  }

  /* compute RSS of model */
  
  // sum( (X[,selector] %*% XtX1Xt %*% y * 1/p - y)^2 )
  res[0] = edf / ((double) dim);
  return res;
}




/* R wrapper for above function */

SEXP R_estmodel(SEXP npXtX, SEXP X, SEXP Amat, SEXP df, SEXP n, SEXP p, SEXP selector, SEXP logrho){

  SEXP res = PROTECT(allocVector(REALSXP, 2));
  double * result;
  
  struct callinfo model;

  model.npXtX = REAL(npXtX);
  model.X = REAL(X);
  model.Amat = REAL(Amat);
  model.df = REAL(df);
  model.n = INTEGER(n);
  model.p = INTEGER(p);
  model.selector = *INTEGER(selector);
    
  result = estmodel(&model, *REAL(logrho));

  REAL(res)[0] = result[0];
  REAL(res)[1] = result[1];
  UNPROTECT(1);
  return res;
}
