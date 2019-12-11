#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include "flrtest.h"


double * estmodel(struct callinfo *model, double logrho, int retbeta){

  double * res;

  double * npXtXplusA, *npXtXplusA1Xt, *npXtXy;
  int i, j, ll, ur, cnt, add, info;
  int dim = *model->dim;
  int selector = model->selector;
  int n= *model->n;
  int p = *model->p;   // p shall be dim -1 if intercept is included in model
  double exprho = exp(logrho);
  double edf = 0, rss=0, yi = 0;

  double alpha = 1/ ((double) *model->n);
  double beta = 0;


  npXtXplusA = (double *) Calloc( ( selector * selector), double);
  npXtXplusA1Xt = (double *) Calloc( (selector * (*model->n)), double);
  npXtXy = (double *) Calloc( selector, double);

  if (retbeta){
    res = (double *) Calloc(selector, double);
  } else {
    res = (double *) Calloc(2, double);
  }

  // select rows and columns according to selector
  if (selector > dim){

    error("selector can be at max as large as dim!");

  } else if (selector == dim){

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

    // compute (1/np XtX + rho A)^(-1) * Xt
    F77_CALL(dgemm)("n","t",&selector, &n, &selector, &alpha, npXtXplusA, &selector,
		    model->X, &n, &beta, npXtXplusA1Xt, &selector);

  }

  if (retbeta){
    // compute beta
    for (i = 0; i < selector; i++){
      res[i] = 0.0;
      for (j = 0; j< n; j++){
        res[i] += npXtXplusA1Xt[j*selector + i] * model->y[j];
      }
    }

    // Free memory
    Free(npXtXplusA1Xt);
    Free(npXtXplusA);
    Free(npXtXy);

    return res;
  }

  // compute trace of X * npXtXplusA1Xt
  for (i=0; i<n; i++){
    for (j=0; j<selector; j++){
      edf += model->X[j*n + i] * npXtXplusA1Xt[i * selector + j];
    }
  }

  /* compute RSS of model */
  // 1. XtX %*% y
  for (i = 0; i < selector; i++){
    npXtXy[i] = 0;
    for (j=0; j < n; j++){
      npXtXy[i] +=( npXtXplusA1Xt[selector * j + i] * ( model->y[j] ));
    }

  }

  // 2. aggregate (X %*% XtX %*%y /p - y)^2
  for (i = 0; i<n; i++){
    yi = 0;
    for (j = 0; j<selector; j++){
      yi += model->X[ j * n + i ] * npXtXy[j];
    }
    rss += pow( yi/ ((double) p) - model->y[i], 2);
  }

  // Free memory
  Free(npXtXplusA1Xt);
  Free(npXtXplusA);
  Free(npXtXy);

  // sum( (X[,selector] %*% XtX1Xt %*% y * 1/p - y)^2 )
  res[1] = rss;
  res[0] = edf / ((double) p);
  return res;
}




/* R wrapper for above function */

SEXP R_estmodel(SEXP npXtX, SEXP X, SEXP y, SEXP Amat, SEXP df, SEXP n, SEXP p, SEXP dim, SEXP selector, SEXP logrho, SEXP retbeta){


  int lres;

  if (*INTEGER(retbeta)){
    lres = *INTEGER(selector);
  } else {
    lres = 2;
  }

  SEXP res = PROTECT(allocVector(REALSXP, lres));

  double * result;

  struct callinfo model;

  model.npXtX = REAL(npXtX);
  model.X = REAL(X);
  model.y = REAL(y);
  model.Amat = REAL(Amat);
  model.df = REAL(df);
  model.n = INTEGER(n);
  model.p = INTEGER(p);
  model.dim = INTEGER(dim);
  model.selector = *INTEGER(selector);


  result = estmodel(&model, *REAL(logrho), *INTEGER(retbeta));


  for (int i=0; i<lres; i++){
    REAL(res)[i] = result[i];
  }

  UNPROTECT(1);
  return res;
}
