#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include "flrtest.h"

/*
  estmodel computes RSS, effective degrees of freedom or coefficient function
  for a smoothing spline estimation of a functional linear model. Via a selector,
  some parts of the coefficient function can be restricted to a specified value.

  Inputs:
   - model: a pointer to a callinfo_smspl stuct (defined in flrtest.h) containing all the
       relevant data and model parameters for the fit.
   - logrho: predefined smoothing parameter
   - retbeta: if zero, rss and edf are computed otherwise the coefficientfunction

  Output:
   A pointer to a double array. Either of length 2 (retbeta == 0), in which case
   rss and edf is returned, or of lenght model->selector (retbeta == 1), then beta
   itself is returned.


  Details:
   The entry model->selector is an integer which in the currect implementation
   marks the end on the right hand side of the unrestricted domain.

   The model estimate is obtained as follows after subsetting entries (1:selector) of
   rows and columns respectively of all matrices:

   - for beta:
     (npXtX + rho * Am)^{-1} Xt y
   - for rss:
     ((npXtX + rho * Am)^{-1} Xt y) - y
   - for edf:
     trace(X  (npXtX + rho * Am)^{-1} Xt)

   Here, every matrix multiplication is a call to the Fortran routine dgemm and
   the matrix inversion is obtained by using dpotrf (cholesky) and
   dpotri (inverse based on cholesky).
*/


double * estmodel_smspl(struct callinfo_smspl *model, double logrho, int retbeta){

  // initialize container objects for intermediate results
  double * npXtXplusA, *npXtXplusA1Xt, *npXtXy;

  // helper variables
  int i, j, ll, ur, cnt, add, info;

  // obtain the information from the callinfo
  int dim = *model->dim;
  int selector = model->selector;
  int n= *model->n;
  int p = *model->p;   // p shall be dim -1 if intercept is included in model
  double exprho = exp(logrho);
  double edf = 0, rss=0, yi = 0;

  // variables used in the fortran routines
  double alpha = 1/ ((double) *model->n);
  double beta = 0;

  // allocate memory for all arrays containing intermediate results
  npXtXplusA = (double *) Calloc( ( selector * selector), double);
  npXtXplusA1Xt = (double *) Calloc( (selector * (*model->n)), double);
  npXtXy = (double *) Calloc( selector, double);

  // initialize return object and allocate memory
  double * res;

  if (retbeta){
    res = (double *) Calloc(selector, double);
  } else {
    res = (double *) Calloc(2, double);
  }

  // compute npXtXplusA by using subsets of npXtX and Amat
  // select rows and columns according to selector
  if (selector > dim){

    error("selector can be at max as large as dim!");

  } else if (selector == dim){

    // then this is the full model
    for (i=0; i< (dim * dim); i++){
      npXtXplusA[i] = model->npXtX[i] + exprho * model->Amat[i];
    }

  } else {

    cnt = 1;  // we use cnt to know in wich row of the matrix we are
    add = 0;  // how many entries we skip whenever we jump from one column to the next

    for (i=0; i< (selector * selector); i++, cnt++){

      // Rprintf("cnt: %i; selector: %i; small: %i; large: %i\n",cnt, selector, i, i+add);
      npXtXplusA[i] = model->npXtX[i + add] +  exprho * model->Amat[i+add];

      // jump to next column and increase add by the number of jumped entries
      if (cnt == selector){
        add += (dim - selector);
        cnt = 0;
      }
    }
  }

  // compute cholesky
  F77_CALL(dpotrf)("U",&selector, npXtXplusA, &selector, &info);

  if (info){
    error("Error ocurred during cholesky decomposition, Info = %i!", info);
  } else {

    F77_CALL(dpotri)("U",&selector, npXtXplusA, &selector, &info);
  }

  if (info){

    error("Error ocurred while inverting after the cholesky decomposition!");
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

  // if edf and rss are to be returned
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

SEXP R_estmodel_smspl(SEXP npXtX, SEXP X, SEXP y, SEXP Amat, SEXP df, SEXP n, SEXP p, SEXP dim, SEXP selector, SEXP logrho, SEXP retbeta){


  int lres;

  if (*INTEGER(retbeta)){
    lres = *INTEGER(selector);
  } else {
    lres = 2;
  }

  SEXP res = PROTECT(allocVector(REALSXP, lres));

  double * result;

  struct callinfo_smspl model;

  model.npXtX = REAL(npXtX);
  model.X = REAL(X);
  model.y = REAL(y);
  model.Amat = REAL(Amat);
  model.df = REAL(df);
  model.n = INTEGER(n);
  model.p = INTEGER(p);
  model.dim = INTEGER(dim);
  model.selector = *INTEGER(selector);


  result = estmodel_smspl(&model, *REAL(logrho), *INTEGER(retbeta));


  for (int i=0; i<lres; i++){
    REAL(res)[i] = result[i];
  }

  UNPROTECT(1);
  return res;
}
