#include "flrtest.h"

/*
  estmodel computes RSS, effective degrees of freedom or coefficient function
  for a smoothing spline estimation of a functional linear model. Via a selector,
  some parts of the coefficient function can be restricted to a specified value.

  Inputs:
   - model: a pointer to a callinfo_spl stuct (defined in flrtest.h) containing all the
       relevant data and model parameters for the fit.
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

double * estmodel_spl(struct callinfo_spl *model, int retbeta){

  /* array containing result of estimation */
  double *res;

  /* allocate space for necessary intermediate results */
  double *pXB;   // 1/p * X * Basis
  pXB = (double *) R_alloc( *model->n * *model->dim, sizeof(double));

  /* this in fact will hold pXB'pXB */
  double *XtX;
  XtX = (double *) R_alloc( *model->dim * *model->dim, sizeof(double));

  double *XtX1Xt;   // 1/p * X * Basis
  XtX1Xt = (double *) R_alloc( *model->n * *model->dim, sizeof(double));

  /* compute pXB  (1/p * X %*% B)*/
  //double alpha = 1.0;
  double alpha = 1/(double) *model->p;
  double beta = 0.0;

  double temp = 0.0;


  F77_CALL(dgemm)("n", "n", model->n, model->dim, model->p, &alpha,
                  model->X, model->n, model->Basis, model->p, &beta,
                  pXB, model->n);

  /* compute XtX*/
  internal_crossprod(pXB, *model->n, *model->dim, pXB, *model->n, *model->dim, XtX);

  /* invert XtX using cholesky decomposition */
  int info, ll, ur;

  F77_CALL(dpotrf)("U", model->dim, XtX, model->dim, &info);

  if (info){
    error("Error ocurred during cholesky decomposition, Info = %i!", info);
  } else {

    F77_CALL(dpotri)("U", model->dim, XtX, model->dim, &info);
  }

  if (info){

    error("Error ocurred while inverting after cholesky decomposition!");
  } else {

    // make matrix symmetric
    for (int j = 0; j<*model->dim; j++){
      for (int i = j+1; i< *model->dim; i++){
        ll = j**model->dim + i;
        ur = i**model->dim + j;
        XtX[ll] = XtX[ur];
      }
    }
  }

  // compute XtX^(-1) * Xt
  alpha = 1.0;
  F77_CALL(dgemm)("n","t", model->dim, model->n,
                  model->dim, &alpha, XtX, model->dim,
                  pXB, model->n, &beta, XtX1Xt, model->dim);


  /* compute splines coefficient: theta = XtX^(-1) %*% Xt %*% y */
  double *theta;
  theta = (double *) R_alloc(*model->dim, sizeof(double));

  for (int i = 0; i < *model->dim; i++){
    temp = 0.0;
    for (int j = 0; j < *model->n; j++){
      temp += XtX1Xt[i + j * *model->dim] * model->y[j];
    }

    theta[i] = temp;
    // Rprintf("%f;;", temp);
  }

  if (retbeta){

    /* compute functional coefficient: beta = basis %*% theta */
    res = (double *) R_alloc(*model->p, sizeof(double));

    for (int i = 0; i < *model->p; i++){
      temp = 0.0;
      for (int j = 0; j < *model->dim; j++){
        temp += model->Basis[i + j * *model->p] * theta[j];
      }

      res[i] = temp;
    }

    return res;
  }

  /* return only rss and df */
  res = (double *) R_alloc(2, sizeof(double));

  /* compute RSS */
  double rss = 0.0;

  for (int i = 0; i < *model->n; i++){
    temp = 0.0;
    for (int j = 0; j < *model->dim; j++){
      temp += pXB[i + j* *model->n] * theta[j];
    }
    rss += pow(model->y[i] - temp, 2);
  }
  res[1] = rss;
  res[0] = *model->dim;

  for (int i = 0; i<2; i++){
    Rprintf("%f\n; ", res[i]);
  }

  for (int i=0; i < 10; i++) {
    Rprintf("%f; ", XtX1Xt[i]);
  }
  return res;
}




SEXP R_estmodel_spl(SEXP y, SEXP X, SEXP basis, SEXP n, SEXP p, SEXP dim, SEXP retbeta) {

  struct callinfo_spl model;

  model.X = REAL(X);
  model.y = REAL(y);
  model.Basis = REAL(basis);
  model.n = INTEGER(n);
  model.p = INTEGER(p);
  model.dim = INTEGER(dim);

  double *val;
  val = estmodel_spl(&model, *INTEGER(retbeta));

  /* copy result into return SEXP */
  SEXP res;
  if (*INTEGER(retbeta)){
    res = PROTECT(allocVector(REALSXP, *model.p));
    for (int i = 0; i < *model.p; i++){
      REAL(res)[i] = val[i];
    }
    UNPROTECT(1);
    return res;
  }

  res = PROTECT(allocVector(REALSXP, 2));

  REAL(res)[0] = val[0];
  REAL(res)[1] = val[1];

  UNPROTECT(1);
  return res;
}
