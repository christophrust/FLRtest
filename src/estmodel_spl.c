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
  res = (double *) R_alloc(1, sizeof(double));
  res[0] = 2.0;

  /* allocate space for necessary intermediate results */
  double *pXB;   // 1/p * X * Basis
  pXB = (double *) R_alloc( *model->n * *model->dim, sizeof(double));

  /* this in fact will hold pXB'pXB */
  double *XtX;
  XtX = (double *) R_alloc( *model->dim * *model->dim, sizeof(double));

  /* compute pXB  (1/p * X %*% B)*/
  //double alpha = 1.0;
  double alpha = 1/(double) *model->p;
  double beta = 0.0;

  Rprintf("%i;--\n", (*model->n * *model->dim));

  F77_CALL(dgemm)("n", "n", model->n, model->dim, model->p, &alpha,
                  model->X, model->n, model->Basis, model->p, &beta,
                  pXB, model->n);

  /* compute matrix product by hand: */
  /* double temp = 0.0; */
  /* for (int i=0; i < *model->n; i++){ */
  /*   for (int j = 0; j < *model->dim; j++){ */
  /*     temp = 0.0; */
  /*     for (int k =0; k < *model->p; k++){ */
  /*       temp += model->X[i + k* *model->n] * model->Basis[j * *model->p + k]; */
  /*     } */
  /*     pXB1[i + j * *model->dim ] = temp; */
  /*   } */
  /* } */

  /* for (int i=0; i < 10; i++){ */
  /*   Rprintf("%f, ",model->Basis[i]); */

  /* } */
  /* Rprintf("\n--------\n"); */


  /* for (int i=0; i < 10; i++){ */
  /*   Rprintf("%f, ",model->X[i]); */
  /* } */
  /* Rprintf("\n--------\n"); */


  /* for (int i=0; i < 10; i++){ */
  /*   Rprintf("%f, ",pXB[i]); */
  /* } */
  /* Rprintf("\n-----a---\n"); */


  /* compute XtX*/
  internal_crossprod(pXB, *model->n, *model->dim, pXB, *model->n, *model->dim, XtX);

  /* for (int i=0; i < 10; i++) { */
  /*   Rprintf("%f; ", XtX[i]); */
  /* } */


  return res;
}




SEXP R_estmodel_spl(SEXP y, SEXP X, SEXP basis, SEXP n, SEXP p, SEXP dim) {

  struct callinfo_spl model;

  model.X = REAL(X);
  model.y = REAL(y);
  model.Basis = REAL(basis);
  model.n = INTEGER(n);
  model.p = INTEGER(p);
  model.dim = INTEGER(dim);

  double *val;
  val = estmodel_spl(&model, 0);

  SEXP res;
  res = PROTECT(allocVector(REALSXP, 1));
  UNPROTECT(1);

  REAL(res)[0] = val[0];

  return res;
}
