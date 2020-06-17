#include "flrtest.h"

/*
  estmodel computes RSS, effective degrees of freedom or coefficient function
  for a smoothing spline estimation of a functional linear model. Via a selector,
  some parts of the coefficient function can be restricted to a specified value.

  Inputs:
   - model: a pointer to a callinfo_spl stuct (defined in flrtest.h) containing all the
       relevant data and model parameters for the fit.
       The member 'selector' plays a crucial role in the fit as it determines,
       whether the full model (selector = dim) or whether the small model
       (in this case, selector must be equal to the selector member of the return of the
       SplitSplineBasis function)
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


  /* dimension of the full model */
  //int dim_full = *model->k + model->intercept;

  /* and dimension of potentially small model */
  int dim = model->selector + model->intercept;

  /* allocate space for necessary intermediate results */
  double *pXB = model->pXB;   // 1/p * X * Basis
  //pXB = (double *) R_alloc( *model->n * dim, sizeof(double));

  /* this in fact will hold pXB'pXB */
  double *XtX = model->XtX;
  // XtX = (double *) R_alloc( dim * dim, sizeof(double));

  double *XtX1Xt = model->XtX1Xt;   // 1/p * X * Basis
  // XtX1Xt = (double *) R_alloc( *model->n * dim, sizeof(double));

  /* compute pXB  (1/p * X %*% B)*/
  //double alpha = 1.0;
  double alpha = 1/(double) *model->p;
  double beta = 0.0;

  double temp = 0.0;

  /* include intercept if necessary */
  if (model->intercept == 1){
    for (int i=0; i < *model->n; i++) {
      pXB[i] = 1.0;
    }
  }


  /*
    DGEMM 'n', 'n' arguments:
    C = alpha * A %*% B + beta C
    rows A/rows C, cols B/cols C, cols A/rows B, alpha, A, rows A, B, cols B,
    beta, rows C
  */

  F77_CALL(dgemm)("n", "n", model->n, &model->selector, model->p, &alpha,
                  model->X, model->n, model->Basis, model->p, &beta,
                  pXB + model->intercept * *model->n, model->n);



  /* compute XtX*/
  /* computes C = A'B with
   A, rows A, cols A, B, rows B, cols B, C */
  internal_crossprod(pXB, *model->n, dim, pXB,
                     *model->n, dim, XtX);


  /* invert XtX using cholesky decomposition */
  int info, ll, ur;

  /*
    DPOTRF call, version 'U':
    rows A, A, rows A, info

   */
  F77_CALL(dpotrf)("U", &dim, XtX, &dim, &info);

  if (info){
    error("Error ocurred during cholesky decomposition, Info = %i!", info);
  } else {

    F77_CALL(dpotri)("U", &dim, XtX, &dim, &info);
  }

  if (info){

    error("Error ocurred while inverting after cholesky decomposition!");
  } else {

    // make matrix symmetric
    for (int j = 0; j < dim; j++){
      for (int i = j+1; i < dim; i++){
        ll = j * dim + i;
        ur = i * dim + j;
        XtX[ll] = XtX[ur];
      }
    }
  }


  // compute XtX^(-1) * Xt
  alpha = 1.0;
  F77_CALL(dgemm)("n","t", &dim, model->n,
                  &dim, &alpha, XtX, &dim,
                  pXB, model->n, &beta, XtX1Xt, &dim);



  /* compute splines coefficient: theta = XtX^(-1) %*% Xt %*% y */
  double *theta;
  theta = (double *) R_alloc(dim, sizeof(double));

    for (int i = 0; i < dim; i++){
    temp = 0.0;

    if (i < dim){
      for (int j = 0; j < *model->n; j++){
        temp += XtX1Xt[i + j * dim] * model->y[j];
      }
    }

    theta[i] = temp;
    // Rprintf("%f;;", temp);
  }


  if (retbeta){

    /* compute functional coefficient: beta = basis %*% theta */
    res = (double *) R_alloc(*model->p + model->intercept, sizeof(double));

    /* first entry of returned beta is the constant */
    if (model->intercept == 1){
      res[0] = theta[0];
    }

    for (int i = 0; i < *model->p; i++){
      temp = 0.0;

      for (int j = 0; j < model->selector; j++){
        temp += model->Basis[i + j * *model->p] * theta[j + model->intercept];
      }

      res[i + model->intercept] = temp;
    }

    return res;
  }

  /* return only rss and df */
  res = (double *) R_alloc(2, sizeof(double));

  /* compute RSS */
  double rss = 0.0;

  for (int i = 0; i < *model->n; i++){
    temp = 0.0;
    for (int j = 0; j < dim; j++){
      temp += pXB[i + j* *model->n] * theta[j];
    }
    rss += pow(model->y[i] - temp, 2);
  }
  res[1] = rss;
  res[0] = dim;

  return res;
}




SEXP R_estmodel_spl(SEXP y, SEXP X, SEXP basis, SEXP n, SEXP p, SEXP k, SEXP selector,
                    SEXP intercept, SEXP retbeta) {

  struct callinfo_spl model;


  /* allocate space for intermediate matrix results */
  double *pXB, *XtX, *XtX1Xt;
  pXB = (double *) R_alloc( *INTEGER(n) * (*INTEGER(k) + *INTEGER(intercept)),
                            sizeof(double));
  XtX = (double *) R_alloc( (*INTEGER(k) + *INTEGER(intercept)) *
                            (*INTEGER(k) + *INTEGER(intercept)),
                            sizeof(double));
  XtX1Xt = (double *) R_alloc( *INTEGER(n) * (*INTEGER(k) + *INTEGER(intercept)),
                               sizeof(double));

  model.X = REAL(X);
  model.y = REAL(y);
  model.Basis = REAL(basis);
  model.n = INTEGER(n);
  model.p = INTEGER(p);
  model.k = INTEGER(k);
  model.intercept = *INTEGER(intercept);
  model.selector = *INTEGER(selector);
  model.pXB = pXB;
  model.XtX = XtX;
  model.XtX1Xt = XtX1Xt;

  double *val;
  val = estmodel_spl(&model, *INTEGER(retbeta));

  /* copy result into return SEXP */

  /* functional coefficient */
  SEXP res;
  if (*INTEGER(retbeta)){
    res = PROTECT(allocVector(REALSXP, *model.p + model.intercept));
    for (int i = 0; i < *model.p + model.intercept; i++){
      REAL(res)[i] = val[i];
    }
    UNPROTECT(1);
    return res;
  }

  /* rss and df */
  res = PROTECT(allocVector(REALSXP, 2));

  REAL(res)[0] = val[0];
  REAL(res)[1] = val[1];

  UNPROTECT(1);
  return res;
}
