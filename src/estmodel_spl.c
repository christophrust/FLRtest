#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
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

  /* return array */
  double *res;

  /* allocate space for necessary intermediate results */
  double *pXB;   // 1/p * X * Basis
  pXB = (double *) R_alloc( *model->n * *model->p, sizeof(double));




  return res;
}
