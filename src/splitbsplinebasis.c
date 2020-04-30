#include "flrtest.h"

/*
  SplitSplineBasis computes the evaluated basis functions of a cubic spline basis with
  one split point (hence the function space also includes not continouos functions) and
  equidistant knots.

  Inputs:
   - grd: a pointer to a double array containing the grid values (spanning [0,1]) where the
          basis is to be evaluated
   - df: (integer) number of basis functions (dimension of the function space)
   - splitpoint: double value indicating where to split. this not necessarily has to be
         an element of grd but should be and  anything else does not make any sense
         nor will it change results.

  Output:
   A pointer to a double array containing the evaluated basis.
*/

double * SplitSplineBasis(double * grd, int df, double splitpoint, int lgrd){

  double * res;

  int order = 4;
  int dderiv = 0;
  int * deriv = &dderiv;
  int nd = 1;

  /* number of knots (=df + 4 as we use cubic splines) */
  int nk = df + 4;
  /* number of inner knots without splitpoint */
  int nik = nk - 12;
  /* distance between two adjacent knots (+2 because we do not count the split point)*/
  double delta = 1.0/(double) (nik + 2);

  /* knot sequence spanning [0,1] */
  double * knots;
  knots = (double *) R_alloc(nk, sizeof(double));

  int j = 4;
  /* four initial and four end knots */
  for (int i = 0; i < 4; i++) knots[i] = 0.0;
  for (int i = df; i < nk; i++) knots[i] = 1.0;

  /* inner knots */
  for (int i=0; i < nik; i++){
    if ((knots[j-1] < splitpoint) && ((knots[j-1] + delta) >= splitpoint)) {

      if ((splitpoint - knots[j-1]) <= ((knots[j-1] + delta) - splitpoint)) {
        j--;
        for (int k=0; k < 4; k++) {
          knots[j] = splitpoint;
          j++;
        }
        knots[j] = knots[j-5] + 2 * delta; j++;
        knots[j] = knots[j-1] + delta; j++;
      } else {
        for (int k=0; k < 4; k++) {
          knots[j] = splitpoint;
          j++;
        }
        knots[j] = knots[j-5] + 2 * delta; j++;
      }
    } else {
      knots[j] = knots[j-1] + delta;
      j++;
    }
  }
  // for (int i=0; i< nk; i++) Rprintf("%f; ", knots[i]);


  /* compute basis */
  // for (int i=0; i < lgrd; i++) Rprintf("%f; ", grd[i]);

  //Rprintf("nk = %i, lgrd = %i", nk, lgrd);
  res = spline_basis(knots, order, grd, deriv, nk, lgrd, nd);

  // for (int i=0; i< lgrd * df; i++) Rprintf("%f; ", res[i]);
  return res;
}


/* R wrapper */
SEXP R_SplitSplineBasis(SEXP grd, SEXP df, SEXP splitpoint){

  double * val;

  val = SplitSplineBasis(REAL(grd), *INTEGER(df), *REAL(splitpoint), length(grd));

  SEXP res = PROTECT(allocMatrix(REALSXP, length(grd), *INTEGER(df)));

  for (int i = 0; i < length(grd)* *INTEGER(df); i++){
    //Rprintf("%f; ", val[i]);
    REAL(res)[i] = val[i];
  }

  UNPROTECT(1);
  return res;
}
