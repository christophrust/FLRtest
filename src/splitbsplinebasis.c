#include "flrtest.h"

/*
  helper function to construct an array containing in the first or last columns
  in the upper (or lower) block an identity matrix.

  Inputs:
   - grd: a pointer to a double array, containing the full grid values
   - startvalidx: integer, specifying, at which grid val to start the spline basis
   - endvalidx: integer, specifying, at which grid val to end the spline basis
   - lgrd: integer, length of the array grd
   - df: dimension of the full basis

  Output:
   - a pointer to a double array of dimension lgrd * df where if startvalidx > 0 the first
     startvalidx columns will contain in the upper block a diagonal matrix and the...
 */
double * SimpleSplineBasis(double *grd, int startvalidx, int endvalidx, int lgrd, int df){

  /* initialize return */
  double * res;
  res = (double *) R_alloc(lgrd * df, sizeof(double));

  /* some objects necessary for the call of spline_basis */
  int order = 4;
  int dderiv = 0;
  int * deriv = &dderiv;
  int nd = 1;

  int nk = df + 4 - startvalidx - lgrd + endvalidx - 1;
  int nik = nk - 8;
  double *knots;
  knots = (double *) R_alloc(nk, sizeof(double));

  double delta = (grd[endvalidx]-grd[startvalidx])/(double) (nik + 1);

  /* start and end points of knots */
  for (int i = 0; i < 4; i++) knots[i] = grd[startvalidx];
  for (int i = (nk-4); i < nk; i++) knots[i] = grd[endvalidx];

  for (int i = 0; i < nik; i++){
    knots[i+4] = knots[i+3] + delta;
  }

}


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

splitsplPTR SplitSplineBasis(double * grd, int df, double splitpoint, int lgrd){

  double * basis;

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

  /*
   selector holds the index of the last non-zero basis
   function of the first part of the domain.
  */
  int selector = 4;

  /* we assume split point is always at least as large as the first inner knots! */
  /* inner knots */
  double nextknot = delta;

  for (int i=0; i < (nik + 1); i++){


    /* check whether nextknot has to be replaced by the 4 splitpoints */
    if (((nextknot - 0.5 * delta) < splitpoint) && ((nextknot + 0.5 * delta) >= splitpoint)){


      /* replace nextknot by splitpoint */
      for (int k=0; k < 4; k++) {
        knots[j] = splitpoint;
        j++;
      }

      nextknot = nextknot + delta;

    } else {
      knots[j] = nextknot;
      j++;
      if (knots[j-1] < splitpoint) {
        selector++;
      }

      nextknot = nextknot + delta;
    }
  }


  /* compute basis */
  basis = spline_basis(knots, order, grd, deriv, nk, lgrd, nd);

  splitsplPTR res = (struct split_spl_struct *) R_alloc(1, sizeof(struct split_spl_struct));

  res->basis = basis;
  res->selector = selector;
  res->knots = knots;

  return res;
}


/* R wrapper */
SEXP R_SplitSplineBasis(SEXP grd, SEXP df, SEXP splitpoint){

  splitsplPTR val;

  const char *names[] = {"basis", "selector", "knots", ""};

  val = SplitSplineBasis(REAL(grd), *INTEGER(df), *REAL(splitpoint), length(grd));

  SEXP res = PROTECT(mkNamed(VECSXP,names));
  SEXP basis = PROTECT(allocMatrix(REALSXP, length(grd), *INTEGER(df)));
  SEXP selector = PROTECT(allocVector(INTSXP, 1));
  SEXP knots = PROTECT(allocVector(REALSXP, *INTEGER(df) + 4));

  /* copy content into SEXPs */
  *INTEGER(selector) = val->selector;

  double * pbasis = REAL(basis);
  for (int i = 0; i < length(grd)* *INTEGER(df); i++){
    //Rprintf("%f; ", val[i]);
    pbasis[i] = val->basis[i];
  }

  for (int i = 0; i < (*INTEGER(df) + 4); i++){
    REAL(knots)[i] = val->knots[i];
  }

  SET_VECTOR_ELT(res, 0 , basis);
  SET_VECTOR_ELT(res, 1 , selector);
  SET_VECTOR_ELT(res, 2 , knots);

  UNPROTECT(4);
  return res;
}
