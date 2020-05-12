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
splitsplPTR SimpleSplineBasis(double *grd, int startvalidx, int endvalidx, int lgrd, int df){

  /* either startval or endval may deviate from beginning/end */
  if (startvalidx * (lgrd - 1 - endvalidx) != 0)
    error("Identity basis either at beginning OR end of domain!");

  /* initialize return */
  double * basis;
  basis = (double *) R_alloc(lgrd * df, sizeof(double));

  /* some objects necessary for the call of spline_basis */
  int order = 4;
  int deriv = 0;
  int nd = 1;
  int lgrd_short = endvalidx - startvalidx + 1;

  int dim_id_block = startvalidx + (lgrd - endvalidx - 1);
  int nk = df + 4 - dim_id_block;
  int nik = nk - 8;

  double *knots;
  knots = (double *) R_alloc(nk, sizeof(double));

  double delta = (grd[endvalidx] - grd[startvalidx])/(double) (nik + 1);

  /* start and end points of knots */
  for (int i = 0; i < 4; i++) knots[i] = grd[startvalidx];
  for (int i = (nk-4); i < nk; i++) knots[i] = grd[endvalidx];

  /* remaining inner knot sequence */
  for (int i = 0; i < nik; i++){
    knots[i+4] = knots[i+3] + delta;
  }

  /* spline basis object */
  double *spl_basis;
  spl_basis = (double *) R_alloc(lgrd_short * (nk-4), sizeof(double));
  spl_basis = spline_basis(knots, order, (grd + startvalidx), &deriv, nk, lgrd_short, nd);

  int basisIdx, selector;
  /* fill final matrix basis */
  if (startvalidx > 0){

    selector = dim_id_block;

    /* fill upper left dim_id_block */
    for (int i=0; i < pow(dim_id_block,2); i++){
      basisIdx = i + lgrd_short * (i/dim_id_block);
      if ((i % (dim_id_block + 1)) != 0){
        basis[basisIdx] = 1.0;
      } else {
        basis[basisIdx] = 0.0;
      }
    }

    /* fill block containing spline basis */
    for (int i=0; i < (lgrd_short * (nk-4)); i++){
      basisIdx = dim_id_block * (lgrd + 1 + i/(lgrd - dim_id_block)) + i;
      basis[basisIdx] = spl_basis[i];
    }

    /* fill remaining zero blocks */
    /* block left of basis */
    for (int i=0; i<(dim_id_block * lgrd_short); i++){
      basisIdx = (i/lgrd_short + 1) * dim_id_block + i;
      basis[basisIdx] = 0.0;
    }
    /* block above basis */
    for (int i=0; i<(dim_id_block * (nk-4)); i++){
      basisIdx = lgrd * dim_id_block + (i/dim_id_block) * lgrd_short + i; // integer division!
      basis[basisIdx] = 0.0;
    }

  } else if ((lgrd - 1 - endvalidx) > 0) {

    selector = lgrd_short;

    /* fill upper left block with spline basis */
    for (int i=0; i < (lgrd_short * (nk-4)); i++){
      basisIdx = (i/lgrd_short) * dim_id_block + i; // integer division!
      basis[basisIdx] = spl_basis[i];
    }


    /* fill lower right block with identity matrix */
    for (int i = 0; i< pow(dim_id_block, 2); i++){
      basisIdx = lgrd * (nk - 4) + (i/dim_id_block + 1) * lgrd_short + i; // integer division!
      if ((i % (dim_id_block + 1)) != 0){
        basis[basisIdx] = 1.0;
      } else {
        basis[basisIdx] = 0.0;
      }
    }

    /* fill remaining blocks with zeros */

    /* block right of basis */
    for (int i=0; i<(dim_id_block * lgrd_short); i++){
      basisIdx = lgrd * (nk - 4) + (i/lgrd_short) * dim_id_block + i; // integer division!
      basis[basisIdx] = 0.0;
    }
    /* block right of basis */
    for (int i=0; i<(dim_id_block * (nk-4)); i++){
      basisIdx = lgrd * (nk - 4) + (i/lgrd_short) * dim_id_block + i; // integer division!
      basis[basisIdx] = 0.0;
    }
  }

  splitsplPTR res = (struct split_spl_struct *) R_alloc(1, sizeof(struct split_spl_struct));

  res->basis = basis;
  res->selector = selector;
  res->knots = knots;
  res->nk = nk;

  return res;
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
         nor will it change results. Splitpoint is still part of the first part
         of the domain.

  Output:
   A pointer to a split_spl_struct (see flrtest.h), inlcuding i.a. a pointer to a
   double array containing the evaluated basis.

   If splitspoint is smaller than grd[3] (or larger than grd[lgrd-4]), the first
   (or last) columns of the resulting array contain in the uppermost (lowermost) block
   a diagonal matrix.

*/

splitsplPTR SplitSplineBasis(double * grd, int df, double splitpoint, int lgrd){

  double * basis;


  /* simple spline basis (diagonal block at beginning or end) */
  if ((splitpoint < grd[4]) | (splitpoint > grd[lgrd-5])) {

    /* determine startvalidx/envalidx */
    int i = 0, j = 0;
    while ((splitpoint < grd[i]) | (splitpoint > grd[lgrd-1-j])){
      j++; i++;
    }
    int startvalidx, endvalidx;
    if (i>j) {
      startvalidx = i;
      endvalidx = lgrd-1;
    } else {
      startvalidx = 0;
      endvalidx = lgrd- 1 - j;
    }

    splitsplPTR res;
    res = SimpleSplineBasis(grd, startvalidx, endvalidx, lgrd, df);
    return res;
  }


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
  res->nk = nk;

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
  SEXP knots = PROTECT(allocVector(REALSXP, val->nk));

  /* copy content into SEXPs */
  *INTEGER(selector) = val->selector;

  double * pbasis = REAL(basis);
  for (int i = 0; i < length(grd)* *INTEGER(df); i++){
    //Rprintf("%f; ", val[i]);
    pbasis[i] = val->basis[i];
  }

  for (int i = 0; i < (val->nk); i++){
    REAL(knots)[i] = val->knots[i];
  }

  SET_VECTOR_ELT(res, 0 , basis);
  SET_VECTOR_ELT(res, 1 , selector);
  SET_VECTOR_ELT(res, 2 , knots);

  UNPROTECT(4);
  return res;
}
