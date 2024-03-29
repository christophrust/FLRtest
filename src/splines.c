/*
  This code is originally taken from the R package splines and adjusted to our needs.
 */


/*  Routines for manipulating B-splines.  These are intended for use with
 *  S or S-PLUS or R.
 *
 *     Copyright (C) 1998 Douglas M. Bates and William N. Venables.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 *
 * These functions are distributed in the hope that they will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the
 * GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 * The routines are loosely based on the pseudo-code in Schumacher (Wiley,
 * 1981) and the CMLIB library DBSPLINES.
 */

#include <R.h>
#include <Rinternals.h>
#include <string.h> // for memcpy


typedef struct spl_struct {
  int order,			/* order of the spline */
    ordm1,			/* order - 1 (3 for cubic splines) */
    nknots,			/* number of knots */
    curs,			/* current position in knots vector */
    boundary;		/* must have knots[curs] <= x < knots[curs+1] */
  /* except for the boundary case */

  double *ldel,		/* differences from knots on the left */
    *rdel,			/* differences from knots on the right */
    *knots,			/* knot vector */
    *coeff,			/* coefficients */
    *a;			/* scratch array */
} *splPTR;


/* set sp->curs to the index of the first knot position > x.
   Special handling for x == sp->knots[sp->nknots - sp-order + 1] */
static int
set_cursor(splPTR sp, double x)
{
  int i;
  /* don't assume x's are sorted */

  sp->curs = -1; /* Wall */
  sp->boundary = 0;
  for (i = 0; i < sp->nknots; i++) {
    if (sp->knots[i] >= x) sp->curs = i;
    if (sp->knots[i] > x) break;
  }
  if (sp->curs > sp->nknots - sp->order) {
    int lastLegit = sp->nknots - sp->order;
    if (x == sp->knots[lastLegit]) {
	    sp->boundary = 1; sp->curs = lastLegit;
    }
  }
  return sp->curs;
}

static void
diff_table(splPTR sp, double x, int ndiff)
{
  int i;
  for (i = 0; i < ndiff; i++) {
    sp->rdel[i] = sp->knots[sp->curs + i] - x;
    sp->ldel[i] = x - sp->knots[sp->curs - (i + 1)];
  }
}

/* fast evaluation of basis functions */
static void
basis_funcs(splPTR sp, double x, double *b)
{
  diff_table(sp, x, sp->ordm1);
  b[0] = 1.;
  for (int j = 1; j <= sp->ordm1; j++) {
    double saved = 0.;
    for (int r = 0; r < j; r++) { // do not divide by zero
	    double den = sp->rdel[r] + sp->ldel[j - 1 - r];
	    if(den != 0) {
        double term = b[r]/den;
        b[r] = saved + sp->rdel[r] * term;
        saved = sp->ldel[j - 1 - r] * term;
	    } else {
        if(r != 0 || sp->rdel[r] != 0.)
          b[r] = saved;
        saved = 0.;
	    }
    }
    b[j] = saved;
  }
}

/* "slow" evaluation of (derivative of) basis functions */
static double
evaluate(splPTR sp, double x, int nder)
{
  register double *lpt, *rpt, *apt, *ti = sp->knots + sp->curs;
  int inner, outer = sp->ordm1;

  if (sp->boundary && nder == sp->ordm1) { /* value is arbitrary */
    return 0.0;
  }
  while(nder--) {  // FIXME: divides by zero
    for(inner = outer, apt = sp->a, lpt = ti - outer; inner--; apt++, lpt++)
	    *apt = outer * (*(apt + 1) - *apt)/(*(lpt + outer) - *lpt);
    outer--;
  }
  diff_table(sp, x, outer);
  while(outer--)
    for(apt = sp->a, lpt = sp->ldel + outer, rpt = sp->rdel, inner = outer + 1;
        inner--; lpt--, rpt++, apt++)
	    // FIXME: divides by zero
	    *apt = (*(apt + 1) * *lpt + *apt * *rpt)/(*rpt + *lpt);
  return sp->a[0];
}

/* called from	predict.bSpline() and predict.pbSpline() : */

/* called from	splineDesign() : */
double *
spline_basis(double *knots, int order, double *xvals, int *derivs,
             int nk, int nx, int nd)
{
  /* evaluate the non-zero B-spline basis functions (or their derivatives)
   * at xvals.  */


  // PROTECT(knots = coerceVector(knots, REALSXP));
  // double *kk = REAL(knots); int nk = length(knots);
  // int ord = asInteger(order);
  // PROTECT(xvals = coerceVector(xvals, REALSXP));
  // double *xx = REAL(xvals); int nx = length(xvals);
  // PROTECT(derivs = coerceVector(derivs, INTSXP));
  // int *ders = INTEGER(derivs), nd = length(derivs);

  splPTR sp = (struct spl_struct *) R_alloc(1, sizeof(struct spl_struct));

  /* fill sp : */
  sp->order = order;
  sp->ordm1 = order - 1;
  sp->rdel = (double *) R_alloc(sp->ordm1, sizeof(double));
  sp->ldel = (double *) R_alloc(sp->ordm1, sizeof(double));
  sp->knots = knots; sp->nknots = nk;
  sp->a = (double *) R_alloc(sp->order, sizeof(double));

  double * val, *res;
  int *offsets;

  val = (double *) R_alloc((sp->order * nx), sizeof(double));
  res = (double *) R_alloc(((nk - sp->order) * nx), sizeof(double));
  offsets = (int *) R_alloc(nx, sizeof(int));

  // initialize array res to zero
  for (int i=0; i < ((nk - sp->order) * nx); i++) res[i] = 0.0;

  /* SEXP val = PROTECT(allocMatrix(REALSXP, sp->order, nx)), */
  /* offsets = PROTECT(allocVector(INTSXP, nx)); */
  /* double *valM = REAL(val); */
  /* int *ioff = INTEGER(offsets); */

  for(int i = 0; i < nx; i++) {
    set_cursor(sp, xvals[i]);
    int io = offsets[i] = sp->curs - sp->order;
    if (io < 0 || io > nk) {
	    for (int j = 0; j < sp->order; j++) {
        val[i * sp->order + j] = R_NaN;
	    }
    } else if (derivs[i % nd] > 0) { /* slow method for derivatives */
	    for(int ii = 0; ii < sp->order; ii++) {
        for(int j = 0; j < sp->order; j++) sp->a[j] = 0;
        sp->a[ii] = 1;
        val[i * sp->order + ii] =
          evaluate(sp, xvals[i], derivs[i % nd]);
	    }
    } else {		/* fast method for value */
	    basis_funcs(sp, xvals[i], val + i * sp->order);
    }
  }
  //Rprintf("test1\n");
  // construct the final matrix according to R function splineDesign assuming
  // no outer knots and dense design
  int cnt=0;
  for (int i = 0; i < nx; i++){
    for (int j=0; j<order; j++){
      /* Rprintf("i=%i, j=%i, arraypos = %i, val = %f\n", */
      /*         i, */
      /*         j + offsets[i], */
      /*         i + (j + offsets[i]) * nx, */
      /*         val[cnt]); */

      res[ i + (j + offsets[i]) * nx ] = val[cnt];
      cnt++;
    }
  }
  // Rprintf("test1\n");
  // setAttrib(val, install("Offsets"), offsets);
  // UNPROTECT(5);
  return res;
}

// Wrapper around above function to make code accessible from R
SEXP R_spline_basis(SEXP knots, SEXP order, SEXP xvals, SEXP derivs){

  // we assume, the arguments are provided properly, therefore no checking
  int nk = length(knots);
  int nx = length(xvals);
  int nd = length(derivs);

  SEXP res = PROTECT(allocMatrix(REALSXP, nx, nk - *INTEGER(order)));
  double *val;

  val = spline_basis(REAL(knots), *INTEGER(order), REAL(xvals), INTEGER(derivs), nk, nx, nd);

  // copy the values as the pointer seems to be difficult to access from outside
  for (int i=0; i < ((nk - *INTEGER(order)) * nx); i++){
    REAL(res)[i] = val[i];
  }

  UNPROTECT(1);
  return res;
}
