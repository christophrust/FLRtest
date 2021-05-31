#include "flrtest.h"


/*
  tstatseq_spl computes the test sequence described in the manuscript
  "directed local tesing in functional regression".

  Inputs:
   - y: numeric vector of length n containing observations of dependent variable
   - X: numeric matrix n times (p + intercept) containing discretized functions
        (every row is one observation)
   - Basis: an (optional) stacked (as.vector) numeric array (p times k times (p-1))
        containing all the basis objects for all splits
   - selectors: an (optional) vector containing the split indices of the basis objects
   - p: integer, number of discretization points.
   - n: number of observations.
   - df: number of degrees of freedom for full model
   - intercept: integer 0 - no constant in model, 1 - constant included

  Output:
   - A numeric matrix with p-1 rows and 5 columns:
     1. rss of the full model
     2. edf of the full model
     3. rss of null model
     4. edf of null model
     5. test statistic (see manuscript)
 */

SEXP tstatseq_spl(SEXP y, SEXP X,  SEXP Basis, SEXP selectors, SEXP p,
                  SEXP n,  SEXP df, SEXP intercept){

  /* obtain input variables and initialize pointers to these objects */
  int np = INTEGER(p)[0];
  int intercpt = INTEGER(intercept)[0];
  int k = *INTEGER(df);


  /* initialize return object */
  SEXP res;
  res = PROTECT(allocMatrix(REALSXP, np-1, 5));

  /* containers for intermediate results */
  //double logrho;
  double *fullmodel, *nullmodel;

  /* allocate space for intermediate matrix results */
  double *pXB, *XtX, *XtX1Xt;
  pXB = (double *) R_alloc( *INTEGER(n) * (k + intercpt), sizeof(double));
  XtX = (double *) R_alloc( (k + intercpt) * (k + intercpt), sizeof(double));
  XtX1Xt = (double *) R_alloc( *INTEGER(n) * (k + intercpt), sizeof(double));

  /* callinfo_spl structure containing all the necessary information for the model fit */
  struct callinfo_spl info;

  info.X = REAL(X);
  info.y = REAL(y);
  info.Basis = REAL(Basis);
  info.n = INTEGER(n);
  info.k = &k;
  info.p = INTEGER(p);
  info.intercept = intercpt;
  info.pXB = pXB;
  info.XtX = XtX;
  info.XtX1Xt = XtX1Xt;

  /* main iteration */
  for (int j=0; j < (np-1); j++){

    /* estimate full model */
    info.selector = k;
    fullmodel = estmodel_spl(&info, 0);

    // Rprintf("Full: [%i] df: %f, rss: %f\n",j,fullmodel[0], fullmodel[1]);

    REAL(res)[j] = fullmodel[0];          // 1st col: rss of full model
    REAL(res)[j + (np-1)] = fullmodel[1]; // 2nd col: df of full model

    /* estimate null model */
    info.selector = INTEGER(selectors)[j];


    nullmodel = estmodel_spl(&info, 0);

    // Rprintf("Null: [%i] df: %f, rss: %f\n",j,nullmodel[0], nullmodel[1]);

    REAL(res)[j + 2 * (np-1)] = nullmodel[0]; // 3rd col: rss of null model
    REAL(res)[j + 3 * (np-1)] = nullmodel[1]; // 4th col: df of null model

    /* compute test statistic: ( (rss0 - rss1)/df1-df0) / (rss1 / (n - df1) */
    REAL(res)[j +4*(np-1)] =
      ( (REAL(res)[j +3*(np-1)] - REAL(res)[j +1*(np-1)]) /
        (REAL(res)[j +0*(np-1)] - REAL(res)[j +2*(np-1)])) /
      (REAL(res)[j +1*(np-1)] /
       ( (double) *INTEGER(n) - REAL(res)[j +0*(np-1)]));



    if (j < (np-2)){
      info.Basis += *INTEGER(p) * k;
    }
  }



  UNPROTECT(1);
  return res;
}
