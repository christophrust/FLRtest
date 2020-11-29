#ifndef FLRTEST_H
#define FLRTEST_H

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>


/* structure containing pointers to the necessary objects for model estimation*/
struct callinfo_smspl {
  double *npXtX; // pointer to matrix 1/(Np) X'X
  double *X;     // pointer to matrix X
  double *y;     // pointer to vector y
  double *Amat;  // pointer to matrix Amat
  double *df;    // pointer to df of full model
  int *n;        // pointer to number of obs
  int *p;        // pointer to dimension of full model (p+1) if intercept is included
  int *dim;      // pointer to dimension of the model (number of colums of X/ npXtX, Amat)
  int selector;  // mutable selector (selecting subset of model)
};


struct callinfo_spl {
  double *X;     // pointer to matrix X
  double *y;     // pointer to vector y
  double *Basis; // pointer to matrix containing evaluated B spline basis. Basis has to be
                 // created using the function SplitSplineBasis()
  double *pXB;   // memory reserved for matrix pXB (n * (k + intercept))
  double *XtX;   // memory reserved for symmetric matrix XtX (k + intercept)
  double *XtX1Xt;// memory reserved for matrix XtX1Xt (n * (k + intercept))
  int *n;        // pointer to number of obs
  int *p;        // pointer to number of grid points
  int *k;        // number of basis functions
  int selector;  // mutable selector (selecting subset of spline basis)
  int intercept; // 0 - no intercept, 1 - intercept in model
};


/* return of function SplitSplineBasis */
typedef struct split_spl_struct {
  double *basis;  // pointer to array containing the evaluated basis
  double *knots;  // full knot sequence used to compute basis
  int selector;   // selector, specifying which basis function is the
                  // last one before split point
  int nk;         // number of knots
} *splitsplPTR;

/* taken from  src/library/stats/src/zeroin.c */
double R_zeroin2(     /* An estimate of the root */
    double ax,        /* Left border | of the range */
    double bx,        /* Right border| the root is seeked*/
    double fa, double fb, /* f(a), f(b) */
    double (*f)(double x, void *info), /* Function under investigation */
    void *info,       /* Add'l info passed on to f */
    double *Tol,      /* Acceptable tolerance */
    int *Maxit);


double dfGivenRho(double x, double *npXtX, double *X, double *Amat, int Nobs, int dim, int p);



SEXP R_dfGivenRho(SEXP rho, SEXP npXtX, SEXP X, SEXP Amat, SEXP Nobs, SEXP dim, SEXP p);


double * estmodel_smspl(struct callinfo_smspl *model, double rho, int retbeta);

SEXP R_estmodel_smspl(SEXP npXtX, SEXP X, SEXP y, SEXP Amat, SEXP df,
                SEXP n, SEXP p, SEXP dim, SEXP selector, SEXP logrho, SEXP retbeta);

void internal_crossprod(double *x, int nrx, int ncx,
                      double *y, int nry, int ncy, double *z);

double * estmodel_spl(struct callinfo_spl *model, int retbeta);
// SEXP matmult(SEXP a, SEXP b);

SEXP R_estmodel_spl(SEXP y, SEXP X, SEXP basis, SEXP n, SEXP p,
                    SEXP k, SEXP selector, SEXP intercept, SEXP retbeta);

double * spline_basis(double *knots, int order, double *xvals, int *derivs,
                      int nk, int nx, int nd);

SEXP tstatseq_smspl(SEXP y, SEXP X,  SEXP Amats, SEXP p, SEXP n,  SEXP df,
                    SEXP npXtX, SEXP tol, SEXP maxit,SEXP intercept);

SEXP R_spline_basis(SEXP knots, SEXP order, SEXP xvals, SEXP derivs);

SEXP R_SplitSplineBasis(SEXP grd, SEXP df, SEXP splitpoint);

SEXP tstatseq_spl(SEXP y, SEXP X,  SEXP Basis, SEXP selectors, SEXP p,
                  SEXP n,  SEXP df, SEXP intercept);
#endif

