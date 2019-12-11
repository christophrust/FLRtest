#ifndef FLRTEST_H
#define FLRTEST_H

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>


struct callinfo {
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




double R_zeroin2(			/* An estimate of the root */
    double ax,				/* Left border | of the range	*/
    double bx,				/* Right border| the root is seeked*/
    double fa, double fb,		/* f(a), f(b) */
    double (*f)(double x, void *info),	/* Function under investigation	*/
    void *info,				/* Add'l info passed on to f	*/
    double *Tol,			/* Acceptable tolerance		*/
    int *Maxit);


double dfGivenRho(double x, double *npXtX, double *X, double *Amat, int Nobs, int dim, int p);



SEXP R_dfGivenRho(SEXP rho, SEXP npXtX, SEXP X, SEXP Amat, SEXP Nobs, SEXP dim, SEXP p);


double * estmodel(struct callinfo *model, double rho, int retbeta);

SEXP R_estmodel(SEXP npXtX, SEXP X, SEXP y, SEXP Amat, SEXP df,
                SEXP n, SEXP p, SEXP dim, SEXP selector, SEXP logrho, SEXP retbeta);

// SEXP matmult(SEXP a, SEXP b);

#endif

