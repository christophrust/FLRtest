#ifndef FLRTEST_H
#define FLRTEST_H

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>


struct callinfo {
  double *npXtX; // pointer to matrix 1/(Np) X'X
  double *X;     // pointer to matrix X
  double *Amat;  // pointer to matrix Amat
  double *df;    // pointer to df of full model
  int *n;        // pointer to number of obs
  int *p;        // pointer to dimension of full model (p+1) if intercept is included
  int selector;  // mutable selector (selection subset of model)
};




double R_zeroin2(			/* An estimate of the root */
    double ax,				/* Left border | of the range	*/
    double bx,				/* Right border| the root is seeked*/
    double fa, double fb,		/* f(a), f(b) */
    double (*f)(double x, void *info),	/* Function under investigation	*/
    void *info,				/* Add'l info passed on to f	*/
    double *Tol,			/* Acceptable tolerance		*/
    int *Maxit);


double dfGivenRho(double x, double *npXtX, double *X, double *Amat, int Nobs, int dim);



SEXP R_dfGivenRho(SEXP rho, SEXP npXtX, SEXP X, SEXP Amat, SEXP Nobs, SEXP dim);


double * estmodel(struct callinfo *model, double rho);

SEXP R_estmodel(SEXP npXtX, SEXP X, SEXP Amat, SEXP df, SEXP n, SEXP p, SEXP selector, SEXP logrho);

// SEXP matmult(SEXP a, SEXP b);

#endif

