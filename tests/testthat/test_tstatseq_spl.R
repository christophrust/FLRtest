if (0) {
    devtools::load_all()
}

set.seed(1234)
k <- 10L
p <- 100L
N <- 1000L
intercept <- 0L

y <- rnorm(N)
X <- matrix(rnorm(N * p), ncol=p)
grd <- seq(0,1,len = p)


BasisAndSelectors <- lapply(grd[-p], function(splitpt){
    .Call("R_SplitSplineBasis",
          grd = grd,
          df = k ,
          splitpoint = splitpt)
})

a <- .Call("R_SplitSplineBasis",
      grd = grd,
      df = k ,
      splitpoint = grd[99])


Basis <- vapply(BasisAndSelectors, function(x) x$basis, matrix(0, nrow = p, ncol = k))

## check that all basis matrices are full-rank:
fullrank <- vapply(BasisAndSelectors, function(x) abs(det(crossprod(x$basis)))>0, TRUE)

context("tstatseq_spl")

test_that("simple call tstatseq_spl",{

    .Call("tstatseq_spl",
          y = y,
          X = X,
          Basis = as.vector(Basis),
          selectors = as.integer(Selectors),
          p = p,
          n = N,
          df = as.integer(k + intercept),
          intercept = intercept)


})
