if (0) {
    devtools::unload()
    devtools::load_all()
}

##set.seed(1234)
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

Basis <- vapply(BasisAndSelectors, function(x) x$basis, matrix(0, nrow = p, ncol = k))
Selectors <- vapply(BasisAndSelectors, function(x) x$selector, 0)
## check that all basis matrices are full-rank:
fullrank <- vapply(BasisAndSelectors, function(x) abs(det(crossprod(x$basis)))>0, TRUE)

context("tstatseq_spl")

test_that("simple call tstatseq_spl",{

    intercept = 0L
    a <- .Call("tstatseq_spl",
               y = y,
               X = X,
               Basis = as.vector(Basis),
               selectors = as.integer(Selectors),
               p = p,
               n = N,
               df = as.integer(k + intercept),
               intercept = intercept)

    intercept = 1L
    b <- .Call("tstatseq_spl",
               y = y,
               X = X,
               Basis = as.vector(Basis),
               selectors = as.integer(Selectors),
               p = p,
               n = N,
               df = as.integer(k + intercept),
               intercept = intercept)

    a_pvals <- pf(a[,5], df1 = a[,1] - a[,3],
                  df2 = N - a[,1], lower.tail = FALSE)
    b_pvals <- pf(b[,5], df1 = b[,1] - b[,3],
                  df2 = N - b[,1], lower.tail = FALSE)


    expect_true(mean(a_pvals < 0.01)<=0.01)
    expect_true(mean(b_pvals < 0.01)<=0.01)

})


if (0) {
    devtools::unload()
    devtools::load_all()
}

k <- 10L
p <- 100L
N <- 1000L
intercept <- 1L
y <- rnorm(N)
X <- matrix(rnorm(N * p), ncol=p)

res <- EstFLM(y = y, X=X, type = "spline", df = 20)
plot(res)

debug(testseq)
testseq(res)
