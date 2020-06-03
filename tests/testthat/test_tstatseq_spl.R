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

    for (i in 1:100){
        expect_equal(.Call("tstatseq_spl",
                           y = y,
                           X = X,
                           Basis = as.vector(Basis),
                           selectors = as.integer(Selectors),
                           p = p,
                           n = N,
                           df = as.integer(k + intercept),
                           intercept = intercept),
                     .Call("tstatseq_spl",
                           y = y,
                           X = X,
                           Basis = as.vector(Basis),
                           selectors = as.integer(Selectors),
                           p = p,
                           n = N,
                           df = as.integer(k + intercept),
                           interceptgeorg
                           ))
    }

})


if (0) {
    devtools::unload()
    devtools::load_all()
}


test_that("Simple Example: testseq",{

    k <- 10L
    p <- 100L
    N <- 1000L
    intercept <- 1L
    y <- rnorm(N)
    X <- matrix(rnorm(N * p), ncol=p)

    res_spl <- EstFLM(y = y, X=X, type = "spline", df = 20)
    res_smspl <- EstFLM(y = y, X=X, type = "smoothspline", df = 20)
    ## plot(res_spl)
    ## plot(res_smspl)
    ##debug(testseq)

    pval_spl <- testseq(res_spl)[,"pval"]

    expect_equal(mean(pval_spl < 0.05) < 0.1, TRUE)

    ## pval_smspl <- testseq(res_smspl)[,"pval"]
    ## expect_equal(mean(pval_smspl < 0.05) < 0.1, TRUE)
})


test_that("Sometimes not reproducible!", {

    Nobs <- 250
    p <- 100
    set.seed(1234)

    beta <- pnorm( seq(6,-3, length = p))
    conf.level <- 0.95


    grd <- seq(0,1,length = p)

    ## same curves for all replications
    X <- t(bs(x= grd, df = p, intercept = TRUE) %*% matrix(rnorm(Nobs*p), ncol = Nobs))

    ## also the same for all replications:
    y <- X %*% beta/p + rnorm(Nobs, sd = sqrt(var(X %*% beta/p)/0.5))

    obj <- EstFLM(y, X, type = "spline", df = 8, intercept = TRUE)

    for (i in 1:100){
        expect_equal(dirTest(obj, conf.level = conf.level),
                     dirTest(obj, conf.level = conf.level))
    }

    for (i in 1:100){
        expect_equal(testseq(obj)[,4],
                     testseq(obj)[,4])
    }

    for (i in 1:100){
        print(i)
        f1 <- capture.output(a <- testseq(obj)[,2])
        f2 <- capture.output(b <- testseq(obj)[,2])
        expect_equal(a, b)

    }
    writeLines(f1, "f1.txt")
    writeLines(f2, "f2.txt")
})

for (i in 1:100){

    ## expect_equal(
    ## lapply(grd[-p], function(splitpt){
    ## .Call("R_SplitSplineBasis",
    ##       grd = grd,
    ##       df = as.integer(k) ,
    ##       splitpoint = splitpt)
    ## }),
    ## lapply(grd[-p], function(splitpt){
    ##     .Call("R_SplitSplineBasis",
    ##           grd = grd,
    ##           df = as.integer(k) ,
    ##           splitpoint = splitpt)
    ## }))

    expect_equal(
        .Call("tstatseq_spl",
          y = y,
          X = obj$data$X,
          Basis = as.vector(Basis),
          selectors = as.integer(Selectors),
          p = as.integer(p),
          n = as.integer(length(y)),
          df = as.integer(k + intercept),
          intercept = intercept,
          PACKAGE = "FLRtest")[,2],
        .Call("tstatseq_spl",
              y = y,
              X = obj$data$X,
              Basis = as.vector(Basis),
              selectors = as.integer(Selectors),
              p = as.integer(p),
              n = as.integer(length(y)),
              df = as.integer(k + intercept),
              intercept = intercept,
              PACKAGE = "FLRtest")[,2])

}

devtools::unload()
devtools::load_all()

load("callInfoLast.RData")

for (i in 1:100){
    expect_equal(
        .Call("tstatseq_spl",
              y = y,
              X = X,
              Basis = as.vector(Basis),
              selectors = as.integer(Selectors),
              p = as.integer(p),
              n = as.integer(length(y)),
              df = as.integer(k + intercept),
              intercept = intercept,
              PACKAGE = "FLRtest"),
        .Call("tstatseq_spl",
              y = y,
              X = X,
              Basis = as.vector(Basis),
              selectors = as.integer(Selectors),
              p = as.integer(p),
              n = as.integer(length(y)),
              df = as.integer(k + intercept),
              intercept = intercept,
              PACKAGE = "FLRtest")
    )
}


rss <- .Call("R_estmodel_spl",
             y = y,
             X = X,
             basis = Basis[,,1],
             n = length(y),
             p = dim(X)[2],
             dim = dim(Basis[,,1])[2],
             selector = dim(Basis[,,1])[2],
             intercept = 1L,
             retbeta = 0L)[2]

expect_equal(rss, sum((y - pXB %*% (solve(crossprod(pXB)) %*% t(pXB)) %*% y)^2))
