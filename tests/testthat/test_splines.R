if (0) {
    devtools::load_all()
}

context("Spline Basis")

test_that("Spline Basis Functions",{
    basis1 <- .Call("R_spline_basis",
                    knots = c(0,0,0,0,0.5,1,1,1,1),
                    order = 4L,
                    xvals = seq(0,1, len = 10),
                    derivs = 0L)
    basis2 <- splines::splineDesign(knots = c(0,0,0,0,0.5,1,1,1,1),
                                    ord = 4L,
                                    x = seq(0,1, len = 10),
                                    derivs = 0L)
    expect_equal(basis1, basis2)

    knots <- c(rep(0,time = 4), sort(runif(10)), rep(1,time = 4))
    xvals <- runif(1000)
    basis1 <- .Call("R_spline_basis",
                    knots = knots,
                    order = 4L,
                    xvals = xvals,
                    derivs = 0L)
    basis2 <- splines::splineDesign(knots = knots,
                                    ord = 4L,
                                    x = xvals,
                                    derivs = 0L)
    expect_equal(basis1, basis2)
})


test_that("SplitSpline - correct basis",{

    grd <- seq(0,1, len = 100)

    a <- .Call("R_SplitSplineBasis",
               grd = grd,
               df = 16L,
               splitpoint = 0.5)

    b <- splineDesign(knots = a$knots, x = grd)
    ##matplot(x = grd, y = a, type = "l")

    expect_equal(a$basis, b)


})


test_that("SplitSpline - sequential split", {

    grd <- seq(0,1, by = 0.05)

    for (splitpoint in grd){

        ##(splitpoint <- runif(1, min = 0.1, max = 0.9))
        ## knots <- sort(c(rep(0,times = 3), rep(1, times = 3),
        ##                 seq(0,1, by = 0.1), rep(splitpoint, times = 3)))
        ## splitpoint <- 0.55
        ## splitpoint <- 0.20
        a <- .Call("R_SplitSplineBasis",
                   grd = grd,
                   df = 16L,
                   splitpoint = splitpoint)

        ## test invertibility
        expect_true((det(t(a$basis) %*% a$basis) != 0))
        ##solve(t(a$basis[1:a$selector,1:a$selector]) %*% a$basis[1:a$selector,1:a$selector])


        if (splitpoint <grd[5]){
            rIdx <- which(grd>splitpoint)
            nrIdx <- which(!(grd>splitpoint))
            cIdx <- (a$selector+1):ncol(a$basis)
            ncIdx <- 1:a$selector
            b <- splineDesign(knots = a$knots, x = grd[rIdx])
            basis <- a$basis[rIdx,cIdx]

            ## identity block
            expect_equal(a$basis[nrIdx, ncIdx, drop = FALSE], diag(length(nrIdx)))

        } else if ((splitpoint >grd[length(grd)-4]) && (splitpoint <grd[length(grd)])){
            rIdx <- which(grd<=splitpoint)
            nrIdx <- which(!(grd<=splitpoint))
            cIdx <- 1:a$selector
            ncIdx <- (a$selector+1):ncol(a$basis)
            b <- splineDesign(knots = a$knots, x = grd[rIdx])

            basis <- a$basis[rIdx,cIdx]
            ## identity block
            expect_equal(a$basis[nrIdx, ncIdx, drop = FALSE], diag(length(nrIdx)))

        } else if (isTRUE(all.equal(splitpoint, grd[length(grd)]))){
            rIdx <- seq_along(grd)
            nrIdx <- NULL
            cIdx <- 1:ncol(a$basis)
            ncIdx <- NULL
            b <- splineDesign(knots = a$knots, x = grd[rIdx])
            basis <- a$basis
        } else {
            rIdx <- which(grd<=splitpoint)
            nrIdx <- which(!(grd<=splitpoint))
            cIdx <- 1:a$selector
            ncIdx <- (a$selector+1):ncol(a$basis)
            b <- splineDesign(knots = a$knots, x = grd)
            basis <- a$basis
        }

        ## equality of spline basis
        expect_equal(sum(abs(basis-b)), 0)


        ## remaining zero blocks
        expect_equal(sum(a$basis[nrIdx, cIdx]), 0)
        expect_equal(sum(a$basis[rIdx, ncIdx]), 0)
    }

})

set.seed(1234)
y <- rnorm(1000)
X <- matrix(rnorm(1e5), ncol=100)
SplitBasis <- .Call("R_SplitSplineBasis",
               grd = seq(0,1,len = 100),
               df = 10L,
               splitpoint = 0.5)

test_that("Spline Model Estimation -- Full model",{

    if (0) {
        devtools::load_all()
    }


    basis <- SplitBasis$basis

    pXB <- 1/dim(X)[2] * X %*% basis

    betahat1 <- .Call("R_estmodel_spl",
                     y = y,
                     X = X,
                     basis = basis,
                     n = 1000L,
                     p = 100L,
                     dim = 10L,
                     selector = 10L,
                     retbeta = 1L)

    betahat2 <- as.vector(basis %*% (solve(crossprod(pXB)) %*% t(pXB)) %*% y)

    expect_equal(betahat1, betahat2)

    rss <- .Call("R_estmodel_spl",
                 y = y,
                 X = X,
                 basis = basis,
                 n = 1000L,
                 p = 100L,
                 dim = 10L,
                 selector = 10L,
                 retbeta = 0L)[2]

    expect_equal(rss, sum((y - pXB %*% (solve(crossprod(pXB)) %*% t(pXB)) %*% y)^2))
})



test_that("Spline Model Estimation -- restricted model",{

    if (0) {
        devtools::load_all()
    }


    basis <- SplitBasis$basis

    pXB <- 1/dim(X)[2] * X %*% basis[,1:SplitBasis$selector]

    betahat1 <- .Call("R_estmodel_spl",
                      y = y,
                      X = X,
                      basis = basis,
                      n = 1000L,
                      p = 100L,
                      dim = 10L,
                      selector = as.integer(SplitBasis$selector),
                      retbeta = 1L)

    (solve(t(pXB) %*% pXB) %*% t(pXB))[1:100]

    betahat2 <- as.vector(basis[,1:SplitBasis$selector] %*% (solve(crossprod(pXB)) %*% t(pXB)) %*% y)

    plot(betahat1~betahat2, t = "l")
    expect_equal(betahat1, betahat2)

    rss <- .Call("R_estmodel_spl",
                 y = y,
                 X = X,
                 basis = basis,
                 n = 1000L,
                 p = 100L,
                 dim = 10L,
                 selector = as.integer(SplitBasis$selector),
                 retbeta = 0L)[2]

    expect_equal(rss, sum((y - pXB %*% (solve(crossprod(pXB)) %*% t(pXB)) %*% y)^2))
})
