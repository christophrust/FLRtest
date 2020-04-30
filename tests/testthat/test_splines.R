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


test_that("SplitSpline",{

    grd <- seq(0,1, len = 100)
    knots <- sort(c(rep(0,times = 3), rep(1, times = 3),
                    seq(0,1, by = 0.1), rep(0.5, times = 3)))
    a <- .Call("R_SplitSplineBasis",
          grd = grd,
          df = 16L,
          splitpoint = 0.5)
    b <- splineDesign(knots = knots, x = grd)
    matplot(x = grd, y = a, type = "l")

    expect_equal(a, b)

})


test_that("Spline Model Estimation",{

    if (0) {
        devtools::load_all()
    }

    basis <- .Call("R_SplitSplineBasis",
                   grd = seq(0,1,len = 100),
                   df = 10L,
                   splitpoint = 0.5)


    y <- rnorm(1000)
    X <- matrix(rnorm(1e5), ncol=100)

    pXB <- 1/dim(X)[2] * X %*% basis

    betahat1 <- .Call("R_estmodel_spl",
                     y = y,
                     X = X,
                     basis = basis,
                     n = 1000L,
                     p = 100L,
                     dim = 10L,
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
                 retbeta = 0L)[2]

    expect_equal(rss, sum((y - pXB %*% (solve(crossprod(pXB)) %*% t(pXB)) %*% y)^2))
})
