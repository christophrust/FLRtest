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
})
