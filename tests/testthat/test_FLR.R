if (0) {library(FLRtest)
    library(splines)
}
if (1){
    ## simulate some data
    set.seed(123)
    p <- 100
    N <- 1000
    grd <- seq(0,1,length = p)
    X <- t(bs(x= grd, df = p, intercept = TRUE) %*% matrix(rnorm(N*p), ncol = N))
    beta <- sin(10*grd)
    y <- X %*% beta/p + rnorm(N, sd=0.1)
    res <- EstFLM(y, X, type = "smoothspline", rho = 0.00000001, intercept = TRUE)
    res1 <- EstFLM(y, X, type = "fpc", cpv = 0.99, intercept = TRUE)

    ##matplot(y = cbind(beta, res$coefficient$beta,res1$coefficient$beta), x = grd, t = "l")

    context("Functional Regression - Error messages")

    test_that("EstFLM", {
        expect_error(EstFLM())
        expect_error(EstFLM(2,5,type = "abc",rho = 3))
    })



    context("Functional Regression - Estimation")

    test_that("Estimation", {

        ## estimate, given rho
        res1 <- EstFLM(y, X, type = "smoothspline", rho = 0.00000001, intercept = TRUE)
        expect_equal(mean(res1$residuals), 0)

        ## estimate, given df
        res2 <- EstFLM(y, X, type = "smoothspline", df = 6, intercept = TRUE)
        expect_equal(mean(res2$residuals), 0)

        ## fpca estimator
        res3 <- EstFLM(y, X, type = "fpc", df = 6, intercept = TRUE)
        expect_equal(mean(res3$residuals), 0)

        res4 <- EstFLM(y, X, type = "fpc", cpv = 0.9, intercept = TRUE)
        expect_equal(mean(res4$residuals), 0)

        res5 <- EstFLM(y, X, type = "spline", df = 7, intercept = TRUE)
        expect_equal(mean(res5$residuals), 0)

    })



    context("Global Test")

    test_that("Global Test: error messages",{
        expect_error(GlobalTest(res, alternative = c(2,3)))
    })

    test_that("Global Test: plausible values",{
        expect_true(GlobalTest(res)$pval < 0.6)
        expect_true(GlobalTest(res, null = beta)$pval>0.0001)
    })

}
