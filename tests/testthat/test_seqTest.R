context("Generating sequence of test statistics")

set.seed(123)
p <- 100
N <- 1000
grd <- seq(0,1,length = p)
X <- t(bs(x= grd, df = p, intercept = TRUE) %*% matrix(rnorm(N*p), ncol = N))
beta <- sin(10*grd)
y <- X %*% beta/p + rnorm(N, sd=0.1)
obj <- EstFLM(y, X, type = "smoothspline", rho = 0.00000001, intercept = TRUE)

if (0) {
    devtools::load_all()
    debug(tStatSequence)
    aa <- tStatSequence(obj)
    undebug(tStatSequence)
}

obj1 <- list()
class(obj1) <- "abc"


test_that("Input checking",{
    
    ## correct class of object
    expect_error(tStatSequence(1:10))
    
    ## correct estimation method
    expect_error(tStatSequence(obj1))
    
    expect_warning(tStatSequence(obj, direction = "left"))
    expect_warning(tStatSequence(obj, startval = 1))
    expect_warning(tStatSequence(obj, gridvals = (0:99)/99))
    
})
