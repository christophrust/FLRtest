context("Generating sequence of test statistics")


if (0){
    devtools::load_all()


set.seed(123)
p <- 100
N <- 200
grd <- seq(0,1,length = p)
X <- t(bs(x= grd, df = p, intercept = TRUE) %*% matrix(rnorm( N * (p) ), ncol = N))
beta <- sin(10*grd)
y <- X %*% beta/p + rnorm(N, sd=0.1)
obj <- EstFLM(y, X, type = "smoothspline", rho = 0.00000001, intercept = TRUE)


AmSeq <- lapply(3:p, function(k){
    natSplBasis( seq(0, 1, length= k))$A_m * k/p
})

intercept <- TRUE
mdim <- if (intercept) p + 1 else p

Amats <- vapply(1:p, function(k){
    Am <- matrix(0, ncol = mdim, nrow = mdim)
    if (k <3){
        Am[ (k+1):p + intercept, (k+1):p + intercept] <- AmSeq[[ p-k-2 ]]
    } else if (k < (p-2)){
        Am[1:k + intercept, 1:k + intercept] <- AmSeq[[ k-2 ]]
        Am[(k+1):p + intercept, (k+1):p + intercept] <- AmSeq[[ p-k-2 ]]
    } else {
        Am[1:k + intercept, 1:k + intercept] <- AmSeq[[ k-2 ]]
    }
    Am
}, matrix(0,nrow = mdim, ncol = mdim))


obj1 <- list()
class(obj1) <- "abc"


if (0){
    test_that("Input checking",{

        ## correct class of object
        expect_error(tStatSequence(1:10))

        ## correct estimation method
        expect_error(tStatSequence(obj1))

        ##expect_warning(tStatSequence(obj, direction = "left"))
        ##expect_warning(tStatSequence(obj, startval = 1))
        ##expect_warning(tStatSequence(obj, gridvals = (0:99)/99))

    })
}

test_that("c-routine",{

    #print(as.vector(obj$model$smspl$npXtX + exp(-2.3) * Amats[,,1])[1:6])
    #cc <- chol(obj$model$smspl$npXtX + exp(-2.3) * Amats[,,1])
    print(N)
    expect_equal(
        .Call("tstatseq_smspl", y=y, X = X,
              Amats = as.vector(Amats),
              p = as.integer(p),
              n = as.integer(N),
              df = obj$model$effDf,
              npXtX = obj$model$smspl$npXtX,
              tol = 1e-8,
              maxit = 1000L,
              intercept = 1L,
              PACKAGE = "FLRtest")[,6]
       ,
        t(tStatSequenceR(obj))[,1])


    ## a <- matrix(rnorm(9), ncol = 3)
    ## b <- matrix(rnorm(9), ncol = 3)
    ## expect_equal(
    ##     .Call("matmult", a = a, b = b,
    ##           PACKAGE = "FLRtest"),
    ##     a %*% b)
})
}

if (0){
    devtools::load_all()
    mm <- matrix(rnorm(100), ncol = 10)
    mm1 <- matrix(rnorm(100), ncol = 10)
    A <- tcrossprod(mm)
    npXtx <- tcrossprod(mm1)
    yy <- rnorm(20)
    x <- matrix(rnorm(200), ncol = 10)


    aaaa <- .Call("tstatseq_smspl", y=yy, X = x,
                  Amats = A,
                  p = 10L,
                  n = 20L,
                  df = 2,
                  npXtX = npXtx,
                  PACKAGE = "FLRtest")



    .Call("tstatseq_smspl", y=y, X = X,
          Amats = Amats[,,1],
          p = as.integer(p),
          n = as.integer(N),
          df = obj$model$effDf,
          npXtX = obj$model$smspl$npXtX,
          tol = 1e-8,
          maxit = 1000L,
          intercept = 1L,
          PACKAGE = "FLRtest")


    dfGivenRho <- function(lrho , X, npXtX, Amat){

        XtX1Xt <- 1/nrow(X) * chol2inv( chol( npXtX + exp(lrho) * Amat)) %*% t(X)
        sum(vapply(1:nrow(X), function(i)  sum(x[i,] * XtX1Xt[,i]),0))/ncol(X)

    }


    dfGivenRho(0, X, obj$model$smspl$npXtX, Amats[,,1])
}
