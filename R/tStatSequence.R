#' Test Sequence for directed Test
#'
#' This function is used by \code{dirTest} and calls the underlying C routine
#' \code{tstatseq}.
#'
#' @param obj A fitted model object of class "flm", as returned by
#' \code{EstFLM}
#'
#' @param null An optional global function to test against. The null
#' can be specified as a vector holding the discretized function object or,
#' alternatively, as an expansion of a given basis system. Defaults to
#' the zero function.
#'
#' @param startval Grad value where the sequential test starts. Currently this
#' has no effect and the procedure starts at the left border of the domain.
#'
#' @param direction One of "right", "left", "both". Specifies into which direction
#' from \code{startval} to procedure is applied.
#'
#' @param gridvals Currently without effect.
#'
#' @param maxit The maximum number of iterations in the root finding routine used
#' to relate the smoothing parameter \equation{\rho} and \code{df}.
#'
#' @param tol Tolerance for the numerical routine used to relate the smoothing
#' parameter \equation{\rho} and \code{df}.
#'
#'
#' @return A matrix with \code{p} rows and 7 columns. The first six columns
#' contain all the information used to calculate the test statistic and the
#' seventh contains the corresponding (approximate) p value.
#'
#'
#'
#' @useDynLib FLRtest
testseq <- function(obj, null, startval, direction, gridvals = NULL,
                    maxit = 1000L, tol = .Machine$double.eps^0.25) {

    ## check that obj is of class flm
    if (!identical(class(obj), "flm")){
        stop("obj must be of class flm!")
    }

    ## check that obj is estimated by smoothing splines
    if (!identical(obj$model$type, "smoothspline")){
        warning("Currently only smoothing spline estimation is supported")
    }

    if (missing(direction))
        direction <- "right"

    if (!(direction %in% c("right","left","both"))) {
        stop("direction must be either of 'right', 'left' or 'both'!")
    }

    if (!is.null(gridvals)) {
        warning("You specified gridvals, this currently does not have an effect!")
    }

    ## temporary warnings for current limitations
    if (!missing(startval) && startval != 0) {
        warning("Currently, only startval = 0 is supported!")
    }

    if (direction != "right") {
        warning("Currently, only direction = 'right' is supported!")
    }

    intercept <- obj$model$intercept

    ## compute models with increasing support
    ## RSSfull <- sum(obj$residuals^2)

    ## obtain p and rho from fitted model
    p <- length(obj$coefficients$beta)
    ## rhoEst <- obj$model$rho
    Nobs <- dim(obj$data$X)[1]
    df <- obj$model$effDf
    X <- if (intercept) cbind(1,obj$data$X) else obj$data$X
    y <- obj$data$y


    ## precompute basis and matrix a for all splits
    ## can be optimized
    AmSeq <- lapply(3:p, function(k){
        natSplBasis( seq(0, 1, length= k))$A_m * k/p
    })

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


    ## compute sequence of test statistic
    tSeq <- .Call("tstatseq",
                  y=y,
                  X = X,
                  Amats = Amats,
                  p = as.integer(p),
                  n = as.integer(Nobs),
                  df = obj$model$effDf,
                  npXtX = obj$model$smspl$npXtX,
                  tol = as.numeric(tol),
                  maxit = as.integer(maxit),
                  intercept = as.integer(intercept),
                  PACKAGE = "FLRtest")



    colnames(tSeq) <- c("edfFull", "rssFull", "edfNull", "rssNull", "logrho", "statistic")


    ## add column p-Value
    tSeq <- cbind(tSeq, pval = pf(tSeq[,"statistic"], df1 = tSeq[,"edfFull"] - tSeq[,"edfNull"],
                   df2 = Nobs - tSeq[,"edfFull"], lower.tail = FALSE))

    ## return
    tSeq
}
