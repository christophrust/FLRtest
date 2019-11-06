tStatSequence <- function(obj, null, startval, direction, gridvals = NULL) {

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

    
    ## function to compute rho from effDf
    dfGivenRho <- function(x, k) {
        XtX1Xt <- 1/Nobs * chol2inv( chol( obj$model$smspl$npXtX + exp(x) * Amats[,,k])) %*% t(X)
        sum(vapply(1:Nobs, function(i)  sum(X[i,] * XtX1Xt[,i]),0))/p
    }
    

    ## currently hard-coded: startval = 0 and direction = right
    tStatSeq <- vapply(1:(p-1), function(k){

        ## obtain rho
        rho <- exp(
            uniroot( function(x) {dfGivenRho(x, k = k) - df},lower  = -100, upper = 100, f.lower = p,
                    extendInt = "downX")$root
        )
        
        ## estimate full model
        XtX1Xt <- 1/Nobs * chol2inv( chol( obj$model$smspl$npXtX + rho * Amats[,,k])) %*% t(X)
        effDfFull <- sum(vapply(1:Nobs, function(i)  sum(X[i,] * XtX1Xt[,i]),0))/p
        RSSfull <- sum( (X %*% XtX1Xt %*% y * 1/p - y)^2)
        
        ## estimate null model
        selector <- if (intercept) 0:k + 1 else 1:k
        XtX1Xt <- 1/Nobs * chol2inv( chol( obj$model$smspl$npXtX[ selector, selector] +
                                           rho * Amats[,,k][selector, selector])) %*% t(X[,selector])
        effDfNull <- sum(vapply(1:Nobs, function(i)  sum(X[i, selector] * XtX1Xt[,i]),0))/p
        RSSnull <- sum( (X[,selector] %*% XtX1Xt %*% y * 1/p - y)^2 )

        ## compute statistic
        Stat <- ( (RSSnull - RSSfull) / (effDfFull - effDfNull) ) /
            (RSSfull /  (Nobs -  effDfFull - intercept))
        
        ## result
        cbind(Stat, effDfFull, effDfNull, rho)
    }, numeric(4))

    
    ## return
    tStatSeq
}
