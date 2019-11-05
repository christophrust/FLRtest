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
    
    
    ## compute models with increasing support
    RSSfull <- sum(obj$residuals^2)

    ## obtain p and rho from fitted model
    p <- length(obj$coefficients$beta)
    rhoEst <- obj$model$rho

    
    AmSeq <- lapply(9:p, function(k){
        natSplBasis( seq(0, 1, length= k))$A_m * k/p
    })
    
    p <- dim(obj$data$X)[2]
    Nobs <- dim(obj$data$X)[1]
    
    ## currently hard-coded: startval = 0 and direction = right
    RSSnullSeq <- vapply(1:p, function(k) {
        
        ## if k is less than 5, use simple linear model to obtain rss,
        ## otherwise the smoothings spline estimator with same rho
        ## as in original estimation
        if (k <= 8){
            
            ## estimate a simple linear model where regressors are columns 1:k
            rss <- sum(
                .lm.fit(
                    x = if (obj$model$intercept) {
                            cbind(1,obj$data$X[,1:k, drop = FALSE])
                        } else {
                            obj$data$X[,1:k, drop = FALSE]
                        },
                    y = obj$data$y
                )$residuals^2
            )

            ## return
            c(rss, k + obj$model$intercept)
            
        } else {

            ## estimate smoothing spline model
            selector <- if (obj$model$intercept) c(1:k,p+1) else 1:k
            npXtX <- obj$model$smspl$npXtX[selector, selector, drop = FALSE]
            Xsub <- if (obj$model$intercept) cbind(obj$data$X[, 1:k],1) else obj$data$X[, 1:k]
            A_m <- AmSeq[[ k - 8 ]]
            
            rss <- sum(
            ( 1/k * Xsub %*% (xtx1xt <- 1/Nobs * chol2inv( chol( npXtX + obj$model$rho * A_m))
                %*% t(Xsub)) %*% obj$data$y -
                obj$data$y
            )^2)

            edf <- sum(vapply(1:Nobs, function(i)  sum(Xsub[i,] * xtx1xt[,i]),0))/k
            c(rss, edf)
        }
    },numeric(2))

    ## TfSeq <- vapply(RSSnullSeq, function(x){
    ##     ( (x[1] - RSSfull) / (df1 <- obj$model$effDf + obj$model$intercept - intercept) ) /
    ##         (RSSfull / (df2 <- Nobs - obj$model$effDf - obj$model$intercept) )
    ## }, numeric(1))
}
