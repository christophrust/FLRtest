#' EstFLM
#'
#' Estimate a functional linear model
#'
#' @param y Dependent variable, a numeric vector of length N.
#'
#' @param X A N times p matrix holding the functional predictors
#' (one curve per row), assuming an equidistant grid of length
#' \code{p}.
#'
#' @param model A list containing entries \code{type} (either
#' \code{"smoothspline"} or \code{"fpc"} and \code{df},
#' the latter controls flexibility of the fit (either number
#' of eigenfunctions or edf).
#'
#' @param intercept Logical, if TRUE, an intercept is included in the model.
#'
#' @return A list containing estimation results with the
#' following entries:
#' \item{coefficients}{Estimated coefficients, a list
#' with entries \code{coefm} (coefficient matrix),
#' \code{beta} (beta function) and \code{beta_sd} (pointwise standard errors).
#' }
#' \item{residuals}{
#'  Residuals of the fit.
#' }
#' \item{fitted}{
#'  Fitted values
#' }
#' \item{model}{
#'  Parameters of the model fit.
#' }
#'
#' @author Christoph Rust
#' @keywords FDA, smoothing splines, regression
#'
#' @export
EstFLM <- function(y, X, intercept = TRUE, type = "smoothspline", df = NULL, rho = NULL, cpv=NULL){


    if (missing(X)) stop("No functional predictors specified.")
    if (missing(y)) stop("No dependent variable specified.")

    Nobs <- length(y)
    p <- dim(X)[2]
    grd <- seq(0, 1, length.out = p)


    if (type == "smoothspline") {

        ## basis and Amat
        splMat <- natSplBasis(grd)

        
        if (intercept){
            X <- cbind(X,1)
            splMat$A_m <- rbind( cbind(splMat$A_m,0),0)
        }
        
        
        if (!is.null(rho)) {
            
            NPXtX <- crossprod(X) * 1/(Nobs * p)
            XtX1Xt <- 1/Nobs * chol2inv( chol( NPXtX + rho * splMat$A_m)) %*% t(X)
            
            ## solution for all parameters
            beta <- XtX1Xt %*% y
            
            ## number of effective degrees of freedom
            effDf <- sum(vapply(1:Nobs, function(i)  sum(X[i,] * XtX1Xt[,i]),0))/p
            
            yHat <- X %*% beta * 1/p
        } else if (!is.null(df)){

            ## find rho s.t. effdf == df
            NPXtX <- crossprod(X) * 1/(Nobs * p)

            dfGivenRho <- function(x) {
                XtX1Xt <- 1/Nobs * chol2inv( chol( NPXtX + exp(x) * splMat$A_m)) %*% t(X)
                sum(vapply(1:Nobs, function(i)  sum(X[i,] * XtX1Xt[,i]),0))/p
            }

            rho <- exp(uniroot( function(x) {dfGivenRho(x) - df},lower  = -100, upper = 100, f.lower = p,
                               extendInt = "downX")$root)
            
            XtX1Xt <- 1/Nobs * chol2inv( chol( NPXtX + rho * splMat$A_m)) %*% t(X)
            
            ## solution for all parameters
            beta <- XtX1Xt %*% y
            
            ## number of effective degrees of freedom
            effDf <- sum(vapply(1:Nobs, function(i)  sum(X[i,] * XtX1Xt[,i]),0))/p
            
            yHat <- X %*% beta * 1/p
            
        } else {
            ## use GCV to find optimal rho
            ## ... tbd
        }
        
        
    } else if (type == "fpc"){

        ## spectral decomposition
        Cov <- cov(X)
        eigendec <- eigen(Cov)
        
        if (!is.null(df)){
            ## scale by sqrt(p) to be of length 1 in L^2 norm
            efuncs <- eigendec$vectors[,1:df, drop = FALSE] * sqrt(p)
            K <- df
            
        } else if (!is.null(cpv)){
            cumvals <- vapply(1:length(eigendec$values),
                              function(k) sum(eigendec$values[1:k])/sum(eigendec$values),0)
            K <- which.max( cumvals > cpv)
            efuncs <- eigendec$vectors[,1:K, drop = FALSE] * sqrt(p)
        } else {
            ## automatic selection tbd
            ## ...
        }

        scores <- X %*% efuncs * 1/p


        ## estimate approx. model
        lmEst <- .lm.fit(x = if(intercept) cbind(1,scores) else scores,
                         y = y)
        xi <- lmEst$coefficients
        yHat <- (if (intercept) cbind(1,scores) else scores) %*% xi

        ## recover curve
        beta <- c( if(intercept) xi[1] else NULL, efuncs %*% (if(intercept) xi[-1] else xi))

        effDf <- ncol(efuncs)
        
        
    } else{
        stop(sprintf("type '%s' not a valid estimation specification!", model$type))
    }

    ## return stuff
    obj <- list(
        coefficients =
            list(beta = if (intercept) beta[-1] else beta,
                 intercept = if(intercept) beta[1] else NULL),
        residuals = y - yHat,
        fitted = yHat,
        model = list(intercept = intercept,
                     type = type,
                     effDf = effDf,
                     rho = rho,
                     cpv = cpv,
                     fpca = if (type=="fpc") list(efuncs = efuncs,
                                 evals = eigendec$values,
                                 K = K) else NULL,
                     smspl = if (type=="smoothspline") list(npXtX = NPXtX,
                                  A_m = splMat$A_m) else NULL
                     ),
        data = list(y=y, X=X[,-ncol(X), drop = FALSE])
    )
    
    class(obj) <- "flm"
    obj
}
