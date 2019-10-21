#' EstFLM
#'
#' Estimate a functional linear model
#'
#' @param y Dependent variable, a numeric vector of length N.
#' @param X A N times p matrix holding the functional predictors (one curve per row)
#' @param model A list containing entries \code{type} (either \code{"smoothspline"} or \code{"fpc"} and \code{df},
#' the latter controls flexibility of the fit (either number of eigenfunctions or edf).
#' @param intercept Logical, if TRUE, an intercept is included in the model.
#'
#' @return A list containing estimation results with the following entries:
#' \item{coefficients}{Estimated coefficients, a list with entries \code{coefm} (coefficient matrix),
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
#' @author Johann Eppelsheimer, Christoph Rust
#' @keywords FDA, smoothing splines, regression
#' @importFrom stats pnorm
#'
#' @export
EstFLM <- function(y,X, model = list(type = "smoothspline", df = NULL){


    if (missing(X)) stop("No functional predictors specified.")
    if (missing(y)) stop("No dependent variable specified.")

    Nobs <- length(y)
    p <- dim(X)[2]

}

oldfun <-function(
    y, Curves = NULL, ScalarPred = NULL, lostDF = 0,
                                clusterVar = NULL, clusterType = NULL, rho = 1e-6 ,intercept = TRUE) {

    
    if (is.null(Curves) && is.null(ScalarPred)) stop("No predictor variables specified.")


    Nobs <- length(y)
    

    
    if (!is.null(Curves)) {
        ## grid where functions are evaluated
        p <- if (is.matrix(Curves)) dim(Curves)[2] else dim(Curves[[1]])[2]
    }
    ## design matrix
    if (intercept || !is.null(ScalarPred)){
        Xvars <- if (intercept) {
                     if (!is.null(ScalarPred)) {
                         as.matrix(cbind(1, as.matrix(ScalarPred)))
                     } else {
                         matrix(1, nrow = Nobs, ncol = 1)
                     }
                 } else {
                     as.matrix(ScalarPred)
                 }
    } else {
        Xvars <- NULL
    }
    if (intercept) colnames(Xvars)[1] <-"Constant"
    Xfull <- cbind(if (is.null(Curves)) NULL else if (is.matrix(Curves)) Curves else do.call(cbind, Curves),
                   p * Xvars)

    
    Yfull <- y
    
    dimen <- ncol(Xfull)
    Afull <- matrix(0, ncol = dimen, nrow = dimen)

    ## spline Basis and Amat
    if (!is.null(Curves)) {

        
        ## check for correct input
        if (is.list(Curves)){
            dims <- vapply(Curves , dim, numeric(2))
            if (any(dims[2,] != p) || any(dims[1,] != Nobs)) stop("Elements in 'Curves' have different dimensions!")
        }
        
        if (!is.null(ScalarPred) && dim(ScalarPred)[1] != Nobs) stop("'ScalarPred' has wrong number of observations")
        
        grd <- seq(0, 1, len = p)

        
        splMat <- calSecDerNatSpline(grd)
        if (is.matrix(Curves)) {
            Afull[1:p,1:p] <- splMat[["A_m"]]
            ScPredStartIdx <- p+1
        } else{
            for (i in 0:(length(Curves)-1)) Afull[i*p + 1:p, i*p + 1:p] <- splMat[["A_m"]]
            ScPredStartIdx <- length(Curves) * p +1
        }
    } else {
        ScPredStartIdx <- 1
    }
    
    NPXtX <- crossprod(Xfull) * 1/(Nobs * p)
    XtX1Xt <- (XtX1 <- 1/Nobs * chol2inv( chol( NPXtX + rho * Afull))) %*% t(Xfull)

    ## solution for all parameters
    beta_full <- XtX1Xt %*% Yfull
    
    
    ## number of effective degrees of freedom
    effDf <- sum(vapply(1:Nobs, function(i)  sum(Xfull[i,] * XtX1Xt[,i]),0))/p
    
    yHat <- Xfull %*% (beta_full * 1/p)
    
    ## variance of estimator
    if (is.null(clusterVar)){
        ssr <- sum( (resid <- Yfull - yHat)^2)
        var_beta <- diag( tcrossprod(XtX1Xt) ) * 1/(Nobs - lostDF - effDf) * ssr   # also substract DF form FE-estimates ('lostDF')
    } else {
        if (is.null(clusterType) || clusterType == "LZ"){
            ## LZ cluster version, see AbadieAtheyImbensWoolridge, eq. 2.3
            resid <- as.vector(Yfull - yHat)
            
            ##XtOmegaX <- rowSums(vapply(unique(clusterVar), function(c){
            ##    idx <- (clusterVar == c)
            ##    tcrossprod(colSums(resid[idx] * Xfull[idx,,drop=FALSE]))
            ##}, matrix(0,dimen,dimen)), dims = 2)
            #XtOmegaX <- RobClusterSEs(resid, Xfull, clusterVar, unique(clusterVar))
            srtIdx <- order(clusterVar)
            uCl <- sort(unique(clusterVar))
            clustersAndTimes <- cbind(as.integer(uCl), as.integer(tabulate(clusterVar))[uCl])
            XtOmegaX <- RobClusterSEs( resid[srtIdx], Xfull[srtIdx,], clustersAndTimes)

            var_beta <- diag(XtX1 %*% XtOmegaX %*% XtX1)
        } else if (clusterType == "EHW"){
            ## Eicker-Huber-White
            resid <- as.vector(Yfull - yHat)
            XtOmegaX <- rowSums(vapply(1:length(resid), function(i){
                tcrossprod(resid[i] * Xfull[i,, drop = FALSE])
            }, matrix(0,dimen,dimen)), dims = 2)
            var_beta <- diag(XtX1 %*% XtOmegaX %*% XtX1)
        } else if (clusterType =="KLOEK"){
            resid <- as.vector(Yfull - yHat)
            
        } else {
            stop("clusterType = '", clusterType, "' is not valid clustering version!")
        }
    }
        
    ## summary table of ScalarPredictors coefficients
    if (!is.null(ScalarPred)){
        coefs <- beta_full[ScPredStartIdx:length(beta_full)]
        coefsVar <- var_beta[ScPredStartIdx:length(beta_full)]
        coefm <- matrix(NA, ncol = 4 , nrow = length(coefs))
        coefm[,1] <- coefs
        coefm[,2] <- sqrt(coefsVar)
        coefm[,3] <- coefm[,1]/coefm[,2]
        coefm[,4] <- 1 - 2 * pnorm( abs(coefm[,3]))
        colnames(coefm) <- c("Estimate", "Std. Error", "t vlaue" , "Pr(>|t|)")
        rownames(coefm) <- colnames(Xvars)
    } else coefm <- NULL
    
    
    ## construct beta
    if (is.null(Curves)){
        beta <- NULL
    } else if (is.matrix(Curves)){
        beta <- list(slope = beta_full[1:p],
                     StdErr = sqrt(var_beta[1:p]))
    } else {
        beta <- lapply(1:length(Curves) , function(i){
            list(slope = beta_full[(i-1)*p + 1:p],
                 StdErr = sqrt(var_beta[(i-1)*p + 1:p]))
        })
        names(beta) <- if (!is.null(names(Curves))) names(Curves) else paste0("Curves.",1:length(Curves))
    }
    
    
    
    ## return summary, yHat and plot
    obj <- list(
        coefficients = list(coefm = coefm,
                            beta = beta),
        residuals = resid,
        fitted = yHat,
        model = list(effDf = effDf - if (!is.null(coefm)) dim(coefm)[1] else 0, rho = rho, intercept = intercept))
                                        #gcvVal = opt$value,
                                        #rho = optPar))
    class(obj) <- "SmSplFLR"
    obj
}

