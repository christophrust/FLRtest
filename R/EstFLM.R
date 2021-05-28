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
#' @param intercept Logical, if TRUE, an intercept is included in the model.
#'
#' @param type A character specifying, which basis to use. Possible options:
#' "spline" (default), "smoothspline" or "fpc".
#'
#' @param df Degrees of freedom of the model fit. This is either the cardinality of
#' the underlying basis or it is the edf computed as the trace of the smoother matrix
#' in the case of smoothing spline estimation.
#'
#' @param rho If df is not specified, in the smoothing spline case, also the smoothing
#' parameter rho can control model flexibility.
#'
#' @param cpv For the fpc basis, the cumulative proportion of the explained variance
#' can be taken to select the number of basis functions.
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
EstFLM <- function(y, X, intercept = TRUE, type = "spline", df = NULL, rho = NULL, cpv=NULL){


    if (missing(X)) stop("No functional predictors specified.")
    if (missing(y)) stop("No dependent variable specified.")
    if (missing(df) & identical(type, "spline"))
        stop("For a predefined spline basis, please define df")

    Nobs <- length(y)
    p <- dim(X)[2]
    grd <- seq(0, 1, length.out = p)


    if (type == "smoothspline") {

        ## basis and Amat
        splMat <- natSplBasis(grd)


        if (intercept){
            X <- cbind(p,X)
            splMat$A_m <- rbind( 0, cbind(0, splMat$A_m))
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

            rho <- exp(stats::uniroot( function(x) {dfGivenRho(x) - df},
                               lower  = -100, upper = 100, f.lower = p,
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
        Cov <- stats::cov(X)
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
        lmEst <- stats::lm.fit(x = if(intercept) cbind(1,scores) else scores,
                         y = y)
        xi <- lmEst$coefficients
        yHat <- (if (intercept) cbind(1,scores) else scores) %*% xi

        ## recover curve
        beta <- c( if(intercept) xi[1] else NULL, efuncs %*% (if(intercept) xi[-1] else xi))

        effDf <- ncol(efuncs)


    } else if (type == "spline") {

        ## create basis
        knots <- c(0,0,0,seq(0,1, length.out = df - 2),1,1,1)
        basis <- splines::splineDesign(knots = knots, x = grd)

        ## product of X and basis
        pXBasis <- 1/p * X %*% basis

        ## include intercept
        if (intercept) pXBasis <- cbind(1,pXBasis)

        ## compute least square estimate
        res <- stats::lm.fit(x = pXBasis, y = y)

        ## store resulst
        beta <- if (intercept) {
                    c(res$coefficients[1], basis %*% res$coefficients[-1])
                } else {
                    basis %*% res$coefficients
                }
        effDf <- length(res$coefficients)
        yHat <- res$fitted

    } else{
        stop(sprintf("type '%s' not a valid estimation specification!", type))
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
                                                            A_m = splMat$A_m) else NULL,
                     spline = if (type == "spline") list(basis = basis,
                                                         lmfit = res) else NULL
                     ),
        data = list(y=y, X=if (intercept && type == "smoothspline") X[,-1, drop = FALSE] else X)
    )

    class(obj) <- "flm"
    obj
}
