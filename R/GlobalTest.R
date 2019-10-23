#' GlobalTest
#' 
#' Test versus global alternative in Functional Linear Regression
#'
#' 
#' @param obj A fitted model object, see also estFLR
#' 
#' @param alternative The global alternative to test against. The alternative
#' can be specified as a vector holding the discretized function object or,
#' alternatively, as an expansion of a given basis system. Defaults to the zero function.
#'
#' @param type Which test to use. Currently, the following alternatives are available: "F" (see Kong et al. 2017)
#'
#' @param intercept TRUE or FALSE, indicating whether intercept to be included in alternative.
#' 
#' @export
GlobalTest <- function(obj, alternative, type, intercept = TRUE){

    if (!identical(class(obj), "flm")){
        stop("obj must be of class flm!")
    }

    p <- length(obj$coefficients$beta)
    Nobs <- length(obj$fitted)

    if (!missing(alternative)) {
        ## check for correctly specified alternative
        if (is.numeric(alternative) && (length(alternative) != p)){

            stop(gettext("Wrongly specified alternative!\nIf specified as numeric vector, it needs to be of same length as the number of discretization points used in '%s': %i", obj, p))
            
        } else if (is.list(alternative) &&
                   (!all(names(alternative) %in% c("basis","coefs")) ||
                    (dim(alternative$basis)[1] != p ||
                     dim(alternative$basis)[2] != length(alternative$coefs)))){

            stop("Wrongly specified alternative!\nIf specified as basis expansion 'alternative' must contain both a basis (matric of dimension p times k) and coefs (vector of length k)!")
            
        } else if (!is.list(alternative) && !is.numeric(alternative)){
            stop("Wrongly specified alternative! Please consult the help file.")
        }
    }
    
    if (missing(type) || type == "F"){
        ## RSS of full model
        RSSfull <- sum(obj$residuals^2)
    
        ## RSS of null model
        beta0 <- if (missing(alternative)){
                     NULL
                 } else if (is.numeric(alternative)) {
                     alternative
                 } else {
                     alternative$basis %*% alternative$coefs
                 }

        RSSnull <- if (missing(alternative)){
                       sum( (obj$data$y - if(intercept) mean(obj$data$y) else 0)^2)
                   } else {
                       sum( (obj$data$y -  obj$data$X[,1:p, drop=FALSE] %*% beta0/p -
                             if(intercept) mean(obj$data$y) else 0)^2)
                   }


        ## Test statistic
        Tf <- ( (RSSnull - RSSfull) / (df1 <- obj$model$effDf + obj$model$intercept - intercept) ) /
            (RSSfull / (df2 <- Nobs - obj$model$effDf - obj$model$intercept) )


        ## Finite sample null distribution, see Remark 1 in Kong et al 2017 (JNonpStat)
        ## if spline basis, replace s_n by effDf

        if (missing(alternative)){
            pval <- 1 - pf(Tf, df1 = df1, df2 = df2)
        } else {
            ncp <- 1/p^2 * beta0 %*% cov(obj$data$X) %*% beta0
            pval <- 1 - pf(Tf, df1 = df1, df2 = df2, ncp = ncp)
        }
    } else {
        stop("Currently only the F test is supported!")
    }

    
    ## collect return
    obj <- list(stat = Tf,
                pval = pval,
                RSS0 = RSSnull,
                RSS1 = RSSfull,
                beta = obj$coefficients$beta,
                alternative = beta0)
    class(obj) <- "flm.test"
    obj
    
}
