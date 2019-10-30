#' SeqTest
#' 
#' Sequential test in the functional linear model
#'
#' 
#' @param obj A fitted model object of class "flm", as returned by
#' \code{EstFLM}
#' 
#' @param startval A value in the domain \eqn{[0,1]} where it is
#' assumed that the effect is strongest. Defaults to 0.
#'
#' @param direction Either of \code{"right","left","both"}, specifying
#' into which direction to test starting from \code{startvalue}.
#' Defaults to \code{"right"}. 
#'
#' @param null An optional global function to test against. The null
#' can be specified as a vector holding the discretized function object or,
#' alternatively, as an expansion of a given basis system. Defaults to
#' the zero function.
#'
#' @details
#' \code{SeqTest} performs a sequential test of no association in the
#' functional linear model in order to detect regions where the functions
#' have a significant influence on the scalar outcome. The procedure can
#' also be used to test whether the estimate differs from a specific
#' functional slope parameter.
#'
#' @return An object of type \code{flm.dirtest} with S3-methods \code{print} and
#' \code{summary}.
#'
#' @references
#' Rust, C. (2019) Directed Local Testing in the Functional Linear Model.
#' Author Manuscript.
#' 
#' 
#' @author Christoph Rust
#'
#' @seealso EstFLM, GlobalTest
#' @export
SeqTest <- function(obj, startval = 0, direction = "right", null){

    if (!identical(class(obj), "flm")){
        stop("obj must be of class flm!")
    }
    
    p <- length(obj$coefficients$beta)
    Nobs <- length(obj$fitted)

    if (!missing(null)) {
        ## check for correctly specified null
        if (is.numeric(null) && (length(null) != p)){

            stop(gettext("Wrongly specified null!\nIf specified as numeric vector, it needs to be of same length as the number of discretization points used in '%s': %i", obj, p))
            
        } else if (is.list(null) &&
                   (!all(names(null) %in% c("basis","coefs")) ||
                    (dim(null$basis)[1] != p ||
                     dim(null$basis)[2] != length(null$coefs)))){

            stop("Wrongly specified null!\nIf specified as basis expansion 'null' must contain both a basis (matric of dimension p times k) and coefs (vector of length k)!")
            
        } else if (!is.list(null) && !is.numeric(null)){
            stop("Wrongly specified null! Please consult the help file.")
        }
    }
    
    ## if (missing(type) || type == "F"){
    ##     ## RSS of full model
    ##     RSSfull <- sum(obj$residuals^2)
        
    ##     ## RSS of null model
    ##     beta0 <- if (missing(null)){
    ##                  NULL
    ##              } else if (is.numeric(null)) {
    ##                  matrix(null, ncol = 1)
    ##              } else {
    ##                  null$basis %*% null$coefs
    ##              }
        
    ##     RSSnull <- if (missing(null)){
    ##                    sum( (obj$data$y - if(intercept) mean(obj$data$y) else 0)^2)
    ##                } else {
    ##                    sum( (obj$data$y -  obj$data$X[,1:p, drop=FALSE] %*% beta0/p -
    ##                          if(intercept) mean(obj$data$y) else 0)^2)
    ##                }
        

    ##     ## Test statistic
    ##     Tf <- ( (RSSnull - RSSfull) / (df1 <- obj$model$effDf + obj$model$intercept - intercept) ) /
    ##         (RSSfull / (df2 <- Nobs - obj$model$effDf - obj$model$intercept) )


    ##     ## Finite sample null distribution, see Remark 1 in Kong et al 2017 (JNonpStat)
    ##     ## if spline basis, replace s_n by effDf

    ##     if (missing(null)){
    ##         pval <- 1 - pf(Tf, df1 = df1, df2 = df2)
    ##     } else {
    ##         ncp <- 1/p^2 * (t(beta0) %*% cov(obj$data$X) %*% beta0)
    ##         pval <- 1 - pf(Tf, df1 = df1, df2 = df2, ncp = ncp)
    ##     }
    ## } else {
    ##     stop("Currently only the F-test is supported!")
    ## }
    
    
    ## collect return
    obj <- list(statistic = Tf,
                df = c(df1,df2),
                pval = pval,
                RSS0 = RSSnull,
                RSS1 = RSSfull,
                beta = obj$coefficients$beta,
                type = if (missing(type)) "F" else type,
                null = beta0)
    class(obj) <- "flm.test"
    obj
}
