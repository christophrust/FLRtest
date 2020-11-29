#' Directed Test
#'
#' Test versus global alternative in Functional Linear Regression
#'
#'
#' @param obj A fitted model object of class "flm", as returned by
#' \code{EstFLM}
#'
#' @param null An optional global function to test against. The null
#' can be specified as a vector holding the discretized function object or,
#' alternatively, as an expansion of a given basis system. Defaults to
#' the zero function.
#'
#' @param conf.level Confidence level (1 - significance level) used to decide
#' whether to reject or not.
#'
#' @param ... Further arguments passed to underlying C-routine. See source code of
#' \code{testseq}.
#'
#' @details
#' \code{dirTest} performs the directed local test in the functional linear
#' regression model. This test can be used to discover a threshold t beyond which
#' the functional covariat have no significant effect on the outcome variable.
#'
#' @return
#' An object of type \code{flm.test} with S3-methods \code{print} and
#' \code{summary}.
#'
#' @references
#' Rust, C. (2020) Directed Local Testing in Functional Linear Models. Available
#' from the author upon request.
#'
#' @author Christoph Rust
#'
#' @seealso EstFLM, GlobalTest
#' @export
dirTest <- function(obj, null, conf.level = 0.95, ...){

    ## compute sequence of test statistics
    tseq <- testseq(obj, null,...)

    ## find rejections
    rej <- tseq[,"pval"] < (1 - conf.level)

    if (any(!rej)) {
        ## replace all TRUE after the first FALSE by FALSE
        rej[which(!rej)[1]:length(rej)] <- FALSE
    }


    ## return obj
    obj <- list(rejections = rej, test.sequence = tseq, conf.level = conf.level, obj = obj)

    class(obj) <- "flm.test"
    attr(obj, which = "test.type") <- "directed"

    obj
}

