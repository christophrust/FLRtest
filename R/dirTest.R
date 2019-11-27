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
#' @param ... Further arguments passed to underlying C-routine. See source code.
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
dirTest <- function(obj, null, ...){

    tseq <- testseq(obj, null,...)
    
}

