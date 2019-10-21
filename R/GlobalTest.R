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
#' @param type Which test to use. Currently, the following alternatives are available: "F" (see Kong et al. 2016)
#' 
#' @export
GlobalTest <- function(obj, alternative, type,...){

    if (!is.identical(class(obj), "FLR")){
        stop("obj must be of class FLR!")
    }

    
    
    
}
