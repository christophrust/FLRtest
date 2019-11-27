## print method for test object
#' @export
print.flm.test <- function(x,...){

    if (attr(x, which = "test.type") == "global"){
        cat("Test of association in the functional linear model\n")
        cat("--------------------------------------------------\n")
        cat(sprintf("Type: %s-test\n",x$type))
        cat(sprintf("Null hypothesis: %s\n",
                    if (is.null(x$alternative))
                        "no effect" else "specified alternative"))
        cat(sprintf("Test statistic: %.2f at %.2f and %.2f degrees of freedom\n",
                    x$statistic, x$df[1], x$df[2]))
        cat(sprintf("RSS0 = %f; RSS1 = %f\n", x$RSS0, x$RSS1))
        cat(sprintf("P-value: %f\n",x$pval))

    } else  if (attr(x, which = "test.type") == "directed"){

        pm <- max(which(x$rej))
        
        cat("Directed test of association in the functional linear model\n")
        cat("-----------------------------------------------------------\n")
        cat(sprintf(
            "Significant influence up to grid index %i.\nThis corresponds to %.2f in normalized domain\n", pm, seq(0,1, len = length(x$rej +1))[pm]
))
        cat(sprintf("Significance level used: %.2f\n", 1-x$conf.level))
    }
}



#' @export
summary.flm.test <- function(object,...){
    print(object)
}



## plot method for flm (estimation object)
#' @export
plot.flm <- function(x,...){
    cat("Yet to be implemented...\n")
}


#' @export
print.flm <- function(x,...){
    cat("Yet to be implemented...\n")
}


#' @export
summary.flm <- function(object,...){
    cat("Yet to be implemented...\n")
}
