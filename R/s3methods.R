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

      if (any(x$rej)){
        pm <- max(which(x$rej))

        cat("Directed test of association in the functional linear model\n")
        cat("-----------------------------------------------------------\n")
        cat(sprintf(
          "Significant influence up to grid index %i.\nThis corresponds to %.2f in normalized domain\n",
          pm, seq(0,1, len = length(x$rej +1))[pm]
        ))
      } else {
          cat("Directed test of association in the functional linear model\n")
        cat("-----------------------------------------------------------\n")
        cat(sprintf(
          "The test did not reject at all!\n"))
      }

        cat(sprintf("Significance level used: %.2f\n", 1-x$conf.level))
    }
}


#' @export
plot.flm.test <- function(x,...){

    p <- length(x$rejections) + 1
    X <- x$obj$data$X
    y <- x$obj$data$y

    grd <- seq(0,1, length.out = p)

    if (any(x$rejections)) {

        pm <- max(which(x$rejections))

        df <- max(x$obj$model$spline$basis,8)

        basis <- .Call("R_SplitSplineBasis",
                       grd = grd,
                       df = as.integer(df) ,
                       splitpoint = grd[pm])$basis
        pXB <- 1/dim(X)[2] * X %*% basis

        betahat <- as.vector(basis %*% (solve(crossprod(pXB)) %*% t(pXB)) %*% y)

        plot(x = NULL, y = NULL, xlim = c(0,1), ylim = range(betahat),
             xlab = "Normalized domain", ylab = expression(paste(beta, "-function")))
        polygon(x = c(0, grd[pm],grd[pm],0), y = c(rep(par("usr")[3:4], each = 2)),
                col = "lightgrey")
        abline(h = 0, lty = 3)
        lines(x = grd[1:pm], y = betahat[1:pm], col = "#444444", lwd = 4)
        lines(x = grd[(pm+1):p], y = betahat[(pm+1):p], col = "#d2d2d2", lwd = 2, lty = 2)

    } else {


        betahat <- x$obj$coefficients$beta

        plot(x = grd, y = betahat, lty = 2, col = "#b2b2b2", lwd = 2, type = "l",
             xlab = "Normalized domain", ylab = expression(paste(beta, "-function")))
    }
}


#' @export
summary.flm.test <- function(object,...){
    print(object)
}



## plot method for flm (estimation object)
#' @export
plot.flm <- function(x,...){

    cat("Plotting method yet not fully implemented...\n")
    graphics::plot(x = seq(0,1,length.out = length(x$coefficients$beta)),
         y = x$coefficients$beta, type = "l")
}


#' @export
print.flm <- function(x,...){
    cat("Yet to be implemented...\n")
}


#' @export
summary.flm <- function(object,...){
    cat("Yet to be implemented...\n")
}
