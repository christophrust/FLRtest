# Directed Local Testing in the Functional Linear Model

This R package implements the testing procedure described in the manuscript 
[Directed Local Testing in the Functional Linear Model](https://christophrust.de/manuscripts/DLTFLR.pdf).

## Installation

```splus
devtools::install_github("https://github.com/christophrust/FLRTest.git")
```


## Example

```splus
set.seed(123)

## number of discretization points
p <- 100

## number of observations
N <- 1000

## draw some curves
grd <- seq(0,1,length = p)
X <- t(bs(x= grd, df = p, intercept = TRUE) %*% matrix(rnorm(N*p), ncol = N))

## generate data via the functional linear model
beta <- pmax(sin(5*grd), 0)
y <- X %*% beta/p + rnorm(N, sd=0.1)

## estimate a FLM and perform the test
est  <- EstFLM(y, X, type = "spline", df=20, intercept = TRUE)
test <- dirTest(est)

## plot results + true functional coefficient
plot(test)
lines(grd, beta, lty = 2, lwd = 2)
```
