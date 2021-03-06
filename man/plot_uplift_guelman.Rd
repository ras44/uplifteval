% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_uplift_guelman.R
\name{plot_uplift_guelman}
\alias{plot_uplift_guelman}
\title{A direct copy of Leo Guelman's uplift::qini function
available in the R uplift package at commit
95965272e71c312623c95c439fb0b84f95c185b7:
https://github.com/cran/uplift/blob/95965272e71c312623c95c439fb0b84f95c185b7/R/qini.R#L5}
\usage{
plot_uplift_guelman(p1, W, Y, p0 = rep(0, length(p1)), plotit = TRUE,
  direction = 1, groups = 10)
}
\arguments{
\item{p1}{vector of numeric uplift predictions.  Some uplift models produce two predictions: if-treated and if-control.  In this case, if-treated predictions can be provided as p1, and if-control predictions can be provided as p0.}

\item{W}{vector of {0,1} treatment indicators}

\item{Y}{vector of {0,1} outcomes}

\item{p0}{vector of numeric control predictions (default 0)}

\item{plotit}{boolean plot the Qini chart}

\item{direction}{1: calculate the differential response as p1-p0, 2: p0-p1}

\item{groups}{5, 10, or 20: the number of quantiles in which to divide the population}
}
\description{
A direct copy of Leo Guelman's uplift::qini function
available in the R uplift package at commit
95965272e71c312623c95c439fb0b84f95c185b7:
https://github.com/cran/uplift/blob/95965272e71c312623c95c439fb0b84f95c185b7/R/qini.R#L5
}
\examples{

set.seed(0)
rl <- function(x){
  round(1/(1+exp(-x)))
}
n <- 2000; p <- 3
beta <- -0.5
X <- matrix(rnorm(n*p), n, p)
W <- rbinom(n, 1, 0.5)
Y <- rl(pmax(beta+X[,1], 0) * W + X[,2])
p1 <- 1/(1+exp(-(beta+X[,1])))
plot_uplift_guelman(p1, W, Y, groups=10, plotit=TRUE)

\donttest{
library(grf)
set.seed(123)

alpha <- 0.1
n <- 1000
W <- rbinom(n, 1, 0.5)
Y <- W
p1 <- Y + alpha*rnorm(n)
plot_uplift_guelman(p1, W, Y, groups=10)


rl <- function(x){
  round(1/(1+exp(-x)))
}
n <- 2000; p = 10
X <- matrix(rnorm(n*p), n, p)
W <- rbinom(n, 1, 0.2)
Y <- rl(rl(X[,1]) * W - rl(X[,3]) * W + rnorm(n))
tau.forest <- causal_forest(X, Y, W)
tau.hat <- predict(tau.forest, X)
plot_uplift_guelman(tau.hat$predictions, W, Y)

}
}
