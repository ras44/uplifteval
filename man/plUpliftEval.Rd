% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plUpliftEval.R
\name{plUpliftEval}
\alias{plUpliftEval}
\title{A helper for the new_PlUpliftEval function that validates the treatment,
outcome, prediction, p, and n_bins arguments.}
\usage{
plUpliftEval(treatment, outcome, prediction, p = "infer", n_bins = 20)
}
\arguments{
\item{treatment}{numeric vector of treatment identifiers}

\item{outcome}{numeric vector of outcomes}

\item{prediction}{numeric vector of uplift predictions}

\item{p}{optional "infer", numeric, numeric vector representing treatment
propensities}

\item{n_bins}{integer number of bins on x-axis; default 20}
}
\value{
a list representing a pylift uplift eval object
}
\description{
A helper for the new_PlUpliftEval function that validates the treatment,
outcome, prediction, p, and n_bins arguments.
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
plUpliftEval(W, Y, p1)

\donttest{
library(grf)
set.seed(123)

rl <- function(x){
  round(1/(1+exp(-x)))
}
n <- 2000; p <- 10
X <- matrix(rnorm(n*p), n, p)
W <- rbinom(n, 1, 0.2)
Y <- rl(rl(X[,1]) * W - rl(X[,3]) * W + rnorm(n))
tau.forest <- causal_forest(X, Y, W)
tau.hat <- predict(tau.forest, X)
plue <- plUpliftEval(W, Y, tau.hat$predictions)
plue
}

}
