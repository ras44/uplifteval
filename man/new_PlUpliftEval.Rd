% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plUpliftEval.R
\name{new_PlUpliftEval}
\alias{new_PlUpliftEval}
\title{A non-S3 "constructor" function that returns a list representing a pylift uplift
eval object, PlUpliftEval.  This object contains metrics and can be used to generate plots.}
\usage{
new_PlUpliftEval(treatment = integer(), outcome = integer(),
  prediction = numeric(), p = "infer", n_bins = 20)
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
A non-S3 "constructor" function that returns a list representing a pylift uplift
eval object, PlUpliftEval.  This object contains metrics and can be used to generate plots.
}
