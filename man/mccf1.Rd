% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mccf1.R
\name{mccf1}
\alias{mccf1}
\title{Perform MCCF1 analysis}
\usage{
mccf1(response, predictor)
}
\arguments{
\item{response}{numeric vector representing ground truth classes (0 or 1).}

\item{predictor}{numeric vector representing prediction scores (in the range [0,1]).}
}
\value{
S3 object of class "mccf1", a list with the following members: `thresholds`: vector of
doubles describing the thresholds; `normalized_mcc`: vector of doubles representing normalized
MCC for each threshold; `f1`: vector of doubles representing F1 for each threshold.
}
\description{
`mccf1()` performs MCC (Matthews correlation coefficient)-F1 analysis for paired vectors
of binary response classes and fractional prediction scores representing the performance of
a binary classification task.
}
\examples{
response <- c(rep(1L, 1000L), rep(0L, 10000L))
set.seed(2017)
predictor <- c(rbeta(300L, 12, 2), rbeta(700L, 3, 4), rbeta(10000L, 2, 3))
x <- mccf1(response, predictor)
head(x$thresholds)
# [1]  Inf 0.9935354 0.9931493 0.9930786 0.9925507 0.9900520
head(x$normalized_mcc)
# [1]  NaN 0.5150763 0.5213220 0.5261152 0.5301566 0.5337177
head(x$f1)
# [1]  NaN 0.001998002 0.003992016 0.005982054 0.007968127 0.009950249
}
