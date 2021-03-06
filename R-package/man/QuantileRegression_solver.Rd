% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QuantileRegression_solver.R
\name{QuantileRegression_solver}
\alias{QuantileRegression_solver}
\title{Solve given quantile regression problem in parametric simplex method}
\usage{
QuantileRegression_solver(X, y, max_it = 50, lambda_threshold = 0.01,
  tau = 0.5)
}
\arguments{
\item{X}{\code{x} is an \code{n} by \code{d} data matrix}

\item{y}{\code{y} is a length \code{n} response vector}

\item{max_it}{This is the number of the maximum path length one would like to achieve. The default length is \code{50}.}

\item{lambda_threshold}{The parametric simplex method will stop when the calculated parameter is smaller than lambda. The default value is \code{0.01}.}

\item{tau}{The quantile number you want. The default quantile is \code{0.5}}
}
\description{
Solve given quantile regression problem in parametric simplex method
}
\note{
The returned value is similar with \code{Dantzig_solver()}
}
\seealso{
\code{\link{primal-package}}, \code{\link{Dantzig_solver}}
}
