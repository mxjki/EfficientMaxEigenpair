% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tridiag.R
\name{tridiag}
\alias{tridiag}
\title{Tridiagonal matrix}
\usage{
tridiag(upper, lower, main)
}
\arguments{
\item{upper}{The upper diagonal vector.}

\item{lower}{The lower diagonal vector.}

\item{main}{The main diagonal vector.}
}
\value{
A tridiagonal matrix is returned.
}
\description{
Generate tridiagonal matrix Q based on three input vectors.
}
\examples{
a = c(1:7)^2
b = c(1:7)^2
c = -c(1:8)^2
tridiag(b, a, c)
}
