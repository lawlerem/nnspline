% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nnodes.R
\name{symsqrt}
\alias{symsqrt}
\title{Compute symmetric square roots of positive definite matrices}
\usage{
symsqrt(m, invert = TRUE)
}
\arguments{
\item{m}{A positive definite matrix}

\item{invert}{If TRUE, the take the symmetric square root of m^-1}
}
\value{
A symmetric matrix x where x \%*\% x == m (or m^-1)
}
\description{
Compute symmetric square roots of positive definite matrices
}
