\name{der_5points_stencil}
\alias{der_5points_stencil}
\title{
Numerical differentiation by five-point stencil method
}
\description{
This module performs numerical differentiation using the five-point stencil method.
}
\usage{
der_5points_stencil(x, y, n)
}
\arguments{
  \item{x}{
a vector of values for x.
}
  \item{y}{
a vector of values for y.
}
  \item{n}{
order of numerical differentiation (n=1-4).
}
}
\value{
A matrix of 2 columns. The first column represents x and the second column represents numerical differentiation values. This matrix has four rows (two rows from the beginning and 2 rows from the end) less than length of x or y.
}
\examples{
data(peak_spline)
rt <- peak_spline[, 1]
int <- peak_spline[, 2]
n <- 2 # second order derivative
der_5points_stencil(rt, int, n)
}
