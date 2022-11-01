\name{peakAreaCalculator}
\alias{peakAreaCalculator}
\title{
peak area
}
\description{
This function calculates area under the curve using a trapezoid method.
}
\usage{
peakAreaCalculator(x, y)
}
\arguments{
  \item{x}{
is a vector of x values.
}
  \item{y}{
is a vector of y values.
}
}
\value{
A number for the integrated peak area.
}
\examples{
data("peak_spline")
rt <- peak_spline[, 1]
int <- peak_spline[, 2]
peakAreaCalculator(rt, int)
}
