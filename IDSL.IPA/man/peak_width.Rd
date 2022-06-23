\name{peak_width}
\alias{peak_width}
\title{
peak width measuement
}
\description{
This function measures peak width at different peak heights.
}
\usage{
peak_width(rt, int, gauge)
}
\arguments{
  \item{rt}{
a vector of retention times of the chromatographic peak.
}
  \item{int}{
a vector of intensities of the chromatographic peak.
}
  \item{gauge}{
a height gauge to measure the peak width. This parameter should be between 0-1.
}
}
\value{
A peak width at the guaged height.
}
\examples{
data("peak_spline")
rt <- peak_spline[, 1]
int <- peak_spline[, 2] - peak_spline[, 3]
gauge <- 0.5
peak_width(rt, int, gauge)
}
