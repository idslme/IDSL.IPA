\name{pseudomoments_symmetry}
\alias{pseudomoments_symmetry}
\title{
pseudomoments symmetry
}
\description{
This function measures peak symmetry and skewness using the inflection points
of the peak on both sides.
}
\usage{
pseudomoments_symmetry(rt, int)
}
\arguments{
  \item{rt}{
 a vector of retention times for the chromatographic peak.
}
  \item{int}{
 a vector of intensities corresponding to the vector of retention times for
 the chromatographic peak.
}
}
\value{
  \item{PeakSymmetry}{
 peak symmetry for the chromatographic peak.
}
  \item{Skewness}{
 skewness for the chromatographic peak.
}
}
\examples{
data("peak_spline")
rt <- peak_spline[, 1]
int <- peak_spline[, 2] - peak_spline[, 3]
pseudomoments_symmetry(rt, int)
}
