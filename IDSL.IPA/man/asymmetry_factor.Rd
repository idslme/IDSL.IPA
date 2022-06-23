\name{asymmetry_factor}
\alias{asymmetry_factor}
\title{
Asymmetry factor for a chromatographic peak
}
\description{
This function calculates an asymmetry factor for a chromatographic peak.
}
\usage{
asymmetry_factor(rt, int)
}
\arguments{
  \item{rt}{
a vector of retention times for the chromatographic peak.
}
  \item{int}{
a vector of intensities corresponding to the vector of retention times for the chromatographic peak.
}
}
\value{
asymmetry of the chromatographic peak. 1 is for very symmetric peak.
}
\examples{
data(peak_spline)
rt <- peak_spline[, 1]
int <- peak_spline[, 2] - peak_spline[, 3]
asymmetry_factor(rt, int)
}
