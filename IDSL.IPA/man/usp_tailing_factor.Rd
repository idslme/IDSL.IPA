\name{usp_tailing_factor}
\alias{usp_tailing_factor}
\title{
USP tailing factor
}
\description{
This function calculates USP tailing factor at above 10 percent of the height.
}
\usage{
usp_tailing_factor(rt, int)
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
USP tailing factor for the chromatographic peak.
}
\examples{
data(peak_spline)
rt <- peak_spline[, 1]
int <- peak_spline[, 2] - peak_spline[, 3]
usp_tailing_factor(rt, int)
}
