\name{SNRbaseline}
\alias{SNRbaseline}
\title{
SNR baseline
}
\description{
This function calculates S/N using local noise levels from baseline,
}
\usage{
SNRbaseline(int, baseline)
}
\arguments{
  \item{int}{
a vector of intensities corresponding to the vector of retention times for the chromatographic peak.
}
  \item{baseline}{
a vector of baseline of the chromatographic peak.
}
}
\value{
S/N value
}
\examples{
data("peak_spline")
int <- peak_spline[, 2]
baseline <- peak_spline[, 3]
SNRbaseline(int, baseline)
}
