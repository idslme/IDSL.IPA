\name{SNRxcms}
\alias{SNRxcms}
\title{
SNR xcms
}
\description{
This function calculates S/N values using a method suggested in the xcms paper (Tautenhahn, 2008).
}
\usage{
SNRxcms(int)
}
\arguments{
  \item{int}{
a vector of intensities corresponding to the vector of retention times for the
chromatographic peak.
}
}
\value{
S/N value
}
\references{
Tautenhahn, R., Böttcher, C. and Neumann, S. (2008). Highly sensitive feature detection for high resolution LC/MS. \emph{BMC bioinformatics}, 9(1), 1-16, \doi{10.1186/1471-2105-9-504}.
}
\examples{
data(peak_spline)
int <- peak_spline[, 2]
SNRxcms(int)
}
