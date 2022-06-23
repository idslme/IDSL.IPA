\name{snr_signal2baseline}
\alias{snr_signal2baseline}
\title{
SNR baseline
}
\description{
This function calculates S/N using local noise levels from baseline,
}
\usage{
snr_signal2baseline(int, baseline)
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
snr_signal2baseline(int, baseline)
}
