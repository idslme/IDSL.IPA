\name{snr_rms}
\alias{snr_rms}
\title{
SNR RMS
}
\description{
This function calculates signal-to-noise ratio using root mean square.
}
\usage{
snr_rms(int, baseline, gauge)
}
\arguments{
  \item{int}{
 is the vector of intensities corresponding to the vector of retention times for the chromatographic peak.
}
  \item{baseline}{
is a vector of baseline of the chromatographic peak.
}
  \item{gauge}{
represents the gauge height of peak for gaussianity measurement.
}
}
\value{
S/N value
}
\examples{
data("peak_spline")
int <- peak_spline[, 2]
baseline <- peak_spline[, 3]
gauge <- 0.8
snr_rms(int, baseline, gauge)
}
