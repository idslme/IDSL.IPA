\name{peakSharpness}
\alias{peakSharpness}
\title{
Peak Sharpness
}
\description{
This function measures sharpness of a chromatographic peak
}
\usage{
peakSharpness(int)
}
\arguments{
  \item{int}{
a vector of intensities of the chromatographic peak.
}
}
\value{
A number representing peak sharpness. The higher values indicate higher
sharpness.
}
\examples{
data("peak_spline")
int <- peak_spline[, 2]
peakSharpness(int)
}
