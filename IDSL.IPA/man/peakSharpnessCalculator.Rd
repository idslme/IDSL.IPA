\name{peakSharpnessCalculator}
\alias{peakSharpnessCalculator}
\title{
Peak Sharpness Calculator
}
\description{
This function measures sharpness of a chromatographic peak
}
\usage{
peakSharpnessCalculator(int)
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
peakSharpnessCalculator(int)
}
