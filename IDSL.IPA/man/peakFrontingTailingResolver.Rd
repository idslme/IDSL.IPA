\name{peakFrontingTailingResolver}
\alias{peakFrontingTailingResolver}
\title{
Fronting and tailing peaks resolver
}
\description{
This function attempts to resolve peak tailings or frontings into the main peak in case they were detected as seperate peaks.
}
\usage{
peakFrontingTailingResolver(segment, int, maxScanDifference, peakResolvingPower = 0.025)
}
\arguments{
  \item{segment}{
a matrix or a vector of peak boundaries.
}
  \item{int}{
a vector of intensities of the entire chromatogram.
}
  \item{maxScanDifference}{
maximum scan number difference between peak tailing or fronting and the main peak.
}
  \item{peakResolvingPower}{
power of peak resolving tool.
}
}
\value{
A matrix of 2 columns. Each row indicates peak boundary indices on the 'int' vector after resolving fronting and tailing peaks.
}
\examples{
data(segment)
data(chromatogramMatrix)
int <- chromatogramMatrix$smoothChromatogram
maxScanDifference <- 7
peakResolvingPower <- 0.2
peakFrontingTailingResolver(segment, int, maxScanDifference, peakResolvingPower)
}
