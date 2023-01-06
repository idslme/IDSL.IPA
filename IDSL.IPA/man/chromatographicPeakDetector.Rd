\name{chromatographicPeakDetector}
\alias{chromatographicPeakDetector}
\title{
peak detection
}
\description{
This function detects separated chromatographic peaks on the chromatogram.
}
\usage{
chromatographicPeakDetector(int)
}
\arguments{
  \item{int}{
a vector of intensities of the chromatogram.
}
}
\value{
A matrix of 2 columns. Each row indicates peak boundary indices on the 'int' vector.
}
\examples{
data(chromatogramMatrix)
int <- chromatogramMatrix$smoothChromatogram
chromatographicPeakDetector(int)
}
