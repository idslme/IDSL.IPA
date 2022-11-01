\name{IPA_baselineDeveloper}
\alias{IPA_baselineDeveloper}
\title{
Develop a baseline for the chromatogram using local minima
}
\description{
This function generates a vector of baselines for the chromatogram using local
minima. It also is capable of excluding outlier local minima to generate a
realistic baseline including true baseline regions. This baseline may represent
the local noise levels for the chromatogram.
}
\usage{
IPA_baselineDeveloper(segment, int)
}
\arguments{
  \item{segment}{
a matrix or a vecotr of adjusted scan number of local minima w/ or w/o redundant
local minima. Adjusted scan numbers are the scan numbers but adjusted to start
at 1.
}
  \item{int}{
a vector of intensities of the chromatogram.
}
}
\value{
A vector of baselines in the same size of the "int" vector.
}
\examples{
data(segment)
data(chromatogram_builder)
int <- chromatogram_builder$SmoothedChromatogram
IPA_baselineDeveloper(segment, int)
}
