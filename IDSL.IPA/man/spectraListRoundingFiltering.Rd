\name{spectraListRoundingFiltering}
\alias{spectraListRoundingFiltering}
\title{
spectraList filtering
}
\description{
This function reduces the size of the spectraList value by removing m/z values with no correspondence to 12C/13C isotopologue pairs.
}
\usage{
spectraListRoundingFiltering(specScanXIC, spectraList, roundingDigit = 1)
}
\arguments{
  \item{specScanXIC}{
a matrix of any size, but the first column containing the m/z of 12C
isotopologues are used.
}
  \item{spectraList}{
a list of mass spectra in each chromatogram scan.
}
  \item{roundingDigit}{
rounding digit to choose power of size reduction.
}
}
\value{
a list of mass spectra
}
