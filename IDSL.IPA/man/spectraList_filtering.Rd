\name{spectraList_filtering}
\alias{spectraList_filtering}
\title{
spectraList filtering
}
\description{
This function reduces the size of the spectraList value by removing m/z values
with no correspondence to 12C/13C isotopologue pairs.
}
\usage{
spectraList_filtering(spec_scan.xic, spectraList, rounding_digit)
}
\arguments{
  \item{spec_scan.xic}{
a matrix of any size, but the first column containing the m/z of 12C
isotopologues are used.
}
  \item{spectraList}{
a list of mass spectra in each chromatogram scan.
}
  \item{rounding_digit}{
rounding digit to choose power of size reduction.
}
}
\value{
a list of mass spectrals
}
