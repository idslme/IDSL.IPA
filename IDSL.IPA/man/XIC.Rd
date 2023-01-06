\name{XIC}
\alias{XIC}
\title{
XIC
}
\description{
XIC
}
\usage{
XIC(aggregatedSpectraList, scanNumberStart, scanNumberEnd, mzTarget, massAccuracy)
}
\arguments{
  \item{aggregatedSpectraList}{
aggregated spectraList and spectra matrix from the `IPA_spectraListAggregator` module
}
  \item{scanNumberStart}{
the first scan number.
}
  \item{scanNumberEnd}{
the last scan number.
}
  \item{mzTarget}{
an m/z value to perform XIC analysis.
}
  \item{massAccuracy}{
a mass error to perform XIC analysis.
}
}
\value{
A matrix of three columns representing scan number, m/z, and intensity.
}
