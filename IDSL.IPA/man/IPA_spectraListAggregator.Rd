\name{IPA_spectraListAggregator}
\alias{IPA_spectraListAggregator}
\title{
spectraList filtering
}
\description{
This module stacks the spectraList object and creates a list of ions for a rapid spectra query.
}
\usage{
IPA_spectraListAggregator(spectraList)
}
\arguments{
  \item{spectraList}{
a list of mass spectra in each chromatogram scan.
}
}
\value{
  \item{aggregatedSpectraList}{
aggregated spectraList
}
  \item{spectraListMatrix}{
matrix of row bounded spectraList
}
}
