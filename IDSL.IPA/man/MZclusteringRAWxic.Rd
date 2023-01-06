\name{mzClusteringRawXIC}
\alias{mzClusteringRawXIC}
\title{
m/z clustering raw XIC
}
\description{
This function clusters related 12C m/z values.
}
\usage{
mzClusteringRawXIC(spectraScan123, massAccuracy, minNIonPair, minPeakHeightXIC)
}
\arguments{
  \item{spectraScan123}{
a matrix consists of 3 columns. The column contents are the m/z of 12C isotopologues, intensity of 12C isotopologues, and scan number (t).
}
  \item{massAccuracy}{
mass accuracy to detect related 12C m/z values.
}
  \item{minNIonPair}{
minimum number of nIsoPair for an individual peak.
}
  \item{minPeakHeightXIC}{
minimum peak height for an individual raw EIC
}
}
\value{
This function returns an list on index numbers of EICs for the "spectraScan" variable.
}
