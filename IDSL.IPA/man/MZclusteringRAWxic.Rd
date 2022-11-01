\name{MZclusteringRAWxic}
\alias{MZclusteringRAWxic}
\title{
m/z clustering XIC
}
\description{
This function clusters related 12C m/z values.
}
\usage{
MZclusteringRAWxic(spectraScan, massAccuracyXIC, minPeakHeight, minNIonPair)
}
\arguments{
  \item{spectraScan}{
a matrix consists of 3 columns. The column contents are the m/z of 12C
isotopologues, intensity of 12C isotopologues, and scan number (t).
}
  \item{massAccuracyXIC}{
mass accuracy to detect related 12C m/z values.
}
  \item{minNIonPair}{
minimum number of nIsoPair for an individual peak.
}
  \item{minPeakHeight}{
minimum peak height for an individual peak.
}
}
\value{
This function returns an list on index numbers of EICs for the "spectraScan" variable.
}
