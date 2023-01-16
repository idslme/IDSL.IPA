\name{IPA_IonPairing}
\alias{IPA_IonPairing}
\title{
IPA Ion Pairing
}
\description{
This function pairs two ions with a fixed distance in high-resolution mass spectral datasets
}
\usage{
IPA_IonPairing(spectraList, minSpectraNoiseLevel, massAccuracyIonPair = 0.015,
ionMassDifference = 1.003354835336)
}
\arguments{
  \item{spectraList}{
list of mass spectra in each chromatogram scan
}
  \item{minSpectraNoiseLevel}{
intensity threshold at each chromatogram scan
}
  \item{massAccuracyIonPair}{
mass error to detect pair ions
}
  \item{ionMassDifference}{
mass difference to pair ions. (Default = DeltaC = 13C - 12C = 1.003354835336), or DeltaS = 34S - 32S = 1.9957958356, or any numerical value.
}
}
\value{
A matrix consists of 5 columns. The column contents are the m/z of 12C isotopologues, intensity of 12C isotopologues, scan number (t), m/z of 13C isotopologues, and intensity of 13C isotopologues, respectively.
}
