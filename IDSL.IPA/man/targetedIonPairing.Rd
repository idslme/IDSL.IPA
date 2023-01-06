\name{targetedIonPairing}
\alias{targetedIonPairing}
\title{
Targeted Ion Pairing
}
\description{
This module only pairs `mzTarget` values across `scanNumberStart` through `scanNumberEnd` scan numbers.
}
\usage{
targetedIonPairing(spectraList, scanNumberStart, scanNumberEnd, mzTarget,
massAccuracy, ionMassDifference = 1.003354835336, massAccuracyIonPair = massAccuracy*1.5)
}
\arguments{
  \item{spectraList}{
spectraList which is a list of mass spectra
}
  \item{scanNumberStart}{
the first scan number.
}
  \item{scanNumberEnd}{
the last scan number.
}
  \item{mzTarget}{
m/z value to perform chromatography analysis
}
  \item{massAccuracy}{
mass accuracy to select the dominant ion
}
  \item{ionMassDifference}{
mass difference to pair ions. (Default = DeltaC = 13C - 12C = 1.003354835336), or DeltaS = 34S - 32S = 1.9957958356, or any numerical value.
}
  \item{massAccuracyIonPair}{
mass accuracy to select the second ion
}
}
\value{
A targeted ion paired spectra and their scan numbers
}
