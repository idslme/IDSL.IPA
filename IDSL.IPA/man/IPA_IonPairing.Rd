\name{IPA_IonPairing}
\alias{IPA_IonPairing}
\title{
IPA Ion Pairing
}
\description{
This function pairs two ions in high-resolution mass spectral datasets
}
\usage{
IPA_IonPairing(spectraList, int_threshold, mass_accuracy_isotope_pair,
massDifferenceIsotopes = 1.003354835336)
}
\arguments{
  \item{spectraList}{
list of mass spectra in each chromatogram scan
}
  \item{int_threshold}{
intensity threshold at each chromatogram scan
}
  \item{mass_accuracy_isotope_pair}{
mass error to detect pair ions
}
  \item{massDifferenceIsotopes}{
mass difference to pair ions. (Default =  \eqn{\Delta C} = \out{<sup>13</sup>C} - \out{<sup>12</sup>C} = 1.003354835336), \eqn{\Delta S} = \out{<sup>34</sup>S} - \out{<sup>32</sup>S} = 1.9957958356, or any numerical value.
}
}
\value{
A matrix consists of 5 columns. The column contents are the m/z of \out{<sup>12</sup>C}
isotopologues, intensity of \out{<sup>12</sup>C} isotopologues, scan number (t), m/z of \out{<sup>13</sup>C}
isotopologues, and intensity of \out{<sup>13</sup>C} isotopologues, respectively.
}
