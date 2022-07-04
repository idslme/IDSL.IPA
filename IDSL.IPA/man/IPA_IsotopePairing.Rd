\name{IPA_IsotopePairing}
\alias{IPA_IsotopePairing}
\title{
IPA Isotope Pairing
}
\description{
This function pairs isotopologue in high-resolution mass spectral datasets
}
\usage{
IPA_IsotopePairing(spectraList, int_threshold, mass_accuracy_isotope_pair,
massDifferenceIsotopes)
}
\arguments{
  \item{spectraList}{
list of mass spectra in each chromatogram scan
}
  \item{int_threshold}{
intensity threshold at each chromatogram scan
}
  \item{mass_accuracy_isotope_pair}{
mass error to detect pair isotopologues
}
  \item{massDifferenceIsotopes}{
mass difference to pair isotopologues. (Default =  \eqn{\Delta C} = \out{<sup>13</sup>C} - \out{<sup>12</sup>C} = 1.003354835336), \eqn{\Delta S} = \out{<sup>34</sup>S} - \out{<sup>32</sup>S} = 1.9957958356, or any numerical value.
}
}
\value{
A matrix consists of 5 columns. The column contents are the m/z of \out{<sup>12</sup>C}
isotopologues, intensity of \out{<sup>12</sup>C} isotopologues, scan number (t), m/z of \out{<sup>13</sup>C}
isotopologues, and intensity of \out{<sup>13</sup>C} isotopologues, respectively.
}
