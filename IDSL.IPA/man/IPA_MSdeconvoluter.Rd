\name{IPA_MSdeconvoluter}
\alias{IPA_MSdeconvoluter}
\title{
MS deconvoluter
}
\description{
This function deconvolutes mass spectrometry files into a list of mass spectra and a vector of retention times.
}
\usage{
IPA_MSdeconvoluter(inputHRMSfolderPath, MSfileName, MSlevel = 1)
}
\arguments{
  \item{inputHRMSfolderPath}{
address of the mass spectrometry file
}
  \item{MSfileName}{
mass spectrometry file.
}
  \item{MSlevel}{
MS level to extract information.
}
}
\value{
  \item{spectraList}{
 a list of mass spectra.
}
  \item{retentionTime}{
 a vector of retention times for scan numbers.
}
  \item{MS_polarity}{
mass spectrometry ionization mode (+/-)
}
}
