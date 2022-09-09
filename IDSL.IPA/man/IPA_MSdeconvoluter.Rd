\name{IPA_MSdeconvoluter}
\alias{IPA_MSdeconvoluter}
\title{
MS deconvoluter
}
\description{
This function deconvolutes mass spectrometry files into a list of mass spectrals and a vector of retention times.
}
\usage{
IPA_MSdeconvoluter(HRMS_path, MSfile, MS_level = 1)
}
\arguments{
  \item{HRMS_path}{
address of the mass spectrometry file
}
  \item{MSfile}{
mass spectrometry file.
}
  \item{MS_level}{
MS level to extract information.
}
}
\value{
  \item{spectraList}{
 a list of mass spectra.
}
  \item{RetentionTime}{
 a vector of retention times for scan numbers.
}
  \item{MS_polarity}{
mass spectromtery ionization mode (+/-)
}
}
