\name{XIC}
\alias{XIC}
\title{
XIC
}
\description{
XIC
}
\usage{
XIC(spectraList.xic, scan_number_start = 1, mz_target, mass_accuracy_xic)
}
\arguments{
  \item{spectraList.xic}{
a list of mass spectra in each chromatogram scan.
}
  \item{scan_number_start}{
the first scan number.
}
  \item{mz_target}{
an m/z value to perform XIC analysis.
}
  \item{mass_accuracy_xic}{
a mass error to perform XIC analysis.
}
}
\value{
A matrix of three columns represpenting scan number, m/z, and intensity.
}
