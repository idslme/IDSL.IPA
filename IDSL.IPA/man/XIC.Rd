\name{XIC}
\alias{XIC}
\title{
XIC
}
\description{
XIC
}
\usage{
XIC(spectraListXIC, scanNumberStart = 1, mzTarget, massAccuracyXIC)
}
\arguments{
  \item{spectraListXIC}{
a list of mass spectra in each chromatogram scan.
}
  \item{scanNumberStart}{
the first scan number.
}
  \item{mzTarget}{
an m/z value to perform XIC analysis.
}
  \item{massAccuracyXIC}{
a mass error to perform XIC analysis.
}
}
\value{
A matrix of three columns representing scan number, m/z, and intensity.
}
