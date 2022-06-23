\name{chromatogram_builder}
\alias{chromatogram_builder}
\docType{data}
\title{
chromatogram builder for m/z = 263.1678 in 003.d from cord blood sample
}
\description{
illusterates a chromatogram and baseline vectors to indicate chromatogram
development.
}
\usage{data("chromatogram_builder")}
\format{
  A data frame with 219 observations on the following 6 variables.
  \describe{
    \item{\code{ScanNumber}}{a numeric vector}
    \item{\code{RetentionTIme}}{a numeric vector}
    \item{\code{SmoothedChromatogram}}{a numeric vector}
    \item{\code{RawChromatogram}}{a numeric vector}
    \item{\samp{12C/13C Isotopologue Pairs}}{a numeric vector}
    \item{\code{Baseline}}{a numeric vector}
  }
}
\examples{
data(chromatogram_builder)
}
\keyword{datasets}
