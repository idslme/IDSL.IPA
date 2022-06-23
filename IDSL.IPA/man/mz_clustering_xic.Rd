\name{mz_clustering_xic}
\alias{mz_clustering_xic}
\title{
mz clustering XIC
}
\description{
This function clusters related 12C m/z values.
}
\usage{
mz_clustering_xic(spec_scan, mass_accuracy_xic, min_peak_height, min_nIsoPair)
}
\arguments{
  \item{spec_scan}{
a matrix consists of 3 columns. The column contents are the m/z of 12C
isotopologues, intensity of 12C isotopologues, and scan number (t).
}
  \item{mass_accuracy_xic}{
mass accuracy to detect related 12C m/z values.
}
  \item{min_nIsoPair}{
minimum number of nIsoPair for an individual peak.
}
  \item{min_peak_height}{
minimum peak height for an individual peak.
}
}
\value{
This function returns an list on index numbers of EICs for the "spec_scan" variable.
}
