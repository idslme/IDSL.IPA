\name{fronting_tailing_resolver}
\alias{fronting_tailing_resolver}
\title{
Fronting and tailing peaks resolver
}
\description{
This function attempts to resolve peak tailings or frontings into the main peak in case they were detected as seperate peaks.
}
\usage{
fronting_tailing_resolver(segment, int, max_space, peak_resolving_power)
}
\arguments{
  \item{segment}{
a matrix or a vector of peak boundaries.
}
  \item{int}{
a vector of intensities of the entire chromatogram.
}
  \item{max_space}{
maximum scan number difference between peak tailing or fronting and the main peak.
}
  \item{peak_resolving_power}{
power of peak resolving tool.
}
}
\value{
A matrix of 2 columns. Each row indicates peak boundary indices on the 'int' vector after resolving fronting and tailing peaks.
}
\examples{
data(segment)
data(chromatogram_builder)
int <- chromatogram_builder$SmoothedChromatogram
max_space <- 7
peak_resolving_power <- 0.2
fronting_tailing_resolver(segment, int, max_space, peak_resolving_power)
}
