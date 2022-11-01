\name{peakXcolFiller}
\alias{peakXcolFiller}
\title{
Peak table producer
}
\description{
This function fills the peak table from individual peaklists.
}
\usage{
peakXcolFiller(peakXcol, inputPathPeaklist)
}
\arguments{
  \item{peakXcol}{
a matrix of index numbers in individual peaklists for each peak (m/z-RT).
}
  \item{inputPathPeaklist}{
address of the peaklists.
}
}
\value{
  \item{peak_height}{
peak table for height values
}
  \item{peak_area}{
peak table for area values
}
  \item{peak_R13C}{
peak table for R13C values
}
}
