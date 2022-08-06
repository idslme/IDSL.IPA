\name{peak_Xcol2}
\alias{peak_Xcol2}
\title{
Peak table producer
}
\description{
This function fills the peak table from individual peaklists.
}
\usage{
peak_Xcol2(input_path_peaklist, file_names_peaklist, peak_Xcol)
}
\arguments{
  \item{input_path_peaklist}{
address of the peaklists.
}
  \item{file_names_peaklist}{
a vector of the peaklist file names.
}
  \item{peak_Xcol}{
a matrix of index numbers in individual peaklists for each peak (m/z-RT).
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
