\name{peak_alignment}
\alias{peak_alignment}
\title{
Peak alignment
}
\description{
This function aligns peaks from multiple peaklists and produce a peak table to
find common peaks among multiple samples.
}
\usage{
peak_alignment(input_path_pl, file_names_pl, RT_pl, mz_error, rt_tol,
n_quantile, number_processing_cores)
}
\arguments{
  \item{input_path_pl}{
path to directory of peaklists.
}
  \item{file_names_pl}{
name of peaklists for peak table production.
}
  \item{RT_pl}{
a list of corrected or uncorrected retention times for each peaklist.
}
  \item{mz_error}{
mass error to detect common peaks.
}
  \item{rt_tol}{
retention time tolerance to detect common peaks.
}
  \item{n_quantile}{
number of total m/z quantiles to split the whole table for faster processing.
}
  \item{number_processing_cores}{
number of processing cores.
}
}
\value{
This function returns an aligned peak table with index numbers from individual peaklists for each peak.
}
