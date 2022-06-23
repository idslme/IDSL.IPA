\name{reference_peaks_detector}
\alias{reference_peaks_detector}
\title{
Reference peaks detector
}
\description{
This function detects reocurring reference peaks (m/z-RT) for retention time
correction.
}
\usage{
reference_peaks_detector(input_path_pl, file_names_ref, min_frequency_ref_peaks,
RT_pl_ref, mz_error, rt_tol, n_quantile, number_processing_cores)
}
\arguments{
  \item{input_path_pl}{
path to directory of peaklists.
}
  \item{file_names_ref}{
name of peaklists files to detect reocurring reference peaks (m/z-RT).
}
  \item{min_frequency_ref_peaks}{
minimum frequency of the reocurring reference peaks (m/z-RT) in the reference
files.
}
  \item{RT_pl_ref}{
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
a matrix of two columns of m/z and RT of common peaks in the reference samples.
}
