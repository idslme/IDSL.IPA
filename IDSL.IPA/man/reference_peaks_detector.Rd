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
reference_peaks_detector(input_path_peaklist, file_names_peaklist_ref,
min_frequency_ref_peaks, mz_error, rt_tol, n_quantile, number_processing_threads = 1)
}
\arguments{
  \item{input_path_peaklist}{
path to directory of peaklists.
}
  \item{file_names_peaklist_ref}{
name of peaklists files to detect reocurring reference peaks (m/z-RT).
}
  \item{min_frequency_ref_peaks}{
minimum frequency of the reocurring reference peaks (m/z-RT) in the reference
files.
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
  \item{number_processing_threads}{
number of processing threads
}
}
\value{
  \item{reference_mz_rt_peaks}{
a matrix of two columns of m/z and RT of common peaks in the reference samples.
}
  \item{listRefRT}{
a list of corrected or uncorrected retention times for each peaklist.
}
}
