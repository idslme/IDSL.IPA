\name{sample_rt_corrector}
\alias{sample_rt_corrector}
\title{
sample retention time corrector
}
\description{
This function calculates corrected retention times for the peaklists.
}
\usage{
sample_rt_corrector(reference_mz_rt_peaks, input_path_peaklist, file_name_peaklist,
mz_error, rt_correction_method, reference_peak_tol = 1, polynomial_degree = 3)
}
\arguments{
  \item{reference_mz_rt_peaks}{
a matrix of reference peaks for retention time correction.
}
  \item{input_path_peaklist}{
input path to peaklist
}
  \item{file_name_peaklist}{
file name peaklist
}
  \item{mz_error}{
mass error to detect common reference peaks.
}
  \item{rt_correction_method}{
c('RetentionIndex','Polynomial')
}
  \item{reference_peak_tol}{
number of reference peaks for retention time correction using the 'RetentionIndex' method.
}
  \item{polynomial_degree}{
polynomial degree for retention time correction using the 'Polynomial' method.
}
}
\value{
a list of corrected retention times for each peaklist.
}
