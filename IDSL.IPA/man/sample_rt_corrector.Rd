\name{sample_rt_corrector}
\alias{sample_rt_corrector}
\title{
sample retention time corrector
}
\description{
This function calculates corrected retention times for the peaklists.
}
\usage{
sample_rt_corrector(reference_mz_rt_peaks, peaklist, mz_error,
rt_correction_method, reference_peak_tol = 1, polynomial_degree = 3)
}
\arguments{
  \item{reference_mz_rt_peaks}{
a matrix of reference peaks for retention time correction.
}
  \item{peaklist}{
an IPA peaklist.
}
  \item{mz_error}{
mass error to detect common reference peaks.
}
  \item{rt_correction_method}{
This parameter can be either 'RetentionIndex' or 'Polynomial'.
}
  \item{reference_peak_tol}{
number of reference peaks for retention time correction using 'RetentionIndex'
method.
}
  \item{polynomial_degree}{
polynomial degree for retention time correction using 'Polynomial' method.
}
}
\value{
a list of corrected retention times for each peaklist.
}
