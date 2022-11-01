\name{reference_peaks_detector}
\alias{reference_peaks_detector}
\title{
Reference peaks detector
}
\description{
This function detects recurring reference peaks (m/z-RT) for retention time
correction.
}
\usage{
reference_peaks_detector(inputPathPeaklist, refPeaklistFileNames,
minFrequencyRefPeaks, massAccuracy, RTtolerance, noQuantile,
number_processing_threads = 1)
}
\arguments{
  \item{inputPathPeaklist}{
path to directory of peaklists.
}
  \item{refPeaklistFileNames}{
name of peaklists files to detect recurring reference peaks (m/z-RT).
}
  \item{minFrequencyRefPeaks}{
minimum frequency of the recurring reference peaks (m/z-RT) in the reference
files.
}
  \item{massAccuracy}{
mass error to detect common peaks.
}
  \item{RTtolerance}{
retention time tolerance to detect common peaks.
}
  \item{noQuantile}{
number of total m/z quantiles to split the whole table for faster processing.
}
  \item{number_processing_threads}{
number of processing threads
}
}
\value{
  \item{referenceMZRTpeaks}{
a matrix of two columns of m/z and RT of common peaks in the reference samples.
}
  \item{listRefRT}{
a list of corrected or uncorrected retention times for each peaklist.
}
}
