\name{peak_alignment}
\alias{peak_alignment}
\title{
Peak alignment
}
\description{
This function aligns peaks from multiple peaklists and produce a peak table to find common peaks among multiple samples.
}
\usage{
peak_alignment(peaklistInputFolderPath, peaklistFileNames, listCorrectedRTpeaklists,
massAccuracy, RTtolerance, noQuantile, number_processing_threads = 1)
}
\arguments{
  \item{peaklistInputFolderPath}{
path to directory of peaklists.
}
  \item{peaklistFileNames}{
name of peaklists for peak table production.
}
  \item{listCorrectedRTpeaklists}{
a list of corrected or uncorrected retention times for each peaklist.
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
This function returns an aligned peak table with index numbers from individual peaklists for each peak.
}
