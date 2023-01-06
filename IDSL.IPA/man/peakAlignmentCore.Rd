\name{peakAlignmentCore}
\alias{peakAlignmentCore}
\title{
Peak Alignment Core
}
\description{
This function aligns peaks from multiple peaklists and produces an aligned table of common peaks among multiple samples.
}
\usage{
peakAlignmentCore(peaklistInputFolderPath, peaklistFileNames, listCorrectedRTpeaklists,
massAccuracy, RTtolerance, number_processing_threads = 1)
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
  \item{number_processing_threads}{
number of processing threads
}
}
\value{
This function returns an aligned peak table with index numbers from individual peaklists for each peak.
}
