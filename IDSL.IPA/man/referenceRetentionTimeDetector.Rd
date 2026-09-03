\name{referenceRetentionTimeDetector}
\alias{referenceRetentionTimeDetector}
\title{
Reference retention time detector
}
\description{
This module detects recurring reference peaks (m/z-RT) for retention time correction.
}
\usage{
referenceRetentionTimeDetector(inputPathPeaklist, inputPathMZML,
refPeaklistFileNames, minFrequencyRefPeaks, massAccuracy, RTtolerance,
RTcorrectionIonSource, number_processing_threads = 1)
}
\arguments{
  \item{inputPathPeaklist}{
path to directory of peaklists.
}
  \item{inputPathMZML}{
path to the directory containing raw mass spectrometry data files (e.g., mzML, mzXML, or CDF).
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
  \item{RTcorrectionIonSource}{
method or ion source mode used for retention time correction.
}
  \item{number_processing_threads}{
number of processing threads.
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
