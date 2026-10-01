\name{gapFillingCore}
\alias{gapFillingCore}
\title{
Gap-Filling Core Function
}
\description{
Gap-Filling Core Function
}
\usage{
gapFillingCore(input_path_hrms, peakXcol, massAccuracy, smoothingWindow,
peakResolvingPower, RTtolerance, scanTolerance, nSpline,
retentionTimeCorrectionCheck = FALSE, listCorrectedRTpeaklists = NULL,
inputPathPeaklist = NULL, ionMassDifference = 1.003354835336,
output_path = NA, number_processing_threads = 1)
}
\arguments{
  \item{input_path_hrms}{
  input_path_hrms
}
  \item{peakXcol}{
peakXcol
}
  \item{massAccuracy}{
massAccuracy
}
  \item{smoothingWindow}{
number of scans for peak smoothing
}
  \item{peakResolvingPower}{
a value to represent peak resolving power
}
  \item{RTtolerance}{
RTtolerance
}
  \item{scanTolerance}{
a scan tolerance to extend the chromatogram for better calculations.
}
  \item{nSpline}{
number of points for further smoothing using a cubic spline smoothing method to calculate ancillary chromatographic parameters
}
  \item{retentionTimeCorrectionCheck}{
retentionTimeCorrectionCheck
}
  \item{listCorrectedRTpeaklists}{
listCorrectedRTpeaklists
}
  \item{inputPathPeaklist}{
inputPathPeaklist
}
  \item{ionMassDifference}{
ionMassDifference
}
  \item{output_path}{
output_path
}
  \item{number_processing_threads}{
number of processing threads
}
}
\value{
A list of gap-filled data
}
