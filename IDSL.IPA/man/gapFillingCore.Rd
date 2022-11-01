\name{gapFillingCore}
\alias{gapFillingCore}
\title{
Gap-Filling Core Function
}
\description{
Gap-Filling Core Function
}
\usage{
gapFillingCore(input_path_hrms, peakXcol, massAccuracy, RTtolerance, scanTolerance,
retentionTimeCorrectionCheck = FALSE, listCorrectedRTpeaklists = NULL,
inputPathPeaklist = NULL, ionMassDifference = 1.003354835336,
number_processing_threads = 1)
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
  \item{RTtolerance}{
RTtolerance
}
  \item{scanTolerance}{
scanTolerance
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
  \item{number_processing_threads}{
number of processing threads
}
}
\value{
A list of gap-filled data
}
