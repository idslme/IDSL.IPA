\name{analyteRetentionTimeCorrector}
\alias{analyteRetentionTimeCorrector}
\title{
analyte retention time corrector
}
\description{
This function calculates corrected retention times for the peaklists.
}
\usage{
analyteRetentionTimeCorrector(referenceMZRTpeaks, inputPathPeaklist, peaklistFileName,
massAccuracy, RTcorrectionMethod, refPeakTolerance = 1, degreePolynomial = 3)
}
\arguments{
  \item{referenceMZRTpeaks}{
a matrix of reference peaks for retention time correction.
}
  \item{inputPathPeaklist}{
input path to peaklist
}
  \item{peaklistFileName}{
file name peaklist
}
  \item{massAccuracy}{
mass error to detect common reference peaks.
}
  \item{RTcorrectionMethod}{
c('RetentionIndex','Polynomial')
}
  \item{refPeakTolerance}{
number of reference peaks for retention time correction using the 'RetentionIndex' method.
}
  \item{degreePolynomial}{
polynomial degree for retention time correction using the 'Polynomial' method.
}
}
\value{
a list of corrected retention times for each peaklist.
}
