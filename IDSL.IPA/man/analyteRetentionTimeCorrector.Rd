\name{analyteRetentionTimeCorrector}
\alias{analyteRetentionTimeCorrector}
\title{
analyte retention time corrector
}
\description{
This function calculates corrected retention times for the peaklists.
}
\usage{
analyteRetentionTimeCorrector(referenceMZRTpeaks, inputPathMZML,
inputPathPeaklist, peaklistFileName, massAccuracy, RTtolerance,
RTcorrectionIonSource, refPeakTolerance = 1, degreePolynomial = 3)
}
\arguments{
  \item{referenceMZRTpeaks}{
a matrix of reference peaks for retention time correction.
}
  \item{inputPathMZML}{
input path to the raw mass spectrometry data files (e.g., mzML, mzXML, or CDF).
}
  \item{inputPathPeaklist}{
input path to peaklist.
}
  \item{peaklistFileName}{
file name of the peaklist.
}
  \item{massAccuracy}{
mass error to detect common reference peaks.
}
  \item{RTtolerance}{
retention time tolerance for matching reference peaks.
}
  \item{RTcorrectionIonSource}{
method for retention time correction, e.g., c('RetentionIndex','Polynomial').
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
