\name{chromatographicPeakAnalysis}
\alias{chromatographicPeakAnalysis}
\title{
Chromatography analysis
}
\description{
This function detects individual chromatographic peaks and measures their peak qualification metrics.
}
\usage{
chromatographicPeakAnalysis(spectraScanXIC, aggregatedSpectraList, retentionTime,
LretentionTime, massAccuracy, mzTarget, rtTarget = NULL, scanNumberStart,
scanNumberEnd, smoothingWindow, peakResolvingPower, minNIonPair, minPeakHeight,
minRatioIonPair, maxRPW, minSNRbaseline, maxR13CcumulatedIntensity,
maxPercentageMissingScans, nSpline, exportEICparameters  = NULL)
}
\arguments{
  \item{spectraScanXIC}{
a matrix consists of 5 columns. The column contents are the m/z of 12C isotopologues, intensity of 12C isotopologues, scan number (t), m/z of 13C isotopologues, and intensity of 13C isotopologues, respectively. Redundant scan numbers are not allowed for this module.
}
  \item{aggregatedSpectraList}{
aggregated spectraList and spectra matrix from the `IPA_spectraListAggregator` module
}
  \item{retentionTime}{
a vector of retention times vs. corresponding scan numbers
}
  \item{LretentionTime}{
length of the retention time vector
}
  \item{massAccuracy}{
mass error to perform chromatography analysis
}
  \item{mzTarget}{
m/z value to perform chromatography analysis
}
  \item{rtTarget}{
retention time value for a targeted peak to calculate the ancillary chromatography parameters. When this parameter set at 0, the ancillary chromatography parameters are calculated for the entire detected peaks.
}
  \item{scanNumberStart}{
the first scan number.
}
  \item{scanNumberEnd}{
the last scan number.
}
  \item{smoothingWindow}{
number of scans for peak smoothing
}
  \item{peakResolvingPower}{
a value to represent peak resolving power
}
  \item{minNIonPair}{
minimum number of nIsoPair for an individual peak
}
  \item{minPeakHeight}{
minimum peak height for an individual peak
}
  \item{minRatioIonPair}{
minimum ratio of nIsoPair per number of available scans within an individual peak
}
  \item{maxRPW}{
maximum allowed value of ratio of peak width at half-height to baseline (RPW) for an individual peak
}
  \item{minSNRbaseline}{
minimum S/N baseline for an individual peak
}
  \item{maxR13CcumulatedIntensity}{
maximum allowed value of average R13C for an individual peak
}
  \item{maxPercentageMissingScans}{
maximum allowed value of percentage missing scans on the raw chromatogram for an individual peak.
}
  \item{nSpline}{
number of points for further smoothing using a cubic spline smoothing method to calculate ancillary chromatographic parameters
}
  \item{exportEICparameters}{
When `NULL`, EICs are not plotted. `exportEICparameters` should contain three variables of 1) an address to save IPA EICs figures, 2) 'HRMS' file name, and 3) a valid string of characters.
}
}
\value{
a data frame consisting of 24 columns representing chromatography and mass spectrometry parameters. Each row represents an individual separated chromatographic peak.
}
