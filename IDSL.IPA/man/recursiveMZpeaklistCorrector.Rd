\name{recursiveMZpeaklistCorrector}
\alias{recursiveMZpeaklistCorrector}
\title{
recursive mass correction
}
\description{
This function performs recursive mass correction.
}
\usage{
recursiveMZpeaklistCorrector(peaklist, spectraScan, scanTolerance,
aggregatedSpectraList, retentionTime, massAccuracy, smoothingWindow,
peakResolvingPower, minNIonPair, minPeakHeight, minRatioIonPair, maxRPW,
minSNRbaseline, maxR13CcumulatedIntensity, maxPercentageMissingScans, nSpline,
exportEICparameters = NULL, number_processing_threads = 1)
}
\arguments{
  \item{peaklist}{
an IPA peaklist from 'primaryXICdeconvoluter' function.
}
  \item{spectraScan}{
a matrix consists of 5 columns. The column contents are the m/z of 12C isotopologues, intensity of 12C isotopologues, scan number (t), m/z of 13C isotopologues, and intensity of 13C isotopologues.
}
  \item{scanTolerance}{
a scan tolerance to extend the chromatogram for better calculations.
}
  \item{aggregatedSpectraList}{
aggregated spectraList and spectra matrix from the `IPA_spectraListAggregator` module
}
  \item{retentionTime}{
a vector of retention times for corresponding scan numbers.
}
  \item{massAccuracy}{
an m/z value to perform chromatography analysis.
}
  \item{smoothingWindow}{
a number of scans for peak smoothing.
}
  \item{peakResolvingPower}{
a value to represent peak resolving power.
}
  \item{minNIonPair}{
minimum number of nIsoPair for an individual peak.
}
  \item{minPeakHeight}{
minimum peak height for an individual peak.
}
  \item{minRatioIonPair}{
minimum ratio of nIsoPair per number of available scans within an individual peak.
}
  \item{maxRPW}{
maximum allowed value of ratio of peak width at half-height to baseline (RPW) for an individual peak.
}
  \item{minSNRbaseline}{
minimum S/N baseline for an individual peak.
}
  \item{maxR13CcumulatedIntensity}{
maximum allowed value of average R13C for an individual peak.
}
  \item{maxPercentageMissingScans}{
maximum allowed value of percentage missing scans on the raw chromatogram for an individual peak.
}
  \item{nSpline}{
number of points for further smoothing using a cubic spline smoothing method to calculate ancillary chromatographic parameters.
}
  \item{exportEICparameters}{
When `NULL`, EICs are not plotted. `exportEICparameters` should contain three variables of 1) an address to save IPA EICs figures, 2) 'HRMS' file name, and 3) a valid string of characters.
}
  \item{number_processing_threads}{
number of processing threads
}
}
\value{
a dataframe consisting of 24 columns representing chromatography and mass spectrometry parameters. Each row represents an individual separated chromatographic peak.
}
