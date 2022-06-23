\name{recursive_mass_correction}
\alias{recursive_mass_correction}
\title{
recursive mass correction
}
\description{
This function performs recursive mass correction.
}
\usage{
recursive_mass_correction(peaklist, spec_scan, scan_tol,
spectraList, RetentionTime, mass_accuracy_xic, smoothing_window,
peak_resolving_power, min_nIsoPair, min_peak_height, min_ratio_IsoPair,
max_rpw, min_snr_baseline, max_R13C_integrated_peak,
max_percentage_missing_scans, n_spline)
}
\arguments{
  \item{peaklist}{
an IPA peaklist from 'primary_peak_analyzer' function.
}
  \item{spec_scan}{
a matrix consists of 5 columns. The column contents are the m/z of 12C
isotopologues, intensity of 12C isotopologues, scan number (t), m/z of 13C
isotopologues, and intensity of 13C isotopologues.
}
  \item{scan_tol}{
a scan tolerance to extend the chromatogram for better calculations.
}
  \item{spectraList}{
a list of mass spectra in each chromatogram scan.
}
  \item{RetentionTime}{
a vector of retention times for corresponding scan numbers.
}
  \item{mass_accuracy_xic}{
an m/z value to perform chromatography analysis.
}
  \item{smoothing_window}{
a number of scans for peak smoothing.
}
  \item{peak_resolving_power}{
a value to represent peak resolving power.
}
  \item{min_nIsoPair}{
minimum number of nIsoPair for an individual peak.
}
  \item{min_peak_height}{
minimum peak height for an individual peak.
}
  \item{min_ratio_IsoPair}{
minimum ratio of nIsoPair per number of available scans within an individual
peak.
}
  \item{max_rpw}{
maximum allowed value of ratio of peak width at half-height to baseline (RPW)
for an individual peak.
}
  \item{min_snr_baseline}{
minimum S/N baseline for an individual peak.
}
  \item{max_R13C_integrated_peak}{
maximum allowed value of average R13C for an individual peak.
}
  \item{max_percentage_missing_scans}{
maximum allowed value of percentage missing scans on the raw chromatogram for an individual peak.
}
  \item{n_spline}{
number of points for further smoothing using a cubic spline smoothing method to
calculate ancillary chromatographical parameters.
}
}
\value{
a data frame consisting of 24 columns representing chromatography and mass
spectrometry parameters. Each row represents an individual seperated
chromatographic peak.
}
