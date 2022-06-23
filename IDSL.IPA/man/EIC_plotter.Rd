\name{EIC_plotter}
\alias{EIC_plotter}
\title{
EIC plotter
}
\description{
This function plots the EIC figure and annex the chromatographic
properties to the EIC figures.
}
\usage{
EIC_plotter(spec_scan_xic, peak_property_xic, smoothing_window,
peak_resolving_power, mass_accuracy_xic, spectraList, RetentionTime,
mz_target, rt_target, file_name, legend_EIC)
}
\arguments{
  \item{spec_scan_xic}{
a matrix consists of 5 columns. The column contents are the m/z
of 12C isotopologues, intensity of 12C isotopologues, scan number (t),
m/z of 13C isotopologues, and intensity of 13C isotopologues.
}
  \item{peak_property_xic}{
a data frame representing chromatographic peak properties.
}
  \item{smoothing_window}{
number of scans for peak smoothing.
}
  \item{peak_resolving_power}{
a value to represent peak resolving power.
}
  \item{mass_accuracy_xic}{
a mass accuracy value to perform chromatography analysis.
}
  \item{spectraList}{
a list of mass spectra in each chromatogram scan.
}
  \item{RetentionTime}{
a vector of retention times vs. corresponding scan numbers.
}
  \item{mz_target}{
an m/z value to perform chromatography analysis.
}
  \item{rt_target}{
the retention time value of the candidate peak.
}
  \item{file_name}{
name of HRMS file used for peak construction.
}
  \item{legend_EIC}{
A file to attach the legends on the EIC figures.
}
}
\value{
A figure to show the EIC and its property table.
}
