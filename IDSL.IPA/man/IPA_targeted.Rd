\name{IPA_targeted}
\alias{IPA_targeted}
\title{
IPA Targeted Analysis
}
\description{
This function plots extracted ion chromatogram (EIC) figures in the targeted mode.
}
\usage{
IPA_targeted(PARAM_targeted, allowedVerbose = TRUE)
}
\arguments{
  \item{PARAM_targeted}{
IPA parameters to feed the `IPA_targeted` module. This variable can be produced using the `IPA_targeted_xlsxAnalyzer` module.
}
  \item{allowedVerbose}{
c(TRUE, FALSE). A `TRUE` allowedVerbose provides messages about the flow of the function.
}
}
\value{
This module saves extracted ion chromatograms (EICs) in .png format in the "Targeted_EICs" folder and saves a table of peak properties.
}
