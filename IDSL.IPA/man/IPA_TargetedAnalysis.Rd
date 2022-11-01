\name{IPA_TargetedAnalysis}
\alias{IPA_TargetedAnalysis}
\title{
IPA Targeted Analysis
}
\description{
This function plots extracted ion chromatogram (EIC) figures in the targeted mode.
}
\usage{
IPA_TargetedAnalysis(spreadsheet, mzCandidate, rtCandidate, exportEIC = TRUE,
exportTable = FALSE)
}
\arguments{
  \item{spreadsheet}{
a spreadsheet containing the parameters.
}
  \item{mzCandidate}{
a vector of candidate m/z values.
}
  \item{rtCandidate}{
a vector of candidate RT values.
}
  \item{exportEIC}{
TRUE by default. To plot and save EICs.
}
  \item{exportTable}{
FALSE by default. select TRUE to return the whole peaklists for the m/z and RT candidate vectors.
}
}
\value{
This function saves extracted ion chromatograms in .png format in the "EICs" folder when "exportEIC = TRUE", and saves a table of peak properties when "exportTable = TRUE".
}
