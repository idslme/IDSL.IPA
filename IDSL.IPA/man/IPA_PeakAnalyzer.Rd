\name{IPA_PeakAnalyzer}
\alias{IPA_PeakAnalyzer}
\title{
IPA Peak Analyzer
}
\description{
This function performs the IPA peak detection module.
}
\usage{
IPA_PeakAnalyzer(PARAM)
}
\arguments{
  \item{PARAM}{
is a data frame from IPA_xlsxAnalyzer function.
}
}
\value{
This function saves individual peaklist files in `.csv` and `.Rdata` formats for HRMS files in the `peaklists` folder.
}
