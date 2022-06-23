\name{IPA_PeakAlignment}
\alias{IPA_PeakAlignment}
\title{
IPA peak alignment
}
\description{
This function produce an aligned peak table from individual peaklists.
}
\usage{
IPA_PeakAlignment(PARAM)
}
\arguments{
  \item{PARAM}{
is a data frame from IPA_xlsxAnalyzer function.
}
}
\value{
This function saves individual .CSV and .Rdata files for the aligned peak tables for peak height, area, and R13C properties in the "peak_alignment" folder.
}
