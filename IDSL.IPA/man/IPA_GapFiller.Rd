\name{IPA_GapFiller}
\alias{IPA_GapFiller}
\title{
IPA GapFiller
}
\description{
This function fills the gaps on the peak table.
}
\usage{
IPA_GapFiller(PARAM)
}
\arguments{
  \item{PARAM}{
a data frame from IPA_xlsxAnalyzer function containing the IPA parameters.
}
}
\value{
This function saves individual .CSV and .Rdata files for the gap-filled peak tables for peak height, area, and R13C properties in the "peak_alignment" folder.
}
