\name{IPA_PeaklistAnnotation}
\alias{IPA_PeaklistAnnotation}
\title{
IPA Peaklist Annotation
}
\description{
This function performs sample-centric peak annotation.
}
\usage{
IPA_PeaklistAnnotation(PARAM)
}
\arguments{
  \item{PARAM}{
a data frame from IPA_xlsxAnalyzer function.
}
}
\value{
This function saves individual .CSV files for peak height, area, and R13C properties in the "sample_centeric_annotation" folder.
}
