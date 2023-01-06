\name{alignedPeakPropertyTableCorrelationListCalculator}
\alias{alignedPeakPropertyTableCorrelationListCalculator}
\title{
Aligned Peak Property Table Correlation List Calculator
}
\description{
Aligned Peak Property Table Correlation List Calculator
}
\usage{
alignedPeakPropertyTableCorrelationListCalculator(peakPropertyTable,
RTtolerance = 0.05, minFreqDetection = 3, minRatioDetection = 0.01,
method = "pearson", minThresholdCorrelation = 0, number_processing_threads = 1)
}
\arguments{
  \item{peakPropertyTable}{
peak property table such as `peak_height`, `peak_area' and  `peak_R13C'
}
  \item{RTtolerance}{
retention time tolerance (min)
}
  \item{minFreqDetection}{
minimum frequency of detection for a (m/z-RT) peak across the peak property table
}
  \item{minRatioDetection}{
minimum ratio of detection for a (m/z-RT) peak across the peak property table. This value should be between (0 - 1).
}
  \item{method}{
a character string indicating which correlation coefficient (or covariance) is to be computed. One of "pearson" (default), or "spearman": can be abbreviated. (from `cor` function of the `stats` package)
}
  \item{minThresholdCorrelation}{
minimum threshold for the correlation method
}
  \item{number_processing_threads}{
number of processing threads
}
}
\value{
A list of related peak IDs for each individual (m/z-RT) pair on the peak property table
}
