IPA_MSdeconvoluter <- function(HRMS_path, MSfile, MS_level = 1) {
  p2l <- IDSL.MXP::peak2list(HRMS_path, MSfile)
  scanTable <- p2l[["scanTable"]]
  spectraList <- p2l[["spectraList"]]
  x_MS <- which(scanTable$peaksCount > 0 & scanTable$msLevel == MS_level) ## some files may not have data in calibration scan acquisitions.
  spectraList <- spectraList[x_MS]
  scanTable <- scanTable[x_MS, ]
  RetentionTime <- as.numeric(scanTable$retentionTime)  # Retention times in minute
  MS_polarity <- ifelse(scanTable$polarity[x_MS[1]] == 1, "+", "-")
  outputer <- list(spectraList, RetentionTime, MS_polarity)
  return(outputer)
}
