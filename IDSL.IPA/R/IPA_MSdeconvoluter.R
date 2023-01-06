IPA_MSdeconvoluter <- function(inputHRMSfolderPath, MSfileName, MSlevel = 1) {
  p2l <- IDSL.MXP::peak2list(inputHRMSfolderPath, MSfileName)
  scanTable <- p2l[["scanTable"]]
  spectraList <- p2l[["spectraList"]]
  x_MS <- which(scanTable$peaksCount > 0 & scanTable$msLevel == MSlevel) ## some files may not have data in calibration scan acquisitions.
  spectraList <- spectraList[x_MS]
  scanTable <- scanTable[x_MS, ]
  retentionTime <- as.numeric(scanTable$retentionTime)  # Retention times in minute
  msPolarity <- if (scanTable$polarity[1] == 1) {"+"} else {"-"}
  outputer <- list(spectraList, retentionTime, msPolarity)
  names(outputer)  <- c("spectraList", "retentionTime", "msPolarity")
  return(outputer)
}
