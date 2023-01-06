peakPropertyTableMedianCalculator <- function(peakPropertyTable, falggingVector = NULL, number_processing_threads = 1, allowedVerbose = TRUE) {
  ##
  dimPeakPropertyTable <- dim(peakPropertyTable)
  nPeaks <- dimPeakPropertyTable[1]
  nSamples3 <- dimPeakPropertyTable[2]
  ##
  ##############################################################################
  ##############################################################################
  ##
  call_peakPropertyTableMedianCalculator <- function(i) {
    x_h <- which(peakPropertyTable[i, 4:nSamples3] != 0)
    if (length(x_h) > 0) {
      median(peakPropertyTable[i, (x_h + 3)])
    } else {
      0
    }
  }
  ##
  ##############################################################################
  ##############################################################################
  ##
  if (number_processing_threads == 1) {
    ##
    if (allowedVerbose) {progressBARboundaries <- txtProgressBar(min = 0, max = nPeaks, initial = 0, style = 3)}
    ##
    medianPeakProperty <- do.call(c, lapply(1:nPeaks, function(i) {
      if (allowedVerbose) {setTxtProgressBar(progressBARboundaries, i)}
      ##
      call_peakPropertyTableMedianCalculator(i)
    }))
    ##
    if (allowedVerbose) {close(progressBARboundaries)}
    ##
  } else {
    ## Processing OS
    osType <- Sys.info()[['sysname']]
    ##
    ############################################################################
    ##
    if (osType == "Linux") {
      ##
      medianPeakProperty <- do.call(c, mclapply(1:nPeaks, function(i) {
        call_peakPropertyTableMedianCalculator(i)
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
      ##
      ##########################################################################
      ##
    } else if (osType == "Windows") {
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      medianPeakProperty <- foreach(i = 1:nPeaks, .combine = 'c', .verbose = FALSE) %dopar% {
        call_peakPropertyTableMedianCalculator(i)
      }
      ##
      stopCluster(clust)
      ##
    }
  }
  ##
  ##############################################################################
  ##############################################################################
  ##
  names(medianPeakProperty) <- "medianNon0"
  ##
  peakPropertyTable <- cbind(peakPropertyTable[, c(1, 2, 3)], medianPeakProperty, falggingVector, peakPropertyTable[, 4:nSamples3])
  rownames(peakPropertyTable) <- NULL
  ##
  return(peakPropertyTable)
}
