peakPropertyTableFreqCalculator <- function(peakPropertyTable, startColumnIndex = 3, number_processing_threads = 1, allowedVerbose = TRUE) {
  ##
  dimPeakPropertyTable <- dim(peakPropertyTable)
  nPeaks <- dimPeakPropertyTable[1]
  endColumnIndex <- dimPeakPropertyTable[2]
  ##
  ##############################################################################
  ##############################################################################
  ##
  if (number_processing_threads == 1) {
    ##
    if (allowedVerbose) {progressBARboundaries <- txtProgressBar(min = 0, max = nPeaks, initial = 0, style = 3)}
    ##
    freqPeakProperty <- do.call(c, lapply(1:nPeaks, function(i) {
      if (allowedVerbose) {setTxtProgressBar(progressBARboundaries, i)}
      ##
      length(which(peakPropertyTable[i, startColumnIndex:endColumnIndex] != 0))
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
      freqPeakProperty <- do.call(c, mclapply(1:nPeaks, function(i) {
        length(which(peakPropertyTable[i, startColumnIndex:endColumnIndex] != 0))
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
      freqPeakProperty <- foreach(i = 1:nPeaks, .combine = 'c', .verbose = FALSE) %dopar% {
        length(which(peakPropertyTable[i, startColumnIndex:endColumnIndex] != 0))
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
  return(freqPeakProperty)
}
