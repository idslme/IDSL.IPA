alignedPeakPropertyTableCorrelationListCalculator <- function(peakPropertyTable, RTtolerance = 0.05, minFreqDetection = 3, minRatioDetection = 0.01,
                                                              method = "pearson", minThresholdCorrelation = 0, number_processing_threads = 1) {
  ##
  dimPeakPropertyTable <- dim(peakPropertyTable)
  nPeaks <- dimPeakPropertyTable[1]
  LpeakPropertyTable <- dimPeakPropertyTable[2]
  nSamples <- LpeakPropertyTable - 5
  peakPropertyTable <- matrix(peakPropertyTable, ncol = LpeakPropertyTable)
  ##
  ##############################################################################
  ##
  call_alignedPeakPropertyTableCorrelationListCalculator <- function(i) {
    ##
    tPeak <- NULL
    ##
    if (xS[i] >= minFreqDetection) {
      if ((xS[i]/nSamples) >= minRatioDetection) {
        x_t <- which(abs(peakPropertyTable[i, 2] - peakPropertyTable[, 2]) <= RTtolerance)
        x_t <- setdiff(x_t, i)
        ##
        if (length(x_t) > 0) {
          ##
          xiPeak <- xSList[[i]]
          ##
          tPeak <- do.call(c, lapply(x_t, function(t) {
            if (xS[t] >= minFreqDetection) {
              if ((xS[t]/nSamples) >= minRatioDetection) {
                ##
                xtPeak <- xSList[[t]]
                xCommonPeak <- xiPeak[(xiPeak %in% xtPeak)]
                ##
                LxCommonPeak <- length(xCommonPeak)
                if (LxCommonPeak >= minFreqDetection) {
                  if ((LxCommonPeak/nSamples) >= minRatioDetection) {
                    xCommonPeak <- 5 + xCommonPeak
                    ##
                    corRho <- suppressWarnings(cor(peakPropertyTable[i, xCommonPeak], peakPropertyTable[t, xCommonPeak], method = method))
                    if (!is.na(corRho)) {
                      if (corRho >= minThresholdCorrelation) {
                        t
                      }
                    }
                  }
                }
              }
            }
          }))
        }
      }
    }
    return(tPeak)
  }
  ##
  ##############################################################################
  ##
  if (number_processing_threads == 1) {
    ##
    xSList <- lapply(1:nPeaks, function(i) {
      which(peakPropertyTable[i, 6:LpeakPropertyTable] > 0)
    })
    ##
    xS <- do.call(c, lapply(1:nPeaks, function(i) {
      length(xSList[[i]])
    }))
    ##
    progressBARboundaries <- txtProgressBar(min = 0, max = nPeaks, initial = 0, style = 3)
    ##
    correlationList <- lapply(1:nPeaks, function(i) {
      setTxtProgressBar(progressBARboundaries, i)
      ##
      tryCatch(call_alignedPeakPropertyTableCorrelationListCalculator(i), error = function(e) {NULL}, warning = function(w) {NULL})
    })
    ##
    close(progressBARboundaries)
    ##
  } else {
    ## Processing OS
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Windows") {
      ##
      clust <- makeCluster(number_processing_threads)
      clusterExport(clust, c("peakPropertyTable", "LpeakPropertyTable"), envir = environment())
      xSList <- parLapply(clust, 1:nPeaks, function(i) {
        which(peakPropertyTable[i, 6:LpeakPropertyTable] > 0)
      })
      stopCluster(clust)
      ##
      clust <- makeCluster(number_processing_threads)
      clusterExport(clust, "xSList", envir = environment())
      xS <- do.call(c, parLapply(clust, 1:nPeaks, function(i) {
        length(xSList[[i]])
      }))
      stopCluster(clust)
      ##
      clust <- makeCluster(number_processing_threads)
      clusterExport(clust, setdiff(ls(), c("clust", "nPeaks")), envir = environment())
      correlationList <- parLapply(clust, 1:nPeaks, function(i) {
        tryCatch(call_alignedPeakPropertyTableCorrelationListCalculator(i), error = function(e) {NULL}, warning = function(w) {NULL})
      })
      stopCluster(clust)
      ##
      ##########################################################################
      ##
    } else {
      ##
      xSList <- mclapply(1:nPeaks, function(i) {
        which(peakPropertyTable[i, 6:LpeakPropertyTable] > 0)
      }, mc.cores = number_processing_threads)
      ##
      xS <- do.call(c, mclapply(1:nPeaks, function(i) {
        length(xSList[[i]])
      }, mc.cores = number_processing_threads))
      ##
      correlationList <- mclapply(1:nPeaks, function(i) {
        tryCatch(call_alignedPeakPropertyTableCorrelationListCalculator(i), error = function(e) {NULL}, warning = function(w) {NULL})
      }, mc.cores = number_processing_threads)
      ##
      closeAllConnections()
      ##
    }
  }
  ##
  return(correlationList)
}
