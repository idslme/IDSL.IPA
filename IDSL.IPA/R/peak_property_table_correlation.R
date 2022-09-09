peak_property_table_correlation <- function(peakPropertyTable, RTtolerance = 0.05, minFreqDetection = 1, method = "pearson", minThresholdCorrelation = 0, number_processing_threads = 1) {
  ##
  dimPeakPropertyTable <- dim(peakPropertyTable)
  peakPropertyTable <- matrix(peakPropertyTable, ncol = dimPeakPropertyTable[2])
  ##
  ##############################################################################
  ##
  call_correlationList <- function(i) {
    ##
    tPeak <- NULL
    ##
    if (xS[i] >= minFreqDetection) {
      x_t <- which((abs(peakPropertyTable[i, 2] - peakPropertyTable[, 2]) <= RTtolerance))
      x_t <- setdiff(x_t, i)
      ##
      if (length(x_t) > 0) {
        ##
        xiPeak <- xSList[[i]]
        ##
        tPeak <- do.call(rbind, lapply(x_t, function(t) {
          if (xS[t] >= minFreqDetection) {
            ##
            xtPeak <- xSList[[t]]
            ##
            xCommonPeak <- xiPeak[xiPeak %in% xtPeak]
            ##
            if (length(xCommonPeak) >= minFreqDetection) {
              xCommonPeak <- 2 + xCommonPeak
              ##
              corRho <- cor(peakPropertyTable[i, xCommonPeak], peakPropertyTable[t, xCommonPeak], method = method)
              if (!is.na(corRho)) {
                if (corRho >= minThresholdCorrelation) {
                  c(i, xS[i], t, xS[t], corRho)
                }
              }
            }
          }
        }))
        ##
        if (!is.null(tPeak)) {
          colnames(tPeak) <- c("PeakID1", "Freq1", "PeakID2", "Freq2", "Rho")
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
    xSList <- lapply(1:dimPeakPropertyTable[1], function(i) {
      which(peakPropertyTable[i, 3:dimPeakPropertyTable[2]] > 0)
    })
    ##
    xS <- do.call(c, lapply(1:dimPeakPropertyTable[1], function(i) {
      length(xSList[[i]])
    }))
    ##
    progressBARboundaries <- txtProgressBar(min = 0, max = dimPeakPropertyTable[1], initial = 0, style = 3)
    ##
    correlationList <- lapply(1:dimPeakPropertyTable[1], function(i) {
      setTxtProgressBar(progressBARboundaries, i)
      ##
      tryCatch(call_correlationList(i), error = function(e) {NULL}, warning = function(w) {NULL})
    })
    ##
    close(progressBARboundaries)
    ##
  } else {
    ## Processing OS
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Linux") {
      ##
      xSList <- mclapply(1:dimPeakPropertyTable[1], function(i) {
        which(peakPropertyTable[i, 3:dimPeakPropertyTable[2]] > 0)
      }, mc.cores = number_processing_threads)
      ##
      xS <- do.call(c, mclapply(1:dimPeakPropertyTable[1], function(i) {
        length(xSList[[i]])
      }, mc.cores = number_processing_threads))
      ##
      correlationList <- mclapply(1:dimPeakPropertyTable[1], function(i) {
        tryCatch(call_correlationList(i), error = function(e) {NULL}, warning = function(w) {NULL})
      }, mc.cores = number_processing_threads)
      ##
      closeAllConnections()
      ##
    } else if (osType == "Windows") {
      ##
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      xSList <- foreach(i = 1:dimPeakPropertyTable[1], .verbose = FALSE) %dopar% {
        which(peakPropertyTable[i, 3:dimPeakPropertyTable[2]] > 0)
      }
      ##
      xS <- foreach(i = 1:dimPeakPropertyTable[1], .combine = 'c', .verbose = FALSE) %dopar% {
        length(xSList[[i]])
      }
      ##
      correlationList <- foreach(i = 1:dimPeakPropertyTable[1], .verbose = FALSE) %dopar% {
        tryCatch(call_correlationList(i), error = function(e) {NULL}, warning = function(w) {NULL})
      }
      ##
      stopCluster(clust)
      ##
    }
  }
  ##
  return(correlationList)
}
