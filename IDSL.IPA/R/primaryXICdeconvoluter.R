primaryXICdeconvoluter <- function(spectraScan, scanTolerance, indexXIC, aggregatedSpectraList, retentionTime,
                                   massAccuracy, smoothingWindow, peakResolvingPower, minNIonPair, minPeakHeight,
                                   minRatioIonPair, maxRPW, minSNRbaseline, maxR13CcumulatedIntensity,
                                   maxPercentageMissingScans, nSpline, exportEICparameters = NULL, number_processing_threads = 1) {
  ##
  LretentionTime <- length(retentionTime)
  ##
  call_primaryXICdeconvoluter <- function(i) {
    ##
    mzTarget <- spectraScan[i[1], 1]
    ##
    Lx <- length(i)
    ##
    spectraScanPrimary <- spectraScan[i, ]
    if (Lx == 1) {
      spectraScanPrimary <- matrix(spectraScan, nrow = 1)
    } else {
      spectraScanPrimary <- spectraScanPrimary[order(spectraScanPrimary[, 3], decreasing = FALSE), ]
    }
    ##
    scanNumberStart <- spectraScanPrimary[1, 3] - scanTolerance
    if (scanNumberStart < 1) {
      scanNumberStart <- 1
    }
    scanNumberEnd <- spectraScanPrimary[Lx, 3] + scanTolerance
    if (scanNumberEnd > LretentionTime) {
      scanNumberEnd <- LretentionTime
    }
    ##
    chromatographicPeakAnalysis(spectraScanPrimary, aggregatedSpectraList, retentionTime, LretentionTime, massAccuracy, mzTarget,
                                rtTarget = NULL, scanNumberStart, scanNumberEnd, smoothingWindow, peakResolvingPower, minNIonPair,
                                minPeakHeight, minRatioIonPair, maxRPW, minSNRbaseline, maxR13CcumulatedIntensity,
                                maxPercentageMissingScans, nSpline, exportEICparameters)
  }
  ##
  ##############################################################################
  ##############################################################################
  ##
  if (number_processing_threads == 1) {
    ##
    peaklist <- do.call(rbind, lapply(indexXIC, function(i) {
      call_primaryXICdeconvoluter(i)
    }))
    ##
  } else {
    ## Processing OS
    osType <- Sys.info()[['sysname']]
    ##
    ############################################################################
    ##
    if (osType == "Windows") {
      ##
      clust <- makeCluster(number_processing_threads)
      clusterExport(clust, setdiff(ls(), c("clust", "indexXIC")), envir = environment())
      ##
      peaklist <- do.call(rbind, parLapply(clust, indexXIC, function(i) {
        call_primaryXICdeconvoluter(i)
      }))
      ##
      stopCluster(clust)
      ##
      ##########################################################################
      ##
    } else {
      ##
      peaklist <- do.call(rbind, mclapply(indexXIC, function(i) {
        call_primaryXICdeconvoluter(i)
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
      ##
    }
  }
  ##
  ##############################################################################
  ##############################################################################
  ##
  return(peaklist)
}
