recursiveMZpeaklistCorrector <- function(peaklist, spectraScan, scanTolerance, aggregatedSpectraList, retentionTime, massAccuracy,
                                         smoothingWindow, peakResolvingPower, minNIonPair, minPeakHeight, minRatioIonPair, maxRPW,
                                         minSNRbaseline, maxR13CcumulatedIntensity, maxPercentageMissingScans, nSpline, exportEICparameters = NULL,
                                         number_processing_threads = 1) {
  ##
  LretentionTime <- length(retentionTime)
  ##
  call_recursiveMZpeaklistCorrector <- function(i) {
    rtTarget <- peaklist[i, 3]
    mzTarget <- peaklist[i, 8]
    xSpec <- which((abs(spectraScan[, 1] - mzTarget) <= massAccuracy) &
                     spectraScan[, 3] >= (peaklist[i, 1] - scanTolerance) &
                     spectraScan[, 3] <= (peaklist[i, 2] + scanTolerance))
    LxSpec <- length(xSpec)
    if (LxSpec >= minNIonPair) {
      if (LxSpec > 1) {
        ## To remove repeated scans
        xSpec <- IPA_aggregate(idVec = spectraScan[xSpec, 3], variableVec = spectraScan[xSpec, 1], indexVec = xSpec, targetVar = mzTarget)
        LxSpec <- length(xSpec)
      }
      if (LxSpec >= minNIonPair) {
        spectraScanRecursive <- spectraScan[xSpec, ]
        if (LxSpec == 1) {
          spectraScanRecursive <- matrix(spectraScan , nrow = 1)
        } else {
          spectraScanRecursive <- spectraScanRecursive[order(spectraScanRecursive[, 3], decreasing = FALSE), ]
        }
        ##
        scanNumberStart <- peaklist[i, 1] - scanTolerance
        if (scanNumberStart < 1) {
          scanNumberStart <- 1
        }
        scanNumberEnd <- peaklist[i, 2] + scanTolerance
        if (scanNumberEnd > LretentionTime) {
          scanNumberEnd <- LretentionTime
        }
        ##
        chromatographicPeakAnalysis(spectraScanRecursive, aggregatedSpectraList, retentionTime, LretentionTime, massAccuracy, mzTarget,
                                    rtTarget, scanNumberStart, scanNumberEnd, smoothingWindow, peakResolvingPower, minNIonPair,
                                    minPeakHeight, minRatioIonPair, maxRPW, minSNRbaseline, maxR13CcumulatedIntensity,
                                    maxPercentageMissingScans, nSpline, exportEICparameters)
      }
    }
  }
  ##
  ##############################################################################
  ##############################################################################
  ##
  if (number_processing_threads == 1) {
    ##
    peaklist_rec <- do.call(rbind, lapply(1:dim(peaklist)[1], function(i) {
      call_recursiveMZpeaklistCorrector(i)
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
      clusterExport(clust, setdiff(ls(), "clust"), envir = environment())
      ##
      peaklist_rec <- do.call(rbind, parLapply(clust, 1:dim(peaklist)[1], function(i) {
        call_recursiveMZpeaklistCorrector(i)
      }))
      ##
      stopCluster(clust)
      ##
      ##########################################################################
      ##
    } else {
      ##
      peaklist_rec <- do.call(rbind, mclapply(1:dim(peaklist)[1], function(i) {
        call_recursiveMZpeaklistCorrector(i)
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
  return(peaklist_rec)
}
