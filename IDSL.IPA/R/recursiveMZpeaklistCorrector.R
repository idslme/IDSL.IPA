recursiveMZpeaklistCorrector <- function(peaklist, spectraScan, scanTolerance, aggregatedSpectraList, retentionTime, massAccuracy,
                                         smoothingWindow, peakResolvingPower, minNIonPair, minPeakHeight, minRatioIonPair, maxRPW,
                                         minSNRbaseline, maxR13CcumulatedIntensity, maxPercentageMissingScans, nSpline, exportEICparameters = NULL) {
  ##
  LretentionTime <- length(retentionTime)
  ##
  peaklist_rec <- do.call(rbind, lapply(1:dim(peaklist)[1], function(k) {
    rtTarget <- peaklist[k, 3]
    mzTarget <- peaklist[k, 8]
    xSpec <- which((abs(spectraScan[, 1] - mzTarget) <= massAccuracy) &
                     spectraScan[, 3] >= (peaklist[k, 1] - scanTolerance) &
                     spectraScan[, 3] <= (peaklist[k, 2] + scanTolerance))
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
        scanNumberStart <- peaklist[k, 1] - scanTolerance
        if (scanNumberStart < 1) {
          scanNumberStart <- 1
        }
        scanNumberEnd <- peaklist[k, 2] + scanTolerance
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
  }))
  ##
  return(peaklist_rec)
}
