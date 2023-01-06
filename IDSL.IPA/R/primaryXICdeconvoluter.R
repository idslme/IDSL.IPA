primaryXICdeconvoluter <- function(spectraScan, scanTolerance, indexXIC, aggregatedSpectraList, retentionTime,
                                   massAccuracy, smoothingWindow, peakResolvingPower, minNIonPair, minPeakHeight,
                                   minRatioIonPair, maxRPW, minSNRbaseline, maxR13CcumulatedIntensity,
                                   maxPercentageMissingScans, nSpline, exportEICparameters = NULL) {
  ##
  LretentionTime <- length(retentionTime)
  ##
  peaklist <- do.call(rbind, lapply(indexXIC, function(x) {
    ##
    mzTarget <- spectraScan[x[1], 1]
    ##
    Lx <- length(x)
    ##
    spectraScanPrimary <- spectraScan[x, ]
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
  }))
  ##
  return(peaklist)
}
