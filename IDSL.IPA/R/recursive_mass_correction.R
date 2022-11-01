recursive_mass_correction <- function(peaklist, spectraScan, scanTolerance, spectraList, RetentionTime, massAccuracyXIC,
                                      smoothingWindow, peakResolvingPower, minNIonPair, minPeakHeight, minRatioIonPair, maxRPW,
                                      minSNRbaseline, maxR13CcumulatedIntensity, maxPercentageMissingScans, nSpline, exportEICparameters = NULL) {
  n_RT <- length(RetentionTime)
  if (minNIonPair > 1) {
    peaklist_rec <- do.call(rbind, lapply(1:dim(peaklist)[1], function(k) {
      cc_k <- NULL
      rtTarget <- peaklist[k, 3]
      mzTarget <- peaklist[k, 8]
      x_Spec <- which((abs(spectraScan[, 1] - mzTarget) <= massAccuracyXIC) &
                        spectraScan[, 3] >= (peaklist[k, 1] - scanTolerance) &
                        spectraScan[, 3] <= (peaklist[k, 2] + scanTolerance))
      if (length(x_Spec) >= minNIonPair) {
        spectraScanRec <- spectraScan[x_Spec, ]
        spectraScanRec <- spectraScanRec[order(spectraScanRec[, 3], decreasing = FALSE), ]
        t1 <- peaklist[k, 1] - scanTolerance
        if (t1 < 1) {
          t1 <- 1
        }
        t2 <- peaklist[k, 2] + scanTolerance
        if (t2 > n_RT) {
          t2 <- n_RT
        }
        if (length(spectraScanRec) == 0 | spectraScanRec[1, 3] != t1) {
          spectraScanRec <- rbind(c(mzTarget, 0, t1, 0, 0), spectraScanRec)
        }
        if (length(spectraScanRec) == 5 | spectraScanRec[nrow(spectraScanRec),3] != t2) {
          spectraScanRec <- rbind(spectraScanRec, c(mzTarget, 0, t2, 0, 0))
        }
        cc_k <- chromatographyPeakAnalysis(spectraScanRec, smoothingWindow, peakResolvingPower,
                                           minNIonPair, minPeakHeight, minRatioIonPair, maxRPW,
                                           minSNRbaseline, maxR13CcumulatedIntensity, maxPercentageMissingScans,
                                           mzTarget, rtTarget, massAccuracyXIC, spectraList, RetentionTime, nSpline,
                                           exportEICparameters)
      }
      cc_k
    }))
  } else {
    peaklist_rec <- do.call(rbind, lapply(1:dim(peaklist)[1], function(k) {
      cc_k <- NULL
      rtTarget <- peaklist[k, 3]
      mzTarget <- peaklist[k, 8]
      x_Spec <- which((abs(spectraScan[, 1] - mzTarget) <= massAccuracyXIC) &
                        spectraScan[, 3] >= (peaklist[k, 1] - scanTolerance) &
                        spectraScan[, 3] <= (peaklist[k, 2] + scanTolerance))
      if (length(x_Spec) >= minNIonPair) {
        spectraScanRec <- matrix(spectraScan[x_Spec, ] , ncol = 5)
        spectraScanRec <- spectraScanRec[order(spectraScanRec[, 3], decreasing = FALSE), ]
        spectraScanRec <- matrix(spectraScanRec, ncol = 5)
        t1 <- peaklist[k, 1] - scanTolerance
        if (t1 < 1) {
          t1 <- 1
        }
        t2 <- peaklist[k, 2] + scanTolerance
        if (t2 > n_RT) {
          t2 <- n_RT
        }
        if (length(spectraScanRec) == 0 | spectraScanRec[1, 3] != t1) {
          spectraScanRec <- rbind(c(mzTarget, 0, t1, 0, 0), spectraScanRec)
        }
        if (length(spectraScanRec) == 5 | spectraScanRec[nrow(spectraScanRec), 3] != t2) {
          spectraScanRec <- rbind(spectraScanRec, c(mzTarget, 0, t2, 0, 0))
          spectraScanRec <- matrix(spectraScanRec, ncol = 5)
        }
        cc_k <- chromatographyPeakAnalysis(spectraScanRec, smoothingWindow, peakResolvingPower,
                                           minNIonPair, minPeakHeight, minRatioIonPair, maxRPW,
                                           minSNRbaseline, maxR13CcumulatedIntensity, maxPercentageMissingScans,
                                           mzTarget, rtTarget, massAccuracyXIC, spectraList, RetentionTime, nSpline,
                                           exportEICparameters)
      }
      cc_k
    }))
  }
  return(peaklist_rec)
}
