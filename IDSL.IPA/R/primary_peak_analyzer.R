primary_peak_analyzer <- function(spectraScan, indexXIC, scanTolerance, spectraList, RetentionTime,
                                  massAccuracyXIC, smoothingWindow, peakResolvingPower, minNIonPair,
                                  minPeakHeight, minRatioIonPair, maxRPW, minSNRbaseline,
                                  maxR13CcumulatedIntensity, maxPercentageMissingScans, nSpline,
                                  exportEICparameters = NULL) {
  n_RT <- length(RetentionTime)
  if (minNIonPair > 1) {
    peaklist <- do.call(rbind, lapply(1:length(indexXIC), function(i) {
      x <- indexXIC[[i]]
      A <- spectraScan[x, ]
      mz_interim <- A[1, 1]
      A <- A[order(A[, 3]), ]
      t1 <- A[1, 3] - scanTolerance
      if (t1 < 1) {
        t1 <- 1
      }
      t2 <- A[nrow(A), 3] + scanTolerance
      if (t2 > n_RT) {
        t2 <- n_RT
      }
      if (A[1, 3] != t1) {
        A <- rbind(c(mz_interim, 0, t1, 0, 0), A)
      }
      if (A[nrow(A), 3] != t2) {
        A <- rbind(A, c(mz_interim, 0, t2, 0, 0))
      }
      chromatographyPeakAnalysis(A, smoothingWindow, peakResolvingPower, minNIonPair, minPeakHeight,
                                 minRatioIonPair, maxRPW, minSNRbaseline, maxR13CcumulatedIntensity,
                                 maxPercentageMissingScans, mz_interim, rtTarget = NULL, massAccuracyXIC,
                                 spectraList, RetentionTime, nSpline, exportEICparameters)
    }))
  } else {
    peaklist <- do.call(rbind, lapply(1:length(indexXIC), function(i) {
      x <- indexXIC[[i]]
      A <- matrix(spectraScan[x, ], ncol = 5)
      mz_interim <- A[1, 1]
      A <- A[order(A[, 3]), ]
      A <- matrix(A, ncol = 5)
      t1 <- A[1, 3] - scanTolerance
      if (t1 < 1) {
        t1 <- 1
      }
      t2 <- A[nrow(A), 3] + scanTolerance
      if (t2 > n_RT) {
        t2 <- n_RT
      }
      if (A[1, 3] != t1) {
        A <- rbind(c(mz_interim, 0, t1, 0, 0), A)
        A <- matrix(A, ncol = 5)
      }
      if (A[nrow(A), 3] != t2) {
        A <- rbind(A, c(mz_interim, 0, t2, 0, 0))
        A <- matrix(A, ncol = 5)
      }
      chromatographyPeakAnalysis(A, smoothingWindow, peakResolvingPower, minNIonPair, minPeakHeight,
                                 minRatioIonPair, maxRPW, minSNRbaseline, maxR13CcumulatedIntensity,
                                 maxPercentageMissingScans, mz_interim, rtTarget = NULL, massAccuracyXIC,
                                 spectraList, RetentionTime, nSpline, exportEICparameters)
    }))
  }
  return(peaklist)
}
