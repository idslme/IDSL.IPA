mzClusteringRawXIC <- function(spectraScan123, massAccuracy, minNIonPair, minPeakHeightXIC) {
  ##
  if (minPeakHeightXIC <= 0) {
    minPeakHeightXIC <- 1e-16
  }
  ##
  nPeaks <- nrow(spectraScan123)
  orderMZ <- order(spectraScan123[, 1], decreasing = FALSE)
  xDIffMZ <- c(0, which(diff(spectraScan123[orderMZ, 1]) > massAccuracy), nPeaks)
  LxDIffMZ <- length(xDIffMZ) - 1
  ##
  maxLengthClusters <- ceiling(nPeaks/max(c(1, minNIonPair)))
  indexXIC <- vector(mode = "list", maxLengthClusters)
  mzxic <- rep(0, maxLengthClusters)
  Counter <- 0
  ##
  for (q in 1:LxDIffMZ) {
    nQ <- xDIffMZ[q + 1] - xDIffMZ[q]
    if (nQ >= minNIonPair) {         # Minimum number of scans
      xQ <- seq((xDIffMZ[q] + 1), xDIffMZ[q + 1], 1)
      xQ <- orderMZ[xQ]
      ##
      xQ <- xQ[order(spectraScan123[xQ, 2], decreasing = TRUE)]
      ##
      for (i in xQ) {
        if (spectraScan123[i, 2] > minPeakHeightXIC) {
          x1 <- which(abs(spectraScan123[xQ, 1] - spectraScan123[i, 1]) <= massAccuracy)
          Lx1 <- length(x1)
          if (Lx1 >= minNIonPair) {         # Minimum number of scans
            x1 <- xQ[x1]
            if (Lx1 > 1) {
              ## To remove repeated scans
              x1 <- IPA_aggregate(idVec = spectraScan123[x1, 3], variableVec = spectraScan123[x1, 1],
                                  indexVec = x1, targetVar = spectraScan123[i, 1])
              Lx1 <- length(x1)
            }
            ##
            if (Lx1 >= minNIonPair) {         # Minimum number of scans
              Counter <- Counter + 1
              indexXIC[[Counter]] <- x1
              mzxic[Counter] <- spectraScan123[i, 1]
              spectraScan123[x1, ] <- 0
            }
          }
        }
      }
    }
  }
  ##
  indexXIC <- indexXIC[1:Counter]
  mzxic <- mzxic[1:Counter]
  names(indexXIC) <- mzxic
  ##
  return(indexXIC)
}
