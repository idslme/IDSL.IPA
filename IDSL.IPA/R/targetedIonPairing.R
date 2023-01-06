targetedIonPairing <- function(spectraList, scanNumberStart, scanNumberEnd, mzTarget, massAccuracy, ionMassDifference = 1.003354835336, massAccuracyIonPair = massAccuracy*1.5) {
  ##
  targetedIonPairs <- do.call(rbind, lapply(scanNumberStart:scanNumberEnd, function(t) {
    Spec <- spectraList[[t]]
    if (length(Spec) > 0) {
      x_mz1 <- which(abs(Spec[, 1] - mzTarget) <= massAccuracy)
      L_x_mz1 <- length(x_mz1)
      if (L_x_mz1 > 0) {
        mzTargetPair <- ionMassDifference + mzTarget
        x_mz2 <- which(abs(Spec[, 1] - mzTargetPair) <= massAccuracyIonPair)
        L_x_mz2 <- length(x_mz2)
        if (L_x_mz2 > 0) {
          if (L_x_mz1 > 1) {
            x_min <- which.min(abs(Spec[x_mz1, 1] - mzTarget))
            x_mz1 <- x_mz1[x_min[1]]
          }
          if (L_x_mz2 > 1) {
            x_min <- which.min(abs(Spec[x_mz2, 1] - mzTargetPair))
            x_mz2 <- x_mz2[x_min[1]]
          }
          if (Spec[x_mz1, 2] >= Spec[x_mz2, 2]) {
            c(Spec[x_mz1, 1], Spec[x_mz1, 2], t, Spec[x_mz2, 1], Spec[x_mz2, 2])
          }
        }
      }
    }
  }))
  ##
  return(targetedIonPairs)
}
