XIC <- function(spectraListXIC, scanNumberStart = 1, mzTarget, massAccuracyXIC) {
  chrom_ScN_Int <- do.call(rbind, lapply(1:length(spectraListXIC), function(t) {
    PEAKS <- spectraListXIC[[t]]
    PEAKSt <- c(t, 0, 0)
    x <- which(abs(PEAKS[ ,1] - mzTarget) <= massAccuracyXIC)
    L_x <- length(x)
    if (L_x > 0) {
      if (L_x > 1) {
        x2 <- which.min(abs(PEAKS[x, 1] - mzTarget))
        x <- x[x2[1]]
      }
      PEAKSt[2] <- PEAKS[x, 1]
      PEAKSt[3] <- PEAKS[x, 2]
    }
    PEAKSt # c(Scan number, m/z, Intensity)
  }))
  chrom_ScN_Int[, 1] <- chrom_ScN_Int[, 1] + scanNumberStart - 1
  return(chrom_ScN_Int)
}
