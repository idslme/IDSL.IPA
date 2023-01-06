peakXcolFlagger <- function(mzPeakXcol, rtPeakXcol, freqPeakXcol, massAccuracy, RTtolerance, maxRedundantPeakFlagging) {
  ##
  nPeaks <- length(mzPeakXcol)
  falggingVector <- rep(1, nPeaks)
  ##
  orderMZ <- order(mzPeakXcol, decreasing = FALSE)
  xDIffMZ <- c(0, which(diff(mzPeakXcol[orderMZ]) > massAccuracy), nPeaks)
  LxDIffMZ <- length(xDIffMZ) - 1
  ##
  for (q in 1:LxDIffMZ) {
    xQ <- seq((xDIffMZ[q] + 1), xDIffMZ[q + 1], 1)
    xQ <- orderMZ[xQ]
    ##
    xQ <- xQ[order(freqPeakXcol[xQ], decreasing = TRUE)]
    ##
    for (i in xQ) {
      if (mzPeakXcol[i] != 0) {
        xC <- which((abs(mzPeakXcol[xQ] - mzPeakXcol[i]) <= massAccuracy) &
                      (abs(rtPeakXcol[xQ] - rtPeakXcol[i]) <= RTtolerance))
        if (length(xC) > 1) {
          xC <- xQ[xC]
          xF <- which(freqPeakXcol[xC]/max(freqPeakXcol[xC]) <= maxRedundantPeakFlagging)
          if (length(xF) > 0) {
            ##
            falggingVector[xC[xF]] <- 0
            ##
          }
          mzPeakXcol[xC] <- 0
          rtPeakXcol[xC] <- 0
        }
      }
    }
  }
  ##
  return(falggingVector)
}
