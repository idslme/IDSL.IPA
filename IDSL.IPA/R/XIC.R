XIC <- function(aggregatedSpectraList, scanNumberStart, scanNumberEnd, mzTarget, massAccuracy) {
  ##
  chromScanMzInt <- NULL
  LchromScanMzInt <- scanNumberEnd - scanNumberStart + 1
  minScanNumbers <- 3
  if (LchromScanMzInt >= minScanNumbers) {
    ##
    roundMZtarget <- round(mzTarget, digits = 2)
    ##
    xS <- do.call(c, lapply(c(-0.01, 0, 0.01), function(i) {
      aggregatedSpectraList[["aggregatedSpectraList"]][[as.character(roundMZtarget + i)]]
    }))
    ##
    if (length(xS) >= minScanNumbers) {
      xS <- sort(xS, decreasing = FALSE)
      ##
      spectraListMatrix.sb <- aggregatedSpectraList[["spectraListMatrix"]][xS, ]
      xS <- which((spectraListMatrix.sb[, 1] >= scanNumberStart) & (spectraListMatrix.sb[, 1] <= scanNumberEnd))
      if (length(xS) >= minScanNumbers) {
        spectraListMatrix.sb <- spectraListMatrix.sb[xS, ]
        ##
        xS <- which(abs(spectraListMatrix.sb[, 2] - mzTarget) <= massAccuracy)
        if (length(xS) >= minScanNumbers) {
          xS <- IPA_aggregate(idVec = spectraListMatrix.sb[xS, 1],
                              variableVec = spectraListMatrix.sb[xS, 2], indexVec = xS, targetVar = mzTarget)
          if (length(xS) >= minScanNumbers) {
            spectraListMatrix.sb <- spectraListMatrix.sb[xS, ]
            ##
            ####################################################################
            ##
            chromScanMzInt <- matrix(rep(0, LchromScanMzInt*3), ncol = 3) # c(Scan number, m/z, Intensity)
            chromScanMzInt[, 1] <- seq(scanNumberStart, scanNumberEnd, 1)
            ##
            ####################################################################
            ##
            xChromScanMzInt <- which(chromScanMzInt[, 1] %in% spectraListMatrix.sb[, 1])
            chromScanMzInt[xChromScanMzInt, 2] <- spectraListMatrix.sb[, 2]
            chromScanMzInt[xChromScanMzInt, 3] <- spectraListMatrix.sb[, 3]
            ##
            ####################################################################
            ##
          }
        }
      }
    }
  }
  ##
  return(chromScanMzInt)
}
