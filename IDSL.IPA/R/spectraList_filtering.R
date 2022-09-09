spectraList_filtering <- function(spec_scan.xic, spectraList, rounding_digit = 1) {
  ## To reduce the number of arrays in the spectraList variable
  N_chromatogramScans <- length(spectraList)
  ##
  spectraList.mat <- do.call(rbind, lapply(1:N_chromatogramScans, function(t) {
    cbind(spectraList[[t]], rep(t, dim(spectraList[[t]])[1]))
  }))
  selectIndex <- round(spectraList.mat[,1], digits = rounding_digit) %in% unique(round(spec_scan.xic[ ,1], digits = rounding_digit))
  spectraList.mat.sb <- spectraList.mat[selectIndex, ]
  ##############################################################################
  ## To aggregate peaks in the same chromatogram scans
  spectraList.mat.sb <- matrix(spectraList.mat.sb[order(spectraList.mat.sb[, 3], decreasing = FALSE), ], ncol = 3)
  ##
  xDiff <- which(diff(spectraList.mat.sb[, 3]) > 0)
  #
  xCSt <- matrix(rep(0, 2*N_chromatogramScans), ncol = 2)
  #
  uCS <- unique(spectraList.mat.sb[, 3])
  xCSt[uCS, 1] <- c(1, (xDiff + 1))
  xCSt[uCS, 2] <- c(xDiff, dim(spectraList.mat.sb)[1])
  ##
  spectra00 <- matrix(rep(0, 2), ncol = 2)
  ##
  spectraList.Reduced <- lapply(1:N_chromatogramScans, function(t) {
    ##
    if (xCSt[t, 1] != 0) {
      reducedSL <- matrix(spectraList.mat.sb[seq(xCSt[t, 1], xCSt[t, 2], 1), 1:2], ncol = 2)
    } else {
      reducedSL <- spectra00 # This must be here
    }
    return(reducedSL)
  })
  ##
  return(spectraList.Reduced)
}
