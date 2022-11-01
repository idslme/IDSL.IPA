peakXcolFiller <- function(peakXcol, inputPathPeaklist) {
  L_2 <- dim(peakXcol)[2]
  peak_height <- peakXcol
  peak_height[, 3:L_2] <- 0
  peak_area <- peak_height
  peak_R13C <- peak_height
  ##
  colnamesPeakXcol <- colnames(peakXcol)
  ##
  progressBARboundaries <- txtProgressBar(min = 3, max = L_2, initial = 3, style = 3)
  for (i in 3:L_2) {
    setTxtProgressBar(progressBARboundaries, i)
    iPeaklistFileName <- paste0("peaklist_", colnamesPeakXcol[i], ".Rdata")
    peaklist <- loadRdata(paste0(inputPathPeaklist, "/", iPeaklistFileName))
    x_non0 <- which(peakXcol[, i] != 0)
    ##
    for (j in x_non0) {
      x_peak <- peakXcol[j, i]
      peak_height[j, i] <- peaklist[x_peak, 4]
      peak_area[j, i] <- peaklist[x_peak, 5]
      peak_R13C[j, i] <- peaklist[x_peak, 11]
    }
  }
  close(progressBARboundaries)
  listHeightAreaR13C <- list(peak_height, peak_area, peak_R13C)
  names(listHeightAreaR13C) <- c("peak_height", "peak_area", "peak_R13C")
  return(listHeightAreaR13C)
}
