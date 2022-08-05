peak_Xcol2 <- function(input_path_peaklist, file_names_peaklist, peak_Xcol) {
  L_2 <- dim(peak_Xcol)[2]
  peak_height <- peak_Xcol
  peak_height[, 3:L_2] <- 0
  peak_area <- peak_height
  peak_R13C <- peak_height
  ##
  progressBARboundaries <- txtProgressBar(min = 3, max = L_2, initial = 3, style = 3)
  for (i in 3:L_2) {
    setTxtProgressBar(progressBARboundaries, i)
    peaklist <- loadRdata(paste0(input_path_peaklist, "/", file_names_peaklist[i - 2]))
    x_non0 <- which(peak_Xcol[, i] != 0)
    ##
    for (j in x_non0) {
      x_peak <- peak_Xcol[j, i]
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
