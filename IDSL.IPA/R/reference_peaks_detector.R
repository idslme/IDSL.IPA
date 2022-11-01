reference_peaks_detector <- function(inputPathPeaklist, refPeaklistFileNames, minFrequencyRefPeaks,
                                     massAccuracy, RTtolerance, noQuantile, number_processing_threads = 1) {
  L_sS <- length(refPeaklistFileNames)
  ##
  listRefRT <- lapply(refPeaklistFileNames, function(i) {
    RTvec <- loadRdata(paste0(inputPathPeaklist, "/", i))[, 3]
    if (is.na(RTvec[1])) {
      stop(IPA_logRecorder(paste0("IMPORTANT: `", i, "' CANNOT be a reference file in PARAM0030!")))
    }
    RTvec
  })
  ##
  names(listRefRT) <- refPeaklistFileNames
  ##
  refPeakXcol <- peak_alignment(inputPathPeaklist, refPeaklistFileNames, listRefRT,
                                massAccuracy, RTtolerance, noQuantile, number_processing_threads)
  mz_rt_Xmed_ref <- do.call(rbind, lapply(1:nrow(refPeakXcol), function(i) {
    x_ref <- which(refPeakXcol[i, 3:(L_sS + 2)] != 0)
    freqRefPeaks <- round(length(x_ref)/L_sS, digits = 2)*100
    if (freqRefPeaks >= (minFrequencyRefPeaks)) {
      c(refPeakXcol[i, 1:2], freqRefPeaks)
    }
  }))
  ## To remove isomeric peaks
  round_mz <- round(mz_rt_Xmed_ref[, 1], digits = 1)
  x_unique <- which(table(round_mz) == 1)
  unique_mz <- as.numeric(names(x_unique))
  select_mz <- round_mz %in% unique_mz
  referenceMZRTpeaks <- matrix(mz_rt_Xmed_ref[select_mz, ], ncol = 3)
  colnames(referenceMZRTpeaks) <- c("m/z", "RT", "freqRefPeaks(%)")
  rownames(referenceMZRTpeaks) <- NULL
  ##
  listReferencePeaks <- list(referenceMZRTpeaks, listRefRT)
  names(listReferencePeaks) <- c("referenceMZRTpeaks", "listRefRT")
  ##
  return(listReferencePeaks)
}
