referenceRetentionTimeDetector <- function(inputPathPeaklist, refPeaklistFileNames, minFrequencyRefPeaks,
                                           massAccuracy, RTtolerance, number_processing_threads = 1) {
  L_sS <- length(refPeaklistFileNames)
  ##
  listRefRT <- lapply(refPeaklistFileNames, function(i) {
    RTvec <- loadRdata(paste0(inputPathPeaklist, "/", i))[, 3]
    if (is.na(RTvec[1])) {
      stop(IPA_logRecorder(paste0("IMPORTANT: `", i, "` CANNOT be a reference file in PARAM0030!")))
    }
    as.numeric(RTvec)
  })
  ##
  names(listRefRT) <- refPeaklistFileNames
  ##
  refPeakXcol <- peakAlignmentCore(inputPathPeaklist, refPeaklistFileNames, listRefRT,
                                   massAccuracy, RTtolerance, number_processing_threads)
  ##
  refPeakXcol[, 3] <- round(refPeakXcol[, 3]/L_sS*100, digits = 0)
  x_ref <- which(refPeakXcol[, 3] >= minFrequencyRefPeaks)
  mz_rt_Xmed_ref <- matrix(refPeakXcol[x_ref, 1:3], ncol = 3)
  ## To remove isomeric peaks
  round_mz <- round(mz_rt_Xmed_ref[, 1], digits = 1)
  x_unique <- which(table(round_mz) == 1)
  unique_mz <- as.numeric(names(x_unique))
  selectedMZ <- which(round_mz %in% unique_mz)
  referenceMZRTpeaks <- matrix(mz_rt_Xmed_ref[selectedMZ, ], ncol = 3)
  colnames(referenceMZRTpeaks) <- c("m/z", "RT", "freqRefPeaks(%)")
  rownames(referenceMZRTpeaks) <- NULL
  ##
  listReferencePeaks <- list(referenceMZRTpeaks, listRefRT)
  names(listReferencePeaks) <- c("referenceMZRTpeaks", "listRefRT")
  ##
  return(listReferencePeaks)
}
