reference_peaks_detector <- function(input_path_peaklist, file_names_peaklist_ref, min_frequency_ref_peaks,
                                     mz_error, rt_tol, n_quantile, number_processing_threads = 1) {
  L_sS <- length(file_names_peaklist_ref)
  ##
  listRefRT <- lapply(file_names_peaklist_ref, function(i) {
    matrix(loadRdata(paste0(input_path_peaklist, "/", i))[, 3], ncol = 1)
  })
  ##
  names(listRefRT) <- file_names_peaklist_ref
  ##
  peak_table_Xcol_ref <- peak_alignment(input_path_peaklist, file_names_peaklist_ref, listRefRT,
                                        mz_error, rt_tol, n_quantile, number_processing_threads)
  mz_rt_Xmed_ref <- do.call(rbind, lapply(1:nrow(peak_table_Xcol_ref), function(i) {
    x_ref <- which(peak_table_Xcol_ref[i, 3:(L_sS + 2)] != 0)
    if ((length(x_ref)/L_sS) >= (min_frequency_ref_peaks/100)) {
      peak_table_Xcol_ref[i, 1:2]
    }
  }))
  ## To remove isomeric peaks
  round_mz <- round(mz_rt_Xmed_ref[, 1], digits = 1)
  x_unique <- which(table(round_mz) == 1)
  unique_mz <- as.numeric(names(x_unique))
  select_mz <- round_mz %in% unique_mz
  reference_mz_rt_peaks <- matrix(mz_rt_Xmed_ref[select_mz, ], ncol = 2)
  ##
  listReferencePeaks <- list(reference_mz_rt_peaks, listRefRT)
  names(listReferencePeaks) <- c("referenceMZRTpeaks", "listRefRT")
  ##
  return(listReferencePeaks)
}
