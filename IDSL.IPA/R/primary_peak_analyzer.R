primary_peak_analyzer <- function (spec_scan, index_xic, scan_tol, spectraList, RetentionTime,
                                   mass_accuracy_xic, smoothing_window, peak_resolving_power, min_nIsoPair,
                                   min_peak_height, min_ratio_IsoPair, max_rpw, min_snr_baseline,
                                   max_R13C_integrated_peak, max_percentage_missing_scans, n_spline) {
  n_RT <- length(RetentionTime)
  if (min_nIsoPair > 1) {
    peaklist <- do.call(rbind, lapply(1:length(index_xic), function(i) {
      x <- index_xic[[i]]
      A <- spec_scan[x, ]
      mz_interim <- A[1, 1]
      A <- A[order(A[, 3]), ]
      t1 <- A[1, 3] - scan_tol
      if (t1 < 1) {
        t1 <- 1
      }
      t2 <- A[nrow(A), 3] + scan_tol
      if (t2 > n_RT) {
        t2 <- n_RT
      }
      if (A[1, 3] != t1) {
        A <- rbind(c(mz_interim, 0, t1, 0, 0), A)
      }
      if (A[nrow(A), 3] != t2) {
        A <- rbind(A, c(mz_interim, 0, t2, 0, 0))
      }
      chromatography_analysis(A, smoothing_window, peak_resolving_power, min_nIsoPair, min_peak_height,
                              min_ratio_IsoPair, max_rpw, min_snr_baseline, max_R13C_integrated_peak,
                              max_percentage_missing_scans, mz_interim, rt_target = 0,
                              mass_accuracy_xic, spectraList, RetentionTime, n_spline)
    }))
  } else {
    peaklist <- do.call(rbind, lapply(1:length(index_xic), function(i) {
      x <- index_xic[[i]]
      A <- matrix(spec_scan[x, ], ncol = 5)
      mz_interim <- A[1, 1]
      A <- A[order(A[, 3]), ]
      A <- matrix(A, ncol = 5)
      t1 <- A[1, 3] - scan_tol
      if (t1 < 1) {
        t1 <- 1
      }
      t2 <- A[nrow(A), 3] + scan_tol
      if (t2 > n_RT) {
        t2 <- n_RT
      }
      if (A[1, 3] != t1) {
        A <- rbind(c(mz_interim, 0, t1, 0, 0), A)
        A <- matrix(A, ncol = 5)
      }
      if (A[nrow(A), 3] != t2) {
        A <- rbind(A, c(mz_interim, 0, t2, 0, 0))
        A <- matrix(A, ncol = 5)
      }
      chromatography_analysis(A, smoothing_window, peak_resolving_power, min_nIsoPair, min_peak_height,
                              min_ratio_IsoPair, max_rpw, min_snr_baseline, max_R13C_integrated_peak,
                              max_percentage_missing_scans, mz_interim, rt_target = 0,
                              mass_accuracy_xic, spectraList, RetentionTime, n_spline)
    }))
  }
  return(peaklist)
}
