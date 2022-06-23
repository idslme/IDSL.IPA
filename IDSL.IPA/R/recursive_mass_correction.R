recursive_mass_correction <- function(peaklist, spec_scan, scan_tol, spectraList, RetentionTime,
                                      mass_accuracy_xic, smoothing_window, peak_resolving_power,
                                      min_nIsoPair, min_peak_height, min_ratio_IsoPair, max_rpw,
                                      min_snr_baseline, max_R13C_integrated_peak, max_percentage_missing_scans, n_spline) {
  n_RT <- length(RetentionTime)
  if (min_nIsoPair > 1) {
    peaklist_rec <- do.call(rbind, lapply(1:dim(peaklist)[1], function(k) {
      cc_k <- c()
      rt_target <- peaklist[k, 3]
      mz_target <- peaklist[k, 8]
      x_Spec <- which(abs(spec_scan[, 1] - mz_target) <= mass_accuracy_xic &
                        spec_scan[, 3] >= (peaklist[k, 1] - scan_tol) &
                        spec_scan[, 3] <= (peaklist[k, 2] + scan_tol))
      if (length(x_Spec) >= min_nIsoPair) {
        spec_scan_rec <- spec_scan[x_Spec, ]
        spec_scan_rec <- spec_scan_rec[order(spec_scan_rec[,3]), ]
        t1 <- peaklist[k, 1] - scan_tol
        if (t1 < 1) {
          t1 <- 1
        }
        t2 <- peaklist[k, 2] + scan_tol
        if (t2 > n_RT) {
          t2 <- n_RT
        }
        if (length(spec_scan_rec) == 0 || spec_scan_rec[1, 3] != t1) {
          spec_scan_rec <- rbind(c(mz_target, 0, t1, 0, 0), spec_scan_rec)
        }
        if (length(spec_scan_rec) == 5 || spec_scan_rec[nrow(spec_scan_rec),3] != t2) {
          spec_scan_rec <- rbind(spec_scan_rec, c(mz_target, 0, t2, 0, 0))
        }
        cc_k <- chromatography_analysis(spec_scan_rec, smoothing_window, peak_resolving_power,
                                        min_nIsoPair, min_peak_height, min_ratio_IsoPair, max_rpw,
                                        min_snr_baseline, max_R13C_integrated_peak, max_percentage_missing_scans,
                                        mz_target, rt_target, mass_accuracy_xic, spectraList, RetentionTime, n_spline)
      }
      cc_k
    }))
  } else {
    peaklist_rec <- do.call(rbind, lapply(1:dim(peaklist)[1], function(k) {
      cc_k <- c()
      rt_target <- peaklist[k, 3]
      mz_target <- peaklist[k, 8]
      x_Spec <- which(abs(spec_scan[, 1] - mz_target) <= mass_accuracy_xic &
                        spec_scan[, 3] >= (peaklist[k, 1] - scan_tol) &
                        spec_scan[, 3] <= (peaklist[k, 2] + scan_tol))
      if (length(x_Spec) >= min_nIsoPair) {
        spec_scan_rec <- matrix(spec_scan[x_Spec, ] , ncol = 5)
        spec_scan_rec <- spec_scan_rec[order(spec_scan_rec[, 3]), ]
        spec_scan_rec <- matrix(spec_scan_rec, ncol = 5)
        t1 <- peaklist[k, 1] - scan_tol
        if (t1 < 1) {
          t1 <- 1
        }
        t2 <- peaklist[k, 2] + scan_tol
        if (t2 > n_RT) {
          t2 <- n_RT
        }
        if (length(spec_scan_rec) == 0 || spec_scan_rec[1, 3] != t1) {
          spec_scan_rec <- rbind(c(mz_target, 0, t1, 0, 0), spec_scan_rec)
        }
        if (length(spec_scan_rec) == 5 || spec_scan_rec[nrow(spec_scan_rec), 3] != t2) {
          spec_scan_rec <- rbind(spec_scan_rec, c(mz_target, 0, t2, 0, 0))
          spec_scan_rec <- matrix(spec_scan_rec, ncol = 5)
        }
        cc_k <- chromatography_analysis(spec_scan_rec, smoothing_window, peak_resolving_power,
                                        min_nIsoPair, min_peak_height, min_ratio_IsoPair, max_rpw,
                                        min_snr_baseline, max_R13C_integrated_peak, max_percentage_missing_scans,
                                        mz_target, rt_target, mass_accuracy_xic, spectraList, RetentionTime, n_spline)
      }
      cc_k
    }))
  }
  return(peaklist_rec)
}
