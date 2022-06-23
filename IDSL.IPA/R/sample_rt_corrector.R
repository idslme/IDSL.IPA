sample_rt_corrector <- function(reference_mz_rt_peaks, peaklist, mz_error,
                                rt_correction_method, reference_peak_tol = 1, polynomial_degree = 3) {
  ##
  sample_reference_mz_rt_peaks <- do.call(rbind, lapply(1:nrow(reference_mz_rt_peaks), function(i) {
    REF <- c()
    x_peak <- which(abs(reference_mz_rt_peaks[i, 1] - peaklist[, 8]) <= mz_error)
    if (length(x_peak) > 0) {
      if (length(x_peak) > 1) {
        x_min <- which.min(abs(reference_mz_rt_peaks[i, 1] - peaklist[x_peak, 8]))
        x_peak <- x_peak[x_min[1]]
      }
      REF <- c(peaklist[x_peak, 8], peaklist[x_peak, 3], reference_mz_rt_peaks[i, 2])
    }
    REF
  }))
  sample_reference_mz_rt_peaks <- sample_reference_mz_rt_peaks[order(sample_reference_mz_rt_peaks[, 2]), ]
  rtdf1 <- abs(sample_reference_mz_rt_peaks[, 2] - sample_reference_mz_rt_peaks[, 3])
  res3 <- boxplot(rtdf1, plot = FALSE)
  sample_reference_mz_rt_peaks <- sample_reference_mz_rt_peaks[!rtdf1 %in% res3$out, ]
  ##################### Retention Time Index ###################################
  if (gsub(" ", "", tolower(rt_correction_method)) == "retentionindex") {
    #
    nROW_sample_reference_mz_rt_peaks <- nrow(sample_reference_mz_rt_peaks)
    #
    RT_following_vector <- sapply((nROW_sample_reference_mz_rt_peaks - reference_peak_tol + 1):nROW_sample_reference_mz_rt_peaks, function(j) {
      RT_sample_following <- sample_reference_mz_rt_peaks[j, 2]
      RT_ref_following <- sample_reference_mz_rt_peaks[j, 3]
      (RT_ref_following - RT_sample_following)
    })
    RT_diff_following <- median(RT_following_vector)
    #
    RT_preceding_vector <- sapply(1:reference_peak_tol, function(j) {
      RT_sample_preceding <- sample_reference_mz_rt_peaks[j, 2]
      RT_ref_preceding <- sample_reference_mz_rt_peaks[j, 3]
      (RT_ref_preceding - RT_sample_preceding)
    })
    RT_diff_preceding <- median(RT_preceding_vector)
    #
    corrected_RT_peaklists <- sapply(1:nrow(peaklist), function(i) {
      RT_uncorrected <- peaklist[i, 3]
      #
      x_sample_preceding_main <- which((sample_reference_mz_rt_peaks[, 2] - RT_uncorrected) <= 0)
      L_x_sample_preceding <- 0
      if (length(x_sample_preceding_main) > 0) {
        x_sample_preceding <- x_sample_preceding_main[length(x_sample_preceding_main)]
        x_sample_preceding <- seq((x_sample_preceding - reference_peak_tol + 1), x_sample_preceding, by = 1)
        x_sample_preceding <- x_sample_preceding[x_sample_preceding > 0]
        L_x_sample_preceding <- length(x_sample_preceding)
      }
      #
      x_sample_following_main <- which((sample_reference_mz_rt_peaks[, 2] - RT_uncorrected) > 0)
      L_x_sample_following <- 0
      if (length(x_sample_following_main) > 0) {
        x_sample_following <- x_sample_following_main[1]
        x_sample_following <- seq(x_sample_following, (x_sample_following + reference_peak_tol - 1), by = 1)
        x_sample_following <- x_sample_following[x_sample_following <= nROW_sample_reference_mz_rt_peaks]
        L_x_sample_following <- length(x_sample_following)
      }
      #
      if (L_x_sample_preceding > 0 & L_x_sample_following > 0) {
        #
        rt_intermediate_vector <- do.call(c, lapply(x_sample_preceding, function(j) {
          #
          RT_sample_preceding <- sample_reference_mz_rt_peaks[j, 2]
          RT_ref_preceding <- sample_reference_mz_rt_peaks[j, 3]
          #
          do.call(c, lapply(x_sample_following, function(k) {
            #
            RT_sample_following <- sample_reference_mz_rt_peaks[k, 2]
            RT_ref_following <- sample_reference_mz_rt_peaks[k, 3]
            #
            exp(((log(RT_uncorrected) - log(RT_sample_preceding))/(log(RT_sample_following) - log(RT_sample_preceding)))*(log(RT_ref_following) - log(RT_ref_preceding)) + log(RT_ref_preceding))
          }))
        }))
        rt_intermediate <- median(rt_intermediate_vector)
        ##
      } else if (L_x_sample_following == 0) {
        rt_intermediate <- RT_uncorrected + RT_diff_following
        ##
      } else if (L_x_sample_preceding == 0) {
        rt_intermediate <- RT_uncorrected + RT_diff_preceding
      }
      ##
      rt_intermediate
    })
  }
  ##################### Polynomial Regression ##################################
  if (gsub(" ", "", tolower(rt_correction_method)) == "polynomial") {
    idf <- data.frame(oRT = sample_reference_mz_rt_peaks[, 2], eRT = sample_reference_mz_rt_peaks[, 3])
    rtmodel <- lm(eRT ~ poly(oRT, polynomial_degree), idf)
    new.df <- data.frame(oRT = peaklist[, 3]) # predict new RTs
    corrected_RT_peaklists <- predict(rtmodel, new.df)
  }
  ##############################################################################
  corrected_RT_peaklists <- matrix(corrected_RT_peaklists, ncol = 1)
  return(corrected_RT_peaklists)
}
