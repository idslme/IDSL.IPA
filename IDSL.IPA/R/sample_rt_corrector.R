sample_rt_corrector <- function(referenceMZRTpeaks, inputPathPeaklist, peaklistFileName, massAccuracy,
                                RTcorrectionMethod, refPeakTolerance = 1, degreePolynomial = 3) {
  ##
  listCorrectedRTpeaklists <- NULL
  ##
  peaklist <- loadRdata(paste0(inputPathPeaklist, "/", peaklistFileName))
  ##
  sample_referenceMZRTpeaks <- do.call(rbind, lapply(1:nrow(referenceMZRTpeaks), function(i) {
    x_peak <- which(abs(referenceMZRTpeaks[i, 1] - peaklist[, 8]) <= massAccuracy)
    if (length(x_peak) > 0) {
      if (length(x_peak) > 1) {
        x_min <- which.min(abs(referenceMZRTpeaks[i, 1] - peaklist[x_peak, 8]))
        x_peak <- x_peak[x_min[1]]
      }
      c(peaklist[x_peak, 8], peaklist[x_peak, 3], referenceMZRTpeaks[i, 2])
    }
  }))
  ##
  if (!is.null(sample_referenceMZRTpeaks)) {
    sample_referenceMZRTpeaks <- matrix(sample_referenceMZRTpeaks[order(sample_referenceMZRTpeaks[, 2]), ], ncol = 3)
    rtdf1 <- abs(sample_referenceMZRTpeaks[, 2] - sample_referenceMZRTpeaks[, 3])
    res <- boxplot(rtdf1, plot = FALSE)
    res3 <- res$out
    if (length(res3) > 2) {
      sample_referenceMZRTpeaks <- matrix(sample_referenceMZRTpeaks[!rtdf1 %in% res3, ], ncol = 3)
      nROW_sample_referenceMZRTpeaks <- nrow(sample_referenceMZRTpeaks)
      ##################### Retention Time Index ###############################
      if (gsub(" ", "", tolower(RTcorrectionMethod)) == "retentionindex") {
        ##
        if (refPeakTolerance > nROW_sample_referenceMZRTpeaks) {
          refPeakTolerance <- nROW_sample_referenceMZRTpeaks
          ##
          IPA_logRecorder(paste0("Number of reference peak tolerance for 'RetentionIndex' [PARAM0033] was reduced to ",
                                 refPeakTolerance, " for `", peaklistFileName, "` due to lack of sufficient reference peaks!"))
        }
        #
        RT_following_vector <- do.call(c, lapply((nROW_sample_referenceMZRTpeaks - refPeakTolerance + 1):nROW_sample_referenceMZRTpeaks, function(j) {
          RT_sample_following <- sample_referenceMZRTpeaks[j, 2]
          RT_ref_following <- sample_referenceMZRTpeaks[j, 3]
          (RT_ref_following - RT_sample_following)
        }))
        RT_diff_following <- median(RT_following_vector)
        #
        RT_preceding_vector <- do.call(c, lapply(1:refPeakTolerance, function(j) {
          RT_sample_preceding <- sample_referenceMZRTpeaks[j, 2]
          RT_ref_preceding <- sample_referenceMZRTpeaks[j, 3]
          (RT_ref_preceding - RT_sample_preceding)
        }))
        RT_diff_preceding <- median(RT_preceding_vector)
        #
        listCorrectedRTpeaklists <- sapply(1:nrow(peaklist), function(i) {
          RT_uncorrected <- peaklist[i, 3]
          #
          x_sample_preceding_main <- which((sample_referenceMZRTpeaks[, 2] - RT_uncorrected) <= 0)
          L_x_sample_preceding <- 0
          if (length(x_sample_preceding_main) > 0) {
            x_sample_preceding <- x_sample_preceding_main[length(x_sample_preceding_main)]
            x_sample_preceding <- seq((x_sample_preceding - refPeakTolerance + 1), x_sample_preceding, by = 1)
            x_sample_preceding <- x_sample_preceding[x_sample_preceding > 0]
            L_x_sample_preceding <- length(x_sample_preceding)
          }
          #
          x_sample_following_main <- which((sample_referenceMZRTpeaks[, 2] - RT_uncorrected) > 0)
          L_x_sample_following <- 0
          if (length(x_sample_following_main) > 0) {
            x_sample_following <- x_sample_following_main[1]
            x_sample_following <- seq(x_sample_following, (x_sample_following + refPeakTolerance - 1), by = 1)
            x_sample_following <- x_sample_following[x_sample_following <= nROW_sample_referenceMZRTpeaks]
            L_x_sample_following <- length(x_sample_following)
          }
          #
          if (L_x_sample_preceding > 0 & L_x_sample_following > 0) {
            #
            rt_intermediate_vector <- do.call(c, lapply(x_sample_preceding, function(j) {
              #
              RT_sample_preceding <- sample_referenceMZRTpeaks[j, 2]
              RT_ref_preceding <- sample_referenceMZRTpeaks[j, 3]
              #
              do.call(c, lapply(x_sample_following, function(k) {
                #
                RT_sample_following <- sample_referenceMZRTpeaks[k, 2]
                RT_ref_following <- sample_referenceMZRTpeaks[k, 3]
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
      ##################### Polynomial Regression ##############################
      if (gsub(" ", "", tolower(RTcorrectionMethod)) == "polynomial") {
        ##
        if (degreePolynomial > nROW_sample_referenceMZRTpeaks) {
          degreePolynomial <- nROW_sample_referenceMZRTpeaks - 1
          ##
          IPA_logRecorder(paste0("Polynomial degree for 'Polynomial' regression [PARAM0034] was reduced to ",
                                 degreePolynomial, " for `", peaklistFileName, "` due to lack of sufficient reference peaks!"))
        }
        ##
        idf <- data.frame(oRT = sample_referenceMZRTpeaks[, 2], eRT = sample_referenceMZRTpeaks[, 3])
        rtmodel <- lm(eRT ~ poly(oRT, degreePolynomial), idf)
        new.df <- data.frame(oRT = peaklist[, 3]) # predict new RTs
        listCorrectedRTpeaklists <- predict(rtmodel, new.df)
      }
      ##########################################################################
      listCorrectedRTpeaklists <- listCorrectedRTpeaklists
    }
  }
  ##
  if (is.null(listCorrectedRTpeaklists)) {
    listCorrectedRTpeaklists <- peaklist[, 3]
    IPA_logRecorder(paste0("Problem with retention time correction with `", peaklistFileName, "`!"))
  }
  ##
  return(listCorrectedRTpeaklists)
}
