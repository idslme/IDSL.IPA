IPA_PeakAnalyzer <- function (PARAM) {
  ##
  number_processing_cores <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0006'), 2])
  input_path_hrms <- PARAM[which(PARAM[, 1] == 'PARAM0007'), 2]
  if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]) == "all") {
    file_name_hrms <- dir(path = input_path_hrms)
    file_name_hrms <- file_name_hrms[grep(paste0(".", tolower(PARAM[which(PARAM[, 1] == 'PARAM0009'), 2]), "$"), file_name_hrms, ignore.case = TRUE)]
  } else {
    samples_string <- PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]
    file_name_hrms <- strsplit(samples_string, ";")[[1]] # files used as reference m/z-RT
  }
  output_path <- PARAM[which(PARAM[, 1] == 'PARAM0010'), 2]
  output_path_peaklist <- paste0(output_path, "/peaklists")
  if (!dir.exists(output_path_peaklist)) {
    dir.create(output_path_peaklist)
  }
  opendir(output_path_peaklist)
  ## To select monoisotopic peaks that have 13C isotopologues in the same scan
  int_threshold <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0011'), 2])     # Intensity threshold in each scan
  massDifferenceIsotopes <- tryCatch(as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0012'), 2]), error = function(e) {1.003354835336}, warning = function(w) {1.003354835336})     # Mass difference for isotopic pairs
  mass_accuracy_xic <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0013'), 2]) # Mass accuracy to cluster m/z in consecutive scans
  mass_accuracy_isotope_pair <- 1.5*mass_accuracy_xic      # Mass accuracy to find 13C isotopologues
  delta_rt <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0014'), 2])           	    # The retention time deviations to detect redundant peaks
  smoothing_window <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0015'), 2])
  peak_resolving_power <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0017'), 2])

  ## HRMS size reduction
  power_mass_spectral_size_reduction <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0018'), 2])

  # Data reduction
  mz_recursion <- PARAM[which(PARAM[, 1] == 'PARAM0019'), 2]
  scan_tol <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0020'), 2])		            # Number of scans to include in the search before and after of boundaries of detected peaks
  min_peak_height <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0021'), 2])	    	  # Intensity threshold for peak height
  max_percentage_missing_scans <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0022'), 2])	  # Maximum percentage of missing scans on the raw chromatogram
  min_nIsoPair <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0023'), 2])          # Minimum number of scans in each peak
  min_ratio_IsoPair <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0024'), 2])   # Ratio of number of detected scans per number of available scans (RCS)
  max_R13C_integrated_peak <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0025'), 2])  	      	  # Max ratio of 13C for the integrated spectra
  max_rpw <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0026'), 2])      	    	  # Peak width at half height to peak width at baseline (RPW)
  min_snr_baseline <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0027'), 2])     		        # S/N threshold
  n_spline <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0028'), 2])         		  # Level of peak smoothing to calculate ancillary chromatography parameters
  ##
  if (tolower(mz_recursion) != "yes") {
    n_spline1 <- n_spline
  } else {
    n_spline1 <- 0
    n_spline2 <- n_spline
  }
  ##
  call_carbon_IPA_parallel <- function(i) {
    ## To convert mzML/mzXML datafiles
    outputer <- IPA_MSdeconvoluter(input_path_hrms, file_name_hrms[i])
    spectraList <- outputer[[1]]
    RetentionTime <- outputer[[2]]
    ## IPA_isotope_pairing
    spec_scan <- IPA_isotope_pairing(spectraList, int_threshold, mass_accuracy_isotope_pair, massDifferenceIsotopes)
    ## m/z clustering
    index_xic <- mz_clustering_xic(spec_scan, mass_accuracy_xic, min_peak_height, min_nIsoPair)
    ## spectraList size reduction
    if (power_mass_spectral_size_reduction > 0) {
      spectraList <- spectraList_filtering(spec_scan, spectraList, rounding_digit = power_mass_spectral_size_reduction)
    }
    ## primary peak analyzer
    peaklist <- primary_peak_analyzer(spec_scan, index_xic, scan_tol, spectraList, RetentionTime, mass_accuracy_xic,
                                      smoothing_window, peak_resolving_power, min_nIsoPair, min_peak_height, min_ratio_IsoPair,
                                      max_rpw, min_snr_baseline, max_R13C_integrated_peak, max_percentage_missing_scans, n_spline1)
    ## Recursive analysis
    if (tolower(mz_recursion) == "yes") {
      if (!is.null(peaklist)) {
        peaklist <- recursive_mass_correction(peaklist, spec_scan, scan_tol, spectraList, RetentionTime, mass_accuracy_xic,
                                              smoothing_window, peak_resolving_power, min_nIsoPair, min_peak_height, min_ratio_IsoPair,
                                              max_rpw, min_snr_baseline, max_R13C_integrated_peak, max_percentage_missing_scans, n_spline2)
      }
    }
    if (!is.null(peaklist)) {
      if (length(peaklist) > 24) {
        ## To remove repeated peaks
        peaklist_subset <- cbind(peaklist[, 8], peaklist[, 3], peaklist[, 6])
        for (k in 1:nrow(peaklist_subset)) {
          x_mz_rt <- which(abs(peaklist_subset[, 1] - peaklist_subset[k, 1]) <= mass_accuracy_xic &
                             abs(peaklist_subset[, 2] - peaklist_subset[k, 2]) <= delta_rt)
          if (length(x_mz_rt) > 1) {
            x_max <- which.max(peaklist_subset[x_mz_rt, 3])
            x_remove <- x_mz_rt[-x_max[1]]
            peaklist_subset[x_remove, ] <- rep(0, 3)
          }
        }
        x_0  <- which(peaklist_subset[, 3] != 0)
        peaklist <- matrix(peaklist[x_0, ], ncol = 24)
        ##
        peaklist <- matrix(peaklist[order(peaklist[, 8]), ], ncol = 24) # Sort candidate m/z values by their masses
      }
      ##
      peaklist[, 1] <- round(peaklist[, 1], 0)
      peaklist[, 2] <- round(peaklist[, 2], 0)
      peaklist[, 3] <- round(peaklist[, 3], 3)
      peaklist[, 4] <- round(peaklist[, 4], 0)
      peaklist[, 5] <- round(peaklist[, 5], 0)
      peaklist[, 6] <- round(peaklist[, 6], 0)
      peaklist[, 7] <- round(peaklist[, 7], 0)
      peaklist[, 8] <- round(peaklist[, 8], 5)
      peaklist[, 9] <- round(peaklist[, 9], 0)
      peaklist[, 10] <- round(peaklist[, 10], 5)
      peaklist[, 11] <- round(peaklist[, 11], 3)
      peaklist[, 12] <- round(peaklist[, 12], 3)
      peaklist[, 13] <- round(peaklist[, 13], 3)
      peaklist[, 14] <- round(peaklist[, 14], 0)
      peaklist[, 15] <- round(peaklist[, 15], 3)
      peaklist[, 16] <- round(peaklist[, 16], 3)
      peaklist[, 17] <- round(peaklist[, 17], 3)
      peaklist[, 18] <- round(peaklist[, 18], 3)
      peaklist[, 19] <- round(peaklist[, 19], 3)
      peaklist[, 20] <- round(peaklist[, 20], 3)
      peaklist[, 21] <- round(peaklist[, 21], 3)
      peaklist[, 22] <- round(peaklist[, 22], 3)
      peaklist[, 23] <- round(peaklist[, 23], 3)
      peaklist[, 24] <- round(peaklist[, 24], 0)
      ##
    } else {
      peaklist <- matrix(rep(0, 24), ncol = 24)
    }
    peaklist <- data.frame(peaklist)
    colnames(peaklist) <- c("ScanNumberStart","ScanNumberEnd","RetentionTimeApex","PeakHeight","PeakArea",
                            "NumberDetectedScans(nIsoPair)","RCS(%)","m/z MonoIsotopic","CumulatedIntensity",
                            "m/z 13C","Ratio 13C CumulatedIntensity","PeakWidthBaseline","Ratio PeakWidth @ 50%",
                            "SeperationTray","AsymmetryFactor @ 10%","USPTailingFactor @ 5%",
                            "Skewness_DerivativeMethod", "Symmetry PseudoMoments","Skewness PseudoMoments",
                            "Gaussianity", "S/N baseline", "S/N xcms method", "S/N RMS", "Sharpness")
    ##
    save(peaklist, file = paste0(output_path_peaklist, "/peaklist_",file_name_hrms[i],".Rdata"))
    write.csv(peaklist, file = paste0(output_path_peaklist, "/peaklist_",file_name_hrms[i],".csv"))
    ##
    return()
  }
  print("Initiated HRMS peak detection!")
  ##
  if (number_processing_cores == 1) {
    ##
    progressBARboundaries <- txtProgressBar(min = 0, max = length(file_name_hrms), initial = 1, style = 3)
    ##
    Null_variable <- do.call(rbind, lapply(1:length(file_name_hrms), function(k) {
      call_carbon_IPA_parallel(k)
      ##
      setTxtProgressBar(progressBARboundaries, k)
    }))
    ##
    close(progressBARboundaries)
    ##
  } else {
    ## Processing OS
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Linux") {
      ##
      Null_variable <- do.call(rbind, mclapply(1:length(file_name_hrms), function(k) {
        call_carbon_IPA_parallel(k)
      }, mc.cores = number_processing_cores))
      ##
      closeAllConnections()
    } else if (osType == "Windows") {
      clust <- makeCluster(number_processing_cores)
      registerDoParallel(clust)
      ##
      Null_variable <- foreach(k = 1:length(file_name_hrms), .combine = 'rbind', .verbose = FALSE) %dopar% {
        call_carbon_IPA_parallel(k)
      }
      stopCluster(clust)
      ##
    }
  }
  ##
  print("Completed HRMS peak detection!")
  return()
}
