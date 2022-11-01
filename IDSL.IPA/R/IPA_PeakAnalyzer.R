IPA_PeakAnalyzer <- function(PARAM) {
  ##
  IPA.dir.create <- function(folder) {
    unlink(folder, recursive = TRUE)
    dir.create(folder, recursive = TRUE)
  }
  ##
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  peaklist024 <- matrix(rep(NA, 24), ncol = 24)
  ##
  number_processing_threads <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0006'), 2])
  input_path_hrms <- PARAM[which(PARAM[, 1] == 'PARAM0007'), 2]
  if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]) == "all") {
    file_name_hrms <- dir(path = input_path_hrms)
    file_name_hrms <- file_name_hrms[grep(paste0(".", tolower(PARAM[which(PARAM[, 1] == 'PARAM0009'), 2]), "$"), file_name_hrms, ignore.case = TRUE)]
  } else {
    samples_string <- PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]
    file_name_hrms <- strsplit(samples_string, ";")[[1]] # files used as reference m/z-RT
  }
  ##
  if (length(file_name_hrms) == 0) {
    stop(IPA_logRecorder("EMPTY HRMS FOLDER!!!"))
  }
  ##
  output_path <- PARAM[which(PARAM[, 1] == 'PARAM0010'), 2]
  output_path_peaklist <- paste0(output_path, "/peaklists")
  if (!dir.exists(output_path_peaklist)) {
    dir.create(output_path_peaklist, recursive = TRUE)
  }
  opendir(output_path_peaklist)
  ##
  exportEICcheck <- if (tolower(PARAM[which(PARAM[, 1] == 'PARAM_EIC'), 2]) == "yes") {TRUE} else {FALSE}
  if (exportEICcheck) {
    output_path_eics <- paste0(output_path, "/IPA_EIC")
    if (!dir.exists(output_path_eics)) {
      dir.create(output_path_eics, recursive = TRUE)
    }
    ##
    tryCatch(dev.off(), error = function(e) {Sys.sleep(0.0001)})
  }
  opendir(output_path_peaklist)
  ## To select monoisotopic peaks that have 13C isotopologues in the same scan
  minSpectraNoiseLevel <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0011'), 2])     # Intensity threshold in each scan
  ionMassDifference <- tryCatch(as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0012'), 2]), error = function(e) {1.003354835336}, warning = function(w) {1.003354835336})     # Mass difference for isotopic pairs
  #
  if (ionMassDifference <= 1.00336 & ionMassDifference >= 1.00335) {
    IPA_logRecorder("Carbon isotopes are selected for ion pairing!")
  } else {
    IPA_logRecorder(paste0("Mass difference to pair isotopes are = '", ionMassDifference, " Da'!"))
  }
  #
  massAccuracyXIC <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0013'), 2]) # Mass accuracy to cluster m/z in consecutive scans
  massAccuracy1.5 <- 1.5*massAccuracyXIC      # Mass accuracy to find 13C isotopologues
  RTtolerance <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0014'), 2])           	    # The retention time deviations to detect redundant peaks
  smoothingWindow <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0015'), 2])
  peakResolvingPower <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0017'), 2])

  ## HRMS size reduction
  powerSpectraReduction <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0018'), 2])

  # Data reduction
  recursiveMZcorrectionCheck <-  if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0019'), 2]) == "yes") {TRUE} else {FALSE}
  scanTolerance <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0020'), 2])		            # Number of scans to include in the search before and after of boundaries of detected peaks
  minPeakHeight <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0021'), 2])	    	  # Intensity threshold for peak height
  maxPercentageMissingScans <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0022'), 2])	  # Maximum percentage of missing scans on the raw chromatogram
  minNIonPair <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0023'), 2])          # Minimum number of scans in each peak
  minRatioIonPair <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0024'), 2])   # Ratio of number of detected scans per number of available scans (RCS)
  maxR13CcumulatedIntensity <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0025'), 2])  	      	  # Max ratio of 13C for the integrated spectra
  maxRPW <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0026'), 2])      	    	  # Peak width at half height to peak width at baseline (RPW)
  minSNRbaseline <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0027'), 2])     		        # S/N threshold
  nSpline <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0028'), 2])         		  # Level of peak smoothing to calculate ancillary chromatography parameters
  ##
  if (!recursiveMZcorrectionCheck) {
    nSpline1 <- nSpline
    exportEICcheck1 <- exportEICcheck
  } else {
    nSpline1 <- 0
    nSpline2 <- nSpline
    ##
    exportEICcheck1 <- FALSE
    exportEICcheck2 <- exportEICcheck
  }
  ##
  ##############################################################################
  ##
  call_carbon_IPA_parallel <- function(k) {
    ## To convert mzML/mzXML/CDF datafiles
    outputer <- IPA_MSdeconvoluter(input_path_hrms, file_name_hrms[k])
    spectraList <- outputer[["spectraList"]]
    RetentionTime <- outputer[["retentionTime"]]
    ## IPA_IonPairing
    spectraScan <- IPA_IonPairing(spectraList, minSpectraNoiseLevel, massAccuracy1.5, ionMassDifference)
    ## m/z clustering
    indexXIC <- MZclusteringRAWxic(spectraScan, massAccuracyXIC, minPeakHeight, minNIonPair)
    ## spectraList size reduction
    if (powerSpectraReduction > 0) {
      spectraList <- spectraListRoundingFiltering(spectraScan, spectraList, roundingDigit = powerSpectraReduction)
    }
    ##
    if (exportEICcheck1) {
      .GlobalEnv$counterEIC <- 0
      outputIPAeic <- paste0(output_path_eics, "/", file_name_hrms[k])
      exportEICparameters1 <- c(outputIPAeic, file_name_hrms[k], "UnTargetedWorkflow")
      ##
      IPA.dir.create(outputIPAeic)
      ##
    } else {
      exportEICparameters1 <- NULL
    }
    ## primary peak analyzer
    peaklist <- primary_peak_analyzer(spectraScan, indexXIC, scanTolerance, spectraList, RetentionTime, massAccuracyXIC,
                                      smoothingWindow, peakResolvingPower, minNIonPair, minPeakHeight, minRatioIonPair,
                                      maxRPW, minSNRbaseline, maxR13CcumulatedIntensity, maxPercentageMissingScans, nSpline1,
                                      exportEICparameters = exportEICparameters1)
    ## Recursive analysis
    if (recursiveMZcorrectionCheck) {
      if (!is.null(peaklist)) {
        ##
        if (exportEICcheck2) {
          .GlobalEnv$counterEIC <- 0
          outputIPAeic <- paste0(output_path_eics, "/", file_name_hrms[k])
          exportEICparameters2 <- c(outputIPAeic, file_name_hrms[k], "UnTargetedWorkflow")
          ##
          IPA.dir.create(outputIPAeic)
          ##
        } else {
          exportEICparameters2 <- NULL
        }
        ##
        peaklist <- recursive_mass_correction(peaklist, spectraScan, scanTolerance, spectraList, RetentionTime, massAccuracyXIC,
                                              smoothingWindow, peakResolvingPower, minNIonPair, minPeakHeight, minRatioIonPair,
                                              maxRPW, minSNRbaseline, maxR13CcumulatedIntensity, maxPercentageMissingScans, nSpline2,
                                              exportEICparameters = exportEICparameters2)
      }
    }
    ##
    if (!is.null(peaklist)) {
      if (length(peaklist) > 24) {
        ## Sort candidate m/z values by intensity
        orderINTpeaklist <- order(peaklist[, 4], decreasing = TRUE)
        ## To remove repeated peaks
        peaklistSubset <- cbind(peaklist[, 8], peaklist[, 3], peaklist[, 6])
        for (j in orderINTpeaklist) {
          x_mz_rt <- which((abs(peaklistSubset[, 1] - peaklistSubset[j, 1]) <= massAccuracyXIC) &
                             (abs(peaklistSubset[, 2] - peaklistSubset[j, 2]) <= RTtolerance))
          if (length(x_mz_rt) > 1) {
            x_max <- which.max(peaklistSubset[x_mz_rt, 3])
            xRemove <- x_mz_rt[-x_max[1]]
            peaklistSubset[xRemove, ] <- rep(0, 3)
          }
        }
        xEICpng  <- which(peaklistSubset[, 3] > 0)
        peaklist <- matrix(peaklist[xEICpng, ], ncol = 24)
        ## Sort candidate m/z values by the mass
        order12MZpeaklist <- order(peaklist[, 8], decreasing = FALSE)
        peaklist <- matrix(peaklist[order12MZpeaklist, ], ncol = 24)
      }
      ##
      peaklist[, 3] <- round(peaklist[, 3], 3) ## Retention Time Apex
      peaklist[, 4] <- round(peaklist[, 4], 0) ## Peak Height
      peaklist[, 5] <- round(peaklist[, 5], 0) ## Peak Area
      peaklist[, 7] <- round(peaklist[, 7], 0) ## RCS(%)
      peaklist[, 8] <- round(peaklist[, 8], 5) ## m/z 12C
      peaklist[, 9] <- round(peaklist[, 9], 0) ## Cumulated Intensity
      peaklist[, 10] <- round(peaklist[, 10], 5) ## m/z 13C
      peaklist[, 11] <- round(peaklist[, 11], 3) ## Ratio 13C Cumulated Intensity
      peaklist[, 12] <- round(peaklist[, 12], 3) ## Peak Width Baseline
      peaklist[, 13] <- round(peaklist[, 13], 3) ## Ratio Peak Width @ 50%
      peaklist[, 14] <- round(peaklist[, 14], 0) ## Separation Tray
      peaklist[, 15] <- round(peaklist[, 15], 3) ## Asymmetry Factor @ 10%
      peaklist[, 16] <- round(peaklist[, 16], 3) ## USP Tailing Factor @ 5%
      peaklist[, 17] <- round(peaklist[, 17], 3) ## Skewness Derivative Method
      peaklist[, 18] <- round(peaklist[, 18], 3) ## Symmetry Pseudo-Moments
      peaklist[, 19] <- round(peaklist[, 19], 3) ## Skewness Pseudo-Moments
      peaklist[, 20] <- round(peaklist[, 20], 3) ## Gaussianity
      peaklist[, 21] <- round(peaklist[, 21], 3) ## S/N baseline
      peaklist[, 22] <- round(peaklist[, 22], 3) ## S/N xcms method
      peaklist[, 23] <- round(peaklist[, 23], 3) ## S/N RMS
      peaklist[, 24] <- round(peaklist[, 24], 0) ## Sharpness
      ##
      if (exportEICcheck) {
        ##
        xEICpng <- xEICpng[order12MZpeaklist]
        ##
        PNGdir <- dir(path = outputIPAeic, pattern = ".png")
        strPNGdir <- strsplit(PNGdir, "_")
        ##
        idPNG <- do.call(c, lapply(strPNGdir, function(j) {
          j[length(j) - 5]
        }))
        idPNG <- as.numeric(idPNG)
        PNGdir <- PNGdir[order(idPNG, decreasing = FALSE)]
        ##
        jCounter <- 0
        for (j in xEICpng) {
          jCounter <- jCounter + 1
          oldEICname <- paste0(outputIPAeic, "/", PNGdir[j])
          newEICname <- paste0(outputIPAeic, "/IPA_EIC_", file_name_hrms[k], "_PeakID_", jCounter, "_MZ_", peaklist[jCounter, 8], "_RT_", peaklist[jCounter, 3], "_.png")
          file.rename(from = oldEICname, to = newEICname)
        }
        ##
        lRAWpeaklist <- length(orderINTpeaklist)
        xRemovePNG <- setdiff(seq(1, lRAWpeaklist, 1), xEICpng)
        if (length(xRemovePNG) > 0) {
          for (j in xRemovePNG) {
            unlink(paste0(outputIPAeic, "/", PNGdir[j]))
          }
        }
      }
      ##
    }else {
      peaklist <- peaklist024
    }
    ##
    colnames(peaklist) <- c("ScanNumberStart","ScanNumberEnd","RetentionTimeApex","PeakHeight","PeakArea",
                            "NumberDetectedScans(nIsoPair)","RCS(%)","m/z 12C","CumulatedIntensity",
                            "m/z 13C","Ratio 13C CumulatedIntensity","PeakWidthBaseline","Ratio PeakWidth @ 50%",
                            "SeparationTray","AsymmetryFactor @ 10%","USPTailingFactor @ 5%",
                            "Skewness_DerivativeMethod", "Symmetry PseudoMoments","Skewness PseudoMoments",
                            "Gaussianity", "S/N baseline", "S/N xcms method", "S/N RMS", "Sharpness")
    ##
    save(peaklist, file = paste0(output_path_peaklist, "/peaklist_", file_name_hrms[k], ".Rdata"))
    write.csv(peaklist, file = paste0(output_path_peaklist, "/peaklist_", file_name_hrms[k], ".csv"), row.names = FALSE)
    ##
    return()
  }
  ##
  IPA_logRecorder("Initiated HRMS peak detection!")
  IPA_logRecorder("Individual peaklists are stored in `.Rdata` and `.csv` formats in the `peaklist` folder!")
  if (exportEICcheck) {
    IPA_logRecorder("Extracted ion chromatogram (EIC) figures for the detected ions are stored in the `IPA_EICs` folder!")
    IPA_logRecorder("NOTICE: Please DO NOT open EIC figures and folders until the run is completed since at the end EIC figures are renamed according to their peak IDs in the peaklist files!")
  }
  ##
  if (number_processing_threads == 1) {
    ##
    progressBARboundaries <- txtProgressBar(min = 0, max = length(file_name_hrms), initial = 0, style = 3)
    ##
    for (k in 1:length(file_name_hrms)) {
      Null_variable <- tryCatch(call_carbon_IPA_parallel(k),
                                error = function(e) {save(peaklist024, file = paste0(output_path_peaklist, "/peaklist_", file_name_hrms[k], ".Rdata"))
                                  IPA_logRecorder(paste0("Problem with `", file_name_hrms[k],"`!"))})
      ##
      setTxtProgressBar(progressBARboundaries, k)
    }
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
        tryCatch(call_carbon_IPA_parallel(k),
                 error = function(e) {save(peaklist024, file = paste0(output_path_peaklist, "/peaklist_", file_name_hrms[k], ".Rdata"))
                   IPA_logRecorder(paste0("Problem with `", file_name_hrms[k],"`!"))})
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
      ##
    } else if (osType == "Windows") {
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      Null_variable <- foreach(k = 1:length(file_name_hrms), .combine = 'rbind', .verbose = FALSE) %dopar% {
        tryCatch(call_carbon_IPA_parallel(k),
                 error = function(e) {save(peaklist024, file = paste0(output_path_peaklist, "/peaklist_", file_name_hrms[k], ".Rdata"))
                   IPA_logRecorder(paste0("Problem with `", file_name_hrms[k],"`!"))})
      }
      stopCluster(clust)
      ##
    }
  }
  ##
  if (exportEICcheck) {
    opendir(output_path_eics)
  }
  ##
  IPA_logRecorder("Completed HRMS peak detection!")
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  return()
}
