IPA_PeakAnalyzer <- function(PARAM) {
  ##
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  peaklistNA24 <- matrix(rep(NA, 24), nrow = 1)
  ##
  number_processing_threads <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0006'), 2])
  input_path_hrms <- PARAM[which(PARAM[, 1] == 'PARAM0007'), 2]
  samples_string <- PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]
  if (tolower(samples_string) == "all") {
    file_name_hrms <- dir(path = input_path_hrms)
    file_name_hrms <- file_name_hrms[grep(pattern = ".mzML$|.mzXML$|.CDF$", file_name_hrms, ignore.case = TRUE)]
  } else {
    file_name_hrms <- strsplit(samples_string, ";")[[1]]
  }
  LHRMS <- length(file_name_hrms)
  if (LHRMS == 0) {
    stop(IPA_logRecorder("ERROR!!! EMPTY HRMS FOLDER!!!"))
  }
  ##
  output_path <- PARAM[which(PARAM[, 1] == 'PARAM0010'), 2]
  output_path_peaklist <- paste0(output_path, "/peaklists")
  if (!dir.exists(output_path_peaklist)) {
    dir.create(output_path_peaklist, recursive = TRUE)
  }
  opendir(output_path_peaklist)
  ##
  exportEICcheck <- if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0009'), 2]) == "yes") {TRUE} else {FALSE}
  if (exportEICcheck) {
    ##
    IPA_dir.create <- function(folder) {
      tryCatch(unlink(folder, recursive = TRUE), error = function(e) {IPA_logRecorder(paste0("Can't delete `", folder, "`!"))})
      tryCatch(dir.create(folder, recursive = TRUE), warning = function(w) {NULL})
    }
    ##
    dev.offCheck <- TRUE
    while (dev.offCheck) {
      dev.offCheck <- tryCatch(dev.off(), error = function(e) {FALSE})
    }
    ##
    outputPathEICs <- paste0(output_path, "/IPA_EIC")
    if (!dir.exists(outputPathEICs)) {
      dir.create(outputPathEICs, recursive = TRUE)
    }
  }

  ## To select monoisotopic peaks that have 13C isotopologues in the same scan
  minSpectraNoiseLevel <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0011'), 2])         # Intensity threshold in each scan
  ionMassDifference <- tryCatch(as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0012'), 2]), error = function(e) {1.003354835336}, warning = function(w) {1.003354835336})     # Mass difference for ion pairs
  #
  if (ionMassDifference <= 1.00336 & ionMassDifference >= 1.00335) {
    IPA_logRecorder("Carbon isotopes are selected for ion pairing!")
  } else {
    IPA_logRecorder(paste0("Mass difference to pair ions is = `", ionMassDifference, " Da`!"))
  }
  #
  massAccuracy <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0013'), 2])                # Mass accuracy to cluster m/z in consecutive scans
  massAccuracy1.5 <- 1.5*massAccuracy                                                   # Mass accuracy to find 13C isotopologues
  RTtolerance <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0014'), 2])           	    # The retention time deviations to detect redundant peaks
  smoothingWindow <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0015'), 2])
  peakResolvingPower <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0017'), 2])

  ## Data reduction
  recursiveMZcorrectionCheck <- if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0019'), 2]) == "yes") {TRUE} else {FALSE}
  scanTolerance <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0020'), 2])		            # Number of scans to include in the search before and after of boundaries of detected peaks
  minPeakHeight <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0021'), 2])	    	        # Intensity threshold for peak height
  minPeakHeightXIC <- 0.20*minPeakHeight                                                # Intensity threshold for raw XIC peak clustering
  maxPercentageMissingScans <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0022'), 2])	  # Maximum percentage of missing scans on the raw chromatogram
  minNIonPair <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0023'), 2])                 # Minimum number of scans in each peak
  minRatioIonPair <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0024'), 2])             # Ratio of number of detected scans per number of available scans (RCS)
  maxR13CcumulatedIntensity <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0025'), 2])   # Max ratio of 13C for the integrated spectra
  maxRPW <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0026'), 2])      	          	  # Peak width at half height to peak width at baseline (RPW)
  minSNRbaseline <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0027'), 2])     		      # S/N threshold
  nSpline <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0028'), 2])         	      	  # Level of peak smoothing to calculate ancillary chromatography parameters
  ##
  if (nSpline == 0) {
    IPA_logRecorder("NOTICE: Ancillary chromatography parameters are not calculated because `PARAM0028` was selected `0`!")
  }
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
  call_IPA_PeakAnalyzer <- function(i) {
    ## To convert mzML/mzXML/CDF datafiles
    outputer <- IPA_MSdeconvoluter(input_path_hrms, file_name_hrms[i])
    spectraList <- outputer[["spectraList"]]
    retentionTime <- outputer[["retentionTime"]]
    outputer <- NULL
    ## Pending: This is where torch based smoothing, resolution and de-noising can happen. This needs a GPU hardware. How fast a CPU, depends
    ## IPA_IonPairing
    spectraScan <- IPA_IonPairing(spectraList, minSpectraNoiseLevel, massAccuracy1.5, ionMassDifference)
    ## m/z clustering
    indexXIC <- mzClusteringRawXIC(spectraScan123 = spectraScan[, 1:3], massAccuracy, minNIonPair, minPeakHeightXIC)
    ##
    aggregatedSpectraList <- IPA_spectraListAggregator(spectraList)
    spectraList <- NULL
    ##
    if (exportEICcheck1) {
      outputIPAeic <- paste0(outputPathEICs, "/", file_name_hrms[i])
      exportEICparameters <- c(outputIPAeic, file_name_hrms[i], "UnTargetedWorkflow")
      ##
      IPA_dir.create(outputIPAeic)
      ##
    } else {
      exportEICparameters <- NULL
    }
    ## primary peak analyzer
    peaklist <- primaryXICdeconvoluter(spectraScan, scanTolerance, indexXIC, aggregatedSpectraList, retentionTime, massAccuracy,
                                       smoothingWindow, peakResolvingPower, minNIonPair, minPeakHeight, minRatioIonPair,
                                       maxRPW, minSNRbaseline, maxR13CcumulatedIntensity, maxPercentageMissingScans, nSpline1,
                                       exportEICparameters)
    indexXIC <- NULL
    ## Recursive analysis
    if (recursiveMZcorrectionCheck) {
      if (!is.null(peaklist)) {
        ##
        if (exportEICcheck2) {
          outputIPAeic <- paste0(outputPathEICs, "/", file_name_hrms[i])
          exportEICparameters <- c(outputIPAeic, file_name_hrms[i], "UnTargetedWorkflow")
          ##
          IPA_dir.create(outputIPAeic)
          ##
        } else {
          exportEICparameters <- NULL
        }
        ##
        peaklist <- recursiveMZpeaklistCorrector(peaklist, spectraScan, scanTolerance, aggregatedSpectraList, retentionTime, massAccuracy,
                                                 smoothingWindow, peakResolvingPower, minNIonPair, minPeakHeight, minRatioIonPair,
                                                 maxRPW, minSNRbaseline, maxR13CcumulatedIntensity, maxPercentageMissingScans, nSpline2,
                                                 exportEICparameters)
      }
    }
    aggregatedSpectraList <- NULL
    spectraScan <- NULL
    ##
    ############################################################################
    ##
    if (!is.null(peaklist)) {
      ##
      peaklist <- matrix(peaklist, ncol = 24)
      nPeaks <- nrow(peaklist)
      if (nPeaks > 1) {
        ##
        ########################################################################
        ## Sort candidate m/z values by the mass
        peaklist <- peaklist[order(peaklist[, 8], decreasing = FALSE), ]
        ## To remove similar or redundant peaks
        xDiffQ <- c(0, which(diff(peaklist[, 8]) > massAccuracy), nPeaks)
        LxDiffQ <- length(xDiffQ) - 1
        ##
        for (q in 1:LxDiffQ) {
          xQ <- seq((xDiffQ[q] + 1), xDiffQ[q + 1], 1)
          ## To sort candidate m/z values by intensity
          xQ <- xQ[order(peaklist[xQ, 4], decreasing = TRUE)]
          ##
          for (j in xQ) {
            if (peaklist[j, 1] != 0) {
              xMZRT <- which((abs(peaklist[xQ, 8] - peaklist[j, 8]) <= massAccuracy) &
                               (abs(peaklist[xQ, 3] - peaklist[j, 3]) <= RTtolerance))
              if (length(xMZRT) > 1) {
                xMZRT <- xQ[xMZRT]
                xMax <- which.max(peaklist[xMZRT, 6])
                xRemove <- xMZRT[-xMax[1]]
                peaklist[xRemove, ] <- 0
              }
            }
          }
        }
        xNon0  <- which(peaklist[, 1] > 0)
        nPeaks <- length(xNon0)
        peaklist <- peaklist[xNon0, ]
        ##
        if (nPeaks == 1) {
          peaklist <- matrix(peaklist, nrow = 1)
        }
        ##
        ########################################################################
        ##
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
    } else {
      ##
      peaklist <- peaklistNA24
      nPeaks <- 1
    }
    ##
    colnames(peaklist) <- c("ScanNumberStart","ScanNumberEnd","retentionTimeApex","PeakHeight","PeakArea",
                            "NumberDetectedScans(nIsoPair)","RCS(%)","m/z 12C","CumulatedIntensity", "m/z 13C",
                            "Ratio13C CumulatedIntensity","PeakWidthBaseline","RatioPeakWidth @ 50%",
                            "SeparationTray","peakAsymmetryFactor @ 10%","peakUSPtailingFactor @ 5%",
                            "Skewness_DerivativeMethod", "Symmetry PseudoMoments","Skewness PseudoMoments",
                            "Gaussianity", "S/N baseline", "S/N xcms method", "S/N RMS", "Sharpness")
    ##
    save(peaklist, file = paste0(output_path_peaklist, "/peaklist_", file_name_hrms[i], ".Rdata"))
    write.csv(peaklist, file = paste0(output_path_peaklist, "/peaklist_", file_name_hrms[i], ".csv"), row.names = TRUE)
    ##
    ############################################################################
    ## To rename EICs according to their peak iDs in individual peaklists
    if (exportEICcheck & (!is.na(peaklist[1, 1]))) {
      ##
      cantRenameFiles <- do.call(c, lapply(1:nPeaks, function(j) {
        untargetedEICfilename <- paste0("IPA_EIC_", exportEICparameters[2], "_", exportEICparameters[3], "_MZ_", peaklist[j, 8], "_RT_", peaklist[j, 3], "_.png")
        oldEICfilename <- paste0(exportEICparameters[1], "/", untargetedEICfilename)
        newEICfilename <- paste0(exportEICparameters[1], "/IPA_EIC_", exportEICparameters[2], "_", j, "_MZ_", peaklist[j, 8], "_RT_", peaklist[j, 3], "_.png")
        renameCheck <- tryCatch(file.rename(from = oldEICfilename, to = newEICfilename), error = function(e) {FALSE}, warning = function(w) {FALSE})
        if (!renameCheck) {
          untargetedEICfilename
        }
      }))
      ##
      PNGdir <- dir(path = exportEICparameters[1], pattern = ".png$")
      ##
      if (!is.null(cantRenameFiles)) {
        PNGdir <- setdiff(PNGdir, cantRenameFiles)
      }
      ##
      xRemovePNG <- grep(exportEICparameters[3], PNGdir)
      if (length(xRemovePNG) > 0) {
        for (j in xRemovePNG) {
          tryCatch(unlink(paste0(exportEICparameters[1], "/", PNGdir[j])), error = function(e){NULL})
        }
      }
    }
    ##
    ############################################################################
    ##
    return()
  }
  ##
  ##############################################################################
  ##
  IPA_logRecorder("Initiated HRMS peak detection!")
  IPA_logRecorder("Individual peaklists are stored in `.Rdata` and `.csv` formats in the `peaklist` folder!")
  if (exportEICcheck) {
    IPA_logRecorder("Extracted ion chromatogram (EIC) figures for the detected ions are stored in the `IPA_EICs` folder!")
    IPA_logRecorder("NOTICE: Please DO NOT open EIC figures until the run is completed because EIC figures will be renamed according to their peak IDs in the peaklist files at the end of the run!")
  }
  ##
  ##############################################################################
  ##
  if (number_processing_threads == 1) {
    ##
    progressBARboundaries <- txtProgressBar(min = 0, max = LHRMS, initial = 0, style = 3)
    ##
    for (i in 1:LHRMS) {
      Null_variable <- tryCatch(call_IPA_PeakAnalyzer(i),
                                error = function(e) {save(peaklistNA24, file = paste0(output_path_peaklist, "/peaklist_", file_name_hrms[i], ".Rdata"))
                                  IPA_logRecorder(paste0("Problem with `", file_name_hrms[i],"`!"))})
      ##
      setTxtProgressBar(progressBARboundaries, i)
    }
    ##
    close(progressBARboundaries)
    ##
    ############################################################################
    ##
  } else {
    ## Processing OS
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Linux") {
      ##
      Null_variable <- do.call(rbind, mclapply(1:LHRMS, function(i) {
        tryCatch(call_IPA_PeakAnalyzer(i),
                 error = function(e) {save(peaklistNA24, file = paste0(output_path_peaklist, "/peaklist_", file_name_hrms[i], ".Rdata"))
                   IPA_logRecorder(paste0("Problem with `", file_name_hrms[i],"`!"))})
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
      ##
      ##########################################################################
      ##
    } else if (osType == "Windows") {
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      Null_variable <- foreach(i = 1:LHRMS, .combine = 'rbind', .verbose = FALSE) %dopar% {
        tryCatch(call_IPA_PeakAnalyzer(i),
                 error = function(e) {save(peaklistNA24, file = paste0(output_path_peaklist, "/peaklist_", file_name_hrms[i], ".Rdata"))
                   IPA_logRecorder(paste0("Problem with `", file_name_hrms[i],"`!"))})
      }
      stopCluster(clust)
      ##
    }
  }
  ##
  ##############################################################################
  ##
  if (exportEICcheck) {
    opendir(outputPathEICs)
  }
  ##
  IPA_logRecorder("Completed HRMS peak detection!")
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  ##
  ##############################################################################
  ##
  return()
}
