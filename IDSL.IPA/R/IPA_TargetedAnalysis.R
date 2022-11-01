IPA_TargetedAnalysis <- function(spreadsheet, mzCandidate, rtCandidate, exportEIC = TRUE, exportTable = FALSE) {
  ##
  cc_table <- NULL
  ##
  lCandidate <- length(mzCandidate)
  ##
  if (lCandidate !=  length(rtCandidate)) {
    stop(IPA_logRecorder("Error!!! `mzCandidate` and `rtCandidate` vectors do not have the same length!"))
  }
  ##
  PARAM <- IPA_targeted_xlsxAnalyzer(spreadsheet)
  if (!is.null(PARAM)) {
    ##
    number_processing_threads <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0006'), 2])
    if (number_processing_threads > 1) {
      para_mode <- gsub(" ", "", tolower(PARAM[which(PARAM[, 1] == 'PARAM_PAR'), 2]))
      ##
      if (para_mode != "samplemode") {
        para_mode = "peakmode"
      }
    }
    input_path_hrms <- PARAM[which(PARAM[, 1] == 'PARAM0007'), 2]
    if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]) == "all") {
      file_name_hrms <- dir(path = input_path_hrms)
      file_name_hrms <- file_name_hrms[grep(paste0(".", tolower(PARAM[which(PARAM[, 1] == 'PARAM0009'), 2]), "$"), file_name_hrms, ignore.case = TRUE)]
    } else {
      samples_string <- PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]
      file_name_hrms <- strsplit(samples_string, ";")[[1]] # files used as reference m/z-RT
    }
    output_path <- PARAM[which(PARAM[, 1] == 'PARAM0010'), 2]
    if (!dir.exists(output_path)) {
      IPA_logRecorder("Created output directory!")
      dir.create(output_path, recursive = TRUE)
    }
    if (exportEIC) {
      output_path_eic <- paste0(output_path, "/Targeted_EICs")
      if (!dir.exists(output_path_eic)) {
        dir.create(output_path_eic, recursive = TRUE)
      }
      IPA_logRecorder("Extracted ion chromatograms (EICs) from targted workflow are stored in the `Targeted_EICs` folder!")
      opendir(output_path_eic)
      ##
      tryCatch(dev.off(), error = function(e) {Sys.sleep(0.0001)})
    }
    ##
    ionMassDifference <- tryCatch(as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0012'), 2]), error = function(e) {1.003354835336}, warning = function(w) {1.003354835336})     # Mass difference for isotopic pairs
    massAccuracyXIC <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0013'), 2]) # Mass accuracy to cluster m/z in consecutive scans
    massAccuracy1.5 <- 1.5*massAccuracyXIC      # Mass accuracy to find 13C isotopologues
    smoothingWindow <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0015'), 2])
    peakResolvingPower <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0017'), 2])
    scanTolerance <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0020'), 2])		            # Number of scans to include in the search before and after of boundaries of detected peaks
    nSpline <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0028'), 2])         		  # Level of peak smoothing to calculate ancillary chromatography parameters
    ##
    ############################################################################
    ##
    chrom_call <- function(j, jMZcandidate, jRTcandidate, spectraList, RetentionTime, nRT,
                           scanTolerance, ionMassDifference, massAccuracyXIC, massAccuracy1.5,
                           output_path_eic, iFileNameHRMS, smoothingWindow, peakResolvingPower, nSpline) {
      ##
      ScanNumberApex <- which.min(abs(RetentionTime - jRTcandidate))
      t1 <- ScanNumberApex - scanTolerance
      if (t1 < 1) {
        t1 <- 1
      }
      t2 <- ScanNumberApex + scanTolerance
      if (t2 > nRT) {
        t2 <- nRT
      }
      ##
      chromatogram_segment <- do.call(rbind, lapply(t1:t2, function(t) {
        Spec_ScN_j <- NULL
        Spec <- spectraList[[t]]
        if (length(Spec) > 0) {
          x_mz1 <- which(abs(Spec[, 1] - jMZcandidate) <= massAccuracyXIC)
          L_x_mz1 <- length(x_mz1)
          if (L_x_mz1 > 0) {
            x_mz2 <- which(abs(Spec[, 1] - (ionMassDifference + jMZcandidate)) <= massAccuracy1.5)
            L_x_mz2 <- length(x_mz2)
            if (L_x_mz2 > 0) {
              if (L_x_mz1 > 1) {
                x_min <- which.min(abs(Spec[x_mz1, 1] - jMZcandidate))
                x_mz1 <- x_mz1[x_min[1]]
              }
              if (L_x_mz2 > 1) {
                x_min <- which.min(abs(Spec[x_mz2, 1] - (ionMassDifference + jMZcandidate)))
                x_mz2 <- x_mz2[x_min[1]]
              }
              if (Spec[x_mz1, 2] >= Spec[x_mz2, 2]) {
                Spec_ScN_j <- c(Spec[x_mz1, 1], Spec[x_mz1, 2], t, Spec[x_mz2, 1], Spec[x_mz2, 2])
              }
            }
          }
        }
        Spec_ScN_j
      }))
      if (length(chromatogram_segment) != 0) {
        chromatogram_segment <- chromatogram_segment[order(chromatogram_segment[, 3]), ]
        chromatogram_segment <- matrix(chromatogram_segment, ncol = 5)
      } else {
        chromatogram_segment <- matrix(c(jMZcandidate, jMZcandidate, 0 , 0 , t1, t2, 0, 0, 0, 0), ncol = 5)
      }
      if (chromatogram_segment[1, 3] != t1) {
        chromatogram_segment <- rbind(c(jMZcandidate, 0, t1, 0, 0), chromatogram_segment)
      }
      if (chromatogram_segment[dim(chromatogram_segment)[1], 3] != t2) {
        chromatogram_segment <- rbind(chromatogram_segment, c(jMZcandidate, 0, t2, 0, 0))
      }
      ##
      if (exportEIC) {
        exportEICparameters <- c(output_path_eic, iFileNameHRMS, j)
      } else {
        exportEICparameters = NULL
      }
      ##
      peak_property <- chromatographyPeakAnalysis(chromatogram_segment, smoothingWindow, peakResolvingPower, minNIonPair = 0, minPeakHeight = 0,
                                                  minRatioIonPair = 0, maxRPW = 1, minSNRbaseline = 0, maxR13CcumulatedIntensity = Inf,
                                                  maxPercentageMissingScans = Inf, mzTarget = jMZcandidate, rtTarget = jRTcandidate,
                                                  massAccuracyXIC, spectraList, RetentionTime, nSpline, exportEICparameters)
      ##
      c(iFileNameHRMS, jMZcandidate, jRTcandidate, peak_property)
    }
    ##
    ############################################################################
    ##
    if (number_processing_threads == 1) {
      ##
      cc_table <- do.call(rbind, lapply(1:length(file_name_hrms), function(i) {
        ##
        outputer <- IPA_MSdeconvoluter(input_path_hrms, file_name_hrms[i])
        spectraList <- outputer[["spectraList"]]
        RetentionTime <- outputer[["retentionTime"]]
        nRT <- length(RetentionTime)
        ##
        do.call(rbind, lapply(1:lCandidate, function(j) {
          ##
          chrom_call(j, jMZcandidate = mzCandidate[j], jRTcandidate = rtCandidate[j], spectraList, RetentionTime, nRT,
                     scanTolerance, ionMassDifference, massAccuracyXIC, massAccuracy1.5, output_path_eic,
                     iFileNameHRMS = file_name_hrms[i], smoothingWindow, peakResolvingPower, nSpline)
        }))
      }))
      ##
    } else {
      ##
      osType <- Sys.info()[['sysname']]
      ##
      if (para_mode == "peakmode") {
        if (osType == "Linux") {
          ##
          cc_table <- do.call(rbind, lapply(1:length(file_name_hrms), function(i) {
            ##
            outputer <- IPA_MSdeconvoluter(input_path_hrms, file_name_hrms[i])
            spectraList <- outputer[["spectraList"]]
            RetentionTime <- outputer[["retentionTime"]]
            nRT <- length(RetentionTime)
            ##
            do.call(rbind, mclapply(1:lCandidate, function(j) {
              ##
              chrom_call(j, jMZcandidate = mzCandidate[j], jRTcandidate = rtCandidate[j], spectraList, RetentionTime, nRT,
                         scanTolerance, ionMassDifference, massAccuracyXIC, massAccuracy1.5, output_path_eic,
                         iFileNameHRMS = file_name_hrms[i], smoothingWindow, peakResolvingPower, nSpline)
            }, mc.cores = number_processing_threads))
            ##
          }))
          ##
          closeAllConnections()
          ##
        } else if (osType == "Windows") {
          ##
          clust <- makeCluster(number_processing_threads)
          registerDoParallel(clust)
          ##
          cc_table <- do.call(rbind, lapply(1:length(file_name_hrms), function(i) {
            ##
            outputer <- IPA_MSdeconvoluter(input_path_hrms, file_name_hrms[i])
            spectraList <- outputer[["spectraList"]]
            RetentionTime <- outputer[["retentionTime"]]
            nRT <- length(RetentionTime)
            ##
            foreach(j = 1:lCandidate, .combine = 'rbind', .verbose = FALSE) %dopar% {
              ##
              chrom_call(j, jMZcandidate = mzCandidate[j], jRTcandidate = rtCandidate[j], spectraList, RetentionTime, nRT,
                         scanTolerance, ionMassDifference, massAccuracyXIC, massAccuracy1.5, output_path_eic,
                         iFileNameHRMS = file_name_hrms[i], smoothingWindow, peakResolvingPower, nSpline)
            }
            ##
          }))
          ##
          stopCluster(clust)
          ##
        }
        ##
      } else if (para_mode == "samplemode") {
        ##
        if (osType == "Linux") {
          ##
          cc_table <- do.call(rbind, mclapply(1:length(file_name_hrms), function(i) {
            ##
            outputer <- IPA_MSdeconvoluter(input_path_hrms, file_name_hrms[i])
            spectraList <- outputer[["spectraList"]]
            RetentionTime <- outputer[["retentionTime"]]
            nRT <- length(RetentionTime)
            ##
            do.call(rbind, lapply(1:lCandidate, function(j) {
              ##
              chrom_call(j, jMZcandidate = mzCandidate[j], jRTcandidate = rtCandidate[j], spectraList, RetentionTime, nRT,
                         scanTolerance, ionMassDifference, massAccuracyXIC, massAccuracy1.5, output_path_eic,
                         iFileNameHRMS = file_name_hrms[i], smoothingWindow, peakResolvingPower, nSpline)
            }))
            ##
          }, mc.cores = number_processing_threads))
          ##
          closeAllConnections()
          ##
        } else if (osType == "Windows") {
          ##
          clust <- makeCluster(number_processing_threads)
          registerDoParallel(clust)
          ##
          cc_table <- foreach(i = 1:length(file_name_hrms), .combine = 'rbind', .verbose = FALSE) %dopar% {
            ##
            outputer <- IPA_MSdeconvoluter(input_path_hrms, file_name_hrms[i])
            spectraList <- outputer[["spectraList"]]
            RetentionTime <- outputer[["retentionTime"]]
            nRT <- length(RetentionTime)
            ##
            do.call(rbind, lapply(1:lCandidate, function(j) {
              ##
              chrom_call(j, jMZcandidate = mzCandidate[j], jRTcandidate = rtCandidate[j], spectraList, RetentionTime, nRT,
                         scanTolerance, ionMassDifference, massAccuracyXIC, massAccuracy1.5, output_path_eic,
                         iFileNameHRMS = file_name_hrms[i], smoothingWindow, peakResolvingPower, nSpline)
            }))
            ##
          }
          ##
          stopCluster(clust)
          ##
        }
      }
    }
    ##
    ############################################################################
    ##
    if (exportTable) {
      if (length(cc_table) > 0) {
        cc_table <- matrix(cc_table, ncol = 27)
        cc_table <- data.frame(cc_table)
        colnames(cc_table) <- c("Name HRMS", "m/z candidate", "RT candidate",
                                "ScanNumberStart","ScanNumberEnd","RetentionTimeApex","PeakHeight","PeakArea",
                                "NumberDetectedScans(nIsoPair)","RCS(%)","m/z 12C","CumulatedIntensity",
                                "m/z 13C","Ratio 13C CumulatedIntensity","PeakWidthBaseline","Ratio PeakWidth @ 50%",
                                "SeparationTray","AsymmetryFactor @ 10%","USPTailingFactor @ 5%",
                                "Skewness_DerivativeMethod", "Symmetry PseudoMoments","Skewness PseudoMoments",
                                "Gaussianity", "S/N baseline", "S/N xcms method", "S/N RMS", "Sharpness")
        rownames(cc_table) <- NULL
      }
    } else {
      cc_table <- NULL
    }
  }
  return(cc_table)
}
