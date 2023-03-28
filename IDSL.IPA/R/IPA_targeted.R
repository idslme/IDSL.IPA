IPA_targeted <- function(PARAM_targeted, allowedVerbose = TRUE) {
  ##
  peakPropertiesTable <- NULL
  ##
  if (!is.null(PARAM_targeted)) {
    ##
    mzCandidate <- tryCatch(eval(parse(text = paste0("c(", PARAM_targeted[which(PARAM_targeted[, 1] == 'PARAM_MZ'), 2], ")"))),
                            error = function(e){stop(IPA_logRecorder("Problem with 'PARAM_MZ'! This parameter can't make a vector of numbers!"))})
    rtCandidate <- tryCatch(eval(parse(text = paste0("c(", PARAM_targeted[which(PARAM_targeted[, 1] == 'PARAM_RT'), 2], ")"))),
                            error = function(e){stop(IPA_logRecorder("Problem with 'PARAM_RT'! This parameter can't make a vector of numbers!"))})
    ##
    lCandidate <- length(mzCandidate)
    if ((lCandidate == 0) | (lCandidate != length(rtCandidate))) {
      stop(IPA_logRecorder("Error!!! `PARAM_MZ` and `PARAM_RT` vectors do not have the same length!"))
    }
    xmzNA <- which(is.na(mzCandidate))
    if (length(xmzNA) > 0) {
      stop(IPA_logRecorder("Error!!! `PARAM_MZ` vector contains `NA` values! Maybe some odd characters are in the vector!"))
    }
    xrtNA <- which(is.na(rtCandidate))
    if (length(xrtNA) > 0) {
      stop(IPA_logRecorder("Error!!! `PARAM_RT` vector contains `NA` values! Maybe some odd characters are in the vector!"))
    }
    ##
    number_processing_threads <- as.numeric(PARAM_targeted[which(PARAM_targeted[, 1] == 'PARAM0006'), 2])
    ##
    if (number_processing_threads > 1) {
      parallelizationMode <- gsub(" ", "", tolower(PARAM_targeted[which(PARAM_targeted[, 1] == 'PARAM_PAR'), 2]))
      ##
      if (parallelizationMode != "peakmode") {
        parallelizationMode = "samplemode"
      }
    }
    ##
    input_path_hrms <- PARAM_targeted[which(PARAM_targeted[, 1] == 'PARAM0007'), 2]
    ##
    samples_string <- PARAM_targeted[which(PARAM_targeted[, 1] == 'PARAM0008'), 2]
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
    output_path <- PARAM_targeted[which(PARAM_targeted[, 1] == 'PARAM0010'), 2]
    if (!dir.exists(output_path)) {
      IPA_logRecorder("Created output directory!")
      dir.create(output_path, recursive = TRUE)
    }
    ##
    PARAM_EIC <- tolower(PARAM_targeted[which(PARAM_targeted[, 1] == 'PARAM0009'), 2])
    exportEICcheck <- if (!(PARAM_EIC == "no") | (PARAM_EIC == "n")) {TRUE} else {FALSE}
    if (exportEICcheck) {
      exportEICcheck <- TRUE
      outputPathEIC <- paste0(output_path, "/Targeted_EICs")
      if (!dir.exists(outputPathEIC)) {
        dir.create(outputPathEIC, recursive = TRUE)
      }
      if (allowedVerbose) {IPA_logRecorder("Extracted ion chromatograms (EICs) from targted workflow are stored in the `Targeted_EICs` folder!")}
      opendir(outputPathEIC)
      ##
      dev.offCheck <- TRUE
      while (dev.offCheck) {
        dev.offCheck <- tryCatch(dev.off(), error = function(e) {FALSE})
      }
      ##
    }
    ##
    PARAM_CCT <- tolower(PARAM_targeted[which(PARAM_targeted[, 1] == 'PARAM_CCT'), 2])
    exportTableCheck <- if (!(PARAM_CCT == "no") | (PARAM_CCT == "n")) {TRUE} else {FALSE}
    ##
    ionMassDifference <- tryCatch(as.numeric(PARAM_targeted[which(PARAM_targeted[, 1] == 'PARAM0012'), 2]), error = function(e) {1.003354835336}, warning = function(w) {1.003354835336})     # Mass difference for isotopic pairs
    massAccuracy <- as.numeric(PARAM_targeted[which(PARAM_targeted[, 1] == 'PARAM0013'), 2]) # Mass accuracy to cluster m/z in consecutive scans
    massAccuracy1.5 <- 1.5*massAccuracy      # Mass accuracy to find 13C isotopologues
    smoothingWindow <- as.numeric(PARAM_targeted[which(PARAM_targeted[, 1] == 'PARAM0015'), 2])
    peakResolvingPower <- as.numeric(PARAM_targeted[which(PARAM_targeted[, 1] == 'PARAM0017'), 2])
    scanTolerance <- as.numeric(PARAM_targeted[which(PARAM_targeted[, 1] == 'PARAM0020'), 2])		            # Number of scans to include in the search before and after of boundaries of detected peaks
    nSpline <- as.numeric(PARAM_targeted[which(PARAM_targeted[, 1] == 'PARAM0028'), 2])         		  # Level of peak smoothing to calculate ancillary chromatography parameters
    ##
    ############################################################################
    ##
    call_IPA_targeted <- function(j, jMZcandidate, jRTcandidate, spectraList, aggregatedSpectraList, retentionTime, LretentionTime,
                                  scanTolerance, ionMassDifference, massAccuracy, massAccuracy1.5, outputPathEIC, iFileNameHRMS,
                                  smoothingWindow, peakResolvingPower, nSpline) {
      ##
      ScanNumberApex <- which.min(abs(retentionTime - jRTcandidate))
      if (length(ScanNumberApex) > 0) {
        scanNumberStart <- ScanNumberApex - scanTolerance
        if (scanNumberStart < 1) {
          scanNumberStart <- 1
        }
        scanNumberEnd <- ScanNumberApex + scanTolerance
        if (scanNumberEnd > LretentionTime) {
          scanNumberEnd <- LretentionTime
        }
        ##
        chromatogramSegment <- targetedIonPairing(spectraList, scanNumberStart, scanNumberEnd, jMZcandidate, massAccuracy, ionMassDifference, massAccuracy1.5)
        ##
        if (is.null(chromatogramSegment)) {
          chromatogramSegment <- matrix(c(jMZcandidate, jMZcandidate, 0 , 0 , scanNumberStart, scanNumberEnd, 0, 0, 0, 0), ncol = 5)
        }
        ##
        if (exportEICcheck) {
          exportEICparameters <- c(outputPathEIC, iFileNameHRMS, j)
        } else {
          exportEICparameters = NULL
        }
        ##
        peak_property <- chromatographicPeakAnalysis(chromatogramSegment, aggregatedSpectraList, retentionTime, LretentionTime, massAccuracy,
                                                     mzTarget = jMZcandidate, rtTarget = jRTcandidate, scanNumberStart, scanNumberEnd, smoothingWindow,
                                                     peakResolvingPower, minNIonPair = 0, minPeakHeight = 0, minRatioIonPair = 0, maxRPW = 1, minSNRbaseline = 0,
                                                     maxR13CcumulatedIntensity = Inf, maxPercentageMissingScans = Inf, nSpline, exportEICparameters)
        ##
        c(iFileNameHRMS, jMZcandidate, jRTcandidate, peak_property)
      }
    }
    ##
    ############################################################################
    ##
    if (number_processing_threads == 1) {
      ##
      peakPropertiesTable <- do.call(rbind, lapply(1:LHRMS, function(i) {
        ##
        outputer <- IPA_MSdeconvoluter(input_path_hrms, file_name_hrms[i])
        spectraList <- outputer[["spectraList"]]
        retentionTime <- outputer[["retentionTime"]]
        outputer <- NULL
        LretentionTime <- length(retentionTime)
        aggregatedSpectraList <- IPA_spectraListAggregator(spectraList)
        ##
        do.call(rbind, lapply(1:lCandidate, function(j) {
          ##
          jMZcandidate <- mzCandidate[j]
          jRTcandidate <- rtCandidate[j]
          tryCatch(call_IPA_targeted(j, jMZcandidate, jRTcandidate, spectraList, aggregatedSpectraList, retentionTime,
                                     LretentionTime, scanTolerance, ionMassDifference, massAccuracy, massAccuracy1.5,
                                     outputPathEIC, iFileNameHRMS = file_name_hrms[i], smoothingWindow, peakResolvingPower, nSpline),
                   error = function(e) {IPA_logRecorder(paste0("Problem with `", file_name_hrms[i],"`!"))})
        }))
      }))
      ##
    } else {
      ##
      osType <- Sys.info()[['sysname']]
      ##
      if (parallelizationMode == "peakmode") {
        if (osType == "Windows") {
          ##
          peakPropertiesTable <- do.call(rbind, lapply(1:LHRMS, function(i) {
            ##
            outputer <- IPA_MSdeconvoluter(input_path_hrms, file_name_hrms[i])
            spectraList <- outputer[["spectraList"]]
            retentionTime <- outputer[["retentionTime"]]
            outputer <- NULL
            LretentionTime <- length(retentionTime)
            aggregatedSpectraList <- IPA_spectraListAggregator(spectraList)
            ##
            on.exit(tryCatch(stopCluster(clust), error = function(e){NULL}, warning = function(w){NULL}))
            ##
            clust <- makeCluster(number_processing_threads)
            clusterExport(clust, c("call_IPA_targeted", "i", "mzCandidate", "rtCandidate", "spectraList", "aggregatedSpectraList",
                                   "retentionTime", "LretentionTime", "scanTolerance", "ionMassDifference", "massAccuracy", "massAccuracy1.5",
                                   "outputPathEIC", "file_name_hrms", "smoothingWindow", "peakResolvingPower", "nSpline"), envir = environment())
            ##
            do.call(rbind, parLapply(clust, 1:lCandidate, function(j) {
              ##
              jMZcandidate <- mzCandidate[j]
              jRTcandidate <- rtCandidate[j]
              tryCatch(call_IPA_targeted(j, jMZcandidate, jRTcandidate, spectraList, aggregatedSpectraList, retentionTime,
                                         LretentionTime, scanTolerance, ionMassDifference, massAccuracy, massAccuracy1.5,
                                         outputPathEIC, iFileNameHRMS = file_name_hrms[i], smoothingWindow, peakResolvingPower, nSpline),
                       error = function(e) {IPA_logRecorder(paste0("Problem with `", file_name_hrms[i],"`!"))})
            }))
          }))
        } else {
          ##
          peakPropertiesTable <- do.call(rbind, lapply(1:LHRMS, function(i) {
            ##
            outputer <- IPA_MSdeconvoluter(input_path_hrms, file_name_hrms[i])
            spectraList <- outputer[["spectraList"]]
            retentionTime <- outputer[["retentionTime"]]
            outputer <- NULL
            LretentionTime <- length(retentionTime)
            aggregatedSpectraList <- IPA_spectraListAggregator(spectraList)
            ##
            do.call(rbind, mclapply(1:lCandidate, function(j) {
              ##
              jMZcandidate <- mzCandidate[j]
              jRTcandidate <- rtCandidate[j]
              tryCatch(call_IPA_targeted(j, jMZcandidate, jRTcandidate, spectraList, aggregatedSpectraList, retentionTime,
                                         LretentionTime, scanTolerance, ionMassDifference, massAccuracy, massAccuracy1.5,
                                         outputPathEIC, iFileNameHRMS = file_name_hrms[i], smoothingWindow, peakResolvingPower, nSpline),
                       error = function(e) {IPA_logRecorder(paste0("Problem with `", file_name_hrms[i],"`!"))})
            }, mc.cores = number_processing_threads))
          }))
          ##
          closeAllConnections()
          ##
        }
        ##
      } else if (parallelizationMode == "samplemode") {
        ##
        if (osType == "Windows") {
          ##
          clust <- makeCluster(number_processing_threads)
          clusterExport(clust, setdiff(ls(), c("clust", "LHRMS")), envir = environment())
          ##
          peakPropertiesTable <- do.call(rbind, parLapply(clust, 1:LHRMS, function(i) {
            ##
            outputer <- IPA_MSdeconvoluter(input_path_hrms, file_name_hrms[i])
            spectraList <- outputer[["spectraList"]]
            retentionTime <- outputer[["retentionTime"]]
            outputer <- NULL
            LretentionTime <- length(retentionTime)
            aggregatedSpectraList <- IPA_spectraListAggregator(spectraList)
            ##
            do.call(rbind, lapply(1:lCandidate, function(j) {
              ##
              jMZcandidate <- mzCandidate[j]
              jRTcandidate <- rtCandidate[j]
              tryCatch(call_IPA_targeted(j, jMZcandidate, jRTcandidate, spectraList, aggregatedSpectraList, retentionTime,
                                         LretentionTime, scanTolerance, ionMassDifference, massAccuracy, massAccuracy1.5,
                                         outputPathEIC, iFileNameHRMS = file_name_hrms[i], smoothingWindow, peakResolvingPower, nSpline),
                       error = function(e) {IPA_logRecorder(paste0("Problem with `", file_name_hrms[i],"`!"))})
            }))
          }))
          ##
          stopCluster(clust)
          ##
        } else {
          ##
          peakPropertiesTable <- do.call(rbind, mclapply(1:LHRMS, function(i) {
            ##
            outputer <- IPA_MSdeconvoluter(input_path_hrms, file_name_hrms[i])
            spectraList <- outputer[["spectraList"]]
            retentionTime <- outputer[["retentionTime"]]
            outputer <- NULL
            LretentionTime <- length(retentionTime)
            aggregatedSpectraList <- IPA_spectraListAggregator(spectraList)
            ##
            do.call(rbind, lapply(1:lCandidate, function(j) {
              ##
              jMZcandidate <- mzCandidate[j]
              jRTcandidate <- rtCandidate[j]
              tryCatch(call_IPA_targeted(j, jMZcandidate, jRTcandidate, spectraList, aggregatedSpectraList, retentionTime,
                                         LretentionTime, scanTolerance, ionMassDifference, massAccuracy, massAccuracy1.5,
                                         outputPathEIC, iFileNameHRMS = file_name_hrms[i], smoothingWindow, peakResolvingPower, nSpline),
                       error = function(e) {IPA_logRecorder(paste0("Problem with `", file_name_hrms[i],"`!"))})
            }))
            ##
          }, mc.cores = number_processing_threads))
          ##
          closeAllConnections()
          ##
        }
      }
    }
    ##
    ############################################################################
    ##
    if (exportTableCheck) {
      if (length(peakPropertiesTable) > 0) {
        peakPropertiesTable <- matrix(peakPropertiesTable, ncol = 27)
        peakPropertiesTable <- data.frame(peakPropertiesTable)
        colnames(peakPropertiesTable) <- c("Name HRMS", "m/z candidate", "RT candidate",
                                           "ScanNumberStart","ScanNumberEnd","retentionTimeApex","PeakHeight","PeakArea",
                                           "NumberDetectedScans(nIsoPair)","RCS(%)","m/z 12C","CumulatedIntensity", "m/z 13C",
                                           "Ratio13C CumulatedIntensity","PeakWidthBaseline","RatioPeakWidth @ 50%",
                                           "SeparationTray","peakAsymmetryFactor @ 10%","peakUSPtailingFactor @ 5%",
                                           "Skewness_DerivativeMethod", "Symmetry PseudoMoments","Skewness PseudoMoments",
                                           "Gaussianity", "S/N baseline", "S/N xcms method", "S/N RMS", "Sharpness")
        rownames(peakPropertiesTable) <- NULL
      }
    } else {
      peakPropertiesTable <- NULL
    }
  }
  return(peakPropertiesTable)
}
