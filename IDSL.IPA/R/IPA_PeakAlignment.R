IPA_PeakAlignment <- function(PARAM) {
  ##
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  IPA_logRecorder("Initiated generating data for peak alignment!")
  ##
  number_processing_threads <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0006'), 2])
  ##
  input_path_hrms <- PARAM[which(PARAM[, 1] == 'PARAM0007'), 2]
  samples_string <- PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]
  if (tolower(samples_string) == "all") {
    file_name_hrms <- dir(path = input_path_hrms)
    file_name_hrms <- file_name_hrms[grep(pattern = ".mzML$|.mzXML$|.CDF$", file_name_hrms, ignore.case = TRUE)]
  } else {
    file_name_hrms <- strsplit(samples_string, ";")[[1]]
  }
  LHRMS <- length(file_name_hrms)
  ##
  output_path <- PARAM[which(PARAM[, 1] == 'PARAM0010'), 2]
  inputPathPeaklist <- paste0(output_path, "/peaklists")
  peaklistFileNames <- dir(path = inputPathPeaklist, pattern = ".Rdata$")
  peaklistFileNames <- peaklistFileNames[grep("^peaklist_", peaklistFileNames)]
  L_PL <- length(peaklistFileNames)
  ##
  if (LHRMS > L_PL) {
    peaklistHRMSfileNames <- paste0("peaklist_", file_name_hrms, ".Rdata")
    ndPeaklists <- setdiff(peaklistHRMSfileNames, peaklistFileNames)
    ndPeaklists <- gsub("^peaklist_|.Rdata$", "", ndPeaklists)
    IPA_logRecorder("Error!!! peaklist files are not available for the following HRMS file(s):")
    for (i in ndPeaklists) {
      IPA_logRecorder(i)
    }
    stop()
  }
  ##
  RTcorrectionCheck <- if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0029'), 2]) == "yes") {TRUE} else {FALSE}
  massAccuracy <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0035'), 2])
  RTtolerance <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0036'), 2])
  ##
  OutputPath_peak_alignment <- paste0(output_path, "/peak_alignment")
  if (!dir.exists(OutputPath_peak_alignment)) {
    dir.create(OutputPath_peak_alignment, recursive = TRUE)
  }
  ##
  IPA_logRecorder("Peak alignment datasets are stored in the `peak_alignment` folder!")
  ##
  if (RTcorrectionCheck) {
    reference_samples_string <- PARAM[which(PARAM[, 1] == 'PARAM0030'), 2]
    Ref_name <- strsplit(reference_samples_string, ";")[[1]] # files used as reference m/z-RT
    refPeaklistFileNames <- paste0("peaklist_", Ref_name, ".Rdata")
    ## To find common peaks in the reference samples
    IPA_logRecorder("Initiated detecting reference peaks for RT correction!")
    minFrequencyRefPeaks <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0031'), 2])
    ##
    listReferencePeaks <- referenceRetentionTimeDetector(inputPathPeaklist, refPeaklistFileNames, minFrequencyRefPeaks,
                                                         massAccuracy, RTtolerance, number_processing_threads)
    referenceMZRTpeaks <- listReferencePeaks[["referenceMZRTpeaks"]]
    listRefRT <- listReferencePeaks[["listRefRT"]]
    ##
    IPA_logRecorder("Reference endogenous peaks are stored for retention time correction! Please see the `endogenousReferenceMZRTpeaks.csv` in the `peak_alignment` folder!")
    write.csv(referenceMZRTpeaks, file = paste0(OutputPath_peak_alignment, "/endogenousReferenceMZRTpeaks.csv"), row.names = TRUE)
    ##
    IPA_logRecorder(paste0("Detected `", dim(referenceMZRTpeaks)[1], "` reference peaks for RT correction!"))
    #
    png(paste0(OutputPath_peak_alignment, "/Ref_peaks_distribution.png"), width = 20, height = 10, units = "in", res = 100)
    Ref_peaks_distribution <- referenceMZRTpeaks[, 2]
    histRTreferencePeaks <- hist(Ref_peaks_distribution, breaks = round(max(unlist(listRefRT))*1), xlab = "Retention time (min)")
    tryCatch(dev.off(), error = function(e) {NULL})
    L_x_regions_rt0 <- length(which(histRTreferencePeaks[["counts"]] == 0))
    if (L_x_regions_rt0 > 0) {
      IPA_logRecorder("WARNING!!! Reference peaks were not detected for the entire range of the retention times! Please see the `Ref_peaks_distribution.png` in the `peak_alignment` folder!")
    }
    ##
    IPA_logRecorder("Initiated retention time correction!")
    ##
    RTcorrectionMethod <- PARAM[which(PARAM[, 1] == 'PARAM0032'), 2]
    refPeakTolerance <- tryCatch(as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0033'), 2]), error = function(e) {5}, warning = function(w) {5})
    degreePolynomial <- tryCatch(as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0034'), 2]), error = function(e) {3}, warning = function(w) {3})
    ##
    peaklistFileNameSamples <- setdiff(peaklistFileNames, refPeaklistFileNames)
    ############################################################################
    if (number_processing_threads == 1) {
      listCorrectedRTsamples <- lapply(peaklistFileNameSamples, function(i) {
        analyteRetentionTimeCorrector(referenceMZRTpeaks, inputPathPeaklist, i, massAccuracy, RTcorrectionMethod, refPeakTolerance, degreePolynomial)
      })
    } else {
      ## Processing OS
      osType <- Sys.info()[['sysname']]
      if (osType == "Windows") {
        ##
        clust <- makeCluster(number_processing_threads)
        clusterExport(clust, c("referenceMZRTpeaks", "inputPathPeaklist", "massAccuracy", "RTcorrectionMethod", "refPeakTolerance", "degreePolynomial"), envir = environment())
        ##
        listCorrectedRTsamples <- parLapply(clust, peaklistFileNameSamples, function(i) {
          analyteRetentionTimeCorrector(referenceMZRTpeaks, inputPathPeaklist, i, massAccuracy, RTcorrectionMethod, refPeakTolerance, degreePolynomial)
        })
        ##
        stopCluster(clust)
        ##
      } else {
        ##
        listCorrectedRTsamples <- mclapply(peaklistFileNameSamples, function(i) {
          analyteRetentionTimeCorrector(referenceMZRTpeaks, inputPathPeaklist, i, massAccuracy, RTcorrectionMethod, refPeakTolerance, degreePolynomial)
        }, mc.cores = number_processing_threads)
        ##
        closeAllConnections()
        ##
      }
    }
    names(listCorrectedRTsamples) <- peaklistFileNameSamples
    ##
    ############################################################################
    ##
    listCorrectedRTpeaklists <- lapply(peaklistFileNames, function(i) {
      if (i %in% refPeaklistFileNames) {
        listRefRT[[i]]
      } else {
        listCorrectedRTsamples[[i]]
      }
    })
    ##
    IPA_logRecorder("Corrected retention times for individual peaklists are stored as `listCorrectedRTpeaklists.Rdata` in the `peak_alignment` folder!")
    IPA_logRecorder("Completed retention time correction!")
    ##
  } else {
    listCorrectedRTpeaklists <- lapply(peaklistFileNames, function(i) {
      as.numeric(loadRdata(paste0(inputPathPeaklist, "/", i))[, 3])
    })
  }
  ##
  names(listCorrectedRTpeaklists) <- peaklistFileNames
  ##
  save(listCorrectedRTpeaklists, file = paste0(OutputPath_peak_alignment, "/listCorrectedRTpeaklists.Rdata"))
  ##
  if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0002'), 2]) == "yes") {
    IPA_logRecorder("Initiated peak alignment on the entire peaklists using peak IDs in each peaklist!")
    peakXcol <- peakAlignmentCore(inputPathPeaklist, peaklistFileNames, listCorrectedRTpeaklists, massAccuracy, RTtolerance, number_processing_threads)
    listCorrectedRTpeaklists <- NULL
    colnames(peakXcol) <-  c("mz", "RT", "frequencyPeakXcol", file_name_hrms)
    ##
    save(peakXcol, file = paste0(OutputPath_peak_alignment, "/peakXcol.Rdata"))
    IPA_logRecorder("Stored aligned indexed table as 'peakXcol.Rdata' in the `peak_alignment` folder!")
    ##
    ############################################################################
    ##
    maxRedundantPeakFlagging <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0037'), 2])/100
    ##
    IPA_logRecorder("Initiated flagging suspicious aligned peaks!")
    ##
    falggingVector <- peakXcolFlagger(peakXcol[, 1], peakXcol[, 2], peakXcol[, 3], massAccuracy, 3*RTtolerance, maxRedundantPeakFlagging)
    ##
    IPA_logRecorder("Completed flagging suspicious aligned peaks!")
    ##
    ############################################################################
    ##
    IPA_logRecorder("Initiated generating and saving aligned peak tables for the peak height, peak area, and R13C values!")
    ##
    listHeightAreaR13C <- peakXcolFiller(peakXcol, inputPathPeaklist)
    peakXcol <- NULL
    ##
    peak_height <- peakPropertyTableMedianCalculator(listHeightAreaR13C[["peak_height"]], falggingVector, number_processing_threads)
    colnames(peak_height)[c(3, 4, 5)] <- c("freqPeakHeight", "medianPeakHeight", "Flag")
    save(peak_height, file = paste0(OutputPath_peak_alignment, "/peak_height.Rdata"))
    write.csv(peak_height, file = paste0(OutputPath_peak_alignment, "/peak_height.csv"), row.names = TRUE)
    ##
    peak_area <- peakPropertyTableMedianCalculator(listHeightAreaR13C[["peak_area"]], falggingVector, number_processing_threads)
    colnames(peak_area)[c(3, 4, 5)] <- c("freqPeakArea", "medianPeakArea", "Flag")
    save(peak_area, file = paste0(OutputPath_peak_alignment, "/peak_area.Rdata"))
    write.csv(peak_area, file = paste0(OutputPath_peak_alignment, "/peak_area.csv"), row.names = TRUE)
    peak_area <- NULL
    ##
    peak_R13C <- peakPropertyTableMedianCalculator(listHeightAreaR13C[["peak_R13C"]], falggingVector, number_processing_threads)
    colnames(peak_R13C)[c(3, 4, 5)] <- c("freqR13C", "medianR13C", "Flag")
    save(peak_R13C, file = paste0(OutputPath_peak_alignment, "/peak_R13C.Rdata"))
    write.csv(peak_R13C, file = paste0(OutputPath_peak_alignment, "/peak_R13C.csv"), row.names = TRUE)
    peak_R13C <- NULL
    ##
    listHeightAreaR13C <- NULL
    ##
    IPA_logRecorder("Aligned peak height, peak area, and R13C tables were stored in `.Rdata` and `.csv` formats in the `peak_alignment` folder!")
    ##
    ############################################################################
    ##
    if (LHRMS > 2) {
      IPA_logRecorder("Initiated detecting correlating peaks on the peak height table!")
      ##
      correlationMethod <- tolower(PARAM[which(PARAM[, 1] == 'PARAM_ALG1'), 2])
      minThresholdCorrelation <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM_ALG2'), 2])
      minFreqDetection <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM_ALG3'), 2])
      minRatioDetection <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM_ALG4'), 2])/100
      ##
      IPA_logRecorder(paste0("Initiated detecting correlating peaks on the aligned peak height table using the `", correlationMethod, "` method with coefficients `>=",
                             minThresholdCorrelation,"` and minimum number of complete observations `>=", minFreqDetection, "` combined with observation percentage `>=", minRatioDetection*100, "%`!"))
      ##
      alignedPeakHeightTableCorrelationList <- alignedPeakPropertyTableCorrelationListCalculator(peakPropertyTable = peak_height, RTtolerance, minFreqDetection, minRatioDetection,
                                                                                                 method = correlationMethod, minThresholdCorrelation, number_processing_threads)
      peak_height <- NULL
      ##
      save(alignedPeakHeightTableCorrelationList, file = paste0(OutputPath_peak_alignment, "/alignedPeakHeightTableCorrelationList.Rdata"))
      ##
      IPA_logRecorder("Stored the correlating peaks on the peak height table as `alignedPeakHeightTableCorrelationList.Rdata` in the `peak_alignement` folder (Only correlating aligned peak IDs are presented)!")
    }
    ##
    ############################################################################
    ##
    IPA_logRecorder("Completed generating data for peak alignment!")
    IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  }
  ##
  return()
}
