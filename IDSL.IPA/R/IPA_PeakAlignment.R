IPA_PeakAlignment <- function(PARAM) {
  ##
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  IPA_logRecorder("Initiated generating data for peak alignment!")
  ##
  number_processing_threads <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0006'), 2])
  output_path <- PARAM[which(PARAM[, 1] == 'PARAM0010'), 2]
  inputPathPeaklist <- paste0(output_path, "/peaklists")
  peaklistFileNames1 <- dir(path = inputPathPeaklist, pattern = ".Rdata")
  peaklistFileNames2 <- dir(path = inputPathPeaklist, pattern = "peaklist_")
  peaklistFileNames <- peaklistFileNames1[peaklistFileNames1 %in% peaklistFileNames2]
  L_PL <- length(peaklistFileNames)
  ##
  input_path_hrms <- PARAM[which(PARAM[, 1] == 'PARAM0007'), 2]
  if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]) == "all") {
    file_names_hrms <- dir(path = input_path_hrms)
    file_names_hrms <- file_names_hrms[grep(paste0(".", tolower(PARAM[which(PARAM[, 1] == 'PARAM0009'), 2]), "$"), file_names_hrms, ignore.case = TRUE)]
  } else {
    samples_string <- PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]
    file_names_hrms <- strsplit(samples_string, ";")[[1]] # files used as reference m/z-RT
  }
  ##
  peaklist2HRMS <- gsub("^peaklist_", "", gsub(".Rdata$", "", peaklistFileNames))
  matchPeaklist2HRMS <- file_names_hrms[(file_names_hrms %in% peaklist2HRMS)]
  if (length(matchPeaklist2HRMS) != L_PL) {
    stop(IPA_logRecorder("Error!!! peaklist files are not available for all selected HRMS files!"))
  }
  ##
  RTcorrectionCheck <- if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0029'), 2]) == "yes") {TRUE} else {FALSE}
  massAccuracy <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0035'), 2])
  RTtolerance <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0036'), 2])
  noQuantile <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0037'), 2])
  ##
  OutputPath_peak_alignment <- paste0(output_path, "/peak_alignment")
  if (!dir.exists(OutputPath_peak_alignment)) {
    dir.create(OutputPath_peak_alignment, recursive = TRUE)
  }
  ##
  IPA_logRecorder("Peak alignment data are stored in the `peak_alignment` folder!")
  ##
  if (RTcorrectionCheck) {
    reference_samples_string <- PARAM[which(PARAM[, 1] == 'PARAM0030'), 2]
    Ref_name <- strsplit(reference_samples_string, ";")[[1]] # files used as reference m/z-RT
    refPeaklistFileNames <- paste0("peaklist_", Ref_name, ".Rdata")
    ## To find common peaks in the reference samples
    IPA_logRecorder("Initiated detecting reference peaks for RT correction!")
    minFrequencyRefPeaks <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0031'), 2])
    ##
    listReferencePeaks <- reference_peaks_detector(inputPathPeaklist, refPeaklistFileNames, minFrequencyRefPeaks,
                                                   massAccuracy, RTtolerance, noQuantile, number_processing_threads)
    referenceMZRTpeaks <- listReferencePeaks[["referenceMZRTpeaks"]]
    listRefRT <- listReferencePeaks[["listRefRT"]]
    ##
    IPA_logRecorder("Reference endogenous peaks are stored for retention time correction! Please see the 'referenceMZRTpeaks.csv' in the `peak_alignment` folder!")
    write.csv(referenceMZRTpeaks, file = paste0(OutputPath_peak_alignment, "/referenceMZRTpeaks.csv"), row.names = FALSE)
    ##
    IPA_logRecorder(paste0("Detected " , dim(referenceMZRTpeaks)[1], " reference peaks for RT correction!"))
    #
    png(paste0(OutputPath_peak_alignment, "/Ref_peaks_distribution.png"), width = 20, height = 10, units = "in", res = 100)
    Ref_peaks_distribution <- referenceMZRTpeaks[, 2]
    hist_rt_reference_peaks <- hist(Ref_peaks_distribution, breaks = round(max(unlist(listRefRT))*1), xlab = "Retention time (min)")
    tryCatch(dev.off(), error = function(e) {Sys.sleep(0.0001)})
    L_x_regions_rt0 <- length(which(hist_rt_reference_peaks[["counts"]] == 0))
    if (L_x_regions_rt0 > 0) {
      IPA_logRecorder("WARNING!!! Reference peaks were not detected for the entire range of the retention times! Please see the 'Ref_peaks_distribution.png' in the `peak_alignment` folder!")
    }
    ##
    IPA_logRecorder("Initiated RT correction!")
    ##
    RTcorrectionMethod <- PARAM[which(PARAM[, 1] == 'PARAM0032'), 2]
    refPeakTolerance <- tryCatch(as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0033'), 2]), error = function(e) {5}, warning = function(w) {5})
    degreePolynomial <- tryCatch(as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0034'), 2]), error = function(e) {3}, warning = function(w) {3})
    ##
    peaklistFileName_samples <- setdiff(peaklistFileNames, refPeaklistFileNames)
    ############################################################################
    if (number_processing_threads == 1) {
      listCorrectedRTsamples <- lapply(peaklistFileName_samples, function(i) {
        sample_rt_corrector(referenceMZRTpeaks, inputPathPeaklist, i, massAccuracy, RTcorrectionMethod, refPeakTolerance, degreePolynomial)
      })
    } else {
      ## Processing OS
      osType <- Sys.info()[['sysname']]
      if (osType == "Windows") {
        ##
        cl <- makeCluster(number_processing_threads)
        registerDoParallel(cl)
        ##
        listCorrectedRTsamples <- foreach(i = peaklistFileName_samples, .verbose = FALSE) %dopar% {
          sample_rt_corrector(referenceMZRTpeaks, inputPathPeaklist, i, massAccuracy, RTcorrectionMethod, refPeakTolerance, degreePolynomial)
        }
        ##
        stopCluster(cl)
        ##
      } else if (osType == "Linux") {
        ##
        listCorrectedRTsamples <- mclapply(peaklistFileName_samples, function(i) {
          sample_rt_corrector(referenceMZRTpeaks, inputPathPeaklist, i, massAccuracy, RTcorrectionMethod, refPeakTolerance, degreePolynomial)
        }, mc.cores = number_processing_threads)
        ##
        closeAllConnections()
      }
    }
    names(listCorrectedRTsamples) <- peaklistFileName_samples
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
    IPA_logRecorder("Completed RT correction!")
    ##
  } else {
    listCorrectedRTpeaklists <- lapply(peaklistFileNames, function(i) {
      loadRdata(paste0(inputPathPeaklist, "/", i))[, 3]
    })
  }
  ##
  names(listCorrectedRTpeaklists) <- peaklistFileNames
  ##
  save(listCorrectedRTpeaklists, file = paste0(OutputPath_peak_alignment, "/listCorrectedRTpeaklists.Rdata"))
  ##
  if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0002'), 2]) == "yes") {
    IPA_logRecorder("Initiated peak alignment on the entire peaklists using peak IDs in each peaklist!")
    peakXcol <- peak_alignment(inputPathPeaklist, peaklistFileNames, listCorrectedRTpeaklists, massAccuracy, RTtolerance, noQuantile, number_processing_threads)
    #
    colnames(peakXcol) <- c("m/z", "RT", file_names_hrms)
    save(peakXcol, file = paste0(OutputPath_peak_alignment, "/peak_Xcol.Rdata"))
    IPA_logRecorder("Stored aligned indexed table as 'peak_Xcol.Rdata' in the `peak_alignment` folder!")
    ##
    IPA_logRecorder("Initiated generating aligned peak tables for the peak height, peak area, and R13C values!")
    ##
    listHeightAreaR13C <- peakXcolFiller(peakXcol, inputPathPeaklist)
    peak_height <- listHeightAreaR13C[["peak_height"]]
    peak_area <- listHeightAreaR13C[["peak_area"]]
    peak_R13C <- listHeightAreaR13C[["peak_R13C"]]
    ##
    IPA_logRecorder("Initiated saving aligned peak tables")
    opendir(OutputPath_peak_alignment)
    ##
    save(peak_height, file = paste0(OutputPath_peak_alignment, "/peak_height.Rdata"))
    write.csv(peak_height, file = paste0(OutputPath_peak_alignment, "/peak_height.csv"), row.names = FALSE)
    ##
    save(peak_area, file = paste0(OutputPath_peak_alignment, "/peak_area.Rdata"))
    write.csv(peak_area, file = paste0(OutputPath_peak_alignment, "/peak_area.csv"), row.names = FALSE)
    ##
    save(peak_R13C, file = paste0(OutputPath_peak_alignment, "/peak_R13C.Rdata"))
    write.csv(peak_R13C, file = paste0(OutputPath_peak_alignment, "/peak_R13C.csv"), row.names = FALSE)
    ##
    IPA_logRecorder("Aligned peak height, peak area, and R13C tables were stored in `.Rdata` and `.csv` formats in the `peak_alignment` folder!")
    ##
    ############################################################################
    ##
    IPA_logRecorder("Initiated detecting correlating peaks on the peak height table!")
    ##
    correlationListHeight <- peakPropertyTableCorrelation(peakPropertyTable = peak_height, RTtolerance = RTtolerance, minFreqDetection = 1, method = "pearson", minThresholdCorrelation = 0.50, number_processing_threads)
    ##
    save(correlationListHeight, file = paste0(OutputPath_peak_alignment, "/correlationListHeight.Rdata"))
    ##
    IPA_logRecorder("Stored the correlating peaks on the peak height table as `correlationListHeight.Rdata` in the `peak_alignement` folder!")
    IPA_logRecorder("Completed generating data for peak alignment!")
    IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
    ##
    ############################################################################
    ##
  }
}
