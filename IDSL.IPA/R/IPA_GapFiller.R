IPA_GapFiller <- function(PARAM) {
  ##
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  IPA_logRecorder("Initiated gap-filling!")
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
  peakXcol <- loadRdata(paste0(output_path, "/peak_alignment/peakXcol.Rdata"))
  colnamesPeakXcol <- colnames(peakXcol)[4:(ncol(peakXcol))]
  xMatchedHRMS <- which(colnamesPeakXcol %in% file_name_hrms)
  if (length(xMatchedHRMS) > LHRMS) {
    IPA_logRecorder("Error!!! `peakXcol.Rdata` file do not contain for the following HRMS file(s):")
    for (i in xMatchedHRMS) {
      IPA_logRecorder(colnamesPeakXcol[i])
    }
    stop()
  }
  ##
  listCorrectedRTpeaklists <- loadRdata(paste0(output_path, "/peak_alignment/listCorrectedRTpeaklists.Rdata"))
  ##
  ionMassDifference <- tryCatch(as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0012'), 2]), error = function(e) {1.003354835336}, warning = function(w) {1.003354835336})     # Mass difference for isotopic pairs
  #
  if (ionMassDifference <= 1.00336 & ionMassDifference >= 1.00335) {
    IPA_logRecorder("Carbon isotopes are selected for ion pairing in the gap-filling step!")
  } else {
    IPA_logRecorder(paste0("Mass difference to pair ions is '", ionMassDifference, " Da' in the gap-filling step!"))
  }
  ##
  massAccuracy <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0038'), 2])   # Mass accuracy to cluster m/z in consecutive scans
  RTtolerance <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0039'), 2])
  scanTolerance <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0040'), 2])
  ##
  chromatography_undetected_list <- gapFillingCore(input_path_hrms, peakXcol, massAccuracy, RTtolerance, scanTolerance,
                                                   retentionTimeCorrectionCheck = TRUE, listCorrectedRTpeaklists,
                                                   inputPathPeaklist, ionMassDifference, number_processing_threads)
  peakXcol <- NULL
  listCorrectedRTpeaklists <- NULL
  ##
  IPA_logRecorder("Initiated filling gaps of the aligned peak height, peak area, and R13C tables!")
  OutputPath_peak_alignment <- paste0(output_path, "/peak_alignment/")
  peak_height_gapfilled <- loadRdata(paste0(OutputPath_peak_alignment, "/peak_height.Rdata"))
  peak_height_gapfilled <- peak_height_gapfilled[, -c(4, 5)]
  peak_area_gapfilled <- loadRdata(paste0(OutputPath_peak_alignment, "/peak_area.Rdata"))
  peak_area_gapfilled <- peak_area_gapfilled[, -c(4, 5)]
  peak_R13C_gapfilled <- loadRdata(paste0(OutputPath_peak_alignment, "/peak_R13C.Rdata"))
  peak_R13C_gapfilled <- peak_R13C_gapfilled[, -c(4, 5)]
  ##
  progressBARboundaries <- txtProgressBar(min = 0, max = LHRMS, initial = 0, style = 3)
  for (i in 1:LHRMS) {
    setTxtProgressBar(progressBARboundaries, i)
    iSample <- chromatography_undetected_list[[i]]
    if (!is.null(iSample)) {
      x_j <- iSample[, 1]
      jCounter <- 0
      for (j in x_j) {
        jCounter <- jCounter + 1
        peak_height_gapfilled[j, (i + 3)] <- iSample[jCounter, 2]
        peak_area_gapfilled[j, (i + 3)] <- iSample[jCounter, 3]
        peak_R13C_gapfilled[j, (i + 3)] <- iSample[jCounter, 4]
      }
    }
  }
  close(progressBARboundaries)
  chromatography_undetected_list <- NULL
  ##
  maxRedundantPeakFlagging <-  tryCatch(as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0037'), 2])/100, error = function(e) {0.25}, warning = function(w) {0.25})
  if (maxRedundantPeakFlagging < 0 | maxRedundantPeakFlagging > 1) {
    maxRedundantPeakFlagging <- 0.25
  }
  ##
  IPA_logRecorder("Initiated flagging suspicious gap-filled aligned peaks, calculating median peak properties and saving aligned gap-filled peak tables!")
  ##
  ##############################################################################
  ##############################################################################
  ##
  peak_height_gapfilled[, 3] <- peakPropertyTableFreqCalculator(peak_height_gapfilled, startColumnIndex = 4, number_processing_threads)
  ##
  falggingVector <- peakXcolFlagger(peak_height_gapfilled[, 1], peak_height_gapfilled[, 2], peak_height_gapfilled[, 3],
                                    massAccuracy, 3*RTtolerance, maxRedundantPeakFlagging)
  ##
  peak_height_gapfilled <- peakPropertyTableMedianCalculator(peak_height_gapfilled, falggingVector, number_processing_threads)
  colnames(peak_height_gapfilled)[c(3, 4, 5)] <- c("freqGapFilledPeakHeight", "medianGapFilledPeakHeight", "Flag")
  ##
  save(peak_height_gapfilled, file = paste0(OutputPath_peak_alignment, "peak_height_gapfilled.Rdata"))
  write.csv(peak_height_gapfilled, file = paste0(OutputPath_peak_alignment, "peak_height_gapfilled.csv"), row.names = TRUE)
  ##
  ##############################################################################
  ##
  peak_area_gapfilled[, 3] <- peakPropertyTableFreqCalculator(peak_area_gapfilled, startColumnIndex = 4, number_processing_threads)
  ##
  falggingVector <- peakXcolFlagger(peak_area_gapfilled[, 1], peak_area_gapfilled[, 2], peak_area_gapfilled[, 3],
                                    massAccuracy, 3*RTtolerance, maxRedundantPeakFlagging)
  ##
  peak_area_gapfilled <- peakPropertyTableMedianCalculator(peak_area_gapfilled, falggingVector, number_processing_threads)
  colnames(peak_area_gapfilled)[c(3, 4, 5)] <- c("freqGapFilledPeakArea", "medianGapFilledPeakArea", "Flag")
  ##
  save(peak_area_gapfilled, file = paste0(OutputPath_peak_alignment, "peak_area_gapfilled.Rdata"))
  write.csv(peak_area_gapfilled, file = paste0(OutputPath_peak_alignment, "peak_area_gapfilled.csv"), row.names = TRUE)
  peak_area_gapfilled <- NULL
  ##
  ##############################################################################
  ##
  peak_R13C_gapfilled[, 3] <- peakPropertyTableFreqCalculator(peak_R13C_gapfilled, startColumnIndex = 4, number_processing_threads)
  ##
  falggingVector <- peakXcolFlagger(peak_R13C_gapfilled[, 1], peak_R13C_gapfilled[, 2], peak_R13C_gapfilled[, 3],
                                    massAccuracy, 3*RTtolerance, maxRedundantPeakFlagging)
  ##
  peak_R13C_gapfilled <- peakPropertyTableMedianCalculator(peak_R13C_gapfilled, falggingVector, number_processing_threads)
  colnames(peak_R13C_gapfilled)[c(3, 4, 5)] <- c("freqGapFilledR13C", "medianGapFilledR13C", "Flag")
  ##
  save(peak_R13C_gapfilled, file = paste0(OutputPath_peak_alignment, "peak_R13C_gapfilled.Rdata"))
  write.csv(peak_R13C_gapfilled, file = paste0(OutputPath_peak_alignment, "peak_R13C_gapfilled.csv"), row.names = TRUE)
  peak_R13C_gapfilled <- NULL
  ##
  ##############################################################################
  ##############################################################################
  ##
  IPA_logRecorder("Completed flagging suspicious gap-filled aligned peaks, calculating median peak properties and saving aligned gap-filled peak tables!")
  ##
  ##############################################################################
  ##############################################################################
  ##
  if (LHRMS > 2) {
    correlationMethod <- tryCatch(tolower(PARAM[which(PARAM[, 1] == 'PARAM_ALG1'), 2]), error = function(e) {"pearson"}, warning = function(w) {"pearson"})
    if (!(correlationMethod == "pearson" | correlationMethod == "spearman")) {
      correlationMethod <- "pearson"
    }
    minThresholdCorrelation <- tryCatch(as.numeric(PARAM[which(PARAM[, 1] == 'PARAM_ALG2'), 2]), error = function(e) {0.75}, warning = function(w) {0.75})
    if (minThresholdCorrelation < 0.5 | minThresholdCorrelation > 1) {
      minThresholdCorrelation <- 0.75
    }
    minFreqDetection <- tryCatch(floor(as.numeric(PARAM[which(PARAM[, 1] == 'PARAM_ALG3'), 2])), error = function(e) {0.75}, warning = function(w) {0.75})
    if (minFreqDetection < 3) {
      minFreqDetection <- 3
    }
    minRatioDetection <- tryCatch(as.numeric(PARAM[which(PARAM[, 1] == 'PARAM_ALG4'), 2])/100, error = function(e) {0.10}, warning = function(w) {0.10})
    if (minRatioDetection < 0 | minRatioDetection > 1) {
      minRatioDetection <- 0.1
    }
    ##
    IPA_logRecorder(paste0("Initiated detecting correlating peaks on the gap-filled aligned peak height table using the `", correlationMethod, "` method with coefficients `>=",
                           minThresholdCorrelation,"` and minimum number of complete observations `>=", minFreqDetection, "` combined with observation percentage `>=", minRatioDetection*100, "%`!"))
    ##
    alignedGapFilledPeakHeightTableCorrelationList <- alignedPeakPropertyTableCorrelationListCalculator(peakPropertyTable = peak_height_gapfilled, RTtolerance, minFreqDetection, minRatioDetection,
                                                                                                        method = correlationMethod, minThresholdCorrelation, number_processing_threads)
    peak_height_gapfilled <- NULL
    ##
    save(alignedGapFilledPeakHeightTableCorrelationList, file = paste0(OutputPath_peak_alignment, "/alignedGapFilledPeakHeightTableCorrelationList.Rdata"))
    ##
    IPA_logRecorder("Stored the correlating peaks on the gap-filled peak height table as `alignedGapFilledPeakHeightTableCorrelationList.Rdata` in the `peak_alignement` folder (Only correlating aligned peak IDs are presented)!")
  }
  ##
  ##############################################################################
  ##############################################################################
  ##
  IPA_logRecorder("Completed gap-filling!")
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  ##
  return()
}
