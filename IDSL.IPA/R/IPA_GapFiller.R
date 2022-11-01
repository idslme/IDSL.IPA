IPA_GapFiller <- function(PARAM) {
  ##
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  IPA_logRecorder("Initiated gap-filling!")
  ##
  number_processing_threads <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0006'), 2])
  ##
  input_path_hrms <- PARAM[which(PARAM[, 1] == 'PARAM0007'), 2]
  if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]) == "all") {
    file_name_hrms <- dir(path = input_path_hrms)
    file_name_hrms <- file_name_hrms[grep(paste0(".", tolower(PARAM[which(PARAM[, 1] == 'PARAM0009'), 2]), "$"), file_name_hrms, ignore.case = TRUE)]
  } else {
    samples_string <- PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]
    file_name_hrms <- strsplit(samples_string, ";")[[1]]
  }
  L_HRMS <- length(file_name_hrms)
  ##
  output_path <- PARAM[which(PARAM[, 1] == 'PARAM0010'), 2]
  inputPathPeaklist <- paste0(output_path, "/peaklists")
  file_names_peaklist1 <- dir(path = inputPathPeaklist, pattern = ".Rdata")
  file_names_peaklist2 <- dir(path = inputPathPeaklist, pattern = "peaklist_")
  file_names_peaklist <- file_names_peaklist1[file_names_peaklist1 %in% file_names_peaklist2]
  L_PL <- length(file_names_peaklist)
  ##
  file_names_peaklist_hrms1 <- gsub(".Rdata", "", file_names_peaklist)
  file_names_peaklist_hrms2 <- gsub("peaklist_", "", file_names_peaklist_hrms1)
  file_names_peaklist_hrms <- which(file_name_hrms %in% file_names_peaklist_hrms2)
  if (length(file_names_peaklist_hrms) != L_PL) {
    stop(IPA_logRecorder("Error!!! peaklist files are not available for all selected HRMS files!"))
  }
  ##
  peakXcol <- loadRdata(paste0(output_path, "/peak_alignment/peak_Xcol.Rdata"))
  colnamesPeakXcol <- colnames(peakXcol)[3:(ncol(peakXcol))]
  xMatchedHRMS <- which(file_name_hrms %in% colnamesPeakXcol)
  if (length(xMatchedHRMS) != L_HRMS) {
    stop(IPA_logRecorder("Error!!! `peak_Xcol.Rdata` file do not contain all selected HRMS files!"))
  }
  ##
  listCorrectedRTpeaklists <- loadRdata(paste0(output_path, "/peak_alignment/listCorrectedRTpeaklists.Rdata"))
  ##
  ionMassDifference <- tryCatch(as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0012'), 2]), error = function(e) {1.003354835336}, warning = function(w) {1.003354835336})     # Mass difference for isotopic pairs
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
  peak_area_gapfilled <- loadRdata(paste0(OutputPath_peak_alignment, "/peak_area.Rdata"))
  peak_R13C_gapfilled <- loadRdata(paste0(OutputPath_peak_alignment, "/peak_R13C.Rdata"))
  ##
  progressBARboundaries <- txtProgressBar(min = 0, max = L_HRMS, initial = 0, style = 3)
  for (i in 1:L_HRMS) {
    setTxtProgressBar(progressBARboundaries, i)
    iSample <- chromatography_undetected_list[[i]]
    if (!is.null(iSample)) {
      x_j <- iSample[, 1]
      jCounter <- 0
      for (j in x_j) {
        jCounter <- jCounter + 1
        peak_height_gapfilled[j, (i + 2)] <- iSample[jCounter, 2]
        peak_area_gapfilled[j, (i + 2)] <- iSample[jCounter, 3]
        peak_R13C_gapfilled[j, (i + 2)] <- iSample[jCounter, 4]
      }
    }
  }
  close(progressBARboundaries)
  opendir(OutputPath_peak_alignment)
  IPA_logRecorder("Initiated saving aligned gap-filled peak tables!")
  save(peak_height_gapfilled, file = paste0(OutputPath_peak_alignment, "peak_height_gapfilled.Rdata"))
  write.csv(peak_height_gapfilled, file = paste0(OutputPath_peak_alignment, "peak_height_gapfilled.csv"), row.names = FALSE)
  save(peak_area_gapfilled, file = paste0(OutputPath_peak_alignment, "peak_area_gapfilled.Rdata"))
  write.csv(peak_area_gapfilled, file = paste0(OutputPath_peak_alignment, "peak_area_gapfilled.csv"), row.names = FALSE)
  save(peak_R13C_gapfilled, file = paste0(OutputPath_peak_alignment, "peak_R13C_gapfilled.Rdata"))
  write.csv(peak_R13C_gapfilled, file = paste0(OutputPath_peak_alignment, "peak_R13C_gapfilled.csv"), row.names = FALSE)
  IPA_logRecorder("Gap-filled aligned peak height, peak area, and R13C tables were stored in `.Rdata` and `.csv` formats in the `peak_alignment` folder!")
  ##
  ##############################################################################
  ##
  IPA_logRecorder("Initiated detecting correlating peaks on the gap-filled peak height table!")
  ##
  correlationListHeight_gapfilled <- peakPropertyTableCorrelation(peakPropertyTable = peak_height_gapfilled, RTtolerance = RTtolerance, minFreqDetection = 1, method = "pearson", minThresholdCorrelation = 0.50, number_processing_threads)
  ##
  save(correlationListHeight_gapfilled, file = paste0(OutputPath_peak_alignment, "/correlationListHeight_gapfilled.Rdata"))
  ##
  IPA_logRecorder("Stored the correlating peaks on the gap-filled peak height table as `correlationListHeight_gapfilled.Rdata` in the `peak_alignement` folder!")
  IPA_logRecorder("Completed gap-filling!")
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  ##
  ##############################################################################
  ##
}
