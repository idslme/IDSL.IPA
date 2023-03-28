IPA_PeaklistAnnotation <- function(PARAM) {
  ##
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  IPA_logRecorder("Initiated sample-centric peak annotation!")
  ##
  number_processing_threads <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0006'), 2])
  output_path <- PARAM[which(PARAM[, 1] == 'PARAM0010'), 2]
  Output_Xcol <- paste0(output_path, "/sample_centric_annotation")
  ##
  IPA_logRecorder("Sample centric annotation data are stored in the `sample_centric_annotation` folder!")
  ##
  ref_table <- readxl::read_xlsx(PARAM[which(PARAM[, 1] == 'PARAM0042'), 2])
  massAccuracy <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0043'), 2])   # Mass accuracy to cluster m/z in consecutive scans
  RTtolerance <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0044'), 2])
  retentionTimeCorrectionCheck <- if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0045'), 2]) == "yes") {TRUE} else {FALSE}
  gapFillingCheck <- if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0048'), 2]) == "yes") {TRUE} else {FALSE}
  compoundNames <- ref_table$name
  nCompoundNames <- length(compoundNames)
  if (nCompoundNames > 0) {
    mz_compounds <- ref_table$`m/z`
    rt_compounds <- ref_table$RT
    #
    inputPathPeaklist <- paste0(output_path, "/peaklists")
    peaklistFileNames <- dir(path = inputPathPeaklist, pattern = ".Rdata$")
    peaklistFileNames <- peaklistFileNames[grep("^peaklist_", peaklistFileNames)]
    L_PL <- length(peaklistFileNames)
    ##
    if (gapFillingCheck) {
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
    } else {
      file_name_hrms <- gsub("^peaklist_|.Rdata$", "", peaklistFileNames)
    }
    ##
    annotated_peak_indices <- matrix(rep(0, nCompoundNames*L_PL), nrow = nCompoundNames)
    progressBARboundaries <- txtProgressBar(min = 1, max = L_PL, initial = 1, style = 3)
    if (retentionTimeCorrectionCheck) {
      listCorrectedRTpeaklists <- loadRdata(paste0(output_path, "/peak_alignment/listCorrectedRTpeaklists.Rdata"))
      ##
      for (i in 1:L_PL) {
        peaklist <- loadRdata(paste0(inputPathPeaklist, "/", peaklistFileNames[i]))
        RT_i <- listCorrectedRTpeaklists[[peaklistFileNames[i]]]
        for (j in 1:nCompoundNames) {
          x_compound <- mzRTindexer(MZvec = peaklist[, 8], RTvec = RT_i,
                                    MZref = mz_compounds[j], RTref = rt_compounds[j], massAccuracy, RTtolerance)
          if (!is.null(x_compound)) {
            annotated_peak_indices[j, i] <- x_compound
          }
        }
        setTxtProgressBar(progressBARboundaries, i)
      }
    } else {
      ##
      listCorrectedRTpeaklists <- NULL
      ##
      for (i in 1:L_PL) {
        peaklist <- loadRdata(paste0(inputPathPeaklist, "/", peaklistFileNames[i]))
        for (j in 1:nCompoundNames) {
          x_compound <- mzRTindexer(MZvec = peaklist[, 8], RTvec = peaklist[, 3],
                                    MZref = mz_compounds[j], RTref = rt_compounds[j], massAccuracy, RTtolerance)
          if (!is.null(x_compound)) {
            annotated_peak_indices[j, i] <- x_compound
          }
        }
        setTxtProgressBar(progressBARboundaries, i)
      }
    }
    close(progressBARboundaries)
    ##
    frequencyDetections <- peakPropertyTableFreqCalculator(annotated_peak_indices, startColumnIndex = 1, number_processing_threads, allowedVerbose = FALSE)
    annotated_peak_indices <- cbind(mz_compounds, rt_compounds, frequencyDetections, annotated_peak_indices)
    colnames(annotated_peak_indices) <- c("mz", "RT", "frequencyPeakXcol", file_name_hrms)
    ##
    if (!dir.exists(Output_Xcol)) {
      dir.create(Output_Xcol, recursive = TRUE)
    }
    opendir(Output_Xcol)
    ##
    save(annotated_peak_indices, file = paste0(Output_Xcol, "/annotated_peak_indices.Rdata"))
    IPA_logRecorder("Aligned indexed table from individual peaklists were stored as `annotated_peak_indices.Rdata` in the `sample_centric_annotation` folder!")
    ##
    listHeightAreaR13C <- peakXcolFiller(annotated_peak_indices, inputPathPeaklist)
    annotated_peak_height <- peakPropertyTableMedianCalculator(listHeightAreaR13C[["peak_height"]], falggingVector = NULL, number_processing_threads)
    annotated_peak_height <- cbind(compoundNames, annotated_peak_height)
    colnames(annotated_peak_height) <- c("name", "m/z", "RT", "freqPeakHeight", "medianPeakHeight", file_name_hrms)
    ##
    annotated_peak_area <- peakPropertyTableMedianCalculator(listHeightAreaR13C[["peak_area"]], falggingVector = NULL, number_processing_threads)
    annotated_peak_area <- cbind(compoundNames, annotated_peak_area)
    colnames(annotated_peak_area) <- c("name", "m/z", "RT", "freqPeakArea", "medianPeakArea", file_name_hrms)
    ##
    annotated_peak_R13C <- peakPropertyTableMedianCalculator(listHeightAreaR13C[["peak_R13C"]], falggingVector = NULL, number_processing_threads)
    annotated_peak_R13C <- cbind(compoundNames, annotated_peak_R13C)
    colnames(annotated_peak_R13C) <- c("name", "m/z", "RT", "freqR13C", "medianR13C", file_name_hrms)
    ##
    write.csv(annotated_peak_height, file = paste0(Output_Xcol, "/annotated_peak_height.csv"), row.names = TRUE)
    write.csv(annotated_peak_area, file = paste0(Output_Xcol, "/annotated_peak_area.csv"), row.names = TRUE)
    write.csv(annotated_peak_R13C, file = paste0(Output_Xcol, "/annotated_peak_R13C.csv"), row.names = TRUE)
    ##
    IPA_logRecorder("Annotated peak height, peak area, and R13C tables were stored in `.Rdata` and `.csv` formats in the `sample_centric_annotation` folder!")
    IPA_logRecorder("Completed sample-centric peak annotations for peak height, peak area, and R13C tables!")
    ##
    if (gapFillingCheck) {
      IPA_logRecorder("Initiated gap-filling for sample-centric peak annotation!")
      ##
      ionMassDifference <- tryCatch(as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0012'), 2]), error = function(e) {1.003354835336}, warning = function(w) {1.003354835336})     # Mass difference for isotopic pairs
      massAccuracy <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0038'), 2])   # Mass accuracy to cluster m/z in consecutive scans
      RTtolerance <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0039'), 2])
      scanTolerance <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0040'), 2])
      ##
      chromatography_undetected_list <- gapFillingCore(input_path_hrms, peakXcol = annotated_peak_indices, massAccuracy, RTtolerance,
                                                       scanTolerance, retentionTimeCorrectionCheck, listCorrectedRTpeaklists,
                                                       inputPathPeaklist, ionMassDifference, number_processing_threads)
      ##
      annotated_peak_height_gapfilled <- matrix(as.numeric(annotated_peak_height[, c(-1, -5)]), nrow = nCompoundNames)
      annotated_peak_height <- NULL
      annotated_peak_area_gapfilled <- matrix(as.numeric(annotated_peak_area[, c(-1, -5)]), nrow = nCompoundNames)
      annotated_peak_area <- NULL
      annotated_peak_R13C_gapfilled <- matrix(as.numeric(annotated_peak_R13C[, c(-1, -5)]), nrow = nCompoundNames)
      annotated_peak_R13C <- NULL
      ##
      progressBARboundaries <- txtProgressBar(min = 1, max = L_PL, initial = 1, style = 3)
      for (i in 1:L_PL) {
        setTxtProgressBar(progressBARboundaries, i)
        iSample <- chromatography_undetected_list[[i]]
        if (!is.null(iSample)) {
          x_j <- iSample[, 1]
          jCounter <- 0
          for (j in x_j) {
            jCounter <- jCounter + 1
            annotated_peak_height_gapfilled[j, (i + 3)] <- iSample[jCounter, 2]
            annotated_peak_area_gapfilled[j, (i + 3)] <- iSample[jCounter, 3]
            annotated_peak_R13C_gapfilled[j, (i + 3)] <- iSample[jCounter, 4]
          }
        }
      }
      close(progressBARboundaries)
      ##
      IPA_logRecorder("Initiated saving gap-filled sample-centric peak annotation!")
      ##
      ##########################################################################
      ##
      annotated_peak_height_gapfilled[, 3] <- peakPropertyTableFreqCalculator(annotated_peak_height_gapfilled, startColumnIndex = 4, number_processing_threads)
      annotated_peak_height_gapfilled <- peakPropertyTableMedianCalculator(annotated_peak_height_gapfilled, falggingVector = NULL, number_processing_threads)
      annotated_peak_height_gapfilled <- cbind(compoundNames, annotated_peak_height_gapfilled)
      colnames(annotated_peak_height_gapfilled) <- c("name", "m/z", "RT", "freqGapFilledPeakHeight", "medianGapfilledPeakHeight", file_name_hrms)
      write.csv(annotated_peak_height_gapfilled, file = paste0(Output_Xcol, "/annotated_peak_height_gapfilled.csv"), row.names = TRUE)
      annotated_peak_height_gapfilled <- NULL
      ##
      ##########################################################################
      ##
      annotated_peak_area_gapfilled[, 3] <- peakPropertyTableFreqCalculator(annotated_peak_area_gapfilled, startColumnIndex = 4, number_processing_threads)
      annotated_peak_area_gapfilled <- peakPropertyTableMedianCalculator(annotated_peak_area_gapfilled, falggingVector = NULL, number_processing_threads)
      annotated_peak_area_gapfilled <- cbind(compoundNames, annotated_peak_area_gapfilled)
      colnames(annotated_peak_area_gapfilled) <- c("name", "m/z", "RT", "freqGapFilledPeakArea", "medianGapfilledPeakArea", file_name_hrms)
      write.csv(annotated_peak_area_gapfilled, file = paste0(Output_Xcol, "/annotated_peak_area_gapfilled.csv"), row.names = TRUE)
      annotated_peak_area_gapfilled <- NULL
      ##
      ##########################################################################
      ##
      annotated_peak_R13C_gapfilled[, 3] <- peakPropertyTableFreqCalculator(annotated_peak_R13C_gapfilled, startColumnIndex = 4, number_processing_threads)
      annotated_peak_R13C_gapfilled <- peakPropertyTableMedianCalculator(annotated_peak_R13C_gapfilled, falggingVector = NULL, number_processing_threads)
      annotated_peak_R13C_gapfilled <- cbind(compoundNames, annotated_peak_R13C_gapfilled)
      colnames(annotated_peak_R13C_gapfilled) <- c("name", "m/z", "RT", "freqGapFilledR13C", "medianGapfilledR13C", file_name_hrms)
      write.csv(annotated_peak_R13C_gapfilled, file = paste0(Output_Xcol, "/annotated_peak_R13C_gapfilled.csv"), row.names = TRUE)
      annotated_peak_R13C_gapfilled <- NULL
      ##
      ##########################################################################
      ##
      IPA_logRecorder("Gap-filled annotated peak height, peak area, and R13C tables were stored in `.csv` formats in the `sample_centric_annotation` folder!")
      IPA_logRecorder("Completed gap-filled sample-centric peak annotations!")
    }
    IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  }
  ##
  ##############################################################################
  ##
  return()
}
