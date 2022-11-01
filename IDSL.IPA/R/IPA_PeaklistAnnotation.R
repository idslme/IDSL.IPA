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
  name_compounds <- ref_table$name
  L_nc <- length(name_compounds)
  if (L_nc > 0) {
    mz_compounds <- ref_table$`m/z`
    rt_compounds <- ref_table$RT
    #
    inputPathPeaklist <- paste0(output_path, "/peaklists")
    peaklistFileNames1 <- dir(path = inputPathPeaklist, pattern = ".Rdata")
    peaklistFileNames2 <- dir(path = inputPathPeaklist, pattern = "peaklist_")
    peaklistFileNames <- peaklistFileNames1[peaklistFileNames1 %in% peaklistFileNames2]
    L_PL <- length(peaklistFileNames)
    ##
    if (gapFillingCheck) {
      input_path_hrms <- PARAM[which(PARAM[, 1] == 'PARAM0007'), 2]
      if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]) == "all") {
        file_names_hrms <- dir(path = input_path_hrms)
        file_names_hrms <- file_names_hrms[grep(paste0(".", tolower(PARAM[which(PARAM[, 1] == 'PARAM0009'), 2]), "$"), file_names_hrms, ignore.case = TRUE)]
      } else {
        samples_string <- PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]
        file_names_hrms <- strsplit(samples_string, ";")[[1]] # files used as reference m/z-RT
      }
      peaklistFileNames_hrms1 <- gsub(".Rdata", "", peaklistFileNames)
      peaklistFileNames_hrms2 <- gsub("peaklist_", "", peaklistFileNames_hrms1)
      peaklistFileNames_hrms <- which(file_names_hrms %in% peaklistFileNames_hrms2)
      if (length(peaklistFileNames_hrms) != L_PL) {
        stop(IPA_logRecorder("Error!!! peaklist files are not available for all selected HRMS files!"))
      }
    }
    ##
    annotated_peak_indices <- matrix(rep(0, L_nc*L_PL), nrow = L_nc)
    progressBARboundaries <- txtProgressBar(min = 1, max = L_PL, initial = 1, style = 3)
    if (retentionTimeCorrectionCheck) {
      listCorrectedRTpeaklists <- loadRdata(paste0(output_path, "/peak_alignment/listCorrectedRTpeaklists.Rdata"))
      ##
      for (i in 1:L_PL) {
        setTxtProgressBar(progressBARboundaries, i)
        peaklist <- loadRdata(paste0(inputPathPeaklist, "/", peaklistFileNames[i]))
        mz_i <- peaklist[, 8]
        RT_i <- listCorrectedRTpeaklists[[peaklistFileNames[i]]]
        for (j in 1:L_nc) {
          x_compound <- which(abs(mz_compounds[j] - mz_i) <= massAccuracy &
                                abs(rt_compounds[j] - RT_i) <= RTtolerance)
          if (length(x_compound) > 0) {
            if (length(x_compound) > 1) {
              x_min <- which.min(abs(rt_compounds[j] - RT_i[x_compound]))
              x_compound <- x_compound[x_min[1]]
            }
            annotated_peak_indices[j, i] <- x_compound
          }
        }
      }
    } else {
      listCorrectedRTpeaklists <- NULL
      for (i in 1:L_PL) {
        setTxtProgressBar(progressBARboundaries, i)
        peaklist <- loadRdata(paste0(inputPathPeaklist, "/", peaklistFileNames[i]))
        mz_i <- peaklist[, 8]
        RT_i <- peaklist[, 3]
        for (j in 1:L_nc) {
          x_compound <- which(abs(mz_compounds[j] - mz_i) <= massAccuracy &
                                abs(rt_compounds[j] - RT_i) <= RTtolerance)
          if (length(x_compound) > 0) {
            if (length(x_compound) > 1) {
              x_min <- which.min(abs(rt_compounds[j] - RT_i[x_compound]))
              x_compound <- x_compound[x_min[1]]
            }
            annotated_peak_indices[j, i] <- x_compound
          }
        }
      }
    }
    close(progressBARboundaries)
    annotated_peak_indices <- cbind(mz_compounds, rt_compounds, annotated_peak_indices)
    colnames(annotated_peak_indices) <- c("m/z", "RT", file_names_hrms)
    ##
    if (!dir.exists(Output_Xcol)) {
      dir.create(Output_Xcol, recursive = TRUE)
    }
    opendir(Output_Xcol)
    ##
    save(annotated_peak_indices, file = paste0(Output_Xcol, "/annotated_peak_indices.Rdata"))
    IPA_logRecorder("Aligned indexed table from individual peaklists were stored as 'annotated_peak_indices.Rdata' in the `sample_centric_annotation` folder!")
    listHeightAreaR13C <- peakXcolFiller(annotated_peak_indices, inputPathPeaklist)
    annotated_peak_height <- cbind(name_compounds, listHeightAreaR13C[["peak_height"]])
    colnames(annotated_peak_height) <- c("name", "m/z", "RT", file_names_hrms)
    annotated_peak_area <- cbind(name_compounds, listHeightAreaR13C[["peak_area"]])
    colnames(annotated_peak_area) <- c("name", "m/z", "RT", file_names_hrms)
    annotated_peak_R13C <- cbind(name_compounds, listHeightAreaR13C[["peak_R13C"]])
    colnames(annotated_peak_R13C) <- c("name", "m/z", "RT", file_names_hrms)
    write.csv(annotated_peak_height, file = paste0(Output_Xcol, "/annotated_peak_height.csv"), row.names = FALSE)
    write.csv(annotated_peak_area, file = paste0(Output_Xcol, "/annotated_peak_area.csv"), row.names = FALSE)
    write.csv(annotated_peak_R13C, file = paste0(Output_Xcol, "/annotated_peak_R13C.csv"), row.names = FALSE)
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
      annotated_peak_height_gapfilled <- annotated_peak_height
      annotated_peak_area_gapfilled <- annotated_peak_area
      annotated_peak_R13C_gapfilled <- annotated_peak_R13C
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
      IPA_logRecorder("Initiated saving gap-filled sample-centric peak annotation!")
      write.csv(annotated_peak_height_gapfilled, file = paste0(Output_Xcol, "/annotated_peak_height_gapfilled.csv"), row.names = FALSE)
      write.csv(annotated_peak_area_gapfilled, file = paste0(Output_Xcol, "/annotated_peak_area_gapfilled.csv"), row.names = FALSE)
      write.csv(annotated_peak_R13C_gapfilled, file = paste0(Output_Xcol, "/annotated_peak_R13C_gapfilled.csv"), row.names = FALSE)
      ##
      IPA_logRecorder("Gap-filled annotated peak height, peak area, and R13C tables were stored in `.csv` formats in the `sample_centric_annotation` folder!")
      IPA_logRecorder("Completed gap-filled sample-centric peak annotations!")
    }
    IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  }
}
