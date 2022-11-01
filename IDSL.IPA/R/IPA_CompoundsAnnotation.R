IPA_CompoundsAnnotation <- function(PARAM) {
  ##
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  IPA_logRecorder("Initiated compound-centric peak annotation!")
  ##
  output_path <- gsub("\\", "/", PARAM[which(PARAM[, 1] == 'PARAM0010'), 2], fixed = TRUE)
  Output_CSV <- paste0(output_path, "/compound_centric_annotation")
  if (!dir.exists(Output_CSV)) {
    dir.create(Output_CSV, recursive = TRUE)
  }
  opendir(Output_CSV)
  ##
  ref_table <- readxl::read_xlsx(PARAM[which(PARAM[, 1] == 'PARAM0042'), 2])
  mass_error <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0043'), 2])   # Mass accuracy to cluster m/z in consecutive scans
  rt_error <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0044'), 2])
  x0045 <- tolower(PARAM[which(PARAM[, 1] == 'PARAM0045'), 2])
  name_compounds <- IPA_gsub(c("/", "\\", ":", "*", "?", '"', "<", ">", "|", "."), "_", ref_table$name, fixed = TRUE)  # To replace invalid characters
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
    HRMSfileNames <- gsub("^peaklist_", "", peaklistFileNames)
    HRMSfileNames <- gsub(".Rdata$", "", HRMSfileNames)
    ##
    annotation_list <- lapply(1:L_nc, function(i) {
      MAT <- matrix(rep(0, L_PL*26), nrow = L_PL) # 24 + 2 = 26
      MAT <- data.frame(MAT)
      colnames(MAT) <- c("SampleID", "PeakID", "ScanNumberStart","ScanNumberEnd","RetentionTimeApex","PeakHeight","PeakArea",
                         "NumberDetectedScans(nIsoPair)","RCS(%)","m/z 12C","CumulatedIntensity",
                         "m/z 13C","Ratio 13C CumulatedIntensity","PeakWidthBaseline","Ratio PeakWidth @ 50%",
                         "SeperationTray","AsymmetryFactor @ 10%","USPTailingFactor @ 5%",
                         "Skewness_DerivativeMethod", "Symmetry PseudoMoments","Skewness PseudoMoments",
                         "Gaussianity", "S/N", "S/N xcms method", "S/N RMS", "Sharpness")
      return(MAT)
    })
    progressBARboundaries <- txtProgressBar(min = 1, max = L_PL, initial = 1, style = 3)
    ##
    if (x0045 == "yes") {
      listCorrectedRTpeaklists <- loadRdata(paste0(output_path, "/peak_alignment/listCorrectedRTpeaklists.Rdata"))
      for (i in 1:L_PL) {
        setTxtProgressBar(progressBARboundaries, i)
        peaklist <- loadRdata(paste0(inputPathPeaklist, "/", peaklistFileNames[i]))
        mz_i <- peaklist[, 8]
        RT_i <- listCorrectedRTpeaklists[[peaklistFileNames[i]]]
        for (j in 1:L_nc) {
          annotation_list[[j]][i, 1] <- HRMSfileNames[i]
          x_compound <- which(abs(mz_compounds[j] - mz_i) <= mass_error &
                                abs(rt_compounds[j] - RT_i) <= rt_error)
          if (length(x_compound) > 0) {
            if (length(x_compound) > 1) {
              x_min <- which.min(abs(rt_compounds[j] - RT_i[x_compound]))
              x_compound <- x_compound[x_min[1]]
            }
            annotation_list[[j]][i, 2:26] <- c(x_compound, peaklist[x_compound, ])
          }
        }
      }
    } else {
      for (i in 1:L_PL) {
        setTxtProgressBar(progressBARboundaries, i)
        peaklist <- loadRdata(paste0(inputPathPeaklist, "/", peaklistFileNames[i]))
        mz_i <- peaklist[, 8]
        RT_i <- peaklist[, 3]
        for (j in 1:L_nc) {
          annotation_list[[j]][i, 1] <- HRMSfileNames[i]
          x_compound <- which(abs(mz_compounds[j] - mz_i) <= mass_error &
                                abs(rt_compounds[j] - RT_i) <= rt_error)
          if (length(x_compound) > 0) {
            if (length(x_compound) > 1) {
              x_min <- which.min(abs(rt_compounds[j] - RT_i[x_compound]))
              x_compound <- x_compound[x_min[1]]
            }
            annotation_list[[j]][i, 2:26] <- c(x_compound, peaklist[x_compound, ])
          }
        }
      }
    }
    close(progressBARboundaries)
    ##
    ##############################################################################
    ##
    IPA_logRecorder("Compound centric annotation data are stored in the `compound_centric_annotation` folder!")
    ##
    output_path_eics <- paste0(output_path, "/IPA_EIC")
    if (dir.exists(output_path_eics)) {
      exportEICcheck <- TRUE
      ##
      eicPNGfolderlist <- lapply(1:L_PL, function(j) {
        eic_HRMS_folder <- paste0(output_path_eics, "/", HRMSfileNames[j])
        if (dir.exists(eic_HRMS_folder)) {
          dirPNG <- dir(path = eic_HRMS_folder, pattern = ".png")
          xPNG <- as.numeric(do.call(c, lapply(strsplit(dirPNG, "_"), function(k) {k[length(k) - 5]})))
          list(eic_HRMS_folder, dirPNG, xPNG)
        }
      })
      ##
      IPA_logRecorder("If EICs exist, EIC folders are generated for individual compound in the `compound_centric_annotation` folder!")
      ##
    } else {
      exportEICcheck <- FALSE
    }
    ##
    progressBARboundaries <- txtProgressBar(min = 1, max = L_nc, initial = 1, style = 3)
    ##
    for (i in 1:L_nc) {
      A <- annotation_list[[i]]
      nameAnnotatedCompound <- paste0(i, "_annotated_compound__", name_compounds[i])
      write.csv(A, file = paste0(Output_CSV, "/", nameAnnotatedCompound, ".csv"), row.names = FALSE)
      ##
      if (exportEICcheck) {
        sampleName <- A[, 1]
        peakIDs <- as.numeric(A[, 2])
        xNon0 <- which(peakIDs > 0)
        ##
        if (length(xNon0) > 0) {
          eic_nameAnnotatedCompound <- paste0(Output_CSV, "/", nameAnnotatedCompound)
          if (!dir.exists(eic_nameAnnotatedCompound)) {
            dir.create(eic_nameAnnotatedCompound, recursive = TRUE)
          }
          for (j in xNon0) {
            eicPNG <- eicPNGfolderlist[[j]]
            if (!is.null(eicPNGfolderlist)) {
              eic_HRMS_folder <-  eicPNG[[1]]
              dirPNG <- eicPNG[[2]]
              xPNG <- eicPNG[[3]]
              xCompoundMatch <- which(xPNG == peakIDs[j])
              sourceFileAddress <- paste0(eic_HRMS_folder, "/", dirPNG[xCompoundMatch])
              destFileAddress <- paste0(eic_nameAnnotatedCompound, "/", dirPNG[xCompoundMatch])
              file.copy(from = sourceFileAddress, to = destFileAddress)
            }
          }
        }
      }
      setTxtProgressBar(progressBARboundaries, i)
    }
    ##
    close(progressBARboundaries)
    ##
    ##############################################################################
    ##
    IPA_logRecorder("Completed compound-centric peak annotation!")
    IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  }
}
