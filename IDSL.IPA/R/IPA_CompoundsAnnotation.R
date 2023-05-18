IPA_CompoundsAnnotation <- function(PARAM) {
  ##
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  IPA_logRecorder("Initiated compound-centric peak annotation!")
  ##
  output_path <- gsub("\\", "/", PARAM[which(PARAM[, 1] == 'PARAM0010'), 2], fixed = TRUE)
  OutputCSV <- paste0(output_path, "/compound_centric_annotation")
  if (!dir.exists(OutputCSV)) {
    dir.create(OutputCSV, recursive = TRUE)
  }
  ##
  ref_table <- readxl::read_xlsx(PARAM[which(PARAM[, 1] == 'PARAM0042'), 2])
  massAccuracy <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0043'), 2])   # Mass accuracy to cluster m/z in consecutive scans
  RTtolerance <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0044'), 2])
  retentionTimeCorrectionCheck <- if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0045'), 2]) == "yes") {TRUE} else {FALSE}
  compoundNames <- gsub('/|[\\]|\t|\n|:|[*]|[?]|<|>|[|]|[.][.][.]|[.][.]|^[.]', '_', ref_table$name, fixed = FALSE)  # To replace invalid characters
  nCompoundNames <- length(compoundNames)
  if (nCompoundNames > 0) {
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
    ############################################################################
    ##
    MAT <- matrix(rep(0, L_PL*26), nrow = L_PL) # 24 + 2 = 26
    MAT <- data.frame(MAT)
    colnames(MAT) <- c("SampleID", "PeakID", "ScanNumberStart","ScanNumberEnd","retentionTimeApex","PeakHeight","PeakArea",
                       "NumberDetectedScans(nIsoPair)","RCS(%)","m/z 12C","CumulatedIntensity", "m/z 13C",
                       "Ratio13C CumulatedIntensity","PeakWidthBaseline","RatioPeakWidth @ 50%",
                       "SeparationTray","peakAsymmetryFactor @ 10%","peakUSPtailingFactor @ 5%",
                       "Skewness_DerivativeMethod", "Symmetry PseudoMoments","Skewness PseudoMoments",
                       "Gaussianity", "S/N baseline", "S/N xcms method", "S/N RMS", "Sharpness")
    ##
    annotation_list <- lapply(1:nCompoundNames, function(i) {MAT})
    ##
    ############################################################################
    ##
    progressBARboundaries <- txtProgressBar(min = 1, max = L_PL, initial = 1, style = 3)
    ##
    if (retentionTimeCorrectionCheck) {
      listCorrectedRTpeaklists <- loadRdata(paste0(output_path, "/peak_alignment/listCorrectedRTpeaklists.Rdata"))
      for (i in 1:L_PL) {
        peaklist <- loadRdata(paste0(inputPathPeaklist, "/", peaklistFileNames[i]))
        RT_i <- listCorrectedRTpeaklists[[peaklistFileNames[i]]]
        for (j in 1:nCompoundNames) {
          annotation_list[[j]][i, 1] <- HRMSfileNames[i]
          x_compound <- mzRTindexer(MZvec = peaklist[, 8], RTvec = RT_i,
                                    MZref = mz_compounds[j], RTref = rt_compounds[j], massAccuracy, RTtolerance)
          if (!is.null(x_compound)) {
            annotation_list[[j]][i, 2:26] <- c(x_compound, peaklist[x_compound, ])
          }
        }
        setTxtProgressBar(progressBARboundaries, i)
      }
      ##
    } else {
      for (i in 1:L_PL) {
        peaklist <- loadRdata(paste0(inputPathPeaklist, "/", peaklistFileNames[i]))
        for (j in 1:nCompoundNames) {
          annotation_list[[j]][i, 1] <- HRMSfileNames[i]
          x_compound <- mzRTindexer(MZvec = peaklist[, 8], RTvec = peaklist[, 3],
                                    MZref = mz_compounds[j], RTref = rt_compounds[j], massAccuracy, RTtolerance)
          if (!is.null(x_compound)) {
            annotation_list[[j]][i, 2:26] <- c(x_compound, peaklist[x_compound, ])
          }
        }
        setTxtProgressBar(progressBARboundaries, i)
      }
    }
    close(progressBARboundaries)
    ##
    ############################################################################
    ##
    IPA_logRecorder("Compound centric annotation data are stored in the `compound_centric_annotation` folder!")
    ##
    outputPathEICs <- paste0(output_path, "/IPA_EIC")
    if (dir.exists(outputPathEICs)) {
      exportEICcheck <- TRUE
      ##
      eicPNGfolderlist <- lapply(1:L_PL, function(j) {
        folderIPAeic <- paste0(outputPathEICs, "/", HRMSfileNames[j])
        if (dir.exists(folderIPAeic)) {
          dirPNG <- dir(path = folderIPAeic, pattern = ".png$")
          xPNG <- as.numeric(do.call(c, lapply(strsplit(dirPNG, "_"), function(k) {k[length(k) - 5]})))
          list(folderIPAeic, dirPNG, xPNG)
        }
      })
      ##
      IPA_logRecorder("If EICs exist, EIC folders are generated for individual compound in the `compound_centric_annotation` folder!")
      ##
    } else {
      exportEICcheck <- FALSE
    }
    ##
    progressBARboundaries <- txtProgressBar(min = 1, max = nCompoundNames, initial = 1, style = 3)
    ##
    for (i in 1:nCompoundNames) {
      A <- annotation_list[[i]]
      nameAnnotatedCompound <- paste0(i, "_cpd__", compoundNames[i])
      write.csv(A, file = paste0(OutputCSV, "/", nameAnnotatedCompound, ".csv"), row.names = TRUE)
      ##
      if (exportEICcheck) {
        peakIDs <- as.numeric(A[, 2])
        xNon0 <- which(peakIDs > 0)
        if (length(xNon0) > 0) {
          sampleName <- A[, 1]
          ##
          folderEICAnnotatedCompound <- paste0(OutputCSV, "/", nameAnnotatedCompound)
          if (!dir.exists(folderEICAnnotatedCompound)) {
            dir.create(folderEICAnnotatedCompound, recursive = TRUE)
          }
          if (dir.exists(folderEICAnnotatedCompound)) {
            for (j in xNon0) {
              eicPNG <- eicPNGfolderlist[[j]]
              if (!is.null(eicPNGfolderlist)) {
                folderIPAeic <-  eicPNG[[1]]
                dirPNG <- eicPNG[[2]]
                xPNG <- eicPNG[[3]]
                xCompoundMatch <- which(xPNG == peakIDs[j])
                sourceFileAddress <- paste0(folderIPAeic, "/", dirPNG[xCompoundMatch])
                destFileAddress <- paste0(folderEICAnnotatedCompound, "/", dirPNG[xCompoundMatch])
                file.copy(from = sourceFileAddress, to = destFileAddress)
              }
            }
          } else {
            IPA_logRecorder(paste0("WARNING: EIC directory can not be created for `", nameAnnotatedCompound,"` compound!!!"))
          }
        }
      }
      setTxtProgressBar(progressBARboundaries, i)
    }
    ##
    close(progressBARboundaries)
    ##
    ############################################################################
    ##
    IPA_logRecorder("Completed compound-centric peak annotation!")
    IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  }
  ##
  ##############################################################################
  ##
  return()
}
