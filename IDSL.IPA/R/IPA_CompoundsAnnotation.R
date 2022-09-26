IPA_CompoundsAnnotation <- function(PARAM) {
  ##
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  IPA_logRecorder("Initiated compound-centric peak annotation!")
  ##
  output_path <- PARAM[which(PARAM[, 1] == 'PARAM0010'), 2]
  Output_CSV <- paste0(output_path, "/compound_centric_annotation")
  if (!dir.exists(Output_CSV)) {
    dir.create(Output_CSV, recursive = TRUE)
  }
  opendir(Output_CSV)
  ##
  IPA_logRecorder("Compound centric annotation data are stored in the `compound_centric_annotation` folder!")
  ##
  ref_table <- readxl::read_xlsx(PARAM[which(PARAM[, 1] == 'PARAM0042'), 2])
  mass_error <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0043'), 2])   # Mass accuracy to cluster m/z in consecutive scans
  rt_error <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0044'), 2])
  x0045 <- PARAM[which(PARAM[, 1] == 'PARAM0045'), 2]
  name_compounds <- ref_table$name
  name_compounds <- IPA_gsub(c("/", "\\", ":", "*", "?", '"', "<", ">", "|"), "_", name_compounds, fixed = TRUE)  # To replace invalid characters
  L_nc <- length(name_compounds)
  if (L_nc > 0) {
    mz_compounds <- ref_table$`m/z`
    rt_compounds <- ref_table$RT
    #
    input_path_peaklist <- paste0(output_path, "/peaklists")
    file_names_peaklist1 <- dir(path = input_path_peaklist, pattern = ".Rdata")
    file_names_peaklist2 <- dir(path = input_path_peaklist, pattern = "peaklist_")
    file_names_peaklist <- file_names_peaklist1[file_names_peaklist1 %in% file_names_peaklist2]
    L_PL <- length(file_names_peaklist)
    ##
    annotation_list <- lapply(1:L_nc, function(i) {
      MAT <- matrix(rep(0, L_PL*25), nrow = L_PL) # 24 + 1 = 25
      MAT <- data.frame(MAT)
      names(MAT) <- c("SampleID","ScanNumberStart","ScanNumberEnd","RetentionTimeApex","PeakHeight","PeakArea",
                      "NumberDetectedScans(nIsoPair)","RCS(%)","m/z MonoIsotopic","CumulatedIntensity",
                      "m/z 13C","Ratio 13C CumulatedIntensity","PeakWidthBaseline","Ratio PeakWidth @ 50%",
                      "SeperationTray","AsymmetryFactor @ 10%","USPTailingFactor @ 5%",
                      "Skewness_DerivativeMethod", "Symmetry PseudoMoments","Skewness PseudoMoments",
                      "Gaussianity", "S/N", "S/N xcms method", "S/N RMS", "Sharpness")
      MAT
    })
    progressBARboundaries <- txtProgressBar(min = 1, max = L_PL, initial = 1, style = 3)
    ##
    if (tolower(x0045) == "yes") {
      corrected_RT_peaklists <- loadRdata(paste0(output_path, "/peak_alignment/corrected_RT_peaklists.Rdata"))
      for (i in 1:L_PL) {
        setTxtProgressBar(progressBARboundaries, i)
        peaklist <- loadRdata(paste(input_path_peaklist, "/", file_names_peaklist[i], sep = ""))
        S_ID <- gsub("peaklist_", "", file_names_peaklist[i])
        S_ID <- gsub(".Rdata", "", S_ID)
        mz_i <- matrix(peaklist[, 8], ncol = 1)
        RT_i <- corrected_RT_peaklists[[i]]
        for (j in 1:L_nc) {
          annotation_list[[j]][i, 1] <- S_ID
          x_compound <- which(abs(mz_compounds[j] - mz_i) <= mass_error &
                                abs(rt_compounds[j] - RT_i) <= rt_error)
          if (length(x_compound) > 0) {
            if (length(x_compound) > 1) {
              x_min <- which.min(abs(rt_compounds[j] - RT_i[x_compound]))
              x_compound <- x_compound[x_min[1]]
            }
            annotation_list[[j]][i, 2:25] <- peaklist[x_compound, ]
          }
        }
      }
    } else {
      for (i in 1:L_PL) {
        setTxtProgressBar(progressBARboundaries, i)
        peaklist <- loadRdata(paste(input_path_peaklist, "/", file_names_peaklist[i], sep = ""))
        S_ID <- gsub("peaklist_", "", file_names_peaklist[i])
        S_ID <- gsub(".Rdata", "", S_ID)
        mz_i <- matrix(peaklist[, 8], ncol = 1)
        RT_i <- matrix(peaklist[, 3], ncol = 1)
        for (j in 1:L_nc) {
          annotation_list[[j]][i, 1] <- S_ID
          x_compound <- which(abs(mz_compounds[j] - mz_i) <= mass_error &
                                abs(rt_compounds[j] - RT_i) <= rt_error)
          if (length(x_compound) > 0) {
            if (length(x_compound) > 1) {
              x_min <- which.min(abs(rt_compounds[j] - RT_i[x_compound]))
              x_compound <- x_compound[x_min[1]]
            }
            annotation_list[[j]][i, 2:25] <- peaklist[x_compound, ]
          }
        }
      }
    }
    close(progressBARboundaries)
    ##
    for (i in 1:L_nc) {
      A <- annotation_list[[i]]
      write.csv(A, file = paste0(Output_CSV, "/", i, "_", "annotated_compound__", name_compounds[i], ".csv"))
    }
    ##
    IPA_logRecorder("Annotated tables for each compound were stored in `.csv` formats in the `compound_centric_annotation` folder!")
    IPA_logRecorder("Completed compound-centric peak annotation!")
    IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  }
}
