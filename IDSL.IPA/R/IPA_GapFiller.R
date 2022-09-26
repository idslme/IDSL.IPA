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
  ## To perform chromatography analysis
  output_path <- PARAM[which(PARAM[, 1] == 'PARAM0010'), 2]
  input_path_peaklist <- paste0(output_path, "/peaklists")
  file_names_peaklist1 <- dir(path = input_path_peaklist, pattern = ".Rdata")
  file_names_peaklist2 <- dir(path = input_path_peaklist, pattern = "peaklist_")
  file_names_peaklist <- file_names_peaklist1[file_names_peaklist1 %in% file_names_peaklist2]
  L_PL <- length(file_names_peaklist)
  ##
  file_names_peaklist_hrms1 <- gsub(".Rdata", "", file_names_peaklist)
  file_names_peaklist_hrms2 <- gsub("peaklist_", "", file_names_peaklist_hrms1)
  file_names_peaklist_hrms <- file_name_hrms %in% file_names_peaklist_hrms2
  if (length(which(file_names_peaklist_hrms == TRUE)) != L_PL) {
    stop(IPA_logRecorder("Error!!! peaklist files are not available for all selected HRMS files!"))
  }
  ##
  peak_Xcol <- loadRdata(paste0(output_path, "/peak_alignment/peak_Xcol.Rdata"))
  corrected_RT_peaklists <- loadRdata(paste0(output_path, "/peak_alignment/corrected_RT_peaklists.Rdata"))
  ##
  massDifferenceIsotopes <- tryCatch(as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0012'), 2]), error = function(e) {1.003354835336}, warning = function(w) {1.003354835336})     # Mass difference for isotopic pairs
  ##
  mass_error <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0038'), 2])   # Mass accuracy to cluster m/z in consecutive scans
  mass_error_13c <- 1.5*mass_error
  delta_rt <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0039'), 2])
  scan_tol <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0040'), 2])
  ##
  osType <- Sys.info()[['sysname']]
  if (osType == "Linux") {
    progressBARboundaries <- txtProgressBar(min = 1, max = L_HRMS, initial = 1, style = 3)
    chromatography_undetected_list <- lapply(1:L_HRMS, function(i) {
      setTxtProgressBar(progressBARboundaries, i)
      ##
      chromatography_undetected <- NULL
      x_0 <- which(peak_Xcol[, (i + 2)] == 0)
      L_x_0 <- length(x_0)
      if (L_x_0 > 0) {
        mz_Xcol <- peak_Xcol[x_0, 1]
        ## To back calculate the RT ##
        undeteced_RT <- peak_Xcol[x_0, 2]
        uncorrected_RTi <- matrix(loadRdata(paste0(input_path_peaklist, "/", file_names_peaklist[i]))[, 3], ncol = 1)
        corrected_RTi <- corrected_RT_peaklists[[i]]
        idf <- data.frame(uncoRT = uncorrected_RTi, coRT = corrected_RTi)
        rtmodel <- lm(coRT ~ poly(uncoRT, 3), idf)
        new.df <- data.frame(uncoRT = undeteced_RT) # predict uncorrected RTs
        RT_uncorrected_undeteced <- matrix(predict(rtmodel, new.df), ncol = 1)
        ##
        outputer <- IPA_MSdeconvoluter(input_path_hrms, file_name_hrms[i])
        spectraList <- outputer[["spectraList"]]
        RetentionTime <- outputer[["retentionTime"]]
        nRT <- length(RetentionTime)
        ##
        chromatography_undetected <- do.call(rbind, mclapply(1:L_x_0, function(j) {
          chromatography_undetected_row <- NULL
          pa <- 0
          R13C <- 0
          mzCandidate <- mz_Xcol[j]
          rtCandidate <- RT_uncorrected_undeteced[j]
          sn_apex <- which.min(abs(rtCandidate - RetentionTime)) # scan number at apex
          ScanNumberStart <- max(c((sn_apex - scan_tol), 1))
          ScanNumberEnd <- min(c(nRT, (sn_apex + scan_tol)))
          spectraList.xic <- spectraList[ScanNumberStart:ScanNumberEnd]
          chrom_builder <- XIC(spectraList.xic, scan_number_start = ScanNumberStart, mzCandidate, mass_error)
          length_chrom <- dim(chrom_builder)[1]
          x_apex <- which(chrom_builder[, 1] == sn_apex)
          rt_loc_min <- islocalminimum(chrom_builder[, 3])
          boundary_left <- which(rt_loc_min[1:x_apex] == -1)
          if (length(boundary_left) > 0) {
            boundary_left <- boundary_left[length(boundary_left)]
          } else {
            boundary_left <- 1
          }
          boundary_right <- which(rt_loc_min[x_apex:length_chrom] == -1)
          if (length(boundary_right) > 0) {
            boundary_right <- boundary_right[1] - 1 + x_apex
          } else {
            boundary_right <- length_chrom
          }
          chrom <- cbind(RetentionTime[chrom_builder[boundary_left:boundary_right, 1]], chrom_builder[boundary_left:boundary_right, 3])
          RT_detected <- chrom[which.min(abs(chrom[, 1] - rtCandidate)), 1]
          if (abs(RT_detected - rtCandidate) <= delta_rt) {
            height <- max(chrom[, 2])
            ## R13C
            t1 <- boundary_left + ScanNumberStart - 1
            t2 <- boundary_right + ScanNumberStart - 1
            chromatogram_segment <- do.call(rbind, lapply(t1:t2, function(t) {
              Spec_ScN_j <- NULL
              Spec <- spectraList[[t]]
              if (length(Spec) > 0) {
                x_mz1 <- which(abs(Spec[, 1] - mzCandidate) <= mass_error)
                if (length(x_mz1) > 0) {
                  x_mz2 <- which(abs(Spec[, 1] - (massDifferenceIsotopes + mzCandidate)) <= mass_error_13c)
                  if (length(x_mz2) > 0) {
                    if (length(x_mz1) > 1) {
                      x_min <- which.min(abs(Spec[x_mz1, 1] - mzCandidate))
                      x_mz1 <- x_mz1[x_min]
                    }
                    if (length(x_mz2) > 1) {
                      x_min <- which.min(abs(Spec[x_mz2, 1] - (massDifferenceIsotopes + mzCandidate)))
                      x_mz2 <- x_mz2[x_min]
                    }
                    if (Spec[x_mz1, 2] >= Spec[x_mz2, 2]) {
                      Spec_ScN_j <- c(Spec[x_mz1, 2], Spec[x_mz2, 2])
                    }
                  }
                }
              }
              Spec_ScN_j
            }))
            if (length(chromatogram_segment) > 0) {
              Int12C <- sum(chromatogram_segment[, 1])
              Int13C <- sum(chromatogram_segment[, 2])
              R13C <- Int13C/Int12C*100
            }
            if (boundary_left != boundary_right) {
              pa <- peak_area(chrom[, 1], chrom[, 2])
            }
            chromatography_undetected_row <- c(x_0[j], height, pa, R13C)
          }
          chromatography_undetected_row
        }, mc.cores = number_processing_threads))
      }
      chromatography_undetected
    })
    closeAllConnections()
    close(progressBARboundaries)
    ##
  } else if (osType == "Windows") {
    ##
    clust <- makeCluster(number_processing_threads)
    registerDoParallel(clust)
    chromatography_undetected_list <- foreach(i=1:L_HRMS, .verbose = FALSE) %dopar% {
      chromatography_undetected <- NULL
      x_0 <- which(peak_Xcol[, (i + 2)] == 0)
      L_x_0 <- length(x_0)
      if (L_x_0 > 0) {
        mz_Xcol <- peak_Xcol[x_0, 1]
        ## To back calculate the RT ##
        undeteced_RT <- peak_Xcol[x_0, 2]
        uncorrected_RTi <- matrix(loadRdata(paste0(input_path_peaklist, "/", file_names_peaklist[i]))[, 3], ncol = 1)
        corrected_RTi <- corrected_RT_peaklists[[i]]
        idf <- data.frame(uncoRT = uncorrected_RTi, coRT = corrected_RTi)
        rtmodel <- lm(coRT ~ poly(uncoRT, 3), idf)
        new.df <- data.frame(uncoRT = undeteced_RT) # predict uncorrected RTs
        RT_uncorrected_undeteced <- matrix(predict(rtmodel, new.df), ncol = 1)
        ##
        outputer <- IPA_MSdeconvoluter(input_path_hrms, file_name_hrms[i])
        spectraList <- outputer[["spectraList"]]
        RetentionTime <- outputer[["retentionTime"]]
        nRT <- length(RetentionTime)
        chromatography_undetected <- do.call(rbind, lapply(1:L_x_0, function(j) {
          chromatography_undetected_row <- NULL
          pa <- 0
          R13C <- 0
          mzCandidate <- mz_Xcol[j]
          rtCandidate <- RT_uncorrected_undeteced[j]
          sn_apex <- which.min(abs(rtCandidate - RetentionTime)) # scan number at apex
          ScanNumberStart <- max(c((sn_apex - scan_tol), 1))
          ScanNumberEnd <- min(c(nRT, (sn_apex + scan_tol)))
          spectraList.xic <- spectraList[ScanNumberStart:ScanNumberEnd]
          chrom_builder <- XIC(spectraList.xic, scan_number_start = ScanNumberStart, mzCandidate, mass_error)
          length_chrom <- dim(chrom_builder)[1]
          x_apex <- which(chrom_builder[,1] == sn_apex)
          rt_loc_min <- islocalminimum(chrom_builder[, 3])
          boundary_left <- which(rt_loc_min[1:x_apex] == -1)
          if (length(boundary_left) > 0) {
            boundary_left <- boundary_left[length(boundary_left)]
          } else {
            boundary_left <- 1
          }
          boundary_right <- which(rt_loc_min[x_apex:length_chrom] == -1)
          if (length(boundary_right) > 0) {
            boundary_right <- boundary_right[1] - 1 + x_apex
          } else {
            boundary_right <- length_chrom
          }
          chrom <- cbind(RetentionTime[chrom_builder[boundary_left:boundary_right, 1]], chrom_builder[boundary_left:boundary_right, 3])
          RT_detected <- chrom[which.min(abs(chrom[, 1] - rtCandidate)), 1]
          if (abs(RT_detected - rtCandidate) <= delta_rt) {
            height <- max(chrom[, 2])
            ## R13C
            t1 <- boundary_left + ScanNumberStart - 1
            t2 <- boundary_right + ScanNumberStart - 1
            chromatogram_segment <- do.call(rbind, lapply(t1:t2, function(t) {
              Spec_ScN_j <- NULL
              Spec <- spectraList[[t]]
              if (length(Spec) > 0) {
                x_mz1 <- which(abs(Spec[, 1] - mzCandidate) <= mass_error)
                if (length(x_mz1) > 0) {
                  x_mz2 <- which(abs(Spec[, 1] - (massDifferenceIsotopes + mzCandidate)) <= mass_error_13c)
                  if (length(x_mz2) > 0) {
                    if (length(x_mz1) > 1) {
                      x_min <- which.min(abs(Spec[x_mz1, 1] - mzCandidate))
                      x_mz1 <- x_mz1[x_min]
                    }
                    if (length(x_mz2) > 1) {
                      x_min <- which.min(abs(Spec[x_mz2, 1] - (massDifferenceIsotopes + mzCandidate)))
                      x_mz2 <- x_mz2[x_min]
                    }
                    if (Spec[x_mz1, 2] >= Spec[x_mz2, 2]) {
                      Spec_ScN_j <- c(Spec[x_mz1, 2], Spec[x_mz2, 2])
                    }
                  }
                }
              }
              Spec_ScN_j
            }))
            if (length(chromatogram_segment) > 0) {
              Int12C <- sum(chromatogram_segment[, 1])
              Int13C <- sum(chromatogram_segment[, 2])
              R13C <- Int13C/Int12C*100
            }
            if (boundary_left != boundary_right) {
              pa <- peak_area(chrom[, 1], chrom[, 2])
            }
            chromatography_undetected_row <- c(x_0[j], height, pa, R13C)
          }
          chromatography_undetected_row
        }))
      }
      chromatography_undetected
    }
    stopCluster(clust)
  }
  ##
  IPA_logRecorder("Initiated filling gaps of the aligned peak height, peak area, and R13C tables!")
  OutputPath_peak_alignment <- paste0(output_path, "/peak_alignment/")
  peak_height <- loadRdata(paste0(OutputPath_peak_alignment, "peak_height.Rdata"))
  peak_area <- loadRdata(paste0(OutputPath_peak_alignment, "peak_area.Rdata"))
  peak_R13C <- loadRdata(paste0(OutputPath_peak_alignment, "peak_R13C.Rdata"))
  peak_height_gapfilled <- peak_height
  peak_area_gapfilled <- peak_area
  peak_R13C_gapfilled <- peak_R13C
  progressBARboundaries <- txtProgressBar(min = 0, max = L_HRMS, initial = 0, style = 3)
  for (i in 1:length(chromatography_undetected_list)) {
    setTxtProgressBar(progressBARboundaries, i)
    iSample <- chromatography_undetected_list[[i]]
    if (length(iSample) > 0) {
      x_j <- iSample[, 1]
      counter_j <- 0
      for (j in x_j) {
        counter_j <- counter_j + 1
        peak_height_gapfilled[j, (i + 2)] <- round(iSample[counter_j, 2], 0)
        peak_area_gapfilled[j, (i + 2)] <- round(iSample[counter_j, 3], 0)
        peak_R13C_gapfilled[j, (i + 2)] <- round(iSample[counter_j, 4], 3)
      }
    }
  }
  close(progressBARboundaries)
  opendir(OutputPath_peak_alignment)
  IPA_logRecorder("Initiated saving aligned gap-filled peak tables!")
  save(peak_height_gapfilled, file = paste0(OutputPath_peak_alignment, "peak_height_gapfilled.Rdata"))
  write.csv(peak_height_gapfilled, file = paste0(OutputPath_peak_alignment, "peak_height_gapfilled.csv"))
  save(peak_area_gapfilled, file = paste0(OutputPath_peak_alignment, "peak_area_gapfilled.Rdata"))
  write.csv(peak_area_gapfilled, file = paste0(OutputPath_peak_alignment, "peak_area_gapfilled.csv"))
  save(peak_R13C_gapfilled, file = paste0(OutputPath_peak_alignment, "peak_R13C_gapfilled.Rdata"))
  write.csv(peak_R13C_gapfilled, file = paste0(OutputPath_peak_alignment, "peak_R13C_gapfilled.csv"))
  IPA_logRecorder("Gap-filled aligned peak height, peak area, and R13C tables were stored in `.Rdata` and `.csv` formats in the `peak_alignment` folder!")
  ##
  ############################################################################
  ##
  IPA_logRecorder("Initiated detecting correlating peaks on the gap-filled peak height table!")
  ##
  correlationListHeight_gapfilled <- peak_property_table_correlation(peakPropertyTable = peak_height_gapfilled, RTtolerance = delta_rt, minFreqDetection = 1, method = "pearson", minThresholdCorrelation = 0.50, number_processing_threads)
  ##
  save(correlationListHeight_gapfilled, file = paste0(OutputPath_peak_alignment, "/correlationListHeight_gapfilled.Rdata"))
  ##
  IPA_logRecorder("Stored the correlating peaks on the gap-filled peak height table as `correlationListHeight_gapfilled.Rdata` in the `peak_alignement` folder!")
  IPA_logRecorder("Completed gap-filling!")
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  ##
  ############################################################################
  ##
}
