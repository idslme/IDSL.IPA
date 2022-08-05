IPA_PeaklistAnnotation <- function(PARAM) {
  print("Initiated sample-centric peak annotation!")
  number_processing_threads <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0006'), 2])
  output_path <- PARAM[which(PARAM[, 1] == 'PARAM0010'), 2]
  Output_Xcol <- paste0(output_path, "/sample_centeric_annotation")
  ##
  massDifferenceIsotopes <- tryCatch(as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0012'), 2]), error = function(e) {1.003354835336}, warning = function(w) {1.003354835336})     # Mass difference for isotopic pairs
  ##
  ref_table <- readxl::read_xlsx(PARAM[which(PARAM[, 1] == 'PARAM0042'), 2])
  mass_error <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0043'), 2])   # Mass accuracy to cluster m/z in consecutive scans
  rt_error <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0044'), 2])
  x0045 <- (tolower(PARAM[which(PARAM[, 1] == 'PARAM0045'), 2]) == "yes")
  x0048 <- (tolower(PARAM[which(PARAM[, 1] == 'PARAM0048'), 2]) == "yes")
  name_compounds <- ref_table$name
  L_nc <- length(name_compounds)
  if (L_nc > 0) {
    mz_compounds <- ref_table$'m/z'
    rt_compounds <- ref_table$RT
    #
    input_path_peaklist <- paste0(output_path, "/peaklists")
    file_names_peaklist1 <- dir(path = input_path_peaklist, pattern = ".Rdata")
    file_names_peaklist2 <- dir(path = input_path_peaklist, pattern = "peaklist_")
    file_names_peaklist <- file_names_peaklist1[file_names_peaklist1 %in% file_names_peaklist2]
    L_PL <- length(file_names_peaklist)
    ##
    input_path_hrms <- PARAM[which(PARAM[, 1] == 'PARAM0007'), 2]
    if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]) == "all") {
      file_name_hrms <- dir(path = input_path_hrms)
      file_name_hrms <- file_name_hrms[grep(paste0(".", tolower(PARAM[which(PARAM[, 1] == 'PARAM0009'), 2]), "$"), file_name_hrms, ignore.case = TRUE)]
    } else {
      samples_string <- PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]
      file_name_hrms <- strsplit(samples_string, ";")[[1]] # files used as reference m/z-RT
    }
    file_names_peaklist_hrms1 <- gsub(".Rdata", "", file_names_peaklist)
    file_names_peaklist_hrms2 <- gsub("peaklist_", "", file_names_peaklist_hrms1)
    file_names_peaklist_hrms <- file_name_hrms %in% file_names_peaklist_hrms2
    if (length(which(file_names_peaklist_hrms == TRUE)) != L_PL) {
      stop("Error!!! peaklist files are not available for all selected HRMS files!")
    }
    ##
    annotated_peak_indices <- matrix(rep(0, L_nc*L_PL), nrow = L_nc)
    progressBARboundaries <- txtProgressBar(min = 1, max = L_PL, initial = 1, style = 3)
    if (x0045 == TRUE) {
      corrected_RT_peaklists <- loadRdata(paste0(output_path, "/peak_alignment/corrected_RT_peaklists.Rdata"))
      ##
      for (i in 1:L_PL) {
        setTxtProgressBar(progressBARboundaries, i)
        peaklist <- loadRdata(paste0(input_path_peaklist, "/", file_names_peaklist[i]))
        mz_i <- matrix(peaklist[, 8], ncol = 1)
        RT_i <- corrected_RT_peaklists[[i]]
        for (j in 1:L_nc) {
          x_compound <- which(abs(mz_compounds[j] - mz_i) <= mass_error &
                                abs(rt_compounds[j] - RT_i) <= rt_error)
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
      for (i in 1:L_PL) {
        setTxtProgressBar(progressBARboundaries, i)
        peaklist <- loadRdata(paste0(input_path_peaklist, "/", file_names_peaklist[i]))
        mz_i <- matrix(peaklist[, 8], ncol = 1)
        RT_i <- matrix(peaklist[, 3], ncol = 1)
        for (j in 1:L_nc) {
          x_compound <- which(abs(mz_compounds[j] - mz_i) <= mass_error &
                                abs(rt_compounds[j] - RT_i) <= rt_error)
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
    annotated_peak_indices <- cbind(mz_compounds, rt_compounds , annotated_peak_indices)
    file_names_hrms <- gsub(".Rdata", "", (gsub("peaklist_", "", file_names_peaklist)))
    colnames(annotated_peak_indices) <- c("m/z", "RT", file_names_hrms)
    ##
    dir.create(Output_Xcol)
    opendir(Output_Xcol)
    save(annotated_peak_indices, file = paste0(Output_Xcol, "/annotated_peak_indices.Rdata"))
    print("Peak index numbers from individual peaklists were stored in 'annotated_peak_indices.Rdata'")
    listHeightAreaR13C <- peak_Xcol2(input_path_peaklist, file_names_peaklist, annotated_peak_indices)
    annotated_peak_height <- cbind(name_compounds, listHeightAreaR13C[["peak_height"]])
    colnames(annotated_peak_height) <- c("name", "m/z", "RT", file_names_hrms)
    annotated_peak_area <- cbind(name_compounds, listHeightAreaR13C[["peak_area"]])
    colnames(annotated_peak_area) <- c("name", "m/z", "RT", file_names_hrms)
    annotated_peak_R13C <- cbind(name_compounds, listHeightAreaR13C[["peak_R13C"]])
    colnames(annotated_peak_R13C) <- c("name", "m/z", "RT", file_names_hrms)
    write.csv(annotated_peak_height, file = paste0(Output_Xcol, "/annotated_peak_height.csv"))
    write.csv(annotated_peak_area, file = paste0(Output_Xcol, "/annotated_peak_area.csv"))
    write.csv(annotated_peak_R13C, file = paste0(Output_Xcol, "/annotated_peak_R13C.csv"))
    print("Completed sample-centric peak annotation for peak height, peak area, and R13C!")
    ##
    if (x0048 == TRUE) {
      print("Initiated gap-filling for sample-centric peak annotation")
      mass_error <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0038'), 2])   # Mass accuracy to cluster m/z in consecutive scans
      mass_error_13c <- 1.5*mass_error
      delta_rt <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0039'), 2])
      scan_tol <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0040'), 2])
      ##
      input_path_hrms <- PARAM[which(PARAM[, 1] == 'PARAM0007'), 2]
      if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]) == "all") {
        file_name_hrms <- dir(path = input_path_hrms)
        file_name_hrms <- file_name_hrms[grep(paste0(".", tolower(PARAM[which(PARAM[, 1] == 'PARAM0009'), 2]), "$"), file_name_hrms, ignore.case = TRUE)]
      } else {
        samples_string <- PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]
        file_name_hrms <- strsplit(samples_string, ";")[[1]] # files used as reference m/z-RT
      }
      L_HRMS <- length(file_name_hrms)
      ##
      osType <- Sys.info()[['sysname']]
      if (osType == "Linux") {
        chromatography_undetected_list <- mclapply(1:L_HRMS, function(i) {
          ##
          chromatography_undetected <- NULL
          x_0 <- which(annotated_peak_indices[, (i + 2)] == 0)
          if (length(x_0) > 0) {
            mz_Xcol <- annotated_peak_indices[x_0, 1]
            undeteced_RT <- annotated_peak_indices[x_0, 2]
            if (x0045 == TRUE) {
              ## To back calculate the RT ##
              uncorrected_RTi <- matrix(loadRdata(paste0(input_path_peaklist, "/", file_names_peaklist[i]))[, 3], ncol = 1)
              corrected_RTi <- corrected_RT_peaklists[[i]]
              idf <- data.frame(uncoRT = uncorrected_RTi, coRT = corrected_RTi)
              rtmodel <- lm(coRT ~ poly(uncoRT, 3), idf)
              new.df <- data.frame(uncoRT = undeteced_RT) # predict uncorrected RTs
              RT_uncorrected_undeteced <- matrix(predict(rtmodel, new.df), ncol = 1)
            } else {
              RT_uncorrected_undeteced <- matrix(undeteced_RT, ncol = 1)
            }
            ##
            outputer <- IPA_MSdeconvoluter(input_path_hrms, file_name_hrms[i])
            spectraList <- outputer[["spectraList"]]
            RetentionTime <- outputer[["retentionTime"]]
            nRT <- length(RetentionTime)
            ##
            chromatography_undetected <- do.call(rbind, lapply(1:length(x_0), function(j) {
              chromatography_undetected_row <- NULL
              pa <- 0
              R13C <- 0
              mzCandidate <- mz_Xcol[j]
              rtCandidate <- RT_uncorrected_undeteced[j]
              sn_apex <- which.min(abs(rtCandidate - RetentionTime)) # scan number at apex
              ScanNumberStart <- max(c((sn_apex - scan_tol), 1))
              ScanNumberEnd <- min(c(nRT, (sn_apex + scan_tol)))
              spectraList.xic <- spectraList[ScanNumberStart:ScanNumberEnd]
              chrom_builder <- XIC(spectraList.xic, ScanNumberStart, mzCandidate, mass_error)
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
            }))
          }
          chromatography_undetected
        }, mc.cores = number_processing_threads)
        closeAllConnections()
        ##
      } else if (osType == "Windows") {
        ##
        clust <- makeCluster(number_processing_threads)
        registerDoParallel(clust)
        chromatography_undetected_list <- foreach(i=1:L_HRMS, .verbose = FALSE) %dopar% {
          ##
          chromatography_undetected <- NULL
          x_0 <- which(annotated_peak_indices[, (i + 2)] == 0)
          if (length(x_0) > 0) {
            mz_Xcol <- annotated_peak_indices[x_0, 1]
            undeteced_RT <- annotated_peak_indices[x_0, 2]
            if (x0045 == TRUE) {
              ## To back calculate the RT ##
              uncorrected_RTi <- matrix(loadRdata(paste0(input_path_peaklist, "/", file_names_peaklist[i]))[, 3], ncol = 1)
              corrected_RTi <- corrected_RT_peaklists[[i]]
              idf <- data.frame(uncoRT = uncorrected_RTi, coRT = corrected_RTi)
              rtmodel <- lm(coRT ~ poly(uncoRT, 3), idf)
              new.df <- data.frame(uncoRT = undeteced_RT) # predict uncorrected RTs
              RT_uncorrected_undeteced <- matrix(predict(rtmodel, new.df), ncol = 1)
            } else {
              RT_uncorrected_undeteced <- matrix(undeteced_RT, ncol = 1)
            }
            ##
            outputer <- IPA_MSdeconvoluter(input_path_hrms, file_name_hrms[i])
            spectraList <- outputer[["spectraList"]]
            RetentionTime <- outputer[["retentionTime"]]
            nRT <- length(RetentionTime)
            ##
            chromatography_undetected <- do.call(rbind, lapply(1:length(x_0), function(j) {
              chromatography_undetected_row <- NULL
              pa <- 0
              R13C <- 0
              mzCandidate <- mz_Xcol[j]
              rtCandidate <- RT_uncorrected_undeteced[j]
              sn_apex <- which.min(abs(rtCandidate - RetentionTime)) # scan number at apex
              ScanNumberStart <- max(c((sn_apex - scan_tol), 1))
              ScanNumberEnd <- min(c(nRT, (sn_apex + scan_tol)))
              spectraList.xic <- spectraList[ScanNumberStart:ScanNumberEnd]
              chrom_builder <- XIC(spectraList.xic, ScanNumberStart, mzCandidate, mass_error)
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
            }))
          }
          chromatography_undetected
        }
        stopCluster(clust)
      }
      ##
      annotated_peak_height_gapfilled <- annotated_peak_height
      annotated_peak_area_gapfilled <- annotated_peak_area
      annotated_peak_R13C_gapfilled <- annotated_peak_R13C
      progressBARboundaries <- txtProgressBar(min = 1, max = L_HRMS, initial = 1, style = 3)
      for (i in 1:length(chromatography_undetected_list)) {
        setTxtProgressBar(progressBARboundaries, i)
        iSample <- chromatography_undetected_list[[i]]
        if (length(iSample) > 0) {
          x_j <- iSample[, 1]
          counter_j <- 0
          for (j in x_j) {
            counter_j <- counter_j + 1
            annotated_peak_height_gapfilled[j, (i + 3)] <- round(iSample[counter_j, 2], 0)
            annotated_peak_area_gapfilled[j, (i + 3)] <- round(iSample[counter_j, 3], 0)
            annotated_peak_R13C_gapfilled[j, (i + 3)] <- round(iSample[counter_j, 4], 3)
          }
        }
      }
      close(progressBARboundaries)
      print("Initiated saving gap-filled sample-centric peak annotation!")
      write.csv(annotated_peak_height_gapfilled, file = paste0(Output_Xcol, "/annotated_peak_height_gapfilled.csv"))
      write.csv(annotated_peak_area_gapfilled, file = paste0(Output_Xcol, "/annotated_peak_area_gapfilled.csv"))
      write.csv(annotated_peak_R13C_gapfilled, file = paste0(Output_Xcol, "/annotated_peak_R13C_gapfilled.csv"))
      print("Completed gap-filled sample-centric peak annotation!")
    }
  }
}
