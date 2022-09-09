IPA_TargetedAnalysis <- function(spreadsheet, mzCandidate, rtCandidate, exportEIC = TRUE, exportTable = FALSE) {
  cc_table <- c()
  if (length(mzCandidate) !=  length(rtCandidate)) {
    stop("Error!!! mz and rt vectors do not have the same length!")
  }
  PARAM <- xlsxAnalyzer_EIC(spreadsheet)
  if (length(PARAM) > 0) {
    ##
    number_processing_threads <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0006'), 2])
    input_path_hrms <- PARAM[which(PARAM[, 1] == 'PARAM0007'), 2]
    if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]) == "all") {
      file_name_hrms <- dir(path = input_path_hrms)
      file_name_hrms <- file_name_hrms[grep(paste0(".", tolower(PARAM[which(PARAM[, 1] == 'PARAM0009'), 2]), "$"), file_name_hrms, ignore.case = TRUE)]
    } else {
      samples_string <- PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]
      file_name_hrms <- strsplit(samples_string, ";")[[1]] # files used as reference m/z-RT
    }
    output_path <- PARAM[which(PARAM[, 1] == 'PARAM0010'), 2]
    if (!dir.exists(output_path)) {
      print("Created output directory!")
      dir.create(output_path, recursive = TRUE)
    }
    if (exportEIC == TRUE) {
      output_path_eic <- paste0(output_path, "/EICs")
      if (!dir.exists(output_path_eic)) {
        dir.create(output_path_eic, recursive = TRUE)
      }
      opendir(output_path_eic)
      ##
      img <- png::readPNG(paste0(system.file("extdata", package = "IDSL.IPA"),"/EIC_legend.png"))
      legend_EIC <- grid::rasterGrob(img, interpolate = TRUE)
    }
    ##
    massDifferenceIsotopes <- tryCatch(as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0012'), 2]), error = function(e) {1.003354835336}, warning = function(w) {1.003354835336})     # Mass difference for isotopic pairs
    mass_accuracy_xic <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0013'), 2]) # Mass accuracy to cluster m/z in consecutive scans
    mass_accuracy_isotope_pair <- 1.5*mass_accuracy_xic      # Mass accuracy to find 13C isotopologues
    smoothing_window <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0015'), 2])
    peak_resolving_power <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0017'), 2])
    scan_tol <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0020'), 2])		            # Number of scans to include in the search before and after of boundaries of detected peaks
    n_spline <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0028'), 2])         		  # Level of peak smoothing to calculate ancillary chromatography parameters
    ##
    R0 <- c(1, 2, 4, 5, 6, 7, 9, 14, 24)
    R3 <- c(3, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23)
    R5 <- c(8, 10)
    ##
    osType <- Sys.info()[['sysname']]
    print("Initiated the targeted analysis!")
    if (osType == "Windows") {
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      cc_table <- foreach(i = 1:length(file_name_hrms), .combine ='rbind', .verbose = FALSE) %dopar% {
        ##
        outputer <- IPA_MSdeconvoluter(input_path_hrms, file_name_hrms[i])
        spectraList <- outputer[["spectraList"]]
        RetentionTime <- outputer[["retentionTime"]]
        nRT <- length(RetentionTime)
        ##
        chrome <- do.call(rbind, lapply(1:length(mzCandidate), function(j) {
          ScanNumberApex <- which.min(abs(RetentionTime - rtCandidate[j]))
          t1 <- ScanNumberApex - scan_tol
          if (t1 < 1) {
            t1 <- 1
          }
          t2 <- ScanNumberApex + scan_tol
          if (t2 > nRT) {
            t2 <- nRT
          }
          ##
          chromatogram_segment <- do.call(rbind, lapply(t1:t2, function(t) {
            Spec_ScN_j <- NULL
            Spec <- spectraList[[t]]
            if (length(Spec) > 0) {
              x_mz1 <- which(abs(Spec[, 1] - mzCandidate[j]) <= mass_accuracy_xic)
              L_x_mz1 <- length(x_mz1)
              if (L_x_mz1 > 0) {
                x_mz2 <- which(abs(Spec[, 1] - (massDifferenceIsotopes + mzCandidate[j])) <= mass_accuracy_isotope_pair)
                L_x_mz2 <- length(x_mz2)
                if (L_x_mz2 > 0) {
                  if (L_x_mz1 > 1) {
                    x_min <- which.min(abs(Spec[x_mz1, 1] - mzCandidate[j]))
                    x_mz1 <- x_mz1[x_min[1]]
                  }
                  if (L_x_mz2 > 1) {
                    x_min <- which.min(abs(Spec[x_mz2, 1] - (massDifferenceIsotopes + mzCandidate[j])))
                    x_mz2 <- x_mz2[x_min[1]]
                  }
                  if (Spec[x_mz1, 2] >= Spec[x_mz2, 2]) {
                    Spec_ScN_j <- c(Spec[x_mz1, 1], Spec[x_mz1, 2], t, Spec[x_mz2, 1], Spec[x_mz2, 2])
                  }
                }
              }
            }
            Spec_ScN_j
          }))
          if (length(chromatogram_segment) != 0) {
            chromatogram_segment <- chromatogram_segment[order(chromatogram_segment[, 3]), ]
            chromatogram_segment <- matrix(chromatogram_segment, ncol = 5)
          } else {
            chromatogram_segment <- matrix(c(mzCandidate[j], mzCandidate[j], 0 , 0 , t1, t2, 0, 0, 0, 0), ncol = 5)
          }
          if (chromatogram_segment[1, 3] != t1) {
            chromatogram_segment <- rbind(c(mzCandidate[j], 0, t1, 0, 0), chromatogram_segment)
          }
          if (chromatogram_segment[dim(chromatogram_segment)[1], 3] != t2) {
            chromatogram_segment <- rbind(chromatogram_segment, c(mzCandidate[j], 0, t2, 0, 0))
          }
          peak_property <- chromatography_analysis(chromatogram_segment, smoothing_window, peak_resolving_power,
                                                   min_nIsoPair = 0, min_peak_height = 0, min_ratio_IsoPair = 0, max_rpw = 1, min_snr_baseline = 0,
                                                   max_R13C_integrated_peak = Inf, max_percentage_missing_scans  = Inf, mz_target = mzCandidate[j],
                                                   rt_target = rtCandidate[j], mass_accuracy_xic, spectraList, RetentionTime, n_spline)
          if (length(peak_property) > 0) {
            if (length(peak_property) > 24) {
              x_max <- which.max(peak_property[, 7])
              peak_property <- peak_property[x_max[1], ]
            }
            peak_property <- matrix(peak_property, ncol = 1)
            peak_property[R0] <- round(peak_property[R0], 0)
            peak_property[R3] <- round(peak_property[R3], 3)
            peak_property[R5] <- round(peak_property[R5], 5)
          } else {
            peak_property <- rep(0, 24)
          }
          ##
          peak_property_xic <- c(peak_property[3], peak_property[8], peak_property[10], mass_accuracy_xic, peak_property[1:2], peak_property[4:7], peak_property[9], peak_property[11:24])
          peak_property_xic <- data.frame(peak_property_xic)
          colnames(peak_property_xic) <- ""
          rownames(peak_property_xic) <- c("Retention time (min)", "m/z (monoisotopic)", "m/z (13C)", "Mass tolerance (Da)", "Scan (start)",
                                           "Scan (end)", "Peak height", "Peak area", "nIsoPair", "RCS(%)", "Cumulative intensity",
                                           "R13C (%)", "Peak width @ baseline", "Ratio peak width @ 50%", "Separation trays", "Asymmetry factor @ 10%",
                                           "USP tailing factor @ 5%", "Skewness derivative method", "Symmetry pseudo-moments", "Skewness pseudo-moments",
                                           "Gaussianity", "S/N baseline", "S/N xcms method", "S/N RMS", "Sharpness")
          if (exportEIC == TRUE) {
            EIC_figure <- EIC_plotter(chromatogram_segment, peak_property_xic, smoothing_window, peak_resolving_power,
                                      mass_accuracy_xic, spectraList, RetentionTime, mz_target = mzCandidate[j], rt_target = rtCandidate[j], file_name = file_name_hrms[i], legend_EIC)
            ggsave(filename=paste0("/IPA_EIC_", file_name_hrms[i], "_", j, "_", round(mzCandidate[j], 5), "_",round(rtCandidate[j], 2), ".png"),
                   plot = EIC_figure,
                   device = "png",
                   path = output_path_eic,
                   scale = 1,
                   width = 16,
                   height = 8,
                   units = "in",
                   dpi = 100)
          }
          cbind(file_name_hrms[i], round(mzCandidate[j], 5), round(rtCandidate[j], 2), data.frame(matrix(peak_property, nrow = 1)))
        }))
        chrome
      }
      stopCluster(clust)
      ##
    } else if (osType == "Linux") {
      cc_table <- do.call(rbind, lapply(1:length(file_name_hrms), function(i) {
        ##
        outputer <- IPA_MSdeconvoluter(input_path_hrms, file_name_hrms[i])
        spectraList <- outputer[["spectraList"]]
        RetentionTime <- outputer[["retentionTime"]]
        nRT <- length(RetentionTime)
        ##
        chrome <- do.call(rbind, mclapply(1:length(mzCandidate), function(j) {
          ScanNumberApex <- which.min(abs(RetentionTime - rtCandidate[j]))
          t1 <- ScanNumberApex - scan_tol
          if (t1 < 1) {
            t1 <- 1
          }
          t2 <- ScanNumberApex + scan_tol
          if (t2 > nRT) {
            t2 <- nRT
          }
          ##
          chromatogram_segment <- do.call(rbind, lapply(t1:t2, function(t) {
            Spec_ScN_j <- NULL
            Spec <- spectraList[[t]]
            if (length(Spec) > 0) {
              x_mz1 <- which(abs(Spec[, 1] - mzCandidate[j]) <= mass_accuracy_xic)
              L_x_mz1 <- length(x_mz1)
              if (L_x_mz1 > 0) {
                x_mz2 <- which(abs(Spec[, 1] - (massDifferenceIsotopes + mzCandidate[j])) <= mass_accuracy_isotope_pair)
                L_x_mz2 <- length(x_mz2)
                if (L_x_mz2 > 0) {
                  if (L_x_mz1 > 1) {
                    x_min <- which.min(abs(Spec[x_mz1, 1] - mzCandidate[j]))
                    x_mz1 <- x_mz1[x_min[1]]
                  }
                  if (L_x_mz2 > 1) {
                    x_min <- which.min(abs(Spec[x_mz2, 1] - (massDifferenceIsotopes + mzCandidate[j])))
                    x_mz2 <- x_mz2[x_min[1]]
                  }
                  if (Spec[x_mz1, 2] >= Spec[x_mz2, 2]) {
                    Spec_ScN_j <- c(Spec[x_mz1, 1], Spec[x_mz1, 2], t, Spec[x_mz2, 1], Spec[x_mz2, 2])
                  }
                }
              }
            }
            Spec_ScN_j
          }))
          if (length(chromatogram_segment) != 0) {
            chromatogram_segment <- chromatogram_segment[order(chromatogram_segment[, 3]), ]
            chromatogram_segment <- matrix(chromatogram_segment, ncol = 5)
          } else {
            chromatogram_segment <- matrix(c(mzCandidate[j], mzCandidate[j], 0 , 0 , t1, t2, 0, 0, 0, 0), ncol = 5)
          }
          if (chromatogram_segment[1, 3] != t1) {
            chromatogram_segment <- rbind(c(mzCandidate[j], 0, t1, 0, 0), chromatogram_segment)
          }
          if (chromatogram_segment[dim(chromatogram_segment)[1], 3] != t2) {
            chromatogram_segment <- rbind(chromatogram_segment, c(mzCandidate[j], 0, t2, 0, 0))
          }
          peak_property <- chromatography_analysis(chromatogram_segment, smoothing_window, peak_resolving_power,
                                                   min_nIsoPair = 0, min_peak_height = 0, min_ratio_IsoPair = 0, max_rpw = 1, min_snr_baseline = 0,
                                                   max_R13C_integrated_peak = Inf, max_percentage_missing_scans  = Inf, mz_target = mzCandidate[j],
                                                   rt_target = rtCandidate[j], mass_accuracy_xic, spectraList, RetentionTime, n_spline)
          if (length(peak_property) > 0) {
            if (length(peak_property) > 24) {
              x_max <- which.max(peak_property[, 7])
              peak_property <- peak_property[x_max[1], ]
            }
            peak_property <- matrix(peak_property, ncol = 1)
            peak_property[R0] <- round(peak_property[R0], 0)
            peak_property[R3] <- round(peak_property[R3], 3)
            peak_property[R5] <- round(peak_property[R5], 5)
          } else {
            peak_property <- rep(0, 24)
          }
          ##
          peak_property_xic <- c(peak_property[3], peak_property[8], peak_property[10], mass_accuracy_xic, peak_property[1:2], peak_property[4:7], peak_property[9], peak_property[11:24])
          peak_property_xic <- data.frame(peak_property_xic)
          colnames(peak_property_xic) <- ""
          rownames(peak_property_xic) <- c("Retention time (min)", "m/z (monoisotopic)", "m/z (13C)", "Mass tolerance (Da)", "Scan (start)",
                                           "Scan (end)", "Peak height", "Peak area", "nIsoPair", "RCS(%)", "Cumulative intensity",
                                           "R13C (%)", "Peak width @ baseline", "Ratio peak width @ 50%", "Separation trays", "Asymmetry factor @ 10%",
                                           "USP tailing factor @ 5%", "Skewness derivative method", "Symmetry pseudo-moments", "Skewness pseudo-moments",
                                           "Gaussianity", "S/N baseline", "S/N xcms method", "S/N RMS", "Sharpness")
          if (exportEIC == TRUE) {
            EIC_figure <- EIC_plotter(chromatogram_segment, peak_property_xic, smoothing_window, peak_resolving_power,
                                      mass_accuracy_xic, spectraList, RetentionTime, mz_target = mzCandidate[j], rt_target = rtCandidate[j], file_name = file_name_hrms[i], legend_EIC)
            ggsave(filename=paste0("/IPA_EIC_", file_name_hrms[i], "_", j, "_", round(mzCandidate[j], 5), "_",round(rtCandidate[j], 2), ".png"),
                   plot = EIC_figure,
                   device = "png",
                   path = output_path_eic,
                   scale = 1,
                   width = 16,
                   height = 8,
                   units = "in",
                   dpi = 100)
          }
          cbind(file_name_hrms[i], round(mzCandidate[j], 5), round(rtCandidate[j], 2), data.frame(matrix(peak_property, nrow = 1)))
        }, mc.cores = number_processing_threads))
        chrome
      }))
      closeAllConnections()
    }
    ####
    print("Completed the targeted analysis!")
  }
  if (exportTable == TRUE) {
    if (length(cc_table) > 0) {
      cc_table <- data.frame(cc_table)
      names(cc_table) <- c("Name HRMS", "m/z candidate", "RT candidate",
                           "ScanNumberStart","ScanNumberEnd","RetentionTimeApex","PeakHeight","PeakArea",
                           "NumberDetectedScans(nIsoPair)","RCS(%)","m/z MonoIsotopic","CumulatedIntensity",
                           "m/z 13C","Ratio 13C CumulatedIntensity","PeakWidthBaseline","Ratio PeakWidth @ 50%",
                           "SeperationTray","AsymmetryFactor @ 10%","USPTailingFactor @ 5%",
                           "Skewness_DerivativeMethod", "Symmetry PseudoMoments","Skewness PseudoMoments",
                           "Gaussianity", "S/N baseline", "S/N xcms method", "S/N RMS", "Sharpness")
      rownames(cc_table) <- c()
    }
  }
  return(cc_table)
}
