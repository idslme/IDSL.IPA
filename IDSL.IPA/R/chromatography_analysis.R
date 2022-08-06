chromatography_analysis <- function(spec_scan_xic, smoothing_window, peak_resolving_power,
                                    min_nIsoPair, min_peak_height, min_ratio_IsoPair, max_rpw, min_snr_baseline,
                                    max_R13C_integrated_peak, max_percentage_missing_scans, mz_target, rt_target = 0,
                                    mass_accuracy_xic, spectraList, RetentionTime, n_spline) {
  chromatography_characteristics <- NULL
  ## To create chromatogram builder
  # chromatogram builder <- (scan number, smoothed chromatogram, unprocessed m/z chromatogram, raw chromatogram from 13C filter)
  spec_scan_xic <- spec_scan_xic[order(spec_scan_xic[, 3]), ]             # To sort the array based on the scan number
  spec_scan_xic <- matrix(spec_scan_xic, ncol = 5)
  n_RT <- length(RetentionTime)     # n_RT is the maximum number of scan number
  ScanNumberStart <- spec_scan_xic[1, 3]
  ScanNumberEnd <- spec_scan_xic[nrow(spec_scan_xic), 3]
  filling_window <- floor(0.05*(ScanNumberEnd - ScanNumberStart)) + 1 # Supposedly minimum space between two peaks
  chrom_builder_temp <- XIC(spectraList[ScanNumberStart:ScanNumberEnd], scan_number_start = ScanNumberStart, mz_target, mass_accuracy_xic)
  Top_ScN <- (ScanNumberStart - filling_window - 1):(ScanNumberStart - 1)
  x_Top <- which(Top_ScN > 0)
  L_Top <- length(x_Top)
  if (L_Top > 0) {
    Top_chrom_builder <- cbind(Top_ScN[x_Top], rep(0, L_Top), rep(0, L_Top))
  } else {
    Top_chrom_builder <- c()
  }
  Bottom_ScN <- (ScanNumberEnd + 1):(ScanNumberEnd + filling_window + 1)
  x_Bottom <- which(Bottom_ScN <= n_RT)
  L_Bottom <- length(x_Bottom)
  if (L_Bottom > 0) {
    Bottom_chrom_builder <- cbind(Bottom_ScN[x_Bottom], rep(0, L_Bottom), rep(0, L_Bottom))
  } else {
    Bottom_chrom_builder <- c()
  }
  chrom_builder_temp <- rbind(Top_chrom_builder, chrom_builder_temp, Bottom_chrom_builder)
  SZC <- nrow(chrom_builder_temp)
  mz12C <- chrom_builder_temp[, 2]
  Chrom_Builder <- chrom_builder_temp[, -2]
  Chrom_Builder <- cbind(Chrom_Builder, Chrom_Builder[, 2], rep(0, SZC))
  # To fill scans with the detected 12C/13C isotopologue pairs
  t_T <- table(spec_scan_xic[, 3])
  x_Tn <- which(t_T > 1)
  L_x_Tn <- length(x_Tn)
  if (L_x_Tn == 0) {
    index1 <- Chrom_Builder[, 1] %in% spec_scan_xic[, 3]
    Chrom_Builder[index1, 4] <- spec_scan_xic[, 2]
  } else {
    x_T1 <- which(t_T == 1)
    T1 <- as.numeric(names(x_T1))
    index_Chrom_T1 <- Chrom_Builder[, 1] %in% T1
    index_xic_T1 <- spec_scan_xic[, 3] %in% T1
    Chrom_Builder[index_Chrom_T1, 4] <- spec_scan_xic[index_xic_T1, 2]
    Tn <- as.numeric(names(x_Tn))
    for (t in 1:L_x_Tn) {
      xT_chrom <- which(Tn[t] == Chrom_Builder[, 1])
      index1 <- which(Tn[t] == spec_scan_xic[, 3])
      Chrom_Builder[xT_chrom, 4] <- mean(spec_scan_xic[index1, 2])
    }
  }
  ## Smoothing the chromatogram trace over a smoothing window
  Chrom_Builder <- data.frame(Chrom_Builder)
  colnames(Chrom_Builder) <- c("scan_number", "smooth_chrom", "raw_chrom", "C_pair")
  loess_SZC <- loess(smooth_chrom ~ scan_number, data = Chrom_Builder, span = smoothing_window/SZC, control = loess.control(surface = "direct"))
  Chrom_Builder[, 2] <- predict(loess_SZC)
  x_neg <- which(Chrom_Builder[, 2] < 0)
  Chrom_Builder[x_neg, 2] <- 0
  ##
  if (Chrom_Builder[1, 1] == 1) {
    Chrom_Builder[1, 2] <- 0
  }
  if (Chrom_Builder[SZC, 1] == n_RT) {
    Chrom_Builder[SZC, 2] <- 0
  }
  ## Peak detection module
  Segment <- peak_detection(Chrom_Builder[, 2])
  if (length(Segment) > 0) {
    BL_Base <- baseline_developer(Segment, Chrom_Builder[, 2])
    ## To merge the fronting or tailing of the peaks
    if (peak_resolving_power > 0) {
      Segment <- fronting_tailing_resolver(Segment, Chrom_Builder[, 2], smoothing_window, peak_resolving_power)
    }
    if (rt_target != 0) {
      x_rt <- which.min(abs(RetentionTime[Chrom_Builder[, 1]] - rt_target))
      x_seg <- which(Segment[, 1] <= x_rt & Segment[, 2] >= x_rt)
      if (length(x_seg) == 0) {
        SeG <- sapply(1:nrow(Segment), function(s) {
          (Segment[s, 1] + Segment[s, 2])/2
        })
        x_seg <- which.min(abs(SeG - x_rt))
      }
      if (length(x_seg) > 0) {
        Segment <- matrix(Segment[x_seg[1], ] , ncol = 2)
      } else {
        Segment <- c()
      }
    }
    if (length(Segment) > 0) {
      ## Height adjustment
      for (i in 1:nrow(Segment)) {
        C_max_RAW <- max(Chrom_Builder[Segment[i, 1]:Segment[i, 2], 3])
        if (C_max_RAW >= min_peak_height) {
          C_max <- max(Chrom_Builder[Segment[i, 1]:Segment[i, 2], 2])
          Chrom_Builder[(Segment[i, 1] + 1):(Segment[i, 2] - 1), 2] <- Chrom_Builder[(Segment[i, 1] + 1):(Segment[i, 2] - 1), 2]*C_max_RAW/C_max
        } else {
          Segment[i, ] <- c(0, 0)
        }
      }
      Segment <- Segment[(Segment[, 1] != 0), ]
      ##
      if (length(Segment) > 0) {
        Segment <- matrix(Segment, ncol = 2)
        No_isomers <- nrow(Segment)
        chromatography_characteristics <- do.call(rbind, lapply(1:No_isomers, function(i) {
          cc <- rep(0, 24)
          X_Seg <- Segment[i, 1]:Segment[i, 2]
          nIsoPair <- length(which(Chrom_Builder[X_Seg, 4] > BL_Base[X_Seg])) # Number of actual scans inside a chromatographic peak
          if (nIsoPair >= min_nIsoPair) {
            C_max_RAW <- max(Chrom_Builder[X_Seg, 3])
            ##
            X_Int_80 <- which(Chrom_Builder[X_Seg, 3]/C_max_RAW >= 0.20)
            L_X_Int_80 <- length(X_Int_80)
            if (L_X_Int_80 > 0) {
              L_AvailableScans <- X_Int_80[L_X_Int_80] - X_Int_80[1] + 1
              X_Int_80 <- X_Int_80 + Segment[i, 1] - 1
              L_missing_scans <- length(which(Chrom_Builder[X_Int_80, 3] == 0))
              chrome_gap <- L_missing_scans/L_AvailableScans*100
              if (chrome_gap <= max_percentage_missing_scans) {
                ## Ratio of number of detected versus  number of available scans
                L_Int_80 <- length(which(Chrom_Builder[X_Int_80, 4] > 0)) # Number of actual scans inside a chromatographic peak
                RCS <- L_Int_80/L_AvailableScans*100
                if (RCS >= min_ratio_IsoPair) {
                  Int_Seg <- Chrom_Builder[X_Seg, 2]
                  BL_Seg <- BL_Base[X_Seg]
                  ## Signal to noise ratio (SNR) using the vertical proportions
                  SNR_bl <- snr_signal2baseline(Int_Seg, BL_Seg)
                  if (SNR_bl >= min_snr_baseline) {
                    ## Correcting RT and Int boundaries
                    RT_Seg <- RetentionTime[Chrom_Builder[X_Seg, 1]]
                    N_Seg <- length(RT_Seg)
                    if (Int_Seg[1] != 0) {
                      Int_Seg <- c(0, Int_Seg)
                      RT_Seg <- c((RT_Seg[1] - (RT_Seg[2] - RT_Seg[1])), RT_Seg)
                      BL_Seg <- c(0, BL_Seg)
                      N_Seg <- N_Seg + 1
                    }
                    if (Int_Seg[N_Seg] != 0) {
                      Int_Seg <- c(Int_Seg, 0)
                      RT_Seg <- c(RT_Seg, (RT_Seg[N_Seg] + (RT_Seg[2] - RT_Seg[1])))
                      BL_Seg <- c(BL_Seg, 0)
                      N_Seg <- N_Seg + 1
                    }
                    ## Ratio of peak width at half-height to the peak width at baseline
                    W0.5 <- peak_width(RT_Seg, Int_Seg, 0.50)
                    PWBL <- (RT_Seg[N_Seg] - RT_Seg[1])
                    RPW0.5 <- W0.5/PWBL   # Peak width at half height to peak width at baseline (RPW)
                    if (RPW0.5 <= max_rpw) {
                      ## m/z for the integrated chromatographic peak
                      INT <- sum(Chrom_Builder[X_Seg, 3])
                      MZ <- sum(mz12C[X_Seg]*Chrom_Builder[X_Seg, 3])/INT
                      # To calculate the ratio of the [M+1] peak
                      X_M <- which(spec_scan_xic[, 3] >= Chrom_Builder[Segment[i, 1], 1] & spec_scan_xic[, 3] <= Chrom_Builder[Segment[i, 2], 1])
                      INT13C <- sum(spec_scan_xic[X_M, 5])
                      MZ13C <- sum(spec_scan_xic[X_M, 4]*spec_scan_xic[X_M, 5])/INT13C
                      Ratio13C_12C <- INT13C/INT*100
                      if (Ratio13C_12C <= max_R13C_integrated_peak) {
                        cc[1] <- Chrom_Builder[Segment[i, 1], 1]
                        cc[2] <- Chrom_Builder[Segment[i, 2], 1]
                        x_C_max_RAW <- which(Chrom_Builder[X_Seg, 3] == C_max_RAW)
                        x_C_max_RAW <- x_C_max_RAW[1] + Segment[i, 1] - 1
                        cc[3] <- RetentionTime[Chrom_Builder[x_C_max_RAW, 1]]
                        cc[4] <- C_max_RAW
                        cc[5] <- peak_area(RT_Seg, Int_Seg) ## Peak area calculation
                        cc[6] <- nIsoPair
                        cc[7] <- RCS
                        cc[8] <- MZ
                        cc[9] <- INT
                        cc[10] <- MZ13C
                        cc[11] <- Ratio13C_12C
                        cc[12] <- PWBL
                        cc[13] <- RPW0.5
                        cc[21] <- SNR_bl
                        ## further smoothing peaks for more detailed computations
                        if (n_spline > 0) {
                          W <- spline(RT_Seg, Int_Seg , n = n_spline, method = "fmm",
                                      xmin = RT_Seg[1], xmax = RT_Seg[N_Seg], ties = mean) # To smooth the curve for derivative calculations
                          RT_spline <- W[[1]]
                          Int_spline <- W[[2]]
                          BL_approx <- approx(RT_Seg, BL_Seg, RT_spline, method = "linear", 0, 0, rule = 2, f = 0, ties = mean)[[2]]
                          ## Number of Separation plates
                          cc[14] <- 5.54*((RT_spline[which.max(Int_spline)[1]]/W0.5)^2)
                          ## Asymmetry factor
                          cc[15] <- asymmetry_factor(RT_spline, Int_spline)
                          ## USP tailing factor
                          cc[16] <- usp_tailing_factor(RT_spline, Int_spline)
                          ## Skewness using a derivative method
                          cc[17] <- derivative_skewness(RT_spline, (Int_spline - BL_approx))
                          ## Symmetry using pseudo-moments
                          output_CS <- pseudomoments_symmetry(RT_spline, (Int_spline - BL_approx))
                          cc[18] <- output_CS[1] # peak symmetry
                          cc[19] <- output_CS[2] # skewness
                          ## Gaussianity test
                          cc[20] <- gaussianity_measurement(RT_spline, Int_spline, BL_approx, gauge = 0.8)
                          ## Signal to noise ratio (SNR) using the xcms approach
                          cc[22] <- snr_xcms(Int_spline)
                          ## Signal to noise ratio (SNR) using the root mean square method
                          cc[23] <- snr_rms(Int_spline, BL_approx, 0.80)
                          ## Peak sharpness
                          cc[24] <- peak_sharpness(Int_spline)
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          cc
        }))
        x_0 <- which(chromatography_characteristics[, 1] != 0)
        chromatography_characteristics <- chromatography_characteristics[x_0, ]
      }
    }
  }
  return(chromatography_characteristics)
}
