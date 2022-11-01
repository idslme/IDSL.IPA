chromatographyPeakAnalysis <- function(spectraScan_xic, smoothingWindow, peakResolvingPower,
                                       minNIonPair, minPeakHeight, minRatioIonPair, maxRPW, minSNRbaseline,
                                       maxR13CcumulatedIntensity, maxPercentageMissingScans, mzTarget, rtTarget = NULL,
                                       massAccuracyXIC, spectraList, RetentionTime, nSpline, exportEICparameters  = NULL) {
  ##
  ##############################################################################
  ##
  chromatography_characteristics <- NULL
  ##
  rtTargetedCheck <- if (!is.null(rtTarget)) {TRUE} else {FALSE}
  ##
  ## To create chromatogram builder
  # chromatogram builder <- (scan number, smoothed chromatogram, unprocessed m/z chromatogram, raw chromatogram from 13C filter)
  spectraScan_xic <- spectraScan_xic[order(spectraScan_xic[, 3]), ]             # To sort the array based on the scan number
  spectraScan_xic <- matrix(spectraScan_xic, ncol = 5)
  n_RT <- length(RetentionTime)     # n_RT is the maximum number of scan number
  ScanNumberStart <- spectraScan_xic[1, 3]
  ScanNumberEnd <- spectraScan_xic[nrow(spectraScan_xic), 3]
  filling_window <- floor(0.05*(ScanNumberEnd - ScanNumberStart)) + 1 # Minimum space between two peaks
  chrom_builder_temp <- XIC(spectraList[ScanNumberStart:ScanNumberEnd], scanNumberStart = ScanNumberStart, mzTarget, massAccuracyXIC)
  Top_ScN <- (ScanNumberStart - filling_window - 1):(ScanNumberStart - 1)
  x_Top <- which(Top_ScN > 0)
  L_Top <- length(x_Top)
  if (L_Top > 0) {
    Top_chrom_builder <- cbind(Top_ScN[x_Top], rep(0, L_Top), rep(0, L_Top))
  } else {
    Top_chrom_builder <- NULL
  }
  Bottom_ScN <- (ScanNumberEnd + 1):(ScanNumberEnd + filling_window + 1)
  x_Bottom <- which(Bottom_ScN <= n_RT)
  L_Bottom <- length(x_Bottom)
  if (L_Bottom > 0) {
    Bottom_chrom_builder <- cbind(Bottom_ScN[x_Bottom], rep(0, L_Bottom), rep(0, L_Bottom))
  } else {
    Bottom_chrom_builder <- NULL
  }
  chrom_builder_temp <- rbind(Top_chrom_builder, chrom_builder_temp, Bottom_chrom_builder)
  SZC <- nrow(chrom_builder_temp)
  mz12C <- chrom_builder_temp[, 2]
  Chrom_Builder <- chrom_builder_temp[, -2]
  Chrom_Builder <- cbind(Chrom_Builder, Chrom_Builder[, 2], rep(0, SZC))
  # To fill scans with confirmed 12C/13C isotopologue pairs
  t_T <- table(spectraScan_xic[, 3])
  x_Tn <- which(t_T > 1)
  if (length(x_Tn) == 0) {
    index1 <- Chrom_Builder[, 1] %in% spectraScan_xic[, 3]
    Chrom_Builder[index1, 4] <- spectraScan_xic[, 2]
  } else {
    x_T1 <- which(t_T == 1)
    T1 <- as.numeric(names(x_T1))
    index_Chrom_T1 <- Chrom_Builder[, 1] %in% T1
    indexXIC_T1 <- spectraScan_xic[, 3] %in% T1
    Chrom_Builder[index_Chrom_T1, 4] <- spectraScan_xic[indexXIC_T1, 2]
    Tn <- as.numeric(names(x_Tn))
    for (t in Tn) {
      xT_chrom <- which(t == Chrom_Builder[, 1])
      index1 <- which(t == spectraScan_xic[, 3])
      Chrom_Builder[xT_chrom, 4] <- mean(spectraScan_xic[index1, 2])
    }
  }
  ## Smoothing the chromatogram trace over a smoothing window
  Chrom_Builder <- data.frame(Chrom_Builder)
  colnames(Chrom_Builder) <- c("scan_number", "smooth_chrom", "raw_chrom", "C_pair")
  loess_SZC <- loess(smooth_chrom ~ scan_number, data = Chrom_Builder, span = smoothingWindow/SZC, control = loess.control(surface = "direct"))
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
  Segment <- chromatographicPeakDetector(Chrom_Builder[, 2])
  if (length(Segment) > 0) {
    BL_Base <- IPA_baselineDeveloper(Segment, Chrom_Builder[, 2])
    ## To merge the fronting or tailing of the peaks
    if (peakResolvingPower > 0) {
      Segment <- peakFrontingTailingResolver(Segment, Chrom_Builder[, 2], smoothingWindow, peakResolvingPower)
    }
    if (rtTargetedCheck) {
      segmentOriginal <- Segment
      ##
      x_rt <- which.min(abs(RetentionTime[Chrom_Builder[, 1]] - rtTarget))
      x_seg <- which(Segment[, 1] <= x_rt & Segment[, 2] >= x_rt)
      if (length(x_seg) == 0) {
        SeG <- do.call(c, lapply(1:nrow(Segment), function(s) {
          (Segment[s, 1] + Segment[s, 2])/2
        }))
        x_seg <- which.min(abs(SeG - x_rt))[1]
      }
      if (length(x_seg) > 0) {
        Segment <- matrix(Segment[x_seg[1], ] , ncol = 2)
      } else {
        Segment <- NULL
      }
    }
    if (length(Segment) > 0) {
      ## Height adjustment
      for (i in 1:nrow(Segment)) {
        C_max_RAW <- max(Chrom_Builder[Segment[i, 1]:Segment[i, 2], 3])
        if (C_max_RAW >= minPeakHeight) {
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
          if (nIsoPair >= minNIonPair) {
            C_max_RAW <- max(Chrom_Builder[X_Seg, 3])
            ##
            X_Int_80 <- which(Chrom_Builder[X_Seg, 3]/C_max_RAW >= 0.20)
            L_X_Int_80 <- length(X_Int_80)
            if (L_X_Int_80 > 0) {
              L_AvailableScans <- X_Int_80[L_X_Int_80] - X_Int_80[1] + 1
              X_Int_80 <- X_Int_80 + Segment[i, 1] - 1
              L_missing_scans <- length(which(Chrom_Builder[X_Int_80, 3] == 0))
              chrome_gap <- L_missing_scans/L_AvailableScans*100
              if (chrome_gap <= maxPercentageMissingScans) {
                ## Ratio of number of detected versus  number of available scans
                L_Int_80 <- length(which(Chrom_Builder[X_Int_80, 4] > 0)) # Number of actual scans inside a chromatographic peak
                RCS <- L_Int_80/L_AvailableScans*100
                if (RCS >= minRatioIonPair) {
                  Int_Seg <- Chrom_Builder[X_Seg, 2]
                  BL_Seg <- BL_Base[X_Seg]
                  ## Signal to noise ratio (SNR) using the vertical proportions
                  SNR_bl <- SNRbaseline(Int_Seg, BL_Seg)
                  if (SNR_bl >= minSNRbaseline) {
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
                    W0.5 <- peakWidthCalculator(RT_Seg, Int_Seg, 0.50)
                    PWBL <- (RT_Seg[N_Seg] - RT_Seg[1])
                    RPW0.5 <- W0.5/PWBL   # Peak width at half height to peak width at baseline (RPW)
                    if (RPW0.5 <= maxRPW) {
                      ## m/z for the integrated chromatographic peak
                      INT <- sum(Chrom_Builder[X_Seg, 3])
                      MZ <- sum(mz12C[X_Seg]*Chrom_Builder[X_Seg, 3])/INT
                      # To calculate the ratio of the [M+1] peak
                      X_M <- which(spectraScan_xic[, 3] >= Chrom_Builder[Segment[i, 1], 1] & spectraScan_xic[, 3] <= Chrom_Builder[Segment[i, 2], 1])
                      INT13C <- sum(spectraScan_xic[X_M, 5])
                      MZ13C <- sum(spectraScan_xic[X_M, 4]*spectraScan_xic[X_M, 5])/INT13C
                      Ratio13C_12C <- INT13C/INT*100
                      if (Ratio13C_12C <= maxR13CcumulatedIntensity) {
                        cc[1] <- Chrom_Builder[Segment[i, 1], 1]
                        cc[2] <- Chrom_Builder[Segment[i, 2], 1]
                        x_C_max_RAW <- which(Chrom_Builder[X_Seg, 3] == C_max_RAW)
                        x_C_max_RAW <- x_C_max_RAW[1] + Segment[i, 1] - 1
                        cc[3] <- RetentionTime[Chrom_Builder[x_C_max_RAW, 1]]
                        cc[4] <- C_max_RAW
                        cc[5] <- peakAreaCalculator(RT_Seg, Int_Seg) ## Peak area calculation
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
                        if (nSpline > 0) {
                          W <- spline(RT_Seg, Int_Seg , n = nSpline, method = "fmm",
                                      xmin = RT_Seg[1], xmax = RT_Seg[N_Seg], ties = mean) # To smooth the curve for derivative calculations
                          RT_spline <- W[[1]]
                          Int_spline <- W[[2]]
                          BL_approx <- approx(RT_Seg, BL_Seg, RT_spline, method = "linear", 0, 0, rule = 2, f = 0, ties = mean)[[2]]
                          ## Number of Separation plates
                          cc[14] <- 5.54*((RT_spline[which.max(Int_spline)[1]]/W0.5)^2)
                          ## Asymmetry factor
                          cc[15] <- asymmetryFactor(RT_spline, Int_spline)
                          ## USP tailing factor
                          cc[16] <- USPtailingFactor(RT_spline, Int_spline)
                          ## Skewness using a derivative method
                          cc[17] <- derivativeSkewnessCalculator(RT_spline, (Int_spline - BL_approx))
                          ## Symmetry using pseudo-moments
                          output_CS <- pseudomomentsSymmetry(RT_spline, (Int_spline - BL_approx))
                          cc[18] <- output_CS[1] # peak symmetry
                          cc[19] <- output_CS[2] # skewness
                          ## Gaussianity test
                          cc[20] <- gaussianityMeasurement(RT_spline, Int_spline, BL_approx, gauge = 0.8)
                          ## Signal to noise ratio (SNR) using the xcms approach
                          cc[22] <- SNRxcms(Int_spline)
                          ## Signal to noise ratio (SNR) using the root mean square method
                          cc[23] <- SNRrms(Int_spline, BL_approx, gauge = 0.80)
                          ## Peak sharpness
                          cc[24] <- peakSharpness(Int_spline)
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
        chromatography_characteristics <- matrix(chromatography_characteristics[x_0, ], ncol = 24)
      }
    }
  }
  ##
  ##############################################################################
  ## EIC plotting
  if (!is.null(exportEICparameters)) {
    ##
    names_peak_property_xic <- c("Retention time (min)", "m/z (12C)", "m/z (13C)", "Scan (start)", "Scan (end)", "Peak height",
                                 "Peak area", "nIsoPair", "RCS(%)", "Cumulative intensity", "R13C(%)", "Peak width @ baseline",
                                 "Ratio peak width @ 50%", "Separation trays", "Asymmetry factor @ 10%", "USP tailing factor @ 5%",
                                 "Skewness derivative method", "Symmetry pseudo-moments", "Skewness pseudo-moments",
                                 "Gaussianity", "S/N baseline", "S/N xcms method", "S/N RMS", "Sharpness")
    ##
    RTCB <- RetentionTime[Chrom_Builder[, 1]]
    ##
    if (length(chromatography_characteristics) > 0) {
      if (rtTargetedCheck) {
        ##
        no.Isomer <- nrow(segmentOriginal)
        ## Local minima
        xLocMin <- unique(as.vector(segmentOriginal))
        ##
        for (i in 1:no.Isomer) {
          C_max_RAW <- max(Chrom_Builder[segmentOriginal[i,1]:segmentOriginal[i, 2], 3])
          if (C_max_RAW != 0) {
            C_max <- max(Chrom_Builder[segmentOriginal[i,1]:segmentOriginal[i, 2], 2])
            Chrom_Builder[(segmentOriginal[i, 1] + 1):(segmentOriginal[i, 2] - 1), 2] <- Chrom_Builder[(segmentOriginal[i, 1] + 1):(segmentOriginal[i, 2] - 1), 2]*C_max_RAW/C_max
          }
        }
        ##
        no.Peaks <- 1
        if (length(chromatography_characteristics) > 24) {
          x_max <- which.max(chromatography_characteristics[, 7])
          chromatography_characteristics <- matrix(chromatography_characteristics[x_max[1], ], nrow = 1)
        }
      } else {
        ## Local minima
        xLocMin <- unique(as.vector(Segment))
        ##
        rtTarget <- chromatography_characteristics[, 3]
        ##
        no.Peaks <- dim(chromatography_characteristics)[1]
      }
      ##
      R0 <- c(4, 5, 7, 9, 14, 24)
      R3 <- c(3, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23)
      R5 <- c(8, 10)
      ##
      chromatography_characteristics[, R0] <- round(chromatography_characteristics[, R0], 0)
      chromatography_characteristics[, R3] <- round(chromatography_characteristics[, R3], 3)
      chromatography_characteristics[, R5] <- round(chromatography_characteristics[, R5], 5)
      ##
      for (i in 1:no.Peaks) {
        peak_property_xic <- c(chromatography_characteristics[i, 3], chromatography_characteristics[i, 8],
                               chromatography_characteristics[i, 10], chromatography_characteristics[i, 1:2],
                               chromatography_characteristics[i, 4:7], chromatography_characteristics[i, 9], chromatography_characteristics[i, 11:24])
        ##
        peak_property_xic <- do.call(c, lapply(1:24, function(i) {paste0(names_peak_property_xic[i], " = ", peak_property_xic[i])}))
        ##
        if (!rtTargetedCheck) {
          x_rt <- which.min(abs(RTCB - rtTarget[i]))
        }
        Int_rt <- max(c(Chrom_Builder[x_rt, 2], Chrom_Builder[x_rt, 3]))
        ## Index numbers of scans with ion pairs
        xC_pair <- which(Chrom_Builder[, 4] > 0)
        ##
        if (exportEICparameters[3] == "UnTargetedWorkflow") {
          .GlobalEnv$counterEIC <- .GlobalEnv$counterEIC + 1
          EICfilename <- paste0(exportEICparameters[1], "/IPA_EIC_", exportEICparameters[2], "_", exportEICparameters[3], "_", .GlobalEnv$counterEIC, "_MZ_", round(mzTarget, 5), "_RT_", round(rtTarget[i], 3), "_.png")
        } else {
          EICfilename <- paste0(exportEICparameters[1], "/IPA_EIC_", exportEICparameters[2], "_", exportEICparameters[3], "_MZ_", round(mzTarget, 5), "_RT_", round(rtTarget[i], 3), "_.png")
        }
        ##
        png(EICfilename, width = 20, height = 10, units = "in", res = 100)
        ##
        layout(matrix(c(1, 2), ncol = 2), widths = c(2, 1))
        ##
        plot(RTCB[xC_pair], Chrom_Builder[xC_pair, 4], type = "h", lwd = 2, col = "#00FF7F",
             xlim = c(RTCB[1], RTCB[SZC]),
             ylim = c(0, max(Chrom_Builder[, 2])*1.2),
             xlab = "", ylab = "", xaxs = "i") ## Ion-paired chromatogram scan
        abline(h = 0, col = "black")
        lines(RTCB, Chrom_Builder[, 3], col = "steelblue", lwd = 3 , lty = 6) ## Raw chromatogram
        lines(RTCB, BL_Base, col = "red", lwd = 3 , lty = 4) ## Baseline
        lines(RTCB, Chrom_Builder[, 2], type = "l", lwd = 3, col = "#6A5ACD", lty = 1) ## LOESS smoothed chromatogram
        points(RTCB[xLocMin], Chrom_Builder[xLocMin, 2], type = "p", col = "black", cex = 2, pch = 18) ## Local minimum
        legend(x = "topleft", legend = c("Ion-paired chromatogram scan", "Raw chromatogram", "Smoothed chromatogram", "Baseline", "Local minimum"),
               lty = c(1, 6, 1, 4, NA),
               col = c("#00FF7F", "steelblue", "#6A5ACD", "red", "black"),
               lwd = c(3, 3, 3, 3, 0),
               pch = c(NA, NA, NA, NA, 18),
               pt.cex = c(NA, NA, NA, NA, 2),
               cex = 1.5, bty = "n", seg.len = 2, x.intersp = 0.4, y.intersp = 0.9)
        mtext(exportEICparameters[2], side = 3, adj = 0, line = 0.25, cex = 1.4)
        mtext("Retention Time (min)", side = 1, adj = 0.5, line = 2, cex = 1.35)
        mtext("Intensity", side = 2, adj = 0.5, line = 2, cex = 1.35)
        mtext(text = paste0("Seed m/z = ", round(mzTarget, digits = 5), " +/- ", massAccuracyXIC, " Da"), side = 3, adj = 1, line = 0.25, cex = 1.0)
        text(x = rtTarget[i] , y = Int_rt*1.015, cex = 2.25, label = "*", col = "red")
        ##
        plot.new()
        legend(x = "center", legend = peak_property_xic,
               cex = 1.30, bty = "n", x.intersp = 0.05, y.intersp = 1.22, seg.len = 0)
        ##
        dev.off()
        ##
      }
      ##
    } else {
      ##
      if (rtTargetedCheck) {
        ## To run this block only for the `IPA_TargetedAnalysis` function
        if ((exportEICparameters[3] != "UnTargetedWorkflow")) {
          chromatography_characteristics <- rep(NA, 24)
          chromatography_characteristics[3] <- rtTarget
          ##
          peak_property_xic <- c(rtTarget, rep(NA, 23))
          ##
          peak_property_xic <- do.call(c, lapply(1:24, function(i) {paste0(names_peak_property_xic[i], " = ", peak_property_xic[i])}))
          ##
          EICfilename <- paste0(exportEICparameters[1], "/IPA_EIC_", exportEICparameters[2], "_", exportEICparameters[3], "_MZ_", round(mzTarget, 5), "_RT_", round(rtTarget, 3), "_.png")
          ##
          png(EICfilename, width = 20, height = 10, units = "in", res = 100)
          ##
          layout(matrix(c(1, 2), ncol = 2), widths = c(2, 1))
          ##
          plot(RTCB[1], Chrom_Builder[1, 4], type = "h", lwd = 2, col = "#00FF7F",
               xlim = c(RTCB[1], RTCB[SZC]),
               ylim = c(0, max(Chrom_Builder[, 2])*1.2),
               xlab = "", ylab = "", xaxs = "i") ## Ion-paired chromatogram scan
          abline(h = 0, col = "black")
          mtext(text = paste0("Seed m/z = ", round(mzTarget, digits = 5), " +/- ", massAccuracyXIC, " Da"), side = 3, adj = 1, line = 0.25, cex = 1.0)
          mtext(exportEICparameters[2], side = 3, adj = 0, line = 0.25, cex = 1.4)
          ##
          plot.new()
          legend(x = "center", legend = peak_property_xic,
                 cex = 1.30, bty = "n", x.intersp = 0.05, y.intersp = 1.22, seg.len = 0)
          ##
          dev.off()
          ##
        }
      }
    }
  }
  ##
  ##############################################################################
  ##
  return(chromatography_characteristics)
}
