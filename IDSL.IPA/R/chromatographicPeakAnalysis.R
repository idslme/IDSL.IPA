chromatographicPeakAnalysis <- function(spectraScanXIC, aggregatedSpectraList, retentionTime, LretentionTime, massAccuracy, mzTarget,
                                        rtTarget = NULL, scanNumberStart, scanNumberEnd, smoothingWindow, peakResolvingPower, minNIonPair,
                                        minPeakHeight, minRatioIonPair, maxRPW, minSNRbaseline, maxR13CcumulatedIntensity,
                                        maxPercentageMissingScans, nSpline, exportEICparameters  = NULL) {
  ##
  chromatographyCharacteristics <- NULL
  ##
  ##############################################################################
  ##
  rtTargetedCheck <- if (!is.null(rtTarget)) {TRUE} else {FALSE}
  ##
  ##############################################################################
  ## To create chromatogramMatrix
  chromatogramMatrix <- XIC(aggregatedSpectraList, scanNumberStart, scanNumberEnd, mzTarget, massAccuracy)
  ##
  if (!is.null(chromatogramMatrix)) {
    ##
    fillingWindow <- floor(0.05*(scanNumberEnd - scanNumberStart)) + 1
    topScans <- (scanNumberStart - fillingWindow - 1):(scanNumberStart - 1)
    xTopScans <- which(topScans >= 1)
    LxTopScans <- length(xTopScans)
    if (LxTopScans > 0) {
      topChromatogramMatrix <- cbind(topScans[xTopScans], rep(0, LxTopScans), rep(0, LxTopScans))
    } else {
      topChromatogramMatrix <- NULL
    }
    bottomScans <- (scanNumberEnd + 1):(scanNumberEnd + fillingWindow + 1)
    xBottomScans <- which(bottomScans <= LretentionTime)
    LxBottomScans <- length(xBottomScans)
    if (LxBottomScans > 0) {
      bottomChromatogramMatrix <- cbind(bottomScans[xBottomScans], rep(0, LxBottomScans), rep(0, LxBottomScans))
    } else {
      bottomChromatogramMatrix <- NULL
    }
    chromatogramMatrix <- rbind(topChromatogramMatrix, chromatogramMatrix, bottomChromatogramMatrix)
    SZC <- nrow(chromatogramMatrix)
    mz12C <- chromatogramMatrix[, 2]
    chromatogramMatrix <- chromatogramMatrix[, -2]
    ## chromatogramMatrix <- (scan number, smoothed chromatogram, unprocessed m/z chromatogram, raw chromatogram from 13C filter)
    chromatogramMatrix <- cbind(chromatogramMatrix, chromatogramMatrix[, 2], rep(0, SZC))
    ############################################################################
    ## To fill scans with confirmed ion pairs and update spectraScanXIC
    LspectraScanXIC <- nrow(spectraScanXIC)
    ##
    index12C <- which(chromatogramMatrix[, 1] %in% spectraScanXIC[, 3])
    xSpectraScanXIC <- do.call(c, lapply(1:LspectraScanXIC, function(i) {
      if (abs(mz12C[index12C[i]] - spectraScanXIC[i, 1]) <= 1e-10 ) {
        i
      }
    }))
    ##
    LxSpectraScanXIC <- length(xSpectraScanXIC)
    if ((LxSpectraScanXIC >= minNIonPair) & (LxSpectraScanXIC > 0)) {
      chromatogramMatrix[index12C[xSpectraScanXIC], 4] <- spectraScanXIC[xSpectraScanXIC, 2]
      if (LxSpectraScanXIC == 1) {
        spectraScanXIC <- matrix(spectraScanXIC[xSpectraScanXIC, ], ncol = 5)
      } else {
        spectraScanXIC <- spectraScanXIC[xSpectraScanXIC, ]
      }
      ##########################################################################
      ## Smoothing the chromatogram trace over a smoothing window
      chromatogramMatrix <- data.frame(chromatogramMatrix)
      colnames(chromatogramMatrix) <- c("scanNumber", "smoothChromatogram", "rawChromatogram", "ionPairScans")
      loessSZC <- loess(smoothChromatogram ~ scanNumber, data = chromatogramMatrix, span = smoothingWindow/SZC, control = loess.control(surface = "direct"))
      chromatogramMatrix[, 2] <- predict(loessSZC)
      xNeg <- which(chromatogramMatrix[, 2] < 0)
      chromatogramMatrix[xNeg, 2] <- 0
      ##
      if (chromatogramMatrix[1, 1] == 1) {
        chromatogramMatrix[1, 2] <- 0
      }
      if (chromatogramMatrix[SZC, 1] == LretentionTime) {
        chromatogramMatrix[SZC, 2] <- 0
      }
      ##########################################################################
      ## Peak detection module
      Segment <- chromatographicPeakDetector(chromatogramMatrix[, 2])
      if (length(Segment) > 0) {
        baselineVector <- IPA_baselineDeveloper(Segment, chromatogramMatrix[, 2])
        ## To merge the fronting or tailing of the peaks
        if (peakResolvingPower > 0) {
          Segment <- peakFrontingTailingResolver(Segment, chromatogramMatrix[, 2], smoothingWindow, peakResolvingPower)
        }
        ##
        if (rtTargetedCheck) {
          segmentOriginal <- Segment
          ##
          x_rt <- which.min(abs(retentionTime[chromatogramMatrix[, 1]] - rtTarget))
          x_seg <- which(Segment[, 1] <= x_rt & Segment[, 2] >= x_rt)
          Lx_seg <- length(x_seg)
          if (Lx_seg == 0) {
            SeG <- do.call(c, lapply(1:nrow(Segment), function(s) {
              (Segment[s, 1] + Segment[s, 2])/2
            }))
            x_seg <- which.min(abs(SeG - x_rt))[1]
            Lx_seg <- 1
          }
          if (Lx_seg > 0) {
            Segment <- matrix(Segment[x_seg[1], ], ncol = 2)
          } else {
            Segment <- NULL
          }
        }
        ##
        nPeakIsomers <- length(Segment)/2
        if (nPeakIsomers > 0) {
          ######################################################################
          ## Height adjustment
          for (i in 1:nPeakIsomers) {
            maxIntRawChromatogram <- max(chromatogramMatrix[Segment[i, 1]:Segment[i, 2], 3])
            if (maxIntRawChromatogram >= minPeakHeight) {
              maxIntSmoothedChromatogram <- max(chromatogramMatrix[Segment[i, 1]:Segment[i, 2], 2])
              chromatogramMatrix[(Segment[i, 1] + 1):(Segment[i, 2] - 1), 2] <- chromatogramMatrix[(Segment[i, 1] + 1):(Segment[i, 2] - 1), 2]*maxIntRawChromatogram/maxIntSmoothedChromatogram
            } else {
              Segment[i, ] <- c(0, 0)
            }
          }
          Segment <- Segment[(Segment[, 1] != 0), ]
          nPeakIsomers <- length(Segment)/2
          if (nPeakIsomers > 0) {
            ##
            ####################################################################
            ##
            if (nPeakIsomers == 1) {
              Segment <- matrix(Segment, ncol = 2)
            }
            ##
            chromatographyCharacteristics <- do.call(rbind, lapply(1:nPeakIsomers, function(i) {
              cc <- rep(0, 24)
              ##
              X_Seg <- Segment[i, 1]:Segment[i, 2]
              nIsoPair <- length(which(chromatogramMatrix[X_Seg, 4] > baselineVector[X_Seg])) # Number of actual scans inside a chromatographic peak
              if (nIsoPair >= minNIonPair) {
                maxIntRawChromatogram <- max(chromatogramMatrix[X_Seg, 3])
                ##
                X_Int_80 <- which(chromatogramMatrix[X_Seg, 3]/maxIntRawChromatogram >= 0.20)
                L_X_Int_80 <- length(X_Int_80)
                if (L_X_Int_80 > 0) {
                  L_AvailableScans <- X_Int_80[L_X_Int_80] - X_Int_80[1] + 1
                  X_Int_80 <- X_Int_80 + Segment[i, 1] - 1
                  L_missing_scans <- length(which(chromatogramMatrix[X_Int_80, 3] == 0))
                  chrome_gap <- L_missing_scans/L_AvailableScans*100
                  if (chrome_gap <= maxPercentageMissingScans) {
                    ## Ratio of number of detected versus  number of available scans
                    L_Int_80 <- length(which(chromatogramMatrix[X_Int_80, 4] > 0)) # Number of actual scans inside a chromatographic peak
                    RCS <- L_Int_80/L_AvailableScans*100
                    if (RCS >= minRatioIonPair) {
                      Int_Seg <- chromatogramMatrix[X_Seg, 2]
                      BL_Seg <- baselineVector[X_Seg]
                      ## Signal to noise ratio (SNR) using the vertical proportions
                      SNR_bl <- SNRbaseline(Int_Seg, BL_Seg)
                      if (SNR_bl >= minSNRbaseline) {
                        ## Correcting RT and Int boundaries
                        RT_Seg <- retentionTime[chromatogramMatrix[X_Seg, 1]]
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
                          INT <- sum(chromatogramMatrix[X_Seg, 3])
                          MZ <- sum(mz12C[X_Seg]*chromatogramMatrix[X_Seg, 3])/INT
                          ## To calculate the ratio of the [M+1] peak
                          X_M <- which(spectraScanXIC[, 3] >= chromatogramMatrix[Segment[i, 1], 1] & spectraScanXIC[, 3] <= chromatogramMatrix[Segment[i, 2], 1])
                          INT13C <- sum(spectraScanXIC[X_M, 5])
                          MZ13C <- sum(spectraScanXIC[X_M, 4]*spectraScanXIC[X_M, 5])/INT13C
                          Ratio13C_12C <- INT13C/INT*100
                          if (Ratio13C_12C <= maxR13CcumulatedIntensity) {
                            cc[1] <- chromatogramMatrix[Segment[i, 1], 1]
                            cc[2] <- chromatogramMatrix[Segment[i, 2], 1]
                            xMaxIntRawChromatogram <- which(chromatogramMatrix[X_Seg, 3] == maxIntRawChromatogram)
                            xMaxIntRawChromatogram <- xMaxIntRawChromatogram[1] + Segment[i, 1] - 1
                            cc[3] <- retentionTime[chromatogramMatrix[xMaxIntRawChromatogram, 1]]
                            cc[4] <- maxIntRawChromatogram
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
                            ##
                            ####################################################
                            ## To further smooth the peak for more detailed computations
                            if (nSpline > 0) {
                              chromSpline <- spline(RT_Seg, Int_Seg , n = nSpline, method = "fmm",
                                                    xmin = RT_Seg[1], xmax = RT_Seg[N_Seg], ties = mean) # To smooth the curve for derivative calculations
                              rtSpline <- chromSpline[[1]]
                              intSpline <- chromSpline[[2]]
                              baselineApprox <- approx(RT_Seg, BL_Seg, rtSpline, method = "linear", 0, 0, rule = 2, f = 0, ties = mean)[[2]]
                              ## Number of Separation plates
                              cc[14] <- 5.54*((rtSpline[which.max(intSpline)[1]]/W0.5)^2)
                              ## Asymmetry factor
                              cc[15] <- peakAsymmetryFactorCalculator(rtSpline, intSpline)
                              ## USP tailing factor
                              cc[16] <- peakUSPtailingFactorCalculator(rtSpline, intSpline)
                              ## Skewness using a derivative method
                              cc[17] <- peakDerivativeSkewnessCalculator(rtSpline, (intSpline - baselineApprox))
                              ## Symmetry using pseudo-moments
                              outputPseudomomentsSymmetry <- peakPseudomomentsSymmetryCalculator(rtSpline, (intSpline - baselineApprox))
                              cc[18] <- outputPseudomomentsSymmetry[1] # peak symmetry
                              cc[19] <- outputPseudomomentsSymmetry[2] # skewness
                              ## Gaussianity test
                              cc[20] <- peakGaussianityCalculator(rtSpline, intSpline, baselineApprox, gauge = 0.8)
                              ## Signal to noise ratio (SNR) using the xcms approach
                              cc[22] <- SNRxcms(intSpline)
                              ## Signal to noise ratio (SNR) using the root mean square method
                              cc[23] <- SNRrms(intSpline, baselineApprox, gauge = 0.80)
                              ## Peak sharpness
                              cc[24] <- peakSharpnessCalculator(intSpline)
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
              return(cc)
            }))
            ##
            ####################################################################
            ##
            xNon <- which(chromatographyCharacteristics[, 1] != 0)
            LxNonCC <- length(xNon)
            if (LxNonCC == 0) {
              chromatographyCharacteristics <- NULL
            } else if (LxNonCC == 1) {
              chromatographyCharacteristics <- matrix(chromatographyCharacteristics[xNon, ], nrow = 1)
            } else {
              chromatographyCharacteristics <- chromatographyCharacteristics[xNon, ]
              ##
              if (rtTargetedCheck) {
                xMaxIntensity <- which.max(chromatographyCharacteristics[, 4])
                chromatographyCharacteristics <- matrix(chromatographyCharacteristics[xMaxIntensity[1], ], nrow = 1)
              }
            }
            ##
            ####################################################################
            ##
          }
        }
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
    if (!is.null(chromatographyCharacteristics)) {
      ##
      RTCB <- retentionTime[chromatogramMatrix[, 1]]
      ##
      if (rtTargetedCheck) {
        ##
        nPeakIsomers <- length(segmentOriginal)/2
        ## Local minima
        xLocMin <- unique(as.vector(segmentOriginal))
        ##
        for (i in 1:nPeakIsomers) {
          maxIntRawChromatogram <- max(chromatogramMatrix[segmentOriginal[i,1]:segmentOriginal[i, 2], 3])
          if (maxIntRawChromatogram != 0) {
            maxIntSmoothedChromatogram <- max(chromatogramMatrix[segmentOriginal[i,1]:segmentOriginal[i, 2], 2])
            chromatogramMatrix[(segmentOriginal[i, 1] + 1):(segmentOriginal[i, 2] - 1), 2] <- chromatogramMatrix[(segmentOriginal[i, 1] + 1):(segmentOriginal[i, 2] - 1), 2]*maxIntRawChromatogram/maxIntSmoothedChromatogram
          }
        }
        ##
        no.Peaks <- 1
        ##
      } else {
        ## Local minima
        xLocMin <- unique(as.vector(Segment))
        ##
        rtTarget <- chromatographyCharacteristics[, 3]
        ##
        no.Peaks <- LxNonCC
      }
      ##
      R0 <- c(4, 5, 7, 9, 14, 24)
      R3 <- c(3, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23)
      R5 <- c(8, 10)
      ##
      chromatographyCharacteristics[, R0] <- round(chromatographyCharacteristics[, R0], 0)
      chromatographyCharacteristics[, R3] <- round(chromatographyCharacteristics[, R3], 3)
      chromatographyCharacteristics[, R5] <- round(chromatographyCharacteristics[, R5], 5)
      ##
      for (i in 1:no.Peaks) {
        peak_property_xic <- c(chromatographyCharacteristics[i, 3], chromatographyCharacteristics[i, 8],
                               chromatographyCharacteristics[i, 10], chromatographyCharacteristics[i, 1:2],
                               chromatographyCharacteristics[i, 4:7], chromatographyCharacteristics[i, 9], chromatographyCharacteristics[i, 11:24])
        ##
        peak_property_xic <- do.call(c, lapply(1:24, function(j) {paste0(names_peak_property_xic[j], " = ", peak_property_xic[j])}))
        ##
        if (!rtTargetedCheck) {
          x_rt <- which.min(abs(RTCB - rtTarget[i]))
        }
        Int_rt <- max(c(chromatogramMatrix[x_rt, 2], chromatogramMatrix[x_rt, 3]))
        ## Index numbers of scans with ion pairs
        xC_pair <- which(chromatogramMatrix[, 4] > 0)
        ##
        if (exportEICparameters[3] == "UnTargetedWorkflow") {
          MZpeaklist <- chromatographyCharacteristics[i, 8]
          RTpeaklits <- chromatographyCharacteristics[i, 3]
          EICfilename <- paste0(exportEICparameters[1], "/IPA_EIC_", exportEICparameters[2], "_", exportEICparameters[3], "_MZ_", MZpeaklist, "_RT_", RTpeaklits, "_.png")
        } else {
          EICfilename <- paste0(exportEICparameters[1], "/IPA_EIC_", exportEICparameters[2], "_", exportEICparameters[3], "_MZ_", round(mzTarget, 5), "_RT_", round(rtTarget[i], 3), "_.png")
        }
        ##
        fileCreateRCheck <- file.create(file = EICfilename, showWarnings = FALSE)
        if (fileCreateRCheck) {
          ##
          png(EICfilename, width = 20, height = 10, units = "in", res = 100)
          ##
          layout(matrix(c(1, 2), ncol = 2), widths = c(2, 1))
          ##
          plot(RTCB[xC_pair], chromatogramMatrix[xC_pair, 4], type = "h", lwd = 2, col = "#00FF7F",
               xlim = c(RTCB[1], RTCB[SZC]),
               ylim = c(0, max(chromatogramMatrix[, 2])*1.2),
               xlab = "", ylab = "", xaxs = "i") ## Ion-paired chromatogram scan
          abline(h = 0, col = "black")
          lines(RTCB, chromatogramMatrix[, 3], col = "steelblue", lwd = 3 , lty = 6) ## Raw chromatogram
          lines(RTCB, baselineVector, col = "red", lwd = 3 , lty = 4) ## Baseline
          lines(RTCB, chromatogramMatrix[, 2], type = "l", lwd = 3, col = "#6A5ACD", lty = 1) ## Smoothed chromatogram
          points(RTCB[xLocMin], chromatogramMatrix[xLocMin, 2], type = "p", col = "black", cex = 2, pch = 18) ## Local minimum
          legend(x = "topleft", legend = c("Ion-paired chromatogram scan", "Raw chromatogram", "Smoothed chromatogram", "Baseline", "Local minimum"),
                 lty = c(1, 6, 1, 4, NA),
                 col = c("#00FF7F", "steelblue", "#6A5ACD", "red", "black"),
                 lwd = c(3, 3, 3, 3, NA),
                 pch = c(NA, NA, NA, NA, 18),
                 pt.cex = c(NA, NA, NA, NA, 2),
                 cex = 1.5, bty = "n", seg.len = 2, x.intersp = 0.4, y.intersp = 0.9)
          mtext(exportEICparameters[2], side = 3, adj = 0, line = 0.25, cex = 1.25)
          mtext("Retention Time (min)", side = 1, adj = 0.5, line = 2, cex = 1.35)
          mtext("Intensity", side = 2, adj = 0.5, line = 2, cex = 1.35)
          mtext(text = paste0("Seed m/z = ", round(mzTarget, digits = 5), " +/- ", massAccuracy, " Da"), side = 3, adj = 1, line = 0.25, cex = 1.0)
          text(x = rtTarget[i] , y = Int_rt*1.015, cex = 2.25, label = "*", col = "red")
          ##
          plot.new()
          legend(x = "center", legend = peak_property_xic, cex = 1.30, bty = "n", x.intersp = 0.05, y.intersp = 1.22, seg.len = 0)
          ##
          dev.off()
        } else {
          IPA_logRecorder(paste0("WARNING!!! EIC figures can not be created for `", exportEICparameters[2], "` due to character length limit!"))
        }
      }
      ##
    } else {
      ##
      if (rtTargetedCheck) {
        ## To run this block only for the `IPA_targeted` function
        if ((exportEICparameters[3] != "UnTargetedWorkflow")) {
          ##
          EICfilename <- paste0(exportEICparameters[1], "/IPA_EIC_", exportEICparameters[2], "_", exportEICparameters[3], "_MZ_", round(mzTarget, 5), "_RT_", round(rtTarget, 3), "_.png")
          ##
          fileCreateRCheck <- file.create(file = EICfilename, showWarnings = FALSE)
          if (fileCreateRCheck) {
            ##
            chromatographyCharacteristics <- rep(NA, 24)
            chromatographyCharacteristics[3] <- rtTarget
            ##
            peak_property_xic <- c(rtTarget, rep(NA, 23))
            ##
            peak_property_xic <- do.call(c, lapply(1:24, function(j) {paste0(names_peak_property_xic[j], " = ", peak_property_xic[j])}))
            ##
            png(EICfilename, width = 20, height = 10, units = "in", res = 100)
            ##
            layout(matrix(c(1, 2), ncol = 2), widths = c(2, 1))
            ##
            plot(retentionTime, rep(0, LretentionTime), type = "l", lwd = 1, col = "black",
                 ylim = c(0, 1),
                 xlab = "", ylab = "", xaxs = "i")
            mtext(exportEICparameters[2], side = 3, adj = 0, line = 0.25, cex = 1.25)
            mtext("Retention Time (min)", side = 1, adj = 0.5, line = 2, cex = 1.35)
            mtext("Intensity", side = 2, adj = 0.5, line = 2, cex = 1.35)
            mtext(text = paste0("Seed m/z = ", round(mzTarget, digits = 5), " +/- ", massAccuracy, " Da"), side = 3, adj = 1, line = 0.25, cex = 1.0)
            text(x = rtTarget , y = 0.5, cex = 2.25, label = "*", col = "red")
            ##
            plot.new()
            legend(x = "center", legend = peak_property_xic, cex = 1.30, bty = "n", x.intersp = 0.05, y.intersp = 1.22, seg.len = 0)
            ##
            dev.off()
          } else {
            IPA_logRecorder(paste0("WARNING!!! EIC figures can not be created for `", exportEICparameters[2], "` due to character length limit!"))
          }
        }
      }
    }
  }
  ##
  ##############################################################################
  ##
  return(chromatographyCharacteristics)
}
