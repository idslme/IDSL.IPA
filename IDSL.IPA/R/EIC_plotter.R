EIC_plotter <- function (spec_scan_xic, peak_property_xic, smoothing_window, peak_resolving_power, mass_accuracy_xic, spectraList, RetentionTime, mz_target, rt_target, file_name='HRMS', legend_EIC) {
  spec_scan_xic <- spec_scan_xic[order(spec_scan_xic[, 3]), ] # To sort the array based on the scan number
  spec_scan_xic <- matrix(spec_scan_xic, ncol = 5)
  n_RT <- length(RetentionTime) # n_RT is the maximum number of scan number
  ScanNumberStart <- spec_scan_xic[1, 3]
  ScanNumberEnd <- spec_scan_xic[nrow(spec_scan_xic), 3]
  filling_window <- floor(0.05*(ScanNumberEnd - ScanNumberStart)) + 1 # Supposedly minimum space between two peaks
  chrom_builder_temp <- XIC(spectraList[ScanNumberStart:ScanNumberEnd], scan_number_start = ScanNumberStart, mz_target, mass_accuracy_xic)
  Top_ScN <- (ScanNumberStart - filling_window - 1) : (ScanNumberStart - 1)
  x_Top <- which(Top_ScN > 0)
  L_Top <- length(x_Top)
  if (L_Top > 0) {
    Top_chrom_buider <- cbind(Top_ScN[x_Top], rep(0, L_Top), rep(0, L_Top))
  } else {
    Top_chrom_buider <- c()
  }
  Bottom_ScN <- (ScanNumberEnd + 1) : (ScanNumberEnd + filling_window + 1)
  x_Bottom <- which(Bottom_ScN <= n_RT)
  L_Bottom <- length(x_Bottom)
  if (L_Bottom > 0) {
    Bottom_chrom_buider <- cbind(Bottom_ScN[x_Bottom], rep(0, L_Bottom), rep(0, L_Bottom))
  } else {
    Bottom_chrom_buider <- c()
  }
  chrom_builder_temp <- rbind(Top_chrom_buider, chrom_builder_temp, Bottom_chrom_buider)
  SZC <- nrow(chrom_builder_temp)
  mz12C <- chrom_builder_temp[, 2]
  Chrom_Builder <- chrom_builder_temp[, -2]
  Chrom_Builder <- cbind(Chrom_Builder, Chrom_Builder[, 2], rep(0, SZC))
  dataraw <- data.frame(x=RetentionTime[Chrom_Builder[,1]], y=Chrom_Builder[,2])
  # To fill scans with the detected 12C/13C isotopologue pairs
  t_T <- table(spec_scan_xic[, 3])
  x_Tn <- which(t_T > 1)
  L_x_Tn <- length(x_Tn)
  if (L_x_Tn == 0) {
    index1 <- Chrom_Builder[, 1]%in%spec_scan_xic[, 3]
    Chrom_Builder[index1, 4] <- spec_scan_xic[, 2]
  } else {
    x_T1 <- which(t_T == 1)
    T1 <- as.numeric(names(x_T1))
    index_Chrom_T1 <- Chrom_Builder[, 1]%in%T1
    index_xic_T1 <- spec_scan_xic[, 3]%in%T1
    Chrom_Builder[index_Chrom_T1, 4] <- spec_scan_xic[index_xic_T1, 2]
    Tn <- as.numeric(names(x_Tn))
    for (t in 1:L_x_Tn) {
      xT_chrom <- which(Tn[t] == Chrom_Builder[, 1])
      index1 <- which(Tn[t] == spec_scan_xic[, 3])
      Chrom_Builder[xT_chrom, 4] <- mean(spec_scan_xic[index1, 2])
    }
  }
  data12C <- data.frame(x = RetentionTime[Chrom_Builder[, 1]], y = Chrom_Builder[, 4])
  ## Smoothing the chromatogram trace over a smoothing window
  Chrom_Builder <- data.frame(Chrom_Builder)
  colnames(Chrom_Builder) <- c("scan_number", "smooth_chrom", "raw_chrom", "C_pair")
  loess_SZC <- loess(smooth_chrom ~ scan_number, data=Chrom_Builder, span=smoothing_window/SZC, control = loess.control(surface = "direct"))
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
    ## To calculate the baseline
    BL_Base <- baseline_developer(Segment, Chrom_Builder[, 2])
    dataBL <- data.frame(x=RetentionTime[Chrom_Builder[, 1]], y = BL_Base)
    ## To merge the fronting or tailing of the peaks
    if (peak_resolving_power > 0) {
      Segment <- fronting_tailing_resolver(Segment, Chrom_Builder[, 2], smoothing_window, peak_resolving_power)
    }
    ##
    datalocmin <- unique(as.vector(Segment))
    L_d4 <- length(datalocmin)
    datalocmin <- cbind(datalocmin, rep(0, L_d4))
    for (i in 1:L_d4) {
      datalocmin[i, 2] <- Chrom_Builder[datalocmin[i, 1], 2]
      datalocmin[i, 1] <- Chrom_Builder[datalocmin[i, 1], 1]
    }
    datalocmin <- data.frame(x = RetentionTime[datalocmin[, 1]], y = datalocmin[, 2])
    ## Height adjustment
    No_isomers <- dim(Segment)[1]
    for (i in 1:No_isomers) {
      C_max_RAW <- max(Chrom_Builder[Segment[i,1]:Segment[i, 2], 3])
      if (C_max_RAW != 0) {
        C_max <- max(Chrom_Builder[Segment[i,1]:Segment[i, 2], 2])
        Chrom_Builder[(Segment[i, 1]+1):(Segment[i, 2]-1), 2] <- Chrom_Builder[(Segment[i, 1]+1):(Segment[i, 2]-1), 2]*C_max_RAW/C_max
      }
    }
    datasmooth <- data.frame(x= RetentionTime[Chrom_Builder[, 1]], y = Chrom_Builder[, 2])
    ##
    x_rtt <- which.min(abs(RetentionTime[Chrom_Builder[, 1]] - rt_target))[1]
    Int_rtt <- max(c(Chrom_Builder[x_rtt, 2], Chrom_Builder[x_rtt, 3]))
    ##
    EIC_figure <- ggplot(data = dataraw, aes(x = x)) +
      geom_segment(data=data12C, aes(x=x, xend=x, y=0, yend=y), color = "#00AFBB", size=1) +
      geom_line (data=datasmooth, aes(y=y), color="#009E73", size=2) +
      geom_line (data=dataraw, aes(y=y), color="steelblue", size=1, linetype="F1") +
      geom_line (data=dataBL, aes(y=y), color="#FC4E07", size=1.5, linetype="twodash") +
      geom_point(data=datalocmin, aes(y=y), color="black", size=5, shape=18) +
      annotation_custom(legend_EIC, xmin=RetentionTime[Chrom_Builder[floor(SZC*0.6), 1]], xmax=RetentionTime[Chrom_Builder[SZC, 1]],
                        ymin=max(Chrom_Builder[, 2])*.82, ymax=max(Chrom_Builder[, 2])*1.07) +
      xlab("Retention Time (min)") + ylab("Intensity") +
      scale_y_continuous(limits = c(0, max(Chrom_Builder[, 2])*1.07), expand = c(0, 0)) +
      scale_x_continuous(limits = c(RetentionTime[Chrom_Builder[1, 1]], RetentionTime[Chrom_Builder[SZC, 1]]), expand = c(0, 0)) +
      annotate("text", x = rt_target, y = Int_rtt*1.01, label = "*", color = "red", size = 10) +
      labs(title = file_name) +
      theme_bw() + theme(legend.position = "top", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         text = element_text(size = 24), plot.title = element_text(size = 14))
    my_table <- tableGrob(peak_property_xic, theme = ttheme_default()) ## create inset table
    EIC_figure <- grid.arrange(EIC_figure, my_table, ncol = 2) ## final result
  } else {
    x_rtt <- which.min(abs(RetentionTime[Chrom_Builder[, 1]] - rt_target))
    Int_rtt <- max(c(Chrom_Builder[x_rtt, 2], Chrom_Builder[x_rtt, 3]))
    if (Int_rtt == 0) {
      Int_rtt <- max(Chrom_Builder[, 2])*0.01
    }
    EIC_figure <- ggplot(data=dataraw, aes(x=x)) +
      geom_segment(data=data12C, aes(x=x, xend=x, y=0, yend=y), color = "#00AFBB", size=1) +
      geom_line (data=dataraw, aes(y=y), color="steelblue", size=1, linetype="F1") +
      xlim(RetentionTime[Chrom_Builder[1, 1]], RetentionTime[Chrom_Builder[SZC, 1]]) +
      ylim(0, max(Chrom_Builder[, 2])*1.05) +
      xlab("Retention Time (min)") + ylab("Intensity") +
      annotate("text", x = RetentionTime[Chrom_Builder[5, 1]], y = max(Chrom_Builder[, 2])*1.05, label=paste0("m/z = ", round(mz_target, 5))) +
      annotate("text", x = RetentionTime[Chrom_Builder[SZC - 5, 1]], y = max(Chrom_Builder[, 2])*1.02, label=paste0("RT = ", round(rt_target, 3))) +
      annotate("text", x = rt_target, y = Int_rtt*1.05, label = "*", color = "red", size = 10) +
      labs(title = file_name) +
      theme_bw() + theme(legend.position = "top", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         text = element_text(size = 24), plot.title = element_text(size = 14))
  }
  return(EIC_figure)
}
