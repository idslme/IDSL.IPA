plot_mz_eic <- function(filelist, filelocation, mzTarget, massAccuracy, number_processing_threads = 1, rtstart = 0, rtend = 0, plotTitle = "") {
  ##
  call_plot_mz_eic <- function(i) {
    p2l <- IDSL.MXP::peak2list(filelocation, i)
    scanTable <- p2l[["scanTable"]]
    spectraList <- p2l[["spectraList"]]
    mzdf <- do.call(rbind, lapply(1:length(spectraList),function(l) { cbind(spectraList[[l]],scanTable$retentionTime[l],l) }))
    df <- data.frame(RT = as.numeric(mzdf[, 3]), Intensity = as.numeric(mzdf[, 2]), mz = as.numeric(mzdf[, 1]), scan = as.numeric(mzdf[, 4]))
    df <- df[which(df$mz < mzTarget + massAccuracy & df$mz > mzTarget - massAccuracy), ]
    return(df)
  }
  ##
  ##############################################################################
  ##
  if (number_processing_threads == 1) {
    ##
    dflist <- lapply(filelist, function(i) {
      call_plot_mz_eic(i)
    })
    ##
  } else {
    osType <- Sys.info()[['sysname']]
    ##
    ############################################################################
    ##
    if (osType == "Windows") {
      ##
      clust <- makeCluster(number_processing_threads)
      clusterExport(clust, setdiff(ls(), c("clust", "filelist")), envir = environment())
      ##
      dflist <- parLapply(clust, filelist, function(i) {
        call_plot_mz_eic(i)
      })
      ##
      stopCluster(clust)
      ##
      ##########################################################################
      ##
    } else {
      ##
      dflist <- mclapply(filelist, function(i) {
        call_plot_mz_eic(i)
      }, mc.cores = number_processing_threads)
      ##
      closeAllConnections()
      ##
    }
  }
  ##
  ##############################################################################
  ##
  inten.maxvec <- sapply(1:length(dflist), function(x){max(dflist[[x]][, 2])})
  inten.minvec <- sapply(1:length(dflist), function(x){min(dflist[[x]][, 2])})
  ##
  rt.maxvec <- sapply(1:length(dflist), function(x){max(dflist[[x]][, 1])})
  rt.minvec <- sapply(1:length(dflist), function(x){min(dflist[[x]][, 1])})
  ##
  rtmax <- max(rt.maxvec)
  rtmin <- min(rt.minvec)
  if ((rtstart != 0) & (rtend != 0)) {
    rtmax = rtend
    rtmin = rtstart
  }
  ##
  ptitle <- plotTitle
  if (plotTitle == "") {
    ptitle <- paste0(mzTarget," (+/- ",massAccuracy," Da)")
  }
  ##
  rainbowColors <- rainbow(length(filelist), alpha = 1)
  png(paste0(filelocation, "/_simple_eic.png"), width = 12, height = 8, units = "in", res = 100)
  ##
  for (kk in 1:length(filelist)) {
    df <- dflist[[kk]]
    if (kk == 1) {
      plot(df$RT, df$Intensity, ylim = c(min(inten.minvec), max(inten.maxvec)), xlim = c(rtmin, rtmax), type = "l", lty = 1, lwd = 2, frame = T, pch = 19, col = "black", ylab = "Intensity", xlab = "RT (min)", cex.lab = 4, cex.axis = 2)
      title(main = ptitle, cex.main = 1.5)
    } else {
      lines(df$RT, df$Intensity, pch = 19, col = rainbowColors[kk], type = "l", lty = 1, lwd = 2)
    }
  }
  ##
  dev.off()
  ##
  IPA_logRecorder('simple EICs have been successfully generated!')
}
