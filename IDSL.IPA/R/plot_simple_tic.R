plot_simple_tic <- function(filelist, filelocation, number_processing_threads = 1, plotTitle = "Total Ion Chromatogram") {
  ##
  call_plot_simple_tic <- function(i) {
    p2l <- IDSL.MXP::peak2list(filelocation, i)
    scanTable <- p2l[["scanTable"]]
    df <- data.frame(RT = scanTable$retentionTime, Intensity = scanTable$totIonCurrent)
    return(df)
  }
  ##
  ##############################################################################
  ##
  if (number_processing_threads == 1) {
    ##
    dflist.tic <- lapply(filelist, function(i) {
      call_plot_simple_tic(i)
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
      dflist.tic <- parLapply(clust, filelist, function(i) {
        call_plot_simple_tic(i)
      })
      ##
      stopCluster(clust)
      ##
      ##########################################################################
      ##
    } else {
      ##
      dflist.tic <- mclapply(filelist, function(i) {
        call_plot_simple_tic(i)
      }, mc.cores = number_processing_threads)
      ##
      closeAllConnections()
      ##
    }
  }
  ##
  ##############################################################################
  ##
  rainbowColors <- rainbow(length(filelist), alpha = 1)
  png(paste0(filelocation, "/_simple_tic.png"), width = 12, height = 8, units = "in", res = 100)
  ##
  for (kk in 1:length(filelist)) {
    df <- dflist.tic[[kk]]
    if (kk == 1) {
      inten.maxvec <- sapply(1:length(dflist.tic), function(x){max(dflist.tic[[x]][, 2])})
      inten.minvec <- sapply(1:length(dflist.tic), function(x){min(dflist.tic[[x]][, 2])})
      ##
      rt.maxvec <- sapply(1:length(dflist.tic), function(x){max(dflist.tic[[x]][, 1])})
      rt.minvec <- sapply(1:length(dflist.tic), function(x){min(dflist.tic[[x]][, 1])})
      ##
      plot(df$RT, df$Intensity, ylim = c(min(inten.minvec), max(inten.maxvec)), xlim = c(min(rt.minvec), max(rt.maxvec)), type = "l", lty = 1, lwd = 2, frame = T, pch = 19, col = "white", ylab = "Intensity",xlab = "RT (min)", cex.lab = 3, cex.axis = 2)
      title(main = plotTitle, cex.main = 1.5)
    } else {
      lines(df$RT, df$Intensity, pch = 19, col = rainbowColors[kk], type = "l", lty = 1, lwd = 2)
    }
  }
  ##
  dev.off()
  ##
  IPA_logRecorder("simple TICs have been generated!")
}
