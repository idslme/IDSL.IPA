plot_simple_tic <- function(filelist,filelocation, number_processing_threads = 1, plotTitle = "Total Ion Chromatogram") {



  clust <- makeCluster(number_processing_threads)
  registerDoParallel(clust)
  dflist.tic <- foreach(mzmlfile = filelist) %dopar% {
    p2l <- IDSL.MXP::peak2list(filelocation, mzmlfile)
    scanTable <- p2l[["scanTable"]]
    data.frame(RT=scanTable$retentionTime, Intensity=scanTable$totIonCurrent)
  }
  stopCluster(clust)



  colfunc <- colorRampPalette(c("yellow", "red"))
  colvec <- colfunc(499)



  png(paste0(filelocation,"/_simple_tic.png"),width=12,height =8,units = "in", res = 200)



  for(kk in 1:length(filelist) ) {
    df <- dflist.tic[[kk]]



    if(kk == 1) {



      inten.maxvec <- sapply(1:length(dflist.tic), function(x){max(dflist.tic[[x]][,2])})
      inten.minvec <- sapply(1:length(dflist), function(x){min(dflist.tic[[x]][,2])})



      rt.maxvec <- sapply(1:length(dflist.tic), function(x){max(dflist.tic[[x]][,1])})
      rt.minvec <- sapply(1:length(dflist.tic), function(x){min(dflist.tic[[x]][,1])})



      plot(df$RT, df$Intensity, ylim = c(min(inten.minvec),max(inten.maxvec)), xlim = c(min(rt.minvec),max(rt.maxvec)), type = "l",lty=1, lwd=2, frame = T, pch = 19,col = "white", ylab="Intensity",xlab="RT (min)", cex.lab=3, cex.axis = 2)
      title(main = plotTitle,cex.main = 1.5)
    } else {
      #lines(df$RT, df$Intensity, pch = 19, col = colvec[kk], type = "l", lty = 1,lwd=2)
      lines(df$RT, df$Intensity, pch = 19, col = sample(rainbow(500),1), type = "l", lty = 1,lwd=2)
    }
  }
  print("simple TICs have been generated!")
}
