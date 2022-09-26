plot_mz_eic <- function(filelist,filelocation,mztarget,mzdelta, number_processing_threads = 1, rtstart=0,rtend=0,plotTitle = "") {


  clust <- makeCluster(number_processing_threads)
  registerDoParallel(clust)



  dflist <- foreach(mzmlfile = filelist) %dopar% {
    p2l <- IDSL.MXP::peak2list(filelocation, mzmlfile)
    scanTable <- p2l[["scanTable"]]
    spectraList <- p2l[["spectraList"]]
    mzdf <- do.call(rbind, lapply(1:length(spectraList),function(l) { cbind(spectraList[[l]],scanTable$retentionTime[l],l) }))
    df <- data.frame(RT = as.numeric(mzdf[, 3]), Intensity = as.numeric(mzdf[, 2]), mz = as.numeric(mzdf[, 1]), scan = as.numeric(mzdf[, 4]))
    df <- df[which(df$mz < mztarget + mzdelta & df$mz > mztarget - mzdelta),]
    df
  }
  stopCluster(clust)



  inten.maxvec <- sapply(1:length(dflist), function(x){max(dflist[[x]][,2])})
  inten.minvec <- sapply(1:length(dflist), function(x){min(dflist[[x]][,2])})



  rt.maxvec <- sapply(1:length(dflist), function(x){max(dflist[[x]][,1])})
  rt.minvec <- sapply(1:length(dflist), function(x){min(dflist[[x]][,1])})



  rtmax <- max(rt.maxvec)
  rtmin <- min(rt.minvec)



  if (rtstart != 0 & rtend != 0) {
    rtmax = rtend
    rtmin = rtstart
  }



  colfunc <- colorRampPalette(c("white", "black"))
  colvec <- colfunc(length(filelist))



  ptitle <- plotTitle



  if (plotTitle == "") {
    ptitle <- paste0(mztarget," (+/- ",mzdelta," Da)")
  }



  png(paste0(filelocation,"/_simple_eic.png"),width=12,height =8,units = "in", res = 200)



  for(kk in 1:length(filelist) ) {
    df <- dflist[[kk]]
    if(kk == 1) {
      plot(df$RT, df$Intensity, ylim = c(min(inten.minvec),max(inten.maxvec)), xlim = c(rtmin,rtmax), type = "l",lty=1, lwd=2, frame = T, pch = 19,col = "black", ylab="Intensity",xlab="RT (min)", cex.lab=4, cex.axis = 2)
      title(main = ptitle,cex.main = 1.5)
    } else {
      lines(df$RT, df$Intensity, pch = 19, col = colvec[kk], type = "l", lty = 1,lwd=2)
    }
  }
  IPA_logRecorder('simple EICs have been successfully generated!')
}
