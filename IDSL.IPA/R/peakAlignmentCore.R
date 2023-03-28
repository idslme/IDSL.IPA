peakAlignmentCore <- function(peaklistInputFolderPath, peaklistFileNames, listCorrectedRTpeaklists, massAccuracy, RTtolerance, number_processing_threads = 1) {
  ##
  nSamples <- length(peaklistFileNames)
  nSamples2 <- nSamples + 2
  ##
  ##############################################################################
  ##
  mainImzRTXcol <- do.call(rbind, lapply(1:nSamples, function(i) {
    iPeaklistFileName <- paste0(peaklistInputFolderPath, "/", peaklistFileNames[i])
    peaklist <- loadRdata(iPeaklistFileName)
    nrowPL <- nrow(peaklist)
    cbind(rep((i + 2), nrowPL), peaklist[, 8], listCorrectedRTpeaklists[[peaklistFileNames[i]]], peaklist[, 4], seq(1, nrowPL, 1))
  }))
  ##
  listCorrectedRTpeaklists <- NULL
  ##
  mainImzRTXcol <- mainImzRTXcol[!is.na(mainImzRTXcol[, 3]), ]
  mainImzRTXcol <- mainImzRTXcol[order(mainImzRTXcol[, 2], decreasing = FALSE), ]
  ##
  xDiffMZ <- c(0, which(diff(mainImzRTXcol[, 2]) > massAccuracy), nrow(mainImzRTXcol))
  LxDiffMZ <- length(xDiffMZ) - 1
  ##
  ##############################################################################
  ##
  call_peakAlignmentCore <- function(q) {
    nImzRTXcol <- xDiffMZ[q + 1] - xDiffMZ[q]
    xQ <- seq((xDiffMZ[q] + 1), xDiffMZ[q + 1], 1)
    imzRTXcol <- mainImzRTXcol[xQ, ]
    if (nImzRTXcol == 1) {
      imzRTXcol <- matrix(imzRTXcol, nrow = 1)
      orderINT <- 1
    } else {
      orderINT <- order(imzRTXcol[, 4], decreasing = TRUE)
    }
    ##
    FeatureTable <- matrix(rep(0, nSamples2*nImzRTXcol), ncol = nSamples2)
    counter <- 0
    for (i in orderINT) {
      ##
      if (imzRTXcol[i, 1] != 0) {
        counter <- counter + 1
        ##
        FeatureTable[counter, 1] <- imzRTXcol[i, 2]
        FeatureTable[counter, 2] <- imzRTXcol[i, 3]
        FeatureTable[counter, imzRTXcol[i, 1]] <- imzRTXcol[i, 5]
        x <- which((imzRTXcol[i, 1] != imzRTXcol[, 1]) &
                     (abs(imzRTXcol[i, 2] - imzRTXcol[, 2]) <= massAccuracy) &
                     (abs(imzRTXcol[i, 3] - imzRTXcol[, 3]) <= RTtolerance))
        Lx <- length(x)
        if (Lx > 0) {
          if (Lx > 1) {
            ## To remove repeated sample numbers
            x <- IPA_aggregate(idVec = imzRTXcol[x, 1], variableVec = imzRTXcol[x, 3], indexVec = x, targetVar = imzRTXcol[i, 3])
            Lx <- length(x)
          }
          ##
          for (j in 1:Lx) {
            FeatureTable[counter, imzRTXcol[x[j], 1]] <- imzRTXcol[x[j], 5]
          }
          imzRTXcol[x, ] <- 0
        }
        ##
        imzRTXcol[i, ] <- 0
      }
    }
    FeatureTable <- FeatureTable[1:counter, ]
    ##
    return(FeatureTable)
  }
  ##
  ##############################################################################
  ##
  if (number_processing_threads == 1) {
    ##
    mainPeakTable <- do.call(rbind, lapply(1:LxDiffMZ, function(q) {
      call_peakAlignmentCore(q)
    }))
    ##
  } else {
    osType <- Sys.info()[['sysname']]
    ##
    ############################################################################
    ##
    if (osType == "Windows") {
      ##
      clust <- makeCluster(number_processing_threads)
      clusterExport(clust, setdiff(ls(), c("clust", "LxDiffMZ")), envir = environment())
      ##
      mainPeakTable <- do.call(rbind, parLapply(clust, 1:LxDiffMZ, function(q) {
        call_peakAlignmentCore(q)
      }))
      ##
      stopCluster(clust)
      ##
      ##########################################################################
      ##
    } else {
      ##
      mainPeakTable <- do.call(rbind, mclapply(1:LxDiffMZ, function(q) {
        call_peakAlignmentCore(q)
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
      ##
    }
  }
  ##
  mainImzRTXcol <- NULL
  ##
  ##############################################################################
  ## To resolve redundant peaks in the peak matrix table
  x_s <- peakPropertyTableFreqCalculator(mainPeakTable, startColumnIndex = 3, number_processing_threads, allowedVerbose = FALSE)
  mainPeakTable <- cbind(mainPeakTable[, 1:2], x_s, mainPeakTable[, 3:nSamples2])
  mainPeakTable <- mainPeakTable[order(mainPeakTable[, 1], decreasing = FALSE), ]
  ##
  xDiffQ <- c(0, which(diff(mainPeakTable[, 1]) > massAccuracy), nrow(mainPeakTable))
  LxDiffQ <- length(xDiffQ) - 1
  ##
  progressBARboundaries <- txtProgressBar(min = 0, max = LxDiffQ, initial = 0, style = 3)
  ##
  for (q in 1:LxDiffQ) {
    xQ <- seq((xDiffQ[q] + 1), xDiffQ[q + 1], 1)
    xQ <- xQ[order(mainPeakTable[xQ, 3], decreasing = TRUE)]
    ##
    for (i in xQ) {
      if (mainPeakTable[i, 1] != 0) {
        xC <- which(abs(mainPeakTable[i, 1] - mainPeakTable[xQ, 1]) <= massAccuracy &
                      abs(mainPeakTable[i, 2] - mainPeakTable[xQ, 2]) <= RTtolerance)
        ##
        if (length(xC) > 1) {
          xC <- xQ[xC]
          ##
          xDiffxC <- setdiff(xC, i)
          if (mainPeakTable[i, 3] < nSamples) {
            table_c <- mainPeakTable[xC, ]
            ##
            x_table_main0 <- which(table_c[1, ] == 0)
            for (j in x_table_main0) {
              xNon0 <- which(table_c[, j] > 0)
              LxNon0 <- length(xNon0)
              if (LxNon0 > 0) {
                if (LxNon0 > 1) {
                  xMin <- mzRTindexer(MZvec = table_c[xNon0, 1], RTvec = table_c[xNon0, 2],
                                      MZref = table_c[1, 1], RTref = table_c[1, 2], massAccuracy, RTtolerance)
                  ##
                  xNon0 <- xNon0[xMin[1]]
                }
                table_c[1, j] <- table_c[xNon0, j]
                table_c[1, 3] <- table_c[1, 3] + 1
              }
            }
            mainPeakTable[i, ] <- table_c[1, ]
          }
          mainPeakTable[xDiffxC, ] <- 0
        }
      }
    }
    ##
    setTxtProgressBar(progressBARboundaries, q)
  }
  close(progressBARboundaries)
  ##
  xNon0 <- which(mainPeakTable[, 1] != 0)
  mainPeakTable <- mainPeakTable[xNon0, ]
  rownames(mainPeakTable) <- NULL
  ##
  mainPeakTable[, 1] <- round(mainPeakTable[, 1], digits = 6)
  mainPeakTable[, 2] <- round(mainPeakTable[, 2], digits = 4)
  ##
  return(mainPeakTable)
}
