peak_alignment <- function(peaklistInputFolderPath, peaklistFileNames, listCorrectedRTpeaklists, massAccuracy, RTtolerance, noQuantile, number_processing_threads = 1) {
  ##
  L_PL <- length(peaklistFileNames)
  L_PL2 <- L_PL + 2
  L_PL3 <- L_PL + 3
  ##
  call_mainImzRTXcol <- function(i) {
    iPeaklistFileName <- paste0(peaklistInputFolderPath, "/", peaklistFileNames[i])
    peaklist <- loadRdata(iPeaklistFileName)
    nrowPL <- nrow(peaklist)
    cbind(rep(i, nrowPL), peaklist[, 8], listCorrectedRTpeaklists[[peaklistFileNames[i]]], peaklist[, 4], seq(1, nrowPL, 1))
  }
  ##
  call_mainPeakTable <- function(q) {
    x_Q <- which(mainImzRTXcol[, 2] >= MZ_Q_boundaries[q, 1] &
                   mainImzRTXcol[, 2] <= MZ_Q_boundaries[q, 2])
    imzRTXcol <- mainImzRTXcol[x_Q, ]
    imzRTXcol <- imzRTXcol[order(imzRTXcol[, 4], decreasing = TRUE), ]
    N_imzRTXcol <- length(x_Q)
    FeatureTable <- matrix(rep(0, L_PL2*N_imzRTXcol), ncol = L_PL2)
    counter <- 0
    for (i in 1:N_imzRTXcol) {
      ##
      if (imzRTXcol[i, 1] != 0) {
        counter <- counter + 1
        ##
        FeatureTable[counter, 1] <- imzRTXcol[i, 2]
        FeatureTable[counter, 2] <- imzRTXcol[i, 3]
        FeatureTable[counter, (imzRTXcol[i, 1] + 2)] <- imzRTXcol[i, 5]
        x <- which((abs(imzRTXcol[i, 2] - imzRTXcol[, 2]) <= massAccuracy) &
                     (abs(imzRTXcol[i, 3] - imzRTXcol[, 3]) <= RTtolerance) &
                     (imzRTXcol[i, 1] != imzRTXcol[, 1]))
        ##
        if (length(x) > 0) {
          iSamples <- imzRTXcol[x, 1]
          if (length(iSamples) != length(unique(iSamples))) {
            xi <- table(iSamples)
            xii <- which(xi > 1)
            Rx <- as.numeric(names(xii))  # Repeated peaks in the same sample
            A <- cbind(imzRTXcol[x, ], x)
            xiii <- do.call(c, lapply(Rx, function(j) {
              xj <- which(A[, 1] == j)
              xjj <- which.min(abs(A[xj, 3] - imzRTXcol[i, 3]))
              A[xj[-xjj], 6]
            }))
            x <- setdiff(x, xiii)
            iSamples <- imzRTXcol[x, 1]
          }
          ##
          for (j in 1:length(iSamples)) {
            FeatureTable[counter, iSamples[j] + 2] <- imzRTXcol[x[j], 5]
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
  ##############################################################################
  if (number_processing_threads == 1) {
    ##
    mainImzRTXcol <- do.call(rbind, lapply(1:L_PL, function(i) {
      call_mainImzRTXcol(i)
    }))
    ##
    listCorrectedRTpeaklists <- NULL
    ##
    mainImzRTXcol <- mainImzRTXcol[!is.na(mainImzRTXcol[, 3]), ]
    mainImzRTXcol <- mainImzRTXcol[order(mainImzRTXcol[, 2], decreasing = TRUE), ]
    ##
    if (noQuantile > 1) {
      MZ_Q <- quantile(mainImzRTXcol[, 2], probs = c(1:noQuantile)/noQuantile)
      MZ_Q_boundaries <- cbind(c(min(mainImzRTXcol[, 2]), MZ_Q[1:(noQuantile - 1)]), MZ_Q)
      MZ_Q_boundaries[, 1] <- MZ_Q_boundaries[, 1] - massAccuracy*1.5
      MZ_Q_boundaries[, 2] <- MZ_Q_boundaries[, 2] + massAccuracy*1.5
    } else {
      MZ_Q_boundaries <- matrix(c(min(mainImzRTXcol[, 2]), max(mainImzRTXcol[, 2])), ncol = 2)
    }
    ##
    mainPeakTable <- do.call(rbind, lapply(1:noQuantile, function(q) {
      call_mainPeakTable(q)
    }))
    ##
    mainImzRTXcol <- NULL
    ##
    L_FTmain <- dim(mainPeakTable)[1]
    x_s <- do.call(c, lapply(1:L_FTmain, function(i) {
      length(which(mainPeakTable[i, 3:L_PL2] > 0))
    }))
    ##
  } else {
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Linux") {
      ##
      mainImzRTXcol <- do.call(rbind, mclapply(1:L_PL, function(i) {
        call_mainImzRTXcol(i)
      }, mc.cores = number_processing_threads))
      ##
      listCorrectedRTpeaklists <- NULL
      ##
      mainImzRTXcol <- mainImzRTXcol[!is.na(mainImzRTXcol[, 3]), ]
      mainImzRTXcol <- mainImzRTXcol[order(mainImzRTXcol[, 2], decreasing = TRUE), ]
      ##
      if (noQuantile > 1) {
        MZ_Q <- quantile(mainImzRTXcol[, 2], probs = c(1:noQuantile)/noQuantile)
        MZ_Q_boundaries <- cbind(c(min(mainImzRTXcol[, 2]), MZ_Q[1:(noQuantile - 1)]), MZ_Q)
        MZ_Q_boundaries[, 1] <- MZ_Q_boundaries[, 1] - massAccuracy*1.5
        MZ_Q_boundaries[, 2] <- MZ_Q_boundaries[, 2] + massAccuracy*1.5
      } else {
        MZ_Q_boundaries <- matrix(c(min(mainImzRTXcol[, 2]), max(mainImzRTXcol[, 2])), ncol = 2)
      }
      ##
      mainPeakTable <- do.call(rbind, mclapply(1:noQuantile, function(q) {
        call_mainPeakTable(q)
      }, mc.cores = number_processing_threads))
      ##
      mainImzRTXcol <- NULL
      ##
      L_FTmain <- dim(mainPeakTable)[1]
      x_s <- do.call(c, mclapply(1:L_FTmain, function(i) {
        length(which(mainPeakTable[i, 3:L_PL2] > 0))
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
      ##
    } else if (osType == "Windows") {
      cl <- makeCluster(number_processing_threads)
      registerDoParallel(cl)
      ##
      mainImzRTXcol <- foreach(i = 1:L_PL, .combine = 'rbind', .verbose = FALSE) %dopar% {
        call_mainImzRTXcol(i)
      }
      listCorrectedRTpeaklists <- NULL
      ##
      mainImzRTXcol <- mainImzRTXcol[!is.na(mainImzRTXcol[, 3]), ]
      mainImzRTXcol <- mainImzRTXcol[(mainImzRTXcol[, 2] > 0), ]
      ##
      mainImzRTXcol <- mainImzRTXcol[order(mainImzRTXcol[, 2], decreasing = TRUE), ]
      ##
      if (noQuantile > 1) {
        MZ_Q <- quantile(mainImzRTXcol[, 2], probs = c(1:noQuantile)/noQuantile)
        MZ_Q_boundaries <- cbind(c(min(mainImzRTXcol[, 2]), MZ_Q[1:(noQuantile - 1)]), MZ_Q)
        MZ_Q_boundaries[, 1] <- MZ_Q_boundaries[, 1] - massAccuracy*1.5
        MZ_Q_boundaries[, 2] <- MZ_Q_boundaries[, 2] + massAccuracy*1.5
      } else {
        MZ_Q_boundaries <- matrix(c(min(mainImzRTXcol[, 2]), max(mainImzRTXcol[, 2])), ncol = 2)
      }
      ##
      mainPeakTable <- foreach(q = 1:noQuantile, .combine = 'rbind', .verbose = FALSE) %dopar% {
        call_mainPeakTable(q)
      }
      ##
      mainImzRTXcol <- NULL
      ##
      L_FTmain <- dim(mainPeakTable)[1]
      x_s <- foreach(i = 1:L_FTmain, .combine = 'c', .verbose = FALSE) %dopar% {
        length(which(mainPeakTable[i, 3:L_PL2] > 0))
      }
      ##
      stopCluster(cl)
    }
  }
  ##############################################################################
  ## To resolve redundant peaks in the peak matrix table
  mainPeakTable <- cbind(x_s, mainPeakTable)
  mainPeakTable <- mainPeakTable[order(mainPeakTable[, 1], decreasing = TRUE), ]
  ##
  progressBARboundaries <- txtProgressBar(min = 0, max = L_FTmain, initial = 0, style = 3)
  for (i in 1:L_FTmain) {
    setTxtProgressBar(progressBARboundaries, i)
    if (mainPeakTable[i, 1] != 0) {
      x_c <- which(abs(mainPeakTable[i, 2] - mainPeakTable[, 2]) <= massAccuracy &
                     abs(mainPeakTable[i, 3] - mainPeakTable[, 3]) <= RTtolerance)
      ##
      if (length(x_c) > 1) {
        x_diff <- setdiff(x_c, i)
        if (mainPeakTable[i, 1] < L_PL) {
          table_c <- do.call(rbind, lapply(x_c, function(j) {
            mainPeakTable[j, 1:L_PL3]
          }))
          x_table_main0 <- which(table_c[1, ] == 0)
          for (j in x_table_main0) {
            x_non0 <- which(table_c[, j] > 0)
            if (length(x_non0) > 0) {
              if (length(x_non0) > 1) {
                x_min <- which.min(abs(table_c[1, 3] - table_c[x_non0, 3]))
                x_non0 <- x_non0[x_min[1]]
              }
              table_c[1, j] <- table_c[x_non0, j]
            }
          }
          mainPeakTable[i, 4:L_PL3] <- table_c[1, 4:L_PL3]
        }
        mainPeakTable[x_diff, ] <- 0
      }
    }
  }
  close(progressBARboundaries)
  ##
  x_non0 <- which(mainPeakTable[, 1] != 0)
  mainPeakTable <- mainPeakTable[x_non0, ]
  mainPeakTable <- mainPeakTable[, -1]
  rownames(mainPeakTable) <- NULL
  mainPeakTable <- mainPeakTable[order(mainPeakTable[, 1], decreasing = FALSE), ]
  ##
  mainPeakTable[, 1] <- round(mainPeakTable[, 1], digits = 6)
  mainPeakTable[, 2] <- round(mainPeakTable[, 2], digits = 4)
  ##
  return(mainPeakTable)
}
