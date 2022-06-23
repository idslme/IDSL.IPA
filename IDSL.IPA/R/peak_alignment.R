peak_alignment <- function(input_path_pl, file_names_pl, RT_pl, mz_error, rt_tol, n_quantile, number_processing_cores) {
  ##
  L_PL <- length(file_names_pl)
  L_PL2 <- L_PL + 2
  L_PL3 <- L_PL + 3
  ##
  imzRTXcol_main_call <- function(i) {
    peaklist <- loadRdata(paste0(input_path_pl, "/", file_names_pl[i]))
    cbind(rep(i, nrow(peaklist)), peaklist[, 8], RT_pl[[i]], peaklist[, 4], 1:nrow(peaklist))
  }
  ##
  FeatureTable_Main_call <- function(q) {
    x_Q <- which(imzRTXcol_main[, 2] >= MZ_Q_boundaries[q, 1] &
                   imzRTXcol_main[, 2] <= MZ_Q_boundaries[q, 2])
    imzRTXcol <- imzRTXcol_main[x_Q, ]
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
        x <- which((abs(imzRTXcol[i, 2] - imzRTXcol[, 2]) <= mz_error) &
                     (abs(imzRTXcol[i, 3] - imzRTXcol[, 3]) <= rt_tol) &
                     (imzRTXcol[i, 1] != imzRTXcol[, 1]) & (imzRTXcol[, 1] != 0))
        ##
        if (length(x) > 0) {
          iSamples <- imzRTXcol[x, 1]
          if (length(iSamples) != length(unique(iSamples))) {
            xi <- table(iSamples)
            xii <- which(xi > 1)
            Rx <- as.numeric(names(xii))  # Repeated peaks in the same sample
            A <- cbind(imzRTXcol[x, ], x)
            xiii <- do.call(c, lapply(1:length(Rx), function(j) {
              xj <- which(A[, 1] == Rx[j])
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
          imzRTXcol[x, 1] <- 0
        }
        ##
        imzRTXcol[i, 1] <- 0
      }
    }
    FeatureTable <- FeatureTable[1:counter, ]
    ##
    return(FeatureTable)
  }
  ##############################################################################
  if (number_processing_cores == 1) {
    ##
    imzRTXcol_main <- do.call(rbind, lapply(1:L_PL, function(i) {
      imzRTXcol_main_call(i)
    }))
    ##
    RT_pl <- c()
    ##
    imzRTXcol_main <- imzRTXcol_main[order(imzRTXcol_main[, 2], decreasing = TRUE), ]
    ##
    if (n_quantile > 1) {
      MZ_Q <- quantile(imzRTXcol_main[, 2], probs = c(1:n_quantile)/n_quantile)
      MZ_Q_boundaries <- cbind(c(min(imzRTXcol_main[, 2]), MZ_Q[1:(n_quantile - 1)]), MZ_Q)
      MZ_Q_boundaries[, 1] <- MZ_Q_boundaries[, 1] - mz_error*1.5
      MZ_Q_boundaries[, 2] <- MZ_Q_boundaries[, 2] + mz_error*1.5
    } else {
      MZ_Q_boundaries <- matrix(c(min(imzRTXcol_main[, 2]), max(imzRTXcol_main[, 2])), ncol = 2)
    }
    ##
    FeatureTable_Main <- do.call(rbind, lapply(1:n_quantile, function(q) {
      FeatureTable_Main_call(q)
    }))
    ##
    imzRTXcol_main <- c()
    ##
    L_FTmain <- dim(FeatureTable_Main)[1]
    x_s <- do.call(c, lapply(1:L_FTmain, function(i) {
      length(which(FeatureTable_Main[i, 3:L_PL2] > 0))
    }))
    ##
  } else {
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Linux") {
      ##
      imzRTXcol_main <- do.call(rbind, mclapply(1:L_PL, function(i) {
        imzRTXcol_main_call(i)
      }, mc.cores = number_processing_cores))
      ##
      RT_pl <- c()
      ##
      imzRTXcol_main <- imzRTXcol_main[order(imzRTXcol_main[, 2], decreasing = TRUE), ]
      ##
      if (n_quantile > 1) {
        MZ_Q <- quantile(imzRTXcol_main[, 2], probs = c(1:n_quantile)/n_quantile)
        MZ_Q_boundaries <- cbind(c(min(imzRTXcol_main[, 2]), MZ_Q[1:(n_quantile - 1)]), MZ_Q)
        MZ_Q_boundaries[, 1] <- MZ_Q_boundaries[, 1] - mz_error*1.5
        MZ_Q_boundaries[, 2] <- MZ_Q_boundaries[, 2] + mz_error*1.5
      } else {
        MZ_Q_boundaries <- matrix(c(min(imzRTXcol_main[, 2]), max(imzRTXcol_main[, 2])), ncol = 2)
      }
      ##
      FeatureTable_Main <- do.call(rbind, mclapply(1:n_quantile, function(q) {
        FeatureTable_Main_call(q)
      }, mc.cores = number_processing_cores))
      ##
      imzRTXcol_main <- c()
      ##
      L_FTmain <- dim(FeatureTable_Main)[1]
      x_s <- do.call(c, mclapply(1:L_FTmain, function(i) {
        length(which(FeatureTable_Main[i, 3:L_PL2] > 0))
      }, mc.cores = number_processing_cores))
      ##
      closeAllConnections()
      ##
    } else if (osType == "Windows") {
      cl <- makeCluster(number_processing_cores)
      registerDoParallel(cl)
      ##
      imzRTXcol_main <- foreach(i = 1:L_PL, .combine="rbind", .verbose = FALSE) %dopar% {
        imzRTXcol_main_call(i)
      }
      RT_pl <- c()
      ##
      imzRTXcol_main <- imzRTXcol_main[(imzRTXcol_main[, 2] > 0), ]
      ##
      imzRTXcol_main <- imzRTXcol_main[order(imzRTXcol_main[, 2], decreasing = TRUE), ]
      ##
      if (n_quantile > 1) {
        MZ_Q <- quantile(imzRTXcol_main[, 2], probs = c(1:n_quantile)/n_quantile)
        MZ_Q_boundaries <- cbind(c(min(imzRTXcol_main[, 2]), MZ_Q[1:(n_quantile - 1)]), MZ_Q)
        MZ_Q_boundaries[, 1] <- MZ_Q_boundaries[, 1] - mz_error*1.5
        MZ_Q_boundaries[, 2] <- MZ_Q_boundaries[, 2] + mz_error*1.5
      } else {
        MZ_Q_boundaries <- matrix(c(min(imzRTXcol_main[, 2]), max(imzRTXcol_main[, 2])), ncol = 2)
      }
      ##
      FeatureTable_Main <- foreach(q = 1:n_quantile, .combine="rbind", .verbose = FALSE) %dopar% {
        FeatureTable_Main_call(q)
      }
      ##
      imzRTXcol_main <- c()
      ##
      L_FTmain <- dim(FeatureTable_Main)[1]
      x_s <- foreach(i = 1:L_FTmain, .combine="c", .verbose = FALSE) %dopar% {
        length(which(FeatureTable_Main[i, 3:L_PL2] > 0))
      }
      ##
      stopCluster(cl)
    }
  }
  ##############################################################################
  ## To resolve redundant peaks in the peak matrix table
  FeatureTable_Main <- cbind(x_s, FeatureTable_Main)
  FeatureTable_Main <- FeatureTable_Main[order(FeatureTable_Main[, 1], decreasing = TRUE), ]
  ##
  progressBARboundaries <- txtProgressBar(min = 1, max = L_FTmain, initial = 1, style = 3)
  for (i in 1:L_FTmain) {
    setTxtProgressBar(progressBARboundaries, i)
    x_c <- which(abs(FeatureTable_Main[i, 2] - FeatureTable_Main[, 2]) <= mz_error &
                   abs(FeatureTable_Main[i, 3] - FeatureTable_Main[, 3]) <= rt_tol &
                   FeatureTable_Main[i, 1] != 0)
    ##
    if (length(x_c) > 1) {
      x_diff <- setdiff(x_c, i)
      if (FeatureTable_Main[i, 1] < L_PL) {
        table_c <- do.call(rbind, lapply(x_c, function(j) {
          FeatureTable_Main[j, 1:L_PL3]
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
        FeatureTable_Main[i, 4:L_PL3] <- table_c[1, 4:L_PL3]
      }
      FeatureTable_Main[x_diff, ] <- 0
    }
  }
  close(progressBARboundaries)
  x_non0 <- which(FeatureTable_Main[, 1] != 0)
  FeatureTable_Main <- FeatureTable_Main[x_non0, ]
  FeatureTable_Main <- FeatureTable_Main[, -1]
  rownames(FeatureTable_Main) <- c()
  FeatureTable_Main <- FeatureTable_Main[order(FeatureTable_Main[, 1], decreasing = FALSE), ]
  return(FeatureTable_Main)
}
