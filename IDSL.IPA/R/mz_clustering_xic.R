mz_clustering_xic <- function(spec_scan, mass_accuracy_xic, min_peak_height, min_nIsoPair) {
  MZ_Int_ScN <- spec_scan[, 1:3]
  ##
  maxLengthClusters <- floor(nrow(MZ_Int_ScN)/max(c(1, min_nIsoPair)))
  index_xic <- vector(mode = "list", maxLengthClusters)
  mzxic <- rep(0, maxLengthClusters)
  Counter <- 0
  L_MinHeight <- length(which(MZ_Int_ScN[, 2] >= min_peak_height))
  i <- 1
  if (min_nIsoPair > 1) {
    for (i in 1:L_MinHeight) {
      if (MZ_Int_ScN[i, 1] != 0) {
        x1 <- which(abs(MZ_Int_ScN[, 1] - MZ_Int_ScN[i, 1]) <= mass_accuracy_xic) # to cluster m/z in consecutive scans
        L_x1 <- length(x1)
        if (L_x1 >= min_nIsoPair) {         # Minimum number of scans
          A <- MZ_Int_ScN[x1, ]
          ## To remove repeated scans
          if (L_x1 != length(unique(A[, 3]))) {
            A <- cbind(A, x1)
            A <- A[order(A[, 3], decreasing = FALSE), ]
            xDiff0 <- which(diff(A[, 3]) == 0)
            xDiff1 <- which(diff(xDiff0) > 1)
            xu <- c(0, xDiff1, length(xDiff0))
            xRemove <- do.call(c, lapply(1:(length(xu) - 1), function(j) {
              xuDiff <- xDiff0[(xu[j] + 1):xu[j + 1]]
              xuR <- c(xuDiff, (xuDiff[length(xuDiff)] + 1))  # Repeated scan numbers
              xMin <- which.min(abs(A[xuR, 1] - MZ_Int_ScN[i, 1]))[1]
              xR <- xuR[-xMin]
              A[xR, 4]
            }))
            x1 <- setdiff(x1, xRemove)
          }
          ##
          if (length(x1) >= min_nIsoPair) {
            Counter <- Counter + 1
            index_xic[[Counter]] <- x1
            mzxic[Counter] <- MZ_Int_ScN[i, 1]
            MZ_Int_ScN[x1, 1] <- 0
          }
        }
      }
    }
  } else {
    for (i in 1:L_MinHeight) {
      if (MZ_Int_ScN[i, 1] != 0) {
        x1 <- which(abs(MZ_Int_ScN[, 1] - MZ_Int_ScN[i, 1]) <= mass_accuracy_xic) # to cluster m/z in consecutive scans
        L_x1 <- length(x1)
        if (L_x1 >= min_nIsoPair) {         # Minimum number of scans
          A <- matrix(MZ_Int_ScN[x1, ], ncol = 3)
          ## To remove repeated scans
          if (L_x1 != length(unique(A[, 3]))) {
            A <- cbind(A, x1)
            A <- A[order(A[, 3], decreasing = FALSE), ]
            xDiff0 <- which(diff(A[, 3]) == 0)
            xDiff1 <- which(diff(xDiff0) > 1)
            xu <- c(0, xDiff1, length(xDiff0))
            xRemove <- do.call(c, lapply(1:(length(xu) - 1), function(j) {
              xuDiff <- xDiff0[(xu[j] + 1):xu[j + 1]]
              xuR <- c(xuDiff, (xuDiff[length(xuDiff)] + 1))  # Repeated scan numbers
              xMin <- which.min(abs(A[xuR, 1] - MZ_Int_ScN[i, 1]))[1]
              xR <- xuR[-xMin]
              A[xR, 4]
            }))
            x1 <- setdiff(x1, xRemove)
          }
          ##
          Counter <- Counter + 1
          index_xic[[Counter]] <- x1
          mzxic[Counter] <- MZ_Int_ScN[i, 1]
          MZ_Int_ScN[x1, 1] <- 0
        }
      }
    }
  }
  index_xic <- index_xic[1:Counter]
  mzxic <- mzxic[1:Counter]
  names(index_xic) <- round(mzxic, 5)
  return(index_xic)
}
