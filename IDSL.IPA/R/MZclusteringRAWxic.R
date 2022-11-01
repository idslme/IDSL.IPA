MZclusteringRAWxic <- function(spectraScan, massAccuracyXIC, minPeakHeight, minNIonPair) {
  MZ_Int_ScN <- spectraScan[, 1:3]
  ##
  maxLengthClusters <- floor(nrow(MZ_Int_ScN)/max(c(1, minNIonPair)))
  indexXIC <- vector(mode = "list", maxLengthClusters)
  mzxic <- rep(0, maxLengthClusters)
  Counter <- 0
  L_MinHeight <- length(which(MZ_Int_ScN[, 2] >= minPeakHeight))
  ##
  if (minNIonPair > 1) {
    for (i in 1:L_MinHeight) {
      if (MZ_Int_ScN[i, 1] != 0) {
        x1 <- which(abs(MZ_Int_ScN[, 1] - MZ_Int_ScN[i, 1]) <= massAccuracyXIC) # to cluster m/z in consecutive scans
        L_x1 <- length(x1)
        if (L_x1 >= minNIonPair) {         # Minimum number of scans
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
              xMin <- which.min(abs(A[xuR, 1] - MZ_Int_ScN[i, 1]))
              xR <- xuR[-xMin[1]]
              A[xR, 4]
            }))
            x1 <- setdiff(x1, xRemove)
          }
          ##
          if (length(x1) >= minNIonPair) {
            Counter <- Counter + 1
            indexXIC[[Counter]] <- x1
            mzxic[Counter] <- MZ_Int_ScN[i, 1]
            MZ_Int_ScN[x1, ] <- 0
          }
        }
      }
    }
  } else {
    for (i in 1:L_MinHeight) {
      if (MZ_Int_ScN[i, 1] != 0) {
        x1 <- which(abs(MZ_Int_ScN[, 1] - MZ_Int_ScN[i, 1]) <= massAccuracyXIC) # to cluster m/z in consecutive scans
        L_x1 <- length(x1)
        if (L_x1 >= minNIonPair) {         # Minimum number of scans
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
              xMin <- which.min(abs(A[xuR, 1] - MZ_Int_ScN[i, 1]))
              xR <- xuR[-xMin[1]]
              A[xR, 4]
            }))
            x1 <- setdiff(x1, xRemove)
          }
          ##
          Counter <- Counter + 1
          indexXIC[[Counter]] <- x1
          mzxic[Counter] <- MZ_Int_ScN[i, 1]
          MZ_Int_ScN[x1, ] <- 0
        }
      }
    }
  }
  indexXIC <- indexXIC[1:Counter]
  mzxic <- mzxic[1:Counter]
  names(indexXIC) <- round(mzxic, 5)
  return(indexXIC)
}
