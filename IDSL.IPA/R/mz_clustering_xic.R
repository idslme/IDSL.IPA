mz_clustering_xic <- function(spec_scan, mass_accuracy_xic, min_peak_height, min_nIsoPair) {
  MZ_Int_ScN <- spec_scan[, 1:3]
  index_xic <- vector(mode = "list", floor(nrow(MZ_Int_ScN)/max(c(1, min_nIsoPair))))
  Counter <- 0
  L_MinHeight <- length(which(MZ_Int_ScN[, 2] >= min_peak_height))
  i <- 1
  if (min_nIsoPair > 1) {
    repeat {
      if (i <= L_MinHeight) {
        x1 <- which(abs(MZ_Int_ScN[, 1] - MZ_Int_ScN[i, 1]) <= mass_accuracy_xic) # to cluster m/z in consecutive scans
        L_x1 <- length(x1)
        if (L_x1 >= min_nIsoPair) {         # Minimum number of scans
          A <- MZ_Int_ScN[x1, ]
          ## To remove repeated scans
          if (L_x1 != length(unique(A[, 3]))) {
            x2 <- table(A[, 3])
            x3 <- which(x2 > 1)
            R_ScN <- as.numeric(names(x3))  # Repeated scan numbers
            A <- cbind(A, x1)
            x7 <- sapply(1:length(R_ScN), function(j) {
              x4 <- which(A[, 3] == R_ScN[j])
              x5 <- which.min(abs(A[x4, 1] - MZ_Int_ScN[i, 1]))[1]
              x6 <- x4[-x5]
              A[x6, 4]
            })
            x8 <- unlist(x7)
            x1 <- setdiff(x1, x8)
          }
          ##
          if (length(x1) >= min_nIsoPair) {
            Counter <- Counter + 1
            index_xic[[Counter]] <- x1
            MZ_Int_ScN[x1, 1] <- 0
          }
        }
        ##
        i <- i + which(MZ_Int_ScN[(i + 1):L_MinHeight, 1] != 0)[1]
        if (is.na(i)) {
          break
        }
      } else {
        break
      }
    }
  } else {
    repeat {
      if (i <= L_MinHeight) {
        x1 <- which(abs(MZ_Int_ScN[, 1] - MZ_Int_ScN[i, 1]) <= mass_accuracy_xic) # to cluster m/z in consecutive scans
        L_x1 <- length(x1)
        if (L_x1 >= min_nIsoPair) {         # Minimum number of scans
          A <- matrix(MZ_Int_ScN[x1, ], ncol = 3)
          ## To remove repeated scans
          if (L_x1 != length(unique(A[, 3]))) {
            x2 <- table(A[, 3])
            x3 <- which(x2 > 1)
            R_ScN <- as.numeric(names(x3))  # Repeated scan numbers
            A <- cbind(A, x1)
            x7 <- sapply(1:length(R_ScN), function(j) {
              x4 <- which(A[, 3] == R_ScN[j])
              x5 <- which.min(abs(A[x4, 1] - MZ_Int_ScN[i, 1]))[1]
              x6 <- x4[-x5]
              A[x6, 4]
            })
            x8 <- unlist(x7)
            x1 <- setdiff(x1, x8)
          }
          ##
          Counter <- Counter + 1
          index_xic[[Counter]] <- x1
          MZ_Int_ScN[x1, 1] <- 0
        }
        ##
        i <- i + which(MZ_Int_ScN[(i + 1):L_MinHeight, 1] != 0)[1]
        if (is.na(i)) {
          break
        }
      } else {
        break
      }
    }
  }
  index_xic <- index_xic[1:Counter]
  return(index_xic)
}
