IPA_baselineDeveloper <- function(segment, int) {
  n_BL <- sort(unique(as.vector(segment)))
  SCAN <- 1:length(int)
  BL_Base <- rep(0, length(SCAN))
  BL_Base[n_BL] <- int[n_BL] + 1e-16
  counter_seq <- 0
  seq_BL <- matrix(rep(0, length(n_BL)*2), ncol = 2)
  x_seq <- c(n_BL[1], 0)
  for (i in 2:length(n_BL)) {
    if (median(int[n_BL[i - 1]:n_BL[i]]) != 0) {
      x_seq[2] <- n_BL[i]
    }
    if (int[n_BL[i]] == 0) {
      if (x_seq[2] != 0) {
        counter_seq <- counter_seq + 1
        seq_BL[counter_seq, ] <- x_seq
        x_seq <- c(0, 0)
      }
      x_seq[1] <- n_BL[i]
    }
  }
  if (x_seq[2] != 0) {
    counter_seq <- counter_seq + 1
    seq_BL[counter_seq, ] <- x_seq
  }
  seq_BL <- seq_BL[1:counter_seq, ]
  seq_BL <- matrix(seq_BL, ncol = 2)
  ##
  for (i in 1:nrow(seq_BL)) {
    SEG00 <- seq_BL[i, 1]:seq_BL[i, 2]
    L00 <- seq_BL[i, 2] - seq_BL[i, 1] + 1
    BL_Seg <- BL_Base[SEG00]
    ## To remove outliers using a Median Absolute Deviation (MAD) method
    # Moving forward to detect local outliers
    x0 <- which(BL_Seg >= 1e-16)
    if (length(x0) > 2) {
      for (k in 2:(length(x0) - 1)) {
        x0k <- x0[(k - 1):(k + 1)]
        MAD <- 1.482602*median(abs(BL_Seg[x0k] - median(BL_Seg[x0k])))  # -1/(sqrt(2)*erfcinv(3/2))=1.482602
        x_BL <- which(BL_Seg[x0k] > (median(BL_Seg[x0k]) + 3*MAD))
        if (length(x_BL) >= 1) {
          x1_BL <- x0k[x_BL]
          BL_Seg[x1_BL] <- (min(BL_Seg[x0k]) + max(BL_Seg[x0k]))/2
        }
      }
    }
    # Moving backward to detect local outliers
    x0 <- which(BL_Seg >= 1e-16)
    if (length(x0) > 2) {
      for (k in (length(x0) - 1):2) {
        x0k <- x0[(k - 1):(k + 1)]
        MAD <- 1.482602*median(abs(BL_Seg[x0k] - median(BL_Seg[x0k])))  # -1/(sqrt(2)*erfcinv(3/2))=1.482602
        x_BL <- which(BL_Seg[x0k] > (median(BL_Seg[x0k]) + 3*MAD))
        if (length(x_BL) >= 1) {
          x1_BL <- x0k[x_BL]
          BL_Seg[x1_BL] <- (min(BL_Seg[x0k]) + max(BL_Seg[x0k]))/2
        }
      }
    }
    ##
    x_non0 <- which(BL_Seg > 1e-16)
    if (length(x_non0) > 0) {
      BL_Seg[1] <- BL_Seg[x_non0[1]]
      BL_Seg[L00] <- BL_Seg[x_non0[length(x_non0)]]
    }
    BL_Base[SEG00] <- BL_Seg
    #
    x00 <- which(BL_Base[SEG00] > 1e-16)
    if (length(x00) > 0) {
      x00 <- x00 + SEG00[1] - 1
      if (x00[1] != SEG00[1]) {
        x00 <- c(SEG00[1], x00)
      }
      if (x00[length(x00)] != seq_BL[i,2]) {
        x00 <- c(x00, seq_BL[i,2])
      }
      W <- approx(SCAN[x00], BL_Base[x00], SCAN[SEG00], method = "linear", 0, 0, rule = 2, f = 0, ties = mean)
      BL_Base[SEG00] <- W[[2]]
    }
  }
  x_BL_Neg <- which(BL_Base > int)
  if (length(x_BL_Neg) > 0) {
    BL_Base[x_BL_Neg] <- int[x_BL_Neg]
  }
  return(BL_Base)
}
