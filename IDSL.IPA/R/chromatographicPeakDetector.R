chromatographicPeakDetector <- function(int) {
  segment <- NULL
  Q <- islocalminimum(int)
  x_Q0 <- which(Q == 0)
  L_x_Q0 <- length(x_Q0)
  if (L_x_Q0 > 0) {
    segment <- matrix(rep(0, 2*(L_x_Q0 + 4)), ncol = 2)
    counter_isomer <- 0
    if (x_Q0[1] == 1) {
      x_Seg <- c(1, 0)
    } else {
      x_Seg <- c(0, 0)
    }
    if (L_x_Q0 != 1) {
      for (i in 1:(L_x_Q0 - 1)) {
        if (x_Seg[1] == 0) {
          x_Seg[1] <- x_Q0[i] - 1
        }
        if ((x_Q0[i + 1] - x_Q0[i]) != 1) {
          x_Seg[2] <- x_Q0[i] + 1
          counter_isomer <- counter_isomer + 1
          segment[counter_isomer, ] <- x_Seg
          x_Seg <- c(0, 0)
        }
      }
    } else {
      i <- 0
    }
    x_Seg[2] <- x_Q0[i + 1] + 1
    SZC <- length(int)
    if (x_Seg[2] > SZC) {
      x_Seg[2] <- SZC
    }
    if (x_Seg[1] == 0) {
      x_Seg[1] <- x_Q0[i + 1] - 1
    }
    counter_isomer <- counter_isomer + 1
    segment[counter_isomer, ] <- x_Seg
    segment <- segment[1:counter_isomer, ]
    segment <- matrix(segment, ncol = 2)
  }
  return(segment)
}
