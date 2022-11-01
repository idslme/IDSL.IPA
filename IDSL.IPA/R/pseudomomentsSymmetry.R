pseudomomentsSymmetry <- function(rt, int) {
  H <- max(int)
  x_H <- which(int == H[1])[1]
  nSpline <- length(rt)
  rt1 <- rt[1:x_H]
  int1 <- int[1:x_H]
  rt2 <- rt[x_H:nSpline]
  int2 <- int[x_H:nSpline]
  PeakSymmetry <- NaN
  Skewness <- NaN
  if (x_H >= 5 & (nSpline - x_H - 1) >= 5) {
    ## left side of the peak
    der2_int1 <- derivative5pointsStencil(rt1, int1, n = 2)[, 2]
    J <- islocalminimum(abs(der2_int1))
    x_j <- which(J == -1) + 2
    if (length(x_j) > 0) {
      if (length(x_j) > 1) {
        x_H_2 <- which.min(abs(int1 - H/2))
        x_j1 <- which.min(abs(x_H_2[1] - x_j))
        x_j <- x_j[x_j1[1]]
      }
      Hf <- int1[x_j]
      rtf <- rt1[x_j]
      t1 <- rtf - rt1[1]
      x_j1 <- which.min(abs(rt1 - rtf))
      a1 <- peakAreaCalculator(rt1[1:x_j1], int1[1:x_j1])
      if (a1 > 0) {
        a2 <- peakAreaCalculator(rt1[x_j1:x_H], int1[x_j1:x_H])
        if (a2 > 0) {
          t2 <- rt1[x_H] - rtf
          ## right side of the peak
          der2_int2 <- derivative5pointsStencil(rt2, int2, n = 2)[, 2]
          J <- islocalminimum(abs(der2_int2))
          x_j <- which(J == -1) + 2
          if (length(x_j) > 0) {
            if (length(x_j) > 1) {
              x_H_2 <- which.min(abs(int2 - H/2))
              x_j1 <- which.min(abs(x_H_2[1] - x_j))
              x_j <- x_j[x_j1[1]]
            }
            Hr <- int2[x_j]
            rtr <- rt2[x_j]
            x_j1 <- which.min(abs(rt2 - rtr))
            a3 <- peakAreaCalculator(rt2[1:x_j1], int2[1:x_j1])
            if (a3 > 0) {
              a4 <- peakAreaCalculator(rt2[x_j1:(nSpline - x_H + 1)], int2[x_j1:(nSpline - x_H + 1)])
              if (a4 > 0) {
                t3 <- rtr - rt2[1]
                ##
                m1 <- a1*(t2+a1/(1.5*Hf))
                m2 <- a2^2/(0.5*Hf + 1.5*H)
                m3 <- a3^2/(0.5*Hr + 1.5*H)
                m4 <- a4*(t3+a4/(1.5*Hr))
                PeakSymmetry <- sqrt((m1 + m2)/(m3 + m4))
                Skewness <- m3/m2
                ##
              }
            }
          }
        }
      }
    }
  }
  if (is.nan(PeakSymmetry)) {
    a12 <- peakAreaCalculator(rt1, int1)
    a34 <- peakAreaCalculator(rt2, int2)
    PeakSymmetry <- a12/a34
  }
  return(c(PeakSymmetry, Skewness))
}
