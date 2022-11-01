USPtailingFactor <- function(rt, int) {
  gauge <- 0.05
  H <- max(int)
  x_H <- which(int == H)[1]
  tR <- rt[x_H]
  ## left side of the peak
  rt1 <- rt[1:x_H]
  int1 <- int[1:x_H]
  W1 <- approx(int1, rt1, H*gauge, method = "linear", 0, 0, rule = 2, f = 0, ties = mean)
  f0.05 <- tR - W1[[2]]
  ## left side of the peak
  W0.05 <- peakWidthCalculator(rt, int, gauge)
  ##
  usp_tf <- W0.05/(2*f0.05)
  return(usp_tf)
}
