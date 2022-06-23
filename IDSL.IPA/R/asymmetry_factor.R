asymmetry_factor <- function(rt, int) {
  gauge <- 0.1
  H <- max(int)
  x_H <- which(int == H)[1]
  tR <- rt[x_H]
  ## left side of the peak
  rt1 <- rt[1:x_H]
  int1 <- int[1:x_H]
  W1 <- approx(int1, rt1, H*gauge, method = "linear", 0, 0, rule = 2, f = 0, ties = mean)
  a <- tR - W1[[2]]
  ## right side of the peak
  N_Seg <- length(rt)
  rt2 <- rt[x_H:N_Seg]
  int2 <- int[x_H:N_Seg]
  W2 <- approx(int2, rt2, H*gauge, method = "linear", 0, 0, rule = 2, f = 0, ties = mean)
  b <- W2[[2]] - tR
  return(b/a)
}