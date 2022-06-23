peak_width <- function(rt, int, gauge) {
  PW <- 0
  L_na <- length(which(is.na(int)))
  if (L_na == 0) {
    H <- max(int)
    x_H <- which(int == H)[1]
    H_gauge <- H*gauge
    ## left side of the peak
    int1 <- int[1:x_H]
    if (length(int1) > 1) {
      if (min(int1) <= H_gauge) {
        ## right side of the peak
        N_Seg <- length(rt)
        int2 <- int[x_H:N_Seg]
        if (length(int2) > 1) {
          if (min(int2) <= H_gauge) {
            ##
            rt1 <- rt[1:x_H]
            W1 <- approx(int1, rt1, H_gauge, method = "linear", 0, 0, rule = 2, f = 0, ties = mean)
            rt2 <- rt[x_H:N_Seg]
            W2 <- approx(int2, rt2, H_gauge, method = "linear", 0, 0, rule = 2, f = 0, ties = mean)
            PW <- W2[[2]] - W1[[2]]
          }
        }
      }
    }
  }
  return(PW)
}
