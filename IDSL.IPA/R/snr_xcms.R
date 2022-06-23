snr_xcms <- function(int) {
  int <- sort(int)
  L_int <- length(int)
  L_int95 <- floor(L_int*0.95)
  L_int5 <- floor(L_int*0.05) + 1
  if (L_int5 == 1) {
    L_int5 <- L_int5 + 1
  }
  x90 <- L_int5:L_int95
  int90 <- int[x90]
  BaseLine <- mean(int90)
  NoiseLevel <- sd(int90)
  SNR <- (max(int) - BaseLine)/NoiseLevel
  return(SNR)
}
