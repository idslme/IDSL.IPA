snr_rms <- function(int, baseline, gauge) {
  SNR <- 0
  int <- sort(int - baseline)
  Signal <- max(int)
  xN <- which(int/Signal < (1 - gauge))
  NVec <- int[xN]
  Noise <- sqrt(sum(NVec^2)/length(NVec)) # root mean square
  if (is.numeric(Noise)) {
    if (Noise != 0) {
      SNR <- Signal/Noise
    } else {
      SNR <- Inf
    }
  }
  return(SNR)
}
