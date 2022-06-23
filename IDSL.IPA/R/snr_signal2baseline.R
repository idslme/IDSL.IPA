snr_signal2baseline <- function(int, baseline) {
  SNR <- -Inf
  x_apex <- which.max(int)[1]
  Max_sig <- int[x_apex]
  if (Max_sig > 0) {
    Median_baseline <- baseline[x_apex]
    if (Median_baseline == 0) {
      Median_baseline <- 1
      Max_sig <- Max_sig/length(int)
    }
    SNR <- Max_sig/Median_baseline
  }
  return(SNR)
}