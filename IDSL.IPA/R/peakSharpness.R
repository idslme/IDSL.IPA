peakSharpness <- function(int) {
  SH <- -Inf
  x_0 <- which(int > 0)
  int <- int[x_0]
  L_int <- length(int)
  if (L_int >= 3) {
    x_H <- which.max(int)[1]
    if (x_H >= 2) {
      SH1 <- do.call(c, lapply(2:x_H, function(i) {
        (int[i] - int[i - 1])/int[i - 1]
      }))
      SH2 <- do.call(c, lapply(x_H:(L_int - 1), function(i) {
        (int[i] - int[i + 1])/int[i + 1]
      }))
      SH <- sum(SH1) + sum(SH2)
    }
  }
  return(SH)
}
