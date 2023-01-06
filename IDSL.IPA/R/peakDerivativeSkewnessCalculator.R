peakDerivativeSkewnessCalculator <- function(rt, int) {
  der1_rt_int <- derivative5pointsStencil(rt, int, n = 1)
  der1_skewness <- abs(log(-max(der1_rt_int)/min(der1_rt_int)))
  return(der1_skewness)
}