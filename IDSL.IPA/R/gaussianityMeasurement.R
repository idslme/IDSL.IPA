gaussianityMeasurement <- function(RT, Int, BL, gauge = 0.8) {
  GM <- 0
  Int <- Int - BL
  W0.5 <- peakWidthCalculator(RT, Int, 0.50)
  if (W0.5 > 0) {
    Int_max <- max(Int)
    Sig <- W0.5/(2*sqrt(2*log(2)))
    rtmax <- RT[which.max(Int)[1]]
    G <- Int_max*exp(-(RT - rtmax)^2/(2*Sig^2))
    x_g <- which(Int/Int_max >= (1 - gauge))
    Int_g <- Int[x_g]
    G_g <- G[x_g]
    GM <- cor(Int_g, G_g, method = "pearson")
  }
  return(GM)
}
