IPA_IonPairing <- function(spectraList, minSpectraNoiseLevel, massAccuracyIonPair = 0.015, ionMassDifference = 1.003354835336) {
  ##
  if (minSpectraNoiseLevel <= 0) {
    minSpectraNoiseLevel <- 1e-16 # This condition must be here to avoid interference
  }
  ##
  NumScans <- length(spectraList)
  ##
  spectraScan <- do.call(rbind, lapply(1:NumScans, function(t) {
    Spec <- spectraList[[t]]
    nSpec <- floor(nrow(Spec)/2)
    if (nSpec > 0) {
      ##
      Spec13C <- Spec[, 1] - ionMassDifference
      ##
      orderIntensity <- order(Spec[, 2], decreasing = TRUE) # To sort spectra list by intensity
      ##
      jSpectraScan <- matrix(rep(0, nSpec*5), ncol = 5)
      jCounter <- 0
      ##
      for (j in orderIntensity) {
        if (Spec[j, 2] >= minSpectraNoiseLevel) {           # Intensity threshold in each scan
          x13C <- which(abs(Spec[j, 1] - Spec13C) <= massAccuracyIonPair)
          Lx13C <- length(x13C)
          if (Lx13C > 0) {
            if (Lx13C > 1) {
              x13CMin <- which.min(abs(Spec[j, 1] - Spec13C[x13C]))
              x13C <- x13C[x13CMin]
              ##
              if (length(x13C) > 1) {
                x13CMax <- which.max(Spec[x13C, 2])
                x13C <- x13C[x13CMax[1]]
              }
            }
            jCounter <- jCounter + 1
            jSpectraScan[jCounter, ] <- c(Spec[j, 1], Spec[j, 2], t, Spec[x13C, 1], Spec[x13C, 2])
            Spec[x13C, ] <- 0
            Spec13C[x13C] <- 0
          }
        }
        Spec[j, ] <- 0
        Spec13C[j] <- 0
      }
      ##
      if (jCounter > 0) {
        jSpectraScan <- jSpectraScan[1:jCounter, ]
      }
    }
  }))
  ##
  rownames(spectraScan) <- NULL
  spectraScan <- spectraScan[order(spectraScan[, 2], decreasing = TRUE), ]   # To sort spectraScan by intensity
  return(spectraScan)
}
