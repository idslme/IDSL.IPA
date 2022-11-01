IPA_IonPairing <- function(spectraList, minSpectraNoiseLevel, massAccuracyIonPair, ionMassDifference = 1.003354835336) {
  ##
  if (minSpectraNoiseLevel <= 0) {
    minSpectraNoiseLevel <- 1e-16 # This condition must be here to avoid interference
  }
  ##
  NumScans <- length(spectraList)
  ##
  spectraScan <- do.call(rbind, lapply(1:NumScans, function(t) {
    Spec <- spectraList[[t]]
    L_Spec <- nrow(Spec)
    if (L_Spec > 0) {
      ##
      Spec13C <- Spec[, 1] - ionMassDifference
      ##
      orderIntensity <- order(Spec[, 2], decreasing = TRUE) # To sort spectra list according to their intensities
      ##
      spectraScan_j <- matrix(rep(0, L_Spec*5), ncol = 5)
      counter_j <- 0
      ##
      for (j in orderIntensity) {
        if (Spec[j, 2] >= minSpectraNoiseLevel) {        # Intensity threshold in each scan
          x13C <- which(abs(Spec[j, 1] - Spec13C) <= massAccuracyIonPair)
          L_x13C <- length(x13C)
          if (L_x13C > 0) {
            if (L_x13C > 1) {
              x13CMin <- which.min(abs(Spec[j, 1] - Spec13C[x13C]))
              x13C <- x13C[x13CMin]
              ##
              if (length(x13C) > 1) {
                x13CMax <- which.max(Spec[x13C, 2])
                x13C <- x13C[x13CMax[1]]
              }
            }
            counter_j <- counter_j + 1
            spectraScan_j[counter_j, ] <- c(Spec[j, 1], Spec[j, 2], t, Spec[x13C, 1], Spec[x13C, 2])
            Spec[x13C, ] <- 0
            Spec13C[x13C] <- 0
          }
        }
        Spec[j, ] <- 0
        Spec13C[j] <- 0
      }
      ##
      if (counter_j > 0) {
        spectraScan_j <- spectraScan_j[1:counter_j, ]
      }
    }
  }))
  ##
  rownames(spectraScan) <- NULL
  spectraScan <- spectraScan[order(spectraScan[, 2], decreasing = TRUE), ]   # Sort spectraScan rows by their intensity
  return(spectraScan)
}
