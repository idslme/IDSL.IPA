IPA_IonPairing <- function(spectraList, int_threshold, mass_accuracy_isotope_pair, massDifferenceIsotopes = 1.003354835336) {
  ##
  if (int_threshold <= 0) {
    int_threshold <- 1e-16 # This condition must be here to avoid interference
  }
  ##
  NumScans <- length(spectraList)
  ##
  spec_scan <- do.call(rbind, lapply(1:NumScans, function(t) {
    Spec <- spectraList[[t]]
    L_Spec <- nrow(Spec)
    if (L_Spec > 0) {
      ##
      Spec13C <- Spec[, 1] - massDifferenceIsotopes
      ##
      orderIntensity <- order(Spec[, 2], decreasing = TRUE) # To sort spectra list according to their intensities
      ##
      spec_scan_j <- matrix(rep(0, L_Spec*5), ncol = 5)
      counter_j <- 0
      ##
      for (j in orderIntensity) {
        if (Spec[j, 2] >= int_threshold) {        # Intensity threshold in each scan
          x13C <- which(abs(Spec[j, 1] - Spec13C) <= mass_accuracy_isotope_pair)
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
            spec_scan_j[counter_j, ] <- c(Spec[j, 1], Spec[j, 2], t, Spec[x13C, 1], Spec[x13C, 2])
            Spec[x13C, ] <- 0
            Spec13C[x13C] <- 0
          }
        }
        Spec[j, ] <- 0
        Spec13C[j] <- 0
      }
      ##
      if (counter_j > 0) {
        spec_scan_j <- spec_scan_j[1:counter_j, ]
      }
    }
  }))
  ##
  rownames(spec_scan) <- NULL
  spec_scan <- spec_scan[order(spec_scan[, 2], decreasing = TRUE), ]   # Sort spec_scan rows by their intensity
  return(spec_scan)
}
