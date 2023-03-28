IPA_spectraListAggregator <- function(spectraList) {
  ## To aggregate spectraList variable for faster data parsing
  nSpectraList <- length(spectraList)
  ##
  spectraListMatrix <- do.call(rbind, lapply(1:nSpectraList, function(t) {
    cbind(rep(t, dim(spectraList[[t]])[1]), spectraList[[t]])
  }))
  spectraListMatrix <- spectraListMatrix[(spectraListMatrix[, 3] > 0), ]
  nSpectraListMatrix <- dim(spectraListMatrix)[1]
  ##
  IDroundMass <- cbind(round(spectraListMatrix[, 2], digits = 2), seq(1, nSpectraListMatrix, 1))
  IDroundMass <- IDroundMass[order(IDroundMass[, 1], decreasing = FALSE), ]
  xDiff <- c(0, which(diff(IDroundMass[, 1]) > 0), nSpectraListMatrix)
  nDiff <- length(xDiff) - 1
  aggregatedSpectraList <- lapply(1:nDiff, function(i) {
    IDroundMass[(xDiff[i] + 1):xDiff[i + 1], 2]
  })
  names(aggregatedSpectraList) <- as.character(IDroundMass[(xDiff[1:nDiff] + 1), 1])
  ##
  aggregatedSpectraListMatrixList <- list(aggregatedSpectraList, spectraListMatrix)
  names(aggregatedSpectraListMatrixList) <- c("aggregatedSpectraList", "spectraListMatrix")
  ##
  return(aggregatedSpectraListMatrixList)
}
