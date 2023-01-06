IPA_aggregate <- function(idVec, variableVec, indexVec, targetVar) {
  if (length(idVec) > length(unique(idVec))) {
    A <- cbind(idVec, variableVec, indexVec)
    A <- A[order(A[, 1], decreasing = FALSE), ]
    xDiff0 <- which(diff(A[, 1]) == 0)
    xu <- c(0, which(diff(xDiff0) > 1), length(xDiff0))
    xRemove <- do.call(c, lapply(1:(length(xu) - 1), function(j) {
      xuDiff <- xDiff0[(xu[j] + 1):xu[j + 1]]
      xuR <- c(xuDiff, (xuDiff[length(xuDiff)] + 1))  # Repeated idVec numbers
      xMin <- which.min(abs(A[xuR, 2] - targetVar))
      xR <- xuR[-xMin[1]]
      A[xR, 3]
    }))
    indexVec <- setdiff(indexVec, xRemove)
  }
  return(indexVec)
}
