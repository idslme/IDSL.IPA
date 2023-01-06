mzRTindexer <- function(MZvec, RTvec, MZref, RTref, massAccuracy, RTtolerance) {
  ##
  measuredRTdiff <- abs(RTvec - RTref)
  measuredMassDiff <- abs(MZvec - MZref)
  indexRefVec <- which((measuredMassDiff <= massAccuracy) & (measuredRTdiff <= RTtolerance))
  ##
  LindexRefVec <- length(indexRefVec)
  if (LindexRefVec > 0) {
    if (LindexRefVec > 1) {
      ##
      if (max(measuredMassDiff[indexRefVec]) == 0) {
        xMin <- which.min(measuredRTdiff[indexRefVec])
      } else if (max(measuredRTdiff[indexRefVec]) == 0) {
        xMin <- which.min(measuredMassDiff[indexRefVec])
      } else {
        xMin <- which.min(measuredMassDiff[indexRefVec]*measuredRTdiff[indexRefVec])
      }
      ##
      indexRefVec <- indexRefVec[xMin[1]]
    }
    ##
  } else {
    indexRefVec <- NULL
  }
  ##
  return(indexRefVec)
}
