mzRTindexer <- function(MZvec, RTvec, MZref, RTref, massAccuracy, RTtolerance) {
  ##
  indexRefVec <- which((abs(MZvec - MZref) <= massAccuracy) & (abs(RTvec - RTref) <= RTtolerance))
  ##
  LindexRefVec <- length(indexRefVec)
  if (LindexRefVec > 0) {
    if (LindexRefVec > 1) {
      ##
      measuredRTerror <- abs(RTvec[indexRefVec] - RTref)
      measuredMassError <- abs(MZvec[indexRefVec] - MZref)
      ##
      if (max(measuredMassError) == 0) {
        x_min <- which.min(measuredRTerror)
      } else if (max(measuredRTerror) == 0) {
        x_min <- which.min(measuredMassError)
      } else {
        x_min <- which.min(measuredMassError*measuredRTerror)
      }
      ##
      indexRefVec <- indexRefVec[x_min[1]]
    }
    ##
  } else {
    indexRefVec <- NULL
  }
  ##
  return(indexRefVec)
}
