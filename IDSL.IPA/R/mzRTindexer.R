mzRTindexer <- function(MZvec, RTvec, MZref, RTref, MZtolerance, RTtolerance) {
  ##
  indexRefVec <- which((abs(MZvec - MZref) <= MZtolerance) & (abs(RTvec - RTref) <= RTtolerance))
  ##
  l_indexRefVec <- length(indexRefVec)
  if (l_indexRefVec > 0) {
    if (l_indexRefVec > 1) {
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
