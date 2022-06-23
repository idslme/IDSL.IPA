spectraList_filtering <- function(spec_scan.xic, spectraList, rounding_digit) {
  ## To reduce the number of arrays in the spectraList variable
  spectraList.mat <- do.call(rbind, lapply(1:length(spectraList), function(i) {
    cbind(spectraList[[i]], rep(i, dim(spectraList[[i]])[1]))
  }))
  selectIndex <- round(spectraList.mat[,1], digits = rounding_digit) %in% unique(round(spec_scan.xic[ ,1], digits = rounding_digit))
  spectraList.mat.sb <- spectraList.mat[selectIndex,]
  spectraList.Reduced <- lapply(1:length(spectraList), function(t) {
    x_t <- which(spectraList.mat.sb[, 3] == t)
    if (length(x_t) > 0) {
      spectraList.mat.sb[x_t, 1:2]
    } else {
      c(0,0) # This must be here
    }
  })
  return(spectraList.Reduced)
}
