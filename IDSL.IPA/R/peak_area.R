peak_area <- function(x, y) {
  PA <- do.call(sum, lapply(1:(length(x) - 1), function(i) {
    (x[i + 1] - x[i])*(y[i + 1] + y[i])
  }))/2
  ##
  return(PA)
}
