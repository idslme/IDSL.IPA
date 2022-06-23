fronting_tailing_resolver <- function(segment, int, max_space, peak_resolving_power) {
  R_segment <- length(segment)/2
  if (R_segment > 1) {
    SCAN <- 1:length(int)
    R <- do.call(rbind, lapply(1:R_segment, function(i) {
      x_max <- which.max(int[segment[i, 1]:segment[i, 2]])
      x_max <- segment[i, 1] + x_max[1] - 1
      c(segment[i, 1], x_max, segment[i, 2])
    }))
    i <- 1
    while (i != (length(R)/3)) {
      if ((R[i + 1, 1] - R[i, 3]) <= max_space) {
        r2 <- SCAN[R[i, 2]]
        r3 <- SCAN[R[i, 3]]
        r5 <- SCAN[R[i + 1, 2]]
        l2 <- int[R[i, 2]]
        l5 <- int[R[i + 1, 2]]
        if (l2 > l5) {
          l3 <- int[R[i + 1, 1]]
          l4 <- l5 - l3
          z <- (l2 - l5)/(r5 - r2)*(r5 - r3)
        } else if (l2 < l5) {
          l3 <- int[R[i, 3]]
          l4 <- l2 - l3
          z <- (l5 - l2)/(r5 - r2)*(r3 - r2)
        } else if (l2 == l5) {
          l3 <- int[R[i, 3]]
          l4 <- l5 - l3
          z <- 1e-5
        }
        if (z == 0 | l4/z <= peak_resolving_power) {
          R[i, 3] <- R[i + 1, 3]
          if (l2 < l5) {
            R[i, 2] <- R[i + 1, 2]
          }
          R <- R[-(i + 1), ]
          i <- i - 1
        }
      }
      i <- i + 1
    }
    R <- matrix(R, ncol = 3)
    segment <- matrix(cbind(R[, 1], R[, 3]), ncol = 2)
  }
  return(segment)
}
