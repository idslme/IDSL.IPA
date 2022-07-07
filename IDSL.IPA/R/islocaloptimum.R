islocaloptimum <- function(y) {
  L_y <- length(y)
  Q <- rep(0, L_y)
  if (L_y > 2) {
    if (y[1] <= y[2]) {
      Q[1] <- -1
    }
    for (i in 2:(L_y - 1)) {
      if (y[i] > y[i - 1] & y[i] > y[i + 1]) {
        Q[i] <- +1
      } else if (y[i] <= y[i - 1] & y[i] <= y[i + 1]) {
        Q[i] <- -1
      }
    }
    if (y[L_y] <= y[L_y - 1]) {
      Q[L_y] <- -1
    }
  }
  return(Q)
}
