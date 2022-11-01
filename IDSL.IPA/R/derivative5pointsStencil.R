derivative5pointsStencil <- function(x, y, n) {
  der_5points <- NULL
  L_y <- length(y)
  if (L_y >= 5) {
    K <- matrix(c(-1, -1, 1, 1, 8, 16, -2, -4, 0, -30, 0, 6, -8, 16, 2, -4, 1, -1, -1, 1, 12, 12, 2, 1, 1, 2, 3, 4) , nrow = 4)
    der.x_y <- do.call(c, lapply(3:(L_y - 2), function(i) {
      h <- mean(diff(x[(i - 2):(i + 2)]))
      (K[n, 1]*y[i + 2] + K[n, 2]*y[i + 1] + K[n, 3]*y[i] + K[n, 4]*y[i - 1] + K[n, 5]*y[i - 2])/(K[n, 6]*(h^(K[n, 7])))
    }))
    x_remove <- c(1, 2, (L_y - 1), L_y)
    x <- x[-x_remove]
    der_5points <- cbind(x, der.x_y)
  }
  return(der_5points)
}
