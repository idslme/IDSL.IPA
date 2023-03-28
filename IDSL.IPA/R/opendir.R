opendir <- function(dir) {
  osType <- Sys.info()[['sysname']]
  if (osType == "Windows") {
    tryCatch(shell.exec(dir), error = function(e) {NULL})
  } else {
    tryCatch(system(paste(Sys.getenv("R_BROWSER"), dir), ignore.stderr = TRUE), error = function(e) {NULL})
  }
}
