opendir <- function(dir){
  osType <- Sys.info()[['sysname']]
  if (osType == "Windows") {
    tryCatch(shell.exec(dir), error = function(e){})
  }
  if (osType == "Linux") {
    tryCatch(system(paste(Sys.getenv("R_BROWSER"), dir), ignore.stderr = TRUE), error = function(e){})
  }
  return()
}
