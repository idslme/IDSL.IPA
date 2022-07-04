loadRdata <- function(fileName) {
  #loads an Rdata file, and returns it
  fileName <- gsub("\\", "/", fileName, fixed = TRUE)
  load(fileName)
  get(ls()[ls() != "fileName"])
}
