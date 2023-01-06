loadRdata <- function(fileName) {
  #loads an Rdata file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
