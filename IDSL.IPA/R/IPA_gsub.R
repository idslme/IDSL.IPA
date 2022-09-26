IPA_gsub <- function(pattern, replacement, x, ignore.case = FALSE,
                     perl = FALSE, fixed = FALSE, useBytes = FALSE)  {
  for (i in pattern) {
    x <- gsub(pattern = i, replacement, x, ignore.case, perl, fixed, useBytes)
  }
  ##
  return(x)
}
