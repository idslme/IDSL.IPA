IPA_message <- function(messageQuote, failedMessage = TRUE) {
  ##
  if (failedMessage) {
    col <- 91 ## Red color message
  } else {
    col <- 92 ## Green color message
  }
  ##
  message(paste0("\033[0;", col, "m", messageQuote, "\033[0m"))
}
