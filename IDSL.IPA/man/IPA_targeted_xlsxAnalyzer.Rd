\name{IPA_targeted_xlsxAnalyzer}
\alias{IPA_targeted_xlsxAnalyzer}
\title{
IPA Targeted xlsxAnalyzer
}
\description{
This function processes the spreadsheet of the IPA parameters to ensure the parameter inputs are in agreement with the `IPA_targeted` function.
}
\usage{
IPA_targeted_xlsxAnalyzer(spreadsheet)
}
\arguments{
  \item{spreadsheet}{
contains the IPA parameters.
}
}
\value{
`PARAM_targeted` which is the IPA parameters to feed the `IPA_targeted` function.
}
\examples{
\donttest{
## To generate the IPA spreadsheet parameters
s_path <- system.file("extdata", package = "IDSL.IPA")
SSh1 <- paste0(s_path,"/IPA_parameters.xlsx")
spreadsheet <- readxl::read_xlsx(SSh1, sheet = 'IPA_targeted')
PARAM_targeted = cbind(spreadsheet[, 2], spreadsheet[, 4])
temp_wd <- tempdir()
temp_wd_zip <- paste0(temp_wd,"/idsl_ipa_test_files.zip")
tryCatch({download.file(paste0("https://github.com/idslme/IDSL.IPA/blob/main/",
                               "IPA_educational_files/idsl_ipa_test_files.zip?raw=true"),
                        destfile = temp_wd_zip)
  unzip(temp_wd_zip, exdir = temp_wd)
  pass_download <- TRUE},
  error = function(e) {pass_download <- FALSE},
  warning = function(w) {pass_download <- FALSE})
if (pass_download) {
  PARAM_targeted[3, 2] <- temp_wd
  PARAM_targeted[7, 2] <- temp_wd
  PARAM_targeted[8, 2] <- "53.01853, 61.00759"
  PARAM_targeted[9, 2] <- "0.951, 0.961"
  ##
  PARAM_targeted <- IPA_targeted_xlsxAnalyzer(PARAM_targeted)
}
}
}
