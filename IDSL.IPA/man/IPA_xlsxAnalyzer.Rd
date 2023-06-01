\name{IPA_xlsxAnalyzer}
\alias{IPA_xlsxAnalyzer}
\title{
IPA xlsx Analyzer
}
\description{
This function processes the spreadsheet of the IPA parameters to ensure the parameter inputs are in agreement with the IPA requirements.
}
\usage{
IPA_xlsxAnalyzer(spreadsheet)
}
\arguments{
  \item{spreadsheet}{
IPA spreadsheet
}
}
\value{
This function returns the IPA parameters to feed the IPA_Workflow, IPA_CompoundsAnnotation, IPA_GapFiller, IPA_PeakAlignment, IPA_PeakAnalyzer, and IPA_PeaklistAnnotation functions.
}
\examples{
\donttest{
s_path <- system.file("extdata", package = "IDSL.IPA")
SSh1 <- paste0(s_path, "/IPA_parameters.xlsx")
temp_wd <- tempdir()
temp_wd_zip <- paste0(temp_wd,"/idsl_ipa_test_files.zip")
spreadsheet <- readxl::read_xlsx(SSh1)
PARAM = cbind(spreadsheet[, 2], spreadsheet[, 4])
tryCatch({download.file(paste0("https://github.com/idslme/IDSL.IPA/blob/main/",
                     "IPA_educational_files/idsl_ipa_test_files.zip?raw=true"),
              destfile = temp_wd_zip)
unzip(temp_wd_zip, exdir = temp_wd)
pass_download <- TRUE},
error = function(e) {pass_download <- FALSE},
warning = function(w) {pass_download <- FALSE})
if (pass_download) {
  PARAM[7, 2] <- temp_wd
  PARAM[10, 2] <- temp_wd # output data location
  PARAM[44, 2] <- s_path # reference file location
  PARAM <- IDSL.IPA::IPA_xlsxAnalyzer(PARAM)
}
}
}
