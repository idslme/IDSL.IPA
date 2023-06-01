\name{IPA_workflow}
\alias{IPA_workflow}
\alias{IPA_Workflow}
\title{
IPA Workflow
}
\description{
This function executes the IPA workflow in order.
}
\usage{
IPA_workflow(spreadsheet)
IPA_Workflow(spreadsheet)
}
\arguments{
  \item{spreadsheet}{
IPA spreadsheet
}
}
\value{
This function organizes the IPA file processing for a better performance using the template spreadsheet.
}
\examples{
\donttest{
s_path <- system.file("extdata", package = "IDSL.IPA")
SSh1 <- paste0(s_path,"/IPA_parameters.xlsx")
## To see the results, use a known folder instead of the `tempdir()` command
temp_wd <- tempdir()
temp_wd_zip <- paste0(temp_wd, "/idsl_ipa_test_files.zip")
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
  PARAM[44, 2] <- s_path
  PARAM[10, 2] <- temp_wd
  IPA_workflow(PARAM)
}
}
}
