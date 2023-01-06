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
## To generate the IPA spreadsheet parameters
s_path <- system.file("extdata", package = "IDSL.IPA")
SSh1 <- paste0(s_path,"/IPA_parameters.xlsx")
spreadsheet <- readxl::read_xlsx(SSh1, sheet = 'IPA_targeted')
temp_wd <- tempdir()
temp_wd_zip <- paste0(temp_wd,"/idsl_ipa_test_files.zip")
download.file(paste0("https://github.com/idslme/IDSL.IPA/blob/main/",
"IPA_educational_files/idsl_ipa_test_files.zip?raw=true"),
destfile = temp_wd_zip, mode = "wb")
unzip(temp_wd_zip, exdir = temp_wd)
spreadsheet[3, 4] <- temp_wd
spreadsheet[7, 4] <- temp_wd
spreadsheet[8, 4] <- "53.01853, 61.00759"
spreadsheet[9, 4] <- "0.951, 0.961"
##
PARAM_targeted <- IPA_targeted_xlsxAnalyzer(spreadsheet)
}
