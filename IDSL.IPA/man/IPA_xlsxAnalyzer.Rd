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
temp_wd_zip <- paste0(temp_wd,"/testfiles.zip")
spreadsheet <- readxl::read_xlsx(SSh1)
download.file(
"https://github.com/idslme/IDSL.IPA/raw/main/idsl_ipa_test_files.zip",
destfile = temp_wd_zip)
unzip(temp_wd_zip, exdir = temp_wd)
spreadsheet[7, 4] <- temp_wd
spreadsheet[40, 4] <- s_path # reference file location
spreadsheet[10, 4] <- temp_wd # output data location
PARAM <- IDSL.IPA::IPA_xlsxAnalyzer(spreadsheet)
}
}
