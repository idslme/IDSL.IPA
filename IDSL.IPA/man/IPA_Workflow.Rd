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
\seealso{
\url{https://ipa.idsl.me/home}
}
\examples{
\donttest{
library(IDSL.IPA)
s_path <- system.file("extdata", package = "IDSL.IPA")
SSh1 <- paste0(s_path,"/IPA_parameters.xlsx")
## To see the results, use a known folder instead of the `tempdir()` command
temp_wd <- tempdir()
temp_wd_zip <- paste0(temp_wd, "/idsl_ipa_test_files.zip")
spreadsheet <- readxl::read_xlsx(SSh1)
download.file(paste0("https://github.com/idslme/IDSL.IPA/blob/main/",
"IPA_educational_files/idsl_ipa_test_files.zip?raw=true"),
destfile = temp_wd_zip, mode = "wb")
unzip(temp_wd_zip, exdir = temp_wd)
spreadsheet[7, 4] <- temp_wd
spreadsheet[43, 4] <- s_path
spreadsheet[10, 4] <- temp_wd
IPA_workflow(spreadsheet)
}
}
