\name{IPA_Workflow}
\alias{IPA_Workflow}
\title{
IPA Workflow
}
\description{
This function executes the IPA workflow in order.
}
\usage{
IPA_Workflow(spreadsheet)
}
\arguments{
  \item{spreadsheet}{
IPA spreadsheet
}
}
\value{
This function orgonizes the IPA file processing for a better perofrmance using the template spreadsheet.
}
\seealso{
\url{https://ipa.idsl.me/home}
}
\examples{
\donttest{
library(IDSL.IPA)
s_path <- system.file("extdata", package = "IDSL.IPA")
SSh1 <- paste0(s_path,"/IPA_parameters.xlsx")
temp_wd <- tempdir()
temp_wd_zip <- paste0(temp_wd,"/idsl_ipa_test_files.zip")
spreadsheet <- readxl::read_xlsx(SSh1)
download.file(paste0("https://github.com/idslme/IDSL.IPA/blob/main/",
"IPA_educational_files/idsl_ipa_test_files.zip?raw=true"),
destfile = temp_wd_zip)
unzip(temp_wd_zip, exdir = temp_wd)
spreadsheet[7, 4] <- temp_wd
spreadsheet[40, 4] <- s_path
spreadsheet[10, 4] <- temp_wd
IPA_Workflow(spreadsheet)
}
}
