\name{IPA_TargetedAnalysis}
\alias{IPA_TargetedAnalysis}
\title{
IPA Targeted Analysis
}
\description{
This function plots extracted ion chromatogram (EIC) figures in the targeted mode.
}
\usage{
IPA_TargetedAnalysis(spreadsheet, mzCandidate, rtCandidate, exportEIC = TRUE,
exportTable = FALSE)
}
\arguments{
  \item{spreadsheet}{
a spreadsheet containing the parameters.
}
  \item{mzCandidate}{
a vector of candidate m/z values.
}
  \item{rtCandidate}{
a vector of candidate RT values.
}
  \item{exportEIC}{
TRUE by defult. To plot and save EICs.
}
  \item{exportTable}{
FALSE by defult. To return the whole peaklists for the m/z and RT vectors,
select TRUE.
}
}
\value{
This function saves extracted ion chromatograms in .png format in the "EICs" folder when "exportEIC = TRUE", and it save a table of peak properties when "exportTable = TRUE".
}
\examples{
\donttest{
s_path <- system.file("extdata", package = "IDSL.IPA")
SSh1 <- paste0(s_path,"/IPA_parameters.xlsx")
spreadsheet <- readxl::read_xlsx(SSh1, sheet = 'IPA_targeted')
temp_wd <- tempdir()
temp_wd_zip <- paste0(temp_wd,"/testfiles.zip")
download.file(
  "https://github.com/idslme/IDSL.IPA/raw/main/idsl_ipa_test_files.zip",
  destfile = temp_wd_zip)
unzip(temp_wd_zip, exdir = temp_wd)
spreadsheet[2, 4] <- temp_wd
spreadsheet[5, 4] <- temp_wd
mzCandidate <- c(53.01853, 61.00759)
rtCandidate <- c(0.951, 0.961)
IDSL.IPA::IPA_TargetedAnalysis(spreadsheet, mzCandidate, rtCandidate)
}
}
