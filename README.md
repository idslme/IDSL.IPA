# IDSL.IPA<img src='IPA_educational_files/Figures/IDSL.IPA-logo.PNG' width="250px" align="right" />

<!-- badges: start -->
[![Maintainer](https://img.shields.io/badge/maintainer-Sadjad_Fakouri_Baygi-blue)](https://github.com/sajfb)
[![CRAN status](https://www.r-pkg.org/badges/version/IDSL.IPA)](https://cran.r-project.org/package=IDSL.IPA)
![](http://cranlogs.r-pkg.org/badges/IDSL.IPA?color=orange)
![](http://cranlogs.r-pkg.org/badges/grand-total/IDSL.IPA?color=brightgreen)
[![Dependencies](https://tinyverse.netlify.com/badge/IDSL.IPA)](https://cran.r-project.org/package=IDSL.IPA)
<!-- badges: end -->

[**Intrinsic Peak Analysis (IPA)**](https://pubs.acs.org/doi/10.1021/acs.jproteome.2c00120) by the Integrated Data Science Laboratory for Metabolomics and Exposomics (IDSL.ME) is an R pipeline that extracts peaks for organic small molecules from untargeted LC/HRMS data in population scale projects. 

	install.packages("IDSL.IPA")

## <img src='IPA_educational_files/Figures/IDSL.IPA-TOC_Art.png' align="right" />

## Workflow
To process your mass spectrometry data (**mzXML**, **mzML**, **netCDF**), download the [IPA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.IPA/main/IPA_parameters.xlsx) and fill out the parameters accordingly and then use this spreadsheet as the input for the IDSL.IPA workflow:

	IPA_Workflow("Address of the IPA parameter spreadsheet")

Visit https://ipa.idsl.me/ for the detailed documentation and tutorial.