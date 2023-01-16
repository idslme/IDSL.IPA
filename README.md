# IDSL.IPA<img src='IPA_educational_files/Figures/IDSL.IPA-logo.PNG' width="250px" align="right" />

<!-- badges: start -->
[![Maintainer](https://img.shields.io/badge/maintainer-Sadjad_Fakouri_Baygi-blue)](https://github.com/sajfb)
[![CRAN status](https://www.r-pkg.org/badges/version/IDSL.IPA)](https://cran.r-project.org/package=IDSL.IPA)
![](http://cranlogs.r-pkg.org/badges/IDSL.IPA?color=orange)
![](http://cranlogs.r-pkg.org/badges/grand-total/IDSL.IPA?color=brightgreen)
[![Dependencies](https://tinyverse.netlify.com/badge/IDSL.IPA)](https://cran.r-project.org/package=IDSL.IPA)
<!-- badges: end -->

**Intrinsic Peak Analysis (IPA)** by the [**Integrated Data Science Laboratory for Metabolomics and Exposomics (IDSL.ME)**](https://www.idsl.me) is a light-weight R package that extracts peaks for organic small molecules from untargeted liquid chromatography high resolution mass spectrometry (LC/HRMS) data in population scale projects. IDSL-IPA generates comprehensive and high-quality datasets from untargeted analysis of organic small molecules in a bio-specimen. It is a suite of new algorithms covering extracted ion chromatogram (EIC) candidate generation, peak detection, peak property evaluation, mass-correction, retention time correction across multiple batches and peak annotation. IDSL.IPA has been optimized and tested for analysis of large sample sizes (n = 499) samples. We have shown that IDSL.IPA in our [publication](https://github.com/idslme/IDSL.IPA#citation) is able to outperform similar peak picking tools such as MZmine 2, *xcms*, and MS-DIAL in terms of sensitivity, specificity and speed.

## <img src='IPA_educational_files/Figures/IDSL.IPA-TOC_Art.png' align="right" />

## Table of Contents

- [Features of IDSL.IPA](https://github.com/idslme/IDSL.IPA#features-of-idsl.ipa)
- [Installation](https://github.com/idslme/IDSL.IPA#installation)
- [Workflow](https://github.com/idslme/IDSL.IPA#workflow)
- [Simple Batch Example](https://github.com/idslme/IDSL.IPA#simple-batch-example)
- [Wiki](https://github.com/idslme/IDSL.IPA#wiki)
- [Citation](https://github.com/idslme/IDSL.IPA#citation)
- [News and Update](https://github.com/idslme/IDSL.IPA/update.md)

## Features of IDSL.IPA

1) Parameter selection through a well-described [parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.IPA/main/IPA_parameters.xlsx)
2) Process high-throughput and population size studies (n > 500)
3) Calculating 17 chromatographic peak properties
4) Mass spectra level [Ion Pairing](https://github.com/idslme/IDSL.IPA/wiki/Ion-Pairing) to remove random noises to accelerate the processing speed
5) Compatibility to screen for any ion mass difference in addition to natural carbon signatures (<sup>12</sup>C/<sup>13</sup>C isotopologues) mass difference
6) [Retention time correction](https://github.com/idslme/IDSL.IPA/wiki/Retention-Index) using endogenous reference markers
7) Generating batch untargeted extracted ion chromatograms (EICs)
8) Generating pairwise correlations list for aligned peak height and its gap-filled tables to detect potential recurring adducts, in-source products and fragment peaks
9) Aggregating untargeted EICs after (m/z-RT) annotation for each compound
10) Compatibility with parallel processing in Windows and Linux environments
11) Compatibility with downstream molecular formula annotation tools such as [IDSL.UFA](https://github.com/idslme/IDSL.UFA) and [IDSL.UFAx](https://github.com/idslme/IDSL.UFAx)
12) Compatibility with [IDSL.CSA](https://github.com/idslme/IDSL.CSA) workflow to cluster recurring ions to generate composite spectra

## Installation

	install.packages("IDSL.IPA")
	
**Note:** In case you want to process **netCDF/CDF** mass spectrometry data by IDSL.IPA, you should also install the [**RnetCDF**](https://CRAN.R-project.org/package=RNetCDF) package separately using the below command.

	install.packages("RNetCDF")

## Workflow
To process your mass spectrometry data (**mzXML**, **mzML**, **netCDF**), download the [IPA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.IPA/main/IPA_parameters.xlsx) and select the parameters accordingly and then use this spreadsheet as the input for the `IPA_workflow` function as shown below:

	library(IDSL.IPA)
	IPA_workflow("Address of the IPA parameter spreadsheet")

## Simple Batch Example

## [**Wiki**](https://github.com/idslme/IDSL.IPA/wiki)

1. [**Example of a population size study with 499 indivdual mass spectrometry file**](https://github.com/idslme/IDSL.IPA/wiki/IDSL.IPA-for-MTBLS1684-study)
2. `IPA_targeted` function for a large number of peaks (***m/z***-**RT** pairs) with an [**example for targeted IDSL.IPA**](https://github.com/idslme/IDSL.IPA/wiki/IPA_targeted)
3. [**Ion Pairing**](https://github.com/idslme/IDSL.IPA/wiki/Ion-Pairing)
4. [**Definition of Signal to Noise ratio (S/N)**](https://github.com/idslme/IDSL.IPA/wiki/Definition-Signal-to-Noise-Ratio)
5. [**nIsoPair/RCS**](https://github.com/idslme/IDSL.IPA/wiki/nIsoPair-RCS)
6. [**Ratio of peak width at half-height to peak width at the baseline (RPW)**](https://github.com/idslme/IDSL.IPA/wiki/RPW)
7. [**Chromatogram gap percentage**](https://github.com/idslme/IDSL.IPA/wiki/Chromatogram-gaps-percentage-(missing-scans))
8. [**Peak tailing fronting resolving method**](https://github.com/idslme/IDSL.IPA/wiki/Peak-tailing-fronting-resolving)
9. [**Peak smoothing**](https://github.com/idslme/IDSL.IPA/wiki/Peak-smoothing)
10. [**Extra scans**](https://github.com/idslme/IDSL.IPA/wiki/Extra-scans)
11. [**Retention time correction**](https://github.com/idslme/IDSL.IPA/wiki/Retention-Index).

## Citation

Fakouri Baygi, S., Kumar, Y. Barupal, D.K. [IDSL. IPA characterizes the organic chemical space in untargeted LC/HRMS datasets](https://pubs.acs.org/doi/10.1021/acs.jproteome.2c00120). *Journal of proteome research*, **2022**, *21(6)*, 1485-1494.