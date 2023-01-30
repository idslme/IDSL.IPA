# IDSL.IPA<img src='IPA_educational_files/Figures/IDSL.IPA-logo.PNG' width="250px" align="right" />

<!-- badges: start -->
[![Maintainer](https://img.shields.io/badge/maintainer-Sadjad_Fakouri_Baygi-blue)](https://github.com/sajfb)
[![CRAN status](https://www.r-pkg.org/badges/version/IDSL.IPA)](https://cran.r-project.org/package=IDSL.IPA)
![](http://cranlogs.r-pkg.org/badges/IDSL.IPA?color=orange)
![](http://cranlogs.r-pkg.org/badges/grand-total/IDSL.IPA?color=brightgreen)
[![Dependencies](https://tinyverse.netlify.com/badge/IDSL.IPA)](https://cran.r-project.org/package=IDSL.IPA)
<!-- badges: end -->

**Intrinsic Peak Analysis (IPA)** by the [**Integrated Data Science Laboratory for Metabolomics and Exposomics (IDSL.ME)**](https://www.idsl.me) is a light-weight R package that extracts peaks for organic small molecules from untargeted liquid chromatography high resolution mass spectrometry (LC/HRMS) data in population scale projects. IDSL.IPA is a suite of new algorithms covering extracted ion chromatogram (EIC) candidate generation, peak detection, peak property evaluation, recursive mass correction, retention time correction across multiple batches and peak annotation. IDSL.IPA generates comprehensive and high-quality datasets from untargeted analysis of organic small molecules for population-size studies. We have shown in our [publication](https://github.com/idslme/IDSL.IPA#citation) that IDSL.IPA is able to outperform similar peak picking tools such as MZmine 2, *xcms*, and MS-DIAL in terms of sensitivity, specificity and speed.

## <img src='IPA_educational_files/Figures/IDSL.IPA-TOC_Art.png' align="right" />

## Table of Contents

- [Features of IDSL.IPA](https://github.com/idslme/IDSL.IPA#features-of-idslipa)
- [Installation](https://github.com/idslme/IDSL.IPA#installation)
- [Workflow](https://github.com/idslme/IDSL.IPA#workflow)
- [Quick Batch Example](https://github.com/idslme/IDSL.IPA#quick-batch-example)
- [Wiki](https://github.com/idslme/IDSL.IPA#wiki)
- [News and Updates](https://github.com/idslme/IDSL.IPA#news-and-updates)
- [Citation](https://github.com/idslme/IDSL.IPA#citation)

## Features of IDSL.IPA

1) Parameter selection through a user-friendly and well-described [parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.IPA/main/IPA_parameters.xlsx)
2) Analyzing population size untargeted studies (n > 500)
3) Calculating 19 chromatographic peak properties such as peak area, nIsoPair, RCS, cumulated intensity, R<sup>13</sup>C, peak width, RPW, number of separation trays, asymmetry factor, USP tailing factor, skewness using derivative method, symmetry using pseudo-moments, skewness using pseudo-moments, gaussianity, S/N using baseline, S/N using the *xcms* method, S/N using the RMS method, and sharpness.
4) Mass spectra scan level [Ion Pairing](https://github.com/idslme/IDSL.IPA/wiki/Ion-Pairing) to only select relevant ions
5) Flexibility to screen for any ion mass difference in addition to natural carbon signatures (<sup>12</sup>C/<sup>13</sup>C isotopologues) mass difference
6) [Retention time correction](https://github.com/idslme/IDSL.IPA/wiki/Retention-Index) using endogenous reference markers for multi-batch large scale studies
7) Generating batch untargeted extracted ion chromatograms (EICs)
8) Generating pairwise correlations list for aligned peak height and its gap-filled tables to detect potential recurring adducts, in-source products and fragment peaks
9) Aggregating untargeted EICs after (*m/z*-RT) annotation for each compound
10) Parallel processing in Windows and Linux environments
11) Integration with molecular formula annotation tools [IDSL.UFA](https://github.com/idslme/IDSL.UFA) and [IDSL.UFAx](https://github.com/idslme/IDSL.UFAx)
12) Integration with [IDSL.CSA](https://github.com/idslme/IDSL.CSA) workflow to cluster recurring ions to generate composite spectra

## Installation

	install.packages("IDSL.IPA")
	
**Note:** In case you want to process **netCDF/CDF** mass spectrometry data by IDSL.IPA, you should also install the [**RnetCDF**](https://CRAN.R-project.org/package=RNetCDF) package separately using the below command.

	install.packages("RNetCDF")

## Workflow

To process your mass spectrometry data (**mzXML**, **mzML**, **netCDF**), download the [IPA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.IPA/main/IPA_parameters.xlsx) and select the parameters accordingly and then use this spreadsheet as the input for the `IPA_workflow` function as shown below:

	library(IDSL.IPA)
	IPA_workflow("Address of the IPA parameter spreadsheet")

## Quick Batch Example

Follow these steps for a quick case study (n = 33) [ST002263](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST002263&DataMode=AllData&ResultType=1) which has Thermo Q Exactive HF hybrid Orbitrap data collected in the HILIC-ESI-POS/NEG modes.

1. Download the file ["ST002263_Rawdata.zip" (1.6G)](https://www.metabolomicsworkbench.org/data/DRCCStudySummary.php?Mode=SetupRawDataDownload&StudyID=ST002263)

2. Separate the positive and negative modes *.mzXML* data into different folders and process each mode separately. Positive and negative modes data must be processed separately.

3. Open the [IPA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.IPA/main/IPA_parameters.xlsx) and use default values for all parameters, except:
	
	3.1. Input data location: In **PARAM0007**, specify the location of the MS1 level HRMS data.
	
	3.2. Output location: In **PARAM0010**, specify the location where you want the processed data to be stored.
		
	3.3. Number of processing threads: In **PARAM0006**, increase the number of processing threads based on your computer's computational power.

4. Open the R/Rstudio console or terminal and run the following command:

```
library(IDSL.IPA)
IPA_workflow("Address of the IPA parameter spreadsheet")
```

5. The results will be available in the output location specified in **PARAM0010** and will include:

	5.1. Individual peaklists for each HRMS file in *.Rdata* and *.csv* formats in the "peaklists" directory.
	
	5.2. Peak alignment tables in the "peak_alignment" directory.
	
	5.3. (Optional) untargeted EICs in the "IPA_EIC" directory for each HRMS file, if selected **YES** in **PARAM0009**.

## [**Wiki**](https://github.com/idslme/IDSL.IPA/wiki)

1. [**Example of a population size study with 499 individual mass spectrometry file**](https://github.com/idslme/IDSL.IPA/wiki/IDSL.IPA-for-MTBLS1684-study)
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

## News and Updates

We post major changes in the IDSL.IPA workflow [here](https://github.com/idslme/IDSL.IPA/blob/main/UPDATE.md).

## Citation

[1] Fakouri Baygi, S., Kumar, Y. Barupal, D.K. [IDSL. IPA characterizes the organic chemical space in untargeted LC/HRMS datasets](https://pubs.acs.org/doi/10.1021/acs.jproteome.2c00120). *Journal of proteome research*, **2022**, *21(6)*, 1485-1494.