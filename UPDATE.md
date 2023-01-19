# News and Updates
We slightly revised the IDSL.IPA algorithm relative to our publication in the *Journal of Proteome Research* that are listed for each section:

## Generating untargeted peaklists for individual HRMS files
- **PARAM_PAR** were added to be consistent with **sample mode** and **peak mode** parallelization modes to generate individual peaklists.
- **PARAM0009** was updated to save untargeted extracted ion chromatogram (EIC) figures for individual HRMS files.
- **PARAM0012** has changed its functionality and allows users to decide the mass difference between the isotopologues. This feature allows to screen for C, Cl, Br, S, Se, ... isotopologues or any constant mass difference. Furthermore, negative values are also allowed for this parameter to search  for lower mass isotopologues.
- Gaussianity of the chromatographic peaks are calculated using the Pearson correlation coefficient instead of Kolmogorov–Smirnov test.

## Generating aligned peak tables and correlation lists for the aligned peak height table
- **PARAM0037** also changed its original functionality to flag suspicious peaks on the aligned tables.
- **PARAM_ALG1-ALG4** parameters were added to calculate correlation of aligned peak heights to detect related peaks using Pearson correlation coefficients. These correlation results are stored as `alignedPeakHeightTableCorrelationList.Rdata` in an R list format. Each (*m/z*-RT) on the aligned table is presented with its correlating peaks.

## Annotating peaks from individual peaklists using a reference database of *m/z*-RT-name
- Compound-centric annotation can aggregate untargeted EICs in separate folders for each compound when **PARAM0009** was selected **YES** in the peak-picking step.