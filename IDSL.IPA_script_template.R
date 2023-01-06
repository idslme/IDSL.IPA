install.packages("IDSL.IPA")
library(IDSL.IPA)
spreadsheet <- data.frame(readxl::read_xlsx(system.file("inst/extdata","IPA_parameters.xlsx",package = "IDSL.IPA")))
#######################################################
######################  Global parameters (required)
#######################################################

# PARAM0001  : Generate peaklist for individual file
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0001")] <- "YES"

# PARAM0002  : Generate the aligned peak table
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0002")] <- "YES"

# PARAM0003  : Perform gap-filling for the aligned peak table
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0003")] <- "YES"

# PARAM0004 : Annotate peaks from individual peaklists using a reference database
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0004")] <- "YES"

# PARAMA005 : Employ IPA for targeted analysis (Applicable only for the 'IPA_TargetedAnalysis' module)
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0005")] <- "NO"

# PARAMA006 : Number of processing threads
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0006")] <- 10


#######################################################
######################  Data import and export (required)
#######################################################

# PARAMA007 : Input data location (MS1 level HRMS data)
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0007")] <- "E:/MTBLS136/rp_early_pos/"

# PARAM008 : List of HRMS files to process
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0008")] <- "All"

# PARAM009 : Input HRMS data format
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0009")] <- "mzML"

# PARAM010 : Output location (MS1 processed data)
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0010")] <- "E:/MTBLS136/rp_early_pos/"

# PARAM_EIC : Save extracted ion chromatogram (EIC) figures for individual files
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM_EIC")] <- "NO"

#######################################################
######################  Pairing isotopologues in the mass spectra level
#######################################################

# PARAM011 : Intensity cutoff for the monoisotopic mass in each scan also known as instrument noise level
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0011")] <- 50000

# PARAM012 : Mass difference to pair isotopologues ((Default = DeltaC = 13C - 12C = 1.003354835336))
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0012")] <- 1.003354835336

#######################################################
######################  Chromatographic peak detection
#######################################################

# PARAM013 : Mass accuracy to cluster m/zs and produce extracted ion chromatograms (XICs) (Da)
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0013")] <- 0.01

# PARAM014 : The maximum absolute retention time deviation to remove redundant peaks (min)
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0014")] <- 0.01

# PARAM015 : Smoothing window (number of scans for smoothing)
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0015")] <- 11

# PARAM017 : Peak tailing/fronting resolving power
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0017")] <- 0.05

# PARAM018 : Power of mass spectral size reduction by removing [non-12C/13C] m/zs to accelerate calculation of peaklists
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0018")] <- 0


#######################################################
######################  Peak analysis and data reduction
#######################################################

# PARAM019 : Perform recursive mass correction
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0019")] <- "No"

# PARAM020 : Number of extra scans on the both sides of detected peak boundaries in the recursive analysis
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0020")] <- 50

# PARAM021 : Intensity height threshold for the chromatographic peak
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0021")] <- 100000

# PARAM022 : Maximum percentage of missing scans (gaps) on the top segment of raw chromatographic peaks (%)
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0022")] <- 30

# PARAM023 : Minimum number of scans having 12C/13C isotopologue pairs within the chromatographic peak (nIsoPair)
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0023")] <- 3

# PARAM024 : Minimum ratio of the nIsoPair per number of available scans above 80% of peak height (%) (RCS)
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0024")] <- 30

# PARAM025 : Maximum ratio of 13C for the integrated spectra across the chromatographic peak (%) (Average R13C)
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0025")] <- 80

# PARAM026 : Maximum ratio of peak width at half-height to peak width at the baseline (RPW)
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0026")] <- 0.8

# PARAM027 : Minimum S/N for chromatographic peak
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0027")] <- 2

# PARAM028 : Number of points to smooth individual chromatographic peaks using cubic spline method to calculate ancillary chromatography parameters
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0028")] <- 100

#######################################################
######################  Retention time correction and peak alignment
#######################################################

# PARAM029 : Perform RT correction
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0029")] <- "Yes"

# PARAM030 : Reference samples to detect reoccurring reference peaks for RT correction
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0030")] <- "QXA03POSEAR20160829_NCIA0316ML_HUMAN_SERUM1_03.mzML;QXA03POSEAR20160829_NCIA0316ML_HUMAN_SERUM1_04.mzML;QXA03POSEAR20160829_NCIA0316ML_HUMAN_SERUM1_05.mzML;QXA03POSEAR20160829_NCIA0316ML_HUMAN_SERUM1_06.mzML;QXA03POSEAR20160829_NCIA0316ML_HUMAN_SERUM1_07.mzML;QXA03POSEAR20160829_NCIA0316ML_HUMAN_SERUM1_08.mzML;QXA03POSEAR20160829_NCIA0316ML_HUMAN_SERUM1_09.mzML;QXA03POSEAR20160829_NCIA0316ML_HUMAN_SERUM1_10.mzML;QXA03POSEAR20160829_NCIA0316ML_HUMAN_SERUM1_12.mzML;QXA03POSEAR20160829_NCIA0316ML_HUMAN_SERUM1_13.mzML"

# PARAM031 : Minimum frequency of the reoccurring reference peaks in the reference samples (%)
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0031")] <- 50

# PARAM032 : Retention time correction method
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0032")] <- "RetentionIndex"

# PARAM033 : Reference Peak tolerance for "RetentionIndex" to minimize local RT errors
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0033")] <- 5

# PARAM034 : Polynomial degree for "Polynomial" regression
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0034")] <- ""

# PARAM035 : Mass accuracy to cluster integrated m/z values among multiple samples (Da)
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0035")] <- 0.01

# PARAM036 : Retention time tolerance to detect common peaks (min)
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0036")] <- 0.05

# PARAM037 : Number of slices for the (m/z-RT) alignment
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0037")] <- 20


#######################################################
######################  Gap-filling on the peak table
#######################################################

# PARAM038 : Mass accuracy to cluster integrated m/z values among multiple samples (Da)
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0038")] <- 0.01

# PARAM039 : Retention time tolerance to detect common peaks (min)
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0039")] <- 0.01

# PARAM040 : Number of extra scans on the both sides of peak boundaries to fill
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0040")] <- 20

#######################################################
######################  Peak annotation on individual peaklists using a reference database of m/z and RT
#######################################################

# PARAM041 : Reference file location
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0041")] <- "E:/MTBLS136/rp_early_pos/"

# PARAM042 : Name of the reference file
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0042")] <- "rp_early_pos_targets.xlsx"

# PARAM043 : Mass error (Da)
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0043")] <- 0.01

# PARAM044 : Maximum retention time deviation (min)
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0044")] <- 0.05

# PARAM045 : Use corrected retention times to minimize RT fluctuations in peak annotation
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0045")] <- "YES"

# PARAM046 : Compound-centric annotation
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0046")] <- "YES"

# PARAM047 : Sample-centric annotation
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0047")] <- "YES"

# PARAM048 :Perform gap-filling for sample-centric annotation
spreadsheet$User.provided.input[which(spreadsheet$Parameter.ID=="PARAM0048")] <- "YES"


#######################################################
######################  RUN IDSL.IPA pipeline
#######################################################

IDSL.IPA::IPA_workflow(spreadsheet)

