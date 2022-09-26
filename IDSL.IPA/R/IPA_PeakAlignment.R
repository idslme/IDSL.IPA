IPA_PeakAlignment <- function(PARAM) {
  ##
  IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
  IPA_logRecorder("Inititated generating data for peak alignment!")
  ##
  number_processing_threads <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0006'), 2])
  output_path <- PARAM[which(PARAM[, 1] == 'PARAM0010'), 2]
  input_path_peaklist <- paste0(output_path, "/peaklists")
  file_names_peaklist1 <- dir(path = input_path_peaklist, pattern = ".Rdata")
  file_names_peaklist2 <- dir(path = input_path_peaklist, pattern = "peaklist_")
  file_names_peaklist <- file_names_peaklist1[file_names_peaklist1 %in% file_names_peaklist2]
  L_PL <- length(file_names_peaklist)
  ##
  input_path_hrms <- PARAM[which(PARAM[, 1] == 'PARAM0007'), 2]
  if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]) == "all") {
    file_names_hrms <- dir(path = input_path_hrms)
    file_names_hrms <- file_names_hrms[grep(paste0(".", tolower(PARAM[which(PARAM[, 1] == 'PARAM0009'), 2]), "$"), file_names_hrms, ignore.case = TRUE)]
  } else {
    samples_string <- PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]
    file_names_hrms <- strsplit(samples_string, ";")[[1]] # files used as reference m/z-RT
  }
  ##
  peaklist2HRMS <- gsub("^peaklist_", "", gsub(".Rdata$", "", file_names_peaklist))
  matchPeaklist2HRMS <- file_names_hrms[(file_names_hrms %in% peaklist2HRMS)]
  if (length(matchPeaklist2HRMS) != L_PL) {
    stop(IPA_logRecorder("Error!!! peaklist files are not available for all selected HRMS files!"))
  }
  ##
  RTcorrectionCheck <- tolower(PARAM[which(PARAM[, 1] == 'PARAM0029'), 2])
  mz_error <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0035'), 2])
  rt_tol <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0036'), 2])
  n_quantile <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0037'), 2])
  ##
  OutputPath_peak_alignment <- paste0(output_path, "/peak_alignment")
  if (!dir.exists(OutputPath_peak_alignment)) {
    dir.create(OutputPath_peak_alignment, recursive = TRUE)
  }
  ##
  IPA_logRecorder("Peak alignment data are stored in the `peak_alignment` folder!")
  ##
  if (RTcorrectionCheck == "yes") {
    reference_samples_string <- PARAM[which(PARAM[, 1] == 'PARAM0030'), 2]
    Ref_name <- strsplit(reference_samples_string, ";")[[1]] # files used as reference m/z-RT
    file_names_peaklist_ref <- paste0("peaklist_", Ref_name, ".Rdata")
    ## To find common peaks in the reference samples
    IPA_logRecorder("Initiated detecting reference peaks for RT correction!")
    min_frequency_ref_peaks <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0031'), 2])
    ##
    listReferencePeaks <- reference_peaks_detector(input_path_peaklist, file_names_peaklist_ref, min_frequency_ref_peaks,
                                                   mz_error, rt_tol, n_quantile, number_processing_threads)
    reference_mz_rt_peaks <- listReferencePeaks[["referenceMZRTpeaks"]]
    listRefRT <- listReferencePeaks[["listRefRT"]]
    ##
    IPA_logRecorder(paste0("Detected " , dim(reference_mz_rt_peaks)[1], " reference peaks for RT correction!"))
    #
    png(paste0(OutputPath_peak_alignment, "/Ref_peaks_distribution.png"))
    Ref_peaks_distribution <- reference_mz_rt_peaks[, 2]
    hist_rt_reference_peaks <- hist(Ref_peaks_distribution, breaks = round(max(unlist(listRefRT))*1), xlab = "Retention time (min)")
    dev.off()
    L_x_regions_rt0 <- length(which(hist_rt_reference_peaks[["counts"]] == 0))
    if (L_x_regions_rt0 > 0) {
      IPA_logRecorder("WARNING!!! Reference peaks were not detected for the entire range of the retention times! Please see the 'Ref_peaks_distribution.png' in the `peak_alignment` folder!")
    }
    ##
    IPA_logRecorder("Initiated RT correction!")
    ##
    rt_correction_method <- PARAM[which(PARAM[, 1] == 'PARAM0032'), 2]
    reference_peak_tol <- tryCatch(as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0033'), 2]), error = function(e) {5}, warning = function(w) {5})
    polynomial_degree <- tryCatch(as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0034'), 2]), error = function(e) {3}, warning = function(w) {3})
    ##
    file_name_peaklist_samples <- setdiff(file_names_peaklist, file_names_peaklist_ref)
    ############################################################################
    if (number_processing_threads == 1) {
      correted_RTs_samples <- lapply(file_name_peaklist_samples, function(i) {
        peaklist <- loadRdata(paste0(input_path_peaklist, "/", i))
        sample_rt_corrector(reference_mz_rt_peaks, peaklist, mz_error, rt_correction_method, reference_peak_tol, polynomial_degree)
      })
    } else {
      ## Processing OS
      osType <- Sys.info()[['sysname']]
      if (osType == "Windows") {
        cl <- makeCluster(number_processing_threads)
        registerDoParallel(cl)
        correted_RTs_samples <- foreach(i = file_name_peaklist_samples, .verbose = FALSE) %dopar% {
          peaklist <- loadRdata(paste0(input_path_peaklist, "/", i))
          sample_rt_corrector(reference_mz_rt_peaks, peaklist, mz_error, rt_correction_method, reference_peak_tol, polynomial_degree)
        }
        stopCluster(cl)
        ##
      } else if (osType == "Linux") {
        correted_RTs_samples <- mclapply(file_name_peaklist_samples, function(i) {
          peaklist <- loadRdata(paste0(input_path_peaklist, "/", i))
          sample_rt_corrector(reference_mz_rt_peaks, peaklist, mz_error, rt_correction_method, reference_peak_tol, polynomial_degree)
        }, mc.cores = number_processing_threads)
        closeAllConnections()
      }
    }
    names(correted_RTs_samples) <- file_name_peaklist_samples
    ##
    ############################################################################
    ##
    corrected_RT_peaklists <- lapply(file_names_peaklist, function(i) {
      if (i %in% file_names_peaklist_ref) {
        listRefRT[[i]]
      } else {
        correted_RTs_samples[[i]]
      }
    })
    ##
    IPA_logRecorder("Corrected retention times for individual peaklists are stored as `corrected_RT_peaklists.Rdata` in the `peak_alignment` folder!")
    IPA_logRecorder("Completed RT correction!")
    ##
  } else {
    corrected_RT_peaklists <- lapply(file_names_peaklist, function(i) {
      loadRdata(paste0(input_path_peaklist, "/", i))[, 3]
    })
  }
  ##
  names(corrected_RT_peaklists) <- file_names_peaklist
  ##
  save(corrected_RT_peaklists, file = paste0(OutputPath_peak_alignment, "/corrected_RT_peaklists.Rdata"))
  ##
  if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0002'), 2]) == "yes") {
    IPA_logRecorder("Initiated peak alignment on the entire peaklists using peak IDs in each peaklist!")
    peak_Xcol <- peak_alignment(input_path_peaklist, file_names_peaklist, corrected_RT_peaklists, mz_error, rt_tol, n_quantile, number_processing_threads)
    #
    colnames(peak_Xcol) <- c("m/z", "RT", file_names_hrms)
    save(peak_Xcol, file = paste0(OutputPath_peak_alignment, "/peak_Xcol.Rdata"))
    IPA_logRecorder("Stored aligned indexed table as 'peak_Xcol.Rdata' in the `peak_alignment` folder!")
    ##
    IPA_logRecorder("Initiated generating aligned peak tables for the peak height, peak area, and R13C values!")
    ##
    listHeightAreaR13C <- peak_Xcol2(input_path_peaklist, file_names_peaklist, peak_Xcol)
    peak_height <- listHeightAreaR13C[["peak_height"]]
    peak_area <- listHeightAreaR13C[["peak_area"]]
    peak_R13C <- listHeightAreaR13C[["peak_R13C"]]
    ##
    IPA_logRecorder("Initiated saving aligned peak tables")
    opendir(OutputPath_peak_alignment)
    ##
    save(peak_height, file = paste0(OutputPath_peak_alignment, "/peak_height.Rdata"))
    write.csv(peak_height, file = paste0(OutputPath_peak_alignment, "/peak_height.csv"))
    ##
    save(peak_area, file = paste0(OutputPath_peak_alignment, "/peak_area.Rdata"))
    write.csv(peak_area, file = paste0(OutputPath_peak_alignment, "/peak_area.csv"))
    ##
    save(peak_R13C, file = paste0(OutputPath_peak_alignment, "/peak_R13C.Rdata"))
    write.csv(peak_R13C, file = paste0(OutputPath_peak_alignment, "/peak_R13C.csv"))
    ##
    IPA_logRecorder("Aligned peak height, peak area, and R13C tables were stored in `.Rdata` and `.csv` formats in the `peak_alignment` folder!")
    ##
    ############################################################################
    ##
    IPA_logRecorder("Initiated detecting correlating peaks on the peak height table!")
    ##
    correlationListHeight <- peak_property_table_correlation(peakPropertyTable = peak_height, RTtolerance = rt_tol, minFreqDetection = 1, method = "pearson", minThresholdCorrelation = 0.50, number_processing_threads)
    ##
    save(correlationListHeight, file = paste0(OutputPath_peak_alignment, "/correlationListHeight.Rdata"))
    ##
    IPA_logRecorder("Stored the correlating peaks on the peak height table as `correlationListHeight.Rdata` in the `peak_alignement` folder!")
    IPA_logRecorder("Completed generating data for peak alignment!")
    IPA_logRecorder(paste0(rep("", 100), collapse = "-"))
    ##
    ############################################################################
    ##
  }
}
