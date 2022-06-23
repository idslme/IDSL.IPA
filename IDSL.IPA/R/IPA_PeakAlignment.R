IPA_PeakAlignment <- function (PARAM) {
  number_processing_cores <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0006'), 2])
  ##
  output_path <- PARAM[which(PARAM[, 1] == 'PARAM0010'), 2]
  input_path_peaklist <- paste0(output_path, "/peaklists")
  file_names_peaklist1 <- dir(path = input_path_peaklist, pattern = ".Rdata")
  file_names_peaklist2 <- dir(path = input_path_peaklist, pattern = "peaklist_")
  file_names_peaklist <- file_names_peaklist1[file_names_peaklist1%in%file_names_peaklist2]
  L_PL <- length(file_names_peaklist)
  file_names_peaklist <- cbind(file_names_peaklist, 1:L_PL)
  ##
  input_path_hrms <- PARAM[which(PARAM[, 1] == 'PARAM0007'), 2]
  if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]) == "all") {
    file_name_hrms <- dir(path = input_path_hrms)
    file_name_hrms <- file_name_hrms[grep(paste0(".", tolower(PARAM[which(PARAM[, 1] == 'PARAM0009'), 2]), "$"), file_name_hrms, ignore.case = TRUE)]
  } else {
    samples_string <- PARAM[which(PARAM[, 1] == 'PARAM0008'), 2]
    file_name_hrms <- strsplit(samples_string, ";")[[1]] # files used as reference m/z-RT
  }
  ##
  file_names_peaklist_hrms1 <- gsub(".Rdata", "", file_names_peaklist[, 1])
  file_names_peaklist_hrms2 <- gsub("peaklist_", "", file_names_peaklist_hrms1)
  file_names_peaklist_hrms <- file_name_hrms%in%file_names_peaklist_hrms2
  if (length(file_names_peaklist_hrms) != L_PL) {
    stop("Error!!! peaklist files are not available for all selected HRMS files!")
  }
  ##
  RT_correction_needed <- PARAM[which(PARAM[, 1] == 'PARAM0029'), 2]
  mz_error <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0035'), 2])
  rt_tol <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0036'), 2])
  n_quantile <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0037'), 2])
  ##
  OutputPath_peak_alignment <- paste0(output_path, "/peak_alignment")
  if (!dir.exists(OutputPath_peak_alignment)) {
    dir.create(OutputPath_peak_alignment)
  }
  ##
  if (tolower(RT_correction_needed) == "yes") {
    reference_samples_string <- PARAM[which(PARAM[, 1] == 'PARAM0030'), 2]
    Ref_name <- strsplit(reference_samples_string, ";")[[1]] # files used as reference m/z-RT
    Ref_ID <- sapply(1:length(Ref_name), function(i) {
      ID_name <- paste0("peaklist_", Ref_name[i], ".Rdata")
      x <- which(file_names_peaklist == ID_name)
      as.numeric(file_names_peaklist[x, 2])
    })
    Ref_ID <- Ref_ID[order(Ref_ID)]
    ## To find common peaks in the reference samples
    file_names_peaklist_ref <- file_names_peaklist[Ref_ID, 1]
    print("Initiated detecting reference peaks for RT correction!")
    RT_peaklists_ref <- lapply(1:length(Ref_ID), function(i) {
      matrix(loadRdata(paste0(input_path_peaklist, "/", file_names_peaklist_ref[i]))[, 3], ncol = 1)
    })
    min_frequency_ref_peaks <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0031'), 2])
    reference_mz_rt_peaks <- reference_peaks_detector(input_path_peaklist, file_names_peaklist_ref, min_frequency_ref_peaks,
                                                      RT_peaklists_ref, mz_error, rt_tol, n_quantile, number_processing_cores)
    print(paste0("Detected " , dim(reference_mz_rt_peaks)[1], " reference peaks for RT correction!"))
    #
    png(paste0(OutputPath_peak_alignment, "/Ref_peaks_distribution.png"))
    Ref_peaks_distribution <- reference_mz_rt_peaks[, 2]
    hist_rt_reference_peaks <- hist(Ref_peaks_distribution, breaks = round(max(unlist(RT_peaklists_ref))*1), xlab = "Retention time (min)")
    dev.off()
    L_x_regions_rt0 <- length(which(hist_rt_reference_peaks[["counts"]] == 0))
    if (L_x_regions_rt0 > 0) {
      print("WARNING!!! Reference peaks were not detected for the entire range of the retention times! Please see the 'Ref_peaks_distribution.png' in the 'peak_alignment' folder!")
    }
    ##
    print("Initiated RT correction!")
    ##
    rt_correction_method <- PARAM[which(PARAM[, 1] == 'PARAM0032'), 2]
    reference_peak_tol <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0033'), 2])
    polynomial_degree <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0034'), 2])
    ##
    Sample_ID <- as.numeric(setdiff(file_names_peaklist[, 2], Ref_ID))
    file_name_peaklist_samples <- file_names_peaklist[Sample_ID, 1]
    ############################################################################
    if (number_processing_cores == 1) {
      Correted_RTs_peaklist <- lapply(1:length(Sample_ID), function(i) {
        peaklist <- loadRdata(paste0(input_path_peaklist, "/", file_name_peaklist_samples[i]))
        sample_rt_corrector(reference_mz_rt_peaks, peaklist, mz_error, rt_correction_method, reference_peak_tol, polynomial_degree)
      })
    } else {
      ## Processing OS
      osType <- Sys.info()[['sysname']]
      if(osType == "Windows") {
        cl <- makeCluster(number_processing_cores)
        registerDoParallel(cl)
        Correted_RTs_peaklist <- foreach(i = 1:length(Sample_ID), .verbose = FALSE) %dopar% {
          peaklist <- loadRdata(paste0(input_path_peaklist, "/", file_name_peaklist_samples[i]))
          sample_rt_corrector(reference_mz_rt_peaks, peaklist, mz_error, rt_correction_method, reference_peak_tol, polynomial_degree)
        }
        stopCluster(cl)
        ##
      } else if (osType == "Linux") {
        Correted_RTs_peaklist <- mclapply(1:length(Sample_ID), function(i) {
          peaklist <- loadRdata(paste0(input_path_peaklist, "/", file_name_peaklist_samples[i]))
          sample_rt_corrector(reference_mz_rt_peaks, peaklist, mz_error, rt_correction_method, reference_peak_tol, polynomial_degree)
        }, mc.cores = number_processing_cores)
        closeAllConnections()
      }
    }
    print("Completed RT correction!")
    ############################################################################
    corrected_RT_peaklists <- lapply(1:max(as.numeric(file_names_peaklist[, 2])), function(i) {
      if (i%in%Ref_ID == TRUE) {
        k <- which(Ref_ID == i)
        RT_peaklists_ref[[k]]
      } else {
        k <- which(Sample_ID == i)
        Correted_RTs_peaklist[[k]]
      }
    })
  } else {
    corrected_RT_peaklists <- lapply(1:max(as.numeric(file_names_peaklist[, 2])), function(i) {
      loadRdata(paste0(input_path_peaklist, "/", file_names_peaklist[i]))[, 3]
    })
  }
  ##
  save(corrected_RT_peaklists, file =paste0(OutputPath_peak_alignment, "/corrected_RT_peaklists.Rdata"))
  ##
  if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0002'), 2]) == "yes") {
    print("Initiated peak alignment among the entire samples using peak index numbers in the peaklists!")
    peak_Xcol <- peak_alignment(input_path_peaklist, file_names_peaklist[, 1], corrected_RT_peaklists, mz_error, rt_tol, n_quantile, number_processing_cores)
    #
    file_names_hrms <- gsub(".Rdata", "", (gsub("peaklist_", "", file_names_peaklist[, 1])))
    colnames(peak_Xcol) <- c("m/z", "RT", file_names_hrms)
    save(peak_Xcol, file =paste0(OutputPath_peak_alignment, "/peak_Xcol.Rdata"))
    print("Peak index numbers from individual peaklists were stored in 'peak_Xcol.Rdata'")
    ##
    print("Producing tables for the peak height, peak area, and R13C")
    H_A_R13C <- peak_Xcol2(input_path_peaklist, file_names_peaklist[, 1], peak_Xcol)
    print("Inititated saving peak tables")
    peak_height <- H_A_R13C[[1]]
    peak_area <- H_A_R13C[[2]]
    peak_R13C <- H_A_R13C[[3]]
    opendir(OutputPath_peak_alignment)
    save(peak_height, file =paste0(OutputPath_peak_alignment, "/peak_height.Rdata"))
    write.csv(peak_height, file =paste0(OutputPath_peak_alignment, "/peak_height.csv"))
    save(peak_area, file =paste0(OutputPath_peak_alignment, "/peak_area.Rdata"))
    write.csv(peak_area, file =paste0(OutputPath_peak_alignment, "/peak_area.csv"))
    save(peak_R13C, file =paste0(OutputPath_peak_alignment, "/peak_R13C.Rdata"))
    write.csv(peak_R13C, file =paste0(OutputPath_peak_alignment, "/peak_R13C.csv"))
    print("Completed peak alignment!!!")
  }
}
