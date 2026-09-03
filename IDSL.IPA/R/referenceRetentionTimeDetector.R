referenceRetentionTimeDetector <- function(inputPathPeaklist, inputPathMZML, refPeaklistFileNames, minFrequencyRefPeaks,
                                           massAccuracy, RTtolerance, RTcorrectionIonSource, number_processing_threads = 1) {
  ##
  ####################### Polynomial Regression ################################
  ##
  if (gsub(" ", "", tolower(RTcorrectionIonSource)) == "ms2") {
    ##
    peaklistSource <- paste0(inputPathPeaklist, "/tempMS2LevelPLFolder")
    tryCatch(dir.create(peaklistSource, recursive = TRUE, showWarnings = FALSE), error = function(e) {stop(paste0("Temporary folder can't be created in the `", peaklistSource, "` folder!"))})
    ##
    refMzMLFileNames <- gsub("^peaklist_|.Rdata$", "", refPeaklistFileNames)
    ##
    MS2RefDetector <- function(mzML) {
      ##
      scanTable <- IDSL.MXP::peak2list(paste0(inputPathMZML, "/", mzML), onlyScanTable = TRUE)
      ##
      msLevel2 <- which(scanTable[["msLevel"]] == 2)
      ##
      dummyPL  <-  matrix(NA, nrow = length(msLevel2), ncol = 8)
      dummyPL[, 3] <- scanTable[msLevel2, "retentionTime"]
      dummyPL[, 4] <- scanTable[msLevel2, "precursorIntensity"]
      dummyPL[, 8] <- scanTable[msLevel2, "precursorMZ"]
      ##
      peaklist_name = paste0("peaklist_", mzML, ".Rdata")
      #
      save(dummyPL, file = paste0(peaklistSource, "/", peaklist_name))
      return(peaklist_name)
    }
    ##
    if (number_processing_threads == 1) {
      ##
      progressBARboundaries <- txtProgressBar(min = 0, max = length(refMzMLFileNames), initial = 0, style = 3)
      ##
      refPeaklistFileNames = rep("", length(refMzMLFileNames))
      for (i in 1:length(refMzMLFileNames)) {
        refPeaklistFileNames[i] <- tryCatch(MS2RefDetector(refMzMLFileNames[i]), error = function(e) {IPA_logRecorder(paste0("Problem with extracting MS2 reference ion markers for`", refMzMLFileNames[i],"`!"))})
        ##
        setTxtProgressBar(progressBARboundaries, i)
      }
      ##
      close(progressBARboundaries)
      ##
      ##########################################################################
      ##
    } else {
      ## Processing OS
      osType <- Sys.info()[['sysname']]
      ##
      if (osType == "Windows") {
        ##
        clust <- makeCluster(number_processing_threads)
        clusterExport(clust, setdiff(ls(), c("clust", "refMzMLFileNames")), envir = environment())
        ##
        refPeaklistFileNames <- do.call(c, parLapplyLB(clust, refMzMLFileNames, function(mzML) {
          tryCatch(MS2RefDetector(mzML), error = function(e) {IPA_logRecorder(paste0("Problem with extracting MS2 reference ion markers for`", mzML,"`!"))})
        }))
        ##
        stopCluster(clust)
        ##
        ########################################################################
        ##
      } else {
        ##
        refPeaklistFileNames <- do.call(c, mclapply(refMzMLFileNames, function(mzML) {
          tryCatch(MS2RefDetector(mzML), error = function(e) {IPA_logRecorder(paste0("Problem with extracting MS2 reference ion markers for`", mzML,"`!"))})
        }, mc.cores = number_processing_threads, mc.preschedule = FALSE))
        ##
        closeAllConnections()
        ##
      }
    }
    ##
    ############################################################################
    ##
  } else {
    peaklistSource <- inputPathPeaklist
  }
  ##
  ##############################################################################
  ##
  L_sS <- length(refPeaklistFileNames)
  ##
  listRefRT <- lapply(refPeaklistFileNames, function(i) {
    RTvec <- loadRdata(paste0(peaklistSource, "/", i))[, 3]
    if (is.na(RTvec[1])) {
      stop(IPA_logRecorder(paste0("IMPORTANT: `", i, "` CANNOT be a reference file in PARAM0030!")))
    }
    as.numeric(RTvec)
  })
  ##
  names(listRefRT) <- refPeaklistFileNames
  ##
  refPeakXcol <- peakAlignmentCore(peaklistSource, refPeaklistFileNames, listRefRT,
                                   massAccuracy, RTtolerance, number_processing_threads)
  ##
  ##############################################################################
  ##############################################################################
  ##
  if (gsub(" ", "", tolower(RTcorrectionIonSource)) == "ms2") {
    ##
    if (!is.null(peaklistSource)) {
      tryCatch(unlink(peaklistSource, recursive = TRUE), error = function(e) {IPA_logRecorder(paste0("Can't delete `", peaklistSource, "`!"))})
    }
    ##
    ############################################################################
    ##
    listRefRT <- lapply(refPeaklistFileNames, function(i) {
      RTvec <- loadRdata(paste0(inputPathPeaklist, "/", i))[, 3]
      if (is.na(RTvec[1])) {
        stop(IPA_logRecorder(paste0("IMPORTANT: `", i, "` CANNOT be a reference file in PARAM0030!")))
      }
      as.numeric(RTvec)
    })
    ##
    names(listRefRT) <- refPeaklistFileNames
  }
  ##
  ##############################################################################
  ##############################################################################
  ##
  refPeakXcol[, 3] <- round(refPeakXcol[, 3]/L_sS*100, digits = 0)
  x_ref <- which(refPeakXcol[, 3] >= minFrequencyRefPeaks)
  mz_rt_Xmed_ref <- matrix(refPeakXcol[x_ref, 1:3], ncol = 3)
  ## To remove isomeric peaks
  round_mz <- round(mz_rt_Xmed_ref[, 1], digits = 1)
  x_unique <- which(table(round_mz) == 1)
  unique_mz <- as.numeric(names(x_unique))
  selectedMZ <- which(round_mz %in% unique_mz)
  referenceMZRTpeaks <- matrix(mz_rt_Xmed_ref[selectedMZ, ], ncol = 3)
  colnames(referenceMZRTpeaks) <- c("m/z", "RT", "freqRefPeaks(%)")
  rownames(referenceMZRTpeaks) <- NULL
  ##
  listReferencePeaks <- list(referenceMZRTpeaks, listRefRT)
  names(listReferencePeaks) <- c("referenceMZRTpeaks", "listRefRT")
  ##
  return(listReferencePeaks)
}
