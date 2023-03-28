gapFillingCore <- function(input_path_hrms, peakXcol, massAccuracy, RTtolerance, scanTolerance,
                           retentionTimeCorrectionCheck = FALSE, listCorrectedRTpeaklists = NULL,
                           inputPathPeaklist = NULL, ionMassDifference = 1.003354835336, number_processing_threads = 1) {
  ##
  call_gapFillingCore <- function(i) {
    ##
    chromatography_undetected <- NULL
    x0 <- which(peakXcol[, (i + 3)] == 0)
    Lx0 <- length(x0)
    if (Lx0 > 0) {
      ##
      iHRMSfileName <- colnamesPeakXcol[i]
      undetecedMZ <- peakXcol[x0, 1]
      undetecedRT <- peakXcol[x0, 2]
      ##
      retentionTimeNAcheck <- TRUE
      if (retentionTimeCorrectionCheck) {
        iPeaklistFileName <- paste0("peaklist_", iHRMSfileName, ".Rdata")
        correctedRTi <- listCorrectedRTpeaklists[[iPeaklistFileName]]
        if (!is.na(correctedRTi[1])) {
          ## To back calculate the RT ##
          uncorrectedRTi <- loadRdata(paste0(inputPathPeaklist, "/", iPeaklistFileName))[, 3]
          idf <- data.frame(uncoRT = uncorrectedRTi, coRT = correctedRTi)
          rtmodel <- lm(coRT ~ poly(uncoRT, 3), idf)
          new.df <- data.frame(uncoRT = undetecedRT) ## predict uncorrected RTs
          undetecedRT <- predict(rtmodel, new.df)
        } else {
          retentionTimeNAcheck <- FALSE
        }
      }
      ##
      if (retentionTimeNAcheck) {
        ##
        outputer <- IPA_MSdeconvoluter(input_path_hrms, iHRMSfileName)
        spectraList <- outputer[["spectraList"]]
        retentionTime <- outputer[["retentionTime"]]
        LretentionTime <- length(retentionTime)
        outputer <- NULL
        aggregatedSpectraList <- IPA_spectraListAggregator(spectraList)
        ##
        chromatography_undetected <- do.call(rbind, lapply(1:Lx0, function(j) {
          ##
          peakArea <- 0
          R13C <- 0
          mzCandidate <- undetecedMZ[j]
          rtCandidate <- undetecedRT[j]
          scanNumberApex <- which.min(abs(rtCandidate - retentionTime)) # scan number at apex
          scanNumberStart <- max(c((scanNumberApex - scanTolerance), 1))
          scanNumberEnd <- min(c(LretentionTime, (scanNumberApex + scanTolerance)))
          chromatogramMatrix <- XIC(aggregatedSpectraList, scanNumberStart, scanNumberEnd, mzCandidate, massAccuracy)
          if (!is.null(chromatogramMatrix)) {
            length_chrom <- dim(chromatogramMatrix)[1]
            x_apex <- which(chromatogramMatrix[, 1] == scanNumberApex)
            rt_loc_min <- islocalminimum(chromatogramMatrix[, 3])
            boundary_left <- which(rt_loc_min[1:x_apex] == -1)
            if (length(boundary_left) > 0) {
              boundary_left <- boundary_left[length(boundary_left)]
            } else {
              boundary_left <- 1
            }
            boundary_right <- which(rt_loc_min[x_apex:length_chrom] == -1)
            if (length(boundary_right) > 0) {
              boundary_right <- boundary_right[1] - 1 + x_apex
            } else {
              boundary_right <- length_chrom
            }
            chrom <- cbind(retentionTime[chromatogramMatrix[boundary_left:boundary_right, 1]], chromatogramMatrix[boundary_left:boundary_right, 3])
            RT_detected <- chrom[which.min(abs(chrom[, 1] - rtCandidate)), 1]
            if (abs(RT_detected - rtCandidate) <= RTtolerance) {
              height <- max(chrom[, 2])
              ## R13C
              t1 <- boundary_left + scanNumberStart - 1
              t2 <- boundary_right + scanNumberStart - 1
              ##
              chromatogram_segment <- targetedIonPairing(spectraList, t1, t2, mzCandidate, massAccuracy, ionMassDifference, massAccuracy1.5)
              ##
              if (length(chromatogram_segment) > 0) {
                Int12C <- sum(chromatogram_segment[, 2])
                Int13C <- sum(chromatogram_segment[, 5])
                R13C <- Int13C/Int12C*100
              }
              if (boundary_left < boundary_right) {
                peakArea <- peakAreaCalculator(chrom[, 1], chrom[, 2])
              }
              ##
              c(x0[j], round(height, 0), round(peakArea, 0), round(R13C, 3))
            }
          }
        }))
      }
    }
    return(chromatography_undetected)
  }
  ##
  ##############################################################################
  ##
  massAccuracy1.5 <- 1.5*massAccuracy
  Lsamples3 <- dim(peakXcol)[2]
  Lsamples <- Lsamples3 - 3
  colnamesPeakXcol <- colnames(peakXcol)[4:Lsamples3]
  ##
  if (number_processing_threads == 1) {
    ##
    progressBARboundaries <- txtProgressBar(min = 0, max = Lsamples, initial = 0, style = 3)
    ##
    chromatography_undetected_list <- lapply(1:Lsamples, function(i) {
      setTxtProgressBar(progressBARboundaries, i)
      ##
      call_gapFillingCore(i)
    })
    ##
    close(progressBARboundaries)
    ##
  } else {
    ## Processing OS
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Windows") {
      ##
      clust <- makeCluster(number_processing_threads)
      clusterExport(clust, setdiff(ls(), c("clust", "Lsamples")), envir = environment())
      ##
      chromatography_undetected_list <- parLapply(clust, 1:Lsamples, function(i) {
        call_gapFillingCore(i)
      })
      ##
      stopCluster(clust)
      ##
      #########################################################################
      ##
    } else {
      ##
      chromatography_undetected_list <- mclapply(1:Lsamples, function(i) {
        call_gapFillingCore(i)
      }, mc.cores = number_processing_threads)
      ##
      closeAllConnections()
      ##
    }
  }
  ##
  names(chromatography_undetected_list) <- colnamesPeakXcol
  ##
  return(chromatography_undetected_list)
}
