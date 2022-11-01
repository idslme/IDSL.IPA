gapFillingCore <- function(input_path_hrms, peakXcol, massAccuracy, RTtolerance, scanTolerance,
                           retentionTimeCorrectionCheck = FALSE, listCorrectedRTpeaklists = NULL,
                           inputPathPeaklist = NULL, ionMassDifference = 1.003354835336, number_processing_threads = 1) {
  ##
  call_chromatography_undetected_list <- function(i) {
    ##
    chromatography_undetected <- NULL
    x_0 <- which(peakXcol[, (i + 2)] == 0)
    L_x_0 <- length(x_0)
    if (L_x_0 > 0) {
      ##
      iHRMSfileName <- colnamesPeakXcol[i]
      mzXcol <- peakXcol[x_0, 1]
      ## To back calculate the RT ##
      retentionTimeNAcheck <- TRUE
      if (retentionTimeCorrectionCheck) {
        iPeaklistFileName <- paste0("peaklist_", iHRMSfileName, ".Rdata")
        corrected_RTi <- listCorrectedRTpeaklists[[iPeaklistFileName]]
        ##
        if (!is.na(corrected_RTi[1])) {
          undeteced_RT <- peakXcol[x_0, 2]
          uncorrected_RTi <- loadRdata(paste0(inputPathPeaklist, "/", iPeaklistFileName))[, 3]
          idf <- data.frame(uncoRT = uncorrected_RTi, coRT = corrected_RTi)
          rtmodel <- lm(coRT ~ poly(uncoRT, 3), idf)
          new.df <- data.frame(uncoRT = undeteced_RT) # predict uncorrected RTs
          RT_uncorrected_undeteced <- predict(rtmodel, new.df)
        } else {
          retentionTimeNAcheck <- FALSE
        }
      } else {
        RT_uncorrected_undeteced <- peakXcol[x_0, 2]
      }
      ##
      if (retentionTimeNAcheck) {
        ##
        outputer <- IPA_MSdeconvoluter(input_path_hrms, iHRMSfileName)
        spectraList <- outputer[["spectraList"]]
        RetentionTime <- outputer[["retentionTime"]]
        nRT <- length(RetentionTime)
        ##
        chromatography_undetected <- do.call(rbind, lapply(1:L_x_0, function(j) {
          chromatography_undetected_row <- NULL
          pa <- 0
          R13C <- 0
          mzCandidate <- mzXcol[j]
          rtCandidate <- RT_uncorrected_undeteced[j]
          sn_apex <- which.min(abs(rtCandidate - RetentionTime)) # scan number at apex
          ScanNumberStart <- max(c((sn_apex - scanTolerance), 1))
          ScanNumberEnd <- min(c(nRT, (sn_apex + scanTolerance)))
          spectraListXIC <- spectraList[ScanNumberStart:ScanNumberEnd]
          chrom_builder <- XIC(spectraListXIC, scanNumberStart = ScanNumberStart, mzCandidate, massAccuracy)
          length_chrom <- dim(chrom_builder)[1]
          x_apex <- which(chrom_builder[, 1] == sn_apex)
          rt_loc_min <- islocalminimum(chrom_builder[, 3])
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
          chrom <- cbind(RetentionTime[chrom_builder[boundary_left:boundary_right, 1]], chrom_builder[boundary_left:boundary_right, 3])
          RT_detected <- chrom[which.min(abs(chrom[, 1] - rtCandidate)), 1]
          if (abs(RT_detected - rtCandidate) <= RTtolerance) {
            height <- max(chrom[, 2])
            ## R13C
            t1 <- boundary_left + ScanNumberStart - 1
            t2 <- boundary_right + ScanNumberStart - 1
            chromatogram_segment <- do.call(rbind, lapply(t1:t2, function(t) {
              Spec_ScN_j <- NULL
              Spec <- spectraList[[t]]
              if (length(Spec) > 0) {
                x_mz1 <- which(abs(Spec[, 1] - mzCandidate) <= massAccuracy)
                if (length(x_mz1) > 0) {
                  x_mz2 <- which(abs(Spec[, 1] - (ionMassDifference + mzCandidate)) <= massAccuracy1.5)
                  if (length(x_mz2) > 0) {
                    if (length(x_mz1) > 1) {
                      x_min <- which.min(abs(Spec[x_mz1, 1] - mzCandidate))
                      x_mz1 <- x_mz1[x_min]
                    }
                    if (length(x_mz2) > 1) {
                      x_min <- which.min(abs(Spec[x_mz2, 1] - (ionMassDifference + mzCandidate)))
                      x_mz2 <- x_mz2[x_min]
                    }
                    if (Spec[x_mz1, 2] >= Spec[x_mz2, 2]) {
                      Spec_ScN_j <- c(Spec[x_mz1, 2], Spec[x_mz2, 2])
                    }
                  }
                }
              }
              Spec_ScN_j
            }))
            if (length(chromatogram_segment) > 0) {
              Int12C <- sum(chromatogram_segment[, 1])
              Int13C <- sum(chromatogram_segment[, 2])
              R13C <- Int13C/Int12C*100
            }
            if (boundary_left != boundary_right) {
              pa <- peakAreaCalculator(chrom[, 1], chrom[, 2])
            }
            chromatography_undetected_row <- c(x_0[j], round(height, 0), round(pa, 0), round(R13C, 3))
          }
          chromatography_undetected_row
        }))
      }
    }
    return(chromatography_undetected)
  }
  ##
  ##############################################################################
  ##
  massAccuracy1.5 <- 1.5*massAccuracy
  L_samples2 <- dim(peakXcol)[2]
  L_samples <- L_samples2 - 2
  colnamesPeakXcol <- colnames(peakXcol)[3:L_samples2]
  ##
  if (number_processing_threads == 1) {
    ##
    progressBARboundaries <- txtProgressBar(min = 0, max = L_samples, initial = 0, style = 3)
    ##
    chromatography_undetected_list <- lapply(1:L_samples, function(i) {
      setTxtProgressBar(progressBARboundaries, i)
      ##
      call_chromatography_undetected_list(i)
    })
    ##
  } else {
    ## Processing OS
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Linux") {
      ##
      chromatography_undetected_list <- mclapply(1:L_samples, function(i) {
        call_chromatography_undetected_list(i)
      }, mc.cores = number_processing_threads)
      ##
      closeAllConnections()
      ##
    } else if (osType == "Windows") {
      ##
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      chromatography_undetected_list <- foreach(i = 1:L_samples, .verbose = FALSE) %dopar% {
        call_chromatography_undetected_list(i)
      }
      stopCluster(clust)
      ##
    }
  }
  ##
  names(chromatography_undetected_list) <- colnamesPeakXcol
  ##
  return(chromatography_undetected_list)
}
