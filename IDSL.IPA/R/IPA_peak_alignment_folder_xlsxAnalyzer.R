IPA_peak_alignment_folder_xlsxAnalyzer <- function(PARAM, PARAM_ID, checkpoint_parameter, correctedRTcheck = FALSE, CSAcheck = FALSE, allowedVerbose = TRUE) {
  ##
  xPARAM_ID <- which(PARAM[, 1] == PARAM_ID)
  if (length(xPARAM_ID) == 0) {
    if (allowedVerbose) {IPA_message(paste0("ERROR!!! Problem with ", PARAM_ID, "!"))}
    checkpoint_parameter <- FALSE
  }
  if (checkpoint_parameter) {
    peak_alignment_folder <- as.character(PARAM[xPARAM_ID, 2])
    if (tolower(peak_alignment_folder) != "na") {
      peak_alignment_folder <- gsub("\\", "/", peak_alignment_folder, fixed = TRUE)
      PARAM[xPARAM_ID, 2] <- peak_alignment_folder
      ##
      if (!dir.exists(peak_alignment_folder)) {
        if (allowedVerbose) {IPA_message(paste0("ERROR!!! Problem with ", PARAM_ID, "! Please make sure the full path is provided!"))}
        checkpoint_parameter <- FALSE
      } else {
        peakXcol_FN <- paste0(peak_alignment_folder, "/peakXcol.Rdata")
        ##
        if (!file.exists(peakXcol_FN)) {
          if (allowedVerbose) {IPA_message(paste0("ERROR!!! Problem with ", PARAM_ID, "! `peakXcol.Rdata` was not detected in the `peak_alignment` folder (", PARAM_ID, ")!"))}
          checkpoint_parameter <- FALSE
        }
        ##
        peak_height_FN <- paste0(peak_alignment_folder, "/peak_height.Rdata")
        ##
        if (!file.exists(peak_height_FN)) {
          if (allowedVerbose) {IPA_message(paste0("ERROR!!! Problem with ", PARAM_ID, "! `peak_height.Rdata` was not detected in the `peak_alignment` folder (", PARAM_ID, ")!"))}
          checkpoint_parameter <- FALSE
        }
        ##
        peak_R13C_FN <- paste0(peak_alignment_folder, "/peak_R13C.Rdata")
        ##
        if (!file.exists(peak_R13C_FN)) {
          if (allowedVerbose) {IPA_message(paste0("ERROR!!! Problem with ", PARAM_ID, "! `peak_R13C.Rdata` was not detected in the `peak_alignment` folder (", PARAM_ID, ")!"))}
          checkpoint_parameter <- FALSE
        }
        ##
        if (CSAcheck) {
          alignedPeakHeightTableCorrelationList_FN <- paste0(peak_alignment_folder, "/alignedPeakHeightTableCorrelationList.Rdata")
          ##
          if (!file.exists(alignedPeakHeightTableCorrelationList_FN)) {
            if (allowedVerbose) {IPA_message(paste0("ERROR!!! Problem with ", PARAM_ID, "! `alignedPeakHeightTableCorrelationList.Rdata` was not detected in the `peak_alignment` folder (", PARAM_ID, ")!"))}
            checkpoint_parameter <- FALSE
          }
        }
        ##
        if (correctedRTcheck) {
          listCorrectedRTpeaklists_FN <- paste0(peak_alignment_folder, "/listCorrectedRTpeaklists.Rdata")
          ##
          if (!file.exists(listCorrectedRTpeaklists_FN)) {
            if (allowedVerbose) {IPA_message(paste0("ERROR!!! Problem with ", PARAM_ID, "! `listCorrectedRTpeaklists.Rdata` was not detected in the `peak_alignment` folder (", PARAM_ID, ")!"))}
            checkpoint_parameter <- FALSE
          }
        }
      }
    } else {
      if (CSAcheck | correctedRTcheck) {
        if (allowedVerbose) {IPA_message(paste0("ERROR!!! Problem with ", PARAM_ID, "! Please make sure the full path is provided!"))}
        checkpoint_parameter <- FALSE
      }
    }
  }
  ##
  listAlignmentFolderCheck <- list(PARAM, checkpoint_parameter)
  ##
  return(listAlignmentFolderCheck)
}
