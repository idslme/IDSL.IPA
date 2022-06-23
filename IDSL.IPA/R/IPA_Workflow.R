IPA_Workflow <- function(spreadsheet) {
  initiation_time <- Sys.time()
  ##
  gc()
  closeAllConnections()
  ##
  PARAM <- IPA_xlsxAnalyzer(spreadsheet)
  if (length(PARAM) > 0) {
    if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0001'), 2]) == "yes") {
      IPA_PeakAnalyzer(PARAM)
    }
    if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0002'), 2]) == "yes") {
      IPA_PeakAlignment(PARAM)
    } else {
      x0004 <- PARAM[which(PARAM[, 1] == 'PARAM0004'), 2]
      x0045 <- PARAM[which(PARAM[, 1] == 'PARAM0045'), 2]
      if (!is.na(x0004) & !is.na(x0045)) {
        if (tolower(x0004) == "yes" & tolower(x0045) == "yes") {
          IPA_PeakAlignment(PARAM)
        }
      }
    }
    if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0003'), 2]) == "yes") {
      IPA_GapFiller(PARAM)
    }
    if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0004'), 2]) == "yes") {
      if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0046'), 2]) == "yes") {
        IPA_CompoundsAnnotation(PARAM)
      }
      if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0047'), 2]) == "yes") {
        IPA_PeaklistAnnotation(PARAM)
      }
    }
    ##
    if (tolower(PARAM[which(PARAM[, 1] == 'PARAM0005'), 2]) == "yes") {
      PARAM_targeted <- xlsxAnalyzer_EIC(spreadsheet)
      ##
      mzCandidate <- tryCatch(eval(parse(text = paste0("c(", PARAM_targeted[which(PARAM_targeted[, 1] == 'PARAM_MZ'), 2], ")"))), error = function(e){NULL})
      rtCandidate <- tryCatch(eval(parse(text = paste0("c(", PARAM_targeted[which(PARAM_targeted[, 1] == 'PARAM_RT'), 2], ")"))), error = function(e){NULL})
      #
      if (is.null(mzCandidate) | is.null(rtCandidate)) {
        stop("ERROR!!! Incorrect 'PARAM_MZ' or 'PARAM_RT'")
      }
      ##
      ipa_eic_tar <- tolower(PARAM_targeted[which(PARAM_targeted[, 1] == 'PARAM_EIC'), 2])
      if (ipa_eic_tar == "y" | ipa_eic_tar == "yes") {
        exportEIC_TorF <- TRUE
      } else if (ipa_eic_tar == "n" | ipa_eic_tar == "no") {
        exportEIC_TorF <- FALSE
      }
      ##
      ipa_tab_tar <- tolower(PARAM_targeted[which(PARAM_targeted[, 1] == 'PARAM_CCT'), 2])
      if (ipa_tab_tar == "y" | ipa_tab_tar == "yes") {
        exportTable_TorF <- TRUE
      } else if (ipa_tab_tar == "n" | ipa_tab_tar == "no") {
        exportTable_TorF <- FALSE
      }
      ##
      IPA_TargetedTable <- IPA_TargetedAnalysis(spreadsheet, mzCandidate, rtCandidate, exportEIC = exportEIC_TorF, exportTable = exportTable_TorF)
      ##
      if (exportTable_TorF) {
        output_path <- PARAM_targeted[which(PARAM_targeted[, 1] == 'PARAM0010'), 2]
        save(IPA_TargetedTable, file = paste0(output_path, "/IPA_TargetedTable.Rdata"))
        write.csv(IPA_TargetedTable, file = paste0(output_path, "/IPA_TargetedTable.csv"))
      }
    }
    ##
    print("Completed IDSL.IPA computations successfully!")
    required_time <- Sys.time() - initiation_time
    print(required_time)
  }
  ##
  gc()
  closeAllConnections()
  ##
}
