IPA_targeted_xlsxAnalyzer <- function(spreadsheet) {
  ##
  checkpoint_parameter <- FALSE
  ##
  if (typeof(spreadsheet) == "list") {
    if (ncol(spreadsheet) >= 4) {
      PARAM <- cbind(spreadsheet[, 2], spreadsheet[, 4])
      checkpoint_parameter <- TRUE
      ##
    } else if (ncol(spreadsheet) == 2) {
      PARAM <- spreadsheet
      checkpoint_parameter <- TRUE
      ##
    } else {
      IPA_message("The `IPA_targeted` spreadsheet tab was not produced properly!")
    }
  } else if (typeof(spreadsheet) == "character") {
    if (length(spreadsheet) == 1) {
      if (file.exists(spreadsheet)) {
        PARAM <- readxl::read_xlsx(spreadsheet, sheet = "IPA_targeted")
        PARAM <- cbind(PARAM[, 2], PARAM[, 4])
        checkpoint_parameter <- TRUE
      } else {
        IPA_message("The `IPA_targeted` spreadsheet tab not found! It should be an Excel file with .xlsx extention!")
      }
    } else {
      IPA_message("The `IPA_targeted` spreadsheet tab was not produced properly!")
    }
  } else {
    IPA_message("The `IPA_targeted` spreadsheet tab was not produced properly!")
  }
  ##############################################################################
  if (checkpoint_parameter) {
    ##
    number_processing_threads <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0006'), 2])
    if (is.na(number_processing_threads)) {
      IPA_message("ERROR!!! Problem with PARAM0006!")
      checkpoint_parameter <- FALSE
    } else {
      if (number_processing_threads >= 1) {
        if ((number_processing_threads %% 1) != 0) {
          IPA_message("ERROR!!! Problem with PARAM0006! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      } else {
        IPA_message("ERROR!!! Problem with PARAM0006! This parameter should be at least 1 !")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    if (number_processing_threads > 1) {
      x_par <- which(PARAM[, 1] == 'PARAM_PAR')
      parallelizationMode <- gsub(" ", "", tolower(PARAM[x_par, 2]))
      ##
      if ((parallelizationMode == "peakmode") | (parallelizationMode == "samplemode")) {
        PARAM[x_par, 2] <- parallelizationMode
      } else {
        IPA_message("ERROR!!! Problem with 'PARAM_PAR'! This parameter should be `Sample Mode` or `Peak Mode`!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0007 <- which(PARAM[, 1] == 'PARAM0007')
    if (length(x0007) == 0) {
      IPA_message("ERROR!!! Problem with PARAM0007!")
      checkpoint_parameter <- FALSE
    } else {
      input_path_hrms <- PARAM[x0007, 2]
      input_path_hrms <- gsub("\\", "/", input_path_hrms, fixed = TRUE)
      PARAM[x0007, 2] <- input_path_hrms
      if (!dir.exists(input_path_hrms)) {
        IPA_message("ERROR!!! Problem with PARAM0007! Please make sure the full path is provided!")
        checkpoint_parameter <- FALSE
      }
      ##
      x0008 <- which(PARAM[, 1] == 'PARAM0008')
      if (is.na(PARAM[x0008, 2])) {
        IPA_message("ERROR!!! Problem with PARAM0008!")
        checkpoint_parameter <- FALSE
      } else {
        if (tolower(PARAM[x0008, 2]) == "all") {
          file_name_hrms <- dir(path = input_path_hrms)
          file_name_hrms <- file_name_hrms[grep(pattern = ".mzML$|.mzXML$|.CDF$", file_name_hrms, ignore.case = TRUE)]
          if (length(file_name_hrms) == 0) {
            IPA_message("ERROR!!! Problem with PARAM0008! No mzML/mzXML/CDF file was detected in the folder!")
          }
        } else {
          samples_string <- PARAM[x0008, 2]
          file_name_hrms <- strsplit(samples_string, ";")[[1]]
          ID <- sapply(1:length(file_name_hrms), function(i) {
            ID_name <- paste0(input_path_hrms, "/", file_name_hrms[i])
            as.numeric(file.exists(ID_name))
          })
          x_ID <- which(ID == 0)
          if (length(x_ID) > 0) {
            IPA_message("ERROR!!! Problem with PARAM0008! not detected the following file(s) (case sensitive even for file extensions):")
            for (i in x_ID) {
              IPA_message(file_name_hrms[i])
            }
            checkpoint_parameter <- FALSE
          }
        }
      }
    }
    ##
    PARAM0009 <- tolower(PARAM[which(PARAM[, 1] == 'PARAM0009'), 2])
    if (!(PARAM0009 == "y" | PARAM0009 == "yes" | PARAM0009 == "n" | PARAM0009 == "no")) {
      checkpoint_parameter <- FALSE
      IPA_message("Error!!! Problems with PARAM0009!")
    }
    ##
    PARAM_CCT <- tolower(PARAM[which(PARAM[, 1] == 'PARAM_CCT'), 2])
    if (!(PARAM_CCT == "y" | PARAM_CCT == "yes" | PARAM_CCT == "n" | PARAM_CCT == "no")) {
      checkpoint_parameter <- FALSE
      IPA_message("Error!!! Problems with 'PARAM_CCT'!")
    }
    ##
    x0010 <- which(PARAM[, 1] == 'PARAM0010')
    if (length(x0010) == 0) {
      IPA_message("ERROR!!! Problem with PARAM0010!")
      checkpoint_parameter <- FALSE
    } else {
      output_path <- gsub("\\", "/", PARAM[x0010, 2], fixed = TRUE)
      PARAM[x0010, 2] <- output_path
      if (!dir.exists(output_path)) {
        tryCatch(dir.create(output_path, recursive = TRUE), warning = function(w){warning("Problem with PARAM0010! R cannot create the folder!")})
        if (!dir.exists(output_path)) {
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    x0012 <- which(PARAM[, 1] == 'PARAM0012')
    if (length(x0012) == 0) {
      IPA_message("ERROR!!! Problem with PARAM0012!")
      checkpoint_parameter <- FALSE
    } else {
      ionMassDifference <- tryCatch(as.numeric(PARAM[x0012, 2]), error = function(e) {1.003354835336}, warning = function(w) {1.003354835336})     # Mass difference for isotopic pairs
      PARAM[x0012, 2] <- ionMassDifference
    }
    ##################### Chromatographic peak detection #######################
    massAccuracy <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0013'), 2])
    if (is.na(massAccuracy)) {
      IPA_message("ERROR!!! Problem with PARAM0013!")
      checkpoint_parameter <- FALSE
    } else {
      if (massAccuracy > 0.01) {
        IPA_message("ERROR!!! Problem with PARAM0013! Mass accuracy must be below `0.01 Da`")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0015 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0015'), 2])
    if (is.na(x0015)) {
      IPA_message("ERROR!!! Problem with PARAM0015! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (x0015 <= 0) {
        IPA_message("ERROR!!! Problem with PARAM0015! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0017 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0017'), 2])
    if (is.na(x0017)) {
      IPA_message("ERROR!!! Problem with PARAM0017! This value should be a positive number between 0-0.05 !")
      checkpoint_parameter <- FALSE
    } else {
      if (x0017 < 0 | x0017 > 0.1) {
        IPA_message("ERROR!!! Problem with PARAM0017! This value should be a positive number between 0-0.05 !")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0020 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0020'), 2])
    if (is.na(x0020)) {
      IPA_message("ERROR!!! Problem with PARAM0020! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (x0020 <= 0) {
        IPA_message("ERROR!!! Problem with PARAM0020! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      } else {
        if ((x0020 %% 1) != 0) {
          IPA_message("ERROR!!! Problem with PARAM0020! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    x0028 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0028'), 2])
    if (is.na(x0028)) {
      IPA_message("ERROR!!! Problem with PARAM0028! This parameter should be a positive integer greater than or equal to 11 !")
      checkpoint_parameter <- FALSE
    } else {
      if (x0028 < 0) {
        IPA_message("ERROR!!! Problem with PARAM0028! This parameter should be a positive integer greater than or equal to 11 !")
        checkpoint_parameter <- FALSE
      } else if (x0028 <= 11 && x0028 >= 1) {
        IPA_message("ERROR!!! Problem with PARAM0028! This parameter should be a positive integer greater than or equal to 11 !")
        checkpoint_parameter <- FALSE
      } else {
        if ((x0028 %% 1) != 0) {
          IPA_message("ERROR!!! Problem with PARAM0028! This parameter should be a positive integer greater than or equal to 11 !")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    mzCandidate <- tryCatch(eval(parse(text = paste0("c(", PARAM[which(PARAM[, 1] == 'PARAM_MZ'), 2], ")"))), error = function(e){NULL})
    rtCandidate <- tryCatch(eval(parse(text = paste0("c(", PARAM[which(PARAM[, 1] == 'PARAM_RT'), 2], ")"))), error = function(e){NULL})
    if ((length(mzCandidate) !=  length(rtCandidate)) | is.null(mzCandidate) | is.null(rtCandidate)) {
      checkpoint_parameter <- FALSE
      IPA_message("Error!!! Problems with `PARAM_MZ` and `PARAM_RT` ! mz and RT vectors do not have the same length!")
    }
  }
  ##############################################################################
  if (!checkpoint_parameter) {
    PARAM <- NULL
  }
  return(PARAM)
}
