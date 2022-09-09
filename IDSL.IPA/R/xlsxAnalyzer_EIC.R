xlsxAnalyzer_EIC <- function (spreadsheet) {
  ##
  checkpoint_parameter <- FALSE
  if (length(spreadsheet) >= 4) {
    if (typeof(spreadsheet) == "list") {
      PARAM <- cbind(spreadsheet[, 2], spreadsheet[, 4])
      checkpoint_parameter <- TRUE
    } else {
      print("The IPA input was not produced properly!")
    }
  } else if (length(spreadsheet) == 1) {
    if (typeof(spreadsheet) == "character") {
      if (file.exists(spreadsheet)) {
        spreadsheet_IPA <- readxl::read_xlsx(spreadsheet, sheet = 'IPA_targeted')
        PARAM <- cbind(spreadsheet_IPA[, 2], spreadsheet_IPA[, 4])
        checkpoint_parameter <- TRUE
      } else {
        print("The IPA spreadsheet not found! It should be an Excel file with .xlsx extention!")
      }
    } else {
      print("The IPA spreadsheet was not produced properly!")
    }
  } else {
    print("The IPA spreadsheet was not produced properly!")
  }
  if (checkpoint_parameter) {
    ############################################################################
    x0006 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0006'), 2])
    if (is.na(x0006)) {
      print("ERROR!!! Problem with PARAM0006!")
      checkpoint_parameter <- FALSE
    } else {
      if (x0006 >= 1) {
        if ((x0006 %% 1) != 0) {
          print("ERROR!!! Problem with PARAM0006! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      } else {
        print("ERROR!!! Problem with PARAM0006! This parameter should be at least 1 !")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0007 <- which(PARAM[, 1] == 'PARAM0007')
    if (length(x0007) == 0) {
      print("ERROR!!! Problem with PARAM0007!")
      checkpoint_parameter <- FALSE
    } else {
      address_hrms <- PARAM[x0007, 2]
      address_hrms <- gsub("\\", "/", address_hrms, fixed = TRUE)
      PARAM[x0007, 2] <- address_hrms
      if (!dir.exists(address_hrms)) {
        print("ERROR!!! Problem with PARAM0007! Please make sure the full path is provided!")
        checkpoint_parameter <- FALSE
      }
      ##
      x0008 <- which(PARAM[, 1] == 'PARAM0008')
      if (is.na(PARAM[x0008, 2])) {
        print("ERROR!!! Problem with PARAM0008!")
        checkpoint_parameter <- FALSE
      } else {
        if (tolower(PARAM[x0008, 2]) != "all") {
          samples_string <- PARAM[x0008, 2]
          name <- strsplit(samples_string, ";")[[1]] # files used as reference m/z-RT
          ID <- sapply(1:length(name), function(i) {
            ID_name <- paste0(address_hrms, "/", name[i])
            as.numeric(file.exists(ID_name))
          })
          x_ID <- which(ID == 0)
          if (length(x_ID) > 0) {
            print("ERROR!!! Problem with PARAM0008! not detected the following file(s) (case sensitive even for file extensions):")
            for (i in 1:length(x_ID)) {
              print(name[x_ID[i]])
            }
            checkpoint_parameter <- FALSE
          }
        }
        ##
        if (tolower(PARAM[x0008, 2]) == "all") {
          x0009 <- PARAM[which(PARAM[, 1] == 'PARAM0009'), 2]
          if (is.na(x0009)) {
            print("ERROR!!! Problem with PARAM0009!")
            checkpoint_parameter <- FALSE
          } else {
            if (tolower(x0009) == "mzml" | tolower(x0009) == "mzxml" | tolower(x0009) == "cdf") {
              cat("")
            } else {
              print("ERROR!!! Problem with PARAM0009! HRMS data are incompatible!")
              checkpoint_parameter <- FALSE
            }
          }
        }
      }
    }
    ##
    x0010 <- which(PARAM[, 1] == 'PARAM0010')
    if (length(x0010) == 0) {
      print("ERROR!!! Problem with PARAM0010!")
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
      print("ERROR!!! Problem with PARAM0012!")
      checkpoint_parameter <- FALSE
    } else {
      massDifferenceIsotopes <- tryCatch(as.numeric(PARAM[x0012, 2]), error = function(e) {1.003354835336}, warning = function(w) {1.003354835336})     # Mass difference for isotopic pairs
      PARAM[x0012, 2] <- massDifferenceIsotopes
    }
    ##################### Chromatographic peak detection #######################
    x0013 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0013'), 2])
    if (is.na(x0013)) {
      print("ERROR!!! Problem with PARAM0013!")
      checkpoint_parameter <- FALSE
    } else {
      if (x0013 <= 0) {
        print("ERROR!!! Problem with PARAM0013! This value should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0015 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0015'), 2])
    if (is.na(x0015)) {
      print("ERROR!!! Problem with PARAM0015! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (x0015 <= 0) {
        print("ERROR!!! Problem with PARAM0015! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0017 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0017'), 2])
    if (is.na(x0017)) {
      print("ERROR!!! Problem with PARAM0017! This value should be a positive number between 0-0.05 !")
      checkpoint_parameter <- FALSE
    } else {
      if (x0017 < 0 | x0017 > 0.1) {
        print("ERROR!!! Problem with PARAM0017! This value should be a positive number between 0-0.05 !")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0020 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0020'), 2])
    if (is.na(x0020)) {
      print("ERROR!!! Problem with PARAM0020! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (x0020 <= 0) {
        print("ERROR!!! Problem with PARAM0020! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      } else {
        if ((x0020 %% 1) != 0) {
          print("ERROR!!! Problem with PARAM0020! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    x0028 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0028'), 2])
    if (is.na(x0028)) {
      print("ERROR!!! Problem with PARAM0028! This parameter should be a positive integer greater than or equal to 11 !")
      checkpoint_parameter <- FALSE
    } else {
      if (x0028 < 0) {
        print("ERROR!!! Problem with PARAM0028! This parameter should be a positive integer greater than or equal to 11 !")
        checkpoint_parameter <- FALSE
      } else if (x0028 <= 11 && x0028 >= 1) {
        print("ERROR!!! Problem with PARAM0028! This parameter should be a positive integer greater than or equal to 11 !")
        checkpoint_parameter <- FALSE
      } else {
        if ((x0028 %% 1) != 0) {
          print("ERROR!!! Problem with PARAM0028! This parameter should be a positive integer greater than or equal to 11 !")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    mzCandidate <- tryCatch(eval(parse(text = paste0("c(", PARAM[which(PARAM[, 1] == 'PARAM_MZ'), 2], ")"))), error = function(e){NULL})
    rtCandidate <- tryCatch(eval(parse(text = paste0("c(", PARAM[which(PARAM[, 1] == 'PARAM_RT'), 2], ")"))), error = function(e){NULL})
    if ((length(mzCandidate) !=  length(rtCandidate)) | is.null(mzCandidate) | is.null(rtCandidate)) {
      checkpoint_parameter <- FALSE
      print("Error!!! Problems with PARAM_MZ and PARAM_RT ! mz and RT vectors do not have the same length!")
    }
    ##
    ipa_eic_tar <- tolower(PARAM[which(PARAM[, 1] == 'PARAM_EIC'), 2])
    if (!(ipa_eic_tar == "y" | ipa_eic_tar == "yes" | ipa_eic_tar == "n" | ipa_eic_tar == "no")) {
      checkpoint_parameter <- FALSE
      print("Error!!! Problems with PARAM_EIC !")
    }
    ##
    ipa_tab_tar <- tolower(PARAM[which(PARAM[, 1] == 'PARAM_CCT'), 2])
    if (!(ipa_tab_tar == "y" | ipa_tab_tar == "yes" | ipa_tab_tar == "n" | ipa_tab_tar == "no")) {
      checkpoint_parameter <- FALSE
      print("Error!!! Problems with PARAM_EIC !")
    }
    ##
  }
  ##############################################################################
  if (checkpoint_parameter == FALSE) {
    print("Please visit   https://ipa.idsl.me    for instructions!")
    PARAM <- c()
  }
  return(PARAM)
}
