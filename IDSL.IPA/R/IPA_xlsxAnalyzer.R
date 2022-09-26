IPA_xlsxAnalyzer <- function(spreadsheet) {
  ##
  print("Initiated testing the IPA spreadsheet consistency!")
  ##
  RT_correction_check <- function(checkpoint_parameter, PARAM) {
    ##
    x0029 <- which(PARAM[, 1] == 'PARAM0029')
    if (is.na(x0029)) {
      print("ERROR!!! Problem with PARAM0029!")
      checkpoint_parameter <- FALSE
    } else {
      FA0029 <- tolower(PARAM[x0029, 2])
      if ((FA0029 == "yes") | (FA0029 == "no")) {
        PARAM[x0029, 2] <- FA0029
      } else {
        print("ERROR!!! Problem with PARAM0029!")
        checkpoint_parameter <- FALSE
      }
    }
    if (FA0029 == "yes") {
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
      }
      ##
      reference_samples_string <- PARAM[which(PARAM[, 1] == 'PARAM0030'), 2] ## x0030
      if (is.na(reference_samples_string)) {
        print("ERROR!!! Problem with PARAM0030!")
        checkpoint_parameter <- FALSE
      } else {
        Ref_name <- strsplit(reference_samples_string, ";")[[1]] # files used as reference m/z-RT
        Ref_ID <- sapply(1:length(Ref_name), function(i) {
          ID_name <- paste0(address_hrms, "/", Ref_name[i])
          as.numeric(file.exists(ID_name))
        })
        x_Ref_ID <- which(Ref_ID == 0)
        if (length(x_Ref_ID) > 0) {
          print("ERROR!!! Problem with PARAM0030! not detected the following file(s) (case sensitive even for file extensions):")
          for (i in 1:length(x_Ref_ID)) {
            print(Ref_name[x_Ref_ID[i]])
          }
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0031 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0031'), 2])
      if (is.na(x0031)) {
        print("ERROR!!! Problem with PARAM0031!")
        checkpoint_parameter <- FALSE
      }
      ##
      x0032 <- PARAM[which(PARAM[, 1] == 'PARAM0032'), 2]
      if (is.na(x0032)) {
        print("ERROR!!! Problem with PARAM0032!")
        checkpoint_parameter <- FALSE
      } else {
        if (!(gsub(" ", "", tolower(x0032)) == "polynomial" | gsub(" ", "", tolower(x0032)) == "retentionindex")) {
          print("ERROR!!! Problem with PARAM0032!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      if (!is.na(x0032)) {
        if (gsub(" ", "", tolower(x0032)) == "retentionindex") {
          ##
          x0033 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0033'), 2])
          if (is.na(x0033)) {
            print("ERROR!!! Problem with PARAM0033! This parameter should be a positive integer smaller than 11 !")
            checkpoint_parameter <- FALSE
          } else {
            if (x0033 <= 0) {
              print("ERROR!!! Problem with PARAM0033! This parameter should be a positive integer smaller than 11 !")
              checkpoint_parameter <- FALSE
            } else if (x0033 > 11) {
              print("ERROR!!! Problem with PARAM0033! This parameter should be a positive integer smaller than 11 !")
              checkpoint_parameter <- FALSE
            } else {
              if ((x0033 %% 1) != 0) {
                print("ERROR!!! Problem with PARAM0033! This parameter should be a positive integer smaller than 11 !")
                checkpoint_parameter <- FALSE
              }
            }
          }
          ##
        }
        ##
        if (gsub(" ", "", tolower(x0032)) == "polynomial") {
          ##
          x0034 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0034'), 2])
          if (is.na(x0034)) {
            print("ERROR!!! Problem with PARAM0034! This parameter should be a positive integer smaller than 11 !")
            checkpoint_parameter <- FALSE
          } else {
            if (x0034 <= 0) {
              print("ERROR!!! Problem with PARAM0034! This parameter should be a positive integer smaller than 11 !")
              checkpoint_parameter <- FALSE
            } else if (x0034 > 11) {
              print("ERROR!!! Problem with PARAM0034! This parameter should be a positive integer smaller than 11 !")
              checkpoint_parameter <- FALSE
            } else {
              if ((x0034 %% 1) != 0) {
                print("ERROR!!! Problem with PARAM0034! This parameter should be a positive integer smaller than 11 !")
                checkpoint_parameter <- FALSE
              }
            }
          }
          ##
        }
      }
    }
    ####
    x0035 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0035'), 2])
    if (is.na(x0035)) {
      print("ERROR!!! Problem with PARAM0035! This parameter should be greater than 0 !")
      checkpoint_parameter <- FALSE
    } else {
      if (x0035 <= 0) {
        print("ERROR!!! Problem with PARAM0035! This parameter should be greater than 0 !")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0036 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0036'), 2])
    if (is.na(x0036)) {
      print("ERROR!!! Problem with PARAM0036! This parameter should be greater than 0 !")
      checkpoint_parameter <- FALSE
    } else {
      if (x0036 <= 0) {
        print("ERROR!!! Problem with PARAM0036! This parameter should be greater than 0 !")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0037 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0037'), 2])
    if (is.na(x0037)) {
      print("ERROR!!! Problem with PARAM0037! This parameter should be a positive integer greater than or equal to 2 !")
      checkpoint_parameter <- FALSE
    } else {
      if (x0037 < 1) {
        print("ERROR!!! Problem with PARAM0037! This parameter should be a positive integer greater than or equal to 2 !")
        checkpoint_parameter <- FALSE
      } else {
        if ((x0037 %% 1) != 0) {
          print("ERROR!!! Problem with PARAM0037! This parameter should be a positive integer greater than or equal to 2 !")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    return(list(checkpoint_parameter, PARAM))
  }
  ##############################################################################
  checkpoint_parameter <- FALSE
  if (length(spreadsheet) >= 4) {
    if (typeof(spreadsheet) == "list") {
      PARAM <- cbind(spreadsheet[, 2], spreadsheet[, 4])
      checkpoint_parameter <- TRUE
    } else {
      print("The IPA spreadsheet was not produced properly!")
    }
  } else if (length(spreadsheet) == 1) {
    if (typeof(spreadsheet) == "character") {
      if (file.exists(spreadsheet)){
        spreadsheet_IPA <- readxl::read_xlsx(spreadsheet, sheet = 'parameters')
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
  if (checkpoint_parameter == TRUE) {
    ##################### Global parameters ######################################
    x0001 <- PARAM[which(PARAM[, 1] == 'PARAM0001'), 2]
    if (length(x0001) == 0) {
      print("ERROR!!! Problem with PARAM0001!")
      checkpoint_parameter <- FALSE
      x0001 <- 0
    } else {
      if (!(tolower(x0001) == "yes" | tolower(x0001) == "no")) {
        print("ERROR!!! Problem with PARAM0001!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0002 <- PARAM[which(PARAM[, 1] == 'PARAM0002'), 2]
    if (length(x0002) == 0) {
      print("ERROR!!! Problem with PARAM0002!")
      checkpoint_parameter <- FALSE
      x0002 <- 0
    } else {
      if (!(tolower(x0002) == "yes" | tolower(x0002) == "no")) {
        print("ERROR!!! Problem with PARAM0002!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0003 <- PARAM[which(PARAM[, 1] == 'PARAM0003'), 2]
    if (length(x0003) == 0) {
      print("ERROR!!! Problem with PARAM0003!")
      checkpoint_parameter <- FALSE
      x0003 <- 0
    } else {
      if (!(tolower(x0003) == "yes" | tolower(x0003) == "no")) {
        print("ERROR!!! Problem with PARAM0003!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0004 <- PARAM[which(PARAM[, 1] == 'PARAM0004'), 2]
    if (length(x0004) == 0) {
      print("ERROR!!! Problem with PARAM0004!")
      checkpoint_parameter <- FALSE
      x0004 <- 0
    } else {
      if (!(tolower(x0004) == "yes" | tolower(x0004) == "no")) {
        print("ERROR!!! Problem with PARAM0004!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0005 <- PARAM[which(PARAM[, 1] == 'PARAM0005'), 2]
    if (length(x0005) == 0) {
      print("ERROR!!! Problem with PARAM0005!")
      checkpoint_parameter <- FALSE
    } else {
      if (tolower(x0005) == "yes" | tolower(x0005) == "no") {
        if (tolower(x0005) == "yes") {
          PARAM_targeted <- xlsxAnalyzer_EIC(spreadsheet)
          if (length(PARAM_targeted) == 0) {
            checkpoint_parameter <- FALSE
          }
        }
      } else {
        print("ERROR!!! Problem with PARAM0005!")
        checkpoint_parameter <- FALSE
      }
    }
    # print("WARNING!!! IPA Targeted Analysis was selected in PARAM0005! You may use the targeted analysis only with the 'IPA_TargetedAnalysis' module!")
    ##
    x0006 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0006'), 2])
    if (length(x0006) == 0) {
      print("ERROR!!! Problem with PARAM0006! This parameter should be a positive integer!")
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
    ################################## Data ####################################
    hrms_address_needed <- 0
    if (tolower(x0001) == "yes" | tolower(x0003) == "yes") {
      hrms_address_needed <- 1
    }
    if (tolower(x0002) == "yes") {
      x0029 <- PARAM[which(PARAM[, 1] == 'PARAM0029'), 2]
      if (length(x0029) > 0) {
        if (tolower(x0029) == "yes") {
          hrms_address_needed <- 1
        }
      }
    }
    if (hrms_address_needed == 1) {
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
      }
      ##
      x0008 <- which(PARAM[, 1] == 'PARAM0008')
      if (is.na(PARAM[x0008, 2])) {
        print("ERROR!!! Problem with PARAM0008!")
        checkpoint_parameter <- FALSE
      } else {
        if (tolower(PARAM[x0008, 2]) != "all") {
          samples_string <- PARAM[x0008, 2]
          name <- strsplit(samples_string, ";")[[1]] # files used to detect reference m/z-RT peaks
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
            if (!(tolower(x0009) == "mzml" | tolower(x0009) == "mzxml" | tolower(x0009) == "cdf")) {
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
    ######################### Peaklist production ################################
    ##############################################################################
    ####### Pairing 12C/13C isotopologues in individual chromatogram scans #######
    if (tolower(x0001) == "yes") {
      x0011 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0011'), 2])
      if (is.na(x0011)) {
        print("ERROR!!! Problem with PARAM0011! This value should be a number greater than or equal to 0 !")
        checkpoint_parameter <- FALSE
      } else {
        if (x0011 < 0) {
          print("ERROR!!! Problem with PARAM0011! This value should be a number greater than or equal to 0 !")
          checkpoint_parameter <- FALSE
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
      #################### Chromatographic peak detection ########################
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
      x0014 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0014'), 2])
      if (is.na(x0014)) {
        print("ERROR!!! Problem with PARAM0014! This value should be a positive number!")
        checkpoint_parameter <- FALSE
      } else {
        if (x0014 < 0) {
          print("ERROR!!! Problem with PARAM0014! This value should be a positive number!")
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
      # x0016 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0016'), 2])
      # if (is.na(x0016)) {
      #   print("ERROR!!! Problem with PARAM0016! This parameter should be a positive integer!")
      #   checkpoint_parameter <- FALSE
      # } else {
      #   if (x0016 <= 0) {
      #     print("ERROR!!! Problem with PARAM0016! This parameter should be a positive integer!")
      #     checkpoint_parameter <- FALSE
      #   }
      # }
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
      x0018 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0018'), 2])
      if (is.na(x0018)) {
        print("ERROR!!! Problem with PARAM0018! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      } else {
        if (x0018 < 0) {
          print("ERROR!!! Problem with PARAM0018! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        } else {
          if ((x0018 %% 1) != 0) {
            print("ERROR!!! Problem with PARAM0018! This parameter should be a positive integer!")
            checkpoint_parameter <- FALSE
          }
        }
      }
      ############### Peak quality measurement and data reduction ###################
      x0019 <- PARAM[which(PARAM[, 1] == 'PARAM0019'), 2]
      if (is.na(x0019)) {
        print("ERROR!!! Problem with PARAM0019!")
        checkpoint_parameter <- FALSE
      } else {
        if (!(tolower(x0019) == "yes" | tolower(x0019) == "no")) {
          print("ERROR!!! Problem with PARAM0019!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      if (!is.na(x0019)) {
        if (tolower(x0019) == "yes") {
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
        }
      }
      ##
      x0021 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0021'), 2])
      if (is.na(x0021)) {
        print("ERROR!!! Problem with PARAM0021! This value should be a number greater than or equal to 0 !")
        checkpoint_parameter <- FALSE
      } else {
        if (x0021 < 0) {
          print("ERROR!!! Problem with PARAM0021! This value should be a number greater than or equal to 0 !")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0022 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0022'), 2])
      if (is.na(x0022)) {
        print("ERROR!!! Problem with PARAM0022!")
        checkpoint_parameter <- FALSE
      }
      ##
      x0023 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0023'), 2])
      if (is.na(x0023)) {
        print("ERROR!!! Problem with PARAM0023! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      } else {
        if (x0023 < 0) {
          print("ERROR!!! Problem with PARAM0023! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        } else {
          if ((x0023 %% 1) != 0) {
            print("ERROR!!! Problem with PARAM0023! This parameter should be a positive integer!")
            checkpoint_parameter <- FALSE
          }
        }
      }
      ##
      x0024 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0024'), 2])
      if (is.na(x0024)) {
        print("ERROR!!! Problem with PARAM0024!")
        checkpoint_parameter <- FALSE
      }
      ##
      x0025 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0025'), 2])
      if (is.na(x0025)) {
        print("ERROR!!! Problem with PARAM0025!")
        checkpoint_parameter <- FALSE
      }
      ##
      x0026 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0026'), 2])
      if (is.na(x0026)) {
        print("ERROR!!! Problem with PARAM0026!")
        checkpoint_parameter <- FALSE
      }
      ##
      x0027 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0027'), 2])
      if (is.na(x0027)) {
        print("ERROR!!! Problem with PARAM0027! This parameter should be greater than 1 !")
        checkpoint_parameter <- FALSE
      } else {
        if (x0027 < 1) {
          print("ERROR!!! Problem with PARAM0027! This parameter should be greater than 1 !")
          checkpoint_parameter <- FALSE
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
    }
    ############## RT correction and peak alignment table production #############
    if (tolower(x0002) == "yes") {
      RT_correction_checkList <- RT_correction_check(checkpoint_parameter, PARAM)
      checkpoint_parameter <- RT_correction_checkList[[1]]
      PARAM <- RT_correction_checkList[[2]]
    }
    ############################## Gap-filling ###################################
    gapfilling_check <- function(checkpoint_parameter, PARAM) {
      x0038 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0038'), 2])
      if (is.na(x0038)) {
        print("ERROR!!! Problem with PARAM0038! This parameter should be greater than 0 !")
        checkpoint_parameter <- FALSE
      } else {
        if (x0038 <= 0) {
          print("ERROR!!! Problem with PARAM0038! This parameter should be greater than 0 !")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0039 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0039'), 2])
      if (is.na(x0039)) {
        print("ERROR!!! Problem with PARAM0039! This parameter should be greater than 0 !")
        checkpoint_parameter <- FALSE
      } else {
        if (x0039 <= 0) {
          print("ERROR!!! Problem with PARAM0039! This parameter should be greater than 0 !")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0040 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0040'), 2])
      if (is.na(x0040)) {
        print("ERROR!!! Problem with PARAM0040! This parameter should be a positive integer!")
        checkpoint_parameter <- FALSE
      } else {
        if (x0040 < 0) {
          print("ERROR!!! Problem with PARAM0040! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        } else {
          if ((x0040 %% 1) != 0) {
            print("ERROR!!! Problem with PARAM0040! This parameter should be a positive integer!")
            checkpoint_parameter <- FALSE
          }
        }
      }
      return(checkpoint_parameter)
    }
    if (tolower(x0003) == "yes") {
      checkpoint_parameter <- gapfilling_check(checkpoint_parameter, PARAM)
    }
    ##############################################################################
    if (tolower(x0004) == "yes") {
      address_ref <- PARAM[which(PARAM[, 1] == 'PARAM0041'), 2] ## x0041
      if (is.na(address_ref)) {
        print("Error!!! PARAM0041 is empty. Please also check PARAM0004!")
        checkpoint_parameter <- FALSE
      } else {
        address_ref <- gsub("\\", "/", address_ref, fixed = TRUE)
        PARAM[which(PARAM[, 1] == 'PARAM0041'), 2] <- address_ref
        if (!dir.exists(address_ref)) {
          print("ERROR!!! Problem with PARAM0041! Directory not found!")
          checkpoint_parameter <- FALSE
        } else {
          ##
          x0042 <- PARAM[which(PARAM[, 1] == 'PARAM0042'), 2]
          if (is.na(x0042)) {
            print("Error!!! PARAM0042 is empty! Please also check PARAM0004!")
            checkpoint_parameter <- FALSE
          } else {
            ref_xlsx_file <- paste0(address_ref, "/", x0042)
            if (file.exists(ref_xlsx_file)) {
              PARAM[which(PARAM[, 1] == 'PARAM0042'), 2] <- ref_xlsx_file
              ref_table <- readxl::read_xlsx(ref_xlsx_file)
              col <- colnames(ref_table)
              x_name <- which(col == 'name')
              x_mz <- which(col == 'm/z')
              x_RT <- which(col == 'RT')
              if (!(length(x_name) > 0 & length(x_mz) > 0 & length(x_RT) > 0)) {
                print("ERROR!!! Problem with PARAM0042! Incorrect column headers in the annotation spreadsheet --> The following columns should be detected in the spreadsheet : 'm/z', 'RT', 'name' - case sensitive")
                checkpoint_parameter <- FALSE
              }
            } else {
              print("ERROR!!! Problem with PARAM0042! The annotation spreadsheet not found!")
              checkpoint_parameter <- FALSE
            }
          }
        }
      }
      ##
      x0043 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0043'), 2])
      if (is.na(x0043)) {
        print("ERROR!!! Problem with PARAM0043! This parameter should be greater than 0 !")
        checkpoint_parameter <- FALSE
      } else {
        if (x0043 <= 0) {
          print("ERROR!!! Problem with PARAM0043! This parameter should be greater than 0 !")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0044 <- as.numeric(PARAM[which(PARAM[, 1] == 'PARAM0044'), 2])
      if (is.na(x0044)) {
        print("ERROR!!! Problem with PARAM0044! This parameter should be greater than 0 !")
        checkpoint_parameter <- FALSE
      } else {
        if (x0044 <= 0) {
          print("ERROR!!! Problem with PARAM0044! This parameter should be greater than 0 !")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0045 <- PARAM[which(PARAM[, 1] == 'PARAM0045'), 2]
      if (is.na(x0045)) {
        print("ERROR!!! Problem with PARAM0045!")
        checkpoint_parameter <- FALSE
      } else {
        if (!(tolower(x0045) == "yes" | tolower(x0045) == "no")) {
          print("ERROR!!! Problem with PARAM0045!")
          checkpoint_parameter <- FALSE
        }
        if (tolower(x0045) == "yes") {
          if (tolower(x0002) != "yes") {
            x0029 <- PARAM[which(PARAM[, 1] == 'PARAM0029'), 2]
            if (is.na(x0029)) {
              print("ERROR!!! Problem with PARAM0029! It should be YES when PARAM0045 is YES")
              checkpoint_parameter <- FALSE
            } else {
              if (tolower(x0029) == "no") {
                print("ERROR!!! Problem with PARAM0029! It should be YES when PARAM0045 is YES")
                checkpoint_parameter <- FALSE
              }
            }
            RT_correction_checkList <- RT_correction_check(checkpoint_parameter, PARAM)
            checkpoint_parameter <- RT_correction_checkList[[1]]
            PARAM <- RT_correction_checkList[[2]]
          }
        }
      }
      ##
      x0046 <- PARAM[which(PARAM[, 1] == 'PARAM0046'), 2]
      if (is.na(x0046)) {
        print("ERROR!!! Problem with PARAM0046!")
        checkpoint_parameter <- FALSE
      } else {
        if (!(tolower(x0046) == "yes" | tolower(x0046) == "no")) {
          print("ERROR!!! Problem with PARAM0046!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0047 <- PARAM[which(PARAM[, 1] == 'PARAM0047'), 2]
      if (is.na(x0047)) {
        print("ERROR!!! Problem with PARAM0047!")
        checkpoint_parameter <- FALSE
      } else {
        if (!(tolower(x0047) == "yes" | tolower(x0047) == "no")) {
          print("ERROR!!! Problem with PARAM0047!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      x0048 <- PARAM[which(PARAM[, 1] == 'PARAM0048'), 2]
      if (is.na(x0048)) {
        print("ERROR!!! Problem with PARAM0048!")
        checkpoint_parameter <- FALSE
      } else {
        if (tolower(x0048) == "yes" | tolower(x0048) == "no") {
          if (tolower(x0048) == "yes") {
            if (tolower(x0003) == "no") {
              checkpoint_parameter <- gapfilling_check(checkpoint_parameter, PARAM)
            }
          }
        } else {
          print("ERROR!!! Problem with PARAM0048!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
  }
  if (checkpoint_parameter == FALSE) {
    print("Please visit   https://ipa.idsl.me    for instructions!")
    PARAM <- NULL
  } else {
    print("The spreadsheet is consistent with the IDSL.IPA workflow!")
  }
  return(PARAM)
}
