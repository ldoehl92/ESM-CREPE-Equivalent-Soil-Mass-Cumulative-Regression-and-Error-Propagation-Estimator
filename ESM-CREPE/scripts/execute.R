### ESM-CREPE package - execute ################################################

# Author:       Leonard Nikolas Konstantin Doehl, PhD, FMCHC
# Affiliation:  Faculty of Environment, School of Geography, 
#               University of Leeds, Leeds LS2 9JT, United Kingdom
# E-mail:       L.N.K.Doehl@leeds.ac.uk

# Package name: ESM-CREPE (Equivalent Soil Mass - Cumulative Regression and 
#               Error Propagation Estimator)
# Version:      v1.1
# Date:         05/01/2026
# Copyright:    (C) 2026 Leonard N. K. Doehl <L.N.K.Doehl@leeds.ac.uk>
#               Open-source license (see "License.docx" for details).

# File name:    execute.R
# File purpose: Import and format data set, select key input parameters to 
#               execute and run the main ESM-CREPE function.
#               Tip: Make a copy of this file and then edit it for use on your
#               own data sets. Edit your copy for additional analyses outside 
#               of the main package.

### Define dataset #############################################################

# define filename of dataset
data_filename <- "./data/case_study_biffi2024.xlsx"
data_sheet <- "Cumbria Hedgerows"

# import dataset
data_import <-
  read_xlsx(path = data_filename, 
            sheet = data_sheet)

# identify the column number of the key inputs in the following sequence:
# 1. Study name:            "study"
# 2. Site name of study:    "site"
# 3. Sample increment:      "increment"
# 4. Core:                  "core"
# 5. Sample thickness (cm): "sample_thickness_cm"
# 6. Fine soil mass (g):    "fine_soil_mass_g"
# 7. SOC content (%):       "soc_content_perc"
# 8. SOM content (%):       "som_content_perc"
# 9. Sample volume (cm3):   "sample_volume_cm3"
data_column_index <- c(1,2,3,4,7,10,6,5,11)

data_reformat <-
  data_renamer(data = data_import, original_label_index = data_column_index)

### Define inputs ##############################################################

# define which is the control site to which soil mineral mass are equated to
input_control_site <- c("control")

# define which content method to calculate SOC stock and soil mineral mass with
# select one value using the [] entry
input_soil_content <- c("soc", "som", "both")[1]

# define which van Bemmelen factor assumption to use
# if guessing the factor value, define that number
# select one value using the [] entry
input_k <- c("conventional", "no intercept", "intercept", "guess", "both")[1]
input_guess_k <- 1.669

# define which regression function(s) to use for ESM correction
# select one, two or all entries; if one or two, use [] and define the entry
input_regression_function <- c("LIN", "MCS", "CEDF")[1]

### Execute ESM correction #####################################################

# run ESM correction
data_esm_results <-
  esm_correction_main(data = data_reformat, 
                      control_site = input_control_site, 
                      soil_content = input_soil_content,
                      k = input_k, 
                      regression_function = input_regression_function,
                      guess_k = input_guess_k)
