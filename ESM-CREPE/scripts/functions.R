### ESM-CREPE package - functions #############################################

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

# File name:    functions.R
# File purpose: All functions to execute the ESM-CREPE correction, calculating
#               the cumulative FD-based and ESM-corrected soil organic carbon
#               stocks as well as associated error propagation analysis and
#               ESM-corrected soil depths.

### Main function to execute run ##############################################

# define main function to ESM correction
esm_correction_main <-
  function(data, control_site, soil_content, k, regression_function, guess_k){
    
    # create empty data frames
    data_final             <- data.frame()
    data_curve             <- data.frame()
    data_calculations      <- data.frame()
    data_error_propagation <- data.frame()
    
    # print start of calculation
    message("Starting ESM correction...")
    message("")
    
    # check that input variables have been defined to accepted inputs
    # terminate run otherwise
    if(!(control_site %in% data$site)){
      
      message("Error: control_site")
      message("   Defined argument incorrectly or missing.")
      message("   Please check the input.")
      message("   Define the control site.")
      message("")
      message("   Terminated ESM-CREPE.")
      
    } else if(!any(soil_content %in% c("soc", "som", "both"))){
      
      message("Error: soil_content")
      message("   Defined argument incorrectly or missing.")
      message("   Please check the input.")
      message("   Define soil_content as any of the following inputs:")
      message("      soc")
      message("      som")
      message("      both")
      message("")
      message("   Terminated ESM-CREPE.")
      
    } else if(!any(k %in% c("conventional", "no intercept", "intercept", "guess", "both"))){
      
      message("Error: k")
      message("   Defined argument incorrectly or missing.")
      message("   Please check the input.")
      message("   Define k as any of the following inputs:")
      message("      conventional")
      message("      no intercept")
      message("      intercept")
      message("      guess")
      message("      both")
      message("")
      message("   Terminated ESM-CREPE.")
      
    } else if(!any(regression_function %in% c("LIN", "MCS", "CEDF"))){
      
      
      message("Error: regression_function")
      message("   Defined argument incorrectly or missing.")
      message("   Please check the input.")
      message("   Define regression_function as any of the following inputs:")
      message("      LIN")
      message("      MCS")
      message("      CEDF")
      message("")
      message("   Terminated ESM-CREPE.")
      
    #} else if(k %in% c("guess") && guess_k %in% NA){
    #  
    #  message("Error: guess_k")
    #  message("   Defined argument incorrectly or missing when.")
    #  message("   k = both")
    #  message("   Please check the input.")
    #  message("   Define guess_k as a number.")
    #  message("")
    #  message("   Terminated ESM-CREPE.")
      
    } else{
      
      message(paste("Import data of study ", unique(data$study), sep = ""))
      message("")
      
      # loop through all inputs
      for(i in 1:length(soil_content)){
        for(j in 1:length(k)){
          for(m in 1:length(regression_function)){
            for(l in 1:length(guess_k)){
              
              # get ESM results
              data_results <-
                esm_calculator(data = data,
                               control_site = control_site, 
                               soil_content = soil_content[i],
                               k = k[j],
                               regression_function = regression_function[m],
                               guess_k = guess_k[l])
              
              # assign and consolidate results
              data_final <- rbind(
                data_final, data.frame(data_results[1]))
              data_curve <- rbind(
                data_curve, data.frame(data_results[2]))
              data_calculations <- rbind(
                data_calculations, data.frame(data_results[3]))
              data_error_propagation <- rbind(
                data_error_propagation, data.frame(data_results[4]))
            }
          }
        }
      }
      
      # remove any duplicated calculations if any
      data_final <-
        data_final %>%
        distinct(`study`, `site`, `increment`, `method`, `soil_content`, 
                 `k`, `regression_function`, `cum_soc_stock_ave_Mgha`, 
                 `cum_soc_stock_err_Mgha`, `cum_soc_stock_ci_Mgha`, 
                 .keep_all = TRUE)
      message("Finished calculations.")
      
      # plot ESM correction curve and uncertainty effects for each treatment site
      message("Plotting visual of ESM correction...")
      visual_analyser(
        data_results = list(data_final, data_curve), 
        control_site = control_site, k = k, soil_content = soil_content)
      message("Plotted ESM-corrected SOC stock. and soil depths.")
      message("")
      message("Plotting uncertainty analysis...")
      uncertainty_analyser(
        data_results = data_error_propagation, 
        control_site = control_site, k = k, soil_content = soil_content)
      message("Plotted uncertainty analysis.")
      message("")
      message("Finished ESM-CREPE run.")
      
      # return results of all systems as list
      return(list(data_final, data_curve, 
                  data_calculations, data_error_propagation)) 
    }
  }

### ESM correction estimator ###################################################

# define ESM calculator
esm_calculator <-
  function(data, control_site, soil_content, k, regression_function, guess_k){
    
    # create empty data frames
    data_final             <- data.frame()
    data_curve             <- data.frame()
    data_calculations      <- data.frame()
    data_error_propagation <- data.frame()
    
    # create TRUE/FALSE binary to run calculations with/without uncertainties 
    # included in calculation
    uncertainty_switch <- 
      c("none", "sample_thickness", "fine_soil_mass", "sample_volume", "soc_content",
        "som_content", "k", "a", "b")
    
    # for loop over all uncertainties included and then excluded
    for(i in 1:length(uncertainty_switch)){
        
      # calculate SOC stock and mineral mass
      data_cum_calculation <-
        stock_and_mass_calculator(
          data, control_site, soil_content, k, guess_k, 
          uncertainty_switch, i)
      
      # create curve of soil mineral mass of fitted regression function
      curve_soil_mineral_mass_g <-
        seq(0, max(data_cum_calculation$cum_soil_mineral_mass_ave2_g) + 50, 10)
      
      # separate data of control and treatment sites
      data_treatment <-
        data_cum_calculation %>%
        dplyr::filter(!(`site` %in% control_site))
      data_control <-
        data_cum_calculation %>%
        dplyr::filter(`site` %in% control_site)
      
      # run ESM correction with all uncertainties included
      if(uncertainty_switch[i] == "none"){
        
        # collect soil mineral masses, soil depth and FD-based SOC stock 
        # of all control and treatment sites
        data_final = 
          rbind(data_final,
                data.frame(
                  study = unique(data$study),
                  site = data_cum_calculation$site,
                  increment = data_cum_calculation$increment,
                  method = "FD",
                  soil_content = soil_content, 
                  k = k, 
                  regression_function = "NA",
                  sample_size = data_cum_calculation$sample_size,
                  cum_soil_mineral_mass_ave_g = data_cum_calculation$cum_soil_mineral_mass_ave2_g,
                  cum_soil_mineral_mass_err_g = data_cum_calculation$cum_soil_mineral_mass_err2_g,
                  cum_soil_mineral_mass_ci_g = data_cum_calculation$cum_soil_mineral_mass_ci2_g,
                  cum_soil_depth_ave_cm = data_cum_calculation$cum_soil_depth_ave2_cm,
                  cum_soil_depth_err_cm = data_cum_calculation$cum_soil_depth_err2_cm,
                  cum_soil_depth_ci_cm = data_cum_calculation$cum_soil_depth_ci2_cm,
                  cum_soc_stock_ave_Mgha = data_cum_calculation$cum_soc_stock_ave2_Mgha,
                  cum_soc_stock_err_Mgha = data_cum_calculation$cum_soc_stock_err2_Mgha,
                  cum_soc_stock_ci_Mgha = data_cum_calculation$cum_soc_stock_ci2_Mgha))
        
        # bind all calculations with no uncertainties switched off
        data_error_propagation <- rbind(
          data_error_propagation, 
          cbind(study = unique(data$study),
                method = "FD",
                soil_content = soil_content, 
                k = k, 
                regression_function = "NA",
                uncertainty_off = "none",
                data.frame(
                  site = data_cum_calculation$site,
                  increment = data_cum_calculation$increment,
                  sample_size = data_cum_calculation$sample_size,
                  cum_soil_mineral_mass_ave_g = data_cum_calculation$cum_soil_mineral_mass_ave2_g,
                  cum_soil_mineral_mass_err_g = data_cum_calculation$cum_soil_mineral_mass_err2_g,
                  cum_soil_mineral_mass_ci_g = data_cum_calculation$cum_soil_mineral_mass_ci2_g,
                  cum_soil_depth_ave_cm = data_cum_calculation$cum_soil_depth_ave2_cm,
                  cum_soil_depth_err_cm = data_cum_calculation$cum_soil_depth_err2_cm,
                  cum_soil_depth_ci_cm = data_cum_calculation$cum_soil_depth_ci2_cm,
                  cum_soc_stock_ave_Mgha = data_cum_calculation$cum_soc_stock_ave2_Mgha,
                  cum_soc_stock_err_Mgha = data_cum_calculation$cum_soc_stock_err2_Mgha,
                  cum_soc_stock_ci_Mgha = data_cum_calculation$cum_soc_stock_ci2_Mgha)))
        
        # for loop to run ESM correction to each treatment site
        for(x in 1:length(unique(data_treatment$site))){
          
          # print which treatment field is used for ESM calculation
          message(paste("Treatment site to fit: ", 
                        unique(data_treatment$site)[x], sep = ""))
          
          # select each individual treatment field
          data_treatment_each <-
            data_treatment %>%
            dplyr::filter(`site` %in% unique(data_treatment$site)[x])
          
          # if regression fitting fails, skip site
          tryCatch({
            
            # fit selected regression function
            data_esm <-
              regression_function_fitting(
                data_treatment = data_treatment_each, 
                data_control = data_control, 
                regression_function = regression_function, 
                curve_soil_mineral_mass = curve_soil_mineral_mass_g)
            
            # get regression function fit and ESM-corrected SOC stocks
            data_regression_curve <- 
              data.frame(data_esm[1]) %>% 
              dplyr::filter(`uncertainty_off` %in% "none") %>%
              dplyr::select(-c("uncertainty_off"))
            data_esm_corrected <- 
              data.frame(data_esm[2]) %>% 
              dplyr::filter(`uncertainty_off` %in% "none") %>%
              dplyr::select(-c("uncertainty_off"))
            
            # relabel and format ESM-corrected stocks
            data_prefinal <-
              data_esm_corrected %>% 
              dplyr::select(
                c("site", "increment", "sample_size",
                  "cum_soil_mineral_mass_ave_g", "cum_soil_mineral_mass_err_g", 
                  "cum_soil_mineral_mass_ci_g", "cum_soil_depth_ave_cm", 
                  "cum_soil_depth_err_cm", "cum_soil_depth_ci_cm",
                  "cum_soc_stock_ave_Mgha", "cum_soc_stock_err_Mgha", 
                  "cum_soc_stock_ci_Mgha"))
            
            # row bind soil mineral mass, soil depth and ESM-corrected SOC 
            # stock to data
            data_final = 
              rbind(data_final,
                    data.frame(
                      study = unique(data$study),
                      site = data_prefinal$site,
                      increment = data_prefinal$increment,
                      method = "ESM",
                      soil_content = soil_content, 
                      k = k, 
                      regression_function = regression_function,
                      sample_size = data_prefinal$sample_size,
                      cum_soil_mineral_mass_ave_g = data_prefinal$cum_soil_mineral_mass_ave_g,
                      cum_soil_mineral_mass_err_g = data_prefinal$cum_soil_mineral_mass_err_g,
                      cum_soil_mineral_mass_ci_g = data_prefinal$cum_soil_mineral_mass_ci_g,
                      cum_soil_depth_ave_cm = data_prefinal$cum_soil_depth_ave_cm,
                      cum_soil_depth_err_cm = data_prefinal$cum_soil_depth_err_cm,
                      cum_soil_depth_ci_cm = data_prefinal$cum_soil_depth_ci_cm,
                      cum_soc_stock_ave_Mgha = data_prefinal$cum_soc_stock_ave_Mgha,
                      cum_soc_stock_err_Mgha = data_prefinal$cum_soc_stock_err_Mgha,
                      cum_soc_stock_ci_Mgha = data_prefinal$cum_soc_stock_ci_Mgha))
            
            # row bind regression function curves to data
            data_curve <-
              rbind(data_curve,
                    data.frame(
                      study = unique(data$study),
                      site = data_regression_curve$site,
                      method = "ESM",
                      soil_content = soil_content, 
                      k = k, 
                      regression_function = regression_function,
                      cum_soil_mineral_mass_ave_g = data_regression_curve$cum_soil_mineral_mass_ave_g,
                      cum_soil_depth_ave_cm = data_regression_curve$cum_soil_depth_ave_cm,
                      cum_soil_depth_err_cm = data_regression_curve$cum_soil_depth_err_cm,
                      cum_soil_depth_ci_cm = data_regression_curve$cum_soil_depth_ci_cm,
                      cum_soc_stock_ave_Mgha = data_regression_curve$cum_soc_stock_ave_Mgha,
                      cum_soc_stock_err_Mgha = data_regression_curve$cum_soc_stock_err_Mgha,
                      cum_soc_stock_ci_Mgha = data_regression_curve$cum_soc_stock_ci_Mgha))
            
            # row bind all raw calculations for reference and reproducibility
            data_calculations <- rbind(
              data_calculations,
              data_cum_calculation
            )
            
            # bind all calculations with no uncertainties switched off
            data_error_propagation <- rbind(
              data_error_propagation, 
              cbind(study = unique(data$study),
                    method = "ESM",
                    soil_content = soil_content, 
                    k = k, 
                    regression_function = regression_function,
                    uncertainty_off = "none",
                    data_prefinal))
          }, 
          error = function(cond) {
            message("")
            NA},
          warning = function(cond) {
            message("")
            NULL},
          finally = {
            message("")
          })
        }
        
        # run ESM correction with one uncertainty excluded  
      } else{
        
        # for loop to run ESM correction to each treatment site
        for(x in 1:length(unique(data_treatment$site))){
          
          # print which treatment field is used for ESM calculation
          message(paste("Treatment site to fit: ", 
                        unique(data_treatment$site)[x], sep = ""))
          
          # select each individual treatment field
          data_treatment_each <-
            data_treatment %>%
            dplyr::filter(`site` %in% unique(data_treatment$site)[x])
          
          # if regression fitting fails, skip site
          tryCatch({
            
            # fit selected regression function
            data_esm <-
              regression_function_fitting(
                data_treatment = data_treatment_each, 
                data_control = data_control, 
                regression_function = regression_function, 
                curve_soil_mineral_mass = curve_soil_mineral_mass_g)
            
            # identify which direct measurement uncertainty is excluded
            if(uncertainty_switch[i] == "sample_thickness" ||
               uncertainty_switch[i] == "fine_soil_mass" || 
               uncertainty_switch[i] == "sample_volume" ||
               uncertainty_switch[i] == "soc_content" ||
               uncertainty_switch[i] == "som_content" ||
               uncertainty_switch[i] == "k"){
              
              # bind all calculations with one uncertainty excluded
              data_error_propagation <- rbind(
                data_error_propagation, 
                cbind(study = unique(data$study),
                      method = "FD",
                      soil_content = soil_content, 
                      k = k, 
                      regression_function = "NA",
                      uncertainty_off = uncertainty_switch[i],
                      data.frame(
                        site = data_cum_calculation$site,
                        increment = data_cum_calculation$increment,
                        sample_size = data_cum_calculation$sample_size,
                        cum_soil_mineral_mass_ave_g = data_cum_calculation$cum_soil_mineral_mass_ave2_g,
                        cum_soil_mineral_mass_err_g = data_cum_calculation$cum_soil_mineral_mass_err2_g,
                        cum_soil_mineral_mass_ci_g = data_cum_calculation$cum_soil_mineral_mass_ci2_g,
                        cum_soil_depth_ave_cm = data_cum_calculation$cum_soil_depth_ave2_cm,
                        cum_soil_depth_err_cm = data_cum_calculation$cum_soil_depth_err2_cm,
                        cum_soil_depth_ci_cm = data_cum_calculation$cum_soil_depth_ci2_cm,
                        cum_soc_stock_ave_Mgha = data_cum_calculation$cum_soc_stock_ave2_Mgha,
                        cum_soc_stock_err_Mgha = data_cum_calculation$cum_soc_stock_err2_Mgha,
                        cum_soc_stock_ci_Mgha = data_cum_calculation$cum_soc_stock_ci2_Mgha)))
              
              # get regression function fit and ESM-corrected SOC stocks
              data_esm_corrected <- 
                data.frame(data_esm[2]) %>% 
                dplyr::filter(`uncertainty_off` %in% "none") %>%
                dplyr::select(-c("uncertainty_off"))
              
              # relabel and format ESM-corrected stocks
              data_prefinal <-
                data_esm_corrected %>% 
                dplyr::select(
                  c("site", "increment", "sample_size",
                    "cum_soil_mineral_mass_ave_g", "cum_soil_mineral_mass_err_g", 
                    "cum_soil_mineral_mass_ci_g", "cum_soil_depth_ave_cm", 
                    "cum_soil_depth_err_cm", "cum_soil_depth_ci_cm",
                    "cum_soc_stock_ave_Mgha", "cum_soc_stock_err_Mgha", 
                    "cum_soc_stock_ci_Mgha"))
              
              # bind all calculations with one uncertainty excluded
              data_error_propagation <- rbind(
                data_error_propagation, 
                cbind(study = unique(data$study),
                      method = "ESM",
                      soil_content = soil_content, 
                      k = k, 
                      regression_function = regression_function,
                      uncertainty_off = uncertainty_switch[i],
                      data_prefinal))
              
              # identify whether CEDF coefficient 'a' uncertainty is excluded
            } else if(uncertainty_switch[i] == "a"){
              
              # bind all calculations with one uncertainty excluded
              data_error_propagation <- rbind(
                data_error_propagation, 
                cbind(study = unique(data$study),
                      method = "FD",
                      soil_content = soil_content, 
                      k = k, 
                      regression_function = "NA",
                      uncertainty_off = uncertainty_switch[i],
                      data.frame(
                        site = data_cum_calculation$site,
                        increment = data_cum_calculation$increment,
                        sample_size = data_cum_calculation$sample_size,
                        cum_soil_mineral_mass_ave_g = data_cum_calculation$cum_soil_mineral_mass_ave2_g,
                        cum_soil_mineral_mass_err_g = data_cum_calculation$cum_soil_mineral_mass_err2_g,
                        cum_soil_mineral_mass_ci_g = data_cum_calculation$cum_soil_mineral_mass_ci2_g,
                        cum_soil_depth_ave_cm = data_cum_calculation$cum_soil_depth_ave2_cm,
                        cum_soil_depth_err_cm = data_cum_calculation$cum_soil_depth_err2_cm,
                        cum_soil_depth_ci_cm = data_cum_calculation$cum_soil_depth_ci2_cm,
                        cum_soc_stock_ave_Mgha = data_cum_calculation$cum_soc_stock_ave2_Mgha,
                        cum_soc_stock_err_Mgha = data_cum_calculation$cum_soc_stock_err2_Mgha,
                        cum_soc_stock_ci_Mgha = data_cum_calculation$cum_soc_stock_ci2_Mgha)))
              
              # get regression function fit and ESM-corrected SOC stocks
              data_esm_corrected <- 
                data.frame(data_esm[2]) %>% 
                dplyr::filter(`uncertainty_off` %in% "a") %>%
                dplyr::select(-c("uncertainty_off"))
              
              # relabel and format ESM-corrected stocks
              data_prefinal <-
                data_esm_corrected %>% 
                dplyr::select(
                  c("site", "increment", "sample_size",
                    "cum_soil_mineral_mass_ave_g", "cum_soil_mineral_mass_err_g", 
                    "cum_soil_mineral_mass_ci_g", "cum_soil_depth_ave_cm", 
                    "cum_soil_depth_err_cm", "cum_soil_depth_ci_cm",
                    "cum_soc_stock_ave_Mgha", "cum_soc_stock_err_Mgha", 
                    "cum_soc_stock_ci_Mgha"))
              
              # bind all calculations with one uncertainty excluded
              data_error_propagation <- rbind(
                data_error_propagation, 
                cbind(study = unique(data$study),
                      method = "ESM",
                      soil_content = soil_content, 
                      k = k, 
                      regression_function = regression_function,
                      uncertainty_off = uncertainty_switch[i],
                      data_prefinal))
              
              # identify whether CEDF coefficient 'b' uncertainty is excluded
            } else if(uncertainty_switch[i] == "b"){
              
              # bind all calculations with one uncertainty excluded
              data_error_propagation <- rbind(
                data_error_propagation, 
                cbind(study = unique(data$study),
                      method = "FD",
                      soil_content = soil_content, 
                      k = k, 
                      regression_function = "NA",
                      uncertainty_off = uncertainty_switch[i],
                      data.frame(
                        site = data_cum_calculation$site,
                        increment = data_cum_calculation$increment,
                        sample_size = data_cum_calculation$sample_size,
                        cum_soil_mineral_mass_ave_g = data_cum_calculation$cum_soil_mineral_mass_ave2_g,
                        cum_soil_mineral_mass_err_g = data_cum_calculation$cum_soil_mineral_mass_err2_g,
                        cum_soil_mineral_mass_ci_g = data_cum_calculation$cum_soil_mineral_mass_ci2_g,
                        cum_soil_depth_ave_cm = data_cum_calculation$cum_soil_depth_ave2_cm,
                        cum_soil_depth_err_cm = data_cum_calculation$cum_soil_depth_err2_cm,
                        cum_soil_depth_ci_cm = data_cum_calculation$cum_soil_depth_ci2_cm,
                        cum_soc_stock_ave_Mgha = data_cum_calculation$cum_soc_stock_ave2_Mgha,
                        cum_soc_stock_err_Mgha = data_cum_calculation$cum_soc_stock_err2_Mgha,
                        cum_soc_stock_ci_Mgha = data_cum_calculation$cum_soc_stock_ci2_Mgha)))
              
              # get regression function fit and ESM-corrected SOC stocks
              data_esm_corrected <- 
                data.frame(data_esm[2]) %>% 
                dplyr::filter(`uncertainty_off` %in% "b") %>%
                dplyr::select(-c("uncertainty_off"))
              
              # relabel and format ESM-corrected stocks
              data_prefinal <-
                data_esm_corrected %>% 
                dplyr::select(
                  c("site", "increment", "sample_size",
                    "cum_soil_mineral_mass_ave_g", "cum_soil_mineral_mass_err_g", 
                    "cum_soil_mineral_mass_ci_g", "cum_soil_depth_ave_cm", 
                    "cum_soil_depth_err_cm", "cum_soil_depth_ci_cm",
                    "cum_soc_stock_ave_Mgha", "cum_soc_stock_err_Mgha", 
                    "cum_soc_stock_ci_Mgha"))
              
              # bind all calculations with one uncertainty excluded
              data_error_propagation <- rbind(
                data_error_propagation, 
                cbind(study = unique(data$study),
                      method = "ESM",
                      soil_content = soil_content, 
                      k = k, 
                      regression_function = regression_function,
                      uncertainty_off = uncertainty_switch[i],
                      data_prefinal))
              
            } else{}
          }, 
          error = function(cond) {
            message("Error in calculation - skipped.")
            message(unique(data_treatment$site)[x])
            message("Here is the original error in method:")
            message(conditionMessage(cond))
            NA},
          warning = function(cond) {
            message("Warning in calculation - skipped.")
            message(unique(data_treatment$site)[x])
            message("Here is the original warning in method:")
            message(conditionMessage(cond))
            NULL},
          finally = {
            message("")
          })
        }
      }
    }
    
    # return data frames
    return(list(data_final, data_curve, data_calculations, data_error_propagation))
  }

### Cumulative FD-based SOC stock and soil mineral mass calculator #############

# define function to calculate SOC stock using the FD method (Poeplau, 2017)
# and soil mineral mass
stock_and_mass_calculator <-
  function(data, control_site, soil_content, k, guess_k, uncertainty_switch, i){
    
    # run van Bemmelen factor calculation
    data_van_bemmelen <- van_bemmelen_selector(data, k, guess_k)
    
    # merge data with van Bemmelen factor
    data_merge1 <-
      merge(data, data_van_bemmelen, by = c("study", "site"))
    
    if(soil_content %in% "both"){
      
      # calculate measurement uncertainties and SOC stocks and soil mineral masses
      # per sampling increment using 'Both':
      # - SOC content for SOC stock
      # - SOM content for soil mineral mass
      data_stock_precalculate <-
        data_merge1 %>%
        dplyr::group_by(`site`, `increment`) %>%
        dplyr::summarise(
          sample_size = n(),
          sample_thickness_ave2_cm = 
            mean(`sample_thickness_cm`, na.rm = TRUE),
          sample_thickness_dev2_cm = 
            if(uncertainty_switch[i] == "sample_thickness"){0
              } else{sd(`sample_thickness_cm`, na.rm = TRUE)},
          sample_thickness_err2_cm = 
            `sample_thickness_dev2_cm` / sqrt(`sample_size`),
          sample_thickness_ci2_cm = 
            qt(0.975, df = `sample_size` - 1) * `sample_thickness_err2_cm`,
          fine_soil_mass_ave2_g = 
            mean(`fine_soil_mass_g`, na.rm = TRUE),
          fine_soil_mass_dev2_g = 
            if(uncertainty_switch[i] == "fine_soil_mass"){0
              } else{sd(`fine_soil_mass_g`, na.rm = TRUE)},
          fine_soil_mass_err2_g = 
            `fine_soil_mass_dev2_g` / sqrt(`sample_size`),
          fine_soil_mass_ci2_g = 
            qt(0.975, df = `sample_size` - 1) * `fine_soil_mass_err2_g`,
          sample_volume_ave2_cm3 = 
            mean(`sample_volume_cm3`, na.rm = TRUE),
          sample_volume_dev2_cm3 = 
            if(uncertainty_switch[i] == "sample_volume"){0
              } else{sd(`sample_volume_cm3`, na.rm = TRUE)},
          sample_volume_err2_cm3 = 
            `sample_volume_dev2_cm3` / sqrt(`sample_size`),
          sample_volume_ci2_cm3 = 
            qt(0.975, df = `sample_size` - 1) * `sample_volume_err2_cm3`,
          soc_content_ave2_perc = 
            mean(`soc_content_perc`, na.rm = TRUE),
          soc_content_dev2_perc = 
            if(uncertainty_switch[i] == "soc_content"){0
              } else{sd(`soc_content_perc`, na.rm = TRUE)},
          soc_content_err2_perc = 
            `soc_content_dev2_perc` / sqrt(`sample_size`),
          soc_content_ci2_perc = 
            qt(0.975, df = `sample_size` - 1) * `soc_content_err2_perc`,
          som_content_ave2_perc = 
            mean(`som_content_perc`, na.rm = TRUE),
          som_content_dev2_perc = 
            if(uncertainty_switch[i] == "som_content"){0
              } else{sd(`som_content_perc`, na.rm = TRUE)},
          som_content_err2_perc = 
            `som_content_dev2_perc` / sqrt(`sample_size`),
          som_content_ci2_perc = 
            qt(0.975, df = `sample_size` - 1) * `som_content_err2_perc`,
          k_ave2 = 
            mean(`k_ave`),
          k_err2 = 
            if(uncertainty_switch[i] == "k"){0
              } else{mean(`k_err`)},
          k_ci2 = 
            if(uncertainty_switch[i] == "k"){0
              } else{mean(`k_ci`)},
          cov_thick_soc = 
            if(uncertainty_switch[i] == "sample_thickness" || 
               uncertainty_switch[i] == "soc_content"){0
              } else{2 * cov(`sample_thickness_cm`, `soc_content_perc`, 
                             use = "complete.obs") /
                  (`sample_thickness_ave2_cm` * `soc_content_ave2_perc`)},
          cov_mass_soc = 
            if(uncertainty_switch[i] == "fine_soil_mass" || 
               uncertainty_switch[i] == "soc_content"){0
              } else{2 * cov(`fine_soil_mass_g`, `soc_content_perc`, 
                             use = "complete.obs") / 
                  (`fine_soil_mass_ave2_g` * `soc_content_ave2_perc`)},
          cov_volume_soc = 
            if(uncertainty_switch[i] == "sample_volume" || 
               uncertainty_switch[i] == "soc_content"){0
              } else{2 * cov(`sample_volume_cm3`, `soc_content_perc`, 
                             use = "complete.obs") /
                  (`sample_volume_ave2_cm3` * `soc_content_ave2_perc`)},
          cov_soc_k = 
            if(uncertainty_switch[i] == "soc_content" || 
               uncertainty_switch[i] == "k"){0
              } else{2 * cov(`soc_content_perc`, `k_ave`,
                             use = "complete.obs") / 
                  (`soc_content_ave2_perc` * `k_ave2`)},
          cov_thick_mass = 
            if(uncertainty_switch[i] == "sample_thickness" || 
               uncertainty_switch[i] == "fine_soil_mass"){0
              } else{2 * cov(`sample_thickness_cm`, `fine_soil_mass_g`, 
                             use = "complete.obs") /
                  (`sample_thickness_ave2_cm` * `fine_soil_mass_ave2_g`)},
          cov_thick_volume = 
            if(uncertainty_switch[i] == "sample_thickness" || 
               uncertainty_switch[i] == "sample_volume"){0
              } else{2 * cov(`sample_thickness_cm`, `sample_volume_cm3`, 
                             use = "complete.obs") / 
                  (`sample_thickness_ave2_cm` * `sample_volume_ave2_cm3`)},
          cov_thick_som = 
            if(uncertainty_switch[i] == "sample_thickness" || 
               uncertainty_switch[i] == "som_content"){0
              } else{2 * cov(`sample_thickness_cm`, `som_content_perc`, 
                             use = "complete.obs") /
                  (`sample_thickness_ave2_cm` * `som_content_ave2_perc`)},
          cov_mass_volume = 
            if(uncertainty_switch[i] == "fine_soil_mass" || 
               uncertainty_switch[i] == "sample_volume"){0
              } else{2 * cov(`fine_soil_mass_g`, `sample_volume_cm3`, 
                             use = "complete.obs") / 
                  (`fine_soil_mass_ave2_g` * `sample_volume_ave2_cm3`)},
          cov_mass_som = 
            if(uncertainty_switch[i] == "fine_soil_mass" || 
               uncertainty_switch[i] == "som_content"){0
              } else{2 * cov(`fine_soil_mass_g`, `som_content_perc`, 
                             use = "complete.obs") / 
                  (`fine_soil_mass_ave2_g` * `som_content_ave2_perc`)},
          cov_volume_som = 
            if(uncertainty_switch[i] == "sample_volume" || 
               uncertainty_switch[i] == "som_content"){0
              } else{2 * cov(`sample_volume_cm3`, `som_content_perc`, 
                             use = "complete.obs") /
                  (`sample_volume_ave2_cm3` * `som_content_ave2_perc`)}
        ) %>%
        dplyr::mutate(
          coe_thick_soc = `cov_thick_soc` / `sample_size`,
          coe_mass_soc = `cov_mass_soc` / `sample_size`,
          coe_volume_soc = `cov_volume_soc` / `sample_size`,
          coe_soc_k = `cov_soc_k` / `sample_size`,
          coe_thick_mass = `cov_thick_mass` / `sample_size`,
          coe_thick_volume = `cov_thick_volume` / `sample_size`,
          coe_thick_som = `cov_thick_som` / `sample_size`,
          coe_mass_volume = `cov_mass_volume` / `sample_size`,
          coe_mass_som = `cov_mass_som` / `sample_size`,
          coe_volume_som = `cov_volume_som` / `sample_size`
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          soc_stock_ave2_Mgha = 
            `fine_soil_mass_ave2_g` * `soc_content_ave2_perc` * 
            `sample_thickness_ave2_cm` / `sample_volume_ave2_cm3`,
          soc_stock_err2_Mgha = `soc_stock_ave2_Mgha` *
            sqrt(((`fine_soil_mass_err2_g` / `fine_soil_mass_ave2_g`) ^ 2) + 
                   ((`soc_content_err2_perc` / `soc_content_ave2_perc`) ^ 2) +
                   ((`sample_thickness_err2_cm` / `sample_thickness_ave2_cm`) ^ 2) +
                   ((`sample_volume_err2_cm3` / `sample_volume_ave2_cm3`) ^ 2) +
                   `coe_thick_mass` - `coe_thick_volume` + `coe_thick_soc` - 
                   `coe_mass_volume` + `coe_mass_soc` - `coe_volume_soc`),
          soc_stock_ci2_Mgha = `soc_stock_ave2_Mgha` *
            sqrt(((`fine_soil_mass_ci2_g` / `fine_soil_mass_ave2_g`) ^ 2) + 
                   ((`soc_content_ci2_perc` / `soc_content_ave2_perc`) ^ 2) +
                   ((`sample_thickness_ci2_cm` / `sample_thickness_ave2_cm`) ^ 2) +
                   ((`sample_volume_ci2_cm3` / `sample_volume_ave2_cm3`) ^ 2) +
                   qt(0.975, df = `sample_size` - 1) *
                   (`coe_thick_mass` - `coe_thick_volume` + `coe_thick_soc` - 
                      `coe_mass_volume` + `coe_mass_soc` - `coe_volume_soc`)),
          soil_mineral_mass_ave2_g =
            (1 - 0.01 * `som_content_ave2_perc`) * `fine_soil_mass_ave2_g`,
          soil_mineral_mass_err2_g =
            sqrt(((`soil_mineral_mass_ave2_g` * `fine_soil_mass_err2_g` / 
                     `fine_soil_mass_ave2_g`) ^ 2) +
                   ((`soil_mineral_mass_ave2_g` * 0.01 * `som_content_err2_perc` / 
                       (1 - 0.01 * `som_content_ave2_perc`)) ^ 2) -
                   (`fine_soil_mass_ave2_g` * `som_content_ave2_perc` * 
                      `soil_mineral_mass_ave2_g` * `coe_mass_som`)),
          soil_mineral_mass_ci2_g =
            sqrt(((`soil_mineral_mass_ave2_g` * `fine_soil_mass_ci2_g` / 
                     `fine_soil_mass_ave2_g`) ^ 2) +
                   ((`soil_mineral_mass_ave2_g` * 0.01 * `som_content_ci2_perc` / 
                       (1 - 0.01 * `som_content_ave2_perc`)) ^ 2) -
                   (qt(0.975, df = `sample_size` - 1) * 
                      `fine_soil_mass_ave2_g` * `som_content_ave2_perc` * 
                      `soil_mineral_mass_ave2_g` * `coe_mass_som`))
        )
      
      # calculate per sample location SOC stock and soil mineral mass
      data_stock_scov <- 
        data_merge1 %>%
        dplyr::mutate(soc_stock_Mgha = 
                        `fine_soil_mass_g` * `soc_content_perc` * 
                        `sample_thickness_cm` / `sample_volume_cm3`,
                      soil_mineral_mass_g = (1 - 0.01 * `som_content_perc`) *
                        `fine_soil_mass_g`)
      
    } else if(soil_content %in% "soc"){
      
      # calculate measurement uncertainties and SOC stocks and soil mineral masses
      # per sampling increment using 'SOC' (SOC content)
      data_stock_precalculate <-
        data_merge1 %>%
        dplyr::group_by(`site`, `increment`) %>%
        dplyr::summarise(
          sample_size = n(),
          sample_thickness_ave2_cm = 
            mean(`sample_thickness_cm`, na.rm = TRUE),
          sample_thickness_dev2_cm = 
            if(uncertainty_switch[i] == "sample_thickness"){0
              } else{sd(`sample_thickness_cm`, na.rm = TRUE)},
          sample_thickness_err2_cm = 
            `sample_thickness_dev2_cm` / sqrt(`sample_size`),
          sample_thickness_ci2_cm = 
            qt(0.975, df = `sample_size` - 1) * `sample_thickness_err2_cm`,
          fine_soil_mass_ave2_g = 
            mean(`fine_soil_mass_g`, na.rm = TRUE),
          fine_soil_mass_dev2_g = 
            if(uncertainty_switch[i] == "fine_soil_mass"){0
              } else{sd(`fine_soil_mass_g`, na.rm = TRUE)},
          fine_soil_mass_err2_g = 
            `fine_soil_mass_dev2_g` / sqrt(`sample_size`),
          fine_soil_mass_ci2_g = 
            qt(0.975, df = `sample_size` - 1) * `fine_soil_mass_err2_g`,
          sample_volume_ave2_cm3 = 
            mean(`sample_volume_cm3`, na.rm = TRUE),
          sample_volume_dev2_cm3 = 
            if(uncertainty_switch[i] == "sample_volume"){0
              } else{sd(`sample_volume_cm3`, na.rm = TRUE)},
          sample_volume_err2_cm3 = 
            `sample_volume_dev2_cm3` / sqrt(`sample_size`),
          sample_volume_ci2_cm3 = 
            qt(0.975, df = `sample_size` - 1) * `sample_volume_err2_cm3`,
          soc_content_ave2_perc = 
            mean(`soc_content_perc`, na.rm = TRUE),
          soc_content_dev2_perc = 
            if(uncertainty_switch[i] == "soc_content"){0
              } else{sd(`soc_content_perc`, na.rm = TRUE)},
          soc_content_err2_perc = 
            `soc_content_dev2_perc` / sqrt(`sample_size`),
          soc_content_ci2_perc = 
            qt(0.975, df = `sample_size` - 1) * `soc_content_err2_perc`,
          k_ave2 = 
            mean(`k_ave`),
          k_err2 = 
            if(uncertainty_switch[i] == "k"){0
              } else{mean(`k_err`)},
          k_ci2 = 
            if(uncertainty_switch[i] == "k"){0
              } else{mean(`k_ci`)},
          cov_depth_mass = 
            if(uncertainty_switch[i] == "sample_thickness" || 
               uncertainty_switch[i] == "fine_soil_mass"){0
              } else{2 * cov(`sample_thickness_cm`, `fine_soil_mass_g`, 
                             use = "complete.obs") /
                  (`sample_thickness_ave2_cm` * `fine_soil_mass_ave2_g`)},
          cov_depth_volume = 
            if(uncertainty_switch[i] == "sample_thickness" || 
               uncertainty_switch[i] == "sample_volume"){0
              } else{2 * cov(`sample_thickness_cm`, `sample_volume_cm3`, 
                             use = "complete.obs") / 
                  (`sample_thickness_ave2_cm` * `sample_volume_ave2_cm3`)},
          cov_depth_soc = 
            if(uncertainty_switch[i] == "sample_thickness" || 
               uncertainty_switch[i] == "soc_content"){0
              } else{2 * cov(`sample_thickness_cm`, `soc_content_perc`, 
                             use = "complete.obs") /
                  (`sample_thickness_ave2_cm` * `soc_content_ave2_perc`)},
          cov_mass_volume = 
            if(uncertainty_switch[i] == "fine_soil_mass" || 
               uncertainty_switch[i] == "sample_volume"){0
              } else{2 * cov(`fine_soil_mass_g`, `sample_volume_cm3`, 
                             use = "complete.obs") / 
                  (`fine_soil_mass_ave2_g` * `sample_volume_ave2_cm3`)},
          cov_mass_soc = 
            if(uncertainty_switch[i] == "fine_soil_mass" || 
               uncertainty_switch[i] == "soc_content"){0
              } else{2 * cov(`fine_soil_mass_g`, `soc_content_perc`, 
                             use = "complete.obs") / 
                  (`fine_soil_mass_ave2_g` * `soc_content_ave2_perc`)},
          cov_volume_soc = 
            if(uncertainty_switch[i] == "sample_volume" || 
               uncertainty_switch[i] == "soc_content"){0
              } else{2 * cov(`sample_volume_cm3`, `soc_content_perc`, 
                             use = "complete.obs") /
                  (`sample_volume_ave2_cm3` * `soc_content_ave2_perc`)},
          cov_mass_k = 
            if(uncertainty_switch[i] == "fine_soil_mass" || 
               uncertainty_switch[i] == "k"){0
              } else{2 * cov(`fine_soil_mass_g`, `k_ave`, 
                             use = "complete.obs") / 
                  (`fine_soil_mass_ave2_g` * `k_ave2`)},
          cov_soc_k = 
            if(uncertainty_switch[i] == "soc_content" || 
               uncertainty_switch[i] == "k"){0
            } else{2 * cov(`soc_content_perc`, `k_ave`, 
                           use = "complete.obs") / 
                  (`soc_content_ave2_perc` * `k_ave2`)}
        ) %>%
        mutate(
          coe_depth_mass = `cov_depth_mass` / `sample_size`,
          coe_depth_volume = `cov_depth_volume` / `sample_size`,
          coe_depth_soc = `cov_depth_soc` / `sample_size`,
          coe_mass_volume = `cov_mass_volume` / `sample_size`,
          coe_mass_soc = `cov_mass_soc` / `sample_size`,
          coe_volume_soc = `cov_volume_soc` / `sample_size`,
          coe_mass_k = `cov_mass_k` / `sample_size`,
          coe_soc_k = `cov_soc_k` / `sample_size`
        ) %>%
        ungroup() %>%
        mutate(
          soc_stock_ave2_Mgha = 
            `fine_soil_mass_ave2_g` * `soc_content_ave2_perc` * 
            `sample_thickness_ave2_cm` / `sample_volume_ave2_cm3`,
          soc_stock_err2_Mgha = `soc_stock_ave2_Mgha` *
            sqrt(((`fine_soil_mass_err2_g` / `fine_soil_mass_ave2_g`) ^ 2) + 
                   ((`soc_content_err2_perc` / `soc_content_ave2_perc`) ^ 2) +
                   ((`sample_thickness_err2_cm` / `sample_thickness_ave2_cm`) ^ 2) +
                   ((`sample_volume_err2_cm3` / `sample_volume_ave2_cm3`) ^ 2) +
                   `coe_depth_mass` - `coe_depth_volume` + `coe_depth_soc` - 
                   `coe_mass_volume` + `coe_mass_soc` - `coe_volume_soc`),
          soc_stock_ci2_Mgha = `soc_stock_ave2_Mgha` *
            sqrt(((`fine_soil_mass_ci2_g` / `fine_soil_mass_ave2_g`) ^ 2) + 
                   ((`soc_content_ci2_perc` / `soc_content_ave2_perc`) ^ 2) +
                   ((`sample_thickness_ci2_cm` / `sample_thickness_ave2_cm`) ^ 2) +
                   ((`sample_volume_ci2_cm3` / `sample_volume_ave2_cm3`) ^ 2) +
                   qt(0.975, df = `sample_size` - 1) * 
                   (`coe_depth_mass` - `coe_depth_volume` + `coe_depth_soc` - 
                      `coe_mass_volume` + `coe_mass_soc` - `coe_volume_soc`)),
          soil_mineral_mass_ave2_g =
            (1 - 0.01 * `soc_content_ave2_perc` * `k_ave2`) * 
            `fine_soil_mass_ave2_g`,
          soil_mineral_mass_err2_g =
            sqrt(((`soil_mineral_mass_ave2_g` * `fine_soil_mass_err2_g` / 
                     `fine_soil_mass_ave2_g`) ^ 2) +
                   ((`fine_soil_mass_ave2_g` - `soil_mineral_mass_ave2_g`) ^ 2) * 
                   (((`k_err2` / `k_ave2`) ^ 2) + 
                      ((`soc_content_err2_perc` / `soc_content_ave2_perc`) ^ 2)) +
                   (((`fine_soil_mass_ave2_g` - `soil_mineral_mass_ave2_g`) ^ 2) * 
                      `coe_soc_k`) +
                   (`soil_mineral_mass_ave2_g` * 
                      (`fine_soil_mass_ave2_g` - `soil_mineral_mass_ave2_g`) * 
                      (`coe_mass_k` + `coe_mass_soc`))),
          soil_mineral_mass_ci2_g =
            sqrt(((`soil_mineral_mass_ave2_g` * `fine_soil_mass_ci2_g` / 
                     `fine_soil_mass_ave2_g`) ^ 2) +
                   ((`fine_soil_mass_ave2_g` - `soil_mineral_mass_ave2_g`) ^ 2) * 
                   (((`k_ci2` / `k_ave2`) ^ 2) + 
                      ((`soc_content_ci2_perc` / `soc_content_ave2_perc`) ^ 2)) +
                   qt(0.975, df = `sample_size` - 1) *
                   ((((`fine_soil_mass_ave2_g` - `soil_mineral_mass_ave2_g`) ^ 2) *
                       `coe_soc_k`) +
                      (`soil_mineral_mass_ave2_g` * 
                         (`fine_soil_mass_ave2_g` - `soil_mineral_mass_ave2_g`) * 
                         (`coe_mass_k` + `coe_mass_soc`)))),
        )
      
      # calculate per composite stock and mineral mass
      data_stock_scov <- 
        data_merge1 %>%
        dplyr::mutate(soc_stock_Mgha = 
                        `fine_soil_mass_g` * `soc_content_perc` * 
                        `sample_thickness_cm` / `sample_volume_cm3`,
                      soil_mineral_mass_g = 
                        (1 - 0.01 * `k_ave` * `soc_content_perc`) * 
                        `fine_soil_mass_g`)
      
    } else if(soil_content %in% "som"){
      
      # calculate measurement uncertainties and SOC stocks and soil mineral masses
      # per sampling increment using 'SOM' (SOM content)
      data_stock_precalculate <-
        data_merge1 %>%
        dplyr::group_by(`site`, `increment`) %>%
        dplyr::summarise(
          sample_size = n(),
          sample_thickness_ave2_cm = 
            mean(`sample_thickness_cm`, na.rm = TRUE),
          sample_thickness_dev2_cm = 
            if(uncertainty_switch[i] == "sample_thickness"){0
              } else{sd(`sample_thickness_cm`, na.rm = TRUE)},
          sample_thickness_err2_cm = 
            `sample_thickness_dev2_cm` / sqrt(`sample_size`),
          sample_thickness_ci2_cm = 
            qt(0.975, df = `sample_size` - 1) * `sample_thickness_err2_cm`,
          fine_soil_mass_ave2_g = 
            mean(`fine_soil_mass_g`, na.rm = TRUE),
          fine_soil_mass_dev2_g = 
            if(uncertainty_switch[i] == "fine_soil_mass"){0
              } else{sd(`fine_soil_mass_g`, na.rm = TRUE)},
          fine_soil_mass_err2_g = 
            `fine_soil_mass_dev2_g` / sqrt(`sample_size`),
          fine_soil_mass_ci2_g = 
            qt(0.975, df = `sample_size` - 1) * `fine_soil_mass_err2_g`,
          sample_volume_ave2_cm3 = 
            mean(`sample_volume_cm3`, na.rm = TRUE),
          sample_volume_dev2_cm3 = 
            if(uncertainty_switch[i] == "sample_volume"){0
              } else{sd(`sample_volume_cm3`, na.rm = TRUE)},
          sample_volume_err2_cm3 = 
            `sample_volume_dev2_cm3` / sqrt(`sample_size`),
          sample_volume_ci2_cm3 = 
            qt(0.975, df = `sample_size` - 1) * `sample_volume_err2_cm3`,
          som_content_ave2_perc = 
            mean(`som_content_perc`, na.rm = TRUE),
          som_content_dev2_perc = 
            if(uncertainty_switch[i] == "som_content"){0
              } else{sd(`som_content_perc`, na.rm = TRUE)},
          som_content_err2_perc = 
            `som_content_dev2_perc` / sqrt(`sample_size`),
          som_content_ci2_perc = 
            qt(0.975, df = `sample_size` - 1) * `som_content_err2_perc`,
          k_ave2 = 
            mean(`k_ave`),
          k_err2 = 
            if(uncertainty_switch[i] == "k"){0
              } else{mean(`k_err`)},
          k_ci2 = 
            if(uncertainty_switch[i] == "k"){0
              } else{mean(`k_ci`)},
          cov_thick_mass = 
            if(uncertainty_switch[i] == "sample_thickness" || 
               uncertainty_switch[i] == "fine_soil_mass"){0
              } else{2 * cov(`sample_thickness_cm`, `fine_soil_mass_g`, 
                             use = "complete.obs") /
                  (`sample_thickness_ave2_cm` * `fine_soil_mass_ave2_g`)},
          cov_thick_volume = 
            if(uncertainty_switch[i] == "sample_thickness" || 
               uncertainty_switch[i] == "sample_volume"){0
              } else{2 * cov(`sample_thickness_cm`, `sample_volume_cm3`, 
                             use = "complete.obs") / 
                  (`sample_thickness_ave2_cm` * `sample_volume_ave2_cm3`)},
          cov_thick_som = 
            if(uncertainty_switch[i] == "sample_thickness" || 
               uncertainty_switch[i] == "som_content"){0
              } else{2 * cov(`sample_thickness_cm`, `som_content_perc`, 
                             use = "complete.obs") /
                  (`sample_thickness_ave2_cm` * `som_content_ave2_perc`)},
          cov_mass_volume = 
            if(uncertainty_switch[i] == "fine_soil_mass" || 
               uncertainty_switch[i] == "sample_volume"){0
              } else{2 * cov(`fine_soil_mass_g`, `sample_volume_cm3`, 
                             use = "complete.obs") / 
                  (`fine_soil_mass_ave2_g` * `sample_volume_ave2_cm3`)},
          cov_mass_som = 
            if(uncertainty_switch[i] == "fine_soil_mass" || 
               uncertainty_switch[i] == "som_content"){0
              } else{2 * cov(`fine_soil_mass_g`, `som_content_perc`,
                             use = "complete.obs") / 
                  (`fine_soil_mass_ave2_g` * `som_content_ave2_perc`)},
          cov_volume_som = 
            if(uncertainty_switch[i] == "sample_volume" || 
               uncertainty_switch[i] == "som_content"){0
              } else{2 * cov(`sample_volume_cm3`, `som_content_perc`, 
                             use = "complete.obs") /
                  (`sample_volume_ave2_cm3` * `som_content_ave2_perc`)},
          cov_mass_k = 
            if(uncertainty_switch[i] == "fine_soil_mass" || 
               uncertainty_switch[i] == "k"){0
              } else{2 * cov(`fine_soil_mass_g`, `k_ave`,
                             use = "complete.obs") / 
                  (`fine_soil_mass_ave2_g` * `k_ave2`)},
          cov_som_k = 
            if(uncertainty_switch[i] == "som_content" || 
               uncertainty_switch[i] == "k"){0
              } else{2 * cov(`som_content_perc`, `k_ave`, 
                             use = "complete.obs") / 
                  (`som_content_ave2_perc` * `k_ave2`)},
          cov_thick_k = 
            if(uncertainty_switch[i] == "sample_thickness" || 
               uncertainty_switch[i] == "k"){0
              } else{2 * cov(`sample_thickness_cm`, `k_ave`, 
                             use = "complete.obs") /
                  (`sample_thickness_ave2_cm` * `k_ave2`)},
          cov_volume_k = 
            if(uncertainty_switch[i] == "sample_volume" || 
               uncertainty_switch[i] == "k"){0
              } else{2 * cov(`sample_volume_cm3`, `k_ave`, 
                             use = "complete.obs") /
                  (`sample_volume_ave2_cm3` * `k_ave2`)},
          soc_content_ave2_perc = 
            `som_content_ave2_perc` / `k_ave2`,
          soc_content_dev2_perc = 
            `soc_content_ave2_perc` * 
            sqrt(((`som_content_dev2_perc` / `som_content_ave2_perc`) ^ 2) +
                   ((`k_err2` / `k_ave2`) ^ 2) -
                   cov_som_k),
          soc_content_err2_perc = 
            `soc_content_ave2_perc` * 
            sqrt(((`som_content_err2_perc` / `som_content_ave2_perc`) ^ 2) +
                   ((`k_err2` / `k_ave2`) ^ 2) -
                   cov_som_k / `sample_size`),
          soc_content_ci2_perc = 
            qt(0.975, df = `sample_size` - 1) * `soc_content_err2_perc`,
        ) %>%
        mutate(
          coe_thick_mass = `cov_thick_mass` / `sample_size`,
          coe_thick_volume = `cov_thick_volume` / `sample_size`,
          coe_thick_som = `cov_thick_som` / `sample_size`,
          coe_mass_volume = `cov_mass_volume` / `sample_size`,
          coe_mass_som = `cov_mass_som` / `sample_size`,
          coe_volume_som = `cov_volume_som` / `sample_size`,
          coe_mass_k = `cov_mass_k` / `sample_size`,
          coe_som_k = `cov_som_k` / `sample_size`,
          coe_thick_k = `cov_thick_k` / `sample_size`,
          coe_volume_k = `cov_volume_k` / `sample_size`,
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          soc_stock_ave2_Mgha = 
            `fine_soil_mass_ave2_g` * `som_content_ave2_perc` * 
            `sample_thickness_ave2_cm` / 
            (`k_ave2` * `sample_volume_ave2_cm3`),
          soc_stock_err2_Mgha = `soc_stock_ave2_Mgha` *
            sqrt(((`fine_soil_mass_err2_g` / `fine_soil_mass_ave2_g`) ^ 2) + 
                   ((`som_content_err2_perc` / `som_content_ave2_perc`) ^ 2) +
                   ((`k_err2` / `k_ave2`) ^ 2) +
                   ((`sample_thickness_err2_cm` / `sample_thickness_ave2_cm`) ^ 2) +
                   ((`sample_volume_err2_cm3` / `sample_volume_ave2_cm3`) ^ 2) +
                   `coe_thick_mass` - `coe_thick_volume` + `coe_thick_som` - 
                   `coe_mass_volume` + `coe_mass_som` - `coe_volume_som` -
                   `coe_mass_k` - `coe_som_k` - `coe_thick_k` + `coe_volume_k`),
          soc_stock_ci2_Mgha = `soc_stock_ave2_Mgha` *
            sqrt(((`fine_soil_mass_ci2_g` / `fine_soil_mass_ave2_g`) ^ 2) + 
                   ((`som_content_ci2_perc` / `som_content_ave2_perc`) ^ 2) +
                   ((`k_ci2` / `k_ave2`) ^ 2) +
                   ((`sample_thickness_ci2_cm` / `sample_thickness_ave2_cm`) ^ 2) +
                   ((`sample_volume_ci2_cm3` / `sample_volume_ave2_cm3`) ^ 2) +
                   qt(0.975, df = `sample_size` - 1) *
                   (`coe_thick_mass` - `coe_thick_volume` + `coe_thick_som` - 
                      `coe_mass_volume` + `coe_mass_som` - `coe_volume_som` -
                      `coe_mass_k` - `coe_som_k` - `coe_thick_k` + `coe_volume_k`)),
          soil_mineral_mass_ave2_g =
            (1 - 0.01 * `som_content_ave2_perc`) * `fine_soil_mass_ave2_g`,
          soil_mineral_mass_err2_g =
            sqrt(((`soil_mineral_mass_ave2_g` * `fine_soil_mass_err2_g` / 
                     `fine_soil_mass_ave2_g`) ^ 2) +
                   ((`soil_mineral_mass_ave2_g` * 0.01 * `som_content_err2_perc` / 
                       (1 - 0.01 * `som_content_ave2_perc`)) ^ 2) -
                   (`fine_soil_mass_ave2_g` * `som_content_ave2_perc` *
                      `soil_mineral_mass_ave2_g` * `coe_mass_som`)),
          soil_mineral_mass_ci2_g =
            sqrt(((`soil_mineral_mass_ave2_g` * `fine_soil_mass_ci2_g` / 
                     `fine_soil_mass_ave2_g`) ^ 2) +
                   ((`soil_mineral_mass_ave2_g` * 0.01 * `som_content_ci2_perc` / 
                       (1 - 0.01 * `som_content_ave2_perc`)) ^ 2) -
                   (qt(0.975, df = `sample_size` - 1) * 
                      `fine_soil_mass_ave2_g` * `som_content_ave2_perc` *
                      `soil_mineral_mass_ave2_g` * `coe_mass_som`)),
        )
      
      # calculate per composite stock and mineral mass
      data_stock_scov <- 
        data_merge1 %>%
        dplyr::mutate(soc_stock_Mgha = 
                        `fine_soil_mass_g` * `som_content_perc` *
                        `sample_thickness_cm` / (`k_ave` * `sample_volume_cm3`),
                      soil_mineral_mass_g = (1 - 0.01 * `som_content_perc`) *
                        `fine_soil_mass_g`)
    }
    
    # calculate double-summed covariances of depths, stocks and mineral masses for cumulative summation
    cov_sum <- data.frame()
    for(o in 1:length(unique(data_stock_scov$site))){
      for(n in 2:length(
        unique(data_stock_scov[data_stock_scov$site %in% 
                               unique(data_stock_scov$site)[o],]$increment))){
        for(m in 2:n-1){
          
          # get all depths at layers n and m
          thickness_cm_n <- data_stock_scov[
            data_stock_scov$site %in% unique(data_stock_scov$site)[o] &
              data_stock_scov$increment %in% unique(data_stock_scov[
                data_stock_scov$site %in% unique(data_stock_scov$site)[o],
              ]$increment)[n],]$sample_thickness_cm 
          thickness_cm_m <- data_stock_scov[
            data_stock_scov$site %in% unique(data_stock_scov$site)[o] &
              data_stock_scov$increment %in% unique(data_stock_scov[
                data_stock_scov$site %in% unique(data_stock_scov$site)[o],
              ]$increment)[m],]$sample_thickness_cm
          
          # get all stocks at layers n and m
          stock_Mgha_n <- data_stock_scov[
            data_stock_scov$site %in% unique(data_stock_scov$site)[o] &
              data_stock_scov$increment %in% unique(data_stock_scov[
                data_stock_scov$site %in% unique(data_stock_scov$site)[o],
              ]$increment)[n],]$soc_stock_Mgha
          stock_Mgha_m <- data_stock_scov[
            data_stock_scov$site %in% unique(data_stock_scov$site)[o] &
              data_stock_scov$increment %in% unique(data_stock_scov[
                data_stock_scov$site %in% unique(data_stock_scov$site)[o],
              ]$increment)[m],]$soc_stock_Mgha
          
          # get all mineral masses at layers n and m
          mineral_g_n <- data_stock_scov[
            data_stock_scov$site %in% unique(data_stock_scov$site)[o] &
              data_stock_scov$increment %in% unique(data_stock_scov[
                data_stock_scov$site %in% unique(data_stock_scov$site)[o],
              ]$increment)[n],]$soil_mineral_mass_g
          mineral_g_m <- data_stock_scov[
            data_stock_scov$site %in% unique(data_stock_scov$site)[o] &
              data_stock_scov$increment %in% unique(data_stock_scov[
                data_stock_scov$site %in% unique(data_stock_scov$site)[o],
              ]$increment)[m],]$soil_mineral_mass_g
          
          # select smallest sample size for covariance if different
          sample_size = pmin(length(stock_Mgha_n), length(stock_Mgha_m))
          
          # calculate covariances of depths, stocks and mienral masses between layers j and k
          cov_thickness <-
            cov(thickness_cm_n[1:sample_size],
                thickness_cm_m[1:sample_size],
                use = "complete.obs")
          cov_stock <-
            cov(stock_Mgha_n[1:sample_size],
                stock_Mgha_m[1:sample_size],
                use = "complete.obs")
          cov_mineral <-
            cov(mineral_g_n[1:sample_size],
                mineral_g_m[1:sample_size],
                use = "complete.obs")
          
          # if first iteration, create zero double-summed covariances for first layer
          if(n == 2 & m == 1){
            cov_sum <- 
              rbind(cov_sum, 
                    data.frame(site = unique(data_stock_scov$site)[o], 
                               increment = unique(data_stock_scov[
                                 data_stock_scov$site %in% 
                                   unique(data_stock_scov$site)[o],
                               ]$increment)[n-1],
                               sample_size = sample_size,
                               cov_thickness_cm2 = 0,
                               coe_thickness_cm2 = 0,
                               coci_thickness_cm2 = 0,
                               cov_stock_Mg2ha2 = 0,
                               coe_stock_Mg2ha2 = 0,
                               coci_stock_Mg2ha2 = 0,
                               cov_mineral_g2 = 0,
                               coe_mineral_g2 = 0,
                               coci_mineral_g2 = 0))
          }
          
          # if first iteration, bind double-summed covariances
          # else add each covariance to summed covariance
          if(m == 1){
            cov_sum <- 
              rbind(cov_sum, 
                    data.frame(
                      site = unique(data_stock_scov$site)[o], 
                      increment = unique(data_stock_scov[
                        data_stock_scov$site %in% 
                          unique(data_stock_scov$site)[o],]$increment)[n], 
                      sample_size = sample_size,
                      cov_thickness_cm2 = cov_thickness,
                      coe_thickness_cm2 = cov_thickness / sample_size,
                      coci_thickness_cm2 = cov_thickness * 
                        qt(0.975, df = sample_size - 1) / sample_size,
                      cov_stock_Mg2ha2 = cov_stock,
                      coe_stock_Mg2ha2 = cov_stock / sample_size,
                      coci_stock_Mg2ha2 = cov_stock * 
                        qt(0.975, df = sample_size - 1) / sample_size,
                      cov_mineral_g2 = cov_mineral,
                      coe_mineral_g2 = cov_mineral / sample_size,
                      coci_mineral_g2 = cov_mineral * 
                        qt(0.975, df = sample_size - 1) / sample_size))
          } else if(n > 2 & m > 2){
            cov_sum[n+(o-1)*5,4] <- cov_sum[n+(o-1)*5,4] + cov_thickness
            cov_sum[n+(o-1)*5,5] <- cov_sum[n+(o-1)*5,5] +
              (cov_thickness / sample_size)
            cov_sum[n+(o-1)*5,6] <- cov_sum[n+(o-1)*5,6] + 
              (cov_thickness * qt(0.975, df = sample_size - 1) / sample_size)
            cov_sum[n+(o-1)*5,7] <- cov_sum[n+(o-1)*5,7] + cov_stock
            cov_sum[n+(o-1)*5,8] <- cov_sum[n+(o-1)*5,8] + 
              (cov_stock / sample_size)
            cov_sum[n+(o-1)*5,9] <- cov_sum[n+(o-1)*5,9] +
              (cov_stock * qt(0.975, df = sample_size - 1) / sample_size)
            cov_sum[n+(o-1)*5,10] <- cov_sum[n+(o-1)*5,10] + cov_mineral
            cov_sum[n+(o-1)*5,11] <- cov_sum[n+(o-1)*5,11] + 
              (cov_mineral / sample_size)
            cov_sum[n+(o-1)*5,12] <- cov_sum[n+(o-1)*5,12] +
              (cov_mineral * qt(0.975, df = sample_size - 1) / sample_size)
          }
        }
      }
    }
    
    # calculate cumulative mineral mass per composite
    data_cum_mm <-
      data_stock_scov %>%
      dplyr::group_by(`site`, `core`) %>%
      dplyr::mutate(cum_soil_mineral_mass_g = cumsum(`soil_mineral_mass_g`))
    
    # get all treatment fields
    data_site_trt <- 
      unique((data_cum_mm %>%
                dplyr::filter(!(`site` %in% control_site)))$site)
    
    # create zero covariances for control field
    cov_mm_mineral_site <- 
      data.frame(site = control_site, 
                 increment = 
                   unique((data_cum_mm %>% 
                             dplyr::filter(!(`site` %in% control_site)))$increment),
                 sample_size = 
                   (data_cum_mm %>% 
                      dplyr::filter(`site` %in% control_site) %>% 
                      group_by(`site`, `increment`) %>% 
                      dplyr::summarise(sample_size = n()))$sample_size,
                 cov_mineral_site_g2 = 0,
                 coe_mineral_site_g2 = 0,
                 coci_mineral_site_g2 = 0)
    
    # loop over all cumulative mineral masses for cumulative mass covariance
    for(n in 1:length(data_site_trt)){
      for(m in 1:length(unique(unique(data_cum_mm[
        data_cum_mm$site %in% data_site_trt[n],]$increment)))){
        
        # calculate covariance of cumulative mineral mass between control and 
        # treatment sites: get all mineral masses at layers j and k
        mineral_g_ctrl <- data_cum_mm[
          data_cum_mm$site %in% control_site &
            data_cum_mm$increment %in% unique(data_cum_mm[
              data_cum_mm$site %in% control_site,]$increment)[m],
        ]$cum_soil_mineral_mass_g
        mineral_g_trt <- data_cum_mm[
          data_cum_mm$site %in% data_site_trt[n] &
            data_cum_mm$increment %in% unique(data_cum_mm[
              data_cum_mm$site %in% data_site_trt[n],]$increment)[m],
        ]$cum_soil_mineral_mass_g
        
        # select smallest sample size for covariance if different
        sample_size_mm = pmin(length(mineral_g_ctrl), length(mineral_g_trt))
        
        # calculate covariances of mineral masses between layers n and m
        cov_mineral_site <-
          cov(mineral_g_ctrl[1:sample_size],
              mineral_g_trt[1:sample_size],
              use = "complete.obs")
        
        # bind covariance to data frame
        cov_mm_mineral_site <- 
          rbind(cov_mm_mineral_site, 
                data.frame(
                  site = data_site_trt[n], 
                  increment = unique(data_cum_mm[
                    data_cum_mm$site %in% data_site_trt[n],]$increment)[m],
                  sample_size = sample_size_mm,
                  cov_mineral_site_g2 = cov_mineral_site,
                  coe_mineral_site_g2 = cov_mineral_site / sample_size_mm,
                  coci_mineral_site_g2 = cov_mineral_site * 
                    qt(0.975, df = sample_size - 1) / sample_size))
      }
    }
    
    # merge depths, stocks and mineral mass per layer by their covariances
    data_consolidate_1 <-
      merge(data_stock_precalculate, subset(cov_sum, select = -c(sample_size)), 
            by = c("site", "increment"))
    data_consolidate_2 <-
      merge(data_consolidate_1, cov_mm_mineral_site[,-c(3)], by = c("site", "increment"))
    
    # calculate cumulative depths, stocks and mineral masses
    data_calculation <-
      data_consolidate_2 %>%
      group_by(`site`) %>%
      mutate(cum_soil_mineral_mass_ave2_g = 
               cumsum(`soil_mineral_mass_ave2_g`),
             cum_soil_mineral_mass_err2_g = 
               sqrt(as.complex(cumsum(`soil_mineral_mass_err2_g` ^ 2) +
                                 2 * cumsum(`coe_mineral_g2`))),
             cum_soil_mineral_mass_ci2_g = 
               sqrt(as.complex(cumsum(`soil_mineral_mass_ci2_g` ^ 2) + 
                                 2 * cumsum(`coci_mineral_g2`))),
             cum_soil_depth_ave2_cm = 
               cumsum(`sample_thickness_ave2_cm`),
             cum_soil_depth_err2_cm = 
               sqrt(as.complex(cumsum(`sample_thickness_err2_cm` ^ 2) +
                                 2 * cumsum(`coe_thickness_cm2`))),
             cum_soil_depth_ci2_cm = 
               sqrt(as.complex(cumsum(`sample_thickness_ci2_cm` ^ 2) +
                                 2 * cumsum(`coci_thickness_cm2`))),
             cum_soc_stock_ave2_Mgha = 
               cumsum(`soc_stock_ave2_Mgha`),
             cum_soc_stock_err2_Mgha = 
               sqrt(as.complex(cumsum(`soc_stock_err2_Mgha` ^ 2) +
                                 2 * cumsum(`coe_stock_Mg2ha2`))),
             cum_soc_stock_ci2_Mgha = 
               sqrt(as.complex(cumsum(`soc_stock_ci2_Mgha` ^ 2) +
                                 2 * cumsum(`coci_stock_Mg2ha2`)))
      ) %>%
      arrange(`site`, `increment`)
    
    # return SOC stock and mineral mass calculation
    return(data_calculation) 
  }

### van Bemmelen factor selector ###############################################

# define van Bemmelen factor from input (if calculating mineral mass using SOC
# content)
van_bemmelen_selector <-
  function(data, k, guess_k){
    
    # create empty data frame to collect data models to
    model_consolidate <- data.frame()
    
    # fit straight line through SOM against SOC
    # run model calculation over all fields
    for(i in 1:length(unique(data$site))){
      
      # decide which approach to calculate/define van Bemmelen factor by
      if(k == "conventional"){
        
        # set conventional van Bemmelen factor value
        k_ave <- 1.724
        
        # create data frame of van Bemmelen model
        model <- 
          data.frame(
            study = unique(data$study), site = unique(data$site)[i],
            van_bemmelen_choice = k,
            k_ave = k_ave, k_err = 0, k_t = 0, k_p = 0, k_ci = 0,
            intercept_ave = 0, intercept_err = 0,
            intercept_t = 0, intercept_p = 0, intercept_ci = 0,
            rse = 0, dof = 0, r2 = 0, f_stat = 0)
        
      } else if(k == "no intercept"){
        
        # fit simple linear regression of straight line without intercept
        # SOM = k * SOC
        van_bemmelen_model <- 
          lm(`som_content_perc` ~ `soc_content_perc` + 0,
             data = data[data$site %in% unique(data$site)[i],])
        
        # get model summary and 95% CI
        van_bemmelen_summary <- summary(van_bemmelen_model)
        ci <- confint2(van_bemmelen_model, interval = "confidence")
        
        # get van Bemmelen factor including t- and p-values and CI
        k_ave = van_bemmelen_summary$coefficients[1,1]
        k_err = van_bemmelen_summary$coefficients[1,2]
        k_t   = van_bemmelen_summary$coefficients[1,3]
        k_p   = van_bemmelen_summary$coefficients[1,4]
        k_ci = (ci[2] - ci[1])
        
        # create data frame of van Bemmelen model
        model <- 
          data.frame(
            study = unique(data$study), site = unique(data$site)[i],
            van_bemmelen_choice = k,
            k_ave = k_ave, k_err = k_err, 
            k_t = k_t, k_p = k_p, k_ci = k_ci,
            intercept_ave = 0, intercept_err = 0,
            intercept_t = 0, intercept_p = 0, intercept_ci = 0,
            rse = van_bemmelen_summary$sigma, 
            dof = van_bemmelen_summary$df[2],
            r2 = van_bemmelen_summary$r.squared, 
            f_stat = van_bemmelen_summary$fstatistic["value"]) 
        
      } else if(k == "intercept"){
        
        # fit simple linear regression of straight line with intercept
        # SOM = k * SOC + SOM_intercept
        van_bemmelen_model <- 
          lm(`som_content_perc` ~ `soc_content_perc`,
             data = data[data$site %in% unique(data$site)[i],])
        
        # get model summary and 95% CI
        van_bemmelen_summary <- summary(van_bemmelen_model)
        ci <- confint2(van_bemmelen_model, interval = "confidence")
        
        # get van Bemmelen factor including t- and p-values and CI
        k_ave = van_bemmelen_summary$coefficients[2,1]
        k_err = van_bemmelen_summary$coefficients[2,2]
        k_t   = van_bemmelen_summary$coefficients[2,3]
        k_p   = van_bemmelen_summary$coefficients[2,4]
        k_ci  = (ci[2,2] - ci[2,1])
        
        # get SOM intercept including t- and p-values and CI
        intercept_ave = van_bemmelen_summary$coefficients[1,1]
        intercept_err = van_bemmelen_summary$coefficients[1,2]
        intercept_t   = van_bemmelen_summary$coefficients[1,3]
        intercept_p   = van_bemmelen_summary$coefficients[1,4]
        intercept_ci  = (ci[1,2] - ci[1,1])
        
        # create data frame of van Bemmelen model
        model <- 
          data.frame(
            study = unique(data$study), site = unique(data$site)[i],
            van_bemmelen_choice = k,
            k_ave = k_ave, k_err = k_err, k_t = k_t, k_p = k_p, k_ci = k_ci,
            intercept_ave = intercept_ave, intercept_err = intercept_err, 
            intercept_t = intercept_t, intercept_p = intercept_p,
            intercept_ci = intercept_ci,
            rse = van_bemmelen_summary$sigma, 
            dof = van_bemmelen_summary$df[2],
            r2 = van_bemmelen_summary$r.squared, 
            f_stat = van_bemmelen_summary$fstatistic["value"])
        
      } else if(k == "guess"){
        
        # set van Bemmelen factor to user-defined guess
        k_ave <- guess_k
        
        # create data frame of van Bemmelen model
        model <- 
          data.frame(
            study = unique(data$study), site = unique(data$site)[i],
            van_bemmelen_choice = k,
            k_ave = k_ave, k_err = 0, k_t = 0, k_p = 0, k_ci = 0,
            intercept_ave = 0, intercept_err = 0,
            intercept_t = 0, intercept_p = 0, intercept_ci = 0,
            rse = 0, dof = 0, r2 = 0, f_stat = 0)
        
      } else if(k == "both"){
        
        # create data frame of van Bemmelen model
        model <- 
          data.frame(
            study = unique(data$study), site = unique(data$site)[i],
            van_bemmelen_choice = k,
            k_ave = 0, k_err = 0, k_t = 0, k_p = 0, k_ci = 0,
            intercept_ave = 0, intercept_err = 0,
            intercept_t = 0, intercept_p = 0, intercept_ci = 0,
            rse = 0, dof = 0, r2 = 0, f_stat = 0)
        
      }
      
      # bind calculated van Bemmelen models to each field
      model_consolidate <- rbind(model_consolidate, model)
      
    }
    
    # return van Bemmelen models and coefficient averages
    return(model_consolidate)
  }

### Cumulative exponential distribution function (CEDF) ########################

# define average of cumulative exponential distribution function
cedf_average <-
  function(cum_soil_mineral_mass_ave_g, a, b){
    
    # calculate SOC stock
    cum_soc_stock_ave_Mgha <-
      a * (1 - exp(-b * cum_soil_mineral_mass_ave_g))
    
    # return average stock values
    return(cum_soc_stock_ave_Mgha)
  }

# define propagation error form of cumulative exponential distribution function
cedf_error <-
  function(cum_soc_stock_ave_Mgha, cum_soil_mineral_mass_ave_g, a, b, 
           cum_soil_mineral_mass_err_g, a_err, b_err, ab_cov){
    
    # calculate propagation deviation/error of SOC stock
    cum_soc_stock_err_Mgha = sqrt(as.complex(
      ((cum_soc_stock_ave_Mgha * a_err / a) ^ 2) +
        ((a - cum_soc_stock_ave_Mgha) ^ 2) * 
        (((cum_soil_mineral_mass_ave_g * b_err) ^ 2) + 
           ((cum_soil_mineral_mass_err_g * b) ^ 2)) -
        (2 * cum_soc_stock_ave_Mgha * (a - cum_soc_stock_ave_Mgha) / a) * 
        cum_soil_mineral_mass_ave_g * ab_cov))
    
    # return deviation/error of stock
    return(cum_soc_stock_err_Mgha)
  }

### Regression function fitting ################################################

# define function to fit selected ESM regression function to selected 
# treatment sites
regression_function_fitting <-
  function(data_treatment, data_control, regression_function, 
           curve_soil_mineral_mass){
    
    # create empty data frames
    curve_esm_consolidate <- data.frame()
    data_esm_consolidate  <- data.frame()
    
    # try core
    # if interpolation fails, skip core
    tryCatch({
      
      # select ESM function
      if(regression_function == "CEDF"){
        
        # create TRUE/FALSE switch data frame of parameters in FD2 method
        uncertainty_switch <- c("none", "a", "b")
        
        # for loop over all uncertainties included/excluded
        for(m in 1:length(uncertainty_switch)){
            
            # define CEDF expression for SOC stock and soil depth profile
            cedf_stock <- 
              cum_soc_stock_ave2_Mgha ~ 
              a * (1 - exp(-b * cum_soil_mineral_mass_ave2_g))
            cedf_depth <- 
              cum_soil_depth_ave2_cm ~ 
              d * (1 - exp(-c * cum_soil_mineral_mass_ave2_g))
            
            # fit function using non-linear least squares of CEDF
            # CEDF
            model_cedf_stock <- 
              nlsLM(cedf_stock, 
                    data = data_treatment, 
                    start = list(a = 180, b = 1e-3),
                    weights = as.numeric(`cum_soil_mineral_mass_err2_g` ^ 2) /
                      as.numeric(`cum_soc_stock_err2_Mgha` ^ 2))
            model_cedf_depth <- 
              nlsLM(cedf_depth, 
                    data = data_treatment, 
                    start = list(d = 100, c = 1e-3))
            
            # get summary, variance-covariance matrix and 95% CI of CEDF fit
            summary_cedf_stock <- summary(model_cedf_stock)
            summary_cedf_depth <- summary(model_cedf_depth)
            variance_covariance_stock <- vcov(model_cedf_stock)
            variance_covariance_depth <- vcov(model_cedf_depth)
            ci_cedf_stock <- confint2(model_cedf_stock, level = 0.95)
            ci_cedf_depth <- confint2(model_cedf_depth, level = 0.95)
            
            # get regression coefficient information
            a_ave = summary_cedf_stock$coefficients[1,1]
            a_err = if(uncertainty_switch[m] == "a"){0
              } else{summary_cedf_stock$coefficients[1,2]}
            a_t = summary_cedf_stock$coefficients[1,3]
            a_p = summary_cedf_stock$coefficients[1,4]
            a_ci = if(uncertainty_switch[m] == "a"){0
            } else{(ci_cedf_stock["a","97.5 %"] - 
                      ci_cedf_stock["a","2.5 %"]) / 2}
            b_ave = summary_cedf_stock$coefficients[2,1]
            b_err = if(uncertainty_switch[m] == "b"){0
              } else{summary_cedf_stock$coefficients[2,2]}
            b_t = summary_cedf_stock$coefficients[2,3]
            b_p = summary_cedf_stock$coefficients[2,4]
            b_ci = if(uncertainty_switch[m] == "b"){0
            } else{(ci_cedf_stock["b","97.5 %"] - 
                      ci_cedf_stock["b","2.5 %"]) / 2}
            ab_cov = 
              if(uncertainty_switch[m] == "a" || 
                 uncertainty_switch[m] == "b"){0
                } else{variance_covariance_stock[1,2]}
            
            d_ave = summary_cedf_depth$coefficients[1,1]
            d_err = summary_cedf_depth$coefficients[1,2]
            d_t = summary_cedf_depth$coefficients[1,3]
            d_p = summary_cedf_depth$coefficients[1,4]
            d_ci = (ci_cedf_depth["d","97.5 %"] - 
                      ci_cedf_depth["d","2.5 %"]) / 2
            c_ave = summary_cedf_depth$coefficients[2,1]
            c_err = summary_cedf_depth$coefficients[2,2]
            c_t = summary_cedf_depth$coefficients[2,3]
            c_p = summary_cedf_depth$coefficients[2,4]
            c_ci = (ci_cedf_depth["c","97.5 %"] - 
                      ci_cedf_depth["c","2.5 %"]) / 2
            dc_cov = variance_covariance_depth[1,2]
            
            # create curve-formatted function SOC stock and soil depth
            curve_soc_stock_ave_Mgha <- 
              cedf_average(
                cum_soil_mineral_mass_ave_g = curve_soil_mineral_mass, 
                a = a_ave, b = b_ave)
            curve_soc_stock_err_Mgha <- 
              cedf_error(
                cum_soc_stock_ave_Mgha = curve_soc_stock_ave_Mgha, 
                cum_soil_mineral_mass_ave_g = curve_soil_mineral_mass, 
                a = a_ave, b = b_ave, cum_soil_mineral_mass_err_g = 0, 
                a_err = a_err, b_err = b_err, ab_cov = ab_cov)
            curve_soc_stock_ci_Mgha <- 
              cedf_error(
                cum_soc_stock_ave_Mgha = curve_soc_stock_ave_Mgha, 
                cum_soil_mineral_mass_ave_g = curve_soil_mineral_mass, 
                a = a_ave, b = b_ave, cum_soil_mineral_mass_err_g = 0, 
                a_err = a_ci, b_err = b_ci, 
                ab_cov = (ab_cov * a_ci * 2 / 
                            if(uncertainty_switch[m] == "a"){1
                              } else{a_err}))
            
            curve_soil_depth_ave_cm <- 
              cedf_average(
                cum_soil_mineral_mass_ave_g = curve_soil_mineral_mass, 
                a = d_ave, b = c_ave)
            curve_soil_depth_err_cm <- 
              cedf_error(
                cum_soc_stock_ave_Mgha = curve_soil_depth_ave_cm, 
                cum_soil_mineral_mass_ave_g = curve_soil_mineral_mass, 
                a = d_ave, b = c_ave, cum_soil_mineral_mass_err_g = 0, 
                a_err = d_err, b_err = c_err, ab_cov = dc_cov)
            curve_soil_depth_ci_cm <- 
              cedf_error(
                cum_soc_stock_ave_Mgha = curve_soil_depth_ave_cm, 
                cum_soil_mineral_mass_ave_g = curve_soil_mineral_mass, 
                a = d_ave, b = c_ave, cum_soil_mineral_mass_err_g = 0, 
                a_err = d_ci, b_err = c_ci,
                ab_cov = (dc_cov * d_ci * 2 / d_err))
            
            # calculate ESM-adjusted C stock from function and control averages
            esm_soc_stock_ave_Mgha <- 
              cedf_average(
                cum_soil_mineral_mass_ave_g = data_control$cum_soil_mineral_mass_ave2_g, 
                a = a_ave, b = b_ave)
            esm_soc_stock_err_Mgha <- 
              cedf_error(
                cum_soc_stock_ave_Mgha = esm_soc_stock_ave_Mgha, 
                cum_soil_mineral_mass_ave_g = data_control$cum_soil_mineral_mass_ave2_g,
                a = a_ave, b = b_ave,
                cum_soil_mineral_mass_err_g = data_control$cum_soil_mineral_mass_err2_g, 
                a_err = a_err, b_err = b_err, ab_cov = ab_cov)
            esm_soc_stock_ci_Mgha <- 
              cedf_error(
                cum_soc_stock_ave_Mgha = esm_soc_stock_ave_Mgha, 
                cum_soil_mineral_mass_ave_g = data_control$cum_soil_mineral_mass_ave2_g,
                a = a_ave, b = b_ave,
                cum_soil_mineral_mass_err_g = data_control$cum_soil_mineral_mass_ci2_g, 
                a_err = a_ci, b_err = b_ci, 
                ab_cov = (ab_cov * a_ci * 2 / 
                            if(uncertainty_switch[m] == "a"){1
                              } else{a_err}))
            
            esm_soil_depth_ave_cm <- 
              cedf_average(
                cum_soil_mineral_mass_ave_g = data_control$cum_soil_mineral_mass_ave2_g, 
                a = d_ave, b = c_ave)
            esm_soil_depth_err_cm <- 
              cedf_error(
                cum_soc_stock_ave_Mgha = esm_soil_depth_ave_cm,
                cum_soil_mineral_mass_ave_g = data_control$cum_soil_mineral_mass_ave2_g,
                a = d_ave, b = c_ave,
                cum_soil_mineral_mass_err_g = data_control$cum_soil_mineral_mass_err2_g, 
                a_err = d_err, b_err = c_err, ab_cov = dc_cov)
            esm_soil_depth_ci_cm <- 
              cedf_error(
                cum_soc_stock_ave_Mgha = esm_soil_depth_ave_cm,
                cum_soil_mineral_mass_ave_g = data_control$cum_soil_mineral_mass_ave2_g,
                a = d_ave, b = c_ave,
                cum_soil_mineral_mass_err_g = data_control$cum_soil_mineral_mass_err2_g, 
                a_err = d_ci, b_err = c_ci, ab_cov = (dc_cov * d_ci * 2 / d_err))
            
            # collect curve data to data frame
            curve_esm <-
              data.frame(uncertainty_off = uncertainty_switch[m],
                         site = unique(data_treatment$site),
                         cum_soil_mineral_mass_ave_g = curve_soil_mineral_mass,
                         cum_soil_depth_ave_cm = curve_soil_depth_ave_cm,
                         cum_soil_depth_err_cm = curve_soil_depth_err_cm,
                         cum_soil_depth_ci_cm = curve_soil_depth_ci_cm,
                         cum_soc_stock_ave_Mgha = curve_soc_stock_ave_Mgha,
                         cum_soc_stock_err_Mgha = curve_soc_stock_err_Mgha,
                         cum_soc_stock_ci_Mgha = curve_soc_stock_ci_Mgha)
            
            # collect ESM-adjusted information to data frame
            data_esm <- 
              data.frame(uncertainty_off = uncertainty_switch[m],
                         site = unique(data_treatment$site),
                         increment = data_control$increment,
                         sample_size = data_control$sample_size,
                         cum_soil_mineral_mass_ave_g = data_control$cum_soil_mineral_mass_ave2_g,
                         cum_soil_mineral_mass_err_g = data_control$cum_soil_mineral_mass_err2_g,
                         cum_soil_mineral_mass_ci_g = data_control$cum_soil_mineral_mass_ci2_g,
                         cum_soil_depth_ave_cm = esm_soil_depth_ave_cm,
                         cum_soil_depth_err_cm = esm_soil_depth_err_cm,
                         cum_soil_depth_ci_cm = esm_soil_depth_ci_cm,
                         cum_soc_stock_ave_Mgha = esm_soc_stock_ave_Mgha,
                         cum_soc_stock_err_Mgha = esm_soc_stock_err_Mgha,
                         cum_soc_stock_ci_Mgha = esm_soc_stock_ci_Mgha)

            # consolidate curves and ESM-adjusted points
            curve_esm_consolidate <- rbind(curve_esm_consolidate, curve_esm)
            data_esm_consolidate <- rbind(data_esm_consolidate, data_esm)
        }
        
      } else if(regression_function == "MCS"){
        
        # add origin to data set of treatment site
        data_treatment_modify <-
          rbind(data.frame(
            site = unique(data_treatment$site),
            cum_soil_mineral_mass_ave2_g = 0,
            cum_soil_mineral_mass_err2_g = as.complex(0),
            cum_soil_mineral_mass_ci2_g = as.complex(0),
            cum_soil_depth_ave2_cm = 0,
            cum_soil_depth_err2_cm = as.complex(0),
            cum_soil_depth_ci2_cm = as.complex(0),
            cum_soc_stock_ave2_Mgha = 0,
            cum_soc_stock_err2_Mgha = as.complex(0),
            cum_soc_stock_ci2_Mgha = as.complex(0)),
            data_treatment[c('site',
                             'cum_soil_mineral_mass_ave2_g',
                             'cum_soil_mineral_mass_err2_g',
                             'cum_soil_mineral_mass_ci2_g',
                             'cum_soil_depth_ave2_cm',
                             'cum_soil_depth_err2_cm',
                             'cum_soil_depth_ci2_cm',
                             'cum_soc_stock_ave2_Mgha',
                             'cum_soc_stock_err2_Mgha',
                             'cum_soc_stock_ci2_Mgha')])
        
        # fit monotonic cubic spline to data
        model_mcs_stock <-
          splinefun(
            x = data_treatment_modify$cum_soil_mineral_mass_ave2_g,
            y = data_treatment_modify$cum_soc_stock_ave2_Mgha,
            method = "hyman")
        model_mcs_depth <-
          splinefun(
            x = data_treatment_modify$cum_soil_mineral_mass_ave2_g,
            y = data_treatment_modify$cum_soil_depth_ave2_cm,
            method = "hyman")
        
        # create curve-formatted function SOC stock
        curve_soc_stock_ave_Mgha <- 
          model_mcs_stock(curve_soil_mineral_mass, deriv = 0L)
        curve_soc_stock_err_Mgha <- 0
        curve_soc_stock_ci_Mgha <- 0
        curve_soil_depth_ave_cm <- 
          model_mcs_depth(curve_soil_mineral_mass, deriv = 0L)
        curve_soil_depth_err_cm <- 0
        curve_soil_depth_ci_cm <- 0
        
        # calculate ESM-adjusted C stock from function and control averages
        esm_soc_stock_ave_Mgha <- 
          model_mcs_stock(data_control$cum_soil_mineral_mass_ave2_g, 
                          deriv = 0L)
        esm_soc_stock_err_Mgha <- sqrt(as.complex(
          ((model_mcs_stock(data_control$cum_soil_mineral_mass_ave2_g,
                            deriv = 1L) *
              data_control$cum_soil_mineral_mass_err2_g) ^ 2) +
            ((model_mcs_stock(data_treatment$cum_soil_mineral_mass_ave2_g,
                              deriv = 1L) *
                data_treatment$cum_soil_mineral_mass_err2_g) ^ 2) +
            (2 * model_mcs_stock(data_control$cum_soil_mineral_mass_ave2_g,
                                 deriv = 1L) *
               model_mcs_stock(data_treatment$cum_soil_mineral_mass_ave2_g,
                               deriv = 1L) *
               data_treatment$coe_mineral_g2)))
        esm_soc_stock_ci_Mgha <- sqrt(as.complex(
          ((model_mcs_stock(data_control$cum_soil_mineral_mass_ave2_g,
                            deriv = 1L) *
              data_control$cum_soil_mineral_mass_ci2_g) ^ 2) +
            ((model_mcs_stock(data_treatment$cum_soil_mineral_mass_ave2_g,
                              deriv = 1L) *
                data_treatment$cum_soil_mineral_mass_ci2_g) ^ 2) +
            (2 * model_mcs_stock(data_control$cum_soil_mineral_mass_ave2_g,
                                 deriv = 1L) *
               model_mcs_stock(data_treatment$cum_soil_mineral_mass_ave2_g,
                               deriv = 1L) *
               data_treatment$coci_mineral_g2)))
        
        esm_soil_depth_ave_cm <- 
          model_mcs_depth(data_control$cum_soil_mineral_mass_ave2_g, 
                          deriv = 0L)
        esm_soil_depth_err_cm <- sqrt(as.complex(
          ((model_mcs_depth(data_control$cum_soil_mineral_mass_ave2_g,
                            deriv = 1L) *
              data_control$cum_soil_mineral_mass_err2_g) ^ 2) +
            ((model_mcs_depth(data_treatment$cum_soil_mineral_mass_ave2_g,
                              deriv = 1L) *
                data_treatment$cum_soil_mineral_mass_err2_g) ^ 2) +
            (2 * model_mcs_depth(data_control$cum_soil_mineral_mass_ave2_g,
                                 deriv = 1L) *
               model_mcs_depth(data_treatment$cum_soil_mineral_mass_ave2_g,
                               deriv = 1L) *
               data_treatment$coe_mineral_g2)))
        esm_soil_depth_ci_cm <- sqrt(as.complex(
          ((model_mcs_depth(data_control$cum_soil_mineral_mass_ave2_g,
                            deriv = 1L) *
              data_control$cum_soil_mineral_mass_ci2_g) ^ 2) +
            ((model_mcs_depth(data_treatment$cum_soil_mineral_mass_ave2_g,
                              deriv = 1L) *
                data_treatment$cum_soil_mineral_mass_ci2_g) ^ 2) +
            (2 * model_mcs_depth(data_control$cum_soil_mineral_mass_ave2_g,
                                 deriv = 1L) *
               model_mcs_depth(data_treatment$cum_soil_mineral_mass_ave2_g,
                               deriv = 1L) *
               data_treatment$coci_mineral_g2)))
        
        # collect curve data to data frame
        curve_esm <-
          data.frame(uncertainty_off = "none",
                     site = unique(data_treatment$site),
                     cum_soil_mineral_mass_ave_g = curve_soil_mineral_mass,
                     cum_soil_depth_ave_cm = curve_soil_depth_ave_cm,
                     cum_soil_depth_err_cm = curve_soil_depth_err_cm,
                     cum_soil_depth_ci_cm = curve_soil_depth_ci_cm,
                     cum_soc_stock_ave_Mgha = curve_soc_stock_ave_Mgha,
                     cum_soc_stock_err_Mgha = curve_soc_stock_err_Mgha,
                     cum_soc_stock_ci_Mgha = curve_soc_stock_ci_Mgha)
        
        # collect ESM-adjusted information to data frame
        data_esm <- 
          data.frame(uncertainty_off = "none",
                     site = unique(data_treatment$site),
                     increment = data_control$increment,
                     sample_size = data_control$sample_size,
                     cum_soil_mineral_mass_ave_g = data_control$cum_soil_mineral_mass_ave2_g,
                     cum_soil_mineral_mass_err_g = data_control$cum_soil_mineral_mass_err2_g,
                     cum_soil_mineral_mass_ci_g = data_control$cum_soil_mineral_mass_ci2_g,
                     cum_soil_depth_ave_cm = esm_soil_depth_ave_cm,
                     cum_soil_depth_err_cm = esm_soil_depth_err_cm,
                     cum_soil_depth_ci_cm = esm_soil_depth_ci_cm,
                     cum_soc_stock_ave_Mgha = esm_soc_stock_ave_Mgha,
                     cum_soc_stock_err_Mgha = esm_soc_stock_err_Mgha,
                     cum_soc_stock_ci_Mgha = esm_soc_stock_ci_Mgha)
        
        # consolidate curves and ESM-adjusted points
        curve_esm_consolidate <- rbind(curve_esm_consolidate, curve_esm)
        data_esm_consolidate <- rbind(data_esm_consolidate, data_esm)
        
      } else if(regression_function == "LIN"){
        
        # add origin to data set of treatment group
        data_treatment_modify <-
          rbind(data.frame(
            site = unique(data_treatment$site),
            cum_soil_mineral_mass_ave2_g = 0,
            cum_soil_mineral_mass_err2_g = as.complex(0),
            cum_soil_mineral_mass_ci2_g = as.complex(0),
            cum_soil_depth_ave2_cm = 0,
            cum_soil_depth_err2_cm = as.complex(0),
            cum_soil_depth_ci2_cm = as.complex(0),
            cum_soc_stock_ave2_Mgha = 0,
            cum_soc_stock_err2_Mgha = as.complex(0),
            cum_soc_stock_ci2_Mgha = as.complex(0)),
            data_treatment[c('site',
                             'cum_soil_mineral_mass_ave2_g',
                             'cum_soil_mineral_mass_err2_g',
                             'cum_soil_mineral_mass_ci2_g',
                             'cum_soil_depth_ave2_cm',
                             'cum_soil_depth_err2_cm',
                             'cum_soil_depth_ci2_cm',
                             'cum_soc_stock_ave2_Mgha',
                             'cum_soc_stock_err2_Mgha',
                             'cum_soc_stock_ci2_Mgha')])
        
        # create linear approximation curve
        for(m in 1:length(curve_soil_mineral_mass)){
          
          # if control mineral mass is beyond interpolation range
          if(curve_soil_mineral_mass[m] > 
             max(data_treatment_modify$cum_soil_mineral_mass_ave2_g)){
            
            # calculate gradient and intercept of straight line between
            # last two data points past approx interpolation to extrapolate
            # C stock
            gradient_stock_curve <- 
              (as.numeric((data_treatment_modify %>% 
                             slice(n()))[1,'cum_soc_stock_ave2_Mgha']) -
                 as.numeric((data_treatment_modify %>%
                               slice(n()-1))[1,'cum_soc_stock_ave2_Mgha'])) /
              (as.numeric((data_treatment_modify %>% 
                             slice(n()))[1,'cum_soil_mineral_mass_ave2_g']) -
                 as.numeric((data_treatment_modify %>% 
                               slice(n()-1))[1,'cum_soil_mineral_mass_ave2_g']))
            
            intercept_stock_curve <-
              as.numeric((data_treatment_modify %>% 
                            slice(n()))[1,'cum_soc_stock_ave2_Mgha']) -
              gradient_stock_curve *
              as.numeric((data_treatment_modify %>% 
                            slice(n()))[1,'cum_soil_mineral_mass_ave2_g'])
            
            model_stock_ave_curve <-
              data.frame(curve_stock_ave =
                           intercept_stock_curve +
                           gradient_stock_curve * 
                           curve_soil_mineral_mass[m])[1,1]
            
            # calculate gradient and intercept of straight line between
            # last two data points past approx interpolation to extrapolate
            # depth
            gradient_depth_curve <- 
              (as.numeric((data_treatment_modify %>%
                             slice(n()))[1,'cum_soil_depth_ave2_cm']) -
                 as.numeric((data_treatment_modify %>% 
                               slice(n()-1))[1,'cum_soil_depth_ave2_cm'])) /
              (as.numeric((data_treatment_modify %>% 
                             slice(n()))[1,'cum_soil_mineral_mass_ave2_g']) -
                 as.numeric((data_treatment_modify %>% 
                               slice(n()-1))[1,'cum_soil_mineral_mass_ave2_g']))
            
            intercept_depth_curve <-
              as.numeric((data_treatment_modify %>% 
                            slice(n()))[1,'cum_soil_depth_ave2_cm']) -
              gradient_depth_curve *
              as.numeric((data_treatment_modify %>% 
                            slice(n()))[1,'cum_soil_mineral_mass_ave2_g'])
            
            model_depth_ave_curve <-
              data.frame(curve_depth_ave =
                           intercept_depth_curve +
                           gradient_depth_curve * 
                           curve_soil_mineral_mass[m])[1,1]
            
          } else{
            
            # interpolate curve using linear approximation
            model_stock_ave_curve <- 
              data.frame(
                curve_stock_ave =
                  approx(data_treatment_modify$cum_soil_mineral_mass_ave2_g,
                         data_treatment_modify$cum_soc_stock_ave2_Mgha,
                         xout = curve_soil_mineral_mass[m])$y)[1,1]
            model_depth_ave_curve <- 
              data.frame(
                curve_depth_ave =
                  approx(data_treatment_modify$cum_soil_mineral_mass_ave2_g,
                         data_treatment_modify$cum_soil_depth_ave2_cm,
                         xout = curve_soil_mineral_mass[m])$y)[1,1]
          }
          
          # stack over all k iterations to create curve data frame
          if(m == 1){
            curve_soc_stock_ave_Mgha <- model_stock_ave_curve
            curve_soc_stock_err_Mgha <- 0
            curve_soc_stock_ci_Mgha <- 0
            curve_soil_depth_ave_cm <- model_depth_ave_curve
            curve_soil_depth_err_cm <- 0
            curve_soil_depth_ci_cm <- 0
          } else {
            curve_soc_stock_ave_Mgha <- rbind(curve_soc_stock_ave_Mgha,
                                              model_stock_ave_curve)
            curve_soc_stock_err_Mgha <- rbind(curve_soc_stock_err_Mgha, 0)
            curve_soc_stock_ci_Mgha <- rbind(curve_soc_stock_ci_Mgha, 0)
            curve_soil_depth_ave_cm <- rbind(curve_soil_depth_ave_cm, 
                                             model_depth_ave_curve)
            curve_soil_depth_err_cm <- rbind(curve_soil_depth_err_cm, 0)
            curve_soil_depth_ci_cm <- rbind(curve_soil_depth_ci_cm, 0)
          }
        }
        
        # fit control mass to treatment using linear interpolation
        for(m in 1:nrow(data_control)){
          
          # if control mineral mass is beyond interpolation range
          if(data_control$cum_soil_mineral_mass_ave2_g[m] > 
             max(data_treatment_modify$cum_soil_mineral_mass_ave2_g)){
            
            # calculate gradient and intercept of straight line between
            # last two data points past approx interpolation to extrapolate
            # SOC stock
            gradient_stock_Mghag <- 
              (as.numeric((data_treatment_modify %>%
                             slice(n()))[1,'cum_soc_stock_ave2_Mgha']) -
                 as.numeric((data_treatment_modify %>% 
                               slice(n()-1))[1,'cum_soc_stock_ave2_Mgha'])) /
              (as.numeric((data_treatment_modify %>% 
                             slice(n()))[1,'cum_soil_mineral_mass_ave2_g']) -
                 as.numeric((data_treatment_modify %>% 
                               slice(n()-1))[1,'cum_soil_mineral_mass_ave2_g']))
            
            intercept_stock_Mgha <-
              as.numeric((data_treatment_modify %>%
                            slice(n()))[1,'cum_soc_stock_ave2_Mgha']) -
              gradient_stock_Mghag *
              as.numeric((data_treatment_modify %>%
                            slice(n()))[1,'cum_soil_mineral_mass_ave2_g'])
            
            lin_soc_stock_ave_Mgha <-
              data.frame(intercept_stock_Mgha +
                           gradient_stock_Mghag * 
                           as.numeric(data_control[
                             m,'cum_soil_mineral_mass_ave2_g']))[1,1]
            
            # calculate gradient and intercept of straight line between
            # last two data points past approx interpolation to extrapolate
            # depth
            gradient_depth_cmg <- 
              (as.numeric((data_treatment_modify %>% 
                             slice(n()))[1,'cum_soil_depth_ave2_cm']) -
                 as.numeric((data_treatment_modify %>%
                               slice(n()-1))[1,'cum_soil_depth_ave2_cm'])) /
              (as.numeric((data_treatment_modify %>% 
                             slice(n()))[1,'cum_soil_mineral_mass_ave2_g']) -
                 as.numeric((data_treatment_modify %>%
                               slice(n()-1))[1,'cum_soil_mineral_mass_ave2_g']))
            
            intercept_depth_cm <-
              as.numeric((data_treatment_modify %>% 
                            slice(n()))[1,'cum_soil_depth_ave2_cm']) -
              gradient_depth_cmg *
              as.numeric((data_treatment_modify %>% 
                            slice(n()))[1,'cum_soil_mineral_mass_ave2_g'])
            
            lin_soil_depth_ave_cm <-
              data.frame(intercept_depth_cm +
                           gradient_depth_cmg * 
                           as.numeric(data_control[
                             m,'cum_soil_mineral_mass_ave2_g']))[1,1]
            
            # find upper index of treatment point with total mineral mass
            j2 <- which(
              max(data_treatment_modify$cum_soil_mineral_mass_ave2_g) == 
                data_treatment_modify$cum_soil_mineral_mass_ave2_g)
            
            # find lower index of treatment point
            j1 = j2 - 1
            
          } else{
            
            # calculate interpolated C stock of linear approximation for control
            # mineral mass point
            lin_soc_stock_ave_Mgha <- 
              data.frame(approx(
                data_treatment_modify$cum_soil_mineral_mass_ave2_g,
                data_treatment_modify$cum_soc_stock_ave2_Mgha,
                xout = data_control$cum_soil_mineral_mass_ave2_g[m])$y)[1,1]
            
            lin_soil_depth_ave_cm <- 
              data.frame(approx(
                data_treatment_modify$cum_soil_mineral_mass_ave2_g,
                data_treatment_modify$cum_soil_depth_ave2_cm,
                xout = data_control$cum_soil_mineral_mass_ave2_g[m])$y)[1,1]
            
            # find lower and upper indices of treatment field relative to 
            # control mineral mass
            j12 <- find_bound_indices(
              vector = data_treatment_modify$cum_soil_mineral_mass_ave2_g, 
              value = data_control$cum_soil_mineral_mass_ave2_g[m])
            j1 <- j12[1]
            j2 <- j12[2]
          }
          
          # calculate linear approximation terms and then PE of ESM-adjusted C stock and depth
          # note: assume that correlation between each interpolated point with approx()
          # is equal to 1:
          #   appropriate as each line intercepts between 2 pts in a cumulative trend 
          #   (2 pts have always abs(corr) = 1 and "cumulative" is never declining to be -1)
          gs = (data_treatment_modify$cum_soc_stock_ave2_Mgha[j2] - 
                  data_treatment_modify$cum_soc_stock_ave2_Mgha[j1]) / 
            (data_treatment_modify$cum_soil_mineral_mass_ave2_g[j2] - 
               data_treatment_modify$cum_soil_mineral_mass_ave2_g[j1])
          
          gd = (data_treatment_modify$cum_soil_depth_ave2_cm[j2] - 
                  data_treatment_modify$cum_soil_depth_ave2_cm[j1]) /
            (data_treatment_modify$cum_soil_mineral_mass_ave2_g[j2] - 
               data_treatment_modify$cum_soil_mineral_mass_ave2_g[j1])
          
          r = (data_control$cum_soil_mineral_mass_ave2_g[m] - 
                 data_treatment_modify$cum_soil_mineral_mass_ave2_g[j1]) / 
            (data_treatment_modify$cum_soil_mineral_mass_ave2_g[j2] - 
               data_treatment_modify$cum_soil_mineral_mass_ave2_g[j1])
          
          lin_soc_stock_err_Mgha <- sqrt(as.complex(
            (gs ^ 2) * 
              ((data_control$cum_soil_mineral_mass_err2_g[m] ^ 2) +
            ((r * data_treatment_modify$cum_soil_mineral_mass_err2_g[j2]) ^ 2) +
              (((1 - r) * 
                data_treatment_modify$cum_soil_mineral_mass_err2_g[j1]) ^ 2)) +
              ((r * data_treatment_modify$cum_soc_stock_err2_Mgha[j2]) ^ 2) +
              (((1 - r) * 
                  data_treatment_modify$cum_soc_stock_err2_Mgha[j1]) ^ 2) +
              2 * (gs ^ 2) * (
                (r * (1 - r) * 
                   data_treatment_modify$cum_soil_mineral_mass_err2_g[j2] *
                   data_treatment_modify$cum_soil_mineral_mass_err2_g[j1]) -
                  (r * data_control$cum_soil_mineral_mass_err2_g[m] * 
                     data_treatment_modify$cum_soil_mineral_mass_err2_g[j2]) -
                  ((1 - r) * data_control$cum_soil_mineral_mass_err2_g[m] *
                     data_treatment_modify$cum_soil_mineral_mass_err2_g[j1])) +
              2 * gs * r * (
                (data_treatment_modify$cum_soc_stock_err2_Mgha[j2] *
                   data_control$cum_soil_mineral_mass_err2_g[m]) -
                  (r * data_treatment_modify$cum_soc_stock_err2_Mgha[j2] *
                     data_treatment_modify$cum_soil_mineral_mass_err2_g[j2]) - 
                  ((1 - r) * data_treatment_modify$cum_soc_stock_err2_Mgha[j2] *
                     data_treatment_modify$cum_soil_mineral_mass_err2_g[j1])) +
              2 * gs * (1 - r) * (
                (data_treatment_modify$cum_soc_stock_err2_Mgha[j1] *
                   data_control$cum_soil_mineral_mass_err2_g[m]) -
                  (r * data_treatment_modify$cum_soc_stock_err2_Mgha[j1] *
                     data_treatment_modify$cum_soil_mineral_mass_err2_g[j2]) - 
                  ((1 - r) * data_treatment_modify$cum_soc_stock_err2_Mgha[j1] *
                     data_treatment_modify$cum_soil_mineral_mass_err2_g[j1])) +
              2 * r * (1 - r) * 
              data_treatment_modify$cum_soc_stock_err2_Mgha[j2] *
              data_treatment_modify$cum_soc_stock_err2_Mgha[j1]))
          
          lin_soc_stock_ci_Mgha <- sqrt(as.complex(
            (gs ^ 2) * 
              ((data_control$cum_soil_mineral_mass_ci2_g[m] ^ 2) +
            ((r * data_treatment_modify$cum_soil_mineral_mass_ci2_g[j2]) ^ 2) +
            (((1 - r) * 
                data_treatment_modify$cum_soil_mineral_mass_ci2_g[j1]) ^ 2)) +
              ((r * data_treatment_modify$cum_soc_stock_ci2_Mgha[j2]) ^ 2) +
              (((1 - r) * 
                  data_treatment_modify$cum_soc_stock_ci2_Mgha[j1]) ^ 2) +
              2 * (gs ^ 2) * (
                (r * (1 - r) * 
                   data_treatment_modify$cum_soil_mineral_mass_ci2_g[j2] *
                   data_treatment_modify$cum_soil_mineral_mass_ci2_g[j1]) -
                  (r * data_control$cum_soil_mineral_mass_ci2_g[m] * 
                     data_treatment_modify$cum_soil_mineral_mass_ci2_g[j2]) -
                  ((1 - r) * data_control$cum_soil_mineral_mass_ci2_g[m] *
                     data_treatment_modify$cum_soil_mineral_mass_ci2_g[j1])) +
              2 * gs * r * (
                (data_treatment_modify$cum_soc_stock_ci2_Mgha[j2] *
                   data_control$cum_soil_mineral_mass_ci2_g[m]) -
                  (r * data_treatment_modify$cum_soc_stock_ci2_Mgha[j2] *
                     data_treatment_modify$cum_soil_mineral_mass_ci2_g[j2]) - 
                  ((1 - r) * data_treatment_modify$cum_soc_stock_ci2_Mgha[j2] *
                     data_treatment_modify$cum_soil_mineral_mass_ci2_g[j1])) +
              2 * gs * (1 - r) * (
                (data_treatment_modify$cum_soc_stock_ci2_Mgha[j1] *
                   data_control$cum_soil_mineral_mass_ci2_g[m]) -
                  (r * data_treatment_modify$cum_soc_stock_ci2_Mgha[j1] *
                     data_treatment_modify$cum_soil_mineral_mass_ci2_g[j2]) - 
                  ((1 - r) * data_treatment_modify$cum_soc_stock_ci2_Mgha[j1] *
                     data_treatment_modify$cum_soil_mineral_mass_ci2_g[j1])) +
              2 * r * (1 - r) * 
              data_treatment_modify$cum_soc_stock_ci2_Mgha[j2] *
              data_treatment_modify$cum_soc_stock_ci2_Mgha[j1]))
          
          lin_soil_depth_err_cm <- sqrt(as.complex(
            (gd ^ 2) * 
              ((data_control$cum_soil_mineral_mass_err2_g[m] ^ 2) +
            ((r * data_treatment_modify$cum_soil_mineral_mass_err2_g[j2]) ^ 2) +
              (((1 - r) * 
                data_treatment_modify$cum_soil_mineral_mass_err2_g[j1]) ^ 2)) +
              ((r * data_treatment_modify$cum_soil_depth_err2_cm[j2]) ^ 2) +
              (((1 - r) * 
                  data_treatment_modify$cum_soil_depth_err2_cm[j1]) ^ 2) +
              2 * (gd ^ 2) * (
                (r * (1 - r) * 
                   data_treatment_modify$cum_soil_mineral_mass_err2_g[j2] * 
                   data_treatment_modify$cum_soil_mineral_mass_err2_g[j1]) -
                  (r * data_control$cum_soil_mineral_mass_err2_g[m] *
                     data_treatment_modify$cum_soil_mineral_mass_err2_g[j2]) -
                  ((1 - r) * data_control$cum_soil_mineral_mass_err2_g[m] *
                     data_treatment_modify$cum_soil_mineral_mass_err2_g[j1])) +
              2 * gd * r * (
                (data_treatment_modify$cum_soil_depth_err2_cm[j2] * 
                   data_control$cum_soil_mineral_mass_err2_g[m]) -
                  (r * data_treatment_modify$cum_soil_depth_err2_cm[j2] *
                     data_treatment_modify$cum_soil_mineral_mass_err2_g[j2]) - 
                  ((1 - r) * data_treatment_modify$cum_soil_depth_err2_cm[j2] *
                     data_treatment_modify$cum_soil_mineral_mass_err2_g[j1])) +
              2 * gd * (1 - r) * (
                (data_treatment_modify$cum_soil_depth_err2_cm[j1] *
                   data_control$cum_soil_mineral_mass_err2_g[m]) -
                  (r * data_treatment_modify$cum_soil_depth_err2_cm[j1] *
                     data_treatment_modify$cum_soil_mineral_mass_err2_g[j2]) - 
                  ((1 - r) * data_treatment_modify$cum_soil_depth_err2_cm[j1] *
                     data_treatment_modify$cum_soil_mineral_mass_err2_g[j1])) +
              2 * r * (1 - r) * 
              data_treatment_modify$cum_soil_depth_err2_cm[j2] *
              data_treatment_modify$cum_soil_depth_err2_cm[j1]))
          
          lin_soil_depth_ci_cm <- sqrt(as.complex(
            (gd ^ 2) * 
              ((data_control$cum_soil_mineral_mass_ci2_g[m] ^ 2) +
            ((r * data_treatment_modify$cum_soil_mineral_mass_ci2_g[j2]) ^ 2) +
              (((1 - r) * 
                  data_treatment_modify$cum_soil_mineral_mass_ci2_g[j1]) ^ 2)) +
              ((r * data_treatment_modify$cum_soil_depth_ci2_cm[j2]) ^ 2) +
              (((1 - r) * 
                  data_treatment_modify$cum_soil_depth_ci2_cm[j1]) ^ 2) +
              2 * (gd ^ 2) * (
                (r * (1 - r) * 
                   data_treatment_modify$cum_soil_mineral_mass_ci2_g[j2] * 
                   data_treatment_modify$cum_soil_mineral_mass_ci2_g[j1]) -
                  (r * data_control$cum_soil_mineral_mass_ci2_g[m] *
                     data_treatment_modify$cum_soil_mineral_mass_ci2_g[j2]) -
                  ((1 - r) * data_control$cum_soil_mineral_mass_ci2_g[m] *
                     data_treatment_modify$cum_soil_mineral_mass_ci2_g[j1])) +
              2 * gd * r * (
                (data_treatment_modify$cum_soil_depth_ci2_cm[j2] * 
                   data_control$cum_soil_mineral_mass_ci2_g[m]) -
                  (r * data_treatment_modify$cum_soil_depth_ci2_cm[j2] *
                     data_treatment_modify$cum_soil_mineral_mass_ci2_g[j2]) - 
                  ((1 - r) * data_treatment_modify$cum_soil_depth_ci2_cm[j2] *
                     data_treatment_modify$cum_soil_mineral_mass_ci2_g[j1])) +
              2 * gd * (1 - r) * (
                (data_treatment_modify$cum_soil_depth_ci2_cm[j1] *
                   data_control$cum_soil_mineral_mass_ci2_g[m]) -
                  (r * data_treatment_modify$cum_soil_depth_ci2_cm[j1] *
                     data_treatment_modify$cum_soil_mineral_mass_ci2_g[j2]) - 
                  ((1 - r) * data_treatment_modify$cum_soil_depth_ci2_cm[j1] *
                     data_treatment_modify$cum_soil_mineral_mass_ci2_g[j1])) +
              2 * r * (1 - r) * 
              data_treatment_modify$cum_soil_depth_ci2_cm[j2] *
              data_treatment_modify$cum_soil_depth_ci2_cm[j1]))
          
          # stack over all m iterations to create curve data frame
          if(m == 1){
            esm_soc_stock_ave_Mgha <- lin_soc_stock_ave_Mgha
            esm_soc_stock_err_Mgha <- lin_soc_stock_err_Mgha
            esm_soc_stock_ci_Mgha <- lin_soc_stock_ci_Mgha
            esm_soil_depth_ave_cm <- lin_soil_depth_ave_cm
            esm_soil_depth_err_cm <- lin_soil_depth_err_cm
            esm_soil_depth_ci_cm <- lin_soil_depth_ci_cm
            
          } else {
            esm_soc_stock_ave_Mgha <- rbind(
              esm_soc_stock_ave_Mgha, lin_soc_stock_ave_Mgha)
            esm_soc_stock_err_Mgha <- rbind(
              esm_soc_stock_err_Mgha, lin_soc_stock_err_Mgha)
            esm_soc_stock_ci_Mgha <- rbind(
              esm_soc_stock_ci_Mgha, lin_soc_stock_ci_Mgha)
            esm_soil_depth_ave_cm <- rbind(
              esm_soil_depth_ave_cm, lin_soil_depth_ave_cm)
            esm_soil_depth_err_cm <- rbind(
              esm_soil_depth_err_cm, lin_soil_depth_err_cm)
            esm_soil_depth_ci_cm <- rbind(
              esm_soil_depth_ci_cm, lin_soil_depth_ci_cm)
          }
          
        }
        
        # collect curve data to data frame
        curve_esm <-
          data.frame(uncertainty_off = "none",
                     site = unique(data_treatment$site),
                     cum_soil_mineral_mass_ave_g = curve_soil_mineral_mass,
                     cum_soil_depth_ave_cm = curve_soil_depth_ave_cm,
                     cum_soil_depth_err_cm = curve_soil_depth_err_cm,
                     cum_soil_depth_ci_cm = curve_soil_depth_ci_cm,
                     cum_soc_stock_ave_Mgha = curve_soc_stock_ave_Mgha,
                     cum_soc_stock_err_Mgha = curve_soc_stock_err_Mgha,
                     cum_soc_stock_ci_Mgha = curve_soc_stock_ci_Mgha)
        
        # collect ESM-adjusted information to data frame
        data_esm <- 
          data.frame(uncertainty_off = "none",
                     site = unique(data_treatment$site),
                     increment = data_control$increment,
                     sample_size = data_control$sample_size,
                     cum_soil_mineral_mass_ave_g = data_control$cum_soil_mineral_mass_ave2_g,
                     cum_soil_mineral_mass_err_g = data_control$cum_soil_mineral_mass_err2_g,
                     cum_soil_mineral_mass_ci_g = data_control$cum_soil_mineral_mass_ci2_g,
                     cum_soil_depth_ave_cm = esm_soil_depth_ave_cm,
                     cum_soil_depth_err_cm = esm_soil_depth_err_cm,
                     cum_soil_depth_ci_cm = esm_soil_depth_ci_cm,
                     cum_soc_stock_ave_Mgha = esm_soc_stock_ave_Mgha,
                     cum_soc_stock_err_Mgha = esm_soc_stock_err_Mgha,
                     cum_soc_stock_ci_Mgha = esm_soc_stock_ci_Mgha)
        
        # consolidate curves and ESM-adjusted points
        curve_esm_consolidate <- rbind(curve_esm_consolidate, curve_esm)
        data_esm_consolidate <- rbind(data_esm_consolidate, data_esm)
        
      } else{
        
        # print error
        message("Warning:")
        message("ESM correction only works for LIN, MCS or CEDF.")
        message("No regression fit performed.")
      }
    }, 
    error = function(cond) {
      message("")
      NA},
    warning = function(cond) {
      message("")
      NULL},
    finally = {
      message("")
    })
    
    return(list(curve_esm_consolidate, data_esm_consolidate))
  }

### Manual ANOVA calculators ###################################################

# define ANOVA calculator comparing control and all treatments
manual_anova_calculator_control <-
  function(data, data_final, control_field, y_ave, y_ci, y_err){
    
    # reduce data to only unique rows and filter to selected farm, field and SOM-LOI assumption
    data_reduced <- 
      data_final %>%
      distinct(`study`, `site`, `increment`, `method`, `spp`, `soil_content`,
               `regression_function`, .keep_all = TRUE) %>%
      arrange(`study`, `site`, `increment`, `method`, `spp`, `soil_content`,
              `regression_function`)
    
    ###
    #
    summary_anova_results <- list()
    label_anova_results <- data.frame()
    
    for(i in 1:length(unique(data_reduced$increment))){
      
      #
      data_anova <-
        data_reduced %>%
        dplyr::mutate(comparison = 
                        paste(`site`, ".", `method`, ".", `regression_function`, 
                              sep = "")) %>%
        dplyr::select(c("study", "site", "increment", "comparison", "method",
                        "soil_content", "k", "regression_function",
                        y_ave, y_ci, y_err, "sample_size")) %>%
        dplyr::filter(`increment` %in% unique(.$increment)[i])
      
      #
      for(o in 2:nrow(data_anova)){
        
        #
        anova_summary <-
          data_anova[c(1,o),] %>%
          dplyr::mutate(
            cum_y_var = (`sample_size`) * .[y_err] ^ 2,
            cum_y_ave_between = sum(`sample_size` * .[y_ave]) / 
              sum(`sample_size`),
            sse = sum(`sample_size` * `cum_y_var`),
            ssr = sum(`sample_size` * (.[y_ave] - `cum_y_ave_between`) ^ 2),
            sst = sse + ssr,
            dof_lower = n() - 1,
            dof_upper = sum(`sample_size`) - n(),
            mse = sse / dof_upper,
            msr = ssr / dof_lower,
            f_observed = msr / mse,
            p_value = pf(f_observed, dof_lower, dof_upper, lower.tail = FALSE)
          )
        
        #
        anova_outcomes <-
          data.frame(
            `Type` = c(anova_summary[2,]$comparison, "Residuals"),
            `Df` = c(anova_summary[2,]$dof_lower, anova_summary[2,]$dof_upper), 
            `Sum Sq` = c(anova_summary[2,]$ssr, anova_summary[2,]$sse),
            `Mean Sq` = c(anova_summary[2,]$msr, anova_summary[2,]$mse), 
            `F value` = c(anova_summary[2,]$f_observed, NA),
            `Pr(>F)` = c(anova_summary[2,]$p_value, NA),
            check.names = FALSE) %>%
          mutate(`significance` = `Pr(>F)`) %>%
          mutate(`significance` = 
                   replace(`significance`, `significance` == 0, 5)) %>%
          mutate(`significance` = 
                   replace(`significance`, 
                           `significance` <= 1 & `significance` > 0.05, 6)) %>%
          mutate(`significance` =
                   replace(`significance`, 
                           `significance` <= 0.05 & 
                             `significance` > 0.01, 7)) %>%
          mutate(`significance` = 
                   replace(`significance`, 
                           `significance` <= 0.01 &
                             `significance` > 0.001, 8)) %>%
          mutate(`significance` =
                   replace(`significance`, `significance` <= 0.001, 9)) %>%
          mutate(`significance` = as.character(`significance`)) %>%
          mutate(`significance` = 
                   replace(`significance`, `significance` %in% "5", "NA")) %>%
          mutate(`significance` =
                   replace(`significance`, `significance` %in% "6", "ns")) %>%
          mutate(`significance` =
                   replace(`significance`, `significance` %in% "7" , "*")) %>%
          mutate(`significance` =
                   replace(`significance`, `significance` %in% "8" , "**")) %>%
          mutate(`significance` =
                   replace(`significance`, `significance` %in% "9" , "***"))
        
        #
        summary_anova_results[[paste0(anova_summary[2,]$comparison)]] <- 
          anova_outcomes
        
        anova_format <- data.frame(anova_summary[2,],
                                   p = anova_outcomes$significance[1])
        
        label_anova_results <- 
          bind_rows(label_anova_results, anova_format)
      }
    }
    
    anova_final_results <-
      rbind(label_anova_results %>%
              dplyr::select(
                c("study", "site", "increment", "method", "regression_function",
                  y_ave, y_ci, y_err, "mse", "msr", "f_observed", "p")),
            data_final %>%
              dplyr::filter(`site` %in% control_site) %>%
              dplyr::select(c("study", "site", "increment", "method", 
                              "regression_function", y_ave, y_ci, y_err)) %>%
              dplyr::mutate(mse = NA, msr = NA, f_observed = NA, p = "")) %>%
      arrange(`study`, `site`, `method`, `regression_function`, `increment`)
    
    return(anova_final_results)
  }

# define ANOVA calculator comparing all ESM corrections to FD calculation
manual_anova_calculator_esm <-
  function(data, data_final, treatment_site, y_ave, y_ci, y_err){
    
    # reduce data to only unique rows and filter to selected farm, field and SOM-LOI assumption
    data_reduced <-
      data_final %>%
      distinct(`study`, `site`, `increment`, `method`, `soil_content`, 
               `regression_function`, .keep_all = TRUE) %>%
      arrange(`study`, `site`, `increment`, `method`, `soil_content`,
              `regression_function`)%>%
      dplyr::filter(`site` %in% c(treatment_site)) %>%
      arrange(desc(`method`))
    
    ###
    #
    summary_anova_results <- list()
    label_anova_results <- data.frame()
    
    for(i in 1:length(unique(data_reduced$increment))){
      
      #
      data_anova <-
        data_reduced %>%
        dplyr::mutate(comparison = 
                        paste(`site`, ".", `method`, ".", `regression_function`,
                              sep = "")) %>%
        dplyr::select(c("study", "site", "increment", "comparison", "method",
                        "soil_content", "k", "regression_function",
                        y_ave, y_err, "sample_size")) %>%
        dplyr::filter(`increment` %in% unique(.$increment)[i]) %>%
        dplyr::arrange(match(`method`, c("FD", "ESM")))
      
      
      #
      for(o in 2:nrow(data_anova)){
        
        #
        anova_summary <-
          data_anova[c(1,o),] %>%
          dplyr::mutate(
            cum_y_var = (`sample_size`) * .[y_err] ^ 2,
            cum_y_ave_between = sum(`sample_size` * .[y_ave]) /
              sum(`sample_size`),
            sse = sum(`sample_size` * `cum_y_var`),
            ssr = sum(`sample_size` * (.[y_ave] - `cum_y_ave_between`) ^ 2),
            sst = sse + ssr,
            dof_lower = n() - 1,
            dof_upper = sum(`sample_size`) - n(),
            mse = sse / dof_upper,
            msr = ssr / dof_lower,
            f_observed = msr / mse,
            p_value = pf(f_observed, dof_lower, dof_upper, lower.tail = FALSE))
        
        #
        anova_outcomes <-
          data.frame(
            `Type` = c(anova_summary[2,]$comparison, "Residuals"),
            `Df` = c(anova_summary[2,]$dof_lower, anova_summary[2,]$dof_upper), 
            `Sum Sq` = c(anova_summary[2,]$ssr, anova_summary[2,]$sse),
            `Mean Sq` = c(anova_summary[2,]$msr, anova_summary[2,]$mse), 
            `F value` = c(anova_summary[2,]$f_observed, NA),
            `Pr(>F)` = c(anova_summary[2,]$p_value, NA), 
            check.names = FALSE) %>%
          mutate(`significance` = `Pr(>F)`) %>%
          mutate(`significance` = 
                   replace(`significance`, `significance` == 0, 5)) %>%
          mutate(`significance` = 
                   replace(`significance`, 
                           `significance` <= 1 & `significance` > 0.05, 6)) %>%
          mutate(`significance` =
                   replace(`significance`, 
                           `significance` <= 0.05 & 
                             `significance` > 0.01, 7)) %>%
          mutate(`significance` = 
                   replace(`significance`, 
                           `significance` <= 0.01 &
                             `significance` > 0.001, 8)) %>%
          mutate(`significance` =
                   replace(`significance`, `significance` <= 0.001, 9)) %>%
          mutate(`significance` = as.character(`significance`)) %>%
          mutate(`significance` = 
                   replace(`significance`, `significance` %in% "5", "NA")) %>%
          mutate(`significance` =
                   replace(`significance`, `significance` %in% "6", "ns")) %>%
          mutate(`significance` =
                   replace(`significance`, `significance` %in% "7" , "*")) %>%
          mutate(`significance` =
                   replace(`significance`, `significance` %in% "8" , "**")) %>%
          mutate(`significance` =
                   replace(`significance`, `significance` %in% "9" , "***"))
        
        summary_anova_results[[paste0(anova_summary[2,]$comparison)]] <- 
          anova_outcomes
        
        anova_format <- data.frame(anova_summary[2,],
                                   p = anova_outcomes$significance[1])
        
        label_anova_results <- 
          bind_rows(label_anova_results, anova_format)
      }
    }
    
    anova_final_results <-
      rbind(label_anova_results %>%
              dplyr::select(
                c("study", "site", "increment", "method", "regression_function",
                  y_ave, y_err, "p")),
            data_final2 %>%
              dplyr::filter(`site` %in% treatment_site &
                              `regression_function` %in% "NA") %>%
              dplyr::select(c("study", "site", "increment", "method",
                              "regression_function", y_ave, y_err,)) %>%
              dplyr::mutate(p = "")) %>%
      arrange(`study`, `site`, `method`, `regression_function`, `increment`)
    
    return(anova_final_results)
  }

### Plotting ESM-corrected graphs ##############################################

# define visual analyser to plot ESM results
visual_analyser <-
  function(data_results, control_site, k, soil_content){
    
    # separate data 
    data_final <- data.frame(data_results[1])
    data_curve <- data.frame(data_results[2])
    
    # convert complex uncertainties to numeric
    data_final[c("cum_soil_mineral_mass_err_g", "cum_soil_mineral_mass_ci_g",
                 "cum_soil_depth_err_cm", "cum_soil_depth_ci_cm",
                 "cum_soc_stock_err_Mgha", "cum_soc_stock_ci_Mgha")] <- 
      apply(data_final[
        c("cum_soil_mineral_mass_err_g", "cum_soil_mineral_mass_ci_g",
          "cum_soil_depth_err_cm", "cum_soil_depth_ci_cm",
          "cum_soc_stock_err_Mgha", "cum_soc_stock_ci_Mgha")], 2,
            function(x) as.numeric(x))
    
    data_curve[c("cum_soil_depth_err_cm", "cum_soil_depth_ci_cm",
                 "cum_soc_stock_err_Mgha", "cum_soc_stock_ci_Mgha")] <- 
      apply(data_curve[
        c("cum_soil_depth_err_cm", "cum_soil_depth_ci_cm",
          "cum_soc_stock_err_Mgha", "cum_soc_stock_ci_Mgha")], 2,
            function(x) as.numeric(x))
    
    # get all treatment field labels
    treatment_sites <- 
      unique(data_final[!(data_final$site %in% c(control_site)),]$site)
    
    # loop analysis through all treatment sites
    for(m in 1:length(treatment_sites)){
      
      # plot FD-based and ESM-corrected SOC stocks and soil depths from LIN
      if("LIN" %in% data_final$regression_function){
        
        plot_esm_lin_stock <-
          ggplot() +
          geom_line(
            data = data_curve %>% 
              dplyr::filter(`site` %in% c(treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("LIN")),
            aes(x = `cum_soil_mineral_mass_ave_g`,
                y = `cum_soc_stock_ave_Mgha`),
            color = "#1b7837", linetype = 1) + 
          geom_errorbar(
            data = data_final %>% 
              dplyr::filter(`site` %in% c(control_site, treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("NA", "LIN")),
            aes(x = `cum_soil_mineral_mass_ave_g`,
                ymin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
                ymax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`),
                color = paste(`site`, ": ", `method`, sep = "")),
                        width = 20) +
          geom_errorbarh(
            data = data_final %>% 
              dplyr::filter(`site` %in% c(control_site, treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("LIN")),
            aes(y = `cum_soc_stock_ave_Mgha`,
                xmin = (`cum_soil_mineral_mass_ave_g` - 
                          `cum_soil_mineral_mass_err_g`),
                xmax = (`cum_soil_mineral_mass_ave_g` + 
                          `cum_soil_mineral_mass_err_g`),
                color = paste(`site`, ": ", `method`, sep = "")), 
            height = 5) +
          geom_point(
            data = data_final %>% 
              dplyr::filter(`site` %in% c(control_site, treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("NA", "LIN")),
            aes(x = `cum_soil_mineral_mass_ave_g`,
                y = `cum_soc_stock_ave_Mgha`,
                color = paste(`site`, ": ", `method`, sep = "")), 
            size = 2, shape = 15) +
          theme_bw() +
          theme(legend.position = "bottom") +
          scale_x_continuous(limits = 
              c(0, max(data_final$cum_soil_mineral_mass_ave_g) * 1.05)) +
          scale_y_continuous(
            limits = c(0, max(data_final$cum_soc_stock_ave_Mgha + 
                                data_final$cum_soc_stock_err_Mgha) * 1.05)) +
          scale_color_manual(
            values = c("#af8dc3", "#7fbf7b", "#762a83", "white"),
            labels = c(paste(control_site, ":", "FD"),
                       paste(treatment_sites[m], ":", "FD"),
                       paste(treatment_sites[m], ":", "ESM"), "")) +
          labs(x = expression(
            paste("Cumulative soil mineral mass, ", m[M], " (g)", sep = "")),
               y = expression(
                 paste("Cumulative SOC stock, ", S[C], " (Mg C/ha)", sep = "")),
               color = "",
               title = "LIN")
        
        plot_esm_lin_depth <-
          ggplot() +
          geom_line(
            data = data_curve %>% 
              dplyr::filter(`site` %in% c(treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("LIN")),
            aes(x = `cum_soil_mineral_mass_ave_g`,
                y = `cum_soil_depth_ave_cm`),
            color = "#1b7837", linetype = 1) + 
          geom_errorbar(
            data = data_final %>% 
              dplyr::filter(`site` %in% c(control_site, treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("NA", "LIN")),
            aes(x = `cum_soil_mineral_mass_ave_g`,
                ymin = (`cum_soil_depth_ave_cm` - `cum_soil_depth_err_cm`),
                ymax = (`cum_soil_depth_ave_cm` + `cum_soil_depth_err_cm`),
                color = paste(`site`, ": ", `method`, sep = "")),
            width = 20) +
          geom_errorbarh(
            data = data_final %>% 
              dplyr::filter(`site` %in% c(control_site, treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("LIN")),
            aes(y = `cum_soil_depth_ave_cm`,
                xmin = (`cum_soil_mineral_mass_ave_g` - 
                          `cum_soil_mineral_mass_err_g`),
                xmax = (`cum_soil_mineral_mass_ave_g` + 
                          `cum_soil_mineral_mass_err_g`),
                color = paste(`site`, ": ", `method`, sep = "")), 
            height = 5) +
          geom_point(
            data = data_final %>% 
              dplyr::filter(`site` %in% c(control_site, treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("NA", "LIN")),
            aes(x = `cum_soil_mineral_mass_ave_g`,
                y = `cum_soil_depth_ave_cm`,
                color = paste(`site`, ": ", `method`, sep = "")), 
            size = 2, shape = 15) +
          theme_bw() +
          theme(legend.position = "bottom") +
          scale_x_continuous(limits = 
              c(0, max(data_final$cum_soil_mineral_mass_ave_g) * 1.05)) +
          scale_y_continuous(
            limits = c(0, 
                       max(data_final$cum_soil_depth_ave_cm + 
                             data_final$cum_soil_depth_err_cm) * 1.05)) +
          scale_color_manual(
            values = c("#af8dc3", "#7fbf7b", "#762a83", "white"),
            labels = c(paste(control_site, ":", "FD"),
                       paste(treatment_sites[m], ":", "FD"),
                       paste(treatment_sites[m], ":", "ESM"), "")) +
          labs(x = expression(
            paste("Cumulative soil mineral mass, ", m[M], " (g)", sep = "")),
               y = expression(
                 paste("Cumulative soil depth, ", z, " (cm)", sep = "")),
               color = "",
               title = "LIN")
      }
      
      # plot FD-based and ESM-corrected SOC stocks and soil depths from MCS
      if("MCS" %in% data_final$regression_function){
        
        plot_esm_mcs_stock <-
          ggplot() +
          geom_line(
            data = data_curve %>% 
              dplyr::filter(`site` %in% c(treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("MCS")),
            aes(x = `cum_soil_mineral_mass_ave_g`,
                y = `cum_soc_stock_ave_Mgha`),
            color = "#1b7837", linetype = 1) + 
          geom_errorbar(
            data = data_final %>% 
              dplyr::filter(`site` %in% c(control_site, treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("NA", "MCS")),
            aes(x = `cum_soil_mineral_mass_ave_g`,
                ymin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
                ymax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`),
                color = paste(`site`, ": ", `method`, sep = "")),
            width = 20) +
          geom_errorbarh(
            data = data_final %>% 
              dplyr::filter(`site` %in% c(control_site, treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("MCS")),
            aes(y = `cum_soc_stock_ave_Mgha`,
                xmin = (`cum_soil_mineral_mass_ave_g` - 
                          `cum_soil_mineral_mass_err_g`),
                xmax = (`cum_soil_mineral_mass_ave_g` + 
                          `cum_soil_mineral_mass_err_g`),
                color = paste(`site`, ": ", `method`, sep = "")), 
            height = 5) +
          geom_point(
            data = data_final %>% 
              dplyr::filter(`site` %in% c(control_site, treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("NA", "MCS")),
            aes(x = `cum_soil_mineral_mass_ave_g`,
                y = `cum_soc_stock_ave_Mgha`,
                color = paste(`site`, ": ", `method`, sep = "")), 
            size = 2, shape = 15) +
          theme_bw() +
          theme(legend.position = "bottom") +
          scale_x_continuous(
            limits = 
              c(0, max(data_final$cum_soil_mineral_mass_ave_g) * 1.05)) +
          scale_y_continuous(
            limits = c(0, max(data_final$cum_soc_stock_ave_Mgha + 
                                data_final$cum_soc_stock_err_Mgha) * 1.05)) +
          scale_color_manual(
            values = c("#af8dc3", "#7fbf7b", "#762a83", "white"),
            labels = c(paste(control_site, ":", "FD"),
                       paste(treatment_sites[m], ":", "FD"),
                       paste(treatment_sites[m], ":", "ESM"), "")) +
          labs(x = expression(
            paste("Cumulative soil mineral mass, ", m[M], " (g)", sep = "")),
            y = expression(
              paste("Cumulative SOC stock, ", S[C], " (Mg C/ha)", sep = "")),
            color = "",
            title = "MCS")
        
        plot_esm_mcs_depth <-
          ggplot() +
          geom_line(
            data = data_curve %>% 
              dplyr::filter(`site` %in% c(treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("MCS")),
            aes(x = `cum_soil_mineral_mass_ave_g`,
                y = `cum_soil_depth_ave_cm`),
            color = "#1b7837", linetype = 1) + 
          geom_errorbar(
            data = data_final %>% 
              dplyr::filter(`site` %in% c(control_site, treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("NA", "MCS")),
            aes(x = `cum_soil_mineral_mass_ave_g`,
                ymin = (`cum_soil_depth_ave_cm` - `cum_soil_depth_err_cm`),
                ymax = (`cum_soil_depth_ave_cm` + `cum_soil_depth_err_cm`),
                color = paste(`site`, ": ", `method`, sep = "")),
            width = 20) +
          geom_errorbarh(
            data = data_final %>% 
              dplyr::filter(`site` %in% c(control_site, treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("MCS")),
            aes(y = `cum_soil_depth_ave_cm`,
                xmin = (`cum_soil_mineral_mass_ave_g` - 
                          `cum_soil_mineral_mass_err_g`),
                xmax = (`cum_soil_mineral_mass_ave_g` + 
                          `cum_soil_mineral_mass_err_g`),
                color = paste(`site`, ": ", `method`, sep = "")), 
            height = 5) +
          geom_point(
            data = data_final %>% 
              dplyr::filter(`site` %in% c(control_site, treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("NA", "MCS")),
            aes(x = `cum_soil_mineral_mass_ave_g`,
                y = `cum_soil_depth_ave_cm`,
                color = paste(`site`, ": ", `method`, sep = "")), 
            size = 2, shape = 15) +
          theme_bw() +
          theme(legend.position = "bottom") +
          scale_x_continuous(limits = 
                               c(0, max(data_final$cum_soil_mineral_mass_ave_g) * 1.05)) +
          scale_y_continuous(
            limits = c(0, 
                       max(data_final$cum_soil_depth_ave_cm + 
                             data_final$cum_soil_depth_err_cm) * 1.05)) +
          scale_color_manual(
            values = c("#af8dc3", "#7fbf7b", "#762a83", "white"),
            labels = c(paste(control_site, ":", "FD"),
                       paste(treatment_sites[m], ":", "FD"),
                       paste(treatment_sites[m], ":", "ESM"), "")) +
          labs(x = expression(
            paste("Cumulative soil mineral mass, ", m[M], " (g)", sep = "")),
            y = expression(
              paste("Cumulative soil depth, ", z, " (cm)", sep = "")),
            color = "",
            title = "MCS")
      }
      
      # plot FD-based and ESM-corrected SOC stocks and soil depths from CEDF
      if("CEDF" %in% data_final$regression_function){
        
        plot_esm_cedf_stock <-
          ggplot() +
          geom_ribbon(
            data = data_curve %>% 
              dplyr::filter(`site` %in% c(treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("CEDF")),
            aes(x = `cum_soil_mineral_mass_ave_g`,
                ymin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
                ymax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`)),
            color = "#1b7837", linetype = 0, alpha = 0.2) +
          geom_line(
            data = data_curve %>% 
              dplyr::filter(`site` %in% c(treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("CEDF")),
            aes(x = `cum_soil_mineral_mass_ave_g`,
                y = `cum_soc_stock_ave_Mgha`),
            color = "#1b7837", linetype = 1) + 
          geom_errorbar(
            data = data_final %>% 
              dplyr::filter(`site` %in% c(control_site, treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("NA", "CEDF")),
            aes(x = `cum_soil_mineral_mass_ave_g`,
                ymin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
                ymax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`),
                color = paste(`site`, ": ", `method`, sep = "")),
            width = 20) +
          geom_errorbarh(
            data = data_final %>% 
              dplyr::filter(`site` %in% c(control_site, treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("CEDF")),
            aes(y = `cum_soc_stock_ave_Mgha`,
                xmin = (`cum_soil_mineral_mass_ave_g` - 
                          `cum_soil_mineral_mass_err_g`),
                xmax = (`cum_soil_mineral_mass_ave_g` + 
                          `cum_soil_mineral_mass_err_g`),
                color = paste(`site`, ": ", `method`, sep = "")), 
            height = 5) +
          geom_point(
            data = data_final %>% 
              dplyr::filter(`site` %in% c(control_site, treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("NA", "CEDF")),
            aes(x = `cum_soil_mineral_mass_ave_g`,
                y = `cum_soc_stock_ave_Mgha`,
                color = paste(`site`, ": ", `method`, sep = "")), 
            size = 2, shape = 15) +
          theme_bw() +
          theme(legend.position = "bottom") +
          scale_x_continuous(
            limits = 
              c(0, max(data_final$cum_soil_mineral_mass_ave_g) * 1.05)) +
          scale_y_continuous(
            limits = c(0, max(data_final$cum_soc_stock_ave_Mgha + 
                                data_final$cum_soc_stock_err_Mgha) * 1.05)) +
          scale_color_manual(
            values = c("#af8dc3", "#7fbf7b", "#762a83", "white"),
            labels = c(paste(control_site, ":", "FD"),
                       paste(treatment_sites[m], ":", "FD"),
                       paste(treatment_sites[m], ":", "ESM"), "")) +
          labs(x = expression(
            paste("Cumulative soil mineral mass, ", m[M], " (g)", sep = "")),
            y = expression(
              paste("Cumulative SOC stock, ", S[C], " (Mg C/ha)", sep = "")),
            color = "",
            title = "CEDF")
        
        plot_esm_cedf_depth <-
          ggplot() +
          geom_ribbon(
            data = data_curve %>% 
              dplyr::filter(`site` %in% c(treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("CEDF")),
            aes(x = `cum_soil_mineral_mass_ave_g`,
                ymin = (`cum_soil_depth_ave_cm` - `cum_soil_depth_err_cm`),
                ymax = (`cum_soil_depth_ave_cm` + `cum_soil_depth_err_cm`)),
            color = "#1b7837", linetype = 0, alpha = 0.2) +
          geom_line(
            data = data_curve %>% 
              dplyr::filter(`site` %in% c(treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("CEDF")),
            aes(x = `cum_soil_mineral_mass_ave_g`,
                y = `cum_soil_depth_ave_cm`),
            color = "#1b7837", linetype = 1) + 
          geom_errorbar(
            data = data_final %>% 
              dplyr::filter(`site` %in% c(control_site, treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("NA", "CEDF")),
            aes(x = `cum_soil_mineral_mass_ave_g`,
                ymin = (`cum_soil_depth_ave_cm` - `cum_soil_depth_err_cm`),
                ymax = (`cum_soil_depth_ave_cm` + `cum_soil_depth_err_cm`),
                color = paste(`site`, ": ", `method`, sep = "")),
            width = 20) +
          geom_errorbarh(
            data = data_final %>% 
              dplyr::filter(`site` %in% c(control_site, treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("CEDF")),
            aes(y = `cum_soil_depth_ave_cm`,
                xmin = (`cum_soil_mineral_mass_ave_g` - 
                          `cum_soil_mineral_mass_err_g`),
                xmax = (`cum_soil_mineral_mass_ave_g` + 
                          `cum_soil_mineral_mass_err_g`),
                color = paste(`site`, ": ", `method`, sep = "")), 
            height = 5) +
          geom_point(
            data = data_final %>% 
              dplyr::filter(`site` %in% c(control_site, treatment_sites[m]),
                            `soil_content` %in% c(soil_content),
                            `k` %in% c(k),
                            `regression_function` %in% c("NA", "CEDF")),
            aes(x = `cum_soil_mineral_mass_ave_g`,
                y = `cum_soil_depth_ave_cm`,
                color = paste(`site`, ": ", `method`, sep = "")), 
            size = 2, shape = 15) +
          theme_bw() +
          theme(legend.position = "bottom") +
          scale_x_continuous(
            limits = c(0, max(data_final$cum_soil_mineral_mass_ave_g) * 1.05)) +
          scale_y_continuous(
            limits = c(0, 
                       max(data_final$cum_soil_depth_ave_cm + 
                             data_final$cum_soil_depth_err_cm) * 1.05)) +
          scale_color_manual(
            values = c("#af8dc3", "#7fbf7b", "#762a83", "white"),
            labels = c(paste(control_site, ":", "FD"),
                       paste(treatment_sites[m], ":", "FD"),
                       paste(treatment_sites[m], ":", "ESM"), "")) +
          labs(x = expression(
            paste("Cumulative soil mineral mass, ", m[M], " (g)", sep = "")),
            y = expression(
              paste("Cumulative soil depth, ", z, " (cm)", sep = "")),
            color = "",
            title = "CEDF")
      }
      
      # combine individual ESM correction plots based on whether LIN, MCS
      # and/or CEDF have been used for ESM correction
      if("LIN" %in% data_final$regression_function &&
         "MCS" %in% data_final$regression_function &&
         "CEDF" %in% data_final$regression_function){
        
        # grid plot of SOC stocks and soil depth
        plot_esm_all_stock <-
          plot_esm_lin_stock +
          plot_esm_mcs_stock +
          plot_esm_cedf_stock +
          plot_layout(axes = "collect",
                      axis_titles = "collect",
                      guides = "collect") &
          theme(legend.position = "bottom")
        plot_esm_all_depth <-
          plot_esm_lin_depth +
          plot_esm_mcs_depth +
          plot_esm_cedf_depth +
          plot_layout(axes = "collect",
                      axis_titles = "collect",
                      guides = "collect") &
          theme(legend.position = "bottom")
        
      } else if("LIN" %in% data_final$regression_function &&
                "MCS" %in% data_final$regression_function &&
                !("CEDF" %in% data_final$regression_function)){
        
        # grid plot of SOC stocks and soil depth
        plot_esm_all_stock <-
          plot_esm_lin_stock +
          plot_esm_mcs_stock +
          plot_layout(axes = "collect",
                      axis_titles = "collect",
                      guides = "collect") &
          theme(legend.position = "bottom")
        plot_esm_all_depth <-
          plot_esm_lin_depth +
          plot_esm_mcs_depth +
          plot_layout(axes = "collect",
                      axis_titles = "collect",
                      guides = "collect") &
          theme(legend.position = "bottom")
        
      } else if("LIN" %in% data_final$regression_function &&
                !("MCS" %in% data_final$regression_function) &&
                "CEDF" %in% data_final$regression_function){
        
        # grid plot of SOC stocks and soil depth
        plot_esm_all_stock <-
          plot_esm_lin_stock +
          plot_esm_cedf_stock +
          plot_layout(axes = "collect",
                      axis_titles = "collect",
                      guides = "collect") &
          theme(legend.position = "bottom")
        plot_esm_all_depth <-
          plot_esm_lin_depth +
          plot_esm_cedf_depth +
          plot_layout(axes = "collect",
                      axis_titles = "collect",
                      guides = "collect") &
          theme(legend.position = "bottom")
        
      } else if(!("LIN" %in% data_final$regression_function) &&
                "MCS" %in% data_final$regression_function &&
                "CEDF" %in% data_final$regression_function){
        
        # grid plot of SOC stocks and soil depth
        plot_esm_all_stock <-
          plot_esm_mcs_stock +
          plot_esm_cedf_stock +
          plot_layout(axes = "collect",
                      axis_titles = "collect",
                      guides = "collect") &
          theme(legend.position = "bottom")
        plot_esm_all_depth <-
          plot_esm_mcs_depth +
          plot_esm_cedf_depth +
          plot_layout(axes = "collect",
                      axis_titles = "collect",
                      guides = "collect") &
          theme(legend.position = "bottom")
        
      } else if("LIN" %in% data_final$regression_function &&
                !("MCS" %in% data_final$regression_function) &&
                !("CEDF" %in% data_final$regression_function)){
        
        # grid plot of SOC stocks and soil depth
        plot_esm_all_stock <- plot_esm_lin_stock
        plot_esm_all_depth <- plot_esm_lin_depth
        
      } else if(!("LIN" %in% data_final$regression_function) &&
                "MCS" %in% data_final$regression_function &&
                !("CEDF" %in% data_final$regression_function)){
        
        # grid plot of SOC stocks and soil depth
        plot_esm_all_stock <- plot_esm_mcs_stock
        plot_esm_all_depth <- plot_esm_mcs_depth
        
      } else if(!("LIN" %in% data_final$regression_function) &&
                !("MCS" %in% data_final$regression_function) &&
                "CEDF" %in% data_final$regression_function){
        
        # grid plot of SOC stocks and soil depth
        plot_esm_all_stock <- plot_esm_cedf_stock
        plot_esm_all_depth <- plot_esm_cedf_depth
        
      } else{
        
        # print statements that no ESM correction was found
        plot_esm_all_stock <- "No ESM correction found"
        plot_esm_all_depth <- "No ESM correction found"
      }
      
      # plot ESM correction curves
      print(plot_esm_all_stock)
      print(plot_esm_all_depth)
      
    }
    
    # plot cumulative FD-based and ESM-corrected SOC stocks by sample increment
    for(m in 1:length(unique(data_final$increment))){
      
      plot_stock_increment <-
        ggplot(data = data_final %>%
                 dplyr::filter(!(`site` %in% control_site),
                               `increment` %in% unique(data_final$increment)[m],
                               `soil_content` %in% soil_content,
                               `k` %in% k) %>%
                 dplyr::mutate(`regression_function` =
                                 replace(`regression_function`,
                                         `regression_function` %in% "NA",
                                         "")) %>%
                 dplyr::mutate(label = 
                                 paste(`method`, " ", `regression_function`,
                                       sep = "")),
               aes(x = `label`, y = `cum_soc_stock_ave_Mgha`, fill = `site`)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        geom_errorbar(aes(x = `label`,
                          ymin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
                          ymax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`)),
                      position = position_dodge()) +
        theme_bw() +
        labs(x = "Calculation",
             y = "Cumulative SOC stock (Mg C/ha)",
             fill = "Sites",
             title = paste("Sample increment = ", unique(data_final$increment)[m], 
                           sep = ""))
      
      print(plot_stock_increment)
    }
}

# define relative effect of uncertainties
uncertainty_analyser <-
  function(data_results, control_site, k, soil_content){
    
    # create empty data frame
    data_uncertainty_effect <- data.frame()
    
    # separate data 
    data_error_propagation <- data_results
    
    # get all treatment field labels
    treatment_sites <- 
      unique(data_error_propagation[
        !(data_error_propagation$site %in% c(control_site)),]$site)
    
    # loop over all treatment sites
    for(m in 1:length(treatment_sites)){
      
      # filter for all unique entries of m'th treatment site
      data_uncertainty <- data_error_propagation %>%
        dplyr::filter(`site` %in% c(treatment_sites[m]),
                      `k` %in% c(k),
                      `soil_content` %in% c(soil_content)) %>%
        distinct(`method`, `regression_function`, `uncertainty_off`,
                 `increment`, .keep_all = TRUE)
      
      # loop over each regression function, sample increment and FD/ESM
      # method
      for(a in 1:length(unique(data_uncertainty$regression_function))){
        for(b in 1:length(unique(data_uncertainty[
          data_uncertainty$regression_function %in%
          unique(data_uncertainty$regression_function)[a],]$increment))){
          for(d in 1:length(unique(data_uncertainty[
            data_uncertainty$increment %in% 
            unique(data_uncertainty[
              data_uncertainty$regression_function %in%
              unique(data_uncertainty$regression_function)[a],]$increment)[b],
          ]$method))){
            
            # get SOC stock errors of all uncertainty included and each
            # uncertainty term excluded as well as 
            data_uncertainty_calculation <-
              data.frame(
                error = as.numeric(
                  (data_uncertainty %>%
                     dplyr::filter(
                       `regression_function` %in% 
                         unique(data_uncertainty$regression_function)[a],
                       `increment` %in% 
                         unique(data_uncertainty[
                           data_uncertainty$regression_function %in%
                             unique(data_uncertainty$regression_function)[a],
                           ]$increment)[b],
                       `method` %in% 
                         unique(data_uncertainty[
                           data_uncertainty$increment %in% 
                             unique(data_uncertainty[
                               data_uncertainty$regression_function %in%
                                 unique(data_uncertainty$regression_function
                                        )[a],]$increment)[b],
                         ]$method)[d]
                     ))$cum_soc_stock_err_Mgha))
            
            # calculate MAE and attach to data frame
            data_uncertainty_calculation_2 <- cbind(
              data_uncertainty_calculation,
              data.frame(
                as.numeric(
                  (data_uncertainty %>%
                     dplyr::filter(
                       `regression_function` %in% 
                         unique(data_uncertainty$regression_function)[a],
                       `increment` %in% 
                         unique(data_uncertainty[
                           data_uncertainty$regression_function %in%
                             unique(data_uncertainty$regression_function)[a],
                         ]$increment)[b],
                       `method` %in% 
                         unique(data_uncertainty[
                           data_uncertainty$increment %in% 
                             unique(data_uncertainty[
                               data_uncertainty$regression_function %in%
                                 unique(data_uncertainty$regression_function
                                 )[a],]$increment)[b],
                         ]$method)[d],
                       `uncertainty_off` %in% c("none")
                     ))$cum_soc_stock_err_Mgha)
                - data_uncertainty_calculation["error"]) %>%
                dplyr::rename(mae = `error`))
            
            # calculate relative effect of uncertainties on error of SOC
            # stocks 
            data_uncertainty_calculation_3 <-
              cbind(data_uncertainty_calculation_2,
                    data.frame(effect = 
                      data_uncertainty_calculation_2$mae * 100 / 
                      (data_uncertainty_calculation_2 %>%
                         dplyr::summarise_all(sum))$mae))
            
            # bind results of relative effect to data frame
            data_uncertainty_effect <-
              rbind(data_uncertainty_effect,
                    cbind(
                      data_uncertainty %>%
                        dplyr::filter(
                          `regression_function` %in% 
                            unique(data_uncertainty$regression_function)[a],
                          `increment` %in% 
                            unique(data_uncertainty[
                              data_uncertainty$regression_function %in%
                                unique(data_uncertainty$regression_function)[a],
                            ]$increment)[b],
                          `method` %in% 
                            unique(data_uncertainty[
                              data_uncertainty$increment %in% 
                                unique(data_uncertainty[
                                  data_uncertainty$regression_function %in%
                                    unique(data_uncertainty$regression_function
                                    )[a],]$increment)[b],
                            ]$method)[d]
                        ),
                    data_uncertainty_calculation_3) %>%
                      dplyr::filter(!(`uncertainty_off` %in% c("none"))))
          }
        }
      }
      
      # check which soil content method is used for uncertainty calculation and
      # assign appropriate plotting labels and configurations
      if(soil_content %in% "soc"){
        soil_content_include <- c("soc_content")
        soil_content_exclude <- c("som_content")
      } else if(soil_content %in% "som"){
        soil_content_include <- c("som_content")
        soil_content_exclude <- c("soc_content")
      } else if(soil_content %in% "both"){
        soil_content_include <- c("soc_content", "som_content")
        soil_content_exclude <- c("") 
      }
      
      # plot relative effect of uncertainties from FD
      if("FD" %in% data_error_propagation$method &&
         treatment_sites[m] %in% data_error_propagation$site){
        
        plot_uncertainty_fd <-
          ggplot(
            data = data_uncertainty_effect %>%
              dplyr::filter(`site` %in% c(treatment_sites[m]),
                            `method` %in% c("FD"),
                            !(`uncertainty_off` %in% soil_content_exclude),
                            `regression_function` %in% c("NA")),
            aes(y = as.character(`increment`),
                x = `effect`,
                fill = factor(`uncertainty_off`, 
                              levels = c("a", "b", "k", "sample_thickness",
                                         "fine_soil_mass", soil_content_include, 
                                         "sample_volume")),
                group = factor(`uncertainty_off`,
                               levels = c(c("a", "b", "k", "sample_thickness",
                                            "fine_soil_mass", soil_content_include, 
                                            "sample_volume"))))) +
          geom_bar(position = "stack", stat = "identity", color = "black") +
          geom_vline(xintercept = c(0, 100), color = "black", linetype = 2) +
          theme_bw() +
          theme(legend.position = "bottom") +
          scale_fill_manual(
            values = if(soil_content %in% "soc"){
              c("#d53e4f", "#fc8d59", "#fee08b", "#ffffbf", "#e6f598", "#99d594",
                "#3288bd")} else if(soil_content %in% "som"){
                  c("#d53e4f", "#fc8d59", "#fee08b", "#ffffbf", "#e6f598", "#99d594",
                    "#3288bd")} else if(soil_content %in% "both"){
                      c("#d73027", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4",
                        "#66c2a5", "#3288bd")}, name = "",
            labels = if(soil_content %in% "soc"){
              c("a", "b", "k", expression(paste(Delta, ~z, sep = "")),
                expression(paste(m[F])), expression(SOC), expression(paste(V)))
            } else if(soil_content %in% "som"){
              c("a", "b", "k", expression(paste(Delta, ~z, sep = "")),
                expression(paste(m[F])), expression(SOM), expression(paste(V)))
            } else if(soil_content %in% "both"){
              c("a", "b", "k", expression(paste(Delta, ~z, sep = "")),
                expression(paste(m[F])), expression(SOC), expression(SOM), 
                expression(paste(V)))}, drop = FALSE) +
          labs(y = "Cumulative sample increments",
               x = expression(
                 paste("Relative effect on cumulative SOC stock error (%)",
                       sep = "")),
               fill = "",
               title = paste("FD     ", treatment_sites[m], sep = ""))
      }
      
      # plot relative effect of uncertainties from LIN
      if("LIN" %in% data_error_propagation$regression_function){
        
        plot_uncertainty_esm_lin <-
          ggplot(
            data = data_uncertainty_effect %>%
              dplyr::filter(`site` %in% c(treatment_sites[m]),
                            !(`uncertainty_off` %in% soil_content_exclude),
                            `regression_function` %in% c("LIN"),
                            `method` %in% c("ESM")),
            aes(y = as.character(`increment`),
                x = `effect`,
                fill = factor(`uncertainty_off`, 
                              levels = c("a", "b", "k", "sample_thickness",
                                         "fine_soil_mass", soil_content_include, 
                                         "sample_volume")),
                group = factor(`uncertainty_off`,
                               levels = c(c("a", "b", "k", "sample_thickness",
                                            "fine_soil_mass", soil_content_include, 
                                            "sample_volume"))))) +
          geom_bar(position = "stack", stat = "identity", color = "black") +
          geom_vline(xintercept = c(0, 100), color = "black", linetype = 2) +
          theme_bw() +
          theme(legend.position = "bottom") +
          scale_fill_manual(
            values = if(soil_content %in% "soc"){
              c("#d53e4f", "#fc8d59", "#fee08b", "#ffffbf", "#e6f598", "#99d594",
                "#3288bd")} else if(soil_content %in% "som"){
                  c("#d53e4f", "#fc8d59", "#fee08b", "#ffffbf", "#e6f598", "#99d594",
                    "#3288bd")} else if(soil_content %in% "both"){
                      c("#d73027", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4",
                        "#66c2a5", "#3288bd")}, name = "",
            labels = if(soil_content %in% "soc"){
              c("a", "b", "k", expression(paste(Delta, ~z, sep = "")),
                expression(paste(m[F])), expression(SOC), expression(paste(V)))
            } else if(soil_content %in% "som"){
              c("a", "b", "k", expression(paste(Delta, ~z, sep = "")),
                expression(paste(m[F])), expression(SOM), expression(paste(V)))
            } else if(soil_content %in% "both"){
              c("a", "b", "k", expression(paste(Delta, ~z, sep = "")),
                expression(paste(m[F])), expression(SOC), expression(SOM), 
                expression(paste(V)))}, drop = FALSE) +
          labs(y = "Cumulative sample increments",
               x = expression(
                 paste("Relative effect on cumulative SOC stock error (%)",
                       sep = "")),
               fill = "",
               title = "ESM: LIN")
      }
      
      # plot relative effect of uncertainties from MCS
      if("MCS" %in% data_error_propagation$regression_function){
        
        plot_uncertainty_esm_mcs <-
          ggplot(
            data = data_uncertainty_effect %>%
              dplyr::filter(`site` %in% c(treatment_sites[m]),
                            !(`uncertainty_off` %in% soil_content_exclude),
                            `regression_function` %in% c("MCS"),
                            `method` %in% c("ESM")),
            aes(y = as.character(`increment`),
                x = `effect`,
                fill = factor(`uncertainty_off`, 
                              levels = c("a", "b", "k", "sample_thickness",
                                         "fine_soil_mass", soil_content_include, 
                                         "sample_volume")),
                group = factor(`uncertainty_off`,
                               levels = c(c("a", "b", "k", "sample_thickness",
                                            "fine_soil_mass", soil_content_include, 
                                            "sample_volume"))))) +
          geom_bar(position = "stack", stat = "identity", color = "black") +
          geom_vline(xintercept = c(0, 100), color = "black", linetype = 2) +
          theme_bw() +
          theme(legend.position = "bottom") +
          scale_fill_manual(
            values = if(soil_content %in% "soc"){
              c("#d53e4f", "#fc8d59", "#fee08b", "#ffffbf", "#e6f598", "#99d594",
                "#3288bd")} else if(soil_content %in% "som"){
                  c("#d53e4f", "#fc8d59", "#fee08b", "#ffffbf", "#e6f598", "#99d594",
                    "#3288bd")} else if(soil_content %in% "both"){
                      c("#d73027", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4",
                        "#66c2a5", "#3288bd")}, name = "",
            labels = if(soil_content %in% "soc"){
              c("a", "b", "k", expression(paste(Delta, ~z, sep = "")),
                expression(paste(m[F])), expression(SOC), expression(paste(V)))
            } else if(soil_content %in% "som"){
              c("a", "b", "k", expression(paste(Delta, ~z, sep = "")),
                expression(paste(m[F])), expression(SOM), expression(paste(V)))
            } else if(soil_content %in% "both"){
              c("a", "b", "k", expression(paste(Delta, ~z, sep = "")),
                expression(paste(m[F])), expression(SOC), expression(SOM), 
                expression(paste(V)))}, drop = FALSE) +
          labs(y = "Cumulative sample increments",
               x = expression(
                 paste("Relative effect on cumulative SOC stock error (%)",
                       sep = "")),
               fill = "",
               title = "ESM: MCS")
      }
      
      # plot relative effect of uncertainties from CEDF
      if("CEDF" %in% data_error_propagation$regression_function){
        
        plot_uncertainty_esm_cedf <-
          ggplot(
            data = data_uncertainty_effect %>%
              dplyr::filter(`site` %in% c(treatment_sites[m]),
                            !(`uncertainty_off` %in% soil_content_exclude),
                            `regression_function` %in% c("CEDF"),
                            `method` %in% c("ESM")),
            aes(y = as.character(`increment`),
                x = `effect`,
                fill = factor(`uncertainty_off`, 
                              levels = c("a", "b", "k", "sample_thickness",
                                         "fine_soil_mass", soil_content_include, 
                                         "sample_volume")),
                group = factor(`uncertainty_off`,
                               levels = c(c("a", "b", "k", "sample_thickness",
                                            "fine_soil_mass", soil_content_include, 
                                            "sample_volume"))))) +
          geom_bar(position = "stack", stat = "identity", color = "black") +
          geom_vline(xintercept = c(0, 100), color = "black", linetype = 2) +
          theme_bw() +
          theme(legend.position = "bottom") +
          scale_fill_manual(
            values = if(soil_content %in% "soc"){
              c("#d53e4f", "#fc8d59", "#fee08b", "#ffffbf", "#e6f598", "#99d594",
                "#3288bd")} else if(soil_content %in% "som"){
                  c("#d53e4f", "#fc8d59", "#fee08b", "#ffffbf", "#e6f598", "#99d594",
                    "#3288bd")} else if(soil_content %in% "both"){
                      c("#d73027", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4",
                        "#66c2a5", "#3288bd")}, name = "",
            labels = if(soil_content %in% "soc"){
              c("a", "b", "k", expression(paste(Delta, ~z, sep = "")),
                expression(paste(m[F])), expression(SOC), expression(paste(V)))
            } else if(soil_content %in% "som"){
              c("a", "b", "k", expression(paste(Delta, ~z, sep = "")),
                expression(paste(m[F])), expression(SOM), expression(paste(V)))
            } else if(soil_content %in% "both"){
              c("a", "b", "k", expression(paste(Delta, ~z, sep = "")),
                expression(paste(m[F])), expression(SOC), expression(SOM), 
                expression(paste(V)))}, drop = FALSE) +
          labs(y = "Cumulative sample increments",
               x = expression(
                 paste("Relative effect on cumulative SOC stock error (%)",
                       sep = "")),
               fill = "",
               title = "ESM: CEDF")
      }
      
      # combine individual ESM correction plots based on whether LIN, MCS
      # and/or CEDF have been used for ESM correction
      if("LIN" %in% data_error_propagation$regression_function &&
         "MCS" %in% data_error_propagation$regression_function &&
         "CEDF" %in% data_error_propagation$regression_function){
        
        # grid plot of uncertainty effect by method
        plot_uncertainty <-
          plot_uncertainty_fd +
          plot_uncertainty_esm_lin +
          plot_uncertainty_esm_mcs +
          plot_uncertainty_esm_cedf +
          plot_layout(axes = "collect",
                      axis_titles = "collect",
                      guides = "collect") &
          theme(legend.position = "bottom")
        
      } else if("LIN" %in% data_error_propagation$regression_function &&
                "MCS" %in% data_error_propagation$regression_function &&
                !("CEDF" %in% data_error_propagation$regression_function)){
        
        # grid plot of uncertainty effect by method
        plot_uncertainty <-
          plot_uncertainty_fd +
          plot_uncertainty_esm_lin +
          plot_uncertainty_esm_mcs +
          plot_layout(axes = "collect",
                      axis_titles = "collect",
                      guides = "collect") &
          theme(legend.position = "bottom")
        
      } else if("LIN" %in% data_error_propagation$regression_function &&
                !("MCS" %in% data_error_propagation$regression_function) &&
                "CEDF" %in% data_error_propagation$regression_function){
        
        # grid plot of uncertainty effect by method
        plot_uncertainty <-
          plot_uncertainty_fd +
          plot_uncertainty_esm_lin +
          plot_uncertainty_esm_cedf +
          plot_layout(axes = "collect",
                      axis_titles = "collect",
                      guides = "collect") &
          theme(legend.position = "bottom")
        
      } else if(!("LIN" %in% data_error_propagation$regression_function) &&
                "MCS" %in% data_error_propagation$regression_function &&
                "CEDF" %in% data_error_propagation$regression_function){
        
        # grid plot of uncertainty effect by method
        plot_uncertainty <-
          plot_uncertainty_fd +
          plot_uncertainty_esm_mcs +
          plot_uncertainty_esm_cedf +
          plot_layout(axes = "collect",
                      axis_titles = "collect",
                      guides = "collect") &
          theme(legend.position = "bottom")
        
      } else if("LIN" %in% data_error_propagation$regression_function &&
                !("MCS" %in% data_error_propagation$regression_function) &&
                !("CEDF" %in% data_error_propagation$regression_function)){
        
        # grid plot of uncertainty effect by method
        plot_uncertainty <-
          plot_uncertainty_fd +
          plot_uncertainty_esm_lin +
          plot_layout(axes = "collect",
                      axis_titles = "collect",
                      guides = "collect") &
          theme(legend.position = "bottom")
        
      } else if(!("LIN" %in% data_error_propagation$regression_function) &&
                "MCS" %in% data_error_propagation$regression_function &&
                !("CEDF" %in% data_error_propagation$regression_function)){
        
        # grid plot of uncertainty effect by method
        plot_uncertainty <-
          plot_uncertainty_fd +
          plot_uncertainty_esm_mcs +
          plot_layout(axes = "collect",
                      axis_titles = "collect",
                      guides = "collect") &
          theme(legend.position = "bottom")
        
      } else if(!("LIN" %in% data_error_propagation$regression_function) &&
                !("MCS" %in% data_error_propagation$regression_function) &&
                "CEDF" %in% data_error_propagation$regression_function){
        
        # grid plot of uncertainty effect by method
        plot_uncertainty <-
          plot_uncertainty_fd +
          plot_uncertainty_esm_cedf +
          plot_layout(axes = "collect",
                      axis_titles = "collect",
                      guides = "collect") &
          theme(legend.position = "bottom")
        
      } else{
        
        # plot of uncertainty effect by FD method
        plot_uncertainty <-
          plot_uncertainty_fd
      }
      
      # plot relative effect of uncertainties on SOC stock error
      print(plot_uncertainty)
      
    }
  }

### Accessory functions ########################################################

# define lower and upper indices finder for fitting LIN
find_bound_indices <- function(vector, value){
  
  if(min(vector[vector >= value]) == Inf){
    upper <- match(max(vector[vector <= value]), vector)
    lower <- upper - 1
  } else if(max(vector[vector <= value]) == Inf){
    lower <- match(min(vector[vector >= value]), vector)
    upper <- lower + 1
  } else{
    lower <- match(max(vector[vector <= value]), vector)
    upper <- match(min(vector[vector >= value]), vector)
  }
  return(c(lower, upper))
}

# define function to rename column names to that used in algorithm
data_renamer <- function(data, original_label_index){
  
  # check whether SOM content is provided
  # if not, relabel without SOM content
  if(is.na(original_label_index[8]) == TRUE){
    
    # define column renames for algorithm
    column_rename <- 
      c("study", "site", "increment", "core", "sample_thickness_cm", 
        "fine_soil_mass_g", "soc_content_perc", 
        "sample_volume_cm3")
    
    original_label_index <- original_label_index[!is.na(original_label_index)]
    
  } else{
  
    # define column renames for algorithm
    column_rename <- 
      c("study", "site", "increment", "core", "sample_thickness_cm", 
        "fine_soil_mass_g", "soc_content_perc", "som_content_perc",
        "sample_volume_cm3")
    
  }
  
  # rename key input data columns according to user input
  data_rename <- data %>%
    dplyr::rename_with(~column_rename,
                       colnames(data[original_label_index]))
  
  # print which column labels were renamed as check
  for(i in 1:length(column_rename)){
    message(paste(column_rename[i], " <- ", 
                colnames(data_import)[data_column_index[i]], sep = ""))
  }
  
  # return renamed data
  return(data_rename)
}
