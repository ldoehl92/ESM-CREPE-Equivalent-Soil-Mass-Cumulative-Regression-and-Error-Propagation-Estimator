### ESM-CREPE package - Doehl et al. (2026) ####################################

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

# File name:    manuscript_doehl2026.R
# File purpose: Modified version of ESM-CREPE to reproduce the analysis and
#               results published in the following manuscript:
#               Doehl, L.N.K., Biffi, S., Black, H.I.J., Chapman, P.J. and
#               Ziv, G. (2026). Improved quantification of changes in soil 
#               organic carbon stocks using equivalent soil mass correction
#               (submitted for publication).

### ESM correction estimator - modified for Doehl et al. (2026) #######

# define main function to ESM correction modified to include the soil sampling 
# protocols used in Doehl et al. (2026)
esm_correction_main_doehl2026 <-
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
                esm_calculator_doehl2026(
                  data = data,
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
                 `k`, `regression_function`, `ssp`, `cum_soc_stock_ave_Mgha`, 
                 `cum_soc_stock_err_Mgha`, `cum_soc_stock_ci_Mgha`, 
                 .keep_all = TRUE)
      message("Finished calculations.")
      
      # plot ESM correction curve and uncertainty effects for each treatment site
      message("Plotting visual of ESM correction...")
      visual_analyser_doehl2026(
        data_results = list(data_final, data_curve), 
        control_site = control_site, k = k, soil_content = soil_content)
      message("Plotted ESM-corrected SOC stock. and soil depths.")
      message("")
      message("Plotting uncertainty analysis...")
      uncertainty_analyser_doehl2026(
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

# define ESM calculator modified to include the soil sampling protocols
# used in Doehl et al. (2026)
esm_calculator_doehl2026 <-
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
    
    # define soil sampling protocols
    ssp <- c("SSP1", "SSP2a", "SSP2b", "SSP3")
    
    # for loop over all uncertainties included and then excluded
    for(i in 1:length(uncertainty_switch)){
      
      # calculate SOC stock and mineral mass
      data_cum_calculation <-
        stock_and_mass_calculator(
          data, control_site, soil_content, k, guess_k, 
          uncertainty_switch, i)
      
      #
      for(j in 1:length(ssp)){
       
        if(ssp[j] == "SSP1"){
          
          # manipulate calculations to match Verra protocol
          sample_increments_mod <- unique(data_cum_calculation$increment)[c(3,5)]
          
        } else if(ssp[j] == "SSP2a"){
          
          # manipulate calculations to match modified protocol
          sample_increments_mod <- unique(data_cum_calculation$increment)[c(3:5)]
          
        } else if(ssp[j] == "SSP2b"){
          
          # manipulate calculations to match modified protocol
          sample_increments_mod <- unique(data_cum_calculation$increment)[c(1,3,5)]
          
        } else if(ssp[j] == "SSP3") {
          
          # use all sampling horizons
          sample_increments_mod <- unique(data_cum_calculation$increment)
          
        }
        
        # manipulate calculations to match selected protocol
        data_cum_calculation_mod <- data_cum_calculation %>%
          dplyr::filter(`increment` %in% sample_increments_mod)
        
        # create mineral mass range for curve to plot
        curve_soil_mineral_mass_g <-
          seq(0, max(data_cum_calculation_mod$cum_soil_mineral_mass_ave2_g) + 50, 10)
        
        # separate data of control and treatment fields
        data_treatment <-
          data_cum_calculation_mod %>%
          dplyr::filter(!(`site` %in% control_site))
        data_control <-
          data_cum_calculation_mod %>%
          dplyr::filter(`site` %in% control_site)
        
        # run ESM correction with all uncertainties included
        if(uncertainty_switch[i] == "none"){
          
          # collect soil mineral masses, soil depth and FD-based SOC stock 
          # of all control and treatment sites
          data_final = 
            rbind(data_final,
                  data.frame(
                    study = unique(data$study),
                    site = data_cum_calculation_mod$site,
                    increment = data_cum_calculation_mod$increment,
                    method = "FD",
                    soil_content = soil_content, 
                    k = k, 
                    regression_function = "NA",
                    ssp = ssp[j],
                    sample_size = data_cum_calculation_mod$sample_size,
                    cum_soil_mineral_mass_ave_g = data_cum_calculation_mod$cum_soil_mineral_mass_ave2_g,
                    cum_soil_mineral_mass_err_g = data_cum_calculation_mod$cum_soil_mineral_mass_err2_g,
                    cum_soil_mineral_mass_ci_g = data_cum_calculation_mod$cum_soil_mineral_mass_ci2_g,
                    cum_soil_depth_ave_cm = data_cum_calculation_mod$cum_soil_depth_ave2_cm,
                    cum_soil_depth_err_cm = data_cum_calculation_mod$cum_soil_depth_err2_cm,
                    cum_soil_depth_ci_cm = data_cum_calculation_mod$cum_soil_depth_ci2_cm,
                    cum_soc_stock_ave_Mgha = data_cum_calculation_mod$cum_soc_stock_ave2_Mgha,
                    cum_soc_stock_err_Mgha = data_cum_calculation_mod$cum_soc_stock_err2_Mgha,
                    cum_soc_stock_ci_Mgha = data_cum_calculation_mod$cum_soc_stock_ci2_Mgha))
          
          # bind all calculations with no uncertainties switched off
          data_error_propagation <- rbind(
            data_error_propagation, 
            cbind(study = unique(data$study),
                  method = "FD",
                  soil_content = soil_content, 
                  k = k, 
                  regression_function = "NA",
                  ssp = ssp[j],
                  uncertainty_off = "none",
                  data.frame(
                    site = data_cum_calculation_mod$site,
                    increment = data_cum_calculation_mod$increment,
                    sample_size = data_cum_calculation_mod$sample_size,
                    cum_soil_mineral_mass_ave_g = data_cum_calculation_mod$cum_soil_mineral_mass_ave2_g,
                    cum_soil_mineral_mass_err_g = data_cum_calculation_mod$cum_soil_mineral_mass_err2_g,
                    cum_soil_mineral_mass_ci_g = data_cum_calculation_mod$cum_soil_mineral_mass_ci2_g,
                    cum_soil_depth_ave_cm = data_cum_calculation_mod$cum_soil_depth_ave2_cm,
                    cum_soil_depth_err_cm = data_cum_calculation_mod$cum_soil_depth_err2_cm,
                    cum_soil_depth_ci_cm = data_cum_calculation_mod$cum_soil_depth_ci2_cm,
                    cum_soc_stock_ave_Mgha = data_cum_calculation_mod$cum_soc_stock_ave2_Mgha,
                    cum_soc_stock_err_Mgha = data_cum_calculation_mod$cum_soc_stock_err2_Mgha,
                    cum_soc_stock_ci_Mgha = data_cum_calculation_mod$cum_soc_stock_ci2_Mgha)))
          
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
                        ssp = ssp[j],
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
                        ssp = ssp[j],
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
                      ssp = ssp[j],
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
                        ssp = ssp[j],
                        uncertainty_off = uncertainty_switch[i],
                        data.frame(
                          site = data_cum_calculation_mod$site,
                          increment = data_cum_calculation_mod$increment,
                          sample_size = data_cum_calculation_mod$sample_size,
                          cum_soil_mineral_mass_ave_g = data_cum_calculation_mod$cum_soil_mineral_mass_ave2_g,
                          cum_soil_mineral_mass_err_g = data_cum_calculation_mod$cum_soil_mineral_mass_err2_g,
                          cum_soil_mineral_mass_ci_g = data_cum_calculation_mod$cum_soil_mineral_mass_ci2_g,
                          cum_soil_depth_ave_cm = data_cum_calculation_mod$cum_soil_depth_ave2_cm,
                          cum_soil_depth_err_cm = data_cum_calculation_mod$cum_soil_depth_err2_cm,
                          cum_soil_depth_ci_cm = data_cum_calculation_mod$cum_soil_depth_ci2_cm,
                          cum_soc_stock_ave_Mgha = data_cum_calculation_mod$cum_soc_stock_ave2_Mgha,
                          cum_soc_stock_err_Mgha = data_cum_calculation_mod$cum_soc_stock_err2_Mgha,
                          cum_soc_stock_ci_Mgha = data_cum_calculation_mod$cum_soc_stock_ci2_Mgha)))
                
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
                        ssp = ssp[j],
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
                        ssp = ssp[j],
                        uncertainty_off = uncertainty_switch[i],
                        data.frame(
                          site = data_cum_calculation_mod$site,
                          increment = data_cum_calculation_mod$increment,
                          sample_size = data_cum_calculation_mod$sample_size,
                          cum_soil_mineral_mass_ave_g = data_cum_calculation_mod$cum_soil_mineral_mass_ave2_g,
                          cum_soil_mineral_mass_err_g = data_cum_calculation_mod$cum_soil_mineral_mass_err2_g,
                          cum_soil_mineral_mass_ci_g = data_cum_calculation_mod$cum_soil_mineral_mass_ci2_g,
                          cum_soil_depth_ave_cm = data_cum_calculation_mod$cum_soil_depth_ave2_cm,
                          cum_soil_depth_err_cm = data_cum_calculation_mod$cum_soil_depth_err2_cm,
                          cum_soil_depth_ci_cm = data_cum_calculation_mod$cum_soil_depth_ci2_cm,
                          cum_soc_stock_ave_Mgha = data_cum_calculation_mod$cum_soc_stock_ave2_Mgha,
                          cum_soc_stock_err_Mgha = data_cum_calculation_mod$cum_soc_stock_err2_Mgha,
                          cum_soc_stock_ci_Mgha = data_cum_calculation_mod$cum_soc_stock_ci2_Mgha)))
                
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
                        ssp = ssp[j],
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
                        ssp = ssp[j],
                        uncertainty_off = uncertainty_switch[i],
                        data.frame(
                          site = data_cum_calculation_mod$site,
                          increment = data_cum_calculation_mod$increment,
                          sample_size = data_cum_calculation_mod$sample_size,
                          cum_soil_mineral_mass_ave_g = data_cum_calculation_mod$cum_soil_mineral_mass_ave2_g,
                          cum_soil_mineral_mass_err_g = data_cum_calculation_mod$cum_soil_mineral_mass_err2_g,
                          cum_soil_mineral_mass_ci_g = data_cum_calculation_mod$cum_soil_mineral_mass_ci2_g,
                          cum_soil_depth_ave_cm = data_cum_calculation_mod$cum_soil_depth_ave2_cm,
                          cum_soil_depth_err_cm = data_cum_calculation_mod$cum_soil_depth_err2_cm,
                          cum_soil_depth_ci_cm = data_cum_calculation_mod$cum_soil_depth_ci2_cm,
                          cum_soc_stock_ave_Mgha = data_cum_calculation_mod$cum_soc_stock_ave2_Mgha,
                          cum_soc_stock_err_Mgha = data_cum_calculation_mod$cum_soc_stock_err2_Mgha,
                          cum_soc_stock_ci_Mgha = data_cum_calculation_mod$cum_soc_stock_ci2_Mgha)))
                
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
                        ssp = ssp[j],
                        uncertainty_off = uncertainty_switch[i],
                        data_prefinal))
                
              } else{}
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
        } 
      }
    }
    
    # return data frames
    return(list(data_final, data_curve, data_calculations, data_error_propagation))
  }

# define visual analyser to plot ESM results to include the soil sampling 
# protocols used in Doehl et al. (2026)
visual_analyser_doehl2026 <-
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
    
    # plot FD-based and ESM-corrected SOC stocks from LIN and SSP1
    plot_esm_lin_stock_ssp1 <-
      ggplot() +
      geom_line(
        data = data_curve %>% 
          dplyr::filter(`site` %in% c(treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("LIN"),
                        `ssp` %in% c("SSP1")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            y = `cum_soc_stock_ave_Mgha`),
        color = "#1b7837", linetype = 1) + 
      geom_errorbar(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "LIN"),
                        `ssp` %in% c("SSP1")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            ymin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
            ymax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`),
            color = paste(`site`, ":", `method`, sep = "")),
        width = 20) +
      geom_errorbarh(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "LIN"),
                        `ssp` %in% c("SSP1")),
        aes(y = `cum_soc_stock_ave_Mgha`,
            xmin = (`cum_soil_mineral_mass_ave_g` - 
                      `cum_soil_mineral_mass_err_g`),
            xmax = (`cum_soil_mineral_mass_ave_g` + 
                      `cum_soil_mineral_mass_err_g`),
            color = paste(`site`, ":", `method`, sep = "")), 
        height = 5) +
      geom_point(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "LIN"),
                        `ssp` %in% c("SSP1")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            y = `cum_soc_stock_ave_Mgha`,
            color = paste(`site`, ":", `method`, sep = "")), 
        size = 1, shape = 15) +
      theme_bw() +
      theme(legend.position = "bottom") +
      scale_x_continuous(limits = 
                           c(0, max(data_final$cum_soil_mineral_mass_ave_g) * 1.05)) +
      scale_y_continuous(
        limits = c(0, max(data_final$cum_soc_stock_ave_Mgha + 
                            data_final$cum_soc_stock_err_Mgha) * 1.05)) +
      scale_color_manual(
        values = c("#af8dc3", "#762a83", "#7fbf7b"),
        labels = c(paste(control_site, ":", "FD"),
                   paste(treatment_sites, ":", "ESM"),
                   paste(treatment_sites, ":", "FD"))) +
      labs(x = expression(
        paste("Cumulative soil mineral mass, ", m[M], " (g)", sep = "")),
        y = expression(
          paste("Cumulative SOC stock, ", S[C], " (Mg C/ha)", sep = "")),
        color = "",
        title = "a   SSP1: LIN")
    
    # plot FD-based and ESM-corrected SOC stocks from LIN and SSP2a
    plot_esm_lin_stock_ssp2a <-
      ggplot() +
      geom_line(
        data = data_curve %>% 
          dplyr::filter(`site` %in% c(treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("LIN"),
                        `ssp` %in% c("SSP2a")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            y = `cum_soc_stock_ave_Mgha`),
        color = "#1b7837", linetype = 1) + 
      geom_errorbar(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "LIN"),
                        `ssp` %in% c("SSP2a")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            ymin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
            ymax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`),
            color = paste(`site`, ": ", `method`, sep = "")),
        width = 20) +
      geom_errorbarh(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "LIN"),
                        `ssp` %in% c("SSP2a")),
        aes(y = `cum_soc_stock_ave_Mgha`,
            xmin = (`cum_soil_mineral_mass_ave_g` - 
                      `cum_soil_mineral_mass_err_g`),
            xmax = (`cum_soil_mineral_mass_ave_g` + 
                      `cum_soil_mineral_mass_err_g`),
            color = paste(`site`, ": ", `method`, sep = "")), 
        height = 5) +
      geom_point(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "LIN"),
                        `ssp` %in% c("SSP2a")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            y = `cum_soc_stock_ave_Mgha`,
            color = paste(`site`, ": ", `method`, sep = "")), 
        size = 1, shape = 15) +
      theme_bw() +
      theme(legend.position = "bottom") +
      scale_x_continuous(limits = 
                           c(0, max(data_final$cum_soil_mineral_mass_ave_g) * 1.05)) +
      scale_y_continuous(
        limits = c(0, max(data_final$cum_soc_stock_ave_Mgha + 
                            data_final$cum_soc_stock_err_Mgha) * 1.05)) +
      scale_color_manual(
        values = c("#af8dc3", "#762a83", "#7fbf7b"),
        labels = c(paste(control_site, ":", "FD"),
                   paste(treatment_sites, ":", "ESM"),
                   paste(treatment_sites, ":", "FD"))) +
      labs(x = expression(
        paste("Cumulative soil mineral mass, ", m[M], " (g)", sep = "")),
        y = expression(
          paste("Cumulative SOC stock, ", S[C], " (Mg C/ha)", sep = "")),
        color = "",
        title = "b   SSP2a: LIN")
    
    # plot FD-based and ESM-corrected SOC stocks from LIN and SSP2b
    plot_esm_lin_stock_ssp2b <-
      ggplot() +
      geom_line(
        data = data_curve %>% 
          dplyr::filter(`site` %in% c(treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("LIN"),
                        `ssp` %in% c("SSP2b")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            y = `cum_soc_stock_ave_Mgha`),
        color = "#1b7837", linetype = 1) + 
      geom_errorbar(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "LIN"),
                        `ssp` %in% c("SSP2b")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            ymin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
            ymax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`),
            color = paste(`site`, ": ", `method`, sep = "")),
        width = 20) +
      geom_errorbarh(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "LIN"),
                        `ssp` %in% c("SSP2b")),
        aes(y = `cum_soc_stock_ave_Mgha`,
            xmin = (`cum_soil_mineral_mass_ave_g` - 
                      `cum_soil_mineral_mass_err_g`),
            xmax = (`cum_soil_mineral_mass_ave_g` + 
                      `cum_soil_mineral_mass_err_g`),
            color = paste(`site`, ": ", `method`, sep = "")), 
        height = 5) +
      geom_point(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "LIN"),
                        `ssp` %in% c("SSP2b")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            y = `cum_soc_stock_ave_Mgha`,
            color = paste(`site`, ": ", `method`, sep = "")), 
        size = 1, shape = 15) +
      theme_bw() +
      theme(legend.position = "bottom") +
      scale_x_continuous(limits = 
                           c(0, max(data_final$cum_soil_mineral_mass_ave_g) * 1.05)) +
      scale_y_continuous(
        limits = c(0, max(data_final$cum_soc_stock_ave_Mgha + 
                            data_final$cum_soc_stock_err_Mgha) * 1.05)) +
      scale_color_manual(
        values = c("#af8dc3", "#762a83", "#7fbf7b"),
        labels = c(paste(control_site, ":", "FD"),
                   paste(treatment_sites, ":", "ESM"),
                   paste(treatment_sites, ":", "FD"))) +
      labs(x = expression(
        paste("Cumulative soil mineral mass, ", m[M], " (g)", sep = "")),
        y = expression(
          paste("Cumulative SOC stock, ", S[C], " (Mg C/ha)", sep = "")),
        color = "",
        title = "c   SSP2b: LIN")
    
    # plot FD-based and ESM-corrected SOC stocks from LIN and SSP3
    plot_esm_lin_stock_ssp3 <-
      ggplot() +
      geom_line(
        data = data_curve %>% 
          dplyr::filter(`site` %in% c(treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("LIN"),
                        `ssp` %in% c("SSP3")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            y = `cum_soc_stock_ave_Mgha`),
        color = "#1b7837", linetype = 1) + 
      geom_errorbar(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "LIN"),
                        `ssp` %in% c("SSP3")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            ymin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
            ymax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`),
            color = paste(`site`, ": ", `method`, sep = "")),
        width = 20) +
      geom_errorbarh(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "LIN"),
                        `ssp` %in% c("SSP3")),
        aes(y = `cum_soc_stock_ave_Mgha`,
            xmin = (`cum_soil_mineral_mass_ave_g` - 
                      `cum_soil_mineral_mass_err_g`),
            xmax = (`cum_soil_mineral_mass_ave_g` + 
                      `cum_soil_mineral_mass_err_g`),
            color = paste(`site`, ": ", `method`, sep = "")), 
        height = 5) +
      geom_point(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "LIN"),
                        `ssp` %in% c("SSP3")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            y = `cum_soc_stock_ave_Mgha`,
            color = paste(`site`, ": ", `method`, sep = "")), 
        size = 1, shape = 15) +
      theme_bw() +
      theme(legend.position = "bottom") +
      scale_x_continuous(limits = 
                           c(0, max(data_final$cum_soil_mineral_mass_ave_g) * 1.05)) +
      scale_y_continuous(
        limits = c(0, max(data_final$cum_soc_stock_ave_Mgha + 
                            data_final$cum_soc_stock_err_Mgha) * 1.05)) +
      scale_color_manual(
        values = c("#af8dc3", "#762a83", "#7fbf7b"),
        labels = c(paste(control_site, ":", "FD"),
                   paste(treatment_sites, ":", "ESM"),
                   paste(treatment_sites, ":", "FD"))) +
      labs(x = expression(
        paste("Cumulative soil mineral mass, ", m[M], " (g)", sep = "")),
        y = expression(
          paste("Cumulative SOC stock, ", S[C], " (Mg C/ha)", sep = "")),
        color = "",
        title = "d   SSP3: LIN")
    
    # plot FD-based and ESM-corrected SOC stocks from MCS and SSP1
    plot_esm_mcs_stock_ssp1 <-
      ggplot() +
      geom_line(
        data = data_curve %>% 
          dplyr::filter(`site` %in% c(treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("MCS"),
                        `ssp` %in% c("SSP1")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            y = `cum_soc_stock_ave_Mgha`),
        color = "#1b7837", linetype = 1) + 
      geom_errorbar(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "MCS"),
                        `ssp` %in% c("SSP1")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            ymin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
            ymax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`),
            color = paste(`site`, ": ", `method`, sep = "")),
        width = 20) +
      geom_errorbarh(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "MCS"),
                        `ssp` %in% c("SSP1")),
        aes(y = `cum_soc_stock_ave_Mgha`,
            xmin = (`cum_soil_mineral_mass_ave_g` - 
                      `cum_soil_mineral_mass_err_g`),
            xmax = (`cum_soil_mineral_mass_ave_g` + 
                      `cum_soil_mineral_mass_err_g`),
            color = paste(`site`, ": ", `method`, sep = "")), 
        height = 5) +
      geom_point(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "MCS"),
                        `ssp` %in% c("SSP1")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            y = `cum_soc_stock_ave_Mgha`,
            color = paste(`site`, ": ", `method`, sep = "")), 
        size = 1, shape = 15) +
      theme_bw() +
      theme(legend.position = "bottom") +
      scale_x_continuous(
        limits = 
          c(0, max(data_final$cum_soil_mineral_mass_ave_g) * 1.05)) +
      scale_y_continuous(
        limits = c(0, max(data_final$cum_soc_stock_ave_Mgha + 
                            data_final$cum_soc_stock_err_Mgha) * 1.05)) +
      scale_color_manual(
        values = c("#af8dc3", "#762a83", "#7fbf7b"),
        labels = c(paste(control_site, ":", "FD"),
                   paste(treatment_sites, ":", "ESM"),
                   paste(treatment_sites, ":", "FD"))) +
      labs(x = expression(
        paste("Cumulative soil mineral mass, ", m[M], " (g)", sep = "")),
        y = expression(
          paste("Cumulative SOC stock, ", S[C], " (Mg C/ha)", sep = "")),
        color = "",
        title = "e   MCS SSP1")
    
    # plot FD-based and ESM-corrected SOC stocks from MCS and SSP2a
    plot_esm_mcs_stock_ssp2a <-
      ggplot() +
      geom_line(
        data = data_curve %>% 
          dplyr::filter(`site` %in% c(treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("MCS"),
                        `ssp` %in% c("SSP2a")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            y = `cum_soc_stock_ave_Mgha`),
        color = "#1b7837", linetype = 1) + 
      geom_errorbar(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "MCS"),
                        `ssp` %in% c("SSP2a")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            ymin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
            ymax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`),
            color = paste(`site`, ": ", `method`, sep = "")),
        width = 20) +
      geom_errorbarh(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "MCS"),
                        `ssp` %in% c("SSP2a")),
        aes(y = `cum_soc_stock_ave_Mgha`,
            xmin = (`cum_soil_mineral_mass_ave_g` - 
                      `cum_soil_mineral_mass_err_g`),
            xmax = (`cum_soil_mineral_mass_ave_g` + 
                      `cum_soil_mineral_mass_err_g`),
            color = paste(`site`, ": ", `method`, sep = "")), 
        height = 5) +
      geom_point(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "MCS"),
                        `ssp` %in% c("SSP2a")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            y = `cum_soc_stock_ave_Mgha`,
            color = paste(`site`, ": ", `method`, sep = "")), 
        size = 1, shape = 15) +
      theme_bw() +
      theme(legend.position = "bottom") +
      scale_x_continuous(
        limits = 
          c(0, max(data_final$cum_soil_mineral_mass_ave_g) * 1.05)) +
      scale_y_continuous(
        limits = c(0, max(data_final$cum_soc_stock_ave_Mgha + 
                            data_final$cum_soc_stock_err_Mgha) * 1.05)) +
      scale_color_manual(
        values = c("#af8dc3", "#762a83", "#7fbf7b"),
        labels = c(paste(control_site, ":", "FD"),
                   paste(treatment_sites, ":", "ESM"),
                   paste(treatment_sites, ":", "FD"))) +
      labs(x = expression(
        paste("Cumulative soil mineral mass, ", m[M], " (g)", sep = "")),
        y = expression(
          paste("Cumulative SOC stock, ", S[C], " (Mg C/ha)", sep = "")),
        color = "",
        title = "f   MCS SSP2a")
    
    # plot FD-based and ESM-corrected SOC stocks from MCS and SSP2b
    plot_esm_mcs_stock_ssp2b <-
      ggplot() +
      geom_line(
        data = data_curve %>% 
          dplyr::filter(`site` %in% c(treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("MCS"),
                        `ssp` %in% c("SSP2b")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            y = `cum_soc_stock_ave_Mgha`),
        color = "#1b7837", linetype = 1) + 
      geom_errorbar(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "MCS"),
                        `ssp` %in% c("SSP2b")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            ymin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
            ymax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`),
            color = paste(`site`, ": ", `method`, sep = "")),
        width = 20) +
      geom_errorbarh(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "MCS"),
                        `ssp` %in% c("SSP2b")),
        aes(y = `cum_soc_stock_ave_Mgha`,
            xmin = (`cum_soil_mineral_mass_ave_g` - 
                      `cum_soil_mineral_mass_err_g`),
            xmax = (`cum_soil_mineral_mass_ave_g` + 
                      `cum_soil_mineral_mass_err_g`),
            color = paste(`site`, ": ", `method`, sep = "")), 
        height = 5) +
      geom_point(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "MCS"),
                        `ssp` %in% c("SSP2b")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            y = `cum_soc_stock_ave_Mgha`,
            color = paste(`site`, ": ", `method`, sep = "")), 
        size = 1, shape = 15) +
      theme_bw() +
      theme(legend.position = "bottom") +
      scale_x_continuous(
        limits = 
          c(0, max(data_final$cum_soil_mineral_mass_ave_g) * 1.05)) +
      scale_y_continuous(
        limits = c(0, max(data_final$cum_soc_stock_ave_Mgha + 
                            data_final$cum_soc_stock_err_Mgha) * 1.05)) +
      scale_color_manual(
        values = c("#af8dc3", "#762a83", "#7fbf7b"),
        labels = c(paste(control_site, ":", "FD"),
                   paste(treatment_sites, ":", "ESM"),
                   paste(treatment_sites, ":", "FD"))) +
      labs(x = expression(
        paste("Cumulative soil mineral mass, ", m[M], " (g)", sep = "")),
        y = expression(
          paste("Cumulative SOC stock, ", S[C], " (Mg C/ha)", sep = "")),
        color = "",
        title = "g   MCS SSP2b")
    
    # plot FD-based and ESM-corrected SOC stocks from MCS and SSP3
    plot_esm_mcs_stock_ssp3 <-
      ggplot() +
      geom_line(
        data = data_curve %>% 
          dplyr::filter(`site` %in% c(treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("MCS"),
                        `ssp` %in% c("SSP3")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            y = `cum_soc_stock_ave_Mgha`),
        color = "#1b7837", linetype = 1) + 
      geom_errorbar(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "MCS"),
                        `ssp` %in% c("SSP3")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            ymin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
            ymax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`),
            color = paste(`site`, ": ", `method`, sep = "")),
        width = 20) +
      geom_errorbarh(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "MCS"),
                        `ssp` %in% c("SSP3")),
        aes(y = `cum_soc_stock_ave_Mgha`,
            xmin = (`cum_soil_mineral_mass_ave_g` - 
                      `cum_soil_mineral_mass_err_g`),
            xmax = (`cum_soil_mineral_mass_ave_g` + 
                      `cum_soil_mineral_mass_err_g`),
            color = paste(`site`, ": ", `method`, sep = "")), 
        height = 5) +
      geom_point(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "MCS"),
                        `ssp` %in% c("SSP3")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            y = `cum_soc_stock_ave_Mgha`,
            color = paste(`site`, ": ", `method`, sep = "")), 
        size = 1, shape = 15) +
      theme_bw() +
      theme(legend.position = "bottom") +
      scale_x_continuous(
        limits = 
          c(0, max(data_final$cum_soil_mineral_mass_ave_g) * 1.05)) +
      scale_y_continuous(
        limits = c(0, max(data_final$cum_soc_stock_ave_Mgha + 
                            data_final$cum_soc_stock_err_Mgha) * 1.05)) +
      scale_color_manual(
        values = c("#af8dc3", "#762a83", "#7fbf7b"),
        labels = c(paste(control_site, ":", "FD"),
                   paste(treatment_sites, ":", "ESM"),
                   paste(treatment_sites, ":", "FD"))) +
      labs(x = expression(
        paste("Cumulative soil mineral mass, ", m[M], " (g)", sep = "")),
        y = expression(
          paste("Cumulative SOC stock, ", S[C], " (Mg C/ha)", sep = "")),
        color = "",
        title = "h   MCS SSP3")
    
    # plot FD-based and ESM-corrected SOC stocks from CEDF and SSP2a
    plot_esm_cedf_stock_ssp2a <-
      ggplot() +
      geom_ribbon(
        data = data_curve %>% 
          dplyr::filter(`site` %in% c(treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("CEDF"),
                        `ssp` %in% c("SSP2a")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            ymin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
            ymax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`)),
        color = "#1b7837", linetype = 0, alpha = 0.2) +
      geom_line(
        data = data_curve %>% 
          dplyr::filter(`site` %in% c(treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("CEDF"),
                        `ssp` %in% c("SSP2a")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            y = `cum_soc_stock_ave_Mgha`),
        color = "#1b7837", linetype = 1) + 
      geom_errorbar(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "CEDF"),
                        `ssp` %in% c("SSP2a")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            ymin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
            ymax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`),
            color = paste(`site`, ": ", `method`, sep = "")),
        width = 20) +
      geom_errorbarh(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "CEDF"),
                        `ssp` %in% c("SSP2a")),
        aes(y = `cum_soc_stock_ave_Mgha`,
            xmin = (`cum_soil_mineral_mass_ave_g` - 
                      `cum_soil_mineral_mass_err_g`),
            xmax = (`cum_soil_mineral_mass_ave_g` + 
                      `cum_soil_mineral_mass_err_g`),
            color = paste(`site`, ": ", `method`, sep = "")), 
        height = 5) +
      geom_point(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "CEDF"),
                        `ssp` %in% c("SSP2a")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            y = `cum_soc_stock_ave_Mgha`,
            color = paste(`site`, ": ", `method`, sep = "")), 
        size = 1, shape = 15) +
      theme_bw() +
      theme(legend.position = "bottom") +
      scale_x_continuous(
        limits = 
          c(0, max(data_final$cum_soil_mineral_mass_ave_g) * 1.05)) +
      scale_y_continuous(
        limits = c(0, max(data_final$cum_soc_stock_ave_Mgha + 
                            data_final$cum_soc_stock_err_Mgha) * 1.05)) +
      scale_color_manual(
        values = c("#af8dc3", "#762a83", "#7fbf7b"),
        labels = c(paste(control_site, ":", "FD"),
                   paste(treatment_sites, ":", "ESM"),
                   paste(treatment_sites, ":", "FD"))) +
      labs(x = expression(
        paste("Cumulative soil mineral mass, ", m[M], " (g)", sep = "")),
        y = expression(
          paste("Cumulative SOC stock, ", S[C], " (Mg C/ha)", sep = "")),
        color = "",
        title = "i   CEDF SSP2a")
    
    # plot FD-based and ESM-corrected SOC stocks from CEDF and SSP2b
    plot_esm_cedf_stock_ssp2b <-
      ggplot() +
      geom_ribbon(
        data = data_curve %>% 
          dplyr::filter(`site` %in% c(treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("CEDF"),
                        `ssp` %in% c("SSP2b")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            ymin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
            ymax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`)),
        color = "#1b7837", linetype = 0, alpha = 0.2) +
      geom_line(
        data = data_curve %>% 
          dplyr::filter(`site` %in% c(treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("CEDF"),
                        `ssp` %in% c("SSP2b")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            y = `cum_soc_stock_ave_Mgha`),
        color = "#1b7837", linetype = 1) + 
      geom_errorbar(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "CEDF"),
                        `ssp` %in% c("SSP2b")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            ymin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
            ymax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`),
            color = paste(`site`, ": ", `method`, sep = "")),
        width = 20) +
      geom_errorbarh(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "CEDF"),
                        `ssp` %in% c("SSP2b")),
        aes(y = `cum_soc_stock_ave_Mgha`,
            xmin = (`cum_soil_mineral_mass_ave_g` - 
                      `cum_soil_mineral_mass_err_g`),
            xmax = (`cum_soil_mineral_mass_ave_g` + 
                      `cum_soil_mineral_mass_err_g`),
            color = paste(`site`, ": ", `method`, sep = "")), 
        height = 5) +
      geom_point(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "CEDF"),
                        `ssp` %in% c("SSP2b")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            y = `cum_soc_stock_ave_Mgha`,
            color = paste(`site`, ": ", `method`, sep = "")), 
        size = 1, shape = 15) +
      theme_bw() +
      theme(legend.position = "bottom") +
      scale_x_continuous(
        limits = 
          c(0, max(data_final$cum_soil_mineral_mass_ave_g) * 1.05)) +
      scale_y_continuous(
        limits = c(0, max(data_final$cum_soc_stock_ave_Mgha + 
                            data_final$cum_soc_stock_err_Mgha) * 1.05)) +
      scale_color_manual(
        values = c("#af8dc3", "#762a83", "#7fbf7b"),
        labels = c(paste(control_site, ":", "FD"),
                   paste(treatment_sites, ":", "ESM"),
                   paste(treatment_sites, ":", "FD"))) +
      labs(x = expression(
        paste("Cumulative soil mineral mass, ", m[M], " (g)", sep = "")),
        y = expression(
          paste("Cumulative SOC stock, ", S[C], " (Mg C/ha)", sep = "")),
        color = "",
        title = "j   CEDF SSP2b")
    
    # plot FD-based and ESM-corrected SOC stocks from CEDF and SSP3
    plot_esm_cedf_stock_ssp3 <-
      ggplot() +
      geom_ribbon(
        data = data_curve %>% 
          dplyr::filter(`site` %in% c(treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("CEDF"),
                        `ssp` %in% c("SSP3")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            ymin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
            ymax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`)),
        color = "#1b7837", linetype = 0, alpha = 0.2) +
      geom_line(
        data = data_curve %>% 
          dplyr::filter(`site` %in% c(treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("CEDF"),
                        `ssp` %in% c("SSP3")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            y = `cum_soc_stock_ave_Mgha`),
        color = "#1b7837", linetype = 1) + 
      geom_errorbar(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "CEDF"),
                        `ssp` %in% c("SSP3")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            ymin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
            ymax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`),
            color = paste(`site`, ": ", `method`, sep = "")),
        width = 20) +
      geom_errorbarh(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "CEDF"),
                        `ssp` %in% c("SSP3")),
        aes(y = `cum_soc_stock_ave_Mgha`,
            xmin = (`cum_soil_mineral_mass_ave_g` - 
                      `cum_soil_mineral_mass_err_g`),
            xmax = (`cum_soil_mineral_mass_ave_g` + 
                      `cum_soil_mineral_mass_err_g`),
            color = paste(`site`, ": ", `method`, sep = "")), 
        height = 5) +
      geom_point(
        data = data_final %>% 
          dplyr::filter(`site` %in% c(control_site, treatment_sites),
                        `soil_content` %in% c(soil_content),
                        `k` %in% c(k),
                        `regression_function` %in% c("NA", "CEDF"),
                        `ssp` %in% c("SSP3")),
        aes(x = `cum_soil_mineral_mass_ave_g`,
            y = `cum_soc_stock_ave_Mgha`,
            color = paste(`site`, ": ", `method`, sep = "")), 
        size = 1, shape = 15) +
      theme_bw() +
      theme(legend.position = "bottom") +
      scale_x_continuous(
        limits = 
          c(0, max(data_final$cum_soil_mineral_mass_ave_g) * 1.05)) +
      scale_y_continuous(
        limits = c(0, max(data_final$cum_soc_stock_ave_Mgha + 
                            data_final$cum_soc_stock_err_Mgha) * 1.05)) +
      scale_color_manual(
        values = c("#af8dc3", "#762a83", "#7fbf7b"),
        labels = c(paste(control_site, ":", "FD"),
                   paste(treatment_sites, ":", "ESM"),
                   paste(treatment_sites, ":", "FD"))) +
      labs(x = expression(
        paste("Cumulative soil mineral mass, ", m[M], " (g)", sep = "")),
        y = expression(
          paste("Cumulative SOC stock, ", S[C], " (Mg C/ha)", sep = "")),
        color = "",
        title = "k   CEDF SSP3")
    
    # design custom figure layout
    custom_layout <- "
    ABCD
    ABCD
    ABCD
    ABCD
    EFGH
    EFGH
    EFGH
    EFGH
    #IJK
    #IJK
    #IJK
    #IJK
    LLLL
    "
    
    # grid plot of SOC stocks and soil depth
    plot_esm_all_stock <-
      plot_esm_lin_stock_ssp1 +
      plot_esm_lin_stock_ssp2a +
      plot_esm_lin_stock_ssp2b +
      plot_esm_lin_stock_ssp3 +
      plot_esm_mcs_stock_ssp1 +
      plot_esm_mcs_stock_ssp2a +
      plot_esm_mcs_stock_ssp2b +
      plot_esm_mcs_stock_ssp3 +
      plot_esm_cedf_stock_ssp2a +
      plot_esm_cedf_stock_ssp2b +
      plot_esm_cedf_stock_ssp3 +
      plot_layout(axes = "collect",
                  axis_titles = "collect",
                  guides = "collect",
                  design = custom_layout) &
      theme(legend.position = "bottom")
    
    # plot ESM correction curves
    print(plot_esm_all_stock)
    
    # export figure as PNG
    ggsave("./results/fig-protocols.png",
           plot_esm_all_stock,
           width = 7,
           height = 7,
           dpi = 1400)
    
  }

# define relative effect of uncertainties to include the soil sampling protocols
# used in Doehl et al. (2026)
uncertainty_analyser_doehl2026 <-
  function(data_results, control_site, k, soil_content){
    
    # create empty data frame
    data_uncertainty_effect <- data.frame()
    
    # separate data 
    data_error_propagation <- data_results
    
    # get all treatment field labels
    treatment_sites <- 
      unique(data_error_propagation[
        !(data_error_propagation$site %in% c(control_site)),]$site)
    
      # filter for all unique entries of m'th treatment site
      data_uncertainty <- data_error_propagation %>%
        dplyr::filter(`site` %in% c(treatment_sites),
                      `k` %in% c(k),
                      `soil_content` %in% c(soil_content)) %>%
        distinct(`method`, `regression_function`, `uncertainty_off`,
                 `increment`, `ssp`, .keep_all = TRUE)
      
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
          ]$ssp))){
            
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
                       `ssp` %in% 
                         unique(data_uncertainty[
                           data_uncertainty$increment %in% 
                             unique(data_uncertainty[
                               data_uncertainty$regression_function %in%
                                 unique(data_uncertainty$regression_function
                                 )[a],]$increment)[b],
                         ]$ssp)[d]
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
                       `ssp` %in% 
                         unique(data_uncertainty[
                           data_uncertainty$increment %in% 
                             unique(data_uncertainty[
                               data_uncertainty$regression_function %in%
                                 unique(data_uncertainty$regression_function
                                 )[a],]$increment)[b],
                         ]$ssp)[d],
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
                          `ssp` %in% 
                            unique(data_uncertainty[
                              data_uncertainty$increment %in% 
                                unique(data_uncertainty[
                                  data_uncertainty$regression_function %in%
                                    unique(data_uncertainty$regression_function
                                    )[a],]$increment)[b],
                            ]$ssp)[d]
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
      
      # plot relative effect of uncertainties from FD and SSP1
      plot_uncertainty_fd_ssp1 <-
        ggplot(
          data = data_uncertainty_effect %>%
            dplyr::filter(`site` %in% c(treatment_sites),
                          `method` %in% c("FD"),
                          `regression_function` %in% c("NA"),
                          !(`uncertainty_off` %in% soil_content_exclude),
                          `ssp` %in% c("SSP1")),
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
        scale_y_discrete(
          labels = c("20-30" = "0-30",
                     "40-50" = "30-50")) +
        labs(y = "Cumulative sample increments",
             x = expression(
               paste("Relative effect on cumulative SOC stock error (%)",
                     sep = "")),
             fill = "",
             title = "a   FD SSP1")
      
      # plot relative effect of uncertainties from FD and SSP2a
      plot_uncertainty_fd_ssp2a <-
        ggplot(
          data = data_uncertainty_effect %>%
            dplyr::filter(`site` %in% c(treatment_sites),
                          `method` %in% c("FD"),
                          `regression_function` %in% c("NA"),
                          !(`uncertainty_off` %in% soil_content_exclude),
                          `ssp` %in% c("SSP2a")),
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
        scale_y_discrete(
          labels = c("20-30" = "0-30",
                     "30-40" = "30-40",
                     "40-50" = "40-50")) +
        labs(y = "Cumulative sample increments",
             x = expression(
               paste("Relative effect on cumulative SOC stock error (%)",
                     sep = "")),
             fill = "",
             title = "b   FD SSP2a")
      
      # plot relative effect of uncertainties from FD and SSP2b
      plot_uncertainty_fd_ssp2b <-
        ggplot(
          data = data_uncertainty_effect %>%
            dplyr::filter(`site` %in% c(treatment_sites),
                          `method` %in% c("FD"),
                          `regression_function` %in% c("NA"),
                          !(`uncertainty_off` %in% soil_content_exclude),
                          `ssp` %in% c("SSP2b")),
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
        scale_y_discrete(
          labels = c("0-10" = "0-10",
                     "20-30" = "10-30",
                     "40-50" = "30-50")) +
        labs(y = "Cumulative sample increments",
             x = expression(
               paste("Relative effect on cumulative SOC stock error (%)",
                     sep = "")),
             fill = "",
             title = "c   FD SSP2b")
      
      # plot relative effect of uncertainties from FD and SSP3
      plot_uncertainty_fd_ssp3 <-
        ggplot(
          data = data_uncertainty_effect %>%
            dplyr::filter(`site` %in% c(treatment_sites),
                          `method` %in% c("FD"),
                          `regression_function` %in% c("NA"),
                          !(`uncertainty_off` %in% soil_content_exclude),
                          `ssp` %in% c("SSP3")),
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
        scale_y_discrete(
          labels = c("0-10" = "0-10",
                     "10-20" = "10-20",
                     "20-30" = "20-30",
                     "30-40" = "30-40",
                     "40-50" = "40-50")) +
        labs(y = "Cumulative sample increments",
             x = expression(
               paste("Relative effect on cumulative SOC stock error (%)",
                     sep = "")),
             fill = "",
             title = "d   FD SSP3")
      
      # plot relative effect of uncertainties from LIN and SSP1
      plot_uncertainty_esm_lin_ssp1 <-
        ggplot(
          data = data_uncertainty_effect %>%
            dplyr::filter(`site` %in% c(treatment_sites),
                          `regression_function` %in% c("LIN"),
                          `method` %in% c("ESM"),
                          !(`uncertainty_off` %in% soil_content_exclude),
                          `ssp` %in% c("SSP1")),
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
        scale_y_discrete(
          labels = c("20-30" = "0-30",
                     "40-50" = "30-50")) +
        labs(y = "Cumulative sample increments",
             x = expression(
               paste("Relative effect on cumulative SOC stock error (%)",
                     sep = "")),
             fill = "",
             title = "e   LIN SSP1")
      
      # plot relative effect of uncertainties from LIN and SSP2a
      plot_uncertainty_esm_lin_ssp2a <-
        ggplot(
          data = data_uncertainty_effect %>%
            dplyr::filter(`site` %in% c(treatment_sites),
                          `regression_function` %in% c("LIN"),
                          `method` %in% c("ESM"),
                          !(`uncertainty_off` %in% soil_content_exclude),
                          `ssp` %in% c("SSP2a")),
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
        scale_y_discrete(
          labels = c("20-30" = "0-30",
                     "30-40" = "30-40",
                     "40-50" = "40-50")) +
        labs(y = "Cumulative sample increments",
             x = expression(
               paste("Relative effect on cumulative SOC stock error (%)",
                     sep = "")),
             fill = "",
             title = "f   LIN SSP2a")
      
      # plot relative effect of uncertainties from LIN and SSP2b
      plot_uncertainty_esm_lin_ssp2b <-
        ggplot(
          data = data_uncertainty_effect %>%
            dplyr::filter(`site` %in% c(treatment_sites),
                          `regression_function` %in% c("LIN"),
                          `method` %in% c("ESM"),
                          !(`uncertainty_off` %in% soil_content_exclude),
                          `ssp` %in% c("SSP2b")),
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
        scale_y_discrete(
          labels = c("0-10" = "0-10",
                     "20-30" = "10-30",
                     "40-50" = "30-50")) +
        labs(y = "Cumulative sample increments",
             x = expression(
               paste("Relative effect on cumulative SOC stock error (%)",
                     sep = "")),
             fill = "",
             title = "g   LIN SSP2b")
      
      # plot relative effect of uncertainties from LIN and SSP3
      plot_uncertainty_esm_lin_ssp3 <-
        ggplot(
          data = data_uncertainty_effect %>%
            dplyr::filter(`site` %in% c(treatment_sites),
                          `regression_function` %in% c("LIN"),
                          `method` %in% c("ESM"),
                          !(`uncertainty_off` %in% soil_content_exclude),
                          `ssp` %in% c("SSP3")),
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
        scale_y_discrete(
          labels = c("0-10" = "0-10",
                     "10-20" = "10-20",
                     "20-30" = "20-30",
                     "30-40" = "30-40",
                     "40-50" = "40-50")) +
        labs(y = "Cumulative sample increments",
             x = expression(
               paste("Relative effect on cumulative SOC stock error (%)",
                     sep = "")),
             fill = "",
             title = "h   LIN SSP3")
      
      # plot relative effect of uncertainties from MCS and SSP1
      plot_uncertainty_esm_mcs_ssp1 <-
        ggplot(
          data = data_uncertainty_effect %>%
            dplyr::filter(`site` %in% c(treatment_sites),
                          `regression_function` %in% c("MCS"),
                          `method` %in% c("ESM"),
                          !(`uncertainty_off` %in% soil_content_exclude),
                          `ssp` %in% c("SSP1")),
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
        scale_y_discrete(
          labels = c("20-30" = "0-30",
                     "40-50" = "30-50")) +
        labs(y = "Cumulative sample increments",
             x = expression(
               paste("Relative effect on cumulative SOC stock error (%)",
                     sep = "")),
             fill = "",
             title = "i   MCS SSP1")
      
      # plot relative effect of uncertainties from MCS and SSP2a
      plot_uncertainty_esm_mcs_ssp2a <-
        ggplot(
          data = data_uncertainty_effect %>%
            dplyr::filter(`site` %in% c(treatment_sites),
                          `regression_function` %in% c("MCS"),
                          `method` %in% c("ESM"),
                          !(`uncertainty_off` %in% soil_content_exclude),
                          `ssp` %in% c("SSP2a")),
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
        scale_y_discrete(
          labels = c("20-30" = "0-30",
                     "30-40" = "30-40",
                     "40-50" = "40-50")) +
        labs(y = "Cumulative sample increments",
             x = expression(
               paste("Relative effect on cumulative SOC stock error (%)",
                     sep = "")),
             fill = "",
             title = "j   MCS SSP2a")
      
      # plot relative effect of uncertainties from MCS and SSP2b
      plot_uncertainty_esm_mcs_ssp2b <-
        ggplot(
          data = data_uncertainty_effect %>%
            dplyr::filter(`site` %in% c(treatment_sites),
                          `regression_function` %in% c("MCS"),
                          `method` %in% c("ESM"),
                          !(`uncertainty_off` %in% soil_content_exclude),
                          `ssp` %in% c("SSP2b")),
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
        scale_y_discrete(
          labels = c("0-10" = "0-10",
                     "20-30" = "10-30",
                     "40-50" = "30-50")) +
        labs(y = "Cumulative sample increments",
             x = expression(
               paste("Relative effect on cumulative SOC stock error (%)",
                     sep = "")),
             fill = "",
             title = "k   MCS SSP2b")
      
      # plot relative effect of uncertainties from MCS and SSP3
      plot_uncertainty_esm_mcs_ssp3 <-
        ggplot(
          data = data_uncertainty_effect %>%
            dplyr::filter(`site` %in% c(treatment_sites),
                          `regression_function` %in% c("MCS"),
                          `method` %in% c("ESM"),
                          !(`uncertainty_off` %in% soil_content_exclude),
                          `ssp` %in% c("SSP3")),
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
        scale_y_discrete(
          labels = c("0-10" = "0-10",
                     "10-20" = "10-20",
                     "20-30" = "30-30",
                     "30-40" = "30-40",
                     "40-50" = "40-50")) +
        labs(y = "Cumulative sample increments",
             x = expression(
               paste("Relative effect on cumulative SOC stock error (%)",
                     sep = "")),
             fill = "",
             title = "l   MCS SSP3")
      
      # plot relative effect of uncertainties from CEDF and SSP2a
      plot_uncertainty_esm_cedf_ssp2a <-
        ggplot(
          data = data_uncertainty_effect %>%
            dplyr::filter(`site` %in% c(treatment_sites),
                          `regression_function` %in% c("CEDF"),
                          `method` %in% c("ESM"),
                          !(`uncertainty_off` %in% soil_content_exclude),
                          `ssp` %in% c("SSP2a")),
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
        scale_y_discrete(
          labels = c("20-30" = "0-30",
                     "30-40" = "30-40",
                     "40-50" = "40-50")) +
        labs(y = "Cumulative sample increments",
             x = expression(
               paste("Relative effect on cumulative SOC stock error (%)",
                     sep = "")),
             fill = "",
             title = "m   CEDF SSP2a")
      
      # plot relative effect of uncertainties from CEDF and SSP2b
      plot_uncertainty_esm_cedf_ssp2b <-
        ggplot(
          data = data_uncertainty_effect %>%
            dplyr::filter(`site` %in% c(treatment_sites),
                          `regression_function` %in% c("CEDF"),
                          `method` %in% c("ESM"),
                          !(`uncertainty_off` %in% soil_content_exclude),
                          `ssp` %in% c("SSP2b")),
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
        scale_y_discrete(
          labels = c("0-10" = "0-10",
                     "20-30" = "10-30",
                     "40-50" = "30-50")) +
        labs(y = "Cumulative sample increments",
             x = expression(
               paste("Relative effect on cumulative SOC stock error (%)",
                     sep = "")),
             fill = "",
             title = "n   CEDF SSP2b")
      
      # plot relative effect of uncertainties from CEDF and SSP3
      plot_uncertainty_esm_cedf_ssp3 <-
        ggplot(
          data = data_uncertainty_effect %>%
            dplyr::filter(`site` %in% c(treatment_sites),
                          `regression_function` %in% c("CEDF"),
                          `method` %in% c("ESM"),
                          !(`uncertainty_off` %in% soil_content_exclude),
                          `ssp` %in% c("SSP3")),
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
        scale_y_discrete(
          labels = c("0-10" = "0-10",
                     "10-20" = "10-20",
                     "20-30" = "20-30",
                     "30-40" = "30-40",
                     "40-50" = "40-50")) +
        labs(y = "Cumulative sample increments",
             x = expression(
               paste("Relative effect on cumulative SOC stock error (%)",
                     sep = "")),
             fill = "",
             title = "o   CEDF SSP3")
      
      # design custom figure layout
      custom_layout <- "
      ABCD
      ABCD
      ABCD
      ABCD
      EFGH
      EFGH
      EFGH
      EFGH
      IJKL
      IJKL
      IJKL
      IJKL
      #MNO
      #MNO
      #MNO
      #MNO
      #PP#
      "
      
      # grid plot of uncertainty effect by method
      plot_uncertainty <-
        plot_uncertainty_fd_ssp1 +
        plot_uncertainty_fd_ssp2a +
        plot_uncertainty_fd_ssp2b +
        plot_uncertainty_fd_ssp3 +
        plot_uncertainty_esm_lin_ssp1 +
        plot_uncertainty_esm_lin_ssp2a +
        plot_uncertainty_esm_lin_ssp2b +
        plot_uncertainty_esm_lin_ssp3 +
        plot_uncertainty_esm_mcs_ssp1 +
        plot_uncertainty_esm_mcs_ssp2a +
        plot_uncertainty_esm_mcs_ssp2b +
        plot_uncertainty_esm_mcs_ssp3 +
        plot_uncertainty_esm_cedf_ssp2a +
        plot_uncertainty_esm_cedf_ssp2b +
        plot_uncertainty_esm_cedf_ssp3 +
        plot_layout(axes = "collect",
                    axis_titles = "collect",
                    guides = "collect",
                    design = custom_layout) &
        theme(legend.position = "bottom")
      
      # plot relative effect of uncertainties on SOC stock error
      print(plot_uncertainty)
      
      # export figure as PNG
      ggsave("./results/fig-uncertainty.png",
             plot_uncertainty,
             width = 7,
             height = 8,
             dpi = 1400)
  }

# define visualiser of SOC vs SOM vs both as soil content methods to include the 
# soil sampling protocols used in Doehl et al. (2026)
method_analyser_doehl2026 <-
  function(data, control_site, soil_content, k, regression_function, guess_k){
    
    # create empty data frames
    data_final_methods <- data.frame()
    data_curve_methods <- data.frame()
    
    # loop over all soil content methods
    for(i in 1:length(input_soil_content)){
      
      # run ESM correction over every soil content method
      data_esm_results <-
        esm_correction_main_doehl2026(
          data = data_reformat,
          control_site = input_control_site, 
          soil_content = input_soil_content[i],
          k = input_k, 
          regression_function = input_regression_function,
          guess_k = input_guess_k)  
      
      # row bind results
      data_final_methods <-
        rbind(data_final_methods,
              data.frame(data_esm_results[1]))
      data_curve_methods <-
        rbind(data_curve_methods,
              data.frame(data_esm_results[2]))
    }
    
    # separate data 
    data_final <- data_final_methods
    data_curve <- data_curve_methods
    
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
    
    # filter results in final ESM correction points and curves by choice of soil
    # content in method, LIN and SSP2b
    data_final_mod <-
      data_final %>%
      dplyr::filter(`soil_content` %in% c("soc", "som", "both"),
                    `regression_function` %in% c("LIN"),
                    `ssp` %in% c("SSP2b"))
    data_curve_mod <-
      data_curve %>%
      dplyr::filter(`soil_content` %in% c("soc", "som", "both"),
                    `regression_function` %in% c("LIN"),
                    `ssp` %in% c("SSP2b"))
    
    # manually calculate Tukey post-hoc differences at each sampling increment of 
    # SSP2
    tukey_results_1 <-
      tukey_manual_calculator_methods_doehl2026(
        data_results = data_final_mod %>% dplyr::filter(`increment` %in% "0-10"),
        x = "soil_content", 
        y = "cum_soc_stock_ave_Mgha", 
        dy = "cum_soc_stock_err_Mgha")
    
    tukey_results_2 <-
      tukey_manual_calculator_methods_doehl2026(
        data_results = data_final_mod %>% dplyr::filter(`increment` %in% "20-30"),
        x = "soil_content", 
        y = "cum_soc_stock_ave_Mgha", 
        dy = "cum_soc_stock_err_Mgha")
    
    tukey_results_3 <-
      tukey_manual_calculator_methods_doehl2026(
        data_results = data_final_mod %>% dplyr::filter(`increment` %in% "40-50"),
        x = "soil_content", 
        y = "cum_soc_stock_ave_Mgha", 
        dy = "cum_soc_stock_err_Mgha")
    
    # plot the ESM-corrected curves and cumulative SOC stocks by modified sampling
    # increment
    plot_methods_curves <-
      ggplot() +
      geom_line(data = data_curve_mod,
                aes(x = `cum_soil_mineral_mass_ave_g`,
                    y = `cum_soc_stock_ave_Mgha`,
                    color = `soil_content`),
                linetype = 1) + 
      geom_errorbar(data = data_final_mod,
                    aes(x = `cum_soil_mineral_mass_ave_g`,
                        ymin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
                        ymax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`),
                        color = `soil_content`),
                    width = 20) +
      geom_errorbarh(data = data_final_mod,
                     aes(y = `cum_soc_stock_ave_Mgha`,
                         xmin = (`cum_soil_mineral_mass_ave_g` - `cum_soil_mineral_mass_err_g`),
                         xmax = (`cum_soil_mineral_mass_ave_g` + `cum_soil_mineral_mass_err_g`),
                         color = `soil_content`),
                     height = 5) +
      geom_point(data = data_final_mod,
                 aes(x = `cum_soil_mineral_mass_ave_g`,
                     y = `cum_soc_stock_ave_Mgha`,
                     color = `soil_content`),
                 size = 1, shape = 15) +
      theme_bw() +
      theme(legend.position = "bottom") +
      scale_x_continuous(
        limits = c(0, max(data_curve_mod$cum_soil_mineral_mass_ave_g))) +
      scale_y_continuous(
        limits = c(0, max(data_final_mod$cum_soc_stock_ave_Mgha + 
                            data_final_mod$cum_soc_stock_err_Mgha))) +
      scale_color_manual(values = c("soc" = "#762a83",
                                    "som" = "#d8b365",
                                    "both" = "#7fbf7b"),
                         labels = c("soc" = "SOC",
                                    "som" = "SOM",
                                    "both" = "Both")) +
      labs(x = expression(paste("Cumulative soil mineral mass, ", m[M], " (g)", sep = "")),
           y = expression(paste("Cumulative SOC stock, ", S[C], " (Mg C/ha)", sep = "")),
           color = "",
           title = "a")
    
    plot_methods_1 <-
      ggplot(data = data_final_mod %>% 
               dplyr::filter(`increment` %in% c("0-10"))) +
      geom_bar(aes(y = `soil_content`,
                   x = `cum_soc_stock_ave_Mgha`,
                   fill = `soil_content`), 
               stat = "identity", position = position_dodge()) +
      geom_errorbarh(aes(y = `soil_content`,
                         xmin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
                         xmax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`)),
                     color = "black", width = 0.5, position = position_dodge()) +
      geom_text(data = tukey_results_1[[2]] %>% 
                  dplyr::rename(soil_content = `treatment`) %>%
                  dplyr::mutate(
                    place = 1.1 * 
                      max(data_final_mod[data_final_mod$increment %in% c("0-10"),
                      ]$cum_soc_stock_ave_Mgha + 
                        data_final_mod[data_final_mod$increment %in% c("0-10"),
                        ]$cum_soc_stock_err_Mgha)),
                aes(x = `place`, y = `soil_content`, label = `group`)) +
      theme_bw() +
      theme(legend.position = "none") +
      scale_fill_manual(values = c("soc" = "#762a83",
                                   "som" = "#d8b365",
                                   "both" = "#7fbf7b"),
                        labels = c("soc" = "SOC",
                                   "som" = "SOM",
                                   "both" = "Both")) +
      scale_y_discrete(labels = c("soc" = "SOC",
                                  "som" = "SOM",
                                  "both" = "Both")) +
      labs(y = "",
           x = expression(paste("Cumulative SOC stock, ", S[C], " (Mg C/ha)", sep = "")),
           fill = "",
           title = "b")
    
    plot_methods_2 <-
      ggplot(data = data_final_mod %>% 
               dplyr::filter(`increment` %in% c("20-30"))) +
      geom_bar(aes(y = `soil_content`,
                   x = `cum_soc_stock_ave_Mgha`,
                   fill = `soil_content`), 
               stat = "identity", position = position_dodge()) +
      geom_errorbarh(aes(y = `soil_content`,
                         xmin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
                         xmax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`)),
                     color = "black", width = 0.5, position = position_dodge()) +
      geom_text(data = tukey_results_1[[2]] %>% 
                  dplyr::rename(soil_content = `treatment`) %>%
                  dplyr::mutate(
                    place = 1.1 * 
                      max(data_final_mod[data_final_mod$increment %in% c("20-30"),
                      ]$cum_soc_stock_ave_Mgha + 
                        data_final_mod[data_final_mod$increment %in% c("20-30"),
                        ]$cum_soc_stock_err_Mgha)),
                aes(x = `place`, y = `soil_content`, label = `group`)) +
      theme_bw() +
      theme(legend.position = "none") +
      scale_fill_manual(values = c("soc" = "#762a83",
                                   "som" = "#d8b365",
                                   "both" = "#7fbf7b"),
                        labels = c("soc" = "SOC",
                                   "som" = "SOM",
                                   "both" = "Both")) +
      scale_y_discrete(labels = c("soc" = "SOC",
                                  "som" = "SOM",
                                  "both" = "Both")) +
      labs(y = "",
           x = expression(paste("Cumulative SOC stock, ", S[C], " (Mg C/ha)", sep = "")),
           fill = "",
           title = "c")
    
    plot_methods_3 <-
      ggplot(data = data_final_mod %>% 
               dplyr::filter(`increment` %in% c("40-50"))) +
      geom_bar(aes(y = `soil_content`,
                   x = `cum_soc_stock_ave_Mgha`,
                   fill = `soil_content`), 
               stat = "identity", position = position_dodge()) +
      geom_errorbarh(aes(y = `soil_content`,
                         xmin = (`cum_soc_stock_ave_Mgha` - `cum_soc_stock_err_Mgha`),
                         xmax = (`cum_soc_stock_ave_Mgha` + `cum_soc_stock_err_Mgha`)),
                     color = "black", width = 0.5, position = position_dodge()) +
      geom_text(data = tukey_results_1[[2]] %>% 
                  dplyr::rename(soil_content = `treatment`) %>%
                  dplyr::mutate(
                    place = 1.1 * 
                      max(data_final_mod[data_final_mod$increment %in% c("40-50"),
                      ]$cum_soc_stock_ave_Mgha + 
                        data_final_mod[data_final_mod$increment %in% c("40-50"),
                        ]$cum_soc_stock_err_Mgha)),
                aes(x = `place`, y = `soil_content`, label = `group`)) +
      theme_bw() +
      theme(legend.position = "none") +
      scale_fill_manual(values = c("soc" = "#762a83",
                                   "som" = "#d8b365",
                                   "both" = "#7fbf7b"),
                        labels = c("soc" = "SOC",
                                   "som" = "SOM",
                                   "both" = "Both")) +
      scale_y_discrete(labels = c("soc" = "SOC",
                                  "som" = "SOM",
                                  "both" = "Both")) +
      labs(y = "",
           x = expression(paste("Cumulative SOC stock, ", S[C], " (Mg C/ha)", sep = "")),
           fill = "",
           title = "d")
    
    custom_layout_figure_methods <- "
    AAABB
    AAACC
    AAADD
    EEE##
    "
    
    # 
    plot_soc_vs_som <-
      plot_methods_curves +
      plot_methods_1 +
      plot_methods_2 +
      plot_methods_3 +
      guide_area() +
      plot_layout(axes = "collect",
                  axis_titles = "collect",
                  guides = "collect",
                  design = custom_layout_figure_methods)
    
    # export figure as PNG
    ggsave("./results/fig-soc_vs_som.png",
           plot_soc_vs_som,
           width = 7,
           height = 7,
           dpi = 1400)
  }

# define manual Tukey post-hoc calculator to include the soil sampling protocols
# used in Doehl et al. (2026)
tukey_manual_calculator_methods_doehl2026 <- 
  function(data_results, x, y, dy){
    
    tryCatch({
      
      # get number of cores per method group, number of method groups and total count
      n = mean(data_results$sample_size)
      k = nrow(data_results)
      N = sum(data_results$sample_size)
      
      # calculate residual sums of squared errors, mean squared errors, Tukey 
      # q-value and HSD value
      sse <- sum((data_results$sample_size - 1) * 
                   data_results$sample_size * (data_results[[dy]] ^ 2), na.rm = TRUE)
      mse <- sse / (N - k)
      q_value <- qtukey(p = 0.95, nmeans = k, df = N - k)
      tukey_hsd <- q_value * sqrt(mse / n)
      
      # calculate mean of all method groups
      means <- tapply(data_results[[y]], data_results[[x]], mean)
      
      # calculate inverse additions of sample sizes into a matrix
      n_inv <- as.data.frame(matrix(0, nrow = k, ncol = k))
      rownames(n_inv) <- data_results[[x]]
      colnames(n_inv) <- data_results[[x]]
      for(i in 1:nrow(n_inv)){
        for(j in 1:ncol(n_inv)){
          n_inv[i,j] <- 1 / data_results$sample_size[i] + 
            1 / data_results$sample_size[j]
        }
      }
      
      # calculate difference between method group means, lower and upper
      # 95%-CIs and labeling which differences are statistically significant
      diff_mean <- as.data.frame(as.matrix(dist(means), labels = TRUE))
      diff_upper <- diff_mean + q_value * sqrt(mse * n_inv / 2)
      diff_lower <- diff_mean - q_value * sqrt(mse * n_inv / 2)
      diff_sig <- data.frame(diff_mean >= tukey_hsd)
      
      # reform data frames from matrices to columns
      diff_tukey_results <-
        as.data.frame(cbind(rep(names(diff_mean), each = nrow(diff_mean)), 
                            rep(names(diff_mean), times = ncol(diff_mean)),
                            unlist(diff_mean, use.names = FALSE), 
                            unlist(diff_lower, use.names = FALSE), 
                            unlist(diff_upper, use.names = FALSE),
                            unlist(diff_sig, use.names = FALSE))) %>%
        rename("appr 1" = `V1`, "appr 2" = `V2`,
               "diff" = `V3`, "lwr" = `V4`, "upr" = `V5`, "sig" = `V6`) %>%
        distinct(`diff`, .keep_all = TRUE) %>%
        mutate_at(c("diff", "lwr", "upr"), as.numeric)
      
      # remove any rows where approach 1 and 2 values are identical
      diff_tukey_results_unique <- 
        diff_tukey_results[diff_tukey_results$`appr 1` != 
                             diff_tukey_results$`appr 2`,]
      
      # calculate adjusted p-values of differences of means
      diff_p_adj <- 
        1 - ptukey(q = as.numeric(diff_tukey_results_unique$diff) / sqrt(mse / n), 
                   nmeans = k, df = N - k)
      
      # bind adjusted p values to data frame
      diff_tukey_final <-
        diff_tukey_results_unique %>% 
        mutate(`p adj` = diff_p_adj, .before = `sig`) %>%
        mutate_at("p adj", as.numeric) %>%
        mutate(`p sig` = `p adj`) %>%
        mutate(`p sig` = replace(`p sig`, `p sig` == 0, 5)) %>%
        mutate(`p sig` = replace(`p sig`, `p sig` <= 1 & `p sig` > 0.05, 6)) %>%
        mutate(`p sig` = replace(`p sig`, `p sig` <= 0.05 & `p sig` > 0.01, 7)) %>%
        mutate(`p sig` = replace(`p sig`, `p sig` <= 0.01 & `p sig` > 0.001, 8)) %>%
        mutate(`p sig` = replace(`p sig`, `p sig` <= 0.001, 9)) %>%
        mutate(`p sig` = as.character(`p sig`)) %>%
        mutate(`p sig` = replace(`p sig`, `p sig` %in% "5", "NA")) %>%
        mutate(`p sig` = replace(`p sig`, `p sig` %in% "6", "ns")) %>%
        mutate(`p sig` = replace(`p sig`, `p sig` %in% "7" , "*")) %>%
        mutate(`p sig` = replace(`p sig`, `p sig` %in% "8" , "**")) %>%
        mutate(`p sig` = replace(`p sig`, `p sig` %in% "9" , "***")) %>%
        mutate(`lab 1` = "ab", `lab 2` = "ab") %>%
        mutate(`lab 1` = replace(`lab 1`, `p sig` %in% c("*", "**", "***"), "a"),
               `lab 2` = replace(`lab 2`, `p sig` %in% c("*", "**", "***"), "b"))
      
      # assign compact letter display manually
      tukey_lab <- data.frame()
      tukey_lab_pre <-
        data.frame(treatment = c(diff_tukey_final$`appr 1`, diff_tukey_final$`appr 2`),
                   group = c(diff_tukey_final$`lab 1`, diff_tukey_final$`lab 2`)) %>%
        distinct(`treatment`, `group`)
      for(i in 1:length(unique(tukey_lab_pre$treatment))){
        if(sum(tukey_lab_pre$treatment %in% unique(tukey_lab_pre$treatment)[i], na.rm = TRUE) %in% 1){
          tukey_lab <- rbind(tukey_lab, tukey_lab_pre[tukey_lab_pre$treatment %in% unique(tukey_lab_pre$treatment)[i],])
        } else{
          tukey_lab <- rbind(tukey_lab, tukey_lab_pre[tukey_lab_pre$treatment %in% unique(tukey_lab_pre$treatment)[i] &
                                                        !(tukey_lab_pre$group %in% c("ab")),])
        }
      }
      
      #print(data_reduced)
      #print(n)
      #print(k)
      #print(sse)
      #print(mse)
      #print(q_value)
      #print(tukey_hsd)
      #print(means)
      
      # print final results
      print("   Manual Tukey multiple comparisons of means")
      print("   95% family-wise confidence level")
      print("")
      print(diff_tukey_final)
      
      # plot Tukey results
      plot_tukey_comparison <-
        ggplot(data = diff_tukey_final) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_point(aes(y = paste(`appr 1`, "-", `appr 2`, sep = ""),
                       x = `diff`), shape = 12) +
        geom_errorbarh(aes(y = paste(`appr 1`, "-", `appr 2`, sep = ""),
                           xmin = `lwr`, xmax = `upr`)) +
        theme_bw() +
        scale_x_continuous(
          breaks = seq(round(min(diff_tukey_final$lwr) - 1), 
                       round(max(diff_tukey_final$upr) + 1), 
                       10)) +
        scale_y_discrete(limits = rev) +
        labs(x = "Differences in mean levels of group",
             y = " ",
             title = paste("95% family-wise confidence level", sep = ""))
      
      print(plot_tukey_comparison)
      
      return(list(diff_tukey_final, tukey_lab))
    }, 
    error = function(cond) {
      message("Error")
      message("Skipping comparison.")
      message("Here is the original error in method:")
      message(conditionMessage(cond))
      NA},
    warning = function(cond) {
      message("Warning")
      message("Skipping comparison.")
      message("Here is the original warning in method:")
      message(conditionMessage(cond))
      NULL},
    finally = {
      message("Completed comparison.")
      message("")
    })
  }

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
data_reformat <- data_reformat %>%
  dplyr::filter(!(`core` %in% c("W13", "GF3")))

### Define inputs ##############################################################

# define which is the control site to which soil mineral mass are equated to
input_control_site <- c("control")

# define which content method to calculate SOC stock and soil mineral mass with
input_soil_content <- c("soc", "som", "both")

# define which van Bemmelen factor assumption to use
# if guesssing the factor value, define that number
input_k <- c("conventional")
input_guess_k <- 1.724

# define which regression function(s) to use for ESM correction
input_regression_function <- c("LIN", "MCS", "CEDF")

### Execute ESM correction #####################################################

# run ESM correction over all soil content methods to include the soil sampling 
# protocols used in Doehl et al. (2026)
method_analyser_doehl2026(
  data = data_reformat,
  control_site = input_control_site,
  soil_content = input_soil_content,
  k = input_k,
  regression_function = input_regression_function,
  guess_k = input_guess_k)