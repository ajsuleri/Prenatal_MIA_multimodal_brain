#############################################################
##########Project: MIA & multimodal brain development########
#############################################################

#Author: Anna Suleri

#In this script all functions are noted that we use in the project. We will then source this script in the other scripts.

###########################################################
#------------Load data based on file extension------------#
###########################################################
# load foreign and haven libraries 

load_files <- function(folder, pattern, read_function) {
  files <- list.files(folder, pattern = pattern, full.names = TRUE)
  
  for (file in files) {
    object_name <- tools::file_path_sans_ext(basename(file))
    assign(object_name, read_function(file), envir = .GlobalEnv)
  }
}

#\ 

################################################################
#------------Calculate between network connectivity------------#
################################################################

calculate_btnetworks <- function(network){
  cols <- grep(network, names(between_con_f9_correct), value = TRUE)
  
  mean_connectivity_df <- between_con_f9_correct %>%
    rowwise() %>%
    mutate(mean_value = mean(c_across(all_of(cols)), na.rm = TRUE)) %>%
    dplyr::select(idc, mean_value)
  
  return(mean_connectivity_df)
}

#\ 

################################################################
#------------Merge sMRI hemispheres------------#
################################################################

merge_smri_hemispheres <- function(structs, data) {
  for (x in structs) {
    for (y in 1:ncol(data)) {
      a <- colnames(data)[y]
      for (z in 1:ncol(data)) {
        b <- colnames(data)[z]
        
        if (a == paste0("Left_", x, "_f09") & b == paste0("Right_", x, "_f09")) {
          newcolname <- paste0(x, "_subcortical_f09")
          data[, newcolname] <- (data[, y] + data[, z]) / 2
        }
        
        if (a == paste0("lh", x, "_f09") & b == paste0("rh", x, "_f09")) {
          newcolname <- paste0(x, "_cortical_f09")
          data[, newcolname] <- (data[, y] + data[, z]) / 2
        }
        
        if (a == paste0("Left_", x, "_f13") & b == paste0("Right_", x, "_f13")) {
          newcolname <- paste0(x, "_subcortical_f13")
          data[, newcolname] <- (data[, y] + data[, z]) / 2
        }
        
        if (a == paste0("lh", x, "_f13") & b == paste0("rh", x, "_f13")) {
          newcolname <- paste0(x, "_cortical_f13")
          data[, newcolname] <- (data[, y] + data[, z]) / 2
        }
      }
    }
  }
  return(data)
}

#\

###################################################
#-----------Convert factor to numeric ------------#
##################################################

transform_to_numeric <- function(x) {
  as.numeric(as.character(x))  
}

#\ 

########################################
#------------Baseline table------------#
########################################

baseline_table_function <- function(baselinevars, df){ #feed baselinevars vector and dataframe 
  
  #create empty dataframe to store results 
  results <- data.frame(Variable = character(0), Output = character(0))
  
  #for each baseline variable as specified in baselinevars decide whether continous or categorical and give output accordingly 
  for(i in baselinevars){ 
    
    #x = i vars that are columns in the dataframe df
    x <- df[, i] 
    
    #show column name as heading per output 
    message(i) 
    
    #function for continuous variables
    summary_continuous <- function(x){
      standev <- sd(x, na.rm = T)
      meanvar <- mean(x, na.rm = T)
      print(paste0(round(meanvar, 1), '(', round(standev, 1), ')'))
    }
    
    #function for categorical variables 
    summary_categorical <- function(x){
      tab1 <- prop.table(table(x, useNA = 'always'))
      tab2 <- table(x, useNA = "always")
      print(paste(round(tab1 * 100, 1), '%', names(tab1), collapse = ','))
      print(paste(tab2, names(tab2)))
    }
    
    #if else to apply correct function for vars type 
    if (class(x) == 'numeric') {
      output <- summary_continuous(x)
    } 
    else 
    {
      output <- summary_categorical(x) 
    }
    results <- rbind(results, data.frame(Variable = i, Output = output))
  }
  write_xlsx(results, 'baseline_table_results.xlsx') #load write_xl package
}

#\ 

########################################
#------Summary of missing values-------#
########################################

miss_values_function <- function(df){
  missvalues <- cbind("# NA" = sort(colSums(is.na(df))),
                      "% NA" = round(sort(colMeans(is.na(df))) * 100, 2))
  print(missvalues)
}

#\ 

########################################
#---------Multiple imputation----------#
########################################
# first load mice library 

imputation_function <- function(data, exclude_imp_vars, exclude_predictors, method = "default") {
  
  # Running setup imputation run
  if (method == "rf") { #e.g. in case of brain to allow any type of combination between variables 
    imp0 <- mice(data, maxit = 0, method = "rf")
  } else {
    imp0 <- mice(data, maxit = 0, defaultMethod = c("norm", "logreg", "polyreg", "polr"))
  }
  
  # Imputation method matrix
  meth <- imp0$method
  meth[exclude_imp_vars] <- ""
  
  # Predictor matrix
  pred <- imp0$predictorMatrix
  pred[exclude_predictors] <- 0
  
  # Visit sequence
  visSeq <- imp0$visitSequence
  
  # Performing the imputation
  imp.test <- mice(data, method = meth, predictorMatrix = pred, visitSequence = visSeq, maxit = 30, m = 30, printFlag = TRUE, seed = 2024)  
  
  #get the dataframe name 
  dataName <- deparse(substitute(data))
  
  #assigning a dynamic name to the imp.test object
  imputedDataName <- paste0("imputedData_", dataName)
  assign(imputedDataName, imp.test)
  
  # Saving the imputed dataset as .RDS file
  saveRDS(imp.test, file = paste0(imputedDataName, ".rds"))
  
  # Output 
  return(list(imp0 = imp0, imp.test = imp.test))
}

#\ 

#############################################
#---------Check model assumptions----------#
############################################

check_model_assumptions <- function(model) {
  # Check multicollinearity using VIF
  vif_values <- car::vif(model)
  print("Variance Inflation Factors (VIF):")
  print(vif_values)
  
  # Distribution of residuals
  res <- resid(model)
  par(mfrow = c(2, 2))  # Create a 2x2 grid of plots
  qqnorm(res)
  qqline(res)
  plot(density(res), main = "Density Plot of Residuals")
  
  # Heteroscedasticity assumption
  plot(fitted(model), res, main = "Residuals vs. Fitted Values")
  abline(h = 0, lty = 2, col = "red")
}

#\

#############################################
#---------Apply FDR-BH correction----------#
############################################

# Apply FDR correction and resave data frames to Excel
add_fdr_pvalue <- function(data_frame, pvalue_column_index) {
  pval <- unlist(data_frame[, pvalue_column_index])
  data_frame$fdr_pvalue <- p.adjust(pval, method = 'fdr')
  return(data_frame)
}

#\ 
