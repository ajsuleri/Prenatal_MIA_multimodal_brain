#############################################################
##########Project: MIA & multimodal brain development########
#############################################################

## Author: Anna Suleri

## General note: After preparing the dataframe in the first script, we will now run the statistical analyses in this script. In brief, we will run mixed-effects models for MIA & CRP for sMRI, DTI and fMRI. We will examine both the main as well as the interaction effect. In  the interaction model we will also check for non-linearity as well as run a model additionally adjusting for icv in the sMRI models. 
## We will further test for if child sex is a moderating factor as well as conducting the analyses for exposure timepoint 1 vs exposure timepoint 2 seperately. And finally we will conduct a sensitivity analysis in which we will additionally adjust for maternal inflammatory conditions, maternal inflammatory factors, and birth complications. 

#----------------------------------------------------#
####--------------Step 0: PREPARATION-------------####
#---------------------------------------------------#

# Clear environment
rm(list = ls()) 

# Set seed 
set.seed(181223)

# Load functions 
setwd("path_to_script")

source('0.Functions.R')

# Load libraries 
libraries <- c('foreign', 'haven', 'tidyverse' ,'mice', 'lattice', 'ggplot2', 'writexl', 'corrplot', 'stringi', 'miceadds', 'mitools', 'CBPS', 'survey', 'survival', 'psych', 'plotly', 'sjPlot', 'RColorBrewer', 'lme4', 'ggpubr', 'broom.mixed', 'ggseg')

invisible(lapply(libraries, require, character.only = T))

# Load imputed data after ipw and further cleaning
setwd("path_to_data")

df_final_imputed <- readRDS('ipw_mids_final_cleaned.rds')

#\

#------------------------------------------------------------#
####-------------Step 1: Mixed-effects models-------------####
#------------------------------------------------------------#

# ========================================================== #
##                    Model assumptions 
# ========================================================== #

# Select last imputed df (that has persumably best convergence) to test assumptions on
assumptions_single_df <- complete(df_final_imputed, 30)

# Specify mixed-effects models for total brain volume and exposure
fit_mia <- lmer(scale(genr_tbv) ~ scale(mia_index_average) + age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average +  BatchDay_average + (1| IDC), weights = ipw_weights ,data = assumptions_single_df)

fit_crp <- lmer(scale(genr_tbv) ~ scale(hs_crp_average) + age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average +  BatchDay_average + (1| IDC), weights = ipw_weights ,data = assumptions_single_df)

# Check model assumptions 
check_model_assumptions(fit_mia) 

check_model_assumptions(fit_crp) 

#\

# ========================================================== #
##                    structural MRI 
# ========================================================== #

## Create vector of exposures 
exposures_mia <- c('mia_index_average')
exposures_crp <- c('hs_crp_average')
exposures_t1_mia <- c( 'mia_index_t1')
exposures_t2_mia <- c( 'mia_index_t2')
exposures_t1_crp <- c('HsCRPmgL_g1')
exposures_t2_crp <- c('HsCRPmgL_g2')

## Create vector of outcomes
smri_outcomes <- c('genr_tbv', 'CortexVol', 'TotalGrayVol', 'CerebralWhiteMatterVol' ,'Cerebellum_Cortex_vol_subcortical', 'Thalamus_Proper_vol_subcortical', 'Hippocampus_vol_subcortical', 'Amygdala_vol_subcortical', 'Lateral_Ventricle_vol_subcortical', 'frontal_vol_cortical', 'cingulate_vol_cortical', 'occipital_vol_cortical', 'temporal_vol_cortical', 'parietal_vol_cortical', 'insula_vol_cortical')

## Create empty dataframes
results_smri_main_mia <- data.frame()
results_smri_main_crp <- data.frame()

results_smri_int_mia <- data.frame()
results_smri_int_crp <- data.frame()

results_smri_int_icv_mia <- data.frame()
results_smri_int_icv_crp <- data.frame()

results_smri_sex_mia <- data.frame()
results_smri_sex_crp <- data.frame()

results_smri_sens_mia <- data.frame()
results_smri_sens_crp <- data.frame()

results_smri_t1_mia <- data.frame()
results_smri_t1_crp <- data.frame()

results_smri_t2_mia <- data.frame()
results_smri_t2_crp <- data.frame()

#fit <- with(df_final_imputed, lmer(scale(genr_tbv) ~ scale(mia_index_average)*age_child_mri + GENDER + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + INCOME + PARITY + GestagePLg_average +  BatchDay_average + pm + im + GESTBIR + BMI_0 + (1| IDC), weights = ipw_weights)) 
#fit2 <- with(df_final_imputed, lmer(scale(genr_tbv) ~ ns(scale(mia_index_average)*age_child_mri, 3) + GENDER + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + INCOME + PARITY + GestagePLg_average +  BatchDay_average + pm + im + GESTBIR + BMI_0 + (1| IDC), weights = ipw_weights)) 
#anova(fit, fit2)

# ========================================================== #
### Maternal immune activation 

## Apply all models and loop over average exposures and outcomes 
for (outcome in smri_outcomes) {
  
  for (exposure in exposures_mia) {
    
    # Model formulas
    a <- paste0('scale(',outcome,') ~ scale(',exposure,') + age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average +  BatchDay_average + (1| IDC)')  # main model
    
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average +  BatchDay_average + (1| IDC)')  # interaction model
   
     c <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + scale(eTIV) + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average +  BatchDay_average + (1| IDC)')  # interaction model + icv 
   
      d <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri*GENDER + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + INCOME + PARITY + GestagePLg_average +  BatchDay_average + (1| IDC)')  # interaction model + sex
  
        e <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average +  BatchDay_average  + pm + im + GESTBIR + BMI_0 + (1| IDC)')  # interaction model + sens covars
    
    # Calculating beta, confidence interval, and p-value for all models
    main_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(a), weights = ipw_weights))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
   
     interaction_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[17, c(2,3, 6, 7, 8)]
   
      icv_interaction_output <- summary(pool(with(df_final_imputed, lmer(as.formula(c), weights = ipw_weights))), conf.int = TRUE)[18, c(2,3, 6, 7, 8)]
   
       sex_interaction_output <- summary(pool(with(df_final_imputed, lmer(as.formula(d), weights = ipw_weights))), conf.int = TRUE)[20, c(2,3, 6, 7, 8)]
    
       sens_covars_interaction_output <- summary(pool(with(df_final_imputed, lmer(as.formula(e), weights = ipw_weights))), conf.int = TRUE)[21, c(2, 3,6, 7, 8)]
    
    # Adding outcome and exposure information
    main_model_output$outcome <- outcome
    main_model_output$exposure <- exposure
    
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    
    icv_interaction_output$outcome <- outcome
    icv_interaction_output$exposure <- exposure
    
    sex_interaction_output$outcome <- outcome
    sex_interaction_output$exposure <- exposure
    
    sens_covars_interaction_output$outcome <- outcome
    sens_covars_interaction_output$exposure <- exposure
    
    # Adding results to data frames
    results_smri_main_mia <- rbind(results_smri_main_mia, main_model_output)
    results_smri_int_mia <- rbind(results_smri_int_mia, interaction_model_output)
    results_smri_int_icv_mia <- rbind(results_smri_int_icv_mia, icv_interaction_output)
    results_smri_sex_mia <- rbind(results_smri_sex_mia, sex_interaction_output)
    results_smri_sens_mia <- rbind(results_smri_sens_mia, sens_covars_interaction_output)
  }
}

# Adding column names to created dataframes 
colnames(results_smri_main_mia) <- c("beta", "SE" , "pval","lowerCI", "upperCI", "Outcome", 'Exposure')
colnames(results_smri_int_mia) <- c("beta", "SE" ,"pval","lowerCI", "upperCI", "Outcome", 'Exposure')
colnames(results_smri_int_icv_mia) <- c("beta", "SE" ,"pval","lowerCI", "upperCI","Outcome", 'Exposure')
colnames(results_smri_sex_mia) <- c("beta", "SE", "pval","lowerCI", "upperCI","Outcome", 'Exposure')
colnames(results_smri_sens_mia) <- c("beta", "SE","pval" ,"lowerCI", "upperCI","Outcome", 'Exposure')

# Clean up tables to 3 decimals
results_smri_main_mia[, 1:5] <- round(results_smri_main_mia[, 1:5], digits = 3)
results_smri_int_mia[, 1:5] <- round(results_smri_int_mia[, 1:5], digits = 3)
results_smri_int_icv_mia[, 1:5] <- round(results_smri_int_icv_mia[, 1:5], digits = 3)
results_smri_sex_mia[, 1:5] <- round(results_smri_sex_mia[, 1:5], digits = 3)
results_smri_sens_mia[, 1:5] <- round(results_smri_sens_mia[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
write_xlsx(results_smri_main_mia, 'results_smri_main_mia.xlsx')
write_xlsx(results_smri_int_mia, 'results_smri_int_mia.xlsx')
write_xlsx(results_smri_int_icv_mia, 'results_smri_int_icv_mia.xlsx')
write_xlsx(results_smri_sex_mia, 'results_smri_sex_mia.xlsx')
write_xlsx(results_smri_sens_mia, 'results_smri_sens_mia.xlsx')

## Apply all models and loop over time specific exposures and outcomes
# T1
for (outcome in smri_outcomes) {
  for (exposure in exposures_t1_mia) {
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg1 +  BatchDayg1 + (1| IDC)')  
    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[17, c(2,3, 6, 7, 8)]
    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    # Adding results to data frames
    results_smri_t1_mia <- rbind(results_smri_t1_mia, interaction_model_output)
  }
}

# T2
for (outcome in smri_outcomes) {
  for (exposure in exposures_t2_mia) {
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg2 +  BatchDayg2 + (1| IDC)')  
    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[17, c(2,3, 6, 7, 8)]
    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    # Adding results to data frames
    results_smri_t2_mia <- rbind(results_smri_t2_mia, interaction_model_output)
  }
}

# Adding column names to created dataframes 
colnames(results_smri_t1_mia) <- c("beta", "SE" , "pval","lowerCI", "upperCI", "Outcome", 'Exposure')
colnames(results_smri_t2_mia) <- c("beta", "SE" , "pval","lowerCI", "upperCI", "Outcome", 'Exposure')

# Clean up tables to 3 decimals
results_smri_t1_mia[, 1:5] <- round(results_smri_t1_mia[, 1:5], digits = 3)
results_smri_t2_mia[, 1:5] <- round(results_smri_t2_mia[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
write_xlsx(results_smri_t1_mia, 'results_smri_t1_mia.xlsx')
write_xlsx(results_smri_t2_mia, 'results_smri_t2_mia.xlsx')

# ========================================================== #

### C-reactive protein 

## Apply all models and loop over average exposures and outcomes 
for (outcome in smri_outcomes) {
  
  for (exposure in exposures_crp) {
    
    # Model formulas
    a <- paste0('scale(',outcome,') ~ scale(',exposure,') + age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average + (1| IDC)')  # main model
    
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average + (1| IDC)')  # interaction model
    
    c <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + scale(eTIV) + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average + (1| IDC)')  # interaction model + icv 
    
    d <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri*GENDER + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + INCOME + PARITY + GestagePLg_average + (1| IDC)')  # interaction model + sex
    
    e <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average + pm + im + GESTBIR + BMI_0 + (1| IDC)')  # interaction model + sens covars
    
    # Calculating beta, confidence interval, and p-value for all models
    main_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(a), weights = ipw_weights))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
    
    interaction_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[16, c(2,3, 6, 7, 8)]
    
    icv_interaction_output <- summary(pool(with(df_final_imputed, lmer(as.formula(c), weights = ipw_weights))), conf.int = TRUE)[17, c(2,3, 6, 7, 8)]
    
    sex_interaction_output <- summary(pool(with(df_final_imputed, lmer(as.formula(d), weights = ipw_weights))), conf.int = TRUE)[19, c(2,3, 6, 7, 8)]
    
    sens_covars_interaction_output <- summary(pool(with(df_final_imputed, lmer(as.formula(e), weights = ipw_weights))), conf.int = TRUE)[20, c(2, 3,6, 7, 8)]
    
    # Adding outcome and exposure information
    main_model_output$outcome <- outcome
    main_model_output$exposure <- exposure
    
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    
    icv_interaction_output$outcome <- outcome
    icv_interaction_output$exposure <- exposure
    
    sex_interaction_output$outcome <- outcome
    sex_interaction_output$exposure <- exposure
    
    sens_covars_interaction_output$outcome <- outcome
    sens_covars_interaction_output$exposure <- exposure
    
    # Adding results to data frames
    results_smri_main_crp <- rbind(results_smri_main_crp, main_model_output)
    results_smri_int_crp <- rbind(results_smri_int_crp, interaction_model_output)
    results_smri_int_icv_crp <- rbind(results_smri_int_icv_crp, icv_interaction_output)
    results_smri_sex_crp <- rbind(results_smri_sex_crp, sex_interaction_output)
    results_smri_sens_crp <- rbind(results_smri_sens_crp, sens_covars_interaction_output)
  }
}

# Adding column names to created dataframes 
colnames(results_smri_main_crp) <- c("beta", "SE" , "pval","lowerCI", "upperCI", "Outcome", 'Exposure')
colnames(results_smri_int_crp) <- c("beta", "SE" ,"pval","lowerCI", "upperCI", "Outcome", 'Exposure')
colnames(results_smri_int_icv_crp) <- c("beta", "SE" ,"pval","lowerCI", "upperCI","Outcome", 'Exposure')
colnames(results_smri_sex_crp) <- c("beta", "SE", "pval","lowerCI", "upperCI","Outcome", 'Exposure')
colnames(results_smri_sens_crp) <- c("beta", "SE","pval" ,"lowerCI", "upperCI","Outcome", 'Exposure')

# Clean up tables to 3 decimals
results_smri_main_crp[, 1:5] <- round(results_smri_main_crp[, 1:5], digits = 3)
results_smri_int_crp[, 1:5] <- round(results_smri_int_crp[, 1:5], digits = 3)
results_smri_int_icv_crp[, 1:5] <- round(results_smri_int_icv_crp[, 1:5], digits = 3)
results_smri_sex_crp[, 1:5] <- round(results_smri_sex_crp[, 1:5], digits = 3)
results_smri_sens_crp[, 1:5] <- round(results_smri_sens_crp[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
write_xlsx(results_smri_main_crp, 'results_smri_main_crp.xlsx')
write_xlsx(results_smri_int_crp, 'results_smri_int_crp.xlsx')
write_xlsx(results_smri_int_icv_crp, 'results_smri_int_icv_crp.xlsx')
results_smri_sex_crp_fdr <- add_fdr_pvalue(results_smri_sex_crp, 3)
write_xlsx(results_smri_sex_crp_fdr, 'results_smri_sex_crp.xlsx')
write_xlsx(results_smri_sens_crp, 'results_smri_sens_crp.xlsx')

## Apply all models and loop over time specific exposures and outcomes
# T1
for (outcome in smri_outcomes) {
  for (exposure in exposures_t1_crp) {
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + gestage_plasma_g1 + (1| IDC)')  
    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[16, c(2,3, 6, 7, 8)]
    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    # Adding results to data frames
    results_smri_t1_crp <- rbind(results_smri_t1_crp, interaction_model_output)
  }
}

# T2
for (outcome in smri_outcomes) {
  for (exposure in exposures_t2_crp) {
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + gestage_plasma_g2 + (1| IDC)')  
    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[16, c(2,3, 6, 7, 8)]
    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    # Adding results to data frames
    results_smri_t2_crp <- rbind(results_smri_t2_crp, interaction_model_output)
  }
}

# Adding column names to created dataframes 
colnames(results_smri_t1_crp) <- c("beta", "SE" , "pval","lowerCI", "upperCI", "Outcome", 'Exposure')
colnames(results_smri_t2_crp) <- c("beta", "SE" , "pval","lowerCI", "upperCI", "Outcome", 'Exposure')

# Clean up tables to 3 decimals
results_smri_t1_crp[, 1:5] <- round(results_smri_t1_crp[, 1:5], digits = 3)
results_smri_t2_crp[, 1:5] <- round(results_smri_t2_crp[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
write_xlsx(results_smri_t1_crp, 'results_smri_t1_crp.xlsx')
write_xlsx(results_smri_t2_crp, 'results_smri_t2_crp.xlsx')

#\

# ========================================================== #
##                          DTI 
# ========================================================== #

## Create vector of outcomes
dti_outcomes <- c('mean_MD_genr', 'mean_FA_genr', 'avg_cgc_FA', 'avg_cgc_MD', 'avg_cst_FA', 'avg_cst_MD', 'avg_unc_FA', 'avg_unc_MD', 'avg_ilf_FA', 'avg_ilf_MD' , 'avg_slf_FA', 'avg_slf_MD')

## Create empty dataframes
results_dti_main_mia <- data.frame()
results_dti_main_crp <- data.frame()

results_dti_int_mia <- data.frame()
results_dti_int_crp <- data.frame()

results_dti_sex_mia <- data.frame()
results_dti_sex_crp <- data.frame()

results_dti_sens_mia <- data.frame()
results_dti_sens_crp <- data.frame()

results_dti_t1_mia <- data.frame()
results_dti_t1_crp <- data.frame()

results_dti_t2_mia <- data.frame()
results_dti_t2_crp <- data.frame()

### Maternal immune activation

## Apply all models and loop over average exposures and outcomes 
for (outcome in dti_outcomes) {
  
  for (exposure in exposures_mia) {
    
    # Model formulas
    a <- paste0('scale(',outcome,') ~ scale(',exposure,') + age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average +  BatchDay_average + (1| IDC)')  # main model
    
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average +  BatchDay_average + (1| IDC)')  # interaction model

    d <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri*GENDER + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + INCOME + PARITY + GestagePLg_average +  BatchDay_average + (1| IDC)')  # interaction model + sex
    
    e <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average +  BatchDay_average  + pm + im + GESTBIR + BMI_0 + (1| IDC)')  # interaction model + sens covars
    
    # Calculating beta, confidence interval, and p-value for all models
    main_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(a), weights = ipw_weights))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
    
    interaction_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[17, c(2,3, 6, 7, 8)]
    
    sex_interaction_output <- summary(pool(with(df_final_imputed, lmer(as.formula(d), weights = ipw_weights))), conf.int = TRUE)[20, c(2,3, 6, 7, 8)]
    
    sens_covars_interaction_output <- summary(pool(with(df_final_imputed, lmer(as.formula(e), weights = ipw_weights))), conf.int = TRUE)[21, c(2, 3,6, 7, 8)]
    
    # Adding outcome and exposure information
    main_model_output$outcome <- outcome
    main_model_output$exposure <- exposure
    
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    
    sex_interaction_output$outcome <- outcome
    sex_interaction_output$exposure <- exposure
    
    sens_covars_interaction_output$outcome <- outcome
    sens_covars_interaction_output$exposure <- exposure
    
    # Adding results to data frames
    results_dti_main_mia <- rbind(results_dti_main_mia, main_model_output)
    results_dti_int_mia <- rbind(results_dti_int_mia, interaction_model_output)
    results_dti_sex_mia <- rbind(results_dti_sex_mia, sex_interaction_output)
    results_dti_sens_mia <- rbind(results_dti_sens_mia, sens_covars_interaction_output)
  }
}

# Adding column names to created dataframes 
colnames(results_dti_main_mia) <- c("beta", "SE" , "pval","lowerCI", "upperCI", "Outcome", 'Exposure')
colnames(results_dti_int_mia) <- c("beta", "SE" ,"pval","lowerCI", "upperCI", "Outcome", 'Exposure')
colnames(results_dti_sex_mia) <- c("beta", "SE", "pval","lowerCI", "upperCI","Outcome", 'Exposure')
colnames(results_dti_sens_mia) <- c("beta", "SE","pval" ,"lowerCI", "upperCI","Outcome", 'Exposure')

# Clean up tables to 3 decimals
results_dti_main_mia[, 1:5] <- round(results_dti_main_mia[, 1:5], digits = 3)
results_dti_int_mia[, 1:5] <- round(results_dti_int_mia[, 1:5], digits = 3)
results_dti_sex_mia[, 1:5] <- round(results_dti_sex_mia[, 1:5], digits = 3)
results_dti_sens_mia[, 1:5] <- round(results_dti_sens_mia[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
results_dti_main_mia_fdr <- add_fdr_pvalue(results_dti_main_mia, 3)
write_xlsx(results_dti_main_mia_fdr, 'results_dti_main_mia.xlsx')
write_xlsx(results_dti_int_mia, 'results_dti_int_mia.xlsx')

results_dti_sex_mia_fdr <- add_fdr_pvalue(results_dti_sex_mia, 3)
write_xlsx(results_dti_sex_mia_fdr, 'results_dti_sex_mia.xlsx')
write_xlsx(results_dti_sens_mia, 'results_dti_sens_mia.xlsx')

## Apply all models and loop over time specific exposures and outcomes
# T1
for (outcome in dti_outcomes) {
  for (exposure in exposures_t1_mia) {
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg1 +  BatchDayg1 + (1| IDC)')  
    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[17, c(2,3, 6, 7, 8)]
    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    # Adding results to data frames
    results_dti_t1_mia <- rbind(results_dti_t1_mia, interaction_model_output)
  }
}

# T2
for (outcome in dti_outcomes) {
  for (exposure in exposures_t2_mia) {
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg2 +  BatchDayg2 + (1| IDC)')  
    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[17, c(2,3, 6, 7, 8)]
    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    # Adding results to data frames
    results_dti_t2_mia <- rbind(results_dti_t2_mia, interaction_model_output)
  }
}

# Adding column names to created dataframes 
colnames(results_dti_t1_mia) <- c("beta", "SE" , "pval","lowerCI", "upperCI", "Outcome", 'Exposure')
colnames(results_dti_t2_mia) <- c("beta", "SE" , "pval","lowerCI", "upperCI", "Outcome", 'Exposure')

# Clean up tables to 3 decimals
results_dti_t1_mia[, 1:5] <- round(results_dti_t1_mia[, 1:5], digits = 3)
results_dti_t2_mia[, 1:5] <- round(results_dti_t2_mia[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
write_xlsx(results_dti_t1_mia, 'results_dti_t1_mia.xlsx')
results_dti_t2_mia_fdr <- add_fdr_pvalue(results_dti_t2_mia, 3)
write_xlsx(results_dti_t2_mia_fdr, 'results_dti_t2_mia.xlsx')

# ========================================================== #

### C-reactive protein 

## Apply all models and loop over average exposures and outcomes 
for (outcome in dti_outcomes) {
  
  for (exposure in exposures_crp) {
    
    # Model formulas
    a <- paste0('scale(',outcome,') ~ scale(',exposure,') + age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average + (1| IDC)')  # main model
    
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average + (1| IDC)')  # interaction model
    
    d <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri*GENDER + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + INCOME + PARITY + GestagePLg_average + (1| IDC)')  # interaction model + sex
    
    e <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average + pm + im + GESTBIR + BMI_0 + (1| IDC)')  # interaction model + sens covars
    
    # Calculating beta, confidence interval, and p-value for all models
    main_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(a), weights = ipw_weights))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
    
    interaction_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[16, c(2,3, 6, 7, 8)]
    
    sex_interaction_output <- summary(pool(with(df_final_imputed, lmer(as.formula(d), weights = ipw_weights))), conf.int = TRUE)[19, c(2,3, 6, 7, 8)]
    
    sens_covars_interaction_output <- summary(pool(with(df_final_imputed, lmer(as.formula(e), weights = ipw_weights))), conf.int = TRUE)[20, c(2, 3,6, 7, 8)]
    
    # Adding outcome and exposure information
    main_model_output$outcome <- outcome
    main_model_output$exposure <- exposure
    
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure

    sex_interaction_output$outcome <- outcome
    sex_interaction_output$exposure <- exposure
    
    sens_covars_interaction_output$outcome <- outcome
    sens_covars_interaction_output$exposure <- exposure
    
    # Adding results to data frames
    results_dti_main_crp <- rbind(results_dti_main_crp, main_model_output)
    results_dti_int_crp <- rbind(results_dti_int_crp, interaction_model_output)
    results_dti_sex_crp <- rbind(results_dti_sex_crp, sex_interaction_output)
    results_dti_sens_crp <- rbind(results_dti_sens_crp, sens_covars_interaction_output)
  }
}

# Adding column names to created dataframes 
colnames(results_dti_main_crp) <- c("beta", "SE" , "pval","lowerCI", "upperCI", "Outcome", 'Exposure')
colnames(results_dti_int_crp) <- c("beta", "SE" ,"pval","lowerCI", "upperCI", "Outcome", 'Exposure')
colnames(results_dti_sex_crp) <- c("beta", "SE", "pval","lowerCI", "upperCI","Outcome", 'Exposure')
colnames(results_dti_sens_crp) <- c("beta", "SE","pval" ,"lowerCI", "upperCI","Outcome", 'Exposure')

# Clean up tables to 3 decimals
results_dti_main_crp[, 1:5] <- round(results_dti_main_crp[, 1:5], digits = 3)
results_dti_int_crp[, 1:5] <- round(results_dti_int_crp[, 1:5], digits = 3)
results_dti_sex_crp[, 1:5] <- round(results_dti_sex_crp[, 1:5], digits = 3)
results_dti_sens_crp[, 1:5] <- round(results_dti_sens_crp[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
write_xlsx(results_dti_main_crp, 'results_dti_main_crp.xlsx')
write_xlsx(results_dti_int_crp, 'results_dti_int_crp.xlsx')
write_xlsx(results_dti_sex_crp, 'results_dti_sex_crp.xlsx')
write_xlsx(results_dti_sens_crp, 'results_dti_sens_crp.xlsx')

## Apply all models and loop over time specific exposures and outcomes
# T1
for (outcome in dti_outcomes) {
  for (exposure in exposures_t1_crp) {
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + gestage_plasma_g1 + (1| IDC)')  
    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[16, c(2,3, 6, 7, 8)]
    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    # Adding results to data frames
    results_dti_t1_crp <- rbind(results_dti_t1_crp, interaction_model_output)
  }
}

# T2
for (outcome in dti_outcomes) {
  for (exposure in exposures_t2_crp) {
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + gestage_plasma_g2 + (1| IDC)')  
    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[16, c(2,3, 6, 7, 8)]
    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    # Adding results to data frames
    results_dti_t2_crp <- rbind(results_dti_t2_crp, interaction_model_output)
  }
}

# Adding column names to created dataframes 
colnames(results_dti_t1_crp) <- c("beta", "SE" , "pval","lowerCI", "upperCI", "Outcome", 'Exposure')
colnames(results_dti_t2_crp) <- c("beta", "SE" , "pval","lowerCI", "upperCI", "Outcome", 'Exposure')

# Clean up tables to 3 decimals
results_dti_t1_crp[, 1:5] <- round(results_dti_t1_crp[, 1:5], digits = 3)
results_dti_t2_crp[, 1:5] <- round(results_dti_t2_crp[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
write_xlsx(results_dti_t1_crp, 'results_dti_t1_crp.xlsx')
write_xlsx(results_dti_t2_crp, 'results_dti_t2_crp.xlsx')

#\

# ========================================================== #
##                     functional MRI 
# ========================================================== #

## Create vector of outcomes
fmri_outcomes <- c('cpl', 'ge', 'mod', 'wiNetwork_none', 'wiNetwork_default', 'wiNetwork_ParietoOccip', 'wiNetwork_FrontoParietal', 'wiNetwork_Salience', 'wiNetwork_CinguloOperc', 'wiNetwork_MedialParietal', 'wiNetwork_DorsalAttn', 'wiNetwork_VentralAttn', 'wiNetwork_Visual', 'wiNetwork_SMhand', 'wiNetwork_SMmouth', 'wiNetwork_Auditory', 'btNetwork_none', 'btNetwork_default', 'btNetwork_ParietoOccip', 'btNetwork_FrontoParietal', 'btNetwork_Salience', 'btNetwork_CinguloOperc', 'btNetwork_MedialParietal', 'btNetwork_DorsalAttn', 'btNetwork_VentralAttn', 'btNetwork_Visual', 'btNetwork_SMhand', 'btNetwork_SMmouth', 'btNetwork_Auditory')

## Create empty dataframes
results_fmri_main_mia <- data.frame()
results_fmri_main_crp <- data.frame()

results_fmri_int_mia <- data.frame()
results_fmri_int_crp <- data.frame()

results_fmri_sex_mia <- data.frame()
results_fmri_sex_crp <- data.frame()

results_fmri_sens_mia <- data.frame()
results_fmri_sens_crp <- data.frame()

results_fmri_t1_mia <- data.frame()
results_fmri_t1_crp <- data.frame()

results_fmri_t2_mia <- data.frame()
results_fmri_t2_crp <- data.frame()

### Maternal immune activation

## Apply all models and loop over average exposures and outcomes 
for (outcome in fmri_outcomes) {
  
  for (exposure in exposures_mia) {
    
    # Model formulas
    a <- paste0('scale(',outcome,') ~ scale(',exposure,') + age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average +  BatchDay_average + (1| IDC)')  # main model
    
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average +  BatchDay_average + (1| IDC)')  # interaction model
    
    d <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri*GENDER + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + INCOME + PARITY + GestagePLg_average +  BatchDay_average + (1| IDC)')  # interaction model + sex
    
    e <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average +  BatchDay_average  + pm + im + GESTBIR + BMI_0 + (1| IDC)')  # interaction model + sens covars
    
    # Calculating beta, confidence interval, and p-value for all models
    main_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(a), weights = ipw_weights))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
    
    interaction_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[17, c(2,3, 6, 7, 8)]
    
    sex_interaction_output <- summary(pool(with(df_final_imputed, lmer(as.formula(d), weights = ipw_weights))), conf.int = TRUE)[20, c(2,3, 6, 7, 8)]
    
    sens_covars_interaction_output <- summary(pool(with(df_final_imputed, lmer(as.formula(e), weights = ipw_weights))), conf.int = TRUE)[21, c(2, 3,6, 7, 8)]
    
    # Adding outcome and exposure information
    main_model_output$outcome <- outcome
    main_model_output$exposure <- exposure
    
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    
    sex_interaction_output$outcome <- outcome
    sex_interaction_output$exposure <- exposure
    
    sens_covars_interaction_output$outcome <- outcome
    sens_covars_interaction_output$exposure <- exposure
    
    # Adding results to data frames
    results_fmri_main_mia <- rbind(results_fmri_main_mia, main_model_output)
    results_fmri_int_mia <- rbind(results_fmri_int_mia, interaction_model_output)
    results_fmri_sex_mia <- rbind(results_fmri_sex_mia, sex_interaction_output)
    results_fmri_sens_mia <- rbind(results_fmri_sens_mia, sens_covars_interaction_output)
  }
}

# Adding column names to created dataframes 
colnames(results_fmri_main_mia) <- c("beta", "SE" , "pval","lowerCI", "upperCI", "Outcome", 'Exposure')
colnames(results_fmri_int_mia) <- c("beta", "SE" ,"pval","lowerCI", "upperCI", "Outcome", 'Exposure')
colnames(results_fmri_sex_mia) <- c("beta", "SE", "pval","lowerCI", "upperCI","Outcome", 'Exposure')
colnames(results_fmri_sens_mia) <- c("beta", "SE","pval" ,"lowerCI", "upperCI","Outcome", 'Exposure')

# Clean up tables to 3 decimals
results_fmri_main_mia[, 1:5] <- round(results_fmri_main_mia[, 1:5], digits = 3)
results_fmri_int_mia[, 1:5] <- round(results_fmri_int_mia[, 1:5], digits = 3)
results_fmri_sex_mia[, 1:5] <- round(results_fmri_sex_mia[, 1:5], digits = 3)
results_fmri_sens_mia[, 1:5] <- round(results_fmri_sens_mia[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
results_fmri_main_mia_fdr <- add_fdr_pvalue(results_fmri_main_mia, 3)
write_xlsx(results_fmri_main_mia_fdr, 'results_fmri_main_mia.xlsx')

results_fmri_int_mia_fdr <- add_fdr_pvalue(results_fmri_int_mia, 3)
write_xlsx(results_fmri_int_mia_fdr, 'results_fmri_int_mia.xlsx')

results_fmri_sex_mia_fdr <- add_fdr_pvalue(results_fmri_sex_mia, 3)
write_xlsx(results_fmri_sex_mia_fdr, 'results_fmri_sex_mia.xlsx')

results_fmri_sens_mia_fdr <- add_fdr_pvalue(results_fmri_sens_mia, 3)
write_xlsx(results_fmri_sens_mia_fdr, 'results_fmri_sens_mia.xlsx')

## Apply all models and loop over time specific exposures and outcomes
# T1
for (outcome in fmri_outcomes) {
  for (exposure in exposures_t1_mia) {
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg1 +  BatchDayg1 + (1| IDC)')  
    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[17, c(2,3, 6, 7, 8)]
    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    # Adding results to data frames
    results_fmri_t1_mia <- rbind(results_fmri_t1_mia, interaction_model_output)
  }
}

# T2
for (outcome in fmri_outcomes) {
  for (exposure in exposures_t2_mia) {
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg2 +  BatchDayg2 + (1| IDC)')  
    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[17, c(2,3, 6, 7, 8)]
    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    # Adding results to data frames
    results_fmri_t2_mia <- rbind(results_fmri_t2_mia, interaction_model_output)
  }
}

# Adding column names to created dataframes 
colnames(results_fmri_t1_mia) <- c("beta", "SE" , "pval","lowerCI", "upperCI", "Outcome", 'Exposure')
colnames(results_fmri_t2_mia) <- c("beta", "SE" , "pval","lowerCI", "upperCI", "Outcome", 'Exposure')

# Clean up tables to 3 decimals
results_fmri_t1_mia[, 1:5] <- round(results_fmri_t1_mia[, 1:5], digits = 3)
results_fmri_t2_mia[, 1:5] <- round(results_fmri_t2_mia[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
write_xlsx(results_fmri_t1_mia, 'results_fmri_t1_mia.xlsx')
results_fmri_t2_mia_fdr <- add_fdr_pvalue(results_fmri_t2_mia, 3)
write_xlsx(results_fmri_t2_mia_fdr, 'results_fmri_t2_mia.xlsx')

# ========================================================== #

### C-reactive protein 

## Apply all models and loop over average exposures and outcomes 
for (outcome in fmri_outcomes) {
  
  for (exposure in exposures_crp) {
    
    # Model formulas
    a <- paste0('scale(',outcome,') ~ scale(',exposure,') + age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average + (1| IDC)')  # main model
    
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average + (1| IDC)')  # interaction model
    
    d <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri*GENDER + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + INCOME + PARITY + GestagePLg_average + (1| IDC)')  # interaction model + sex
    
    e <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average + pm + im + GESTBIR + BMI_0 + (1| IDC)')  # interaction model + sens covars
    
    # Calculating beta, confidence interval, and p-value for all models
    main_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(a), weights = ipw_weights))), conf.int = TRUE)[2, c(2, 3, 6, 7, 8)]
    
    interaction_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[16, c(2,3, 6, 7, 8)]
    
    sex_interaction_output <- summary(pool(with(df_final_imputed, lmer(as.formula(d), weights = ipw_weights))), conf.int = TRUE)[19, c(2,3, 6, 7, 8)]
    
    sens_covars_interaction_output <- summary(pool(with(df_final_imputed, lmer(as.formula(e), weights = ipw_weights))), conf.int = TRUE)[20, c(2, 3,6, 7, 8)]
    
    # Adding outcome and exposure information
    main_model_output$outcome <- outcome
    main_model_output$exposure <- exposure
    
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    
    sex_interaction_output$outcome <- outcome
    sex_interaction_output$exposure <- exposure
    
    sens_covars_interaction_output$outcome <- outcome
    sens_covars_interaction_output$exposure <- exposure
    
    # Adding results to data frames
    results_fmri_main_crp <- rbind(results_fmri_main_crp, main_model_output)
    results_fmri_int_crp <- rbind(results_fmri_int_crp, interaction_model_output)
    results_fmri_sex_crp <- rbind(results_fmri_sex_crp, sex_interaction_output)
    results_fmri_sens_crp <- rbind(results_fmri_sens_crp, sens_covars_interaction_output)
  }
}

# Adding column names to created dataframes 
colnames(results_fmri_main_crp) <- c("beta", "SE" , "pval","lowerCI", "upperCI", "Outcome", 'Exposure')
colnames(results_fmri_int_crp) <- c("beta", "SE" ,"pval","lowerCI", "upperCI", "Outcome", 'Exposure')
colnames(results_fmri_sex_crp) <- c("beta", "SE", "pval","lowerCI", "upperCI","Outcome", 'Exposure')
colnames(results_fmri_sens_crp) <- c("beta", "SE","pval" ,"lowerCI", "upperCI","Outcome", 'Exposure')

# Clean up tables to 3 decimals
results_fmri_main_crp[, 1:5] <- round(results_fmri_main_crp[, 1:5], digits = 3)
results_fmri_int_crp[, 1:5] <- round(results_fmri_int_crp[, 1:5], digits = 3)
results_fmri_sex_crp[, 1:5] <- round(results_fmri_sex_crp[, 1:5], digits = 3)
results_fmri_sens_crp[, 1:5] <- round(results_fmri_sens_crp[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
results_fmri_main_crp_fdr <- add_fdr_pvalue(results_fmri_main_crp, 3)
write_xlsx(results_fmri_main_crp_fdr, 'results_fmri_main_crp.xlsx')
write_xlsx(results_fmri_int_crp, 'results_fmri_int_crp.xlsx')
write_xlsx(results_fmri_sex_crp, 'results_fmri_sex_crp.xlsx')
write_xlsx(results_fmri_sens_crp, 'results_fmri_sens_crp.xlsx')

## Apply all models and loop over time specific exposures and outcomes
# T1
for (outcome in fmri_outcomes) {
  for (exposure in exposures_t1_crp) {
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + gestage_plasma_g1 + (1| IDC)')  
    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[16, c(2,3, 6, 7, 8)]
    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    # Adding results to data frames
    results_fmri_t1_crp <- rbind(results_fmri_t1_crp, interaction_model_output)
  }
}

# T2
for (outcome in fmri_outcomes) {
  for (exposure in exposures_t2_crp) {
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + gestage_plasma_g2 + (1| IDC)')  
    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[16, c(2,3, 6, 7, 8)]
    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    # Adding results to data frames
    results_fmri_t2_crp <- rbind(results_fmri_t2_crp, interaction_model_output)
  }
}

# Adding column names to created dataframes 
colnames(results_fmri_t1_crp) <- c("beta", "SE" , "pval","lowerCI", "upperCI", "Outcome", 'Exposure')
colnames(results_fmri_t2_crp) <- c("beta", "SE" , "pval","lowerCI", "upperCI", "Outcome", 'Exposure')

# Clean up tables to 3 decimals
results_fmri_t1_crp[, 1:5] <- round(results_fmri_t1_crp[, 1:5], digits = 3)
results_fmri_t2_crp[, 1:5] <- round(results_fmri_t2_crp[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
write_xlsx(results_fmri_t1_crp, 'results_fmri_t1_crp.xlsx')
write_xlsx(results_fmri_t2_crp, 'results_fmri_t2_crp.xlsx')

#\

# ========================================================== #
##                individual cytokines - posthoc
# ========================================================== #

## Calculate mean cytokine value across trimesters
# transform mids to lon
long <- complete(df_final_imputed, action = "long", include =TRUE)

# calculate mean
cytokines <- c('IL1ß', 'IL6', 'IL23', 'IFN..', 'IL17A')

long[, paste0(cytokines, '_average')] <- lapply(cytokines, function(col_prefix) {
  cols <- grep(paste0(col_prefix, '_ObsConcg'), names(long), value = TRUE)
  rowMeans(long[, cols, drop = FALSE], na.rm = TRUE)
})

long <- rename(long, 'IFNy_average' = IFN.._average)

# transform back to mids 
df_final_imputed2 <- as.mids(long)

## Perform individual mixed-effects regressions for all cytokines
# Specify exposure/outcomes of interest 
all_cytokines <- c('IL1ß_average', 'IL6_average', 'IFNy_average', 'IL23_average', 'IL17A_average')
all_outcomes <- c(smri_outcomes, dti_outcomes, fmri_outcomes)

# Specify empty dataframe
results_individual_cytokines_interaction_model <- data.frame()

# Specify loop 
for (outcome in all_outcomes) {
  
  for (exposure in all_cytokines) {
    
    # Model formula for interaction model 
    int_model <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average + (1| IDC)') 
    
    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(df_final_imputed2, lmer(as.formula(int_model), weights = ipw_weights))), conf.int = TRUE)[16, c(2,3, 6, 7, 8)]

    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure

    # Adding results to data frames
    results_individual_cytokines_interaction_model <- rbind(results_individual_cytokines_interaction_model, interaction_model_output)
  }
}

# Adding column names to created dataframes 
colnames(results_individual_cytokines_interaction_model) <- c("beta", "SE" , "pval","lowerCI", "upperCI", "Outcome", 'Exposure')

# Clean up tables to 3 decimals
results_individual_cytokines_interaction_model[, 1:5] <- round(results_individual_cytokines_interaction_model[, 1:5], digits = 3)

# Apply FDR correction if needed a
results_individual_cytokines_interaction_model_fdr <- add_fdr_pvalue(results_individual_cytokines_interaction_model, 3)

# Create cytokine specific dataframes
posthoc_il6 <- subset(results_individual_cytokines_interaction_model_fdr, Exposure == 'IL6_average')
posthoc_il23 <- subset(results_individual_cytokines_interaction_model_fdr, Exposure == 'IL23_average')
posthoc_17a <- subset(results_individual_cytokines_interaction_model_fdr, Exposure == 'IL17A_average')
posthoc_ifny <- subset(results_individual_cytokines_interaction_model_fdr, Exposure == 'IFNy_average')
posthoc_ilb <- subset(results_individual_cytokines_interaction_model_fdr, Exposure == 'IL1ß_average')

# Reorder columns
posthoc_il6 <- dplyr::select(posthoc_il6, c('beta', 'lowerCI', 'upperCI', 'pval', 'SE', 'fdr_pvalue', 'Outcome', 'Exposure'))
posthoc_il23 <- dplyr::select(posthoc_il23, c('beta', 'lowerCI', 'upperCI', 'pval', 'SE', 'fdr_pvalue', 'Outcome', 'Exposure'))
posthoc_17a <- dplyr::select(posthoc_17a, c('beta', 'lowerCI', 'upperCI', 'pval', 'SE', 'fdr_pvalue', 'Outcome', 'Exposure'))
posthoc_ifny <- dplyr::select(posthoc_ifny, c('beta', 'lowerCI', 'upperCI', 'pval', 'SE', 'fdr_pvalue', 'Outcome', 'Exposure'))
posthoc_ilb <- dplyr::select(posthoc_ilb, c('beta', 'lowerCI', 'upperCI', 'pval', 'SE', 'fdr_pvalue', 'Outcome', 'Exposure'))

# Save dataframes to excel  
write_xlsx(posthoc_il6, 'posthoc_il6.xlsx')
write_xlsx(posthoc_il23, 'posthoc_il23.xlsx')
write_xlsx(posthoc_17a, 'posthoc_il17a.xlsx')
write_xlsx(posthoc_ifny, 'posthoc_ifn_y.xlsx')
write_xlsx(posthoc_ilb, 'posthoc_il1_b.xlsx')

#\ 

#------------------------------------------------#
####-------------Step 2: Figures-------------####
#-----------------------------------------------#

## Visualize interaction plots before FDR
# Select last imputed dataframe (with best convergence) for figures 
df_figure <- complete(df_final_imputed, 30)

# Rename levels within child sex 
df_figure$GENDER <- factor(df_figure$GENDER, levels = c("boy", "girl"), labels = c("Male", "Female"))

# ========================================================== #
##                        sMRI 
# ========================================================== #

# Specify models 
lm <- lmer(Thalamus_Proper_vol_subcortical ~ hs_crp_average*age_child_mri*GENDER + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + INCOME + PARITY + gestage_plasma_average + (1| IDC), weights = ipw_weights, data = df_figure)

lm2 <- lmer(Lateral_Ventricle_vol_subcortical ~ hs_crp_average*age_child_mri*GENDER + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + INCOME + PARITY + gestage_plasma_average + (1| IDC), weights = ipw_weights, data = df_figure)

# Create longitudinal plots 
a <- plot_model(lm, type = "pred", terms = c("age_child_mri", "hs_crp_average [1.263, 2.104, 2.849]" ,"GENDER")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Thalamus volume (mm3)") +
  ggplot2::ggtitle("") + 
  ggplot2::labs(color = 'Log2 transformed CRP') +
  scale_color_manual(values = c("1.263" = "purple", "2.104" = "orange", "2.849" = "green")) + #plotting 1st q, median, 3rd q 
  theme(strip.text = element_text(face = "bold")) + theme_bw()

b <- plot_model(lm2, type = "pred", terms = c("age_child_mri", "hs_crp_average [1.263, 2.104, 2.849]" ,"GENDER")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Lateral ventricles volume (mm3)") +
  ggplot2::ggtitle("") + 
  ggplot2::labs(color = 'Log2 transformed CRP') +
  scale_color_manual(values = c("1.263" = "purple", "2.104" = "orange", "2.849" = "green")) + #plotting 1st q, median, 3rd q 
  theme(strip.text = element_text(face = "bold")) + theme_bw()

# Create anatomy plots (three way interaction models)
suggestive_sign_brain_findings <- data.frame(
  region = c('lateral ventricle', 'thalamus proper'),
  threeway_interaction_beta = c(-0.017, -0.033)
)

c <- suggestive_sign_brain_findings %>%
  ggplot() +
  geom_brain(atlas = aseg,
             aes(fill = threeway_interaction_beta)) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()   
  ) + 
  theme(
    legend.position = 'bottom',
    legend.text = element_text(size = 8, angle = 25)
  ) +
  scale_fill_gradient(low = 'aquamarine', high = 'aquamarine4') + 
  guides(fill = guide_colorbar(barheight = 1)) + 
  labs(fill = 'Standardized three-way interaction beta coefficient')

# combine plots
combined <- ggarrange(a, b, ncol = 2)

ggarrange(c, combined, labels = 'auto', ncol = 1, nrow = 2)

# ========================================================== #
##                          DTI 
# ========================================================== #

# Specify models 
lm3 <- lmer(avg_cgc_MD ~mia_index_average +age_child_mri + GENDER + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + INCOME + PARITY + BatchDay_average + GestagePLg_average + (1| IDC), weights = ipw_weights, data = df_figure)
lm4 <- lmer(avg_ilf_MD ~mia_index_average*age_child_mri*GENDER + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + INCOME + PARITY + BatchDay_average + GestagePLg_average + (1| IDC), weights = ipw_weights, data = df_figure)
lm5 <- lmer(avg_ilf_MD ~ mia_index_t2*age_child_mri + GENDER + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + INCOME + PARITY + BatchDayg2 + GestagePLg2 + (1| IDC), weights = ipw_weights, data = df_figure)

# Create plots 
c <- plot_model(lm3, type = "pred", terms = c("mia_index_average")) + ggplot2::xlab("Log2 transformed MIA index") +
  ggplot2::ylab("Cingulum tract (MD)") +
  ggplot2::ggtitle("") + theme_bw() 

d <- plot_model(lm4, type = "pred", terms = c("age_child_mri", "mia_index_average[-0.579, 0.0188, 0.5918]" ,"GENDER")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Inferior longitudinal fasciculus (MD)") +
  ggplot2::ggtitle("") + 
  ggplot2::labs(color = 'Log2 transformed MIA index') +
  scale_color_manual(values = c("-0.579" = "purple", "0.0188" = "orange", "0.5918" = "green")) + #plotting 1st q, median, 3rd q 
  theme(strip.text = element_text(face = "bold")) + theme_bw()

e <- plot_model(lm5, type = "pred", terms = c("age_child_mri", "mia_index_t2[-0.6081, 0.6234]")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Inferior longitudinal fasciculus (MD)") +
  ggplot2::ggtitle("") + 
  ggplot2::labs(color = 'Log2 transformed MIA index (late preg)') +
  scale_color_manual(values = c("-0.6081" = "purple", "0.6234" = "orange"))+ theme_bw() #plotting 1st quartile and 3rd quartile 

ggarrange(c,e,d, labels = 'auto', ncol = 2, nrow = 3)

# ========================================================== #
##                         fMRI 
# ========================================================== #

# Define models
lm6 <- lmer(wiNetwork_VentralAttn ~ mia_index_average*age_child_mri + GENDER + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + INCOME + PARITY + BatchDay_average + GestagePLg_average + (1| IDC), weights = ipw_weights, data = df_figure)
lm7 <- lmer(wiNetwork_Visual ~ mia_index_average*age_child_mri*GENDER + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + INCOME + PARITY + BatchDay_average + GestagePLg_average + (1| IDC), weights = ipw_weights, data = df_figure)
lm8 <- lmer(wiNetwork_SMhand ~ mia_index_average*age_child_mri*GENDER + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + INCOME + PARITY + BatchDay_average + GestagePLg_average + (1| IDC), weights = ipw_weights, data = df_figure)
lm9 <- lmer(wiNetwork_Auditory ~ mia_index_average*age_child_mri*GENDER + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + INCOME + PARITY + BatchDay_average + GestagePLg_average + (1| IDC), weights = ipw_weights, data = df_figure)
lm10 <- lmer(wiNetwork_VentralAttn ~ mia_index_t2*age_child_mri + GENDER + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + INCOME + PARITY + BatchDayg2 + GestagePLg2 + (1| IDC), weights = ipw_weights, data = df_figure)
lm11 <- lmer(wiNetwork_MedialParietal ~ hs_crp_average + age_child_mri + GENDER + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + INCOME + PARITY + gestage_plasma_average + (1| IDC), weights = ipw_weights, data = df_figure)
lm12 <- lmer(btNetwork_SMhand ~ mia_index_average*age_child_mri*GENDER + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + INCOME + PARITY + BatchDay_average + GestagePLg_average + (1| IDC), weights = ipw_weights, data = df_figure)

# Create plots interaction + main
f <- plot_model(lm6, type = "pred", terms = c("age_child_mri", "mia_index_average[-0.5790, 0.5918]")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Within network functional connectivity - ventral attention") +
  ggplot2::ggtitle("") + 
  ggplot2::labs(color = 'Log2 transformed MIA index') +
  scale_color_manual(values = c("-0.579" = "purple", "0.5918" = "orange"))+ theme_bw() #plotting 1st quartile and 3rd quartile 

g <- plot_model(lm10, type = "pred", terms = c("age_child_mri", "mia_index_t2[-0.6081, 0.6234]")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Within network functional connectivity - ventral attention") +
  ggplot2::ggtitle("") + 
  ggplot2::labs(color = 'Log2 transformed MIA index (late preg)') +
  scale_color_manual(values = c("-0.6081" = "purple", "0.6234" = "orange"))+ theme_bw() #plotting 1st quartile and 3rd quartile 

h <- plot_model(lm11, type = "pred", terms = c("hs_crp_average")) + ggplot2::xlab("Log2 transformed CRP") +
  ggplot2::ylab("Within network functional connectivity - medial parietal") +
  ggplot2::ggtitle("") + theme_bw() 

ggarrange(f,g,h, labels = 'auto', ncol = 2, nrow = 2)

# Create plots sex three way interaction
i <- plot_model(lm7, type = "pred", terms = c("age_child_mri", "mia_index_average[-0.579, 0.0188, 0.5918]" ,"GENDER")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Within network FC - visual") +
  ggplot2::ggtitle("") + 
  ggplot2::labs(color = 'Log2 transformed MIA index') +
  scale_color_manual(values = c("-0.579" = "purple", "0.0188" = "orange", "0.5918" = "green")) + #plotting min, median, max 
  theme(strip.text = element_text(face = "bold")) + theme_bw()

j <- plot_model(lm8, type = "pred", terms = c("age_child_mri", "mia_index_average[-0.579, 0.0188, 0.5918]" ,"GENDER")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Within network FC - somatomotor hand") +
  ggplot2::ggtitle("") + 
  ggplot2::labs(color = 'Log2 transformed MIA index') +
  scale_color_manual(values = c("-0.579" = "purple", "0.0188" = "orange", "0.5918" = "green")) + #plotting min, median, max 
  theme(strip.text = element_text(face = "bold")) + theme_bw()

k <- plot_model(lm9, type = "pred", terms = c("age_child_mri", "mia_index_average[-0.579, 0.0188, 0.5918]" ,"GENDER")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Within network FC - auditory") +
  ggplot2::ggtitle("") + 
  ggplot2::labs(color = 'Log2 transformed MIA index') +
  scale_color_manual(values = c("-0.579" = "purple", "0.0188" = "orange", "0.5918" = "green")) + #plotting min, median, max 
  theme(strip.text = element_text(face = "bold")) + theme_bw()

l <- plot_model(lm12, type = "pred", terms = c("age_child_mri", "mia_index_average[-0.579, 0.0188, 0.5918]" ,"GENDER")) + ggplot2::xlab("Child Age (years)") +
  ggplot2::ylab("Between network FC - somatomotor hand") +
  ggplot2::ggtitle("") + 
  ggplot2::labs(color = 'Log2 transformed MIA index') +
  scale_color_manual(values = c("-0.579" = "purple", "0.0188" = "orange", "0.5918" = "green")) + #plotting min, median, max 
  theme(strip.text = element_text(face = "bold")) + theme_bw()

ggarrange(i, j, k, l, labels = 'auto', ncol = 2, nrow = 2)

#\

# ========================================================== #
##         Extra sensitivity analyses for revisions
# ========================================================== #

# ========================================================== #
### Sensitivity analysis 1: lateralization effects global brain regions

## Create vector of exposures 
exposures_mia <- c('mia_index_average')
exposures_crp <- c('hs_crp_average')

## Create vector of outcomes
df <- complete(df_final_imputed, 30)
smri_outcomes_lat <- colnames(dplyr::select(df, starts_with("lh"), starts_with("rh"), starts_with("Left"), starts_with("Right")))

## Create empty dataframes
lat_results_smri_int_mia <- data.frame()
lat_results_smri_int_crp <- data.frame()

## Cytokine index 
for (outcome in smri_outcomes_lat) {
  
  for (exposure in exposures_mia) {
    
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average +  BatchDay_average + (1| IDC)')  # interaction model
    
    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[17, c(2,3, 6, 7, 8)]

    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure

    # Adding results to data frames
    lat_results_smri_int_mia <- rbind(lat_results_smri_int_mia, interaction_model_output)
  }
}

# Adding column names to created dataframes 
colnames(lat_results_smri_int_mia) <- c("beta", "SE" ,"pval","lowerCI", "upperCI", "Outcome", 'Exposure')

# Clean up tables to 3 decimals
lat_results_smri_int_mia[, 1:5] <- round(lat_results_smri_int_mia[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
lat_results_smri_int_mia <- add_fdr_pvalue(lat_results_smri_int_mia, 3)
write_xlsx(lat_results_smri_int_mia, 'lat_results_smri_int_mia.xlsx')

# Create table in word
#library(rempsyc)
#my_table <- nice_table(lat_results_smri_int_mia2)
#print(my_table, preview = "docx")

## C-reactive protein 
for (outcome in smri_outcomes_lat) {
  
  for (exposure in exposures_crp) {
    
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + gestage_plasma_average + (1| IDC)')  # interaction model

    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[16, c(2,3, 6, 7, 8)]
    
    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    
    # Adding results to data frames
    lat_results_smri_int_crp <- rbind(lat_results_smri_int_crp, interaction_model_output)
 
  }
}

# Adding column names to created dataframes 
colnames(lat_results_smri_int_crp) <- c("beta", "SE" ,"pval","lowerCI", "upperCI", "Outcome", 'Exposure')

# Clean up tables to 3 decimals
lat_results_smri_int_crp[, 1:5] <- round(lat_results_smri_int_crp[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
lat_results_smri_int_crp <- add_fdr_pvalue(lat_results_smri_int_crp, 3)
write_xlsx(lat_results_smri_int_crp, 'lat_results_smri_int_crp.xlsx')

#\

# ========================================================== #
### Sensitivity analysis 2: rerun main analysis, additionally adjusting for puberty
combined_main_outcomes <- c(smri_outcomes, dti_outcomes, fmri_outcomes)

## Create empty dataframes
puberty_results_smri_int_mia <- data.frame()
puberty_results_smri_int_crp <- data.frame()

## Cytokine index 
for (outcome in combined_main_outcomes) {
  for (exposure in exposures_mia) {
  
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average +  BatchDay_average + puberty_score + (1| IDC)')  
    
    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[18, c(2, 3, 6, 7, 8)]
    
    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    
    # Adding results to data frames
    puberty_results_smri_int_mia <- rbind(puberty_results_smri_int_mia, interaction_model_output)
  }
}

# Adding column names to created dataframes 
colnames(puberty_results_smri_int_mia) <- c("beta", "SE", "pval","lowerCI", "upperCI", "Outcome", 'Exposure')

# Clean up tables to 3 decimals
puberty_results_smri_int_mia[, 1:5] <- round(puberty_results_smri_int_mia[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
puberty_results_smri_int_mia <- add_fdr_pvalue(puberty_results_smri_int_mia, 3)
write_xlsx(puberty_results_smri_int_mia, 'puberty_results_smri_int_mia.xlsx')

## CRP
for (outcome in combined_main_outcomes) {
  for (exposure in exposures_crp) {
    
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + gestage_plasma_average + puberty_score + (1| IDC)')  
    
    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(df_final_imputed, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[17, c(2, 3, 6, 7, 8)]
    
    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    
    # Adding results to data frames
    puberty_results_smri_int_crp <- rbind(puberty_results_smri_int_crp, interaction_model_output)
  }
}

# Adding column names to created dataframes 
colnames(puberty_results_smri_int_crp) <- c("beta", "SE", "pval","lowerCI", "upperCI", "Outcome", 'Exposure')

# Clean up tables to 3 decimals
puberty_results_smri_int_crp[, 1:5] <- round(puberty_results_smri_int_crp[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
puberty_results_smri_int_crp <- add_fdr_pvalue(puberty_results_smri_int_crp, 3)
write_xlsx(puberty_results_smri_int_crp, 'puberty_results_smri_int_crp.xlsx')

#\

# ========================================================== #
### Sensitivity analysis 3: rerun main analyses with more narrow time interval for cytokines

## Create new exposure indices 
# Transform mids object to long dataframe
imp_long <- complete(df_final_imputed, action = 'long', include = T)

# Select 4 week time frame during which most data is collected 
# df <- complete(df_final_imputed, 30)
# hist(df$gestage_plasma_g1)
# hist(df$gestage_plasma_g2)

imp_long_t1 <- subset(imp_long, gestage_plasma_g1 > 10 & gestage_plasma_g1 < 14)
imp_long_t2 <- subset(imp_long, gestage_plasma_g2 > 18 & gestage_plasma_g2 < 22)

# Convert back to mids 
df_final_imputed_t1_narrow <- as.mids(imp_long_t1)
df_final_imputed_t2_narrow <- as.mids(imp_long_t2)

## Run analysis for cytokines - T1
# Specify exposure
exposures_t1_cytokines <- c('mia_index_t1')

# Specify empty dataframe 
narrow_results_cytokines_t1 <- data.frame()

# Run model
for (outcome in combined_main_outcomes) {
  for (exposure in exposures_t1_cytokines) {
    
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg1 +  BatchDayg1 + (1| IDC)')  
    
    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(df_final_imputed_t1_narrow, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[17, c(2, 3, 6, 7, 8)]
    
    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    
    # Adding results to data frames
    narrow_results_cytokines_t1 <- rbind(narrow_results_cytokines_t1, interaction_model_output)
  }
}

# Adding column names to created dataframes 
colnames(narrow_results_cytokines_t1) <- c("beta", "SE", "pval","lowerCI", "upperCI", "Outcome", 'Exposure')

# Clean up tables to 3 decimals
narrow_results_cytokines_t1[, 1:5] <- round(narrow_results_cytokines_t1[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
narrow_results_cytokines_t1 <- add_fdr_pvalue(narrow_results_cytokines_t1, 3)
write_xlsx(narrow_results_cytokines_t1, 'narrow_results_cytokines_t1.xlsx')

## Run analysis for cytokines - T2
# Specify exposure
exposures_t2_cytokines <- c('mia_index_t2')

# Specify empty dataframe 
narrow_results_cytokines_t2 <- data.frame()

# Run model
for (outcome in combined_main_outcomes) {
  for (exposure in exposures_t2_cytokines) {
    
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg2 +  BatchDayg2 + (1| IDC)')  
    
    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(df_final_imputed_t2_narrow, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[17, c(2, 3, 6, 7, 8)]
    
    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    
    # Adding results to data frames
    narrow_results_cytokines_t2 <- rbind(narrow_results_cytokines_t2, interaction_model_output)
  }
}

# Adding column names to created dataframes 
colnames(narrow_results_cytokines_t2) <- c("beta", "SE", "pval","lowerCI", "upperCI", "Outcome", 'Exposure')

# Clean up tables to 3 decimals
narrow_results_cytokines_t2[, 1:5] <- round(narrow_results_cytokines_t2[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
narrow_results_cytokines_t2 <- add_fdr_pvalue(narrow_results_cytokines_t2, 3)
write_xlsx(narrow_results_cytokines_t2, 'narrow_results_cytokines_t2.xlsx')

## Run analysis for CRP - T1 (error ... )
# Specify exposure
exposures_t1_crp <- c('HsCRPmgL_g1')

# Specify empty dataframe 
narrow_results_crp_t1 <- data.frame()

# Run model
for (outcome in combined_main_outcomes) {
  for (exposure in exposures_t1_crp) {
    
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg1 + (1| IDC)')  
    
    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(df_final_imputed_t1_narrow, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[16, c(2, 3, 6, 7, 8)]
    
    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    
    # Adding results to data frames
    narrow_results_crp_t1 <- rbind(narrow_results_crp_t1, interaction_model_output)
  }
}

# Adding column names to created dataframes 
colnames(narrow_results_crp_t1) <- c("beta", "SE", "pval","lowerCI", "upperCI", "Outcome", 'Exposure')

# Clean up tables to 3 decimals
narrow_results_crp_t1[, 1:5] <- round(narrow_results_crp_t1[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
narrow_results_crp_t1 <- add_fdr_pvalue(narrow_results_crp_t1, 3)
write_xlsx(narrow_results_crp_t1, 'narrow_results_crp_t1.xlsx')

## Run analysis for CRP - T2
# Specify exposure
exposures_t2_crp <- c('HsCRPmgL_g2')

# Specify empty dataframe 
narrow_results_crp_t2 <- data.frame()

# Run model
for (outcome in combined_main_outcomes) {
  for (exposure in exposures_t2_crp) {
    
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg2 + (1| IDC)')  
    
    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(df_final_imputed_t2_narrow, lmer(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[16, c(2, 3, 6, 7, 8)]
    
    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    
    # Adding results to data frames
    narrow_results_crp_t2 <- rbind(narrow_results_crp_t2, interaction_model_output)
  }
}

# Adding column names to created dataframes 
colnames(narrow_results_crp_t2) <- c("beta", "SE", "pval","lowerCI", "upperCI", "Outcome", 'Exposure')

# Clean up tables to 3 decimals
narrow_results_crp_t2[, 1:5] <- round(narrow_results_crp_t2[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
narrow_results_crp_t2 <- add_fdr_pvalue(narrow_results_crp_t2, 3)
write_xlsx(narrow_results_crp_t2, 'narrow_results_crp_t2.xlsx')

# ========================================================== #
### Sensitivity analysis 4: rerun main analysis with linear regressions at F9 and F13

### Transform from long to wide format for linear regressions
# Transform mids object to long dataframe
imp_long <- complete(df_final_imputed, action = 'long', include = T)
imp_long$timepoint <- as.factor(imp_long$timepoint)

# Create timepoint specific dataframes 
imp_long_f09 <- subset(imp_long, timepoint == 'f09')
imp_long_f13 <- subset(imp_long, timepoint == 'f13')

# Transform back to mids
imp_df_f09 <- as.mids(imp_long_f09)
imp_df_f13 <- as.mids(imp_long_f13)

### Run linear regressions for f09
# Create vector of exposures 
exposures_mia <- c('mia_index_average')
exposures_crp <- c('hs_crp_average')

# Create empty dataframes
lm_results_smri_mia_f09 <- data.frame()
lm_results_smri_crp_f09 <- data.frame()

## Cytokine index 
for (outcome in combined_main_outcomes) {
  
  for (exposure in exposures_mia) {
    
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,') + age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average +  BatchDay_average')  
    
    # Calculating beta, confidence interval, and p-value for all models
    model <- summary(pool(with(imp_df_f09, lm(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[2, c(2,3, 6, 7, 8)]
    
    # Adding outcome and exposure information
    model$outcome <- outcome
    model$exposure <- exposure
    
    # Adding results to data frames
    lm_results_smri_mia_f09 <- rbind(lm_results_smri_mia_f09, model)
  }
}

# Adding column names to created dataframes 
colnames(lm_results_smri_mia_f09) <- c("beta", "SE" ,"pval","lowerCI", "upperCI", "Outcome", 'Exposure')

# Clean up tables to 3 decimals
lm_results_smri_mia_f09[, 1:5] <- round(lm_results_smri_mia_f09[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
lm_results_smri_mia_f09 <- add_fdr_pvalue(lm_results_smri_mia_f09, 3)
write_xlsx(lm_results_smri_mia_f09, 'lm_results_smri_mia_f09.xlsx')

## C-reactive protein 
for (outcome in combined_main_outcomes) {
  
  for (exposure in exposures_crp) {
    
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,') + age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + gestage_plasma_average')  
    
    # Calculating beta, confidence interval, and p-value for all models
    model <- summary(pool(with(imp_df_f09, lm(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[2, c(2,3, 6, 7, 8)]
    
    # Adding outcome and exposure information
    model$outcome <- outcome
    model$exposure <- exposure
    
    # Adding results to data frames
    lm_results_smri_crp_f09 <- rbind(lm_results_smri_crp_f09, model)
  }
}

# Adding column names to created dataframes 
colnames(lm_results_smri_crp_f09) <- c("beta", "SE" ,"pval","lowerCI", "upperCI", "Outcome", 'Exposure')

# Clean up tables to 3 decimals
lm_results_smri_crp_f09[, 1:5] <- round(lm_results_smri_crp_f09[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
lm_results_smri_crp_f09 <- add_fdr_pvalue(lm_results_smri_crp_f09, 3)
write_xlsx(lm_results_smri_crp_f09, 'lm_results_smri_crp_f09.xlsx')

#\ 

### Run linear regressions for f13
## Create empty dataframes
lm_results_smri_mia_f13 <- data.frame()
lm_results_smri_crp_f13 <- data.frame()

## Cytokine index 
for (outcome in combined_main_outcomes) {
  
  for (exposure in exposures_mia) {
    
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,') + age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average +  BatchDay_average')  
    
    # Calculating beta, confidence interval, and p-value for all models
    model <- summary(pool(with(imp_df_f13, lm(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[2, c(2,3, 6, 7, 8)]
    
    # Adding outcome and exposure information
    model$outcome <- outcome
    model$exposure <- exposure
    
    # Adding results to data frames
    lm_results_smri_mia_f13 <- rbind(lm_results_smri_mia_f13, model)
  }
}

# Adding column names to created dataframes 
colnames(lm_results_smri_mia_f13) <- c("beta", "SE" ,"pval","lowerCI", "upperCI", "Outcome", 'Exposure')

# Clean up tables to 3 decimals
lm_results_smri_mia_f13[, 1:5] <- round(lm_results_smri_mia_f13[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
lm_results_smri_mia_f13 <- add_fdr_pvalue(lm_results_smri_mia_f13, 3)
write_xlsx(lm_results_smri_mia_f13, 'lm_results_smri_mia_f13.xlsx')

## C-reactive protein 
for (outcome in combined_main_outcomes) {
  
  for (exposure in exposures_crp) {
    
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,') + age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + gestage_plasma_average')  
    
    # Calculating beta, confidence interval, and p-value for all models
    model <- summary(pool(with(imp_df_f13, lm(as.formula(b), weights = ipw_weights))), conf.int = TRUE)[2, c(2,3, 6, 7, 8)]
    
    # Adding outcome and exposure information
    model$outcome <- outcome
    model$exposure <- exposure
    
    # Adding results to data frames
    lm_results_smri_crp_f13 <- rbind(lm_results_smri_crp_f13, model)
  }
}

# Adding column names to created dataframes 
colnames(lm_results_smri_crp_f13) <- c("beta", "SE" ,"pval","lowerCI", "upperCI", "Outcome", 'Exposure')

# Clean up tables to 3 decimals
lm_results_smri_crp_f13[, 1:5] <- round(lm_results_smri_crp_f13[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
lm_results_smri_crp_f13 <- add_fdr_pvalue(lm_results_smri_crp_f13, 3)
write_xlsx(lm_results_smri_crp_f13, 'lm_results_smri_crp_f13.xlsx')

# ========================================================== #
### Sensitivity analysis 5: rerun main analysis with between Fc for all outcomes 
# Load data 
between_fc_f09 <- readRDS("fc_fMRI_f09.rds")
col_names <- names(between_fc_f09)
new_col_names <- c(col_names[1], paste0(col_names[-1], "_f09"))
names(between_fc_f09) <- new_col_names
#print(between_fc_f09)

between_fc_f13 <- readRDS("fc_fMRI_f13.rds")
col_names2 <- names(between_fc_f13)
new_col_names2 <- c(col_names2[1], paste0(col_names2[-1], "_f13"))
names(between_fc_f13) <- new_col_names2
#print(between_fc_f13)

between_fc_combined <- merge(between_fc_f09, between_fc_f13, by = 'IDC', all = T)

df_long <- between_fc_combined %>% pivot_longer(-c(IDC), names_to = c('variable', 'timepoint'), values_to = c('value'), names_pattern = "(.*)_(.*)")
df_long_ordered <- pivot_wider(df_long, names_from = c("variable"), values_from = c("value"))

names(df_long_ordered) <- gsub("-", "_", names(df_long_ordered))
#print(names(df_long_ordered))

# Merge with other data
# Transform to long datasets for MRI (smri + dti)
imp.data_combined <- list()

for (i in 0:30) {
  imp.data <- complete(df_final_imputed, i)
  reshape_imp.data <- merge(imp.data, df_long_ordered, by = 'IDC', all.x = T)
  imp.data_combined[[i+1]] <- reshape_imp.data
}

reshape_imp_combined <- datalist2mids(imp.data_combined)

# Define outcomes
bt_fmri_outcomes <- colnames(dplyr::select(df_long_ordered, -c('IDC', 'timepoint')))

# Define exposures
exposures_mia <- c('mia_index_average')
exposures_crp <- c('hs_crp_average')

# Create empty dataframes
fmri_extra_mia <- data.frame()
fmri_extra_crp <- data.frame()

# Run analysis for MIA
for (outcome in bt_fmri_outcomes) {
  
  for (exposure in exposures_mia) {
    
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + GestagePLg_average + BatchDay_average + (1| IDC)')  
    
    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(reshape_imp_combined, lmer(as.formula(b)))), conf.int = TRUE)[17, c(2,3, 6, 7, 8)]
    
    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    
    # Adding results to data frames
    fmri_extra_mia <- rbind(fmri_extra_mia, interaction_model_output)
  }
}

# Adding column names to created dataframes 
colnames(fmri_extra_mia) <- c("beta", "SE" ,"pval","lowerCI", "upperCI", "Outcome", 'Exposure')

# Clean up tables to 3 decimals
fmri_extra_mia[, 1:5] <- round(fmri_extra_mia[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
fmri_extra_mia <- add_fdr_pvalue(fmri_extra_mia, 3)
write_xlsx(fmri_extra_mia, 'fmri_extra_mia.xlsx')

# Create table in word
#library(rempsyc)
#my_table <- nice_table(fmri_extra_mia)
#print(my_table, preview = "docx")

# Run analysis for CRP
for (outcome in bt_fmri_outcomes) {
  
  for (exposure in exposures_crp) {
    
    # Model formulas
    b <- paste0('scale(',outcome,') ~ scale(',exposure,')*age_child_mri + AGE_M_v2 + ETHNMv2 + mdrink_updated + SMOKE_ALL + GSI + GENDER + INCOME + PARITY + gestage_plasma_average + (1| IDC)')  
    
    # Calculating beta, confidence interval, and p-value for all models
    interaction_model_output <- summary(pool(with(reshape_imp_combined, lmer(as.formula(b)))), conf.int = TRUE)[16, c(2,3, 6, 7, 8)]
    
    # Adding outcome and exposure information
    interaction_model_output$outcome <- outcome
    interaction_model_output$exposure <- exposure
    
    # Adding results to data frames
    fmri_extra_crp <- rbind(fmri_extra_crp, interaction_model_output)
  }
}

# Adding column names to created dataframes 
colnames(fmri_extra_crp) <- c("beta", "SE" ,"pval","lowerCI", "upperCI", "Outcome", 'Exposure')

# Clean up tables to 3 decimals
fmri_extra_crp[, 1:5] <- round(fmri_extra_crp[, 1:5], digits = 3)

# Apply FDR correction if needed and save data frames to Excel
fmri_extra_crp <- add_fdr_pvalue(fmri_extra_crp, 3)
write_xlsx(fmri_extra_crp, 'fmri_extra_crp.xlsx')

#\ END OF SCRIPT
