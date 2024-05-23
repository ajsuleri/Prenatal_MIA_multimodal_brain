#############################################################
##########Project: MIA & multimodal brain development########
#############################################################

## Author: Anna Suleri

## General note: In this script we will prepare the data for the analysis; we will select the data we need and create a final combined dataframe, make a baseline table and check the correlations between the variables, create descriptive figures for exposure and outcome, we will transform the dataframe to long format, we will then impute the data using mice and apply inverse probability weighting, then we will apply the inclusion and exclusion criteria, after which we will end with our final multiply imputed dataset for further analyses (script 2). 

#------------------------------------------------------------#
####---------------Step 1: Create dataframe---------------####
#------------------------------------------------------------#

## Clear environment and set seed
rm(list = ls()) 

set.seed(181223)

setwd("set_path_to_scripts")

source('0.Functions.R')

## Set working directory and load libraries
libraries <- c('foreign', 'haven', 'tidyverse' ,'mice', 'lattice', 'ggplot2', 'writexl', 'corrplot', 'stringi', 'miceadds', 'mitools', 'CBPS', 'survey', 'survival', 'psych', 'plotly', 'sjPlot', 'RColorBrewer', 'lme4', 'ggpubr')

invisible(lapply(libraries, require, character.only = T))

## Load data and select variables of interest
dataframes_folder <- 'path_to_data'

load_files(dataframes_folder, "\\.csv$", read.csv)
load_files(dataframes_folder, "\\.sav$", function(file) read.spss(file, to.data.frame = TRUE))
load_files(dataframes_folder, "\\.rds$", readRDS)

# ============================================================================ #

# Load covariate data
covars_m1 <- merge(Second_hits_DF[, c('IDC', 'IDM', 'MOTHER', 'AGE_M_v2', 'BMI_0', 'ETHNMv2', 'mdrink_updated', 'DIAB_GRA', 'PIH_v1', 'preeclampsia', 'SMOKE_ALL', 'GSI', 'GENDER', 'GESTBIR')], `CHILD-ALLGENERALDATA_24102022`[, c('IDC', 'EDUCM', 'INCOME', 'PARITY')], by = 'IDC', all.x = T)
covars_m2 <- merge(covars_m1, MEDICATIONSELFREPORTPREGNANCY_30112017[, c('IDM', 'SSRITOT', 'TRIPTOT', 'PSYTOT', 'TCATOT', 'NSAIDTOT', 'ABIOTOT', 'PMOLTOT', 'CORTTOT', 'MUCOTOT' ,'COUGHTOT')], by = 'IDM', all.x = T)
covars <- merge(covars_m2, `GR1001-D1-37_16072021`[, c('IDM','d2300101', 'd2400101', 'd2500101', 'd1100101')], by = 'IDM', all.x = T)

# Create binary scores for maternal inflammatory conditions, immune medication use and psychotropic medication use
covars$preg_ic <- as.factor(ifelse(covars$DIAB_GRA == 'Yes' | covars$preeclampsia == 'Yes' | covars$PIH_v == 'Yes', 1, ifelse(covars$DIAB_GRA == "No" & covars$preeclampsia == 'No' & covars$PIH_v1 == 'No', 0, NA)))

covars$ic <- as.factor(ifelse(covars$d2300101 == 'Yes' | covars$d2400101 == 'Yes' | covars$d2500101 == 'Yes' | covars$d1100101 == 'Yes', 1, ifelse(covars$d2300101 == 'No' & covars$d2400101 == 'No' & covars$d2500101 == 'No' & covars$d1100101 == 'No', 0, NA)))

covars$mat_inflam <- as.factor(ifelse(covars$preg_ic == 1 | covars$ic == 1, 1, ifelse(covars$preg_ic == 0 & covars$ic == 0, 0, NA)))

covars$im <- as.factor(ifelse(covars$NSAIDTOT == 'no use' & covars$ABIOTOT == 'no use' & covars$PMOLTOT == 'no use' & covars$CORTTOT == 'no use' & covars$MUCOTOT == 'no use' & covars$COUGHTOT == 'no use', 0, ifelse(is.na(covars$NSAIDTOT) | is.na(covars$ABIOTOT) | is.na(covars$PMOLTOT) | is.na(covars$CORTTOT) | is.na(covars$MUCOTOT) | is.na(covars$COUGHTOT), NA, 1)))

covars$pm <- as.factor(ifelse(covars$SSRITOT == 'no use' & covars$TRIPTOT == 'no use' & covars$PSYTOT == 'no use' & covars$TCATOT == 'no use', 0, ifelse(is.na(covars$SSRITOT) | is.na(covars$TRIPTOT) | is.na(covars$PSYTOT) | is.na(covars$PSYTOT), NA, 1)))

# Check structure of covars and relevel if needed 
str(covars)

covars$EDUCM <- as.factor(ifelse(covars$EDUCM == "no education finished" | covars$EDUCM == "primary", "Low", ifelse(covars$EDUCM == "secondary, phase 1" | covars$EDUCM == "secondary, phase 2", "Middle", ifelse(covars$EDUCM == "higher, phase 1" | covars$EDUCM =="higher, phase 2", "High", NA))))

covars$INCOME <- as.factor(ifelse(covars$INCOME == 'less than 450' | covars$INCOME == '450-600 euro' | covars$INCOME == '600-700 euro' | covars$INCOME == '700-800 euro' | covars$INCOME == '800-900 euro' | covars$INCOME == '900-1200 euro' | covars$INCOME == '1200-1400 euro' | covars$INCOME == '1400-1600 euro' | covars$INCOME == '1600-1800 euro' | covars$INCOME == '1800-2000 euro' | covars$INCOME == '2000-2200 euro', '<2000', ifelse(is.na(covars$INCOME), NA, '>2000')))

# Dump columns we don't need anymore 
covars <- dplyr::select(covars, -c('SSRITOT', 'TRIPTOT', 'PSYTOT', 'TCATOT', 'NSAIDTOT', 'ABIOTOT', 'PMOLTOT', 'CORTTOT', 'MUCOTOT' ,'COUGHTOT','d2300101', 'd2400101', 'd2500101', 'd1100101', 'DIAB_GRA', 'PIH_v1', 'preeclampsia', 'preg_ic', 'ic'))

# ============================================================================ #

# Select outcome data for f9 
core <- genr_mri_core_data_20220311[, -grep('f05', names(genr_mri_core_data_20220311))] 
smri.vars_f09 <- merge(f09_freesurfer_v6_09dec2016_aseg_stats_pull06june2017_v1[,c('idc','Left_Cerebellum_Cortex_vol_f09', 'Right_Cerebellum_Cortex_vol_f09', 'Left_Thalamus_Proper_vol_f09', 'Right_Thalamus_Proper_vol_f09', 'Left_Hippocampus_vol_f09', 'Right_Hippocampus_vol_f09', 'Left_Amygdala_vol_f09', 'Right_Amygdala_vol_f09', 'Left_Lateral_Ventricle_vol_f09', 'Right_Lateral_Ventricle_vol_f09')],f09_freesurfer_v6_09dec2016_tbv_stats_pull20june2017_v2[,c('idc', 'genr_tbv_f09', 'eTIV_f09', 'CortexVol_f09', 'TotalGrayVol_f09', 'CerebralWhiteMatterVol_f09')] ,by = 'idc', all.x = T)
smri.vars_f09 <- merge(smri.vars_f09, f09_freesurfer_v6_09dec2016_lobes_stats_pull06june2017, by = 'idc', all.x = T)
smri.vars_f09 <- dplyr::select(smri.vars_f09, -c(ends_with('_thickavg_f09'), ends_with("_surfarea_f09")))

dti.vars_f09 <- dplyr::select(f09_GenR_MRI_eddy_dipy_wls_14Feb2022_autoPtx_dti_stats_inc_glob_measV1, ends_with('_wavg_FA_f09'),ends_with('_wavg_MD_f09'),'mean_MD_genr_f09', 'mean_FA_genr_f09', 'missingness_genr_f09', 'idc')
dti.vars_f09 <- dplyr::select(dti.vars_f09, starts_with('cgc_'), starts_with('cst_'), starts_with('unc_'), starts_with('ilf_'), starts_with('slf_'), starts_with('fma_'), starts_with('fmi_'), 'mean_MD_genr_f09', 'mean_FA_genr_f09', 'missingness_genr_f09', 'idc')

fmri.vars_f09_wiN <- within_con_f9_correct[, c(1:14)]
fmri.vars_f09_gt <- bold_gift_sfnc_f09_graphtheory_20220615[, c('idc', 'cpl' ,'ge', 'mod')]

names <- names(fmri.vars_f09_gt)
names[names != "idc"] <- paste0(names[names != "idc"], "_f09")
names(fmri.vars_f09_gt) <- names

# calculate between connectivity 
networks <- c("None", "Default", "ParietoOccip", "FrontoParietal", "Salience", "CinguloOperc", "MedialParietal", "DorsalAttn", "VentralAttn", "Visual", "SMhand", "SMmouth", "Auditory")

result_list <- lapply(networks, calculate_btnetworks)
result_df <- do.call(dplyr::bind_cols, result_list)
colnames(result_df) <- c("idc", paste0("mean_", networks))
colnames(result_df) <- c("IDC","btNetwork_none_f09", "IDC" ,"btNetwork_default_f09","IDC","btNetwork_ParietoOccip_f09","IDC" ,"btNetwork_FrontoParietal_f09","IDC" ,"btNetwork_Salience_f09","IDC" ,"btNetwork_CinguloOperc_f09","IDC" ,"btNetwork_MedialParietal_f09","IDC" ,"btNetwork_DorsalAttn_f09","IDC" ,"btNetwork_VentralAttn_f09","IDC" ,"btNetwork_Visual_f09","IDC" ,"btNetwork_SMhand_f09","IDC" ,"btNetwork_SMmouth_f09","IDC" ,"btNetwork_Auditory_f09")
fmri.vars_f09_btN <- dplyr::select(result_df, -c(3, 5,7,9,11,13,15,17,19,21,23,25))

# Adapt the names including timepoint to make it easier to transform dataset to long format later 
fmri.vars_f09_wiN <- rename(fmri.vars_f09_wiN, "wiNetwork_none_f09" = None, 'wiNetwork_default_f09' = Default, 'wiNetwork_ParietoOccip_f09' = ParietoOccip, 'wiNetwork_FrontoParietal_f09' = FrontoParietal, 'wiNetwork_Salience_f09' = Salience, 'wiNetwork_CinguloOperc_f09' = CinguloOperc, 'wiNetwork_MedialParietal_f09' = MedialParietal, 'wiNetwork_DorsalAttn_f09' = DorsalAttn, 'wiNetwork_VentralAttn_f09' = VentralAttn, 'wiNetwork_Visual_f09' = Visual, 'wiNetwork_SMhand_f09' = SMhand, 'wiNetwork_SMmouth_f09' = SMmouth, 'wiNetwork_Auditory_f09' = Auditory) 

# Merge all brain vars at f09
brain.vars_f09 <- merge(core, smri.vars_f09, by = 'idc', all.x = T)
brain.vars_f09 <- merge(brain.vars_f09, dti.vars_f09, by = 'idc', all.x = T)
brain.vars_f09 <- merge(brain.vars_f09, fmri.vars_f09_wiN, by = 'idc', all.x = T)
brain.vars_f09 <- merge(brain.vars_f09, fmri.vars_f09_gt, by = 'idc', all.x = T)
brain.vars_f09 <- merge(brain.vars_f09, fmri.vars_f09_btN, by.x = 'idc', by.y = 'IDC', all.x = T)
brain.vars_f09 <- rename(brain.vars_f09, 'IDC' = idc)

# ============================================================================ #

# Select outcome data for f13
smri.vars_f13 <- merge(f13_freesurfer_v6_14oct2020_aseg_stats_pull23Nov2020_v1[,c('idc', 'Left_Cerebellum_White_Matter_vol_f13', 'Right_Cerebellum_Cortex_vol_f13', 'Left_Thalamus_Proper_vol_f13', 'Right_Thalamus_Proper_vol_f13', 'Left_Hippocampus_vol_f13', 'Right_Hippocampus_vol_f13', 'Left_Amygdala_vol_f13', 'Right_Amygdala_vol_f13', 'Left_Lateral_Ventricle_vol_f13', 'Right_Lateral_Ventricle_vol_f13')],f13_freesurfer_v6_14oct2020_tbv_stats_pull23Nov2020_v2[,c('idc', 'genr_tbv_f13', 'eTIV_f13', 'CortexVol_f13', 'TotalGrayVol_f13', 'CerebralWhiteMatterVol_f13')] ,by = 'idc', all.x = T)
smri.vars_f13 <- merge(smri.vars_f13, f13_freesurfer_v6_29may2021_lobes_stats_pull23Nov2020, by = 'idc', all.x = T)
smri.vars_f13 <- dplyr::select(smri.vars_f13, -c(ends_with('_thickavg_f13'), ends_with("_surfarea_f13")))

dti.vars_f13 <- dplyr::select(f13_GenR_MRI_eddy_dipy_wls_14Feb2022_autoPtx_dti_stats_inc_glob_meas, ends_with('_wavg_FA_f13'),ends_with('_wavg_MD_f13'),'mean_MD_genr_f13', 'mean_FA_genr_f13', 'missingness_genr_f13', 'idc')
dti.vars_f13 <- dplyr::select(dti.vars_f13, starts_with('cgc_'), starts_with('cst_'), starts_with('unc_'), starts_with('ilf_'), starts_with('slf_'), starts_with('fma_'), starts_with('fmi_'), 'mean_MD_genr_f13', 'mean_FA_genr_f13', 'missingness_genr_f13', 'idc')

fmri.vars_f13 <- merge(fmri_graphtheory_static_f13[, -which(names(fmri_graphtheory_static_f13) == "X")], fmri_fc_with_and_bt_updated_f13[, -which(names(fmri_fc_with_and_bt_updated_f13) == "X")], by = 'IDC', all.x = TRUE)

# Adapt the names including timepoint to make it easier to transform dataset to long format later 
new_column_names <- names(fmri.vars_f13)
new_column_names[new_column_names != "IDC"] <- paste0(new_column_names[new_column_names != "IDC"], "_f13")
names(fmri.vars_f13) <- new_column_names

# Merge all brain vars at f13
brain.vars_f13 <- merge(smri.vars_f13, dti.vars_f13, by = 'idc', all.x = T)
brain.vars_f13 <- rename(brain.vars_f13, 'IDC' = idc)
brain.vars_f13 <- merge(brain.vars_f13, fmri.vars_f13, by = 'IDC', all.x = T)

# ============================================================================ #

# Merge f9 and f13 brain data
brain.vars <- merge(brain.vars_f09, brain.vars_f13, by = 'IDC', all.x = T)

# Merge brain & covars data 
brain_covars <- merge(covars, brain.vars, by = 'IDC', all.x = T)

# Lateralize brain outcomes across hemispheres for sMRI 
aseg_cols <- colnames(dplyr::select(brain_covars, starts_with("Left_"), starts_with("Right_"))) 
tbv_cols <- colnames(dplyr::select(brain_covars, starts_with("lh"), starts_with('rh'))) 
aseg_structs1 <-  stri_remove_empty(unlist(strsplit(aseg_cols, "Left_"))) 
aseg_structs2 <-  stri_remove_empty(unlist(strsplit(aseg_structs1, "Right_")))
aseg_structs3 <- stri_remove_empty(unlist(strsplit(aseg_structs2, "_f09")))
aseg_structs4 <- stri_remove_empty(unlist(strsplit(aseg_structs3, "_f13")))

tbv_structs1 <- stri_remove_empty(unlist(strsplit(tbv_cols, "lh")))
tbv_structs2 <- stri_remove_empty(unlist(strsplit(tbv_structs1, "rh")))
tbv_structs3 <- stri_remove_empty(unlist(strsplit(tbv_structs2, "_f09")))
tbv_structs4 <- stri_remove_empty(unlist(strsplit(tbv_structs3, "_f13")))

structs <- c(aseg_structs4, tbv_structs4) 

brain_covars2 <- merge_smri_hemispheres(structs, brain_covars)

df_tidy <- dplyr::select(brain_covars2, -c(starts_with("Left_"), starts_with("Right_"), starts_with("lh"), starts_with("rh")))

# Lateralize brain outcomes across hemispheres for DTI 
fa_suffix <- c("FA_f09", "FA_f13")
md_suffix <- c("MD_f09", "MD_f13")
tracts <- c("cgc", "cst", "unc", "ilf", "slf")

for (tract in tracts) {
  for (suffix in fa_suffix) {
    cols <- grep(paste(tract, "_[lr]_dti_dipy_wls_wavg_", suffix, sep = ""), colnames(df_tidy), value = TRUE)
    avg_col <- paste("avg", tract, suffix, sep = "_")
    df_tidy[, avg_col] <- rowMeans(df_tidy[, cols], na.rm = TRUE)
  }
  
  for (suffix in md_suffix) {
    cols <- grep(paste(tract, "_[lr]_dti_dipy_wls_wavg_", suffix, sep = ""), colnames(df_tidy), value = TRUE)
    avg_col <- paste("avg", tract, suffix, sep = "_")
    df_tidy[, avg_col] <- rowMeans(df_tidy[, cols], na.rm = TRUE)
  }
}

df_tidy2 <- dplyr::select(df_tidy, -c(ends_with('_wavg_FA_f13'), ends_with('_wavg_MD_f13'), ends_with('_wavg_FA_f09'), ends_with('_wavg_MD_f09')))

# Clean up names + check if all vars are in the correct structure 
str(df_tidy2)

current_colnames <- colnames(df_tidy2)
new_col <- sub("^_", "", current_colnames)
colnames(df_tidy2) <- new_col

# ============================================================================ #

# Select exposure data
cytokine.vars <- dplyr::select(MOTHERPREGNANCY_Cytokinen_03082023, -c('DDMMG1', 'DDMMG2', ends_with('_FIg1'), ends_with('_FIg2')))
crp.vars <- dplyr::select(`MOTHERPREGNANCY-CRP_06062023`, -c('suffix_HsCRPmgL_g1','suffix_HsCRPmgL_g2'))

# Check structure of exposure data
str(cytokine.vars)
cytokine.vars2 <- mutate_at(cytokine.vars, vars(ends_with("_ObsConcg1") | ends_with("_ObsConcg2")), transform_to_numeric)

str(crp.vars)
crp.vars2 <- mutate_at(crp.vars, vars(ends_with("mgl_g1") | ends_with("mgl_g2")), transform_to_numeric)

# Remove bad quality exposure data for cytokines - empty tubes, low amount, low bead count & select only complete cases 
cytokine.vars_t1 <- subset(cytokine.vars2, complete.cases(IL23_ObsConcg1) & (!remark_cat_1 == 'empty tube - exclude' & !remark_cat_1 == 'low bead count or no beads - exclude')) #6087

cytokine.vars_t2 <- subset(cytokine.vars2, complete.cases(IL23_ObsConcg2) & (!remark_cat_2 == 'empty tube - exclude' & !remark_cat_2 == 'low bead count or no beads - exclude' & !remark_cat_2 == 'low amount - exclude'))  #7427

# Drop columns we don't need
cytokine.vars_t1 <- dplyr::select(cytokine.vars_t1, c('IDM', 'GestagePLg1', 'remark_cat_1', ends_with('_ObsConcg1'))) 

cytokine.vars_t2 <- dplyr::select(cytokine.vars_t2, c('IDM', 'GestagePLg2', 'remark_cat_2', ends_with('_ObsConcg2')))

# Log2 transform exposure data (of note, we add a constant because 1 participant has a zero value in IL-6)
crp.vars3 <- mutate_at(crp.vars2, vars(ends_with("mgl_g1") | ends_with("mgl_g2")), log2)

constant <- 0.1

cytokine.vars_t1 <- mutate_at(cytokine.vars_t1, vars(ends_with("_ObsConcg1")), function(x) {
  ifelse(x == 0, 0, log2(x + constant))
})

cytokine.vars_t2 <- mutate_at(cytokine.vars_t2, vars(ends_with("_ObsConcg2")), function(y) {
  ifelse(y == 0, 0, log2(y + constant))
})

# ============================================================================ #

# Now, we will create the maternal immune activation index based on cytokines 

# Select only cytokine vars for T1 and T2 
cvars.t1 <- dplyr::select(cytokine.vars_t1, ends_with('ObsConcg1'))

cvars.t2 <- dplyr::select(cytokine.vars_t2, ends_with('ObsConcg2'))

# Create correlation matrix for all cytokines per trimester 
cvars.t1 <- rename(cvars.t1, 'IL-1β' = IL1β_ObsConcg1, 'IL-6' = IL6_ObsConcg1, 'IL-17a' = IL17A_ObsConcg1, 'IFN-γ' = IFNγ_ObsConcg1, 'IL-23' = IL23_ObsConcg1)
cvars.t1_cor <- cor(cvars.t1) 
corrplot::corrplot(cvars.t1_cor, method = 'color', addCoef.col = "black", number.cex=0.6, order = 'FPC', type = 'lower', diag = F,  tl.col = 'black', tl.cex = 0.8, sig.level = 0.05, col = colorRampPalette(c("midnightblue","white","darkred"))(100))

cvars.t2 <- rename(cvars.t2, 'IL-1β' = IL1β_ObsConcg2, 'IL-6' = IL6_ObsConcg2, 'IL-17a' = IL17A_ObsConcg2, 'IFN-γ' = IFNγ_ObsConcg2, 'IL-23' = IL23_ObsConcg2)
cvars.t2_cor <- cor(cvars.t2) 
corrplot::corrplot(cvars.t2_cor, method = 'color',addCoef.col = "black",number.cex=0.6,order = 'FPC', type = 'lower',diag = F, tl.col = 'black',tl.cex = 0.8,  sig.level = 0.05, col = colorRampPalette(c("midnightblue","white","darkred"))(100))

# Create correlation matrix for all cytokines and for both trimesters
cor_df <- merge(cytokine.vars_t1, cytokine.vars_t2, by = 'IDM', all = T)
cor_df2 <- merge(cor_df, crp.vars3, by = 'IDM', all = T)
cor_df3 <- dplyr::select(cor_df2, ends_with('ObsConcg1'), ends_with('ObsConcg2'), 'HsCRPmgL_g1', 'HsCRPmgL_g2')
complete_cases_df <- na.omit(cor_df3)
complete_cases_df2 <- rename(complete_cases_df, 'IL-1β (early pregnancy)' = IL1β_ObsConcg1, 'IL-6 (early pregnancy)' = IL6_ObsConcg1, 'IL-17a (early pregnancy)' = IL17A_ObsConcg1, 'IFN-γ (early pregnancy)' = IFNγ_ObsConcg1, 'IL-23 (early pregnancy)' = IL23_ObsConcg1, 'IL-1β (mid pregnancy)' = IL1β_ObsConcg2, 'IL-6 (mid pregnancy)' = IL6_ObsConcg2, 'IL-17a (mid pregnancy)' = IL17A_ObsConcg2, 'IFN-γ (mid pregnancy)' = IFNγ_ObsConcg2, 'IL-23 (mid pregnancy)' = IL23_ObsConcg2, 'hs-CRP (early pregnancy)' = HsCRPmgL_g1, 'hs-CRP (mid pregnancy)' = HsCRPmgL_g2)
cor_df3.cor <- cor(complete_cases_df2)

corrplot::corrplot(cor_df3.cor, method = 'color',addCoef.col = "black",number.cex=0.6,order = 'FPC', type = 'lower',diag = F, tl.col = 'black',tl.cex = 0.8,  sig.level = 0.05, col = colorRampPalette(c("midnightblue","white","darkred"))(100))

# Calculate Kaiser-Meyer-Olkin factor adequacy (>0.5 = good)
psych::KMO(cvars.t1_cor) 

psych::KMO(cvars.t2_cor)

# Calculate multicolinearity between vars (>0.00001 means no problem)
det(cvars.t1_cor) 

det(cvars.t2_cor)

# Calculate principal component
cvars.t1.pca <- principal(cvars.t1, nfactors =1, rotate = 'none')

cvars.t2.pca <- principal(cvars.t2, nfactors =1, rotate = 'none')

# Elbow plot to confirm number of pc's based on point of inflexion (=1 in our case)
plot(cvars.t1.pca$values, type = 'b') 

plot(cvars.t2.pca$values, type = 'b')

# Extract pc1 and add to exposure dataframe
pc_value.t1 <- cvars.t1.pca$scores[, 1]
cvars.t1$index_t1 <- pc_value.t1
mia_index_t1 <- cvars.t1[, c('index_t1')]

pc_value.t2 <- cvars.t2.pca$scores[, 1]
cvars.t2$index_t2 <- pc_value.t2
mia_index_t2 <- cvars.t2[, c('index_t2')]

# Add ID to mia_index dataframes
id.t1 <- subset(cytokine.vars_t1, complete.cases(IL23_ObsConcg1)) 
id.t1 <- cbind(id.t1, mia_index_t1)

id.t2 <- subset(cytokine.vars_t2, complete.cases(IL23_ObsConcg2)) 
id.t2 <- cbind(id.t2, mia_index_t2)

# ============================================================================ #

# Merge data with covars and brain
exp.vars <- merge(id.t1, id.t2, by = c(1:19), all.y = T)
exp.vars2 <- merge(exp.vars, crp.vars3, by = 'IDM', all.y = T)

full_df <- merge(df_tidy2, exp.vars2, by = 'IDM', all.x = T)

# Merge batch effects covars with dataframe
full_df <- merge(full_df, cytokine.vars[, c('IDM', 'BatchDayg1', 'BatchDayg2', 'RackIDg1', 'RackIDg2')], by = 'IDM', all.x = T)

## Save final full dataframe
setwd("path_to_store_final_df")

saveRDS(full_df, 'full_df.rds')

#\

#------------------------------------------------------------------------#
####------------------Step 2: Create inclusion variable---------------####
#------------------------------------------------------------------------#

full_df <- readRDS("full_df.rds")

## Create inclusion variable based on inclusion and exclusion criteria 
# inclusion criterium 1 (complete case exposure at t1 or t2) (n=8223)
inclu1 <- subset(full_df, complete.cases(IL23_ObsConcg1) | complete.cases(IL23_ObsConcg2) | complete.cases(HsCRPmgL_g1) | complete.cases(HsCRPmgL_g2)) 

# inclusion criterium 2 (complete case one of the smri/dti/fmri outcome) (n=4294)
inclu2 <- subset(inclu1, complete.cases(genr_tbv_f09) | complete.cases(genr_tbv_f13) | complete.cases(mean_MD_genr_f09) | complete.cases(mean_MD_genr_f13) | complete.cases(wiNetwork_none_f09) | complete.cases(wiNetwork_none_f13))

# exclusion criterium 1 (exclude twin or sibling with least available data) (n=4022)
exclu1 <- inclu2[sample(nrow(inclu2)),]
exclu1$na_count <- apply(exclu1, 1, function(x) sum(is.na(x)))
exclu1 <- exclu1[order(exclu1$na_count),]
exclu1 <- exclu1[!duplicated(exclu1$MOTHER, fromLast = T),]

# exclusion criterium 2a (unusable scan, poor quality scan scan with IF or braces) - sMRI
exclu1$inclu_f09_smri <- as.factor(ifelse(exclu1$mri_consent_f09 == 'yes' & exclu1$t1_has_nii_f09 == 'yes' & exclu1$t1_asset_has_nii_f09 != 'exclude' & exclu1$has_braces_mri_f09 == 'no' &  exclu1$exclude_incidental_f09 == 'include' & exclu1$freesurfer_qc_f09 == 'usable', 'include', 'exclude'))
exclu1$inclu_f13_smri <- as.factor(ifelse(exclu1$mri_consent_f13 == 'yes' & exclu1$t1_has_nii_f13 == 'yes' & exclu1$has_braces_mri_f13 == 'no' & exclu1$exclude_incidental_f13 == 'include' &exclu1$freesurfer_qc_f13 == 'usable', 'include', 'exclude'))

df_f09_smri <- subset(exclu1, complete.cases(genr_tbv_f09) & inclu_f09_smri == 'include')
df_f13_smri <- subset(exclu1, complete.cases(genr_tbv_f13) & inclu_f13_smri == 'include')

# list of ID with all good quality sMRI children for f9 or f13 (n=3296; n=2513 at f09 and n=1744 at f13)
df_smri <- merge(df_f09_smri, df_f13_smri, by = c(1:212), all = T)
df_smri <- dplyr::select(df_smri, "IDC")

# exclusion criterium 2b (unusable scan, poor quality scan scan with IF or braces) - DTI
exclu1$inclu_f09_dti <- as.factor(ifelse(exclu1$mri_consent_f09 == 'yes' & exclu1$has_braces_mri_f09 == 'no' &  exclu1$exclude_incidental_f09 == 'include' &  exclu1$dti_man_qc_f09 == 'usable' & complete.cases(exclu1$missingness_genr_f09), 'include', 'exclude'))
exclu1$inclu_f13_dti <- as.factor(ifelse(exclu1$mri_consent_f13 == 'yes' & exclu1$has_braces_mri_f13 == 'no' &  exclu1$exclude_incidental_f13 == 'include' &  exclu1$dti_man_qc_f13 == 'usable' & complete.cases(exclu1$missingness_genr_f13), 'include', 'exclude'))

df_f09_dti <- subset(exclu1, complete.cases(mean_MD_genr_f09) & inclu_f09_dti == 'include')
df_f13_dti <- subset(exclu1, complete.cases(mean_MD_genr_f13) & inclu_f13_dti == 'include')

# list of ID with all good quality DTI children for f9 or f13 (n=3267; n=2359 at f09 and n=1878 at f13)
df_dti <- merge(df_f09_dti, df_f13_dti, by = c(1:214), all = T)
df_dti <- dplyr::select(df_dti, "IDC")

# exclusion criterium 2c (unusable scan, poor quality scan scan with IF or braces) - fMRI
exclu1$inclu_f09_fmri <- as.factor(ifelse(exclu1$mri_consent_f09 == 'yes' & exclu1$has_braces_mri_f09 == 'no' &  exclu1$exclude_incidental_f09 == 'include' & exclu1$rsfmri_has_nii_f09 == 'yes' & exclu1$num_vols_bold_f09 == 200 & exclu1$mean_bold_rms_f09 <= 0.25 & exclu1$exclude_bold_f09 == 'include' & is.na(exclu1$exclude_reg_prob_bold_f09), 'include', 'exclude'))
exclu1$inclu_f13_fmri <- as.factor(ifelse(exclu1$mri_consent_f13 == 'yes' & exclu1$has_braces_mri_f13 == 'no' &  exclu1$exclude_incidental_f13 == 'include' & exclu1$rsfmri_has_nii_f13 == 'yes' & exclu1$num_vols_bold_f13 == 200 & exclu1$mean_bold_rms_f13 <= 0.25 & exclu1$exclude_bold_f13 == 'include' & is.na(exclu1$exclude_reg_prob_bold_f13), 'include', 'exclude'))

df_f09_fmri <- subset(exclu1, complete.cases(wiNetwork_none_f09) & inclu_f09_fmri == 'include')
df_f13_fmri <- subset(exclu1, complete.cases(wiNetwork_none_f13) & inclu_f13_fmri == 'include')

# list of ID with all good quality fMRI children for f9 or f13 (n=2914; n=1954 at f09 and n=1697 at f13)
df_fmri <- merge(df_f09_fmri, df_f13_fmri, by = c(1:216), all = T)
df_fmri <- dplyr::select(df_fmri, "IDC")

# merge brain dataframes together with all exclusion criteria (n total=3648)
df_brain <- merge(df_smri, df_dti, by = 'IDC', all = T)
df_brain2 <- merge(df_brain, df_fmri, by = 'IDC', all = T)

# Create inclusion/exclusion variable in full dataset (full_df) based on the children who are left after applying above inclusion/exclusion criteria 
full_df$include <- as.factor(ifelse(full_df$IDC %in% df_brain2$IDC, 1, 0))

#\

#-----------------------------------------------------------------------#
####---------------Step 3: Create descriptive figures---------------####
#----------------------------------------------------------------------#

# Select study participants out of whole dataframe
df_descriptives <- subset(full_df, include == 1)

# ============================================================================ #

# Baseline table (for included sample; n=3648)
baselinevars <- c('AGE_M_v2', 'BMI_0', 'ETHNMv2', 'mdrink_updated', 'SMOKE_ALL', 'GSI', 'GENDER', 'GESTBIR', 'EDUCM', 'PARITY', 'INCOME', 'age_child_mri_f09', 'age_child_mri_f13')

for(i in baselinevars){ 
  #x = i vars that are columns in the dataframe df
  x <- df_descriptives[, i] 
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
    summary_continuous(x)
  } 
  else 
  {
    summary_categorical(x) 
  }
}

# ============================================================================ #

# Make a correlation plot (for included sample)
corr_vars <- dplyr::select(df_descriptives, 
                           c('AGE_M_v2', 'BMI_0', 'ETHNMv2', 'mdrink_updated', 'SMOKE_ALL', 
                             'GSI', 
                             'GENDER', 'GESTBIR', 'EDUCM', 'PARITY', 'INCOME', 'genr_tbv_f09', 
                             'genr_tbv_f13', 'mean_MD_genr_f09', 'mean_MD_genr_f13', 
                             'mean_FA_genr_f09'
                             , 'mean_FA_genr_f13', 'cpl_f09', 'cpl_f13', 'ge_f09', 'ge_f13', 
                             'mod_f09'
                             , 'mod_f13', 'mia_index_t1', 'mia_index_t2', 'IL1ß_ObsConcg1', 
                             'IL1ß_ObsConcg2', 'IL6_ObsConcg1', 'IL6_ObsConcg2', 
                             'IL17A_ObsConcg1', 
                             'IL17A_ObsConcg2', 'IFN.._ObsConcg1', 'IFN.._ObsConcg2', 
                             'IL23_ObsConcg1'
                             , 'IL23_ObsConcg2', 'HsCRPmgL_g1', 'HsCRPmgL_g2', 
                             'frontal_vol_cortical_f09', 'frontal_vol_cortical_f13', 
                             'cingulate_vol_cortical_f09', 'cingulate_vol_cortical_f13', 
                             'occipital_vol_cortical_f09', 'occipital_vol_cortical_f13', 
                             'temporal_vol_cortical_f09', 'temporal_vol_cortical_f13', 
                             'parietal_vol_cortical_f09', 'parietal_vol_cortical_f13', 
                             'insula_vol_cortical_f09', 'insula_vol_cortical_f13'))

corr_vars2 <- rename(corr_vars, 
                     'Maternal age' = AGE_M_v2, 
                     'Maternal pre-pregnancy BMI' = BMI_0, 
                     'Maternal national background' =ETHNMv2, 
                     'Maternal alcohol use' = mdrink_updated, 
                     'Maternal tobacco use' = SMOKE_ALL, 
                     'Maternal psychopathology' = GSI, 
                     'Child sex' = GENDER, 
                     'Gestational age at birth' = GESTBIR, 
                     'Maternal education' = EDUCM, 
                     'Parity' = PARITY, 
                     'Household income' = INCOME, 
                     'Total brain volume (T1)' = genr_tbv_f09, 
                     'Total brain volume (T2)' = genr_tbv_f13, 
                     'Global mean diffusivity (T1)' = mean_MD_genr_f09, 
                     'Global mean diffusivity (T2)' = mean_MD_genr_f13, 
                     'Global fractional anisotropy (T1)' = mean_FA_genr_f09, 
                     'Global fractional anisotropy (T2)' = mean_FA_genr_f13, 
                     'Characteristic path length (T1)' = cpl_f09, 
                     'Characteristic path length (T2)' = cpl_f13, 
                     'Global efficiency (T1)' = ge_f09, 
                     'Global efficiency (T2)' = ge_f13, 
                     'Modularity (T1)' = mod_f09, 
                     'Modularity (T2)' = mod_f13, 
                     'MIA index (T1)' = mia_index_t1, 
                     'MIA index (T2)' = mia_index_t2, 
                     'IL-1ß (early preg)' = IL1ß_ObsConcg1, 
                     'IL-1ß (late preg)' = IL1ß_ObsConcg2, 
                     'IL-6 (early preg)' = IL6_ObsConcg1, 
                     'IL-6 (late preg)' = IL6_ObsConcg2, 
                     'IL-17 (early preg)' = IL17A_ObsConcg1, 
                     'IL-17 (late preg)' = IL17A_ObsConcg2, 
                     'IFN-y (early preg)' = IFN.._ObsConcg1, 
                     'IFN-y (late preg)' = IFN.._ObsConcg2, 
                     'IL-23 (early preg)' = IL23_ObsConcg1, 
                     'IL-23 (late preg)' = IL23_ObsConcg2, 
                     'hs-CRP (early preg)' = HsCRPmgL_g1, 
                     'hs-CRP (late preg)' = HsCRPmgL_g2, 
                     'Frontal lobe volume (T1)' = frontal_vol_cortical_f09, 
                     'Frontal lobe volume (T2)' = frontal_vol_cortical_f13, 
                     'Cingulate volume (T1)' = cingulate_vol_cortical_f09, 
                     'Cingulate volume (T2)' = cingulate_vol_cortical_f13, 
                     'Occipital volume (T1)' = occipital_vol_cortical_f09, 
                     'Occipital volume (T2)' = occipital_vol_cortical_f13, 
                     'Temporal volume (T1)' = temporal_vol_cortical_f09, 
                     'Temporal volume (T2)' = temporal_vol_cortical_f13, 
                     'Parietal volume (T1)' = parietal_vol_cortical_f09, 
                     'Parietal volume (T2)' = parietal_vol_cortical_f13, 
                     'Insula volume (T1)' = insula_vol_cortical_f09, 
                     'Insule volume (T2)' = insula_vol_cortical_f13)  

corr_vars2[] <- lapply(corr_vars2, as.numeric)

correlation <- cor(corr_vars2, use="pairwise.complete.obs")

corrplot(correlation, method = 'color', addCoef.col = "black", number.cex=0.3, order = 'FPC', type = 'lower', diag = F, tl.col = 'black', tl.cex = 0.4, col = colorRampPalette(c("midnightblue","white","darkred"))(100),tl.srt = 45)

# ============================================================================ #

# Descriptive figure MIA index & CRP exposure
p1 <- ggplot(df_descriptives, aes(x = mia_index_t1)) +
  geom_histogram(binwidth = 1, fill = "maroon4", color = "white", alpha = 0.5) +
  labs(x = "Log2 transformed MIA index (early pregnancy)", y = "Frequency") + theme_minimal()
p2 <- ggplot(df_descriptives, aes(x = mia_index_t2)) +
  geom_histogram(binwidth = 1, fill = "orange", color = "white", alpha = 0.5) +
  labs(x = "Log2 transformed MIA index (late pregnancy)", y = "Frequency") + theme_minimal()
p3 <- ggplot(df_descriptives, aes(x = HsCRPmgL_g1)) +
  geom_histogram(binwidth = 1, fill = "maroon4", color = "white", alpha = 0.5) +
  labs(x = "Log2 transformed hs-CRP (early pregnancy)", y = "Frequency") + theme_minimal()
p4 <- ggplot(df_descriptives, aes(x = HsCRPmgL_g2)) +
  geom_histogram(binwidth = 1, fill = "orange", color = "white", alpha = 0.5) +
  labs(x = "Log2 transformed hs-CRP (late pregnancy)", y = "Frequency") + theme_minimal()

ggarrange(p1, p2, p3, p4, labels = 'auto')

df_exposure <- dplyr::select(df_descriptives, c('IDC', 'mia_index_t1', 'mia_index_t2', 'IL1ß_ObsConcg1', 'IL1ß_ObsConcg2', 'IL6_ObsConcg1', 'IL6_ObsConcg2', 'IL17A_ObsConcg1', 'IL17A_ObsConcg2', 'IFN.._ObsConcg1', 'IFN.._ObsConcg2', 'IL23_ObsConcg1', 'IL23_ObsConcg2', 'HsCRPmgL_g1', 'HsCRPmgL_g2',  'GestagePLg1', 'GestagePLg2'))
df_exposure <- rename(df_exposure, 'IFNy_ObsConcg1' = IFN.._ObsConcg1, 'IFNy_ObsConcg2' = IFN.._ObsConcg2)

df_exposure_long <- reshape(df_exposure, idvar = 'IDC', varying = list(c('GestagePLg1', 'GestagePLg2'), c("mia_index_t1", "mia_index_t2"), c('IL1ß_ObsConcg1', 'IL1ß_ObsConcg2'), c('IL6_ObsConcg1', 'IL6_ObsConcg2'), c('IL17A_ObsConcg1', 'IL17A_ObsConcg2'), c('IFNy_ObsConcg1', 'IFNy_ObsConcg2'), c('IL23_ObsConcg1', 'IL23_ObsConcg2'), c('HsCRPmgL_g1', 'HsCRPmgL_g2')), v.names = c("time", 'MIA_index', 'IL-1ß', 'IL-6', 'IL-17A', 'IFN-y', 'IL-23', 'hs_CRP'), direction ='long')

p5 <- ggplot(df_exposure_long, aes(x = time, y = MIA_index)) +
  geom_point(color = '#bc5090', size = 0.9, alpha = 0.2) +
  geom_smooth(method = "gam", se = F, color = '#ffa600', linetype = "solid", size = 0.7) +
  labs(x = "Time (gestational age in weeks)",
       y = "Log2 transformed MIA index") +
  theme_minimal()

p6 <- ggplot(df_exposure_long, aes(x = time, y = hs_CRP)) +
  geom_point(color = '#bc5090', size = 0.9, alpha = 0.2) +
  geom_smooth(method = "gam", se = F, color = '#ffa600', linetype = "solid", size = 0.7) +
  labs(x = "Time (gestational age in weeks)",
       y = "Log2 transformed hs-C-reactive protein") +
  theme_minimal()

ggarrange(p5, p6, labels = 'auto')

# ============================================================================ #

# Descriptive figure global sMRI/DTI/fMRI outcomes over time
df_outcome <- dplyr::select(df_descriptives, c('IDC', 'genr_tbv_f09', 'genr_tbv_f13', 'TotalGrayVol_f09', 'TotalGrayVol_f13', 'CerebralWhiteMatterVol_f09','CerebralWhiteMatterVol_f13', 'Cerebellum_Cortex_vol_subcortical_f09', 'Cerebellum_Cortex_vol_subcortical_f13','mean_MD_genr_f09', 'mean_MD_genr_f13', 'mean_FA_genr_f09', 'mean_FA_genr_f13','cpl_f09', 'cpl_f13', 'ge_f09', 'ge_f13', 'mod_f09', 'mod_f13', 'age_child_mri_f09', 'age_child_mri_f13'))

df_outcome_long <- reshape(df_outcome, 
                           idvar = 'IDC', 
                           varying = list(c('genr_tbv_f09', 'genr_tbv_f13'), 
                                          c("TotalGrayVol_f09", "TotalGrayVol_f13"), 
                                          c('CerebralWhiteMatterVol_f09', 
                                            'CerebralWhiteMatterVol_f13'), 
                                          c('Cerebellum_Cortex_vol_subcortical_f09', 
                                            'Cerebellum_Cortex_vol_subcortical_f13'), 
                                          c('mean_MD_genr_f09', 'mean_MD_genr_f13'),
                                          c('mean_FA_genr_f09', 'mean_FA_genr_f13'),
                                          c('cpl_f09', 'cpl_f13'), 
                                          c('ge_f09', 'ge_f13'),
                                          c('mod_f09', 'mod_f13'),
                                          c('age_child_mri_f09', 'age_child_mri_f13')), 
                           v.names = c('tbv', 'gm', 'wm', 'cc', 'md', 'fa', 'cpl', 'ge', 'mod'
                                       , 'time'), 
                           direction ='long')

b1 <- ggplot(df_outcome_long, aes(x = time, y = tbv)) +
  geom_point(color = '#bc5090', size = 0.9, alpha = 0.2) +
  geom_smooth(method = "gam", se = F, color = '#ffa600', linetype = "solid", size = 0.7) +
  labs(x = "Time (child age in years)",
       y = "Total brain volume (mm^3)") +
  theme_minimal()

b2 <- ggplot(df_outcome_long, aes(x = time, y = gm)) +
  geom_point(color = '#bc5090', size = 0.9, alpha = 0.2) +
  geom_smooth(method = "gam", se = F, color = '#ffa600', linetype = "solid", size = 0.7) +
  labs(x = "Time (child age in years)",
       y = "Gray matter volume (mm^3)") +
  theme_minimal()

b3 <- ggplot(df_outcome_long, aes(x = time, y = wm)) +
  geom_point(color = '#bc5090', size = 0.9, alpha = 0.2) +
  geom_smooth(method = "gam", se = F, color = '#ffa600', linetype = "solid", size = 0.7) +
  labs(x = "Time (child age in years)",
       y = "White matter volume (mm^3)") +
  theme_minimal()

b4 <- ggplot(df_outcome_long, aes(x = time, y = cc)) +
  geom_point(color = '#bc5090', size = 0.9, alpha = 0.2) +
  geom_smooth(method = "gam", se = F, color = '#ffa600', linetype = "solid", size = 0.7) +
  labs(x = "Time (child age in years)",
       y = "Cerebellum volume (mm^3)") +
  theme_minimal()

b5 <- ggplot(df_outcome_long, aes(x = time, y = md)) +
  geom_point(color = '#bc5090', size = 0.9, alpha = 0.2) +
  geom_smooth(method = "gam", se = F, color = '#ffa600', linetype = "solid", size = 0.7) +
  labs(x = "Time (child age in years)",
       y = "Global mean diffusivity") +
  theme_minimal()

b6 <- ggplot(df_outcome_long, aes(x = time, y = fa)) +
  geom_point(color = '#bc5090', size = 0.9, alpha = 0.2) +
  geom_smooth(method = "gam", se = F, color = '#ffa600', linetype = "solid", size = 0.7) +
  labs(x = "Time (child age in years)",
       y = "Global fractional anisotropy") +
  theme_minimal()

b7 <- ggplot(df_outcome_long, aes(x = time, y = ge)) +
  geom_point(color = '#bc5090', size = 0.9, alpha = 0.2) +
  geom_smooth(method = "gam", se = F, color = '#ffa600', linetype = "solid", size = 0.7) +
  labs(x = "Time (child age in years)",
       y = "Global efficiency") +
  theme_minimal()

ggarrange(b1, b2, b3, b4, b5, b6, b7, labels = 'auto', nrow = 2, ncol = 4)

# ============================================================================ #

# Sample size repeated measures per modality per exposure
df_descriptives$mia_index_average <- rowMeans(df_descriptives[, c('mia_index_t1', 'mia_index_t2')], na.rm = TRUE)
df_descriptives$hs_crp_average <- rowMeans(df_descriptives[, c('HsCRPmgL_g1', 'HsCRPmgL_g2')], na.rm = TRUE)

smri_twoscans_crp <- subset(df_descriptives, complete.cases(hs_crp_average) & complete.cases(genr_tbv_f09) & complete.cases(genr_tbv_f13))
smri_onescan_crp <- subset(df_descriptives, complete.cases(hs_crp_average) & (complete.cases(genr_tbv_f09) & is.na(genr_tbv_f13)) | (complete.cases(genr_tbv_f13) & is.na(genr_tbv_f09)))
smri_twoscans_mia <- subset(df_descriptives, complete.cases(mia_index_average) & complete.cases(genr_tbv_f09) & complete.cases(genr_tbv_f13))
smri_onescan_mia <- subset(df_descriptives, complete.cases(mia_index_average) & (complete.cases(genr_tbv_f09) & is.na(genr_tbv_f13)) | (complete.cases(genr_tbv_f13) & is.na(genr_tbv_f09)))

dti_twoscans_crp <- subset(df_descriptives, complete.cases(hs_crp_average) & complete.cases(mean_MD_genr_f09) & complete.cases(mean_MD_genr_f13))
dti_onescan_crp <- subset(df_descriptives, complete.cases(hs_crp_average) & (complete.cases(mean_MD_genr_f09) & is.na(mean_MD_genr_f13)) | (complete.cases(mean_MD_genr_f13) & is.na(mean_MD_genr_f09)))
dti_twoscans_mia <- subset(df_descriptives, complete.cases(mia_index_average) & complete.cases(mean_MD_genr_f09) & complete.cases(mean_MD_genr_f13))
dti_onescan_mia <- subset(df_descriptives, complete.cases(mia_index_average) & (complete.cases(mean_MD_genr_f09) & is.na(mean_MD_genr_f13)) | (complete.cases(mean_MD_genr_f13) & is.na(mean_MD_genr_f09)))

fmri_twoscans_crp <- subset(df_descriptives, complete.cases(hs_crp_average) & complete.cases(wiNetwork_none_f09) & complete.cases(wiNetwork_none_f13))
fmri_onescan_crp <- subset(df_descriptives, complete.cases(hs_crp_average) & (complete.cases(wiNetwork_none_f09) & is.na(wiNetwork_none_f13)) | (complete.cases(wiNetwork_none_f13) & is.na(wiNetwork_none_f09)))
fmri_twoscans_mia <- subset(df_descriptives, complete.cases(mia_index_average) & complete.cases(wiNetwork_none_f09) & complete.cases(wiNetwork_none_f13))
fmri_onescan_mia <- subset(df_descriptives, complete.cases(mia_index_average) & (complete.cases(wiNetwork_none_f09) & is.na(wiNetwork_none_f13)) | (complete.cases(wiNetwork_none_f13) & is.na(wiNetwork_none_f09)))

#\

#------------------------------------------------------------------------#
####---------------Step 4: Wide to long transformation---------------####
#-----------------------------------------------------------------------#

# Drop columns we do not need
full_df_cleaned <- dplyr::select(full_df, -c(18:19, 22:69, 77, 114))

# Transform to long dataset 
df_long <- full_df_cleaned %>% pivot_longer(-c(IDC, IDM, MOTHER, AGE_M_v2, BMI_0, ETHNMv2,mdrink_updated,SMOKE_ALL,GSI, GENDER, GESTBIR, EDUCM, INCOME, PARITY, mat_inflam, im, pm, include, GestagePLg1, GestagePLg2, remark_cat_1, remark_cat_2, IL1ß_ObsConcg1, IL1ß_ObsConcg2, IL6_ObsConcg1, IL6_ObsConcg2, IL17A_ObsConcg1, IL17A_ObsConcg2, IFN.._ObsConcg1, IFN.._ObsConcg2, IL23_ObsConcg1, IL23_ObsConcg2, mia_index_t1, mia_index_t2, gestage_plasma_g1, gestage_plasma_g2, HsCRPmgL_g1, HsCRPmgL_g2, BatchDayg1, BatchDayg2,RackIDg1,RackIDg2, include), names_to = c('variable', 'timepoint'), values_to = c('value'), names_pattern = "(.*)_(.*)")

df_long_ordered <- pivot_wider(df_long, names_from = c("variable"), values_from = c("value"))

#---------------------------------------------------------------#
####---------------Step 5: Multiple imputation---------------####
#---------------------------------------------------------------#

# Check missing values per variable
miss_values_function(df_long_ordered)

## Impute missing data in covariates
imputation_function(data = df_long_ordered,  method = "rf",exclude_imp_vars = c(1:3, 19:25, 27:33, 35, 37:101), exclude_predictors = c(1:3, 19, 27, 38:43, 46:101))

## Load imputed data and check convergence 
imputed_df <- readRDS('imputedData_df_long_ordered.rds')

plot(imputed_df)

#\

#-------------------------------------------------------------------------#
####---------------Step 6: Inverse probability weighting---------------####
#-------------------------------------------------------------------------#

# Specify inclusion variable
incl <- "include" 

# Specify predictors 
pred_baseline <- c("AGE_M_v2", "BMI_0" ,"ETHNMv2","SMOKE_ALL","GSI","GENDER","EDUCM","PARITY", "mia_index_t1", "GESTBIR")  

# Use the reformulate function to create a formula
fit_eq <- reformulate(termlabels = pred_baseline, response = incl)
print(fit_eq)

# Create empty list for storing propensity scores after running model on each imputed set
ps <- list()

# For loop to run CBPS model over each imputed set and get ps's
for (imp in 1:imputed_df$m){
  imp_ipw.imp <- complete(imputed_df, imp)
  #get weights for weighting back to baseline
  fit.imp <-  CBPS((fit_eq), data=imp_ipw.imp) 
  #propensity scores 
  ps[[imp]] <- fit.imp$fitted.values 
}

# Use long fomat to paste ps's in imputed sets
long <- complete(imputed_df, action = "long", include =TRUE)
len <- lengths(ps)

# Make a new df with a vector .imp en .id and each propensity score so we can merge it with our long df 
df <- data.frame(.imp = rep(seq_along(len), len), .id=sequence(len), ps=unlist(ps))

merge <- merge(long, df, by=c(".imp", ".id"), all.x=TRUE) 

# Calculate weights from ps's (individuals with a lower probability of being included would get higher weights; hence we take the inverse)
merge2 <- mutate(merge, ipw_weights = ifelse(get(incl) == 1, (1 / ps), (1 / (1-ps))))

# Convert to mids
ipw_mids <- as.mids(merge2)

# Create final analysis set with weights by selecting only in observations with inclusionvar = 1
ipw_mids_final <- filter(ipw_mids, get(incl) == 1)

# Save imputed dataset with final participants and including the weights
saveRDS(ipw_mids_final, 'ipw_mids_final.RDS')

#\

#----------------------------------------------------------------------#
####--------------Step 7: Post imputation manipulations--------------####
#----------------------------------------------------------------------#

# Load the mids object
ipw_mids_final <- readRDS('ipw_mids_final.rds')

# Transform mids object to long dataframe
imp_long <- complete(ipw_mids_final, action = 'long', include = T)

# Create average MIA value
imp_long$mia_index_average <- rowMeans(imp_long[, c('mia_index_t1', 'mia_index_t2')], na.rm = TRUE)

# Create average CRP value
imp_long$hs_crp_average <- rowMeans(imp_long[, c('HsCRPmgL_g1', 'HsCRPmgL_g2')], na.rm = TRUE)

# Create average gestational age value 
imp_long$gestage_plasma_average <- rowMeans(imp_long[, c('gestage_plasma_g1', 'gestage_plasma_g1')], na.rm = TRUE)
imp_long$GestagePLg_average <- rowMeans(imp_long[, c('GestagePLg1', 'GestagePLg2')], na.rm = TRUE)

# Create average batch effects value 
imp_long$BatchDay_average <- rowMeans(imp_long[, c('BatchDayg1', 'BatchDayg2')], na.rm = TRUE)

# Removing columns we don't need
imp_long <- dplyr::select(imp_long, -c('remark_cat_1', 'remark_cat_2'))

# Save back to mids object 
imp.mids_final2 <- as.mids(imp_long)

saveRDS(imp.mids_final2,'ipw_mids_final_cleaned.rds')

#\ END OF SCRIPT, GO TO PART 2 FOR STATISTICAL ANALYSES. 
