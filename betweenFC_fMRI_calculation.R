#### Calculate between network connectivity for 12 x 12 networks (fMRI) #### 

## Step 1: Open libraries
library(foreign) 
library(haven)
library(dplyr)

##----------------------------

## Step 2: Load data
matrix_f09 <- readRDS('within_between_fullCon_genr_f9_correct.rds')
df1 <- matrix_f09
matrix_f13 <- readRDS('within_between_fullCon_genr_f13_correct.rds')
df2 <- matrix_f13

df_idc <- readRDS("sub-id_has_conMat.rds")

# Extracting characters from matched string in order to get just numbers
df_idc <- as.integer(gsub('sub-','', df_idc)) 

##----------------------------

## Step 3: Extract between connectivity for F09

# Loop to create df with between connectivity for 12 x 12 networks 

# Define the network labels
networks <- c('Default', 'ParietoOccip', 'FrontoParietal', 'Salience', 'CinguloOperc', 'MedialParietal', 'DorsalAttn', 'VentralAttn', 'Visual', 'SMhand', 'SMmouth', 'Auditory')

# initialize dataframe 
fc_fMRI_df <- data.frame()

# for every subject id 
for (i in 1:length(df_idc)) { 
  # for every network as specified in vector 'networks'
  for (n in networks) { 
    # defining subject ID as everyone in df_idc
    subID = df_idc[i] 
    # extracting between fc
    btNetwork_default<- (df1[[i]]['Default',networks])
    btNetwork_ParietoOccip <- (df1[[i]]['ParietoOccip',networks])
    btNetwork_FrontoParietal<- (df1[[i]]['FrontoParietal',networks])
    btNetwork_Salience <- (df1[[i]]['Salience',networks])
    btNetwork_CinguloOperc<- (df1[[i]]['CinguloOperc',networks])
    btNetwork_MedialParietal<- (df1[[i]]['MedialParietal',networks])
    btNetwork_DorsalAttn<- (df1[[i]]['DorsalAttn',networks])
    btNetwork_VentralAttn<- (df1[[i]]['VentralAttn',networks])
    btNetwork_Visual<- (df1[[i]]['Visual',networks])
    btNetwork_SMhand<- (df1[[i]]['SMhand',networks])
    btNetwork_SMmouth<- (df1[[i]]['SMmouth',networks])
    btNetwork_Auditory<- (df1[[i]]['Auditory',networks])
    # save results 
    fc_fMRI_df[i,c(1:145)] <- c(subID, btNetwork_default, btNetwork_ParietoOccip,
                             btNetwork_FrontoParietal, btNetwork_Salience,btNetwork_CinguloOperc,
                             btNetwork_MedialParietal, btNetwork_DorsalAttn,btNetwork_VentralAttn
                             ,btNetwork_Visual,btNetwork_SMhand,btNetwork_SMmouth,btNetwork_Auditory)
  }
}

colnames(fc_fMRI_df) <- c("IDC",
                          
                          "Default-Default",
                          "Default-ParietoOccip", 
                          "Default-FrontoParietal", 
                          "Default-Salience", 
                          "Default-CinguloOperc", 
                          "Default-MedialParietal", 
                          "Default-DorsalAttn", 
                          "Default-VentralAttn", 
                          "Default-Visual", 
                          "Default-SMhand", 
                          "Default-SMmouth", 
                          'Default-Auditory',
                          
                          "ParietoOccip-Default",
                          "ParietoOccip-ParietoOccip", 
                          "ParietoOccip-FrontoParietal", 
                          "ParietoOccip-Salience", 
                          "ParietoOccip-CinguloOperc", 
                          "ParietoOccip-MedialParietal", 
                          "ParietoOccip-DorsalAttn", 
                          "ParietoOccip-VentralAttn", 
                          "ParietoOccip-Visual", 
                          "ParietoOccip-SMhand", 
                          "ParietoOccip-SMmouth", 
                          'ParietoOccip-Auditory',
                          
                          "FrontoParietal-Default",
                          "FrontoParietal-ParietoOccip", 
                          "FrontoParietal-FrontoParietal", 
                          "FrontoParietal-Salience", 
                          "FrontoParietal-CinguloOperc", 
                          "FrontoParietal-MedialParietal", 
                          "FrontoParietal-DorsalAttn", 
                          "FrontoParietal-VentralAttn", 
                          "FrontoParietal-Visual", 
                          "FrontoParietal-SMhand", 
                          "FrontoParietal-SMmouth", 
                          'FrontoParietal-Auditory',
                          
                          "Salience-Default",
                          "Salience-ParietoOccip", 
                          "Salience-FrontoParietal", 
                          "Salience-Salience", 
                          "Salience-CinguloOperc", 
                          "Salience-MedialParietal", 
                          "Salience-DorsalAttn", 
                          "Salience-VentralAttn", 
                          "Salience-Visual", 
                          "Salience-SMhand", 
                          "Salience-SMmouth", 
                          'Salience-Auditory',
                          
                          "CinguloOperc-Default",
                          "CinguloOperc-ParietoOccip", 
                          "CinguloOperc-FrontoParietal", 
                          "CinguloOperc-Salience", 
                          "CinguloOperc-CinguloOperc", 
                          "CinguloOperc-MedialParietal", 
                          "CinguloOperc-DorsalAttn", 
                          "CinguloOperc-VentralAttn", 
                          "CinguloOperc-Visual", 
                          "CinguloOperc-SMhand", 
                          "CinguloOperc-SMmouth", 
                          'CinguloOperc-Auditory',
                          
                          "MedialParietal-Default",
                          "MedialParietal-ParietoOccip", 
                          "MedialParietal-FrontoParietal", 
                          "MedialParietal-Salience", 
                          "MedialParietal-CinguloOperc", 
                          "MedialParietal-MedialParietal", 
                          "MedialParietal-DorsalAttn", 
                          "MedialParietal-VentralAttn", 
                          "MedialParietal-Visual", 
                          "MedialParietal-SMhand", 
                          "MedialParietal-SMmouth", 
                          'MedialParietal-Auditory',

                          "DorsalAttn-Default",
                          "DorsalAttn-ParietoOccip", 
                          "DorsalAttn-FrontoParietal", 
                          "DorsalAttn-Salience", 
                          "DorsalAttn-CinguloOperc", 
                          "DorsalAttn-MedialParietal", 
                          "DorsalAttn-DorsalAttn", 
                          "DorsalAttn-VentralAttn", 
                          "DorsalAttn-Visual", 
                          "DorsalAttn-SMhand", 
                          "DorsalAttn-SMmouth", 
                          'DorsalAttn-Auditory',
                          
                          "VentralAttn-Default",
                          "VentralAttn-ParietoOccip", 
                          "VentralAttn-FrontoParietal", 
                          "VentralAttn-Salience", 
                          "VentralAttn-CinguloOperc", 
                          "VentralAttn-MedialParietal", 
                          "VentralAttn-DorsalAttn", 
                          "VentralAttn-VentralAttn", 
                          "VentralAttn-Visual", 
                          "VentralAttn-SMhand", 
                          "VentralAttn-SMmouth", 
                          'VentralAttn-Auditory',
                          
                          "Visual-Default",
                          "Visual-ParietoOccip", 
                          "Visual-FrontoParietal", 
                          "Visual-Salience", 
                          "Visual-CinguloOperc", 
                          "Visual-MedialParietal", 
                          "Visual-DorsalAttn", 
                          "Visual-VentralAttn", 
                          "Visual-Visual", 
                          "Visual-SMhand", 
                          "Visual-SMmouth", 
                          'Visual-Auditory',
                          
                          "SMhand-Default",
                          "SMhand-ParietoOccip", 
                          "SMhand-FrontoParietal", 
                          "SMhand-Salience", 
                          "SMhand-CinguloOperc", 
                          "SMhand-MedialParietal", 
                          "SMhand-DorsalAttn", 
                          "SMhand-VentralAttn", 
                          "SMhand-Visual", 
                          "SMhand-SMhand", 
                          "SMhand-SMmouth", 
                          'SMhand-Auditory',
                          
                          "SMmouth-Default",
                          "SMmouth-ParietoOccip", 
                          "SMmouth-FrontoParietal", 
                          "SMmouth-Salience", 
                          "SMmouth-CinguloOperc", 
                          "SMmouth-MedialParietal", 
                          "SMmouth-DorsalAttn", 
                          "SMmouth-VentralAttn", 
                          "SMmouth-Visual", 
                          "SMmouth-SMhand", 
                          "SMmouth-SMmouth", 
                          'SMmouth-Auditory',
                          
                          "Auditory-Default",
                          "Auditory-ParietoOccip", 
                          "Auditory-FrontoParietal", 
                          "Auditory-Salience", 
                          "Auditory-CinguloOperc", 
                          "Auditory-MedialParietal", 
                          "Auditory-DorsalAttn", 
                          "Auditory-VentralAttn", 
                          "Auditory-Visual", 
                          "Auditory-SMhand", 
                          "Auditory-SMmouth", 
                          'Auditory-Auditory'
                          )

saveRDS(fc_fMRI_df, 'fc_fMRI_f09.rds')

#\

##----------------------------

## Step 4: Extract between connectivity for F13

# Loop to create df with between connectivity for 12 x 12 networks 

# Define the network labels
networks <- c('Default', 'ParietoOccip', 'FrontoParietal', 'Salience', 'CinguloOperc', 'MedialParietal', 'DorsalAttn', 'VentralAttn', 'Visual', 'SMhand', 'SMmouth', 'Auditory')

# initialize dataframe 
fc_fMRI_df2 <- data.frame()

# for every subject id 
for (i in 1:length(df_idc)) { 
  # for every network as specified in vector 'networks'
  for (n in networks) { 
    # defining subject ID as everyone in df_idc
    subID = df_idc[i] 
    # extracting between fc
    btNetwork_default<- (df2[[i]]['Default',networks])
    btNetwork_ParietoOccip <- (df2[[i]]['ParietoOccip',networks])
    btNetwork_FrontoParietal<- (df2[[i]]['FrontoParietal',networks])
    btNetwork_Salience <- (df2[[i]]['Salience',networks])
    btNetwork_CinguloOperc<- (df2[[i]]['CinguloOperc',networks])
    btNetwork_MedialParietal<- (df2[[i]]['MedialParietal',networks])
    btNetwork_DorsalAttn<- (df2[[i]]['DorsalAttn',networks])
    btNetwork_VentralAttn<- (df2[[i]]['VentralAttn',networks])
    btNetwork_Visual<- (df2[[i]]['Visual',networks])
    btNetwork_SMhand<- (df2[[i]]['SMhand',networks])
    btNetwork_SMmouth<- (df2[[i]]['SMmouth',networks])
    btNetwork_Auditory<- (df2[[i]]['Auditory',networks])
    # save results 
    fc_fMRI_df2[i,c(1:145)] <- c(subID, btNetwork_default, btNetwork_ParietoOccip,
                                btNetwork_FrontoParietal, btNetwork_Salience,btNetwork_CinguloOperc,
                                btNetwork_MedialParietal, btNetwork_DorsalAttn,btNetwork_VentralAttn
                                ,btNetwork_Visual,btNetwork_SMhand,btNetwork_SMmouth
                                ,btNetwork_Auditory)
  }
}

colnames(fc_fMRI_df2) <- c("IDC",
                          
                          "Default-Default",
                          "Default-ParietoOccip", 
                          "Default-FrontoParietal", 
                          "Default-Salience", 
                          "Default-CinguloOperc", 
                          "Default-MedialParietal", 
                          "Default-DorsalAttn", 
                          "Default-VentralAttn", 
                          "Default-Visual", 
                          "Default-SMhand", 
                          "Default-SMmouth", 
                          'Default-Auditory',
                          
                          "ParietoOccip-Default",
                          "ParietoOccip-ParietoOccip", 
                          "ParietoOccip-FrontoParietal", 
                          "ParietoOccip-Salience", 
                          "ParietoOccip-CinguloOperc", 
                          "ParietoOccip-MedialParietal", 
                          "ParietoOccip-DorsalAttn", 
                          "ParietoOccip-VentralAttn", 
                          "ParietoOccip-Visual", 
                          "ParietoOccip-SMhand", 
                          "ParietoOccip-SMmouth", 
                          'ParietoOccip-Auditory',
                          
                          "FrontoParietal-Default",
                          "FrontoParietal-ParietoOccip", 
                          "FrontoParietal-FrontoParietal", 
                          "FrontoParietal-Salience", 
                          "FrontoParietal-CinguloOperc", 
                          "FrontoParietal-MedialParietal", 
                          "FrontoParietal-DorsalAttn", 
                          "FrontoParietal-VentralAttn", 
                          "FrontoParietal-Visual", 
                          "FrontoParietal-SMhand", 
                          "FrontoParietal-SMmouth", 
                          'FrontoParietal-Auditory',
                          
                          "Salience-Default",
                          "Salience-ParietoOccip", 
                          "Salience-FrontoParietal", 
                          "Salience-Salience", 
                          "Salience-CinguloOperc", 
                          "Salience-MedialParietal", 
                          "Salience-DorsalAttn", 
                          "Salience-VentralAttn", 
                          "Salience-Visual", 
                          "Salience-SMhand", 
                          "Salience-SMmouth", 
                          'Salience-Auditory',
                          
                          "CinguloOperc-Default",
                          "CinguloOperc-ParietoOccip", 
                          "CinguloOperc-FrontoParietal", 
                          "CinguloOperc-Salience", 
                          "CinguloOperc-CinguloOperc", 
                          "CinguloOperc-MedialParietal", 
                          "CinguloOperc-DorsalAttn", 
                          "CinguloOperc-VentralAttn", 
                          "CinguloOperc-Visual", 
                          "CinguloOperc-SMhand", 
                          "CinguloOperc-SMmouth", 
                          'CinguloOperc-Auditory',
                          
                          "MedialParietal-Default",
                          "MedialParietal-ParietoOccip", 
                          "MedialParietal-FrontoParietal", 
                          "MedialParietal-Salience", 
                          "MedialParietal-CinguloOperc", 
                          "MedialParietal-MedialParietal", 
                          "MedialParietal-DorsalAttn", 
                          "MedialParietal-VentralAttn", 
                          "MedialParietal-Visual", 
                          "MedialParietal-SMhand", 
                          "MedialParietal-SMmouth", 
                          'MedialParietal-Auditory',
                          
                          "DorsalAttn-Default",
                          "DorsalAttn-ParietoOccip", 
                          "DorsalAttn-FrontoParietal", 
                          "DorsalAttn-Salience", 
                          "DorsalAttn-CinguloOperc", 
                          "DorsalAttn-MedialParietal", 
                          "DorsalAttn-DorsalAttn", 
                          "DorsalAttn-VentralAttn", 
                          "DorsalAttn-Visual", 
                          "DorsalAttn-SMhand", 
                          "DorsalAttn-SMmouth", 
                          'DorsalAttn-Auditory',
                          
                          "VentralAttn-Default",
                          "VentralAttn-ParietoOccip", 
                          "VentralAttn-FrontoParietal", 
                          "VentralAttn-Salience", 
                          "VentralAttn-CinguloOperc", 
                          "VentralAttn-MedialParietal", 
                          "VentralAttn-DorsalAttn", 
                          "VentralAttn-VentralAttn", 
                          "VentralAttn-Visual", 
                          "VentralAttn-SMhand", 
                          "VentralAttn-SMmouth", 
                          'VentralAttn-Auditory',
                          
                          "Visual-Default",
                          "Visual-ParietoOccip", 
                          "Visual-FrontoParietal", 
                          "Visual-Salience", 
                          "Visual-CinguloOperc", 
                          "Visual-MedialParietal", 
                          "Visual-DorsalAttn", 
                          "Visual-VentralAttn", 
                          "Visual-Visual", 
                          "Visual-SMhand", 
                          "Visual-SMmouth", 
                          'Visual-Auditory',
                          
                          "SMhand-Default",
                          "SMhand-ParietoOccip", 
                          "SMhand-FrontoParietal", 
                          "SMhand-Salience", 
                          "SMhand-CinguloOperc", 
                          "SMhand-MedialParietal", 
                          "SMhand-DorsalAttn", 
                          "SMhand-VentralAttn", 
                          "SMhand-Visual", 
                          "SMhand-SMhand", 
                          "SMhand-SMmouth", 
                          'SMhand-Auditory',
                          
                          "SMmouth-Default",
                          "SMmouth-ParietoOccip", 
                          "SMmouth-FrontoParietal", 
                          "SMmouth-Salience", 
                          "SMmouth-CinguloOperc", 
                          "SMmouth-MedialParietal", 
                          "SMmouth-DorsalAttn", 
                          "SMmouth-VentralAttn", 
                          "SMmouth-Visual", 
                          "SMmouth-SMhand", 
                          "SMmouth-SMmouth", 
                          'SMmouth-Auditory',
                          
                          "Auditory-Default",
                          "Auditory-ParietoOccip", 
                          "Auditory-FrontoParietal", 
                          "Auditory-Salience", 
                          "Auditory-CinguloOperc", 
                          "Auditory-MedialParietal", 
                          "Auditory-DorsalAttn", 
                          "Auditory-VentralAttn", 
                          "Auditory-Visual", 
                          "Auditory-SMhand", 
                          "Auditory-SMmouth", 
                          'Auditory-Auditory'
)

saveRDS(fc_fMRI_df2, 'fc_fMRI_f13.rds')

#\

#\ END OF SCRIPT 
