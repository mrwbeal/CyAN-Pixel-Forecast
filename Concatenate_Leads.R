#This notebook supports the pixel level cyanobacteria forecast development test case for Florida. This script concatenates the train/test folds from 1, 2, and 4
# week lead times.
# Author: Max Beal
# Date Created: 2/5/2025

library(readr)
library(ranger)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(data.table)
library(caret)
library(reshape2)
library(lubridate)
library(sf)
library(tigris)


#### Random Forest ####
for (j in c(2,4)) {
  
    fileresults=list()
    for (i in c(3:8)) {
      fileresults[[i-2]]=paste0("/work/HAB4CAST/max_beal/CyANPixelForecast/output/FL_RF_fold_",i,"_lag",j,".csv")
    }
    
    
    
    rf_folds = lapply(fileresults,fread) %>% bind_rows()
    #rf_folds$chl_LTA = exp(rf_folds$logchl_LTA)
    
    
    #Get years
    rf_folds$year <- as.numeric(substr(rf_folds$Time, nchar(rf_folds$Time) - 3, nchar(rf_folds$Time)))
    #Separate LTA to merge 
    #lta = rf_folds %>% dplyr::select(Time,x,y,chl_LTA)
    rf_folds=rf_folds %>% mutate(Bloom_observed = as.numeric(chl_observed>=12),
                                 Bloom_predicted = as.numeric(chl_predicted>=12),
                                 Bloom_observed_al2 = as.numeric(chl_observed>=24),
                                 Bloom_predicted_al2 = as.numeric(chl_predicted>=24),
                                 Bloom_LTA_al1 = as.numeric(chl_LTA>=12),
                                 Bloom_LTA_al2 = as.numeric(chl_LTA>=24)) 
    
    #Get the observed and predicted category by ensemble vote
    ensemble=rf_folds %>% dplyr::select(x,y,Time,contains("enschl")) %>%
      melt(id.vars=c("x","y","Time")) %>% mutate(Bloom_predicted_al1_ens = as.numeric(value>=12 & value<24),Bloom_predicted_al2_ens = as.numeric(value>=24)) %>%
      group_by(Time,x,y) %>% summarise(ensemble_AL1_vote = sum(Bloom_predicted_al1_ens),
                                       ensemble_AL2_vote = sum(Bloom_predicted_al2_ens),
                                       ensemble_nc_vote = 100-(ensemble_AL1_vote+ensemble_AL2_vote)) %>% mutate(Bloom_predicted_al1_ens = as.numeric(ensemble_AL1_vote>ensemble_AL2_vote & ensemble_AL1_vote>ensemble_nc_vote),
                                                                                                                Bloom_predicted_al2_ens = as.numeric(ensemble_AL2_vote>ensemble_AL1_vote & ensemble_AL2_vote>ensemble_nc_vote))
    
    rf_folds = merge(rf_folds,ensemble,by=c("x","y","Time"))
    
    savename = paste0("/work/HAB4CAST/max_beal/CyANPixelForecast/output/RF_folds",j,"lead.csv")
    
    write_csv(rf_folds,savename)
    
    
}

