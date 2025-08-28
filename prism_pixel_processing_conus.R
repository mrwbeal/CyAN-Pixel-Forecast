# This script is for processing PRISM data 
# Author: Natalie Von Tress
# Date Created: 05/12/2022

# Reference:
# https://cran.r-project.org/web/packages/prism/vignettes/prism.html#download-30-year-normal-data

# Set up workspace ----

# Clear environment
rm(list = ls(all = T))

# Load packages
library(tidyverse)
library(lubridate)
library(terra)
library(exactextractr)
library(sf)
library(prism)
library(tigris)
library(feather)
library(ggplot2)
library(dplyr)
library(parallel)
library(future)
library(future.apply)

# Set working directory
# 
# wd <- '/work/HAB4CAST/data_processing'
# 
# setwd(wd)


setwd("/work/HAB4CAST/max_beal/CyANPixelForecast/")

temp_dir <- './temp_files'
temp_create <- ifelse(!dir.exists(temp_dir),
                      dir.create(file.path(temp_dir),recursive = T),
                      F)
terraOptions(tempdir = temp_dir)

# Read in data ----

# Original csv to add rows to
# old_csv <- read_csv('mean_prism.csv')

# Week assignments
week_assignments <- read_csv('data/week_assignments.csv')

# Projection - from CyAN
albers_terra <- crs(rast(x = list.files(path = 'data/OLCI_preprocessing/output_masked_tif', full.names = TRUE, recursive = T)[1]))

# Lakes
conus_lakes <- read_sf('./data/OLCI_resolvable_lakes_2022_09_08/OLCI_resolvable_lakes_2022_09_08.shp') %>% st_transform(albers_terra)


# CONUS bound
states <- states(cb = TRUE)
states = states %>% filter(STUSPS!="AK"&STUSPS!="HI"&STUSPS!="PR"&STUSPS!="VI"&STUSPS!="GU"&STUSPS!="MP"&STUSPS!="AS")
conus_bound = states %>% st_transform(albers_terra)

# PRISM air temp 
prism_set_dl_dir('/work/HAB4CAST/max_beal/CyANPixelForecast/data/prism_atemp_data')
atemp_files <- file.path('/work/HAB4CAST/max_beal/CyANPixelForecast/data','prism_atemp_data',prism_archive_ls(),paste0(prism_archive_ls(),'.bil'))

#atemp_spat <- rast(atemp_files) %>% project(albers_terra)


# PRISM precip
prism_set_dl_dir('/work/HAB4CAST/max_beal/CyANPixelForecast/data/prism_precip_data')
precip_files <- file.path('/work/HAB4CAST/max_beal/CyANPixelForecast/data','prism_precip_data',prism_archive_ls(),paste0(prism_archive_ls(),'.bil'))

#precip_spat <- rast(precip_files) %>% project(albers_terra)

completed_prism_files = list.files("/work/HAB4CAST/max_beal/CyANPixelForecast/data/prism_pixel_conus_lakes/")
precip_files

# Pull img dates
precip_dates<-precip_files %>% str_extract('[:digit:]{8}') %>% ymd()
atemp_dates<-atemp_files %>% str_extract('[:digit:]{8}') %>% ymd()

completed_dates = as.numeric(gsub("\\D", "", completed_prism_files)) %>% ymd()

atemp_files = atemp_files[!atemp_dates %in% completed_dates]
precip_files=precip_files[!precip_dates %in% completed_dates]

# Extract values ----

prism_pixel_byrast = function(sequence, conus_lakes, file_data){
  
  atemp_spat <- rast(file_data$temp[sequence]) %>% project(albers_terra)
  precip_spat <- rast(file_data$precip[sequence]) %>% project(albers_terra)
  

  atemp_daily_mean_per_lake <- exact_extract(atemp_spat,conus_lakes,include_xy=TRUE,include_cell=TRUE,include_cols = c("COMID"))%>% bind_rows()
  
  precip_daily_mean_per_lake<- exact_extract(precip_spat,conus_lakes,include_xy=TRUE,include_cell=TRUE,include_cols = c("COMID"))%>% bind_rows()
  
  # Reshape data ----
  # Need to reshape the data so there are four columns: 
  ## 1. COMID
  ## 2. Cell
  ## 2. Image date
  ## 3. Average daily value
  ## 4. Geometry
  
  mean_atemp <- atemp_daily_mean_per_lake %>% mutate(date = names(atemp_spat) %>% str_extract('[:digit:]{8}') %>% ymd(),
                                                     variable="temperature")
  
  mean_precip<-precip_daily_mean_per_lake %>% mutate(date = names(precip_spat) %>% str_extract('[:digit:]{8}') %>% ymd(),
                                                     variable="precipitation")
  
  #mean_prism<-merge(mean_atemp,mean_precip)

  mean_prism<-full_join(mean_atemp,mean_precip)
    
  dt =names(precip_spat) %>% str_extract('[:digit:]{8}') %>% ymd()
  
  filename = paste0("/work/HAB4CAST/max_beal/CyANPixelForecast/data/prism_pixel_conus_lakes/daily_prism_",dt,".csv")  
  # Save the data to a CSV file
  write_csv(mean_prism, filename)
  
  return(mean_prism)
  
  
}

#test = file_data[1:10,]
#sequence = seq_along(test$temp)
#r = prism_pixel_byrast(2,conus_lakes = conus_lakes, file_data = file_data)

file_data = data.frame("temp"=atemp_files,"precip"=precip_files)
sequence = seq_along(file_data$temp)

plan(multisession, workers = (32))
daily_prism_output = future_lapply(sequence, prism_pixel_byrast, conus_lakes = conus_lakes, file_data = file_data) %>% bind_rows() #Save as separate files
plan(sequential)

write_csv(daily_prism_output,"/work/HAB4CAST/max_beal/CyANPixelForecast/data/daily_prism_output.csv")

# Delete terra temp files
unlink(temp_dir, recursive = TRUE, force = TRUE) 
