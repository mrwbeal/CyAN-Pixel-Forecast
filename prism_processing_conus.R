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

# Set working directory
# 
# wd <- '/work/HAB4CAST/data_processing'
# 
# setwd(wd)

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
albers_terra <- crs(rast(x = list.files(path = './OLCI_preprocessing/output_masked_tif', full.names = TRUE, recursive = T)[1]))

# Lakes
conus_lakes <- read_sf('./data/OLCI_resolvable_lakes_2022_09_08/OLCI_resolvable_lakes_2022_09_08.shp') %>% st_transform(albers_terra)

# CONUS bound
conus_bound <- read_sf('./data/conus_boundary/conus_boundary.shp')  %>% st_transform(albers_terra)

# PRISM air temp 
prism_set_dl_dir('./data/prism_atemp_data')
atemp_files <- file.path('data','prism_atemp_data',prism_archive_ls(),paste0(prism_archive_ls(),'.bil'))
length(atemp_files)
atemp_spat <- rast(atemp_files) %>% project(albers_terra)

# PRISM precip
prism_set_dl_dir('./data/prism_precip_data')
precip_files <- file.path('data','prism_precip_data',prism_archive_ls(),paste0(prism_archive_ls(),'.bil'))
head(precip_files)
precip_spat <- rast(precip_files) %>% project(albers_terra)

# Pull img dates
names(atemp_spat) %>% str_extract('[:digit:]{8}') %>% ymd() -> prism_dates

# Extract values ----

exact_extract(atemp_spat, # This should be your SpatRaster
              conus_lakes, # This should be your vector
              fun = 'mean') %>% # This tells the function you want to calculate a weighted mean for each image and each lake
  as_tibble() %>%
  bind_cols(dplyr::select(conus_lakes,COMID)) -> atemp_daily_mean_per_lake

exact_extract(precip_spat, 
              conus_lakes,
              fun = 'mean') %>%
  as_tibble() %>%
  bind_cols(dplyr::select(conus_lakes,COMID)) -> precip_daily_mean_per_lake

# Reshape data ----

# Need to reshape the data so there are four columns: 
## 1. COMID
## 2. Image date
## 3. Average daily value
## 4. Geometry

atemp_daily_mean_per_lake %>%
  # Reshape data
  pivot_longer(cols = colnames(atemp_daily_mean_per_lake)[1:(ncol(atemp_daily_mean_per_lake)-2)], 
               values_to = 'daily_atemp') %>%
  # Pull image dates from layer names and create a year column
  mutate(date = str_extract(name,'[:digit:]{8}') %>% ymd()) %>%
  # Assign weeks
  inner_join(week_assignments) %>%
  # Group data by COMID, week, and year
  group_by(COMID,year,week) %>%
  # Compute average value
  summarise(COMID = COMID,
            year = year,
            week = week,
            mean_atemp = mean(daily_atemp)) %>%
  distinct() -> mean_atemp

precip_daily_mean_per_lake %>%
  # Reshape data
  pivot_longer(cols = colnames(precip_daily_mean_per_lake)[1:(ncol(precip_daily_mean_per_lake)-2)], 
               values_to = 'daily_precip') %>%
  # Pull image dates from layer names and create a year column
  mutate(date = str_extract(name,'[:digit:]{8}') %>% ymd()) %>%
  # Assign weeks
  inner_join(week_assignments) %>%
  # Group data by COMID, week, and year
  group_by(COMID,year,week) %>%
  # Compute average value
  summarise(COMID = COMID,
            year = year,
            week = week,
            mean_precip = mean(daily_precip)) %>%
  distinct() -> mean_precip

full_join(mean_atemp,mean_precip) -> mean_prism

# unique(full_join(mean_prism, old_csv)) -> new_csv

# Save workspace - last saved 5/13/22
# save.image("~/cyan_forecasting/prism_processing_fl_workspace.RData")

write_csv(mean_prism,'./data/mean_prism_conus.csv')

# Delete terra temp files
unlink(temp_dir, recursive = TRUE, force = TRUE) 
