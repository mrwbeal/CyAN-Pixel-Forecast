# This script is for downloading PRISM data (hopefully :) )
# Author: Natalie Von Tress
# Date Created: 05/10/2022

# Reference:
# https://cran.r-project.org/web/packages/prism/vignettes/prism.html#download-30-year-normal-data

# We want:
## Daily PRISM data 
## from May 2016-May 2019
## Variables: Air temperature (tmean; ATEMP) and preciptiation (ppt; PRECIP)

# Note: Files cannot be directly uploaded to the network drive; you must download locally, then copy the files over to the network

# Set up workspace ----

# Clear environment
rm(list = ls(all = T))

# Load packages
library(tidyverse)
library(raster)
library(prism)

# Download air temperature data ----
# 2017-2024 data

setwd("/work/HAB4CAST/max_beal/CyANPixelForecast")
## Create and set directory
atemp_dir <- './data/prism_atemp_data'

if(!file.exists(atemp_dir)){
  warning("Creating directory for air temp prism data download in ",atemp_dir,"\n")
  dir.create(atemp_dir)
}

prism_set_dl_dir(atemp_dir)

## Download data
get_prism_dailys(
  type = 'tmean',
  minDate = '2016-01-01',
  maxDate = '2024-08-01',
  keepZip = F
)

# Convert to raster files
# atemp_stack <- pd_stack(prism_archive_ls())
#writeRaster(atemp_stack,filename = list.files('/work/HAB4CAST/data_processing/data/prism_atemp_data',full.names = T),format = 'raster',bylayer = T)

# Download precipitation data ----
# 2017-2021 data

## Create and set directory
precip_dir <- './data/prism_precip_data'

if(!file.exists(precip_dir)){
  warning("Creating directory for air temp prism data download in ",precip_dir,"\n")
  dir.create(precip_dir)
}

prism_set_dl_dir(precip_dir)

## Download data
get_prism_dailys(
  type = 'ppt',
  minDate = '2016-01-01',
  maxDate = '2024-08-01',
  keepZip = F
)

# precip_stack <- pd_stack(prism_archive_ls())
#writeRaster(precip_stack,filename = list.files('/work/HAB4CAST/data_processing/data/prism_precip_data',full.names = T),format = 'raster',bylayer = T)

