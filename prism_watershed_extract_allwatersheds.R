# This code is to investigate extraction of prism precipitation data from lake watersheds
# Created 11/22/2024
# Author: Max Beal

setwd("/work/HAB4CAST/max_beal/CyANPixelForecast/")

library(tidyverse)
library(lubridate)
library(terra)
library(exactextractr)
library(sf)
library(prism)
library(parallel)
library(future)
library(future.apply)


#Watershed File
file = "/work/HAB4CAST/max_beal/CyANPixelForecast/data/CyAN_watersheds.gpkg"

#st_layers(file)

ws = st_read(file)


# PRISM precip
prism_set_dl_dir('./data/prism_precip_data')
precip_files <- file.path('data','prism_precip_data',prism_archive_ls(),paste0(prism_archive_ls(),'.bil'))
head(precip_files)

#Cyan projection
albers = "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"


# Pull img dates
precip_files %>% str_extract('[:digit:]{8}') %>% ymd() -> prism_dates
time = paste0(isoweek(prism_dates),year(prism_dates))
unique_times = unique(time)


#Processed files
processed = list.files("/work/HAB4CAST/max_beal/CyANPixelForecast/data/prism_watershed_weeksum_csv")
processed = gsub("\\D", "", processed)

#subset those weeks not yet processed
unique_times = unique_times[!unique_times %in% processed]

wshed_extract = function(t) {
  library(exactextractr)
  #Read watersheds
  ws = st_read("/work/HAB4CAST/max_beal/CyANPixelForecast/data/CyAN_watersheds.gpkg")
  albers = "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  ws = st_transform(ws,albers)
  
  #Get dates and week-years
  prism_set_dl_dir('./data/prism_precip_data')
  precip_files <- file.path('data','prism_precip_data',prism_archive_ls(),paste0(prism_archive_ls(),'.bil'))
  prism_dates = precip_files %>% str_extract('[:digit:]{8}') %>% ymd()
  time = paste0(isoweek(prism_dates),year(prism_dates))
  
  print("getting files")
  
  precip_spat <- rast(precip_files[time==t]) %>% project(albers)
  precip_spat=sum(precip_spat,na.rm=T)
  
  print(precip_spat)
  
  week = exact_extract(precip_spat,ws,fun="sum",append_cols="COMID")
  week$Time = t
  
  
  filename=paste0("/work/HAB4CAST/max_beal/CyANPixelForecast/data/prism_watershed_weeksum_csv/prism_watershedSum_",t,".csv")
  write_csv(week,filename)
  
  gc()
  
  
}



plan(multisession, workers = (availableCores()-1))
print(availableCores())
future_lapply(unique_times,wshed_extract) #Save apparently needs to be within the function
plan(sequential)


#Save all files and make plots
ws = list.files("/work/HAB4CAST/max_beal/CyANPixelForecast/data/prism_watershed_weeksum_csv",full.names=TRUE)

wshed_df = read_csv(ws)

#write_csv(wshed_df,"/work/HAB4CAST/max_beal/CyANPixelForecast/data/precip_watershed_data.csv")

wshed_df$timepad = sprintf("%06d",wshed_df$Time) # fix to 2 characters 
week = substr(wshed_df$timepad,1,2)
year = substr(wshed_df$timepad,3,6)

wshed_df$date <- as.Date(paste(year, week, 1, sep = "-"), format = "%Y-%U-%u")


test = wshed_df %>% filter(COMID==10380929)

ggplot(data=test) + geom_line(aes(x=date,y=meantotalp))

barplot(test$meantotalp,names.arg=test$date)


