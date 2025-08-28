# This script is for scaling and extracting CI data
# Authors: Natalie Von Tress and Hannah Ferriby
# Date Created: 05/16/2022

# Set up workspace ----

# Clear environment
rm(list = ls(all = T))

library(tidyverse)
library(lubridate)
library(exactextractr)
library(sf)
library(terra)
library(future)
library(future.apply)
library(parallel)
library(feather)
library(ggplot2)

t_init = Sys.time()


temp_dir <- './temp_files'
temp_create <- ifelse(!dir.exists(temp_dir),
                      dir.create(file.path(temp_dir),recursive = T),
                      F)
terraOptions(tempdir = temp_dir)

# Read in data ----

setwd("/work/HAB4CAST/max_beal/CyANPixelForecast/")
#Week Assignments
week_assignments <- read_csv('data/week_assignments.csv')
#OLCI Resolvable LA
conus_lakes <- read_sf('data/OLCI_resolvable_lakes_2022_09_08/OLCI_resolvable_lakes_2022_09_08.shp')
#CyAN Files
files <- list.files(path = 'data/OLCI_preprocessing/output_masked_tif', full.names = TRUE, recursive = T)
ci_crs = crs(rast(x = files[1]))

tifNames <- files %>% str_extract('CI_[:digit:]{14}_out')


# # Check for pre-existing mean_cyan_conus and subset files and tifnames accordingly
# if(file.exists('data/mean_cyan_conus.csv')){
#   # Read in existing data
#   mean_cyan_preexisting <- read_csv('data/mean_cyan_conus.csv')
# 
#   # Subset for images that have yet to be processed
#   processed_scene_names <- mean_cyan_preexisting$scene_name
#   unprocessed_sub <- (tifNames %in% processed_scene_names) == F
# 
#   files <- files[unprocessed_sub]
#   tifNames <- tifNames[unprocessed_sub]
# 
# 
#   if(length(files) == 0) {
#     stopifnot(length(files) == 0)
#   } #else {
#   #   ci_brick <- ci_brick[[unprocessed_sub]]
#   # }
# }


if(length(files) == 0) {

  # Stop the code if there's nothing new to process ----

  ot <- Sys.time() - t_init
  print(ot)
  cat("No new images to process -- exiting the script")
  stopifnot(length(files) == 0)

} else {
  
  # Run the code if there are unprocessed images ----

  #Function to get week of year from the date instead of using "week assignments"
  get_week_of_year <- function(date_str) {
    # Convert the string to a Date object
    date_obj <- as.Date(date_str, format = "%Y-%m-%d")
    # Extract the week number
    week_number <- format(date_obj, "%U")
    return(as.numeric(week_number))
  }
  
  #Transform CRS of the shapefile
  ci_crs = crs(rast(x = files[1]))
  conus_lakes <- st_transform(conus_lakes,ci_crs)
  
  #Create Raster Brick
  ci_brick <- rast(x = files)
  
  #Remove duplicates
  ci_brick = ci_brick[[!duplicated(names(ci_brick))]]
  
  
make_extracted_tibble <- function(conus_lake, ci_brick) {
  
  #Create Raster Brick
  ci_brick <- rast(x = files) #Function only works when files are within function for some reason. Look into.

  #Remove duplicates
  ci_brick = ci_brick[[!duplicated(names(ci_brick))]]
  
  lake <- conus_lakes %>% filter(COMID == conus_lake)
  
  
  data<-exact_extract(ci_brick,lake,include_xy=TRUE,include_cell=TRUE)[[1]] %>%
    filter(coverage_fraction == 1) %>%
    dplyr::select(-coverage_fraction) %>%
    mutate(COMID = conus_lake) %>%
    pivot_longer(cols = starts_with("CI")) %>% #pivot dataset so it's indexed by time
    mutate(date = as.Date(str_extract(name,'[:digit:]{7}') %>% as.Date(format = '%Y%j'))) %>% #Get date from image name
    mutate(year = format(date,"%Y")) %>%
    mutate(Week = get_week_of_year(date)) %>%
    mutate(bloom = ifelse(value>=130,T,F)) # 97 for 3 ug/L, 130 for 12 ug/L, 151 for 24 ug/L - subbing u for mu, so ug = micrograms
    
  
    filename = paste0("/work/HAB4CAST/max_beal/CyANPixelForecast/data/cyan_pixel_conus_lakes_current/COMID_",conus_lake,".feather")  
    # Save the data to a CSV file
    write_feather(data, filename)
    
    return("Data saved successfully")
    
}
  

  
  
  #Parallell processing
  plan(multisession, workers = (16))

  #Save as separate files to avoid huge data
  future_lapply(conus_lakes$COMID ,make_extracted_tibble, ci_brick=ci_brick) #Try without reading raster in function
  
  #extracted_tibble <- future_lapply(conus_lakes$COMID,make_extracted_tibble,future.seed = TRUE) %>% bind_rows()
  
  plan(sequential)
  
  cat(str_c('\nFinished extracting values!'))
  
  
  cell_cyan_extract = extracted_tibble
  
  #write_csv(cell_cyan_extract,file = 'data/cell_cyan_conus.csv')
  #write_feather(cell_cyan,"data/cell_cyan_conus.feather")
  
  cat('Wrote mean_cyan_conus.csv\n')
  ot <- Sys.time() - t_init
  print(ot)
}

# Delete terra temp files
unlink(temp_dir, recursive = TRUE, force = TRUE)


