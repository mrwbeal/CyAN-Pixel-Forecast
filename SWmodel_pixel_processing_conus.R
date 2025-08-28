#This code is adapted from the surface water model by Hannah Ferriby
#The purpose of this code is to automate the surface water predictions and allow for pixel-level extraction
#Authors: Max Beal
#Date: 10/4/24

#Load libraries and set seed
library(tidyverse)
library(dplyr)
library(sf)
library(lubridate)
library(randomForest)
library(arrow)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(readr)
library(feather)
library(arrow)
library(stringr)
library(tidyr)
library(sf)
library(zoo)
library(lubridate)
library(future)
library(future.apply)
library(parallel)
library(rlang)
library(stats)
set.seed(42)


setwd('/./work/HAB4CAST/max_beal/SW_model')

#Settings
runKfold = TRUE #Set to true to run 5-fold cross validation on training set


### Build and Validate Surface Water Model ###
#Data Input ----
#Load in lake shapefiles, lake morpho data, PRISM, ice presence, and ARD water temp data
lakes <- st_read("data/OLCI_resolvable_lakes_2022_09_08/OLCI_resolvable_lakes_2022_09_08.shp")
lake_crs = st_crs(lakes)

#Read in the insitu data
in_situ_all <- read_feather("/work/HAB4CAST/max_beal/SW_model/data/all_insitu_2007_2022.feather")

training <- in_situ_all %>% filter(subset == 'Training')
validation <- in_situ_all %>% filter(subset == 'Validation')

predict_2022 <- read_feather("data/prediction_2022.feather")

#Random Forest Model ----
#Train random forest model with date (day of year), elevation (m), latitude, longitude, lake shoreline length (km),
#lake area(sq km), day of air temp (C), and average air temp 30 day prior (C)

formula <- TEMPERATURE ~ LAT + LONG + day_of_year + ElevWs + daily_atemp + mean_30day + lake_sa + lake_shoreline 


k_fold_cv_rf <- function(data,formula, k = 5, ntree = 100, mtry = NULL) {
  # Shuffle the data
  set.seed(123)
  data <- data[sample(nrow(data)), ]
  
  # Create k folds, equal size
  folds <- cut(seq(1, nrow(data)), breaks = k, labels = FALSE)
  
  # Observations and Predictions
  yhat <- c()
  y <- c()
  
  # Perform k-fold cross-validation
  for (i in 1:k) {
    # Split the data into training and testing sets
    test_indices <- which(folds == i, arr.ind = TRUE)
    test_data <- data[test_indices, ]
    train_data <- data[-test_indices, ]
    
    # Train the randomForest model
    #rf_model <- randomForest(formula,data = train_data,ntree = ntree,importance = T,keep.inbag = T,keep.forest = T,
    #                                     na.action=na.exclude)
    
    rf_model <- randomForest(formula,data = in_situ_all,ntree = 100,na.action=na.exclude)
    
    # Make predictions on the test set
    predictions <- predict(rf_model, test_data)
    
    yhat <- c(yhat, predictions)
    y <- c(y, test_data$TEMPERATURE)
  }
  
  outputs = list()
  df <-data.frame("Observations"=y,"Predictions"=yhat)
  
  outputs[[1]] = df
  outputs[[2]] = rf_model

  return(outputs)
}

if (runKfold) {
  #Kfold cross validation
  outputs = k_fold_cv_rf(training,formula)
  df = outputs[[1]]
  rf_model=outputs[[2]]
  plot(df$Observations,df$Predictions,main="K-fold CV output")
  
  #Use model on validation data
  validation$apply_rf <- predict(rf_model,
                                 newdata = validation,
                                 na.rm = T)
  
  plot(validation$TEMPERATURE,validation$apply_rf,main="Validation set")
  
}




### Apply Model to PRISM Data ###

setwd("/work/HAB4CAST/max_beal/CyANPixelForecast")
week_assignments <- read_csv('data/week_assignments.csv')

prism_files = list.files("/work/HAB4CAST/max_beal/CyANPixelForecast/data/prism_pixel_conus_lakes",full.names = TRUE)
cyan_files = list.files("/work/HAB4CAST/max_beal/CyANPixelForecast/data/cyan_pixel_conus_lakes_current",full.names = TRUE)

read_feather(cyan_files[2])$date %>% tail()

#Need lakes for morphological characteristics
lakes <- st_read("data/OLCI_resolvable_lakes_2022_09_08/OLCI_resolvable_lakes_2022_09_08.shp")
morpho = read_feather("/work/HAB4CAST/max_beal/SW_model/data/conus_lake_morpho.feather") 
elevation = read_feather("/work/HAB4CAST/max_beal/SW_model/data/Elevation.feather")

morpho = merge(morpho,elevation,by="COMID")

lakes = merge(lakes,morpho,by="COMID")

COMIDs = as.numeric(str_extract(substr(cyan_files,15,150),"\\d+"))

#check already processed files
processed = list.files("/work/HAB4CAST/max_beal/CyANPixelForecast/data/swtemp_pixel_conus_lakes_current")
COMID_processed = as.numeric(str_extract(processed,"\\d+"))

#Save COMIDs that do not yet have a file
COMIDs = COMIDs[!COMIDs %in% COMID_processed]

#Function to get all prism temperature data, need daily temp, doy, and 30_day rolling mean
#Will also need the lat/lon from the cells
subset_prism_files_sw_model = function(prism_file, COMID){
  
  prism_data = read_csv(prism_file,show_col_types = FALSE)
  
  if (COMID %in% prism_data$COMID) {
    prism_subset=prism_data[prism_data$COMID==COMID,] #Subset COMID
    
    prism_subset$x = round(prism_subset$x) #tiny differences if not rounded dont allow pivot
    prism_subset$y = round(prism_subset$y)
    
    prism_subset = prism_subset %>% pivot_wider(names_from = variable,values_from = value,
                                                id_cols = c("COMID","date","x","y","cell"))
    
    prism_subset = merge(prism_subset,week_assignments)
    

    
    
    return(prism_subset)
  }else{
    
  }
}

#Function to help join prism and cyan by week and year
filter_and_join_date <- function(df1, df2, date) {
  # Filter both dataframes on year and week
  df1_filtered <- df1 %>% filter(date == !!date)
  df2_filtered <- df2 %>% filter(date == !!date)
  
  # Perform nearest feature join
  joined <- st_join(df1_filtered, df2_filtered, join = st_nearest_feature)
  
  return(joined)
}


join_cyan_sw_model = function(COMID, cyan_files, prism_files){
  
  #Load in CyAN File 
  cfile = cyan_files[as.numeric(str_extract(substr(cyan_files,15,150),"\\d+")) %in% COMID]  #Get CyAN file of interest

  cdf = read_feather(cfile)%>%
    st_as_sf(coords=c("x","y"),remove=FALSE)
  
  prism_df = future_lapply(prism_files,subset_prism_files_sw_model,COMID=COMID) %>% bind_rows() %>%
    st_as_sf(coords=c("x","y"))
  
  
  prism_df<-prism_df %>% mutate(atemp_lag = lag(temperature,1),
                                        mean_30day = rollmean(atemp_lag,30, fill = NA, align = 'right'),
                                        day_of_year = yday(date))
  
  
  cdf = cdf %>% select(geometry,x,y) %>% distinct()
  cdf = cdf %>% mutate(id = 1:nrow(cdf))
  
  date_sequence <- unique(prism_df$date)
  
  expanded_data <- expand.grid(
    id = cdf$id,
    date = date_sequence
  )
  
  cdf = merge(expanded_data,cdf,by="id")
  cdf <- st_as_sf(cdf)
  
  cyan_prism_df = st_join(cdf,prism_df,join=st_nearest_feature) %>% filter(date.y==date.x)
  
  cyan_prism_df = lapply(1:length(date_sequence), function(i) {
    filter_and_join_date(cdf, prism_df, date_sequence[i])
  }) %>% bind_rows()
  
  return(cyan_prism_df)
}


cdf = read_feather(cyan_files[2])%>%
  st_as_sf(coords=c("x","y"),remove=FALSE)




apply_sw_model = function(COMID, cyan_files, prism_files, morpho, formula,in_situ_all, crs){
  
  #df = join_cyan_sw_model(COMID,cyan_files = cyan_files, prism_files = prism_files)
  
  df = future_lapply(COMID, join_cyan_sw_model,cyan_files = cyan_files, prism_files = prism_files) %>% bind_rows()
  
  #Join morpho characteristics
  df = df %>% rename(daily_atemp = temperature)
  df = merge(df,morpho, by="COMID")
  
  #Get LAT/LONG predictors
  df <- df %>% dplyr::mutate(LONG = sf::st_coordinates(geometry)[,1],
                             LAT = sf::st_coordinates(geometry)[,2])
  

  #Have to define the RF model within the future lapply, randomly sample to speed processing
  set.seed(123)  # Setting seed for reproducibility
  sample_size <- floor(0.4 * nrow(in_situ_all))  # Calculate 40% of the dataset size
  sampled_data <- in_situ_all[sample(seq_len(nrow(in_situ_all)), size = sample_size), ]
  rf_model <- randomForest(formula,data = sampled_data,ntree = 100,na.action=na.exclude)
  
  #Apply SW model
  df$sw <- predict(rf_model,newdata = df,na.rm = T)


  #Summarise to week
  weekly_output<-df %>% group_by(COMID, geometry, year, week) %>% summarise(mean_sw_temp = mean(sw,na.rm=T),
                                                       ElevWs = mean(ElevWs),
                                                       lake_shoreline = mean(lake_shoreline),
                                                       lake_sa = mean(lake_sa))
  #st_crs(weekly_output) = lake_crs


  weekly_output= weekly_output %>% dplyr::mutate(x = sf::st_coordinates(geometry)[,1],
                                                 y = sf::st_coordinates(geometry)[,2]) %>% st_drop_geometry()
  
  filename = paste0("/work/HAB4CAST/max_beal/CyANPixelForecast/data/swtemp_pixel_conus_lakes_current/COMID_swtemp_",COMID,".csv")
  #Save the data to a CSV file
  write_csv(weekly_output, filename)

  rm(list = setdiff(ls(), c("COMIDs", "cyan_files","prism_files","morpho","formula","in_situ_all","lake_crs","week_assignments",
                            "apply_sw_model","filter_and_join_date","join_cyan_sw_model","subset_prism_files_sw_model")))
  
  gc()
  return()
}


#Run the script in parallel and save as separate COMIDs
plan(multisession, workers = (16))
future_lapply(COMIDs, apply_sw_model, cyan_files = cyan_files, prism_files = prism_files, morpho=morpho,
                    formula=formula,in_situ_all=in_situ_all, crs=lake_crs)
plan(sequential)



prism_files %>% tail()

