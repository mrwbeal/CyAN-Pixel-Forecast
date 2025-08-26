# This code is created to look into our ability to predict the continuous log chl a value derived from CyAN
# Created 12/11/2024
# Author Max Beal



library(readr)
library(ranger)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(data.table)
library(caret)
library(keras)
library(tensorflow)
library(reshape2)
library(lubridate)
library(tidymodels)
library(rhdf5)
library(lme4)
library(terra)
#### data ####

CA = F
FL = T
if (CA) {
  df <- read_csv("/work/HAB4CAST/max_beal/CyANPixelForecast/data/california_data.csv")
}

if (FL) {
  df <- fread("/work/HAB4CAST/max_beal/CyANPixelForecast/data/florida_data_spatial.csv")
  
}


#leadtimes = c(2,4,8)
leadtimes = c(4)

#### Run Model ####
#Create the week-year column, make sure weeks have 2 digits

for (lead in leadtimes) {
  
  print(paste0("Building RF model for lead time of: ",lead," weeks"))

  df$Time = paste0(df$Year,str_pad(df$Week,2,pad=0))
  morpho=read.csv( "/work/HAB4CAST/max_beal/CyANPixelForecast/data/morphological_characteristics.csv")
  df=merge(df,morpho, by="COMID")
  
  # Make the example reproducible
  set.seed(1)
  
  input <- df %>%
    arrange(COMID,x,y,Time) %>%
    # Fill Wtemp data
    group_by(COMID,x,y) %>%
    mutate(chl_lead = lead(chl, n = lead), #SET LEAD TIMES HERE
           chl = chl,
           chl_lag2 = lag(chl, n = 1), # Keep lags and other data the same
           chl_lag3 = lag(chl, n = 2),
           chl_lag4 = lag(chl, n = 3)
    ) %>%
    fill(PRISMARDSW,PRISMSW) %>%
    ungroup() %>%   # # Pre-process Bloom data
    filter(!is.na(chl_lead))
  
  
  #Adjust sw temp to use PRISMARDSW when available, PRISMARD as a backup
  #Note PRISMARDSW now represents all surface water predictors
  input = input %>% filter(input$PRISMARDSW %>% is.na()) %>% mutate(PRISMARDSW = PRISMSW) 
  
  padded_numbers <- sprintf("%06d", input$time)
  input$week = substr(padded_numbers,1,2)
  input$year = substr(padded_numbers,3,6)
  
  
  dft=input %>%
    mutate(beginning = ymd(str_c(year, "-01-01")), date = beginning + weeks(week)) %>% 
    dplyr::select(COMID,date,time,x,y,chl_lead,chl,chl_lag2,chl_lag3,chl_lag4,PRISMARDSW,precip,surfaceRadiation,windSpeed,dMean,Week,mean_neighbor_chl,year)
  
  
  
  #### Try subsetting by year ####

  dft$subset = as.numeric(factor(dft$year))
  
  #Verify time blocking
  verify = dft %>% group_by(subset) %>% summarise(start_time=min(date),
                                          end_time=max(date),
                                          n_obs = n())
  
  
  #Number of observations per pixel
  numobs = dft %>% group_by(x,y) %>% summarise(n=length(COMID))
  numobs %>% arrange(desc(n))
  hist(numobs$n)
  
  
  subset_index = dft %>% group_by(date) %>% summarise(subset = mean(subset)) %>% pull(subset)
  
  
  coef_save = data.frame()
  comids = dft$COMID %>% unique()
  
  #### Choose model type ####
  RF = TRUE
  
  #Start at 104 to train on two years to start
  for (i in c(8:max(verify$subset))) {
    
    
    #Subset and scale/center training set
    train = dft %>% filter(subset <= i-1) %>% na.omit()
    
    #Try sequential training
    # train = dft %>% filter(subset==(i-1)) %>% na.omit() 
    train_predictor_scaled = train %>% dplyr::select(chl,chl_lag2,chl_lag3,chl_lag4, PRISMARDSW, precip, surfaceRadiation, windSpeed,mean_neighbor_chl) %>% scale(center=TRUE,scale=TRUE)
    
    
    #Rejoin scaled predictors and predictand in training set
    train = data.frame("chl_lead"=train$chl_lead,"Time"=train$time,"x"=train$x,"y"=train$y,"COMID"=train$COMID,"dMean"=train$dMean,"Week"=train$Week,train_predictor_scaled)
    
    #Subset and scale/center test set
    test = dft %>% filter(subset==i) %>% na.omit()
    test_predictor_scaled = test %>% dplyr::select(chl,chl_lag2,chl_lag3,chl_lag4, PRISMARDSW, precip, surfaceRadiation, windSpeed,mean_neighbor_chl) %>% scale(center=TRUE,scale=TRUE)
    
    #Rejoin scaled predictors and predictand in testing set
    test = data.frame("chl_lead"=test$chl_lead,"Time"=test$time,"x"=test$x,"y"=test$y,"COMID"=test$COMID,"dMean"=test$dMean,"Week"=test$Week,test_predictor_scaled)
    
    
    #Calculate the null model now based on the train informaton we've "seen" before the prediction takes place
    lta = dft %>% filter(subset<=i-1) %>% group_by(Week,x,y) %>% summarise(chl_LTA = mean(chl_lead,na.rm=T))
    lta_xy = dft %>% filter(subset<=i-1) %>% group_by(x,y) %>% summarise(chl_LTA_xy = mean(chl_lead,na.rm=T))
    
  
  
    #Join with the test set
    test = left_join(test,lta,by=c("Week","x","y"))
    test = left_join(test,lta_xy,by=c("x","y"))
    test=test %>% mutate(chl_LTA = ifelse(is.na(chl_LTA),chl_LTA_xy,chl_LTA)) %>% dplyr::select(-chl_LTA_xy)
    
  
    
    
    if (RF) {
      formula = chl_lead ~ chl + chl_lag2 + chl_lag3 + chl_lag4 + PRISMARDSW + precip + surfaceRadiation + windSpeed + dMean + mean_neighbor_chl + Week + x + y  #removed area, need to put back in
      model <- ranger(formula, data = train,num.trees=500,quantreg = TRUE,importance = "impurity",num.threads = 16)
      coefs=data.frame("kfold"=i,t(model$variable.importance))
      coef_save = rbind(coef_save,coefs)
      
      # get Predictions
      ensemble <- predict(model, data=test, type="quantiles", what = function(x) sample(x, 100, replace = TRUE))$predictions #get ensemble prediction
      
      numbers = c(1:100)
      strings <- str_c("enschl_member_", numbers)
      colnames(ensemble) = strings
    
      pred = rowMeans(ensemble)
      
    }
    
    test$chl_predicted = pred
    
    test$error = test$chl_lead - test$chl_predicted
    
    test = test %>% mutate(chl_observed = chl_lead,
                           Bloom_predicted = chl_predicted>=12,
                            Bloom_observed = chl_observed>=12)
    
    test = cbind(test,ensemble) #Return to chlorophyll from log chl
    
  
    res <- caret::postResample(test$chl_lead, test$chl_predicted)
    print(res[2])
    
    
    tbl = table(test$Bloom_predicted,test$Bloom_observed)
    print(tbl)
    accuracy = (tbl[1,1]+tbl[2,2])/(tbl[1,1]+tbl[1,2]+tbl[2,1]+tbl[2,2])
    print(accuracy)
    
    test$lead = lead
    
  
    stop()
    if (FL) {
      #Save kfold results
      write_csv(test,paste0("/work/HAB4CAST/max_beal/CyANPixelForecast/output/FL_RF_fold_",i,"_lag",lead,".csv"))
    }
  
      # #Save coefficients
      write_csv(coef_save,"/work/HAB4CAST/max_beal/CyANPixelForecast/output/FL_RF_stdcoefs.csv")
      # 
      coef_save=read_csv("/work/HAB4CAST/max_beal/CyANPixelForecast/output/FL_RF_stdcoefs.csv")
      flcoef = melt(coef_save,id.vars="kfold") %>% filter(variable!="X.Intercept." & variable!="bootstrap")
      ggplot(flcoef) + geom_boxplot(aes(x=variable,y=value,fill=kfold))
  
  }
  gc()
}

