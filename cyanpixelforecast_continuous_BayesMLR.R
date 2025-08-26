#This script tests the use of brms MCMC Bayesian models for California and florida pixel level forecasts
#This is a forecast validation script, not a prdiction file
#Created: 11/4/24
#Author: Max Beal

#library(brms, lib.loc = "/home/local-rhel8/apps/R-4.4.0/intel-23.2/lib64/R/library")
library(brms)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
#library(rstan, lib.loc = "/home/local-rhel8/apps/R-4.4.0/intel-23.2/lib64/R/library")
library(future)
library(data.table)
library(lubridate)

k_fold <- function(data, k) {
  # Create folds
  set.seed(123)  # For reproducibility
  folds <- sample(rep(1:k, length.out = nrow(data)))
  return(folds)
}




#### Read in the data frame ####

df <- fread("/work/HAB4CAST/max_beal/CyANPixelForecast/data/florida_data_spatial.csv")

df$chl = df$chl + 0.0001 #add to remove infinity values when logging
df$logchl = log(df$chl)

#Create the week-year column, make sure weeks have 2 digits
df$Time = paste0(df$Year,str_pad(df$Week,2,pad=0))
morpho=read.csv( "/work/HAB4CAST/max_beal/CyANPixelForecast/data/static_morphological_characteristics.csv") %>% dplyr::select(-X,-Area)
df=merge(df,morpho, by="COMID")

# Make the example reproducible
set.seed(1)



input <- df %>%
  arrange(COMID,x,y,Time) %>%
  # Fill Wtemp data
  group_by(COMID,x,y) %>%
  mutate(logchl_lead = lead(logchl, n = 1), #SET LEAD TIME HERE FOR CYAN
         logchl_lag2 = lag(logchl, n = 1),
         logchl_lag3 = lag(logchl, n = 2),
         logchl_lag4 = lag(logchl, n = 3)
  ) %>%
  fill(PRISMARDSW,PRISMSW) %>%
  ungroup() %>%   # # Pre-process Bloom data
  filter(!is.na(logchl_lead))


#Adjust sw temp to use PRISMARDSW when available, PRISMARD as a backup
#Note PRISMARDSW now represents all surface water predictors
input = input %>% filter(input$PRISMARDSW %>% is.na()) %>% mutate(PRISMARDSW = PRISMSW) 

padded_numbers <- sprintf("%06d", input$time)
input$week = substr(padded_numbers,1,2)
input$year = substr(padded_numbers,3,6)
dft=input %>%
  mutate(beginning = ymd(str_c(year, "-01-01")), date = beginning + weeks(week)) %>% 
  dplyr::select(COMID,date,time,x,y,logchl_lead,logchl,PRISMARDSW,precip,surfaceRadiation,windSpeed,dMean,Week,mean_neighbor_chl)

dft = input
rm(input)

#### Perform time blocking ####
dft$subset = year(dft$date) - (min(year(dft$date)) - 1)

#Verify time blocking
verify = dft %>% group_by(subset) %>% summarise(start_time=min(date),
                                                end_time=max(date),
                                                n_obs = n())


#### Define priors ####

imu=mean(dft$logchl,na.rm=T)
isigma=sd(dft$logchl,na.rm=T)

priors = c(prior(normal(mu,sigma),class="Intercept"))

#### Bayesian MCMC model with BRMS ####
library(future)

# Function to perform k-fold cross-validation

#plan(multicore,workers = availableCores())
print(availableCores())
cores = availableCores()

for (i in c(3:max(verify$subset))) {
  
  
  #Subset and scale/center training set
  train = dft %>% filter(subset <= i-1) %>% na.omit()
  
  #Try sequential training
  # train = dft %>% filter(subset==(i-1)) %>% na.omit() 
  train_predictor_scaled = train %>% dplyr::select(logchl,logchl_lag2,logchl_lag3,logchl_lag4,PRISMARDSW, precip, surfaceRadiation, windSpeed,mean_neighbor_chl) %>% scale(center=TRUE,scale=TRUE)
  
  #Rejoin scaled predictors and predictand in training set
  train = data.frame("logchl_lead"=train$logchl_lead,"Time"=train$time,"x"=train$x,"y"=train$y,"COMID"=train$COMID,"dMean"=train$dMean,"Week"=train$Week,train_predictor_scaled)
  
  #Subset and scale/center test set
  test = dft %>% filter(subset==i) %>% na.omit()
  test_predictor_scaled = test %>% dplyr::select(logchl,logchl_lag2,logchl_lag3,logchl_lag4,PRISMARDSW, precip, surfaceRadiation, windSpeed,mean_neighbor_chl) %>% scale(center=TRUE,scale=TRUE)
  
  #Rejoin scaled predictors and predictand in testing set
  test = data.frame("logchl_lead"=test$logchl_lead,"Time"=test$time,"x"=test$x,"y"=test$y,"COMID"=test$COMID,"dMean"=test$dMean,"Week"=test$Week,test_predictor_scaled)
  
  #Calculate the null model now based on the train informaton we've "seen" before the prediction takes place
  lta = dft %>% filter(subset<=i-1) %>% group_by(Week,x,y) %>% summarise(logchl_LTA = mean(logchl_lead,na.rm=T))
  lta_xy = dft %>% filter(subset<=i) %>% group_by(x,y) %>% summarise(logchl_LTA_xy = mean(logchl_lead,na.rm=T))
  
  
  #Join with the test set
  test = left_join(test,lta,by=c("Week","x","y"))
  test = left_join(test,lta_xy,by=c("x","y"))
  test=test %>% mutate(logchl_LTA = ifelse(is.na(logchl_LTA),logchl_LTA_xy,logchl_LTA)) %>% dplyr::select(-logchl_LTA_xy)
  
  
  #Threads
  threads_per_chain = 32
  Sys.setenv(STAN_NUM_THREADS = threads_per_chain)
  
  #### Define formula and run model ####
  #Note: can't run autoregression from brms on a thredded model, so implementing mannually
  formula = logchl_lead ~  PRISMARDSW + precip + surfaceRadiation + windSpeed + dMean + mean_neighbor_chl + logchl + logchl_lag2 + logchl_lag3 + logchl_lag4
  
  brm.1 <- brm(formula, 
               
               family= gaussian(), 
               
               data = train, 
               
               chains = 4, #specify the number of Markov chains,
               
               cores = 4,
               
               threads = threads_per_chain,
               
               iter = 3000, warmup = 1000) 
  
  
  #post <- as_draws(brm.1)
  
  ensemble=posterior_predict(brm.1,newdata=test,ndraws=100,allow_new_levels=TRUE,cores=4)
  
  
  numbers = c(1:100)
  strings <- str_c("enschl_member_", numbers)


  ensemble = t(ensemble)
  colnames(ensemble) = strings
  
  test=cbind(test,exp(ensemble))
  
  test$logchl_predicted = rowMeans(ensemble)
  
  test$error = test$logchl_lead - test$logchl_predicted
  
  test = test %>% mutate( chl_predicted = exp(logchl_predicted),
                          chl_observed = exp(logchl_lead),
                          Bloom_predicted = chl_predicted>=12,
                          Bloom_observed = chl_observed>=12)
  
  
  write_csv(test,paste0("/work/HAB4CAST/max_beal/CyANPixelForecast/output/FL_BRMS_fold_",i,"_lag1.csv"))
  
  
  # fixed_ef = fixef(brm.1)
  # rand_ef = ranef(brm.1)
  # 
  # week_rf=as.data.frame(rand_ef$Week)
  # week_rf$variable = "Week"
  # 
  # write_csv(week_rf,paste0("/work/HAB4CAST/max_beal/CyANPixelForecast/output/FL_BRMS_random_effect_fold_",i,"_lag1.csv"))
  # 
  # fixed_ef = as.data.frame(fixed_ef)
  # write_csv(fixed_ef,paste0("/work/HAB4CAST/max_beal/CyANPixelForecast/output/FL_BRMS_fixed_effect_fold_",i,"_lag1.csv"))
  
  gc()
  

}



