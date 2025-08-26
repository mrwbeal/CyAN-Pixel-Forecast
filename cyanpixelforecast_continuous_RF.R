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
  
  #df = df %>% rename(Cyano = value, Year = year, Bloom=bloom)
  
}

# test = df %>% filter(COMID==166757656) %>% select(PRISMSW,x,y)
# 
# ggplot(test) + geom_tile(aes(x=x,y=y,fill=PRISMARDSW))


#### Functions ####
# Function to perform k-fold cross-validation
k_fold <- function(data, k) {
  # Create folds
  folds <- sample(rep(1:k, length.out = nrow(data)))
  return(folds)
}



time_blocking <- function(data,time_column,block_size){
  
  data= data %>% arrange(date) %>% mutate(subset = as.numeric(cut(as.numeric(factor(.[[time_column]])),
                                                                      breaks=seq(1,length(unique(.[[time_column]])) + block_size, by=block_size),
                                                                      labels=FALSE,
                                                                      include.lowest = TRUE)))
}




# cicyano = 10^(df$Cyano*0.011714-4.1870866)
# df$chl = 6620*cicyano - 3.07
# df = df %>% filter(chl>0)
# df$logchl = log(df$chl) #Make it normal



#### Run Model ####
#Create the week-year column, make sure weeks have 2 digits
df$Time = paste0(df$Year,str_pad(df$Week,2,pad=0))
morpho=read.csv( "/work/HAB4CAST/max_beal/CyANPixelForecast/data/morphological_characteristics.csv")
df=merge(df,morpho, by="COMID")

# Make the example reproducible
set.seed(1)


input <- df %>%
  arrange(COMID,x,y,Time) %>%
  # Fill Wtemp data
  group_by(COMID,x,y) %>%
  mutate(chl_lead = lead(chl, n = 1), #SET LEAD TIME HERE FOR CYAN, try
         chl = chl,
         chl_lag2 = lag(chl, n = 1),
         chl_lag3 = lag(chl, n = 2),
         chl_lag4 = lag(chl, n = 3)
  ) %>%
  fill(PRISMARDSW,PRISMSW) %>%
  ungroup() %>%   # # Pre-process Bloom data
  filter(!is.na(chl_lead))

hist(input$chl)

#Adjust sw temp to use PRISMARDSW when available, PRISMARD as a backup
#Note PRISMARDSW now represents all surface water predictors
input = input %>% filter(input$PRISMARDSW %>% is.na()) %>% mutate(PRISMARDSW = PRISMSW) 



padded_numbers <- sprintf("%06d", input$time)
input$week = substr(padded_numbers,1,2)
input$year = substr(padded_numbers,3,6)


dft=input %>%
  mutate(beginning = ymd(str_c(year, "-01-01")), date = beginning + weeks(week)) %>% 
  dplyr::select(COMID,date,time,x,y,chl_lead,chl,chl_lag2,chl_lag3,chl_lag4,PRISMARDSW,precip,surfaceRadiation,windSpeed,dMean,Week,mean_neighbor_chl,year)



#### Try subsetting buy year ####
#dft$subset = year(dft$date) - (min(year(dft$year)) - 1)

dft$subset = as.numeric(factor(dft$year))

#Verify time blocking
verify = dft %>% group_by(subset) %>% summarise(start_time=min(date),
                                        end_time=max(date),
                                        n_obs = n())


verify
#Number of observations per pixel
numobs = dft %>% group_by(x,y) %>% summarise(n=length(COMID))
numobs %>% arrange(desc(n))
hist(numobs$n)

#verify$year = verify$start_time %>% year()



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
  
  
  # one_to_one<-ggplot(test) + geom_point(aes(x=logchl_lead,y=logchl_predicted)) + geom_abline(slope=1,intercept=0,linetype="dashed",color="red") + theme_bw()
  # print(one_to_one)
  
  res <- caret::postResample(test$chl_lead, test$chl_predicted)
  print(res[2])
  
  
  tbl = table(test$Bloom_predicted,test$Bloom_observed)
  print(tbl)
  accuracy = (tbl[1,1]+tbl[2,2])/(tbl[1,1]+tbl[1,2]+tbl[2,1]+tbl[2,2])
  print(accuracy)
  

  if (FL) {
    #Save kfold results
    write_csv(test,paste0("/work/HAB4CAST/max_beal/CyANPixelForecast/output/FL_RF_fold_",i,"_lag1.csv"))
  }

    # #Save coefficients
    write_csv(coef_save,"/work/HAB4CAST/max_beal/CyANPixelForecast/output/FL_RF_stdcoefs.csv")
    # 
    coef_save=read_csv("/work/HAB4CAST/max_beal/CyANPixelForecast/output/FL_RF_stdcoefs.csv")
    flcoef = melt(coef_save,id.vars="kfold") %>% filter(variable!="X.Intercept." & variable!="bootstrap")
    ggplot(flcoef) + geom_boxplot(aes(x=variable,y=value,fill=kfold))

}

fileresults=list()
for (i in c(1:7)) {
  fileresults[[i]]=paste0("/work/HAB4CAST/max_beal/CyANPixelForecast/output/RF_folds/FL_RF_fold_",i,"_lag1.csv")
}




folds = lapply(fileresults,fread) %>% bind_rows()


# Evaluation metrics
r2 = function(preds,obs){
  ssr = sum((obs-preds)^2)
  sst = sum((obs - mean(obs,na.rm=T))^2)
  R2 = 1 - (ssr/sst)
  return(R2)
}




evaluate_continuous = function(df){
  
  
  clta = df %>% group_by(Week, x, y) %>% mutate(chl_LTA = mean(chl_observed,na.rm=T),
                                                Bloom_LTA = as.numeric(chl_LTA>=12),
                                                Bloom_observed = as.numeric(chl_observed>=12),
                                                Bloom_predicted = as.numeric(chl_predicted>=12),
                                                PctBloom=mean(Bloom_observed,na.rm=T),
                                                variability = PctBloom>0.1 & PctBloom<0.9) %>% na.omit()
  
  #Create factors from categorical data (for comparison with original model)
  clta$Bloom_observed=factor(clta$Bloom_observed, levels = c(0,1),ordered=T)
  clta$Bloom_predicted=factor(clta$Bloom_predicted, levels = c(0,1),ordered=T)
  clta$Bloom_LTA=factor(clta$Bloom_LTA, levels = c(0,1),ordered=T)
  
  # Define Metrics
  metrics = clta %>% group_by(Week,variability) %>% summarise('Accuracy'=confusionMatrix(Bloom_observed,Bloom_predicted)$overall['Accuracy'],
                                                                   'Sensitivity'=confusionMatrix(Bloom_observed,Bloom_predicted)$byClass['Sensitivity'],
                                                                   'Specificity'=confusionMatrix(Bloom_observed,Bloom_predicted)$byClass['Specificity'],
                                                                   'AccuracyNull'=confusionMatrix(Bloom_observed,Bloom_LTA)$overall['Accuracy'],
                                                                   'SensitivityNull'=confusionMatrix(Bloom_observed,Bloom_LTA)$byClass['Sensitivity'],
                                                                   'SpecificityNull'=confusionMatrix(Bloom_observed,Bloom_LTA)$byClass['Specificity'],
                                                                   'Rsquared' = r2(chl_predicted,chl_observed),
                                                                   'MAE'=MAE(chl_predicted,chl_observed),
                                                                   'MAENull'=MAE(chl_LTA,chl_observed))
  
  
  metrics = metrics %>% mutate(Acc_imp = Accuracy-AccuracyNull,
                               Sen_imp = Sensitivity-SensitivityNull,
                               Spe_imp = Specificity-SpecificityNull,
                               MAE_imp = MAE-MAENull)
  
  return(metrics)
}

plots = function(metrics){
  
  #Difference
  a =ggplot()+geom_line(data=metrics, aes(x=Week,y=Acc_imp,color=variability)) + theme_classic() + ylab("Accuracy Improvement") + 
    geom_hline(yintercept = 0, col="red",linetype="dashed")+ ggtitle("Forecast Improvement \nover Climatological Null")+ ylim(-0.5,0.5)
  b = ggplot()+geom_line(data=metrics, aes(x=Week,y=Sen_imp,color=variability)) + theme_classic() + ylab("Sensitivity Improvement")+ 
    geom_hline(yintercept = 0, col="red",linetype="dashed")+ ylim(-0.5,0.5)
  c = ggplot()+geom_line(data=metrics, aes(x=Week,y=Spe_imp,color=variability)) + theme_classic() + ylab("Specificity Improvement")+ 
    geom_hline(yintercept = 0, col="red",linetype="dashed") + ylim(-1,0.5)
  e = ggplot()+geom_line(data=metrics, aes(x=Week,y=MAE_imp,color=variability)) + theme_classic() + ylab("MAE Improvement") + 
    geom_hline(yintercept = 0, col="red",linetype="dashed")
  
  #Forecast Results
  af =ggplot()+geom_line(data=metrics, aes(x=Week,y=Accuracy,color=variability)) + theme_classic() + ylab("Accuracy") + ggtitle("Forecast Performance") +
    ylim(0,1)
  bf = ggplot()+geom_line(data=metrics, aes(x=Week,y=Sensitivity,color=variability)) + theme_classic() + ylab("Sensitivity")+
    ylim(0,1)
  cf = ggplot()+geom_line(data=metrics, aes(x=Week,y=Specificity,color=variability)) + theme_classic() + ylab("Specificity")+
    ylim(0,1)

    #Continuous
  df = ggplot()+geom_line(data=metrics, aes(x=Week,y=Rsquared,color=variability)) + theme_classic() + ylab("R-squared")+
    ylim(0,1)
  ef = ggplot()+geom_line(data=metrics, aes(x=Week,y=MAE,color=variability)) + theme_classic() + ylab("MAE")
  
  #Null Results  
  an =ggplot()+geom_line(data=metrics, aes(x=Week,y=AccuracyNull,color=variability)) + theme_classic() + ylab("Accuracy") + ggtitle("Null (LTA) Performance")+
    ylim(0,1)
  bn = ggplot()+geom_line(data=metrics, aes(x=Week,y=SensitivityNull,color=variability)) + theme_classic() + ylab("Sensitivity")+
    ylim(0,1)
  cn = ggplot()+geom_line(data=metrics, aes(x=Week,y=SpecificityNull,color=variability)) + theme_classic() + ylab("Specificity")+
    ylim(0,1)
  
    #Continuous
  en = ggplot()+geom_line(data=metrics, aes(x=Week,y=MAENull,color=variability)) + theme_classic() + ylab("MAE Null")
  

  
  library(ggpubr)
  categorical = ggarrange(a,af,an,b,bf,bn,c,cf,cn,align = "v")
  continuous = ggarrange(e,ef,en,df)

  print(categorical)
  print(continuous)
}

metrics=evaluate_continuous(folds)

plots(metrics)

metrics


ggplot(folds, aes(x=chl_observed, y=chl_predicted) ) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()








# ## Sample from error distribution
# library(MASS)
# 
# norm_error_distributions = test %>% group_by(x,y) %>% summarise("mean"=fitdistr(error,"normal")$estimate['mean'],
#                                                                  "sd"=fitdistr(error,"normal")$estimate['sd'])
# for (i in c(1:nrow(norm_error_distributions))) {
#   x = norm_error_distributions$x[i]
#   y = norm_error_distributions$y[i]
# 
#   #Pull from cell-specific error ensemble
#   create_ensemble = function (preds) {preds + rnorm(100,mean=norm_error_distributions$mean[i],sd=norm_error_distributions$sd[i])}
# 
#   #Create ensemble members
#   preds = test %>% filter(x==x, y==y) %>% dplyr::select(logchl_predicted)
# 
#   ensemble = apply(preds, 1, create_ensemble) %>% t()
# 
# 
#   boxplot(t(ensemble[c(1:10),]))
#   lines(c(1:10),preds$logchl_predicted[1:10],type="l")
# 
#   dim(ensemble)
# 
# 
# }

# 
# folds = folds %>% group_by(COMID) %>% mutate(extprecip = as.numeric(precip>quantile(precip,0.9)),
#                                              extsw = as.numeric(PRISMARDSW>quantile(PRISMARDSW,0.9))) %>% ungroup()
# 
# 
# #### Correlations ####
# correlations = folds %>%
#   group_by(x,y,Week) %>%
#   summarize(swtempcor = cor(logchl_lead, PRISMARDSW),
#             pcor = cor(logchl_lead, precip),
#             srcor = cor(logchl_lead, surfaceRadiation),
#             wscor = cor(logchl_lead, windSpeed),
#             autocor = cor(logchl_lead, logchl),
#             extpcor = cor(logchl_lead,extprecip),
#             extsw = cor(logchl_lead,extsw))
# 
# cors = melt(correlations,id=c("x","y","Week"))
# 
# mcor=cors %>% group_by(variable,Week) %>% summarise(meancor=mean(value,na.rm=T))
# 
# ggplot(data=mcor) + geom_line(aes(x=Week,y=meancor,group=Week),group=1) + geom_abline(slope=0,intercept=0,linetype="dashed",color="red")+
#   facet_wrap(~variable) + theme_bw()
# 
# ggplot(data=cors) + geom_boxplot(data=cors,aes(x=Week,y=value,group=Week)) + geom_abline(slope=0,intercept=0,linetype="dashed",color="red")+
#   facet_wrap(~variable) + theme_bw()
# 
# ggplot(data=cors) + geom_density(data=cors,aes(x=value,group=variable,fill=variable),alpha=0.5) + geom_vline(xintercept=0,linetype="dashed",color="red") + theme_bw()
# 
# 
# correlations = folds %>%
#   summarize(swtempcor = cor(logchl_lead, PRISMARDSW),
#             pcor = cor(logchl_lead, precip),
#             srcor = cor(logchl_lead, surfaceRadiation),
#             wscor = cor(logchl_lead, windSpeed),
#             autocor = cor(logchl_lead, logchl),
#             extpcor = cor(logchl_lead,extprecip))
# correlations


