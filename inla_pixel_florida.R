#This notebook brings in the Florida COMIDs for INLA construction
#Attempts ot build a mesh for lakes in Florida
#This should not be used for operational prediction

#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)

#### Libraries ####
library(readr)
library(INLA)
library(dplyr)
library(ggplot2)
library(sf)
library(tidyverse)
library(inlabru)
library(ROCR)
library(caret)
library(RColorBrewer)
library(ggpubr)
library(colorspace)
library(lubridate)
library(grid)
library(gridExtra)
library(cowplot)
library(stringr)
library(tigris)
library(future)
library(future.apply)



#### Read in the data frame ####

cyan_files = list.files("/work/HAB4CAST/max_beal/CyANPixelForecast/data/cyan_pixel_conus_lakes",full.names = TRUE)
COMIDs = as.numeric(str_extract(substr(cyan_files,15,150),"\\d+"))
#check already processed files
prediction_files = list.files("/work/HAB4CAST/max_beal/CyANPixelForecast/data/compiled_pixel_conus_lakes_current",full.names = TRUE)

#Read in lakes and subset florida
conus_lakes <- st_cast(st_read("/work/HAB4CAST/max_beal/SW_model/data/OLCI_resolvable_lakes_2022_09_08/OLCI_resolvable_lakes_2022_09_08.shp"), "MULTIPOLYGON")

states <- states(cb = TRUE)
fl = states %>% filter(STUSPS =="FL")
fl = st_transform(fl,st_crs(conus_lakes))
fl_lakes = st_crop(conus_lakes,fl)


#Plot Florida lakes
plot(fl_lakes$geometry) + plot(fl,add=T,col=NA)

#Concatenate florida data files
fl_comid = as.numeric(fl_lakes$COMID)
prediction_names = list.files("/work/HAB4CAST/max_beal/CyANPixelForecast/data/compiled_pixel_conus_lakes_current",full.names = F)
fl_files = prediction_files[parse_number(prediction_names) %in% fl_comid]



#Okechobee
#prediction_files[parse_number(prediction_names) %in% 166757656]



df = lapply(fl_files,read_csv) %>% bind_rows()

write_csv(df,"/work/HAB4CAST/max_beal/CyANPixelForecast/data/florida_data.csv")

#Create the week-year column, make sure weeks have 2 digits
df$Time = paste0(df$Year,str_pad(df$Week,2,pad=0))
morpho=read.csv( "/work/HAB4CAST/max_beal/CyANPixelForecast/data/static_morphological_characteristics.csv") %>% select(-X,-Area)
df=merge(df,morpho, by="COMID")


#### Define Functions ####

# Function to perform k-fold cross-validation
k_fold <- function(data, k) {
  # Create folds
  set.seed(123)  # For reproducibility
  folds <- sample(rep(1:k, length.out = nrow(data)))
  return(folds)
}



#Got rid of train/test scheme in favor of a k-fold cross validation approach for now. Need to set up a different pipeline for prediction that will use all data for training.
create_input_data = function(input, fold){

  input$group<- ifelse(input$subset==fold, 'test', 'train')
  
  predictors = input  %>% group_by(group) %>% select(Wtemp_adjusted,Precip) %>% mutate(Wtemp_adjusted = scale(Wtemp_adjusted,center=T),
                                                                           Precip = scale(Precip,center = T))
  
  input = input %>% select(COMID,Time,Year,Week,date.x,Bloom_adjusted,Cyano,Area,dMax,dMean,LAT,LONG,subset) %>% cbind(predictors)
  
  return(input)
}




#Create Mesh
create_mesh = function(input, conus_lakes){
  loc = input %>% select(LONG, LAT)
  
  #lake = conus_lakes %>% filter(COMID==unique(input$COMID))
  lake = conus_lakes
  
  #Crop lake to point extent
  pts = st_as_sf(loc,coords=c("LONG","LAT"),crs=st_crs(lake))
  lake = st_crop(lake,st_bbox(pts))
  
  #Create a triangle mesh based on initial point locations, specified or automatic boundaries, and mesh quality parameters.
  #Set cutoff to 300m for S3?
  mesh <- inla.mesh.2d(loc = loc, boundary = lake, max.edge=c(300,100000), cutoff=300)
  
  return(mesh)
}



create_inla_stack = function(input,mesh){
  
  #Create SPDE model  
  spde<- inla.spde2.matern(mesh)   
  
  #Constructs observation/prediction weight matrices for models based on inla.mesh and inla.mesh.1d objects
  loc = input %>% select(LONG, LAT)
  A <- inla.spde.make.A(mesh, loc=as.matrix(loc))
  
  #Predictors
  effects <- c("Wtemp_adjusted", "Precip", "Area", "dMean", "dMax", "Year", "Week", "LONG")
  
  #Replace test points with NAs so INLA can fill
  input$Bloom_adjusted[input$group=="test"] = NA
  
  #Generates a list of named index vectors for an SPDE model.
  idx<-inla.spde.make.index('s',mesh$n)
  stk <- inla.stack(data=list(Bloom=input$Bloom_adjusted),
                    A=list(A,1), #spatial matrix
                    effects=list(idx, list(input[,effects]))) #predictors
  return(stk)
}


fit_inla = function(data_stack,mesh) {
  
  #Create SPDE model  
  spde<- inla.spde2.matern(mesh)   
  
  hyprior=list(theta1 = list(prior="pc.prec", param=c(0.5, 0.05)),
               theta2 = list(prior="pc.cor1", param=c(0.1, 0.9)))
  
  #Model with a spatial component, temporal component by week, and partitioned into training and validation data
  formula = Bloom ~ -1 +Wtemp_adjusted + Precip + Area + dMean+ f(Week, model="ar1", hyper=hyprior) + f(s, model=spde)
  
  #Some INLA variable help: https://rdrr.io/github/andrewzm/INLA/man/inla.stack.html
  results <- inla(formula, 
                  family='binomial', #Likelihood family
                  control.family=list(link="logit"),
                  data=inla.stack.data(data_stack), #Extract data for inla call, and optionally join with other variables
                  control.predictor=list(A=inla.stack.A(data_stack), compute=TRUE, link=1), 
                  #Compute = indicates whether the marginal densities for the linear predictor should be computed
                  #Link = link function of the model
                  control.inla=list(int.strategy='eb'), 
                  #permits to specify a list of variables to obtain more accurate approximations or reduce the computational time
                  #'eb' is empirical Bayes -> uses only one integration point equal to the posterior mode of the hyperparameters
                  verbose=TRUE,
                  control.compute=list(dic=TRUE, cpo=TRUE,return.marginals.predictor=TRUE)) #Compute both the DIC and the Cross-validated predictive measures
  
  return(results)
}


transform_results = function(results,train,test){
  
  input = rbind(train,test)
  
  #Get marginal predicted values for each time step
  list_marginals <- results$marginals.fitted.values
  list_marginals = list_marginals[1:nrow(input)]
  marginals <- data.frame(do.call(rbind, list_marginals))
  marginals$Time <- rep(names(list_marginals)[1:nrow(input)],
                        times = sapply(list_marginals, nrow))
  marginals = list_marginals[1:nrow(input)]
  
  #Extract the linear predictor (log-odds), the emarginal function takes the expected value of the posterior distribution.
  #log_odds <- sapply(marginals, function(x) inla.emarginal(function(y) y, x))
  
  log_odds <- sapply(marginals, function(x) {
    tryCatch(
      inla.emarginal(function(y) y, x),
      error = function(e) NA  # Return NA or handle the error as needed
    )
  })
  
  # Transforming log-odds to probabilities
  probabilities <- exp(log_odds) / (1 + exp(log_odds))
  
  # Display the probabilities
  input$probability = probabilities
  
  return(input)
  
}

performance_metrics = function(Observed, Predicted){
  data.frame("Observed"=Observed,"Predicted"=Predicted) %>% na.omit() %>% summarise(
    TP = sum(Observed == 1 & Predicted == 1),
    TN = sum(Observed == 0 & Predicted == 0),
    FP = sum(Observed == 0 & Predicted == 1),
    FN = sum(Observed == 1 & Predicted == 0),
    proportion_negative_observed = sum(Observed==0)/sum(Observed==1|Observed==0),
    proportion_positive_observed = sum(Observed==1)/sum(Observed==1|Observed==0),
    sensitivity = TP / (TP + FN),
    specificity = TN / (TN + FP),
    accuracy = (TP + TN) / (TP + TN + FP + FN))
}



optimize_cutoff = function(data, optimize_metric="accuracy",plots=TRUE){
  cutoffs = seq(0,1,0.01)
  output = data.frame("TP"=NA,"TN"=NA,"FP"=NA,"FN"=NA,
                      "proportion_negative_observed"=NA,
                      "proportion_positive_observed"=NA,
                      "sensitivity"=NA,"specificity"=NA,"accuracy"=NA)
  for (cutoff in cutoffs) {
    Predicted = as.numeric(data$probability>=cutoff)
    Observed = data$Bloom_adjusted
    output=rbind(output,performance_metrics(Observed,Predicted))
  }
  output = output[-1,]
  output$cutoffs = cutoffs
  
  
  #Get cutoff for chosen metric, if multiple cutoffs with same score, choose the lowest.
  optimal_cutoff = output[output[,optimize_metric] == max(output[,optimize_metric]),] %>% select(cutoffs) %>% min()
  
  if (plots) {
    p<-ggplot() + 
      geom_line(data=output,aes_string(x="cutoffs",y=optimize_metric)) + 
      geom_vline(xintercept = optimal_cutoff,color="red")
    print(p)
  }
  
  best=output[output[,optimize_metric] == max(output[,optimize_metric]),]
  print(best[1,])
  
  return(optimal_cutoff)
}

#Generate Forecast performance stats
generate_stats <- function(test, cutoff) {
  
  test = test %>% filter(!is.na(probability)&!is.na(Bloom_adjusted))
  Observations = factor(test$Bloom_adjusted, levels=c(1,0),ordered=TRUE)
  Predictions = factor(as.numeric(test$probability>=cutoff),levels=c(1,0),ordered=TRUE)
  
  roc_instance <- prediction(test$probability,Observations)
  auc <- performance(roc_instance, measure = 'auc')
  cm = table(Predictions,Observations)
  stats <- confusionMatrix(Predictions,Observations)
  
  
  tibble(
    AUC = auc@y.values[[1]],
    Sensitivity = stats$byClass[1],
    Specificity = stats$byClass[2],
    Accuracy = stats$overall[1],
    Precision = cm[1,1]/(cm[1,1] + cm[1,2]),
    Prevalance = (cm[1,1] + cm[2,1])/sum(cm),
    `False Omission Rate` = cm[2,1]/(cm[2,1] + cm[2,2]),
    `F1 Score` = (2*cm[1,1])/(2*cm[1,1] + cm[1,2] + cm[2,1]),
    Kappa = stats$overall[2],
    `Brier Score` = mean((Predictions - Observations)^2)
  ) %>% t() %>% return()
}




#### Prep the input data for florida ####
## Weeks are defined to match the Time column in compiled data, YYYYWW
first_week <- 201701
mid_week <- 202101
final_week <- 202152

input <- df %>%
  # Fill Wtemp data
  group_by(COMID) %>%
  arrange(COMID,Time) %>%
  mutate(Cyano_adjusted = lead(Cyano),
         Bloom_adjusted = as.numeric(lead(Bloom)), #Lead pushes the predictand forward a week
         Wtemp_adjusted = ifelse(Wtemp < 0, NA, Wtemp)) %>%
  fill(Wtemp_adjusted) %>%
  ungroup() %>%
  
  # Keep only data that falls within the wtemp data
  filter(Time >= first_week & Time <= final_week) %>%
  
  # Pre-process Bloom data
  filter(!is.na(Bloom_adjusted))


#Assign groups for K fold cross validation
input$subset = k_fold(input,5) #Can maybe make this work better by assigning Nas instead of breaking data up?

#table(input$subset,input$Bloom_adjusted)

#### Run Functions ####

inla_kfold = function(foldNum,input,fl_lakes) {
  
  #Partition and scale train/test sets
  data = create_input_data(input,foldNum)

  #Create the spatial mesh (Needs whole data frame)
  mesh = create_mesh(data,fl_lakes)
  #plot(mesh)
  
  #Create Data stack
  data_stack = create_inla_stack(data,mesh)
  
  #Run INLA
  results = fit_inla(data_stack,mesh)
  
  #Next get probabilites, optimize cutoff, compare results, etc.
  out = transform_results(results,train,test)
  
  #Look at performance and find cutoff value
  test = out %>% filter(subset==foldNum)
  cutoff = optimize_cutoff(test)
  
  stats = generate_stats(test,cutoff)
  
  return(stats)
  #ggplot() + geom_boxplot(data=test, aes(x=as.factor(Bloom_adjusted),y=probability))
  
}



test = input %>% filter(COMID==84918 | COMID==166758657)
lapply(c(1:2),inla_kfold,input=test,fl_lakes=fl_lakes)

# plan(multisession, workers = (32))
# out = future_lapply(c(1:5),inla_kfold,input=input,fl_lakes=fl_lakes) %>% bind_rows()
# write_csv(out,"/work/HAB4CAST/max_beal/CyANPixelForecast/FL_results.csv")



#### Plots ####

# pr_ts = test %>% group_by(Time) %>% summarise(date = unique(date.x),
#                                               meanprob=mean(probability,na.rm=T),
#                                               maxprob=max(probability,na.rm=T),
#                                               minprob=min(probability),na.rm=T)
# 
# ggplot(data=pr_ts) + geom_line(aes(x=date,y=meanprob))+
#   geom_ribbon(data=pr_ts,aes(x=date,ymin=minprob,ymax=maxprob), fill="blue", alpha=0.1) + theme_classic()

# ggplot() + geom_tile(data=test,aes(x=LONG,y=LAT, fill=probability)) + 
#   facet_wrap(~Week)+ theme(axis.text.x = element_text(angle = 75, vjust = 0.5, hjust=1)) + xlab("") + ylab("") +
#   scale_fill_viridis_c()






