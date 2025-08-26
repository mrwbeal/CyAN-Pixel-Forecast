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
library(abind)


#### data ####

CA = F
FL = T
if (CA) {
  df <- read_csv("/work/HAB4CAST/max_beal/CyANPixelForecast/data/california_data.csv")
}

if (FL) {
  df <- fread("/work/HAB4CAST/max_beal/CyANPixelForecast/data/florida_data_spatial.csv")
  
}



#### Functions ####
# Function to perform k-fold cross-validation
k_fold <- function(data, k) {
  # Create folds
  folds <- sample(rep(1:k, length.out = nrow(data)))
  return(folds)
}




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

hist(input$logchl)

#Adjust sw temp to use PRISMARDSW when available, PRISMARD as a backup
#Note PRISMARDSW now represents all surface water predictors
input = input %>% filter(input$PRISMARDSW %>% is.na()) %>% mutate(PRISMARDSW = PRISMSW) 



padded_numbers <- sprintf("%06d", input$time)
input$week = substr(padded_numbers,1,2)
input$year = substr(padded_numbers,3,6)


dft=input %>%
  mutate(beginning = ymd(str_c(year, "-01-01")), date = beginning + weeks(week)) %>% 
  dplyr::select(COMID,date,time,x,y,chl_lead,chl,chl_lag2,chl_lag3,chl_lag4,PRISMARDSW,precip,surfaceRadiation,windSpeed,dMean,Week,mean_neighbor_chl)



#### Perform time blocking ####
#### Try subsetting buy year ####
dft$subset = year(dft$date) - (min(year(dft$date)) - 1)

#Verify time blocking
verify = dft %>% group_by(subset) %>% summarise(start_time=min(date),
                                        end_time=max(date),
                                        n_obs = n())

verify


#Number of observations per pixel
numobs = dft %>% group_by(x,y) %>% summarise(n=length(COMID))
numobs %>% arrange(desc(n))
hist(numobs$n)



subset_index = dft %>% group_by(date) %>% summarise(subset = mean(subset)) %>% pull(subset)


t = dft%>% dplyr::select(date) %>% unique() %>% arrange(date)
coef_save = data.frame()
comids = dft$COMID %>% unique()

#### Choose model type ####
ME = FALSE
RF = FALSE
LSTM=TRUE #Need to run LSTM with GPU, with a high memory request

#Get X Y data if LSTM
if (LSTM) {
  x=readRDS("/work/HAB4CAST/max_beal/CyANPixelForecast/data/convLSTM_train_test/convlstm_x_dim.rds")
  y=readRDS("/work/HAB4CAST/max_beal/CyANPixelForecast/data/convLSTM_train_test/convlstm_y_dim.rds")
  lta=readRDS("/work/HAB4CAST/max_beal/CyANPixelForecast/data/convLSTM_train_test/convlstm_lta_dim.rds")
  

  
  x=aperm(x,perm=NULL,c(dim(x)[3],dim(x)[2],dim(x)[1]))
  y=aperm(y,perm=NULL,c(dim(y)[3],dim(y)[2],dim(y)[1]))
  lta=aperm(lta,perm=NULL,c(dim(lta)[3],dim(lta)[2],dim(lta)[1]))
  x=flip(rast(x))
  y=flip(rast(y))
  lta=flip(rast(lta))
  
  
  folddates = t$date
  names(x) = folddates
  names(y) = folddates
  names(lta) = folddates
  
  xdf = as.data.frame(x,xy=T)
  xdf = melt(xdf,id.vars=c("x","y")) %>% na.omit()
  
  ydf = as.data.frame(y,xy=T)
  ydf = melt(ydf,id.vars=c("x","y")) %>% na.omit()
  
  ltadf = as.data.frame(lta,xy=T)
  ltadf = melt(ltadf,id.vars=c("x","y")) %>% na.omit() #LTA will have NAN up  until subset 3 which starts in 2019
  
  xydf = merge(xdf,ydf,by=c("x","y","variable")) # maybe this is where data is going missing?
  
  xyltadf = left_join(xydf,ltadf,by=c("x","y","variable"))

  colnames(xyltadf) = c("x","y","date","realx","realy","chl_LTA")
  #rm(x,y,xdf,ydf)
  gc()
}

for (i in c(3:max(verify$subset))) {
  


  #### Define formula and train model ####
 
  
  if (LSTM) {
    
    #Create data generator to feed ConvLSTM
    data_generator <- function(file_path_x,file_path_y, time_steps, fold, subset_index = subset_index, NAvalue=NAvalue){
      #Open file connection
      xdata = readRDS(file_path_x)
      ydata = readRDS(file_path_y)
      
      print(dim(xdata))
      #Subset the correct training data
      mask = subset_index<=fold-1
      
      xdata = xdata[mask,,,]
      ydata = ydata[mask,,]
      
      #dim(xdata) = c(t_dim,x_dim,y_dim,p_dim)
      
      total_time <- dim(xdata)[1]
      current_time = 1 #initialize time
      
      print(dim(ydata))
      
      function() {
        if(current_time>total_time - time_steps + 1){
          current_time<<-1 #Restart for new epoch
        }
      
        
      if (time_steps==1) {
          #deterine the range of the time steps to read in at a time
          end_time = min(current_time, total_time)
          
          
          #Extract block of time steps
          x_train_matrix <- xdata[current_time, , ,] 
          y_train_matrix <- ydata[current_time, , ] 
          
          #Reshape and establish NA value
          #Make sure these dimensions are in the right order
          t_dim = 1
          x_dim = dim(x_train_matrix)[1]
          y_dim = dim(x_train_matrix)[2]
          p_dim = dim(x_train_matrix)[3]
          
      }
      
        if (time_steps>1) {
          #deterine the range of the time steps to read in at a time
          end_time = min(current_time + time_steps-1, total_time)
          
          
          #Extract block of time steps
          x_train_matrix <- xdata[current_time:end_time, , ,] #exclude 7 which is the subset layer
          y_train_matrix <- ydata[current_time:end_time, , ] #Exclude 2 which is the subset layer
          
          #Reshape and establish NA value
          #Make sure these dimensions are in the right order
          t_dim = dim(x_train_matrix)[1]
          x_dim = dim(x_train_matrix)[2]
          y_dim = dim(x_train_matrix)[3]
          p_dim = dim(x_train_matrix)[4]
          
        }

      
      dim(x_train_matrix) = c(1,t_dim, x_dim, y_dim, p_dim)
      dim(y_train_matrix) = c(1,t_dim,x_dim,y_dim,1)
      
      x_train_matrix[is.na(x_train_matrix)] = NAvalue
      y_train_matrix[is.na(y_train_matrix)] = NAvalue
      
      print(dim(x_train_matrix))
      
      list(x_train_matrix,y_train_matrix)
      
      #any(is.nan(x_train_matrix)) #Make sure no NaNs
      }
    }
    
    #Need to reshape data in format time, x, y, preds

    preds = c("chl","PRISMARDSW","precip","surfaceRadiation","windSpeed","dMean")
    
    n_features = length(preds)
    
    names = list.files("/work/HAB4CAST/max_beal/CyANPixelForecast/data/convLSTM_train_test/",full.names = T)
    
    file_x = names[grepl("train_x",names)]
    file_y = names[grepl("train_y",names)]
    

    
    NAvalue = -99
    
    # Define the custom loss function
    mse_custom_loss <- function(y_true, y_pred) {
      mask = k_cast(k_not_equal(y_true, y_true[1,1,1,1,1]),"float32") ## Grabs a value from the corner
      square_error = k_cast(k_square(y_true-y_pred),"float32")
      masked_square_error=square_error*mask
      loss = k_sum(masked_square_error)/k_sum(mask) # Mean Squared Error
      return(loss)
    }
    
    
    gen<-data_generator(file_path_x=file_x,file_y,time_steps = 4,fold=i,subset_index = subset_index,NAvalue=NAvalue)
    gc()
    #batch = gen() #Test the generator

    
    t_dim = 4 #Time chunks
    x_dim = 775
    y_dim = 1000
    p_dim = 6
    
    total_time = 388 #Total number of time steps in the fold
    
    
    model <- keras_model_sequential() %>%
      layer_masking(mask_value=NAvalue) %>%
      layer_conv_lstm_2d(filters=32, kernel_size=c(4,4), #4x4 based on autocor analysis
                         input_shape = c(NULL,x_dim,y_dim,p_dim),padding="same",activation="tanh",
                         return_sequences = T) %>%
      layer_dropout(0.3) %>%
      layer_conv_lstm_2d(filters=32, kernel_size=c(4,4), #4x4 based on autocor analysis
                         input_shape = c(NULL,x_dim,y_dim,p_dim),padding="same", activation="tanh",
                         return_sequences = T) %>%
      layer_conv_2d(
        filters=1,kernel_size=c(3,3),
        padding="same", activation="softplus"
      )

      

    
    # Compile the model
    model %>% compile(
      optimizer = optimizer_adam(learning_rate=0.001),
      loss = mse_custom_loss,
      metrics = c("mae")
    )
    
    #Train the model
    history <- model %>% fit(
      x=gen,
      steps_per_epoch = ceiling(total_time/t_dim),
      epochs = 25
    )
    

    gc()
    
    x_test_matrix = readRDS(file_x)
    y_test_matrix = readRDS(file_y)
    
    
    #Subset the correct testing data
    mask = subset_index==i
    
    #Variables
    x_test_matrix= x_test_matrix[mask,,,]
    y_test_matrix = y_test_matrix[mask,,]
    
    
    
    dim(x_test_matrix) = c(1,dim(x_test_matrix)[1], x_dim, y_dim, p_dim)
    dim(y_test_matrix) = c(1,dim(y_test_matrix)[1],x_dim,y_dim,1)
    
    x_test_matrix[is.na(x_test_matrix)] = NAvalue
    y_test_matrix[is.na(y_test_matrix)] = NAvalue
    
    gc()
    
    partition_predictions=TRUE
    if (partition_predictions) {

      #Necessary to keep within memory limits for wider (more filter) models
      tmat1 = x_test_matrix[,1:26,,,]
      dim(tmat1) = c(1,dim(tmat1)[1], x_dim, y_dim, p_dim)
      predictions1 <- model %>% predict(tmat1)
      
      tmat2 = x_test_matrix[,27:52,,,]
      dim(tmat2) = c(1,dim(tmat2)[1], x_dim, y_dim, p_dim)
      predictions2 <- model %>% predict(tmat2)
    
      predictions=abind(predictions1,predictions2,along=2)
      print(dim(predictions))
      
    }
    
    if (!partition_predictions) {
      # Generate predictions
      predictions <- model %>% predict(x_test_matrix)
    }
    
    
    
    #Get predictions and observations into raster format and compare
    allp = predictions[1,,,,1]
    allo = y_test_matrix[1,,,,1]
    
    allp=aperm(allp,perm=NULL,c(dim(allp)[3],dim(allp)[2],dim(allp)[1]))
    allo=aperm(allo,perm=NULL,c(dim(allo)[3],dim(allo)[2],dim(allo)[1]))
    gc()
    
    rp = flip(rast(allp))
    ro = flip(rast(allo))
    
    # plot(ro[[8]])
    # plot(rp[[8]])
    #ext(r) = c(min(x,na.rm=T),max(x,na.rm=T),min(y,na.rm=T),max(y,na.rm=T))

    #Remove NA values
    rp[rp<=as.numeric(rp[[1]][1])] = NA
    ro[ro==NAvalue] = NA
    
    
    folddates = t$date[mask]
    names(rp) = folddates
    names(ro) = folddates
    
    pdf = as.data.frame(rp,xy=T)
    pdf = melt(pdf,id.vars=c("x","y"))
    colnames(pdf) = c("x","y","date","chl_predicted")
    
    odf = as.data.frame(ro,xy=T)
    odf = melt(odf,id.vars=c("x","y"))
    colnames(odf) = c("x","y","date","chl_lead")
    

    test = merge(odf,pdf,on=c("x","y","date"))
    test = left_join(test,xyltadf,by=c("x","y","date"))

    test = test %>% na.omit() #is this where data is going missing?
    rm(pdf,odf)
    

    
    #Get dates
    test$date = test$date %>% as.Date()
    test$Week= test$date %>% week()
    
    gc()
    
    
    library(MASS)
    test$error = test$chl_lead - test$chl_predicted
    
    # Get ensemble members from normal distribution prediction
    norm_error_distributions = test  %>% summarise("mean"=fitdistr(error,"normal")$estimate['mean'],
                                                   "sd"=fitdistr(error,"normal")$estimate['sd'])
    #Pull from cell-specific error ensemble
    create_ensemble = function (preds) {preds + rnorm(100,mean=norm_error_distributions$mean,sd=norm_error_distributions$sd)}
    
    #Create ensemble members
    preds = test %>% dplyr::select(chl_predicted)
    ensemble = apply(preds, 1, create_ensemble) %>% t()
    
    numbers = c(1:100)
    strings <- str_c("enschl_member_", numbers)
    colnames(ensemble) = strings
    
  }
  
  
  
  if (!LSTM) {
    test$logchl_predicted = pred
  }
  
  
  #test$error = test$chl_lead - test$chl_predicted
  
  test = test %>% mutate(chl_predicted = chl_predicted,
                            chl_observed = chl_lead,
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
  
  #ggplot(test) + geom_tile(aes(x=x,y=y,fill=logchl_predicted))
  a=ggplot(test) + geom_point(aes(x=chl_predicted,y=chl_observed)) + coord_equal() + geom_abline(slope = 1,intercept=0)
  ggsave("clstm.png",plot=a,width=8,height=8,dpi=300)

  
  if (FL) {
    #Save kfold results
    write_csv(test,paste0("/work/HAB4CAST/max_beal/CyANPixelForecast/output/FL_ConvLSTM_fold_",i,"_lag1.csv"))
  }
  
  gc()
  
}



