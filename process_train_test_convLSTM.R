#I am experiencing issues with reformatting the dataset to train the ConvLSTM using the GPU node
# The GPU node is necessary to run deep learning models like ConvLSTM
# I am creating this notebook to reformat the dataset into a 4d (x,y,time,predictor) train and test set
# This needs to be run on a large memory node and can then be added into the continuous training dataset

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


SelectByLake = FALSE


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



#Adjust sw temp to use PRISMARDSW when available, PRISMARD as a backup
#Note PRISMARDSW now represents all surface water predictors
input = input %>% filter(input$PRISMARDSW %>% is.na()) %>% mutate(PRISMARDSW = PRISMSW) 


# numFolds = 10
# input$subset = k_fold(input,numFolds)

padded_numbers <- sprintf("%06d", input$time)
input$week = substr(padded_numbers,1,2)
input$year = substr(padded_numbers,3,6)

dft=input %>%
  mutate(beginning = ymd(str_c(year, "-01-01")), date = beginning + weeks(week)) %>% 
  dplyr::select(COMID,date,Time,x,y,chl_lead,chl,chl_lag2,chl_lag3,chl_lag4,PRISMARDSW,precip,surfaceRadiation,windSpeed,dMean,Week)


plot(dft$chl[1:100],type="l")

#### Time, X, Y, Variables ####

#### Try subsetting be year ####
dft$subset = year(dft$date) - (min(year(dft$date)) - 1)



#Verify time blocking
verify = dft %>% group_by(subset) %>% summarise(start_time=min(date),
                                                end_time=max(date),
                                                n_obs = n())
lta %>% filter(x==799653, y==833076)

create_lta = function(i,input){
  #Calculate the null model now based on the train informaton we've "seen" before the prediction takes place


  lta = input %>% filter(subset<=i-1) %>% group_by(Week,x,y) %>% summarise(chl_LTA = mean(chl,na.rm=T))
  lta_xy = input %>% filter(subset<=i-1) %>% group_by(x,y) %>% summarise(chl_LTA_xy = mean(chl,na.rm=T))
  
  #If pixel-week data not available default to pixel level
  lta = left_join(lta_xy,lta,by=c("x","y"))
  lta$subset = i
  
  return(lta)
}

lta[is.na(lta$chl_LTA),]


lta=lapply(c(3:8),create_lta,input=dft) %>% bind_rows()
dft = left_join(dft,lta,by=c("Week","subset","x","y"))
dft = dft %>% mutate(chl_LTA = ifelse(is.na(chl_LTA),chl_LTA_xy,chl_LTA)) %>% dplyr::select(-chl_LTA_xy)

dft %>% select(Week,x,y,chl_LTA, subset) %>% filter(subset>2)


#### try exporting all as one matrix ####
##### Create and save the data in train test files for the convLSTM format #####
preds = c("chl","PRISMARDSW","precip","surfaceRadiation","windSpeed","dMean")
n_features = length(preds)



print(paste0("Number of timesteps: ",dft$time %>% unique() %>% length))
num_time_steps = dft$time %>% unique() %>% length
time_unique = sort(unique(dft$time))
print(paste0("Number of x: ",dft %>% select(x) %>% unique() %>% nrow))
x_dim = dft %>% select(x) %>% unique() %>% nrow
x_unique = sort(unique(dft$x))
print(paste0("Number of y: ",dft %>% select(y) %>% unique() %>% nrow))
y_dim = dft %>% select(y) %>% unique() %>% nrow
y_unique = sort(unique(dft$y))

#### save predictor data ####

dft = dft %>% arrange(x,y,Time)
head(dft)
Xlist = list()
for (pred in seq_along(preds)) {
  m = acast(dft, Time~x~y, value.var=c(preds[pred]))
  Xlist[[pred]] = m
}


x_train_matrix = abind(Xlist,along = 4) #In form  Time, X, Y, Predictors, Subset

dim(x_train_matrix)

savename = paste0("/work/HAB4CAST/max_beal/CyANPixelForecast/data/convLSTM_train_test/convlstm_train_x_all.rds")
saveRDS(x_train_matrix,file=savename)
print(dim(x_train_matrix))



#### Save chl data ####


y_train_matrix = acast(dft, Time~x~y, value.var="chl_lead")
y_train_matrix=array(y_train_matrix,dim=dim(y_train_matrix))
savename = paste0("/work/HAB4CAST/max_beal/CyANPixelForecast/data/convLSTM_train_test/convlstm_train_y_all.rds")
saveRDS(y_train_matrix,file=savename)
print(dim(y_train_matrix))

rm(y_train_matrix)

dft = dft %>% arrange(x,y,Time)

#### Save x and y data ####
x = acast(dft, Time~x~y, value.var="x")
x=array(x,dim=dim(x))
savename = paste0("/work/HAB4CAST/max_beal/CyANPixelForecast/data/convLSTM_train_test/convlstm_x_dim.rds")
saveRDS(x,file=savename)
rm(x)

dft = dft %>% arrange(x,y,Time)

y = acast(dft, Time~x~y, value.var="y")
y=array(y,dim=dim(y))
savename = paste0("/work/HAB4CAST/max_beal/CyANPixelForecast/data/convLSTM_train_test/convlstm_y_dim.rds")
saveRDS(y,file=savename)


dft = dft %>% arrange(x,y,Time)


lta = acast(dft, Time~x~y, value.var="chl_LTA")
lta=array(lta,dim=dim(lta))
savename = paste0("/work/HAB4CAST/max_beal/CyANPixelForecast/data/convLSTM_train_test/convlstm_lta_dim.rds")
saveRDS(lta,file=savename)
dim(lta)


