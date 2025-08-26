#This code is written to develop CRPSS and RPSS verification metrics for evaluation of the pixel forecast in FL
#RPSS and CRPSS measure the performance of probabilistic forecasts compared to a reference forecast (climatology)
#Code should be run in parallel using RPSS_CRPSS.sbatch in /work/HAB4CAST/max_beal/CyANPixelForecast/code_testing/
#Author: Max Beal

library(scoringRules)
library(verification)
library(dplyr)
library(data.table)
library(stringr)
library(future)
library(future.apply)
library(parallel)
library(readr)
library(lubridate)
library(tidyr)

#### Settings ####
#Choose metric, only one can run at a time (RPSS: Ranked Probability Skill Score, CRPSS: Continuous Ranked Probability Skill Score)
RPSS=F
CRPSS=T

#Choose model output, only one at a time
RF=FALSE
CLSTM=TRUE
BRMS=FALSE

#Choose lead time (1, 2, or 4). Note: this will only work for random forest (RF=TRUE, above)
lead = 1 


#Output format: when lead=1 outputs will have the form Model_Score.csv ("RF_RPSS.csv")
#When lead > 1, outputs have form ModelLead_Score ("RFL4_CRPSS.csv")

#Load results
if (RF) {
  if (lead==1) {
    folds = read_csv("/work/HAB4CAST/max_beal/CyANPixelForecast/output/RF_folds_all.csv")
  }
  
  if (lead==2) {
    folds = read_csv("/work/HAB4CAST/max_beal/CyANPixelForecast/output/RF_folds2lead.csv")
    print("lead 2 week")
  }
  
  if (lead==4) {
    folds = read_csv("/work/HAB4CAST/max_beal/CyANPixelForecast/output/RF_folds4lead.csv")
    print("lead 4 week")
  }
  
}

if (CLSTM) {
  folds = read_csv("/work/HAB4CAST/max_beal/CyANPixelForecast/output/CLSTM_folds_all.csv")
}

if (BRMS) {
  folds = read_csv("/work/HAB4CAST/max_beal/CyANPixelForecast/output/BRMS_folds_all.csv")
}


#### RPSS ####
if (RPSS) {
  #Get Ensemble Predictions
  ens_select = grepl("enschl",colnames(folds))
  ens_cols = folds[,ens_select,with=FALSE] #Data.table weirdness
  boxplot(t(log(ens_cols[1:52,])),outline=FALSE)
  abline(h=c(12,24), col=c("orange", "red"), lty=c(2,2), lwd=c(3, 3))
  
  
  
  #Cut into categories
  breaks = c(-Inf,3,12,24,Inf)
  labels = c("Non Detect","Vigilance","Alert Level 1","Alert Level 2")
  ens_mat = as.matrix(ens_cols,byrow=T)
  ens_cat = cut(ens_mat,breaks,labels)
  
  head(ens_cat)
  
  dim(ens_cat) = dim(ens_mat)
  
  #Bind together
  tbl = apply(ens_cat,1,table)
  tbl=do.call(rbind,tbl)
  tbl=as.data.frame(tbl)
  colnames(tbl) = labels
  
  
  #convert to probs
  forecast_probs = tbl/100
  
  
  #Get observations
  obs = folds$chl_observed
  obs = cut(obs,breaks,labels)
  obs = as.numeric(obs)
  
  ind = folds %>% dplyr::select(x,y,Time, year)
  obs=cbind(ind,obs)
  
  forecast_probs = cbind(ind,forecast_probs)
  
  
  xy = folds %>% dplyr::select(x,y) %>% unique()
  
  #### Create observed probability (reference forecast) ####
  #This only needs to be done once
  get_observed_probs = FALSE
  if (get_observed_probs) {
    
    
    df <- fread("/work/HAB4CAST/max_beal/CyANPixelForecast/data/florida_data_spatial.csv")
    df$subset = year(df$date) - (min(year(df$date)) - 1)
    
    get_ob_probs = function(index, xy=xy, df=df, i){
      pixel = df %>% filter(x==xy$x[index],y==xy$y[index], subset<=i-1)
      pr_ob = cut(pixel$chl,breaks,labels) %>% table() / nrow(pixel)
      return(pr_ob)
    }
    
    
    plan(multisession, workers = (availableCores()))
    tbllist=list()
    for (i in c(3:8)) {
      
      print(availableCores())
  
      prob = future_lapply(c(1:nrow(xy)),xy=xy,df=df,i=i,get_ob_probs)
      
      prob = lapply(c(1:nrow(xy)),xy=xy,df=df,i=i,get_ob_probs)
      tbl=do.call(rbind,prob)
      tbl=cbind(xy,tbl)
      
      
      tbl$subset = i
      tbl$year_to_predict = (2017+(i-1))
  
      
      tbllist[[i-2]]=tbl
      
    }
    plan(sequential)
  
    tbl_out = tbllist %>% bind_rows
    
    write_csv(tbl_out,"/work/HAB4CAST/max_beal/CyANPixelForecast/data/chlorophyll_pixel_climatology_prob.csv")
  
  }
  
  
  if (!get_observed_probs) {
    observed_probs = read_csv("/work/HAB4CAST/max_beal/CyANPixelForecast/data/chlorophyll_pixel_climatology_prob.csv")
  }
  
  
  #Now subset and calculate for each time period, need to calculate for each fold because null model gets updated
  #Will still calculate one rpss for each timestep but the null will be updated
  
  index=1
  calculate_rpss = function(index, xy=xy, obs=obs,forecast_probs=forecast_probs,observed_probs=observed_probs){
    observation = obs %>% dplyr::filter(x==xy$x[index],y==xy$y[index])
    pixel_o = observed_probs %>% filter(x==xy$x[index],y==xy$y[index])
    pixel_f = forecast_probs %>% filter(x==xy$x[index],y==xy$y[index])
    
    ind =observation %>% dplyr::select(x,y,Time,year)
    ob = observation %>% dplyr::select(-x,-y,-Time,-year)
    fctpr = pixel_f %>% dplyr::select(-x,-y,-Time,-year)
    obspr = pixel_o %>% dplyr::select(-x,-y)
    
    
    rpss=c()
    for (i in 1:nrow(ob)) {
      
      
      f = matrix(as.numeric(fctpr[i,]),nrow = 1)
      o = matrix(as.numeric(obspr),nrow=1)
      
      
      rpss = c(rpss,rps(obs = ob$obs[i], pred = f, baseline=o)$rpss)
      
    }
    
    rpss=cbind(rpss,ind) #cbind index
    return(rpss)
  }
  
  #out = calculate_rpss(index=1,xy=xy, obs=obs,forecast_probs=forecast_probs,observed_probs=observed_probs)
  
  #Write for loop to do this for several years
  
  rpss_list = list()
  for (yr in unique(observed_probs$year_to_predict)) {
    
    #Filter observations
    obs_sub = obs %>% filter(year==yr)
    
    #Filter the observed probs to an individual year to get correct null
    observed_probs_sub = observed_probs %>% filter(year_to_predict==yr)
    
    #Filter forecast probs to get correct year of predictions to compare
    forecast_probs_sub = forecast_probs %>% dplyr::filter(year==yr)
    
    #Fitler xy to observed cells
    xy = obs_sub %>% dplyr::select(x,y) %>% unique()
    
    print(yr)
    
    plan(multisession, workers = (availableCores()))
    print(availableCores())
    out = future_lapply(c(1:nrow(xy)),calculate_rpss,xy=xy, obs=obs_sub,forecast_probs=forecast_probs_sub,observed_probs=observed_probs_sub)
    plan(sequential)
    
    print(paste0("RPSS calculation finished for ",yr))
    
    #concatenate rpss values
    #rpss = do.call(c,out)
    rpss = out %>% bind_rows()
    #rpss = cbind(ind,rpss)
    padded_numbers <- sprintf("%06d", rpss$Time)
    rpss$week = substr(padded_numbers,1,2)
    rpss$year = substr(padded_numbers,3,6)
    rpss=rpss %>%
      mutate(beginning = ymd(str_c(year, "-01-01")), date = beginning + weeks(week)) %>% select(-beginning)
  
    rpss$year = yr
    
    #Put values into list
    rpss_list[[yr-min(yr)+1]] = rpss
  }
  
  
  rpss_out = rpss_list %>% bind_rows
  
  
  
  
  
  if (RF) {
    if (lead==1) {
      write_csv(rpss_out,"/work/HAB4CAST/max_beal/CyANPixelForecast/output/RF_RPSS.csv")
    }
    if (lead==2) {
      write_csv(rpss_out,"/work/HAB4CAST/max_beal/CyANPixelForecast/output/RFL2_RPSS.csv")
    }
    if (lead==4) {
      write_csv(rpss_out,"/work/HAB4CAST/max_beal/CyANPixelForecast/output/RFL4_RPSS.csv")
    }
  }
  
  if (BRMS) {
    write_csv(rpss_out,"/work/HAB4CAST/max_beal/CyANPixelForecast/output/BRMS_RPSS.csv")
  }
  
  if (CLSTM) {
    write_csv(rpss_out,"/work/HAB4CAST/max_beal/CyANPixelForecast/output/CLSTM_RPSS.csv")
  }


}



#### CRPSS ####
if (CRPSS) {
  
  if (BRMS) {
    folds = folds %>% rename(chl_lead=logchl_lead,chl=logchl)
  }
  
library(scoringRules)

folds$subset = (folds$year - (min(folds$year) - 3)) #Recreate fold values, need -3 here b/c starts in 2019


df <- fread("/work/HAB4CAST/max_beal/CyANPixelForecast/data/florida_data_spatial.csv") #Read in df to construct observed forecast distribution
df$subset = year(df$date) - (min(year(df$date)) - 1)
df$Time = paste0(df$Year,str_pad(df$Week,2,pad=0))

xy = folds %>% dplyr::select(x,y) %>% unique()
crpss_list = list()

calculate_crpss = function(fold,df,folds,xy){
      
  for (j in c(1:nrow(xy))) {
      #print(head(df))
      pixel_o = df %>% filter(x==xy$x[j],y==xy$y[j],subset<fold)
      pixel = folds %>% filter(x==xy$x[j],y==xy$y[j],subset==fold)
      
      if (nrow(pixel_o)==0) {next}
      
      
      pixel = pixel %>% arrange(x,y,Time)
      pixel_o$year = year(pixel_o$date)
      pixel_o = pixel_o %>% arrange(x,y,Time) %>% select(Week,year,chl)
      
      
      
      pixel_forecast = pixel %>% arrange(Week) %>% select(contains("enschl"))
      pixel_observation = pixel %>% arrange(Week) %>% select(chl_lead)
      W = pixel %>% arrange(Week) %>% select(Week)
      Time = pixel %>% arrange(Week) %>% select(Time)
      pixel_null_forecast =pixel_o %>% pivot_wider(id_cols=Week,names_from=year,values_from=chl,values_fn = mean) %>% arrange(Week) #values_fn=mean takes the mean if there are somehow two values per year-week (more of an issue with the fake data)
      pixel_null_forecast = suppressMessages(left_join(W,pixel_null_forecast))
      
    
      #forecast ensemble
      forecast = pixel_forecast %>% mutate(across(everything(), as.numeric)) %>%
        as.matrix()
      #forecast = as.matrix(pixel_forecast)
      
      #Null forecast, fill na with annual average
      null_forecast = as.matrix(pixel_null_forecast)
      null_forecast = pixel_null_forecast
      
      fill_na_with_mean <- function(df) {
        overall_mean <- mean(unlist(df), na.rm = TRUE)  # Mean of all numeric values
        
        df[] <- lapply(df, function(col) {
          if (all(is.na(col))) {
            return(if (is.numeric(col)) rep(overall_mean, length(col)) else col)
          } else if (is.numeric(col)) {
            col[is.na(col)] <- mean(col, na.rm = TRUE)
            return(col)
          } else {
            return(col)  # Leave non-numeric columns unchanged
          }
        })
        
        return(df)
      }
      
      
      #Fill NA values, first use the annual mean if there's no value for the given week
      #Second, use all time value if there's no data to average over the year
      null_forecast=null_forecast %>% select(-Week)
      null_forecast = fill_na_with_mean(null_forecast)
      
      
      if (nrow(na.omit(null_forecast))==0) {next}

      
      #Observation
      observation = pixel_observation$chl_lead
      
      #CRPSS = 1 - CRPS_exp / CRPS_ref
      
      null_forecast=null_forecast %>% as.matrix()
      # Calculate CRPS
      crps_value <- crps_sample(y = observation, dat = forecast)
      crps_ref <- crps_sample(y = observation, dat = null_forecast) #null model here will just be the observed values of chla before the year of interest
      
      crpss = 1 - crps_value/crps_ref
      
      
      if (CLSTM) {
        crpss_df=data.frame("x"=xy$x[j],"y"=xy$y[j],"crpss"=crpss,"Time"=Time,"subset"=fold)
      }
      
      if (!CLSTM) {
        crpss_df=data.frame("COMID"=pixel$COMID,"x"=xy$x[j],"y"=xy$y[j],"crpss"=crpss,"Time"=Time,"subset"=fold)
      }
      
      crpss_list[[length(crpss_list) + 1]] = crpss_df
      }

  return(crpss_list)
}

#out = calculate_crpss(1,df=df,folds = folds,xy=xy) %>% bind_rows() #For testing a single output

plan(multisession, workers = (availableCores()))
print(availableCores())
crpss_out = future_lapply(c(3:8),calculate_crpss,df=df,folds=folds,xy=xy) %>% bind_rows()
plan(sequential)

hist(crpss_out$crpss)
mean(crpss_out$crpss>0)

if (RF) {
  if (lead==1) {
    write_csv(rpss_out,"/work/HAB4CAST/max_beal/CyANPixelForecast/output/RF_CRPSS.csv")
  }
  if (lead==2) {
    write_csv(rpss_out,"/work/HAB4CAST/max_beal/CyANPixelForecast/output/RFL2_CRPSS.csv")
  }
  if (lead==4) {
    write_csv(rpss_out,"/work/HAB4CAST/max_beal/CyANPixelForecast/output/RFL4_CRPSS.csv")
  }
}

if (BRMS) {
  write_csv(crpss_out,"/work/HAB4CAST/max_beal/CyANPixelForecast/output/BRMS_CRPSS.csv")
}

if (CLSTM) {
  write_csv(crpss_out,"/work/HAB4CAST/max_beal/CyANPixelForecast/output/CLSTM_CRPSS.csv")
}

}
