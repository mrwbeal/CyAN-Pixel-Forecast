# This code is for joining the predictor datasets for Florida pixel level modeling
# This code is the updated version of the cyan_pixel_data_join 
# Author: Max Beal
# Date created: 12/2/2024

library(readr)
library(feather)
library(data.table)
library(dplyr)
library(parallel)
library(future)
library(future.apply)
library(lubridate)
library(stringr)
library(sf)
library(tigris)
library(data.table)
library(RSQLite)
library(DBI)
library(ggplot2)






#Cyanobacteria files
cyan_files = list.files("/work/HAB4CAST/max_beal/CyANPixelForecast/data/cyan_pixel_conus_lakes_current",full.names = TRUE)

dft = read_feather(cyan_files[1])
dft$subset = year(dft$date) - (min(year(dft$date)) - 1)
#Verify time blocking
verify = dft %>% group_by(subset) %>% summarise(start_time=min(date),
                                                end_time=max(date),
                                                n_obs = n())


#Watershed scale precipitation
p = read_csv("/work/HAB4CAST/max_beal/CyANPixelForecast/data/precip_watershed_data.csv")
p=p %>% dplyr::select(COMID,sum,Time) %>% rename(precip= sum, time=Time)

#Wind speed (m/s)
ws = read_csv("/work/HAB4CAST/max_beal/CyANPixelForecast/data/wind_speed_fl.csv")

#Surface radiance (w/m^2)
sr = read_csv("/work/HAB4CAST/max_beal/CyANPixelForecast/data/surface_radiation_fl.csv")


#Surface water temp files
swfn = list.files("/work/HAB4CAST/max_beal/CyANPixelForecast/data/swtemp_pixel_conus_combined_datasets",full.names = T)

# oke = read_csv(swfn[100]) %>% filter(COMID==166757656)
# 
# 
# ggplot(oke%>% na.omit()) + geom_tile(aes(x=x,y=y,fill=`202016PRISMARD`))
# ggplot(oke) + geom_tile(aes(x=x,y=y,fill=`202016PRISM`))




combine_all_sw = TRUE
if (combine_all_sw) {
  read_sw = function(file){
    df = read_csv(file)
    time = gsub("\\D","",substr(file,15,150))
    if (ncol(df)==6) {
      colnames(df)[3] = "PRISMSW"
      colnames(df)[2] = "PRISMARDSW"
    }
    if (ncol(df)==5) {
      colnames(df)[2] = "PRISMSW"
      df = df%>% mutate("PRISMARDSW"=NA)
    }
    
    df <- df[, c("COMID", "PRISMARDSW", "PRISMSW","x","y")]
    df$time=time
    return(df)
  }
  plan(multisession, workers = (availableCores()-1))
  print(availableCores())
  sw = future_lapply(swfn, read_sw) %>% bind_rows()
  write_csv(sw,"/work/HAB4CAST/max_beal/CyANPixelForecast/data/allsw.csv")
  plan(sequential)
}

# Write the full SW file to an SQL database
sw = fread("/work/HAB4CAST/max_beal/CyANPixelForecast/data/allsw.csv")
conn <- dbConnect(RSQLite::SQLite(), "/work/HAB4CAST/max_beal/CyANPixelForecast/data/prism_sw_database.db")
dbWriteTable(conn, "sw_temp_all_prism_ard", sw,overwrite=TRUE)

#query = paste0("SELECT * FROM PRISM_sw_temp_all WHERE week = ",1," AND year = ",2016)

COMIDs = as.numeric(str_extract(substr(cyan_files,15,150),"\\d+"))

#check already processed files
#processed = list.files("/work/HAB4CAST/max_beal/CyANPixelForecast/data/compiled_pixel_conus_lakes_current")
#COMID_processed = as.numeric(str_extract(processed,"\\d+"))

#Save COMIDs that do not yet have a file
#COMIDs = COMIDs[!COMIDs %in% COMID_processed]


#### FOR NOW JUST DO FLORIDA ####
process_fl=TRUE
if (process_fl) {
  
  #Read in lakes and subset florida
  conus_lakes <- st_cast(st_read("/work/HAB4CAST/max_beal/SW_model/data/OLCI_resolvable_lakes_2022_09_08/OLCI_resolvable_lakes_2022_09_08.shp"), "MULTIPOLYGON")
  states <- states(cb = TRUE)
  fl = states %>% filter(STUSPS =="FL")
  fl = st_transform(fl,st_crs(conus_lakes))
  fl_lakes = st_crop(conus_lakes,fl)
  #Concatenate florida data files
  fl_comid = as.numeric(fl_lakes$COMID)
  
 
  
  fl_files = cyan_files[COMIDs %in% fl_comid]
  
  COMIDs = COMIDs[COMIDs %in% fl_comid]
  
}

#Pull florida surface water and
value_string <- paste(shQuote(COMIDs), collapse = ", ")
conn <- dbConnect(RSQLite::SQLite(), "/work/HAB4CAST/max_beal/CyANPixelForecast/data/prism_sw_database.db")
query = paste("SELECT * FROM sw_temp_all_prism_ard WHERE COMID IN (", value_string, ")", sep = "")
sw= dbGetQuery(conn, query)


#ggplot(sub) + geom_point(aes(x=PRISMARDSW,y=PRISMSW))  + geom_abline(intercept = 0, slope = 1) + theme_classic()


#166757656 Okechobee
#Join all together

merge_datasets = function(i, fl_files, COMIDs, sw, sr, ws, p){

  #CyAN file
  flcyan = read_feather(fl_files[i])
  flcyan$time = paste0(flcyan$Week,flcyan$year)
  flcyan = flcyan %>% mutate(x=round(x),y=round(y))
  
  #Merge SW
  sub = sw %>% filter(COMID==COMIDs[i])
  sub = sub %>% mutate(x=round(x),y=round(y))
  flmerged = merge(flcyan,sub, on=c("COMID","x","y","time"))
  
  #Merge Surface Radiation
  srsub = sr %>% filter(COMID==COMIDs[i])
  srsub = srsub %>% mutate(x=round(x),y=round(y))
  flmerged = merge(flmerged,srsub, on=c("COMID","x","y","time"))
  
  #Merge Wind Speed
  wssub = ws %>% filter(COMID==COMIDs[i])
  wssub = wssub %>% mutate(x=round(x),y=round(y))
  flmerged = merge(flmerged,wssub, on=c("COMID","x","y","time"))
  
  #Merge Precip
  psub = p %>% filter(COMID==COMIDs[i])
  flmerged=merge(flmerged,psub,on=c("COMID","time"))
  
  flmerged = flmerged %>% mutate(bloom= as.numeric(bloom))
  
  #Return merged dataset
  return(flmerged)

}



# test = merge_datasets(1,fl_files,COMIDs,sw,sr,ws,p)
# 
# plot(test$value,test$windSpeed)
# 
# ggplot() + geom_point(data=srsub,aes(x=x,y=y))+ geom_point(data=swpts,aes(x=x,y=y),color="blue")


plan(multisession, workers = (availableCores()))
print(availableCores())
cb = lapply(seq_along(fl_files), merge_datasets,fl_files=fl_files,COMIDs=COMIDs,sw=sw,sr=sr,ws=ws,p=p) %>% bind_rows()
write_csv(cb,"/work/HAB4CAST/max_beal/CyANPixelForecast/data/florida_data.csv")
plan(sequential)

#for some reason, read_csv thinks PRISMARDSW is a lgl
df = fread("/work/HAB4CAST/max_beal/CyANPixelForecast/data/florida_data.csv")


df %>% select(bloom,PRISMSW) %>% na.omit() %>% summarise(bis_cor = cor(bloom,PRISMSW,method="pearson"))
df %>% select(bloom,precip) %>% na.omit() %>% summarise(bis_cor = cor(bloom,precip,method="pearson"))
df %>% select(bloom,PRISMARDSW) %>% na.omit() %>% summarise(bis_cor = cor(bloom,PRISMARDSW,method="pearson"))
df %>% select(bloom,surfaceRadiation) %>% na.omit() %>% summarise(bis_cor = cor(bloom,surfaceRadiation,method="pearson"))
df %>% select(bloom,windSpeed) %>% na.omit() %>% summarise(bis_cor = cor(bloom,windSpeed,method="pearson"))

# 
# ggplot(df %>% na.omit()) + geom_tile(aes(x=x,y=y,fill=PRISMARDSW))
# ggplot(df %>% na.omit()) + geom_tile(aes(x=x,y=y,fill=PRISMSW))


df %>% select(value,PRISMSW) %>% na.omit() %>% summarise(bis_cor = cor(value,PRISMSW,method="pearson"))
df %>% select(value,precip) %>% na.omit() %>% summarise(bis_cor = cor(value,precip,method="pearson"))
df %>% select(value,PRISMARDSW) %>% na.omit() %>% summarise(bis_cor = cor(value,PRISMARDSW,method="pearson"))
df %>% select(value,surfaceRadiation) %>% na.omit() %>% summarise(bis_cor = cor(value,surfaceRadiation,method="pearson"))
df %>% select(value,windSpeed) %>% na.omit() %>% summarise(bis_cor = cor(value,windSpeed,method="pearson"))


ggplot(df, aes(x=PRISMARDSW, y=PRISMSW) ) +
  geom_bin2d() + scale_fill_viridis_c() +
  theme_bw()


test= df %>% filter(!is.na(PRISMARDSW))

ggplot(test, aes(x=x,y=y)) + geom_tile()
