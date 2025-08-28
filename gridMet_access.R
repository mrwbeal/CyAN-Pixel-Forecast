library(terra)
library(raster)
library(ncdf4)
library(sf)
library(tidync)
library(exactextractr)
library(dplyr)
library(tidyr)
library(reshape2)
library(lubridate)
library(data.table)

albers = "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"



#### Choose dataset of interest and years ####
dataset = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/MET/srad/srad_' #Surface Radiation
varname = 'surface_downwelling_shortwave_flux_in_air'
  
# dataset = 'http://thredds.northwestknowledge.net:8080/thredds/dodsC/MET/vs/vs_' #Wind Speed
# varname = 'wind_speed'


start_year = 2016
end_year = 2024

years = c(start_year:end_year)

get_gridMEt = function(dataset,varname, year){
  
  url = paste0(dataset,year,".nc")
  nd = nc_open(url)
  
  y = ncvar_get(nd,"lat")
  x = ncvar_get(nd,"lon")
  t = ncvar_get(nd,"day") #days since 1900-01-01
  
  
  #Load in necessary SF data to crop Florida
  states = st_read("/work/HAB4CAST/cyano_forecast/data/cb_2019_us_state_500k/cb_2019_us_state_500k.shp")
  not_conus <- c("VI","HI","AK","MP","PR","GU","AS")
  states = states %>% dplyr::filter(!STUSPS %in% not_conus)
  states=st_transform(states,'epsg:4326')
  fl = states %>% dplyr::filter(STUSPS=="FL")
  
  #Get florida extents
  xmin = ext(fl)[1]
  xmax= ext(fl)[2]
  ymin = ext(fl)[3]
  ymax = ext(fl)[4]
  
  #Define cropping region
  lat_idx = which(y >= ymin & y <= ymax)
  lon_idx = which(x >= xmin & x <= xmax)
  
  min_lon = min(lon_idx)
  max_lon = max(lon_idx)
  
  min_lat = min(lat_idx)
  max_lat = max(lat_idx)
  
  #Get annual dataset and crop to florida
  sr_data = ncvar_get(nd, varname,
            start=c(min_lon,min_lat,1),
            count = c(length(lon_idx),length(lat_idx),-1))
  
  #retrieve indices
  sr_lat = ncvar_get(nd, 'lat')[lat_idx]
  sr_lon = ncvar_get(nd, 'lon')[lon_idx]
  sr_day = ncvar_get(nd, 'day')
  sr_date <- as.Date(t, origin="1900-01-01")
  
  #Reshape data
  sr_reshape = aperm(sr_data,c(2,1,3))
  sr_reshape = melt(sr_data)
  colnames(sr_reshape) = c("lon_idx","lat_idx","date","value")
  
  #Create Rast Object
  sr_rast = rast(sr_reshape,type="xylz")
  ext(sr_rast) <- c(min(sr_lon),max(sr_lon),min(sr_lat),max(sr_lat))
  crs(sr_rast) = 'epsg:4326'
  sr_rast = flip(sr_rast,direction="vertical")
  
  names(sr_rast) = sr_date
  
  #Transform SF
  fl = st_transform(fl,albers)
  
  fl_cells= st_read("/work/HAB4CAST/max_beal/CyANPixelForecast/data/fl_cells_shapefile.shp") %>% st_transform(albers)
  
  fl_cells = fl_cells %>%
    dplyr::mutate(x = sf::st_coordinates(.)[,1],
                  y = sf::st_coordinates(.)[,2])
  
  test = fl_cells %>% filter(COMID==10361596)
  plot(test$geometry)
  
  #Extract data and aggregate to weekly level
  sr_rast_albers = project(sr_rast,albers)
  
  data_extract = terra::extract(sr_rast_albers,fl_cells,ID=FALSE)

  data_extract = cbind(fl_cells %>% dplyr::select(COMID,x,y),data_extract) %>% st_drop_geometry()

  
  
  data = melt(data_extract,id=c("COMID","x","y"))
  data$variable = as.Date(sub("^X","",data$variable),format="%Y.%m.%d")
  data = data %>% mutate(time = paste0(week(variable),year(variable)))
  
  #Final Dataset
  sr_final = data %>% group_by(COMID,x,y,time) %>% summarise('surfaceRadiation'=mean(value,na.rm=T)) #Change value
  
  return(sr_final)
  
}



data = lapply(years,get_gridMEt,dataset=dataset,varname=varname) %>% bind_rows()

fwrite(data,"/work/HAB4CAST/max_beal/CyANPixelForecast/data/surface_radiation_fl.csv")


data





library(ggplot2)

#Load in necessary SF data to crop Florida
lakes <- st_read("/work/HAB4CAST/max_beal/SW_model/data/OLCI_resolvable_lakes_2022_09_08/OLCI_resolvable_lakes_2022_09_08.shp")
states = st_read("/work/HAB4CAST/cyano_forecast/data/cb_2019_us_state_500k/cb_2019_us_state_500k.shp")
not_conus <- c("VI","HI","AK","MP","PR","GU","AS")
states = states %>% dplyr::filter(!STUSPS %in% not_conus)
states=st_transform(states,'epsg:4326')
lakes=st_transform(lakes,'epsg:4326')
fl = states %>% dplyr::filter(STUSPS=="FL")
sf_use_s2(FALSE)
fl_lakes = st_crop(lakes,fl)


ggplot() + geom_tile(data=data,aes(x=x,y=y,fill=windSpeed)) + geom_sf(data=fl_lakes,fill=NA)



srclim=data %>% group_by(time) %>% summarise(msr = mean(windSpeed,na.rm=T),
                                                 q25=quantile(windSpeed,0.25,na.rm=T),
                                                 q75=quantile(windSpeed,0.75,na.rm=T))

srclim=melt(srclim)

ggplot(srclim) + geom_line(aes(x=time,y=value,color=variable,group=variable))


