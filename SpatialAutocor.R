#This notebook looks at pairwise comparisons of chl correlations at the pixel level
#Spatial autocorrelation should be explicity modeled if we're not using a structure that does so
#The idea here is to look into the average region of influence and include the spatial feature as a predictor
#Author: Max Beal


#### Libraries ####
library(readr)
library(INLA)
library(dplyr)
library(ggplot2)
library(sf)
library(tidyverse)
library(INLA) 
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
library(spdep)
library(data.table)


#### Read in the data frame ####
df <- fread("/work/HAB4CAST/max_beal/CyANPixelForecast/data/florida_data.csv")

df = df %>% rename(Cyano = value, Year = year, Bloom=bloom)

cicyano = 10^(df$Cyano*0.011714-4.1870866)
#df$chl = 6620*cicyano - 3.07 
df$chl = 6620*cicyano  # Try changing the intercept to make sure nothing goes below 0
df$chl = df$chl - min(df$chl,na.rm=T)
df = df %>% filter(chl>=0)
df$logchl = log(df$chl) #Make it normal

print(min(df$chl))

rm(cicyano)

comids = df$COMID %>% unique()

#comids = comids[comids!=166757656] #Remove okechobee for now
#Do only Okechobee
#comids = comids[comids==166757656]



cordistance = data.frame()

#### Spatial Correlation Analysis ####
for (comid in comids) {
  
  #Subset data
  lake = df %>% filter(COMID==comid) %>% arrange(date) %>% select(date,x,y,chl)
  wide_data <- lake %>% pivot_wider(names_from=c(x,y), values_from = chl)
  wide_data = wide_data %>% select(-date)
  
  #Create matrix of all grid combination correlations and get coordiantes
  cor_matrix = cor(wide_data,use="pairwise.complete.obs")
  coords = lake %>% select(x,y) %>% unique()
  rownames(coords) = paste(coords$x,coords$y,sep="_")
  
  
  #Merge coordiantes for grid 1
  cor_df = as.data.frame(as.table(cor_matrix)) %>% rename(grid1=Var1, grid2=Var2, correlation=Freq)
  
  cor_df = cor_df %>% left_join(coords %>% mutate(grid1=rownames(.)),by="grid1") %>% 
    rename(x1=x,y1=y)
  
  #Merge coordiantes for grid 1
  cor_df = cor_df %>% left_join(coords %>% mutate(grid2=rownames(.)),by="grid2") %>% 
    rename(x2=x,y2=y)
  
  #Resulting dataframe has head to head match ups of grids and their correlations
  grids = cor_df$grid1 %>% unique()
  
  #subset target grid and plot
  #grid_cor = cor_df %>% filter(grid1==grids[1])
  
  
  
  #Plot cell and target variables
  #ggplot() + geom_tile(data=grid_cor, aes(x=x2,y=y2,fill=correlation)) + geom_point(data=grid_cor,aes(x=x1,y=y1)) + scale_fill_viridis_c()
  print(length(grids))
  if (length(grids)>4) {
    
  
  #Use Pythagorean theorem to get distance
  distance_across_grid = sqrt((300^2) + (300^2))
  cor_df = cor_df %>% mutate(distance_m = sqrt((x2-x1)^2+(y2-y1)^2),
                             distance_bin=cut(
                               distance_m,
                               breaks=seq(0,max(distance_m),by=distance_across_grid),
                               labels=FALSE))
  
  
  cor_by_dist = cor_df %>% group_by(distance_bin) %>% summarise(ave_correlation = mean(correlation,na.rm=T))
  
  
  cor_by_dist$COMID = comid
  
  cordistance=rbind(cordistance,cor_by_dist)
  }
  
  #Plots
  # barplot(names.arg=cor_by_dist$distance_bin,cor_by_dist$ave_correlation,axis.lty=1)
  # grid_cor = cor_df %>% filter(grid1==grids[2])
  # ggplot() + geom_tile(data=grid_cor, aes(x=x2,y=y2,fill=distance_bin)) + geom_point(data=grid_cor,aes(x=x1,y=y1)) + scale_fill_viridis_c()

}


#write.csv(cordistance,"/work/HAB4CAST/max_beal/CyANPixelForecast/output/spatial_autocor.csv")
# cordistance$COMID %>% unique() %>% length()
# 
a=ggplot(data=cordistance) + geom_boxplot(aes(x=distance_bin,y=ave_correlation,group=distance_bin)) + 
  theme_minimal() + geom_abline(slope=0,intercept=0,linetype="dashed",color="red") + ylab("Pearson Correlation") + xlab("Distance (pixels)")
ggsave("SpatialAutoCor.png",plot=a,width=12,height=6,dpi=300)
# 
# 
# ggplot(data=cordistance,aes(x=distance_bin,y=ave_correlation)) + geom_point() + theme_minimal() + 
#   geom_smooth() + geom_abline(slope=0,intercept=0,linetype="dashed",color="red")
# 
# 
# avecor = cordistance %>% group_by(distance_bin) %>% summarise(avecor=mean(ave_correlation,na.rm=T))


# plot(avecor,type="b")
# avecor$avecor %>% min(na.rm=T)

#### Generate 4x4 Spatial Autocor ####
albers = "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
df_sf=st_as_sf(df, coords=c("x","y"),crs=albers)

for (comid in comids) {
  
  lake = df_sf %>% filter(COMID==comid) %>% arrange(date)
  
  
  #df_sf = df_sf %>% filter(COMID!=166757656)
  
  distance_across_grid = sqrt((300^2) + (300^2))
  neighborchl = lake %>% mutate(geometry_buffer=st_buffer(geometry,dist=distance_across_grid*4)) %>%
    rowwise() %>% group_by(date) %>%
    mutate(mean_neighbor_chl = mean(lake$chl[st_intersects(geometry,geometry_buffer)[[1]]],na.rm=T)) %>% ungroup()
  
  neighborchl = neighborchl%>%
    dplyr::mutate(x = sf::st_coordinates(.)[,1],
                  y = sf::st_coordinates(.)[,2]) %>% st_drop_geometry() %>% select(-geometry_buffer)
  
  fp = "/work/HAB4CAST/max_beal/CyANPixelForecast/data/neighbor_chl/"
  fn = paste0("neighbor_chlorophyll_",unique(neighborchl$COMID),".csv")
  save = paste0(fp,fn)
  write_csv(neighborchl,save)

}



#### Generate Lake wide log chlorophyll-a ####
lakechl = df %>% group_by(COMID,date) %>% summarise(whole_lake_chl = mean(chl))
lakechl

spac = list.files("/work/HAB4CAST/max_beal/CyANPixelForecast/data/neighbor_chl/",full.names=TRUE)

sp4=read_csv(spac) %>% bind_rows()

sp4 = merge(sp4,lakechl,by=c("COMID","date"))

sp4 = sp4 %>% select(-geometry, -geometry_buffer)

write_csv(sp4, "/work/HAB4CAST/max_beal/CyANPixelForecast/data/florida_data_spatial.csv")


sp4 %>% select(x,y,date,chl)
min(sp4$chl)


