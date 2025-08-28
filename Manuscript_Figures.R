# This file is for creating and outputting the manuscript quality figures for the Florida CyANPixelForecast
#Author: Max Beal
#Created 12/13/24

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
library(ggspatial)
library(rnaturalearth)
library(cowplot)
library(raster)

#Read in lakes and subset florida
conus_lakes <- st_cast(st_read("/work/HAB4CAST/max_beal/SW_model/data/OLCI_resolvable_lakes_2022_09_08/OLCI_resolvable_lakes_2022_09_08.shp"), "MULTIPOLYGON")
states <- states(cb = TRUE)
fl = states %>% filter(STUSPS =="FL")
fl = st_transform(fl,st_crs(conus_lakes))
fl_lakes = st_crop(conus_lakes,fl)
#Concatenate florida data files
fl_comid = as.numeric(fl_lakes$COMID)

states=st_transform(states,st_crs(conus_lakes))
not_conus <- c("VI","HI","AK","MP","PR","GU","AS")

states = states %>% filter(!STUSPS%in%not_conus)

cells= rast("/work/HAB4CAST/max_beal/CyANPixelForecast/data/cellsRaster.tif")

flcells = crop(cells,fl_lakes)

flcelldf = as.data.frame(flcells,xy=TRUE)

flcelldf %>% group_by(COMID) %>% summarise(n()) %>% filter(COMID==16630838)

#### Figure 1. Map of Study Locations ####
main_map<-ggplot() + geom_sf(data=states)+ geom_sf(data=fl,fill="lightblue",color="black") + geom_sf(data=fl_lakes,fill="white",color="black") +
  annotation_scale(location="bl",width_hint=0.4) +
  annotation_north_arrow(location="tr",which_north="true",style=north_arrow_fancy_orienteering()) + theme_bw()+
  coord_sf(xlim = c(796752, 1600602), ylim = c(269573.6 , 961154.4))

ggplot() + geom_sf(data=fl_lakes,aes(fill=COMID==	16630838))
fl_lakes %>% arrange(desc(LakeArea_R))

e = fl_lakes %>% filter(COMID==16630838) %>% extent()

cyan = list.files("/work/HAB4CAST/max_beal/CyANPixelForecast/data/cyan_pixel_conus_lakes_current")

df = read_feather("/work/HAB4CAST/max_beal/CyANPixelForecast/data/cyan_pixel_conus_lakes_current/COMID_16630838.feather")

df = df %>% filter(date=="2017-01-01")

inset_map <- ggplot() +
  geom_sf(data=fl_lakes %>% filter(COMID==16630838)) +
  geom_tile(data=df, aes(x=x,y=y),fill="lightgreen",color="darkgreen")  + theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_sf(xlim = c(e[1] , e[2]), ylim = c(e[3], e[4])) + xlab("") + ylab("")

#Use theme_void if you want just the lake

combined_map <- ggdraw() +
  draw_plot(main_map) +
  draw_plot(inset_map, x = 0.3, y = 0.3, width = 0.25, height = 0.25)

ggsave("florida_map.png",plot=combined_map,width=8,height=6,dpi=300)


#### Figure 2. Spatial and Temporal Autocor ####
cordistance = read.csv("/work/HAB4CAST/max_beal/CyANPixelForecast/output/spatial_autocor.csv")
acm = read.csv("/work/HAB4CAST/max_beal/CyANPixelForecast/output/temporal_autocor.csv")

acm$variable = as.numeric(substr(acm$variable,3,4))

a=ggplot(data=cordistance %>% filter(distance_bin<50)) + geom_boxplot(aes(x=distance_bin,y=ave_correlation,group=distance_bin)) + 
  theme_minimal() + geom_abline(slope=0,intercept=0,linetype="dashed",color="red") + ylab("Pearson Correlation") + xlab("Distance (pixels)")

b=ggplot(data=acm) + geom_boxplot(aes(x=variable,y=value,group=variable)) + 
  theme_minimal() + geom_abline(slope=0,intercept=0,linetype="dashed",color="red") + ylab("Pearson Correlation") + xlab("Lag (weeks)")

fig2 = ggarrange(a,b,ncol=1)

ggsave("SpatialTemporalAutoCor_abbreviated.png",plot=fig2,width=6,height=6,dpi=300)


#### Figure 5. Map of Forecast Accuracies #### 

map = read_csv("/work/HAB4CAST/max_beal/CyANPixelForecast/output/map_accuracies_rf.csv")

main_map<-ggplot()+geom_sf(data=states)+ geom_sf(data=fl,fill="lightblue",color="black") + geom_sf(data=fl_lakes,fill="white",color="black") +
  annotation_scale(location="bl",width_hint=0.4) +
  annotation_north_arrow(location="tr",which_north="true",style=north_arrow_fancy_orienteering()) + theme_bw()+
  coord_sf(xlim = c(796752, 1600602), ylim = c(269573.6 , 961154.4)) + geom_tile(data=map,aes(x=x,y=y,fill=Accuracy_AL2)) + scale_fill_viridis_c() + 
  ggtitle("Accuracy: Alert Level 2") + labs(fill="Accuracy")

e = fl_lakes %>% filter(COMID==16630838) %>% extent()

df = map %>% filter(COMID==16630838)


inset_map <- ggplot() +
  geom_sf(data=fl_lakes %>% filter(COMID==16630838)) +
  geom_tile(data=df, aes(x=x,y=y,fill=Accuracy_AL2))  + theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_sf(xlim = c(e[1] , e[2]), ylim = c(e[3], e[4])) + xlab("") + ylab("") + scale_fill_viridis_c()+ theme(legend.position="none")



combined_map <- ggdraw() +
  draw_plot(main_map) +
  draw_plot(inset_map, x = 0.3, y = 0.3, width = 0.25, height = 0.25)

ggsave("florida_map_AL1_RF.png",plot=combined_map,width=8,height=6,dpi=300)


#### Figure 6. Correlation and Error plot ####
trm = read_csv("/work/HAB4CAST/max_beal/CyANPixelForecast/output/correlation_errors_95MAE.csv")
xlabs <- c("Chl-a 1 week lag", "Chl-a 2 week lag", "Chl-a 3 week lag", "Chl-a 4 week lag", "Neighboring Pixels Chl-a","Precipitation","Surface Radiation", "Water temperature",
           "Wind Speed")
corplot = ggplot(trm) + geom_boxplot(aes(y=value,x=variable, fill=higherror),outlier.size = 0.5) + theme_classic() + scale_fill_grey(start = 0.8, end = 0.4) + 
  geom_abline(slope=0,intercept=0, linetype="dashed") + scale_x_discrete(labels= xlabs) + ylab("Spearman Rank Correlation") + 
  theme(axis.text.x=element_text(color = "black", size=10, angle=30, vjust=.8, hjust=0.8)) +
  labs(fill = "MAE > 95%") + xlab("")

ggsave("correlation_error_plot.png",plot=corplot,width=8,height=6,dpi=300)

#### Figure 6. Comparison of Pixel and Lake level forecast ####
fl_lakes$COMID = as.numeric(fl_lakes$COMID)


pi=read.csv("/work/HAB4CAST/max_beal/CyANPixelForecast/output/lake_pixel_comparison.csv")
pi=left_join(pi,fl_lakes,by=c("COMID"))

mean(pi$PixelImprovement)

ann_text<-data.frame( 
  x = c(-0.6,0.6), y = c(65,65), 
  label = c("Lake Model More Accurate",
            "Pixel Model More Accurate")
  
) 

lake_pix_compare = ggplot(pi, aes(y = reorder(LAKE_, -PixelImprovement), x = PixelImprovement)) +
  geom_bar(stat = "identity",fill="gray",color="black") +
  labs(title="Lake and Pixel Resolved Exceedances",x = "Accuracy Difference", y = "Lake Name") + xlim(-1,1) + theme_minimal() +
  geom_text(data=ann_text,aes(x=x,y=y,label=label),size=4)+
  theme(axis.text.y = element_text(size=8))

ggsave("lake_pixel_comparison.png",plot=lake_pix_compare,width=8,height=11,dpi=300)






