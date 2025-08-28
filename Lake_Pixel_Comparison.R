# This code is created to look at differences in how many events a pixel-level CyAN metric catches compared to the aggregates lake value
# This is not bringing in forecast results currently, just identifying where events are happening that are not captured by the Lake level

library(readr)
library(data.table)
library(dplyr)
library(ggplot2)
library(lubridate)
library(stringr)

#read datasets
lake=read_csv("/work/HAB4CAST/max_beal/CyANPixelForecast/data/mean_cyan_lake_conus.csv")
pixelfl = fread("/work/HAB4CAST/max_beal/CyANPixelForecast/data/florida_data.csv")
pixelfl = pixelfl %>% select(COMID,time,x,y,value,bloom)                   



flcomid = unique(pixelfl$COMID)

fllake = lake %>% filter(COMID %in% flcomid) %>% mutate(time = paste0(week,year))


cicyano = 10^(fllake$median*0.011714-4.1870866)
fllake$chl = 6620*cicyano
fllake$chl = fllake$chl - min(fllake$chl,na.rm=T)


lakebloom = fllake %>% mutate(AL1_lake = as.numeric(chl>=12 & chl<24),
                              AL2_lake = as.numeric(chl>=24))%>% select(COMID,time,AL1_lake,AL2_lake)


lakebloom$COMID=as.numeric(lakebloom$COMID)
lakebloom$time=as.numeric(lakebloom$time)


#Some notes:
#A lake is categorized as bloom if the median value of the lake exceeds 130
#A pixel is categorized as bloom if it exceeds 230
#This means that just by going to the bloom level we are capturing more events

#So what do we want to know?
#1. how many more bloom events we find at the pixel level compared to the lake level
#2. how often a lake level metric misses a pixel level bloom event
#3. When the lake level and pixel level both show bloom events (these should be 1:1), how well does each forecast perform?

#There are two Alert Levels that we can work with
#CyAN has a lower detection limit of ~7 ug/L, so we want to focus on Alert level 1 (12-24 ug/L) and Alert level 2 (24+ ug/L)


#Get any bloom pixel fl
cicyano = 10^(pixelfl$value*0.011714-4.1870866)
pixelfl$chl = 6620*cicyano
pixelfl$chl = pixelfl$chl - min(pixelfl$chl,na.rm=T)



pixelbloom = pixelfl %>% mutate(AL1 = as.numeric(chl>=12 & chl<24),
                                AL2 = as.numeric(chl>=24)) 


pixelbloom = pixelbloom %>% group_by(COMID,time) %>% summarise(AL1_pixel=sum(AL1,na.rm=T),AL2_pixel=sum(AL2,na.rm=T)) %>% select(COMID,time,AL1_pixel,AL2_pixel)



#Join the lakes and times in which any pixel has a bloom and any lake has a bloom
#Note, the pixel comparison at the lake scale can have both AL1 and AL2
bloom_comp = left_join(pixelbloom,lakebloom)



#When pixel-level 1 occurs,
bloom_comp %>% filter(AL1_lake==1) %>% summarise(AL1=sum(AL1_pixel>0))



bloom_comp2 = bloom_comp %>% mutate(AL1_both = as.numeric(AL1_pixel>=1 & AL1_lake==1),
                                                        AL2_both = as.numeric(AL2_pixel>=1 & AL2_lake==1),
                                                        NC_both = as.numeric(AL2_pixel==0 & AL2_lake==0 & AL1_pixel==0 & AL1_lake==0),
                                                        AL1_pixel = as.numeric(AL1_pixel>=1 & AL1_lake==0),
                                                        AL2_pixel = as.numeric(AL2_pixel>=1 & AL2_lake==0)) 




#How often does both the lake and at least one pixel record AL1, AL2, and None?
#Note that pixel can have both AL1, AL2 and None 

bloom_comp_bylake = bloom_comp2 %>% group_by(COMID) %>% summarise(AL1_both_prop=mean(AL1_both,na.rm=T),
                                             AL2_both_prop=mean(AL2_both,na.rm=T),
                                             NC_both_prop=mean(NC_both,na.rm=T),
                                             AL1_pixel_prop=mean(AL1_pixel,na.rm=T),
                                             AL2_pixel_prop=mean(AL2_pixel,na.rm=T))



hist(bloom_comp_bylake$AL1_both_prop)

bloom_comp_bylake %>% select(AL1_both_prop,AL2_both_prop,NC_both_prop,AL1_pixel_prop,AL2_pixel_prop)


#By definition, all of the lake AL1 and AL2 events should have a corresponding pixel level alert. The question is, how many pixel-level events go unnoticed?

weekswithpixelexceedance = bloom_comp2 %>% filter(AL1_pixel>0 | AL2_pixel>0) %>% nrow()
print(paste0("At least 1 pixel in a lake records an Alert level 1 or 2 in ",weekswithpixelexceedance," lake-weeks of ",nrow(bloom_comp2), " total lake-weeks"))

weekswithlakeexceedance = bloom_comp2 %>% filter(AL1_lake>0 | AL2_lake>0) %>% nrow()
print(paste0("Lakes record an Alert level 1 or 2 in ",weekswithlakeexceedance," lake-weeks of ",nrow(bloom_comp2), " total lake-weeks"))

print(paste0("There are ",weekswithpixelexceedance-weekswithlakeexceedance, " lake-weeks where a pixel records an exceedance that is not captured by the lake metric"))

pctmissed = round(((weekswithpixelexceedance-weekswithlakeexceedance)/weekswithpixelexceedance)*100)

print(paste0("This makes up ",pctmissed,"% of exceedance events recorded at the pixel-level" ))




#Now, what are the characteristics of these missed events?
WeeksWithMissedEvent = bloom_comp2 %>% filter(AL1_pixel>0 | AL2_pixel>0) %>% filter(AL1_lake==0 & AL2_lake==0) %>% select(COMID,time,AL1_pixel,AL2_pixel)

missed = left_join(WeeksWithMissedEvent,pixelfl)

#### How many pixels/how much of the lake is blooming on average in these events? ####
size_missed = missed %>% group_by(COMID, time) %>% summarise(AL1pixels=sum(chl>=12,na.rm=T),
                                               totalpixels=sum(!is.na(x)),
                                               PctLake = AL1pixels/totalpixels,
                                               group="Excluded by Lake Level Forecast")    

#Using the median should mean that at least half the lake has to be above 12 ug/L
#so misses shouldn't have large bloom events
hist(size_missed$PctLake)



#Look at all pixel events
size_all=pixelfl %>% group_by(COMID, time) %>% summarise(AL1pixels=sum(chl>=12,na.rm=T),
                                                totalpixels=sum(!is.na(x)),
                                                PctLake = AL1pixels/totalpixels,
                                                group="All Events") %>% filter(AL1pixels>0)



size=rbind(size_missed,size_all)

#### plot a ####
a = ggplot(data=size)+geom_density(aes(x=PctLake,group=group,fill=group),alpha=0.5) + theme_bw() + xlab("Percent of Lake Area") + ylab("Density") +
  theme(legend.title=element_blank())

missed %>% na.omit()

#### What is the magnitude of these missed events compared to all events? ####
#We want to know, among the pixels that bloomed, what was the average chl for missed events vs all events?
magnitude_missed = missed %>% group_by(COMID, time,x,y)%>% filter(chl>=12) %>% summarise(avchl=mean(chl,na.rm=T)) %>% mutate(group="Excluded by Lake Level Forecast")
magnitude_all = pixelfl %>% group_by(COMID, time,x,y) %>% filter(chl>=12) %>% summarise(avchl=mean(chl,na.rm=T))  %>% mutate(group="All Events")
magnitude=rbind(magnitude_missed,magnitude_all)


#### plot b ####
b = ggplot(data=magnitude)+geom_density(aes(x=log(avchl),group=group,fill=group),alpha=0.5) + theme_bw() + xlab("log(Chlorophyll-a) (μg/L)") + ylab("")+
  theme(legend.title=element_blank())

ggplot(data=magnitude)+geom_boxplot(aes(y=log(avchl),x=group,fill=group)) + theme_bw()

#### What is the average duration of these missed events? ####
padded_numbers <- sprintf("%06d", missed$time)
missed$week = substr(padded_numbers,1,2)
missed$year = substr(padded_numbers,3,6)
duration_missed=missed %>%
  mutate(beginning = ymd(str_c(year, "-01-01")), date = beginning + weeks(week)) %>% mutate(AL1 = as.numeric(chl>=12))

duration_missed = duration_missed %>% group_by(x,y) %>% arrange(date) %>% mutate(grp=cumsum(c(0,diff(AL1))<0)) %>% filter(AL1==1) %>% group_by(x,y,grp) %>% 
  summarise(weeks_above_12=n(),.groups = "drop_last",start_week=min(as.numeric(week)),end_week=max(as.numeric(week))) %>% mutate(group="Excluded by Lake Level Forecast")

pixelfl=pixelfl %>% mutate(AL1 = as.numeric(chl>=12))

padded_numbers <- sprintf("%06d", pixelfl$time)
pixelfl$week = substr(padded_numbers,1,2)
pixelfl$year = substr(padded_numbers,3,6)
duration_pixelfl=pixelfl %>%
  mutate(beginning = ymd(str_c(year, "-01-01")), date = beginning + weeks(week))

duration_pixelfl = duration_pixelfl %>% group_by(x,y) %>% arrange(date) %>% 
  mutate(grp=cumsum(c(0,diff(AL1))<0)) %>% #grp is an identifier for each event
  filter(AL1==1) %>% group_by(x,y,grp) %>% 
  summarise(weeks_above_12=n(),.groups = "drop_last",start_week=min(as.numeric(week)),end_week=max(as.numeric(week))) %>% 
  mutate(group="All Events")




duration_pixelfl %>% na.omit()

duration=rbind(duration_missed,duration_pixelfl)


#### plot c ####
c = ggplot(data=duration)+geom_density(aes(x=weeks_above_12,group=group,fill=group),alpha=0.5) + theme_bw() + xlim(0,100) + xlab("Weeks above Alert Level 1")+
  theme(legend.title=element_blank())


duration$weeks_above_12 %>% mean()


library(ggpubr)

# It looks like the missed events tend to be 1) smaller 2) of less magnitude and 3) more transient.
a = ggplot(data=size)+geom_density(aes(x=PctLake*100,group=group,fill=group),alpha=0.5) + theme_bw() + xlab("Percent of Lake Area") + ylab("Density") +
  theme(legend.position='bottom',legend.title=element_blank(), legend.text = element_text(size=9.5),legend.spacing.x = unit(2,"cm"))+ ggtitle("Size")+ scale_fill_grey() +
  labs(title="Size",tag="A")
b = ggplot(data=magnitude)+geom_density(aes(x=log(avchl),group=group,fill=group),alpha=0.5) + theme_bw() + xlab("log(Chlorophyll-a) (μg/L)") + ylab("")+
  theme(legend.position='bottom',legend.title=element_blank(), legend.text = element_text(size=9.5),legend.spacing.x = unit(2,"cm"))+ ggtitle("Magnitude")+ scale_fill_grey() +
  labs(title="Magnitude",tag="B")
c = ggplot(data=duration)+geom_density(aes(x=weeks_above_12,group=group,fill=group),alpha=0.5) + theme_bw() + xlab("Weeks above Alert Level 1")+
  theme(legend.position='bottom',legend.title=element_blank(), legend.text = element_text(size=9.5),legend.spacing.x = unit(2,"cm"))+ ylab("") + ggtitle("Duration") + scale_fill_grey() + xlim(0,53) +
  labs(title="Duration",tag="C")

grid = ggarrange(a,b,c,nrow=1,common.legend=T,legend="bottom")

ggsave("Size_Mag_Dur.png",plot=grid,width=13.5,height=4.5,dpi=600)

a=ggplot(data=size,aes(y=PctLake,x=group)) + theme_bw() + ylab("Percent of Lake Area")+
  geom_violin(position=position_dodge()) +
  geom_boxplot(width=0.05,color="black",position = position_dodge(width =0.5)) + xlab("")+ ggtitle("Size")


b=ggplot(data=magnitude,aes(y=log(avchl),x=group)) + theme_bw()+ ylab("log(Chlorophyll-a) (μg/L)")+
  geom_violin(position=position_dodge()) +
  geom_boxplot(width=0.05,color="black",position = position_dodge(width =0.5)) + xlab("")+ ggtitle("Magnitude")


c=ggplot(data=duration,aes(y=weeks_above_12,x=group))+ theme_bw() + ylab("Weeks above Alert Level 1")+
  geom_violin(position=position_dodge()) +
  geom_boxplot(width=0.05,color="black",position = position_dodge(width =0.5)) + xlab("") + ggtitle("Duration") + ylim(0,53)


ggarrange(a,b,c,nrow=1)


duration %>% filter(weeks_above_12>53) %>% nrow()/nrow(duration)

size %>% group_by(group) %>% summarise(mean(PctLake), median(AL1pixels))

magnitude %>% group_by(group) %>% summarise(mean(avchl))

duration %>% group_by(group) %>% summarise(median(weeks_above_12), median(as.numeric(start_week),na.rm=T),median(as.numeric(end_week),na.rm=T))

duration$start_end_week = paste0(duration$start_week,"-",duration$end_week)


start_end = duration %>% group_by(group,start_end_week) %>% summarise(n())

alperweek = pixelfl %>% group_by(week) %>% summarise(totalal=sum(AL1,na.rm=T))

plot(alperweek$week,alperweek$totalal)

