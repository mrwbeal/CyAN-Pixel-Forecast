# This script is for playing with ways to visualize the continuous results in terms of different categories

library(dplyr)
library(data.table)
library(sf)

fileresults=list()
for (i in c(1:7)) {
  fileresults[[i]]=paste0("/work/HAB4CAST/max_beal/CyANPixelForecast/output/FL_RF_fold_",i,"_lag1.csv")
}




folds = lapply(fileresults,fread) %>% bind_rows()



lakes <- st_read("/work/HAB4CAST/max_beal/SW_model/data/OLCI_resolvable_lakes_2022_09_08/OLCI_resolvable_lakes_2022_09_08.shp")
lake_crs = st_crs(lakes)

plot(lakes)


l=99
lake = folds %>% filter(COMID==unique(folds$COMID)[l])
lakeshp = lakes %>% filter(COMID==unique(folds$COMID)[l])
ggplot() + geom_tile(data=lake,aes(x=x,y=y),color="lightgray") + geom_sf(data=lakeshp,fill=NA) + 
  geom_point(data=lake %>% filter(x== 1402353,y== 816276),aes(x=x,y=y),color="green")

lake$x %>% unique()
lake$y %>% unique()

cell = lake %>% filter(x== 1402353,y== 815376)


#### Ensemble Prediction ####
edf = cell %>% dplyr::select(x,y,Time,contains("member"))


edf=melt(edf,id.vars = c("x","y","Time"))


preds = edf %>% group_by(x,y,Time) %>% summarise(minpred = min(value),
                                                 maxpred = max(value),
                                                 predicted_chl = mean(value))

obs = cell %>% dplyr::select(x,y,Time,chl_observed)
cellp = merge(obs,preds)
padded_numbers <- sprintf("%06d", cellp$Time)

cellp$week = substr(padded_numbers,1,2)
cellp$year = substr(padded_numbers,3,6)

cellp=cellp %>%
  mutate(beginning = ymd(str_c(year, "-01-01")),
         date = beginning + weeks(week))

ggplot(data=cellp) + 
  geom_ribbon(aes(x=date,ymin=0,ymax=12), fill="lightgreen",alpha=0.8) + 
  geom_ribbon(aes(x=date,ymin=12,ymax=50), fill="lightyellow",alpha=0.8) + 
  geom_ribbon(aes(x=date,ymin=50,ymax=350), fill="indianred",alpha=0.8) + 
  geom_ribbon(aes(x=date,ymin=minpred, ymax=maxpred),color="black",fill=NA,linetype="dashed", linewidth=0.5) + 
  geom_line(aes(x=date,y=predicted_chl),linewidth=0.8) + 
  geom_line(aes(x=date,y=chl_observed),color="green",linewidth=0.8) +
  ylab("Predicted Chloropyll-a") + xlab("Time") +
  theme_bw()


#### Boxplot Example ####
edf = cell %>% dplyr::select(x,y,Time,chl_observed,contains("member"))
cellp=melt(edf,id.vars = c("x","y","Time","chl_observed"))

cellp$Time

padded_numbers <- sprintf("%06d", cellp$Time)
cellp$week = substr(padded_numbers,1,2)
cellp$year = substr(padded_numbers,3,6)
cellp=cellp %>%
  mutate(beginning = ymd(str_c(year, "-01-01")),
         date = beginning + weeks(week))

ggplot(data=cellp %>% filter(date>"2023-01-01" & date<"2024-01-01")) + 
  geom_ribbon(aes(x=date,ymin=0,ymax=12), fill="lightgreen",alpha=0.8) + 
  geom_ribbon(aes(x=date,ymin=12,ymax=50), fill="lightyellow",alpha=0.8) + 
  geom_ribbon(aes(x=date,ymin=50,ymax=350), fill="indianred",alpha=0.8) + 
  geom_boxplot(aes(x=date,y=value,group=date),linewidth=1) + 
  geom_line(aes(x=date,y=chl_observed),color="cyan",linewidth=1) +
  ylab("Predicted Chloropyll-a") + xlab("Time") +
  theme_bw() + ggtitle("2022 Forecast Results")


#breaks = c(0,12,50,100, 1000)

breaks = c(0,1,12,24,1000)
labels = c("nc","who1","who2","who3")
cellp$categories = cut(cellp$value,breaks,labels)
percents = cellp %>% group_by(date,categories) %>% tally %>% 
  spread(categories, n, fill = 0)

percents

cellp$categories_obs = cut(cellp$chl_observed,breaks,labels)

(table(cellp$categories,cellp$categories_obs)/nrow(cellp))*100


percents = cellp %>% group_by(date,categories) %>% tally %>% complete(n)

percents$correct = percents$categories==percents$categories_obs

ggplot(data=percents%>% filter(date>"2022-01-01" & date<"2023-01-01")) + 
  geom_bar(aes(x=date,y=n,fill=forcats::fct_rev(categories),color=correct),position = "stack", stat = "identity") +
  labs(x = "Date", y = "Percent", title = "Hindcast of Percent Chance WHO Cyanobacteria Categories: 2022") +
  theme_minimal() + scale_fill_manual(values=c("indianred","khaki","lightgreen","skyblue")) + 
  scale_color_manual(values=c("black","white"),na.value = NA)


ggplot(data=percents%>% filter(date=="2024-06-24")) + 
  geom_bar(aes(x="",y=n,fill=forcats::fct_rev(categories),color=correct),position = "stack", stat = "identity") +
  labs(title = "Percent change WHO cyano exceedance: 2024-06-24") +
  theme_minimal() + scale_fill_manual(values=c("indianred","khaki","lightgreen","skyblue")) + 
  scale_color_manual(values=c("black"),na.value = NA) + coord_polar(theta="y")


#### try lake vis

edf = lake %>% dplyr::select(x,y,Time,contains("member"))
edf=melt(edf,id.vars = c("x","y","Time"))

breaks = c(0,1,12,24,100)
labels = c("nc","who1","who2","who3")
edf$categories = cut(edf$value,breaks,labels)

padded_numbers <- sprintf("%06d", edf$Time)
edf$week = substr(padded_numbers,1,2)
edf$year = substr(padded_numbers,3,6)
edf=edf %>%
  mutate(beginning = ymd(str_c(year, "-01-01")),
         date = beginning + weeks(week))

percents = edf %>% group_by(date,x,y,categories) %>% tally %>% complete(n)

library(ggforce)


percents$date %>% unique()
day = percents %>% filter(date=="2024-06-17")


ggplot() + geom_sf(data=lakeshp,fill=NA) + geom_tile(data=day, aes(x=x,y=y),fill="lightgray",color="darkgray") + 
  geom_arc_bar(data=day,aes(x0=x,y0=y,amount=n,r0=0,r=100,fill=categories),stat='pie') + theme_bw() +
  scale_fill_manual(values=c("skyblue","khaki","orange","indianred")) + ggtitle("2024-06-17, COMID 21483782")


padded_numbers <- sprintf("%06d", lake$Time)
lake$week = substr(padded_numbers,1,2)
lake$year = substr(padded_numbers,3,6)
lake=lake %>%
  mutate(beginning = ymd(str_c(year, "-01-01")),
         date = beginning + weeks(week))

ggplot(data=lake %>% filter(date>="2024-06-17")) +
  geom_ribbon(aes(x=date,ymin=0,ymax=12), fill="skyblue",alpha=0.6) + 
  geom_ribbon(aes(x=date,ymin=12,ymax=50), fill="khaki",alpha=0.6) + 
  geom_ribbon(aes(x=date,ymin=50,ymax=100), fill="orange",alpha=0.6) + 
  geom_ribbon(aes(x=date,ymin=100,ymax=120), fill="indianred",alpha=0.6)+ 
    geom_point(aes(x=date,y=chl_observed,group=date)) + ylab("Observed Chlorophyll-a") + theme_bw()

breaks = c(0,1,12,24,100)
labels = c("nc","who1","who2","who3")
lake$categories_obs = cut(lake$chl_observed,breaks,labels)

ggplot() + geom_sf(data=lakeshp,fill=NA) + geom_tile(data=lake %>% filter(date=="2024-07-01"), aes(x=x,y=y,fill=categories_obs),color="gray")+
  scale_fill_manual(values=c("skyblue","khaki","orange","indianred")) +theme_bw()+ ggtitle("2024-07-01, COMID 21483782, Observed Chl categories")




