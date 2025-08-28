# This notebook is for exploring the source of errors in the pixel level forecast, supporting the manuscript
#Author Max Beal
#Date created: 2/12/2025

library(readr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(data.table)

dfm=read_csv("/work/HAB4CAST/max_beal/CyANPixelForecast/output/RFMetricsLead1.csv")

df = read_csv("/work/HAB4CAST/max_beal/CyANPixelForecast/output/RF_folds_all.csv")


#First, find the distribution of errors
hist(df$error)
summary(df$error)

err90 = quantile(abs(df$error),0.9)



err90 = quantile(dfm$MAE,0.95)

#Really only shows a huge difference at the largest errors
dfm$ErrorCat <- cut(dfm$MAE, quantile(dfm$MAE, 0:8/8))



#Start by looking at absolute error above the 90th percentile

dfm$higherror = dfm$MAE>err90
df = left_join(df,dfm,by=c("x","y"))

tr = df %>% group_by(x,y,higherror) %>% summarise(chlcor = cor(chl_lead,chl,method="spearman"),
                                                  chllag2cor = cor(chl_lead,chl_lag2,method="spearman"),
                                                  chllag3cor = cor(chl_lead,chl_lag3,method="spearman"),
                                                  chllag4cor = cor(chl_lead,chl_lag4,method="spearman"),
                                                   tempcor = cor(chl_lead,PRISMARDSW,method="spearman"),
                                                   wscor = cor(chl_lead,windSpeed,method="spearman"),
                                                   srcor = cor(chl_lead,surfaceRadiation,method="spearman"),
                                                   neighchlcor=cor(chl_lead,mean_neighbor_chl,method="spearman"),
                                                  prcor=cor(chl_lead,precip,method="spearman"))



trm = tr %>% dplyr::select(higherror,chlcor,chllag2cor,chllag3cor,chllag4cor,neighchlcor,tempcor,wscor,srcor,prcor) %>% melt(id.vars=c("x","y","higherror"))


xlabs <- c("Chl-a 1 week lag", "Chl-a 2 week lag", "Chl-a 3 week lag", "Chl-a 4 week lag", "Neighboring Pixels Chl-a", "Water temperature",
                    "Wind Speed"," Surface Radiation","Precipitation")
corplot = ggplot(trm) + geom_boxplot(aes(y=value,x=variable, fill=higherror)) + theme_classic() + scale_fill_grey(start = 0.8, end = 0.4) + 
  geom_abline(slope=0,intercept=0, linetype="dashed") + scale_x_discrete(labels= xlabs) + ylab("Spearman Rank Correlation") + 
  theme(axis.text.x=element_text(color = "black", size=11, angle=30, vjust=.8, hjust=0.8)) +
  labs(fill = "MAE > 95%")

#write_csv(trm,"/work/HAB4CAST/max_beal/CyANPixelForecast/output/correlation_errors_95MAE.csv")


trm %>% group_by(variable,higherror) %>% summarise(avcor=mean(value,na.rm=T))

he = df %>% group_by(higherror) %>% summarise(sum(error>0)/n(),
                                              max(chl_observed),
                                              mean(Week))

clim = df %>% group_by(higherror, Week) %>% summarise(mchl = mean(chl_observed),
                                                      vchl = var(chl_observed),
                                                      maxchl=max(chl_observed),
                                                      mchlp=mean(chl_predicted),
                                                      mchllta=mean(chl_LTA,na.rm=T))


ggplot(clim) + geom_line(aes(x=Week,y=mchl,color=higherror)) + geom_line(aes(x=Week,y=mchlp,color=higherror),linetype="dashed") + 
  geom_point(aes(x=Week,y=mchllta,color=higherror))



ggplot(clim) + geom_line(aes(x=Week,y=vchl,color=higherror))

sum(he$error<0)/nrow(he)

ggplot(df) + geom_density(aes(x=error,fill=higherror))


trw = df %>% group_by(Week,higherror) %>% summarise(chlcor = cor(chl_lead,chl,method="spearman"),
                                                  chllag2cor = cor(chl_lead,chl_lag2,method="spearman"),
                                                  chllag3cor = cor(chl_lead,chl_lag3,method="spearman"),
                                                  chllag4cor = cor(chl_lead,chl_lag4,method="spearman"),
                                                  tempcor = cor(chl_lead,PRISMARDSW,method="spearman"),
                                                  wscor = cor(chl_lead,windSpeed,method="spearman"),
                                                  srcor = cor(chl_lead,surfaceRadiation,method="spearman"),
                                                  neighchlcor=cor(chl_lead,mean_neighbor_chl,method="spearman"),
                                                  prcor=cor(chl_lead,precip,method="spearman"))


ggplot(trw) + geom_line(aes(x=Week,y=chlcor,color=higherror)) 

trw %>% filter(Week>=16 & Week<=26 & higherror==T) %>% pull(chlcor) %>% mean()

trw %>% filter(higherror==T) %>% pull(chlcor) %>% mean()

cplot = df %>% filter(higherror==TRUE)

ggplot(cplot) + geom_point(aes(x=chl_lead,y=chl),size=0.5) + facet_wrap(~Week)

climcor=merge(clim,trw)

higherrorlakes = cplot$COMID %>% unique()

#### data ####
library(tidyr)
library(stringr)
df <- fread("/work/HAB4CAST/max_beal/CyANPixelForecast/data/florida_data_spatial.csv")
df$Time = paste0(df$Year,str_pad(df$Week,2,pad=0))
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
input = input %>% filter(input$PRISMARDSW %>% is.na()) %>% mutate(PRISMARDSW = PRISMSW) 


he = left_join(df,dfm,by=c("x","y"))

he=input %>% mutate(higherror=COMID %in% higherrorlakes)

ggplot(he) + geom_boxplot(aes(x=higherror,y=PRISMARDSW)) 

he %>% group_by(higherror) %>% summarise(min(PRISMARDSW,na.rm=T),mean(chl,na.rm=T))

ggplot(he) + geom_boxplot(aes(x=higherror,y=precip)) 

ggplot(he, aes(x=chl,y=PRISMSW)) + geom_hex() + facet_wrap(~higherror)

ggplot(he, aes(x=surfaceRadiation,y=chl)) + geom_hex(bins=60,stat='binhex') + scale_fill_viridis_c(trans="log",name="Density")+ facet_wrap(~higherror)


ggplot(he, aes(x=windSpeed,y=chl_lead)) + geom_hex(bins=60,stat='binhex') + scale_fill_viridis_c(trans="log",name="Density")

ggplot(he, aes(x=PRISMARDSW,y=chl_lead)) + geom_hex(bins=60,stat='binhex') + scale_fill_gradient2(low="red",high="blue",name="Density")

he


