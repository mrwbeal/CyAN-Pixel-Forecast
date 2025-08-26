#This notebook looks at temporal autocorrelation at the pixel level across Florida
#We currently have 1 lag in the model, but more may be beneficial
#This will also generate a figure that can go along with spatial correlaiton
#Author: Max Beal
#Date Created: 1/8/2025
setwd("~/Desktop/pixel_forecast_ms")
#### Libraries ####
library(readr)
library(dplyr)
library(ggplot2)
library(sf)
library(tidyverse)
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
library(ggplot2)


#### Read in the data frame ####
df=read_csv("florida_data_spatial.csv")



comids = df$COMID %>% unique()

test = df[1:10000,]


getacf = function(chl,lag){
  
  if (length(chl)>2) {
    ac = acf(chl,plot=FALSE)[[1]][lag]
  }
  if (length(chl)<=2) {
    ac=NA
  }
  
  return(ac)
  
}

#Get observed probabilities
ac = df %>% group_by(x,y) %>% summarise(ac1 = getacf(chl,lag=2), #2 is lag 1, etc.
                                        ac2 = getacf(chl,lag=3),
                                        ac3 = getacf(chl,lag=4),
                                        ac4 = getacf(chl,lag=5),
                                        ac5 = getacf(chl,lag=6),
                                        ac6 = getacf(chl,lag=7),
                                        ac7 = getacf(chl,lag=8),
                                        ac8 = getacf(chl,lag=9),
                                        ac9 = getacf(chl,lag=10),
                                        ac10 = getacf(chl,lag=11),
                                        ac11 = getacf(chl,lag=12),
                                        ac12 = getacf(chl,lag=13),
                                        ac13 = getacf(chl,lag=14),
                                        ac14 = getacf(chl,lag=15),
                                        ac15 = getacf(chl,lag=16)) 


acm = melt(ac,id.vars=c("x","y"))


a=ggplot(data=acm) + geom_boxplot(aes(x=variable,y=value,group=variable)) + 
  theme_minimal() + geom_abline(slope=0,intercept=0,linetype="dashed",color="red") + ylab("Pearson Correlation") + xlab("Lag (weeks)")


#write.csv(acm,"/work/HAB4CAST/max_beal/CyANPixelForecast/output/temporal_autocor.csv")

#ggsave("TemporalAutoCor.png",plot=a,width=12,height=6,dpi=300)


colMeans(ac[,-c(1,2)],na.rm = T)

sapply(ac[,-c(1,2)], function(col) {
  if (is.numeric(col)) median(col, na.rm = TRUE) else NA
})



write_csv(ac,"temporal_autocor.csv")




