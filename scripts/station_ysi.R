# A script for calculating running linear models for satellite and in situ data
# Summer 2022 reservoirs. Sample~YSI

# clear workspace
rm(list=ls())

# 1. load packages
library(janitor)
library(lubridate)
library(readxl)
library(plyr)
library(tidyverse)
library(reshape2)
library(data.table)
library(stringr)
library(ggplot2)
library(rgdal)
library(grid)
library(gridExtra)
library(lattice)
library(rioja)
library(zoo)
library(ggpubr)
library(ggpmisc)


# 2. Read in data files ------------------------------------------

# Exo 2 data

table0<-fread("FLAMe_EXO.dat",skip=1)
table<-table0 %>% as.data.frame()
table_clean<-clean_names(table)
table_clean<-table_clean[-c(1:2),]

exo<-table_clean
cols_numeric<-names(exo)[5:21]
exo<-exo %>% mutate_at(cols_numeric, as.numeric)

exo %>% mutate_if(is.character,as.numeric)

# Correct for lag time between intake and sensor using an approximation
# the Exo 2 flow cell has a volume of 1L
exo$exo_datetime_utc<-as.POSIXct(exo$timestamp,tz="UTC")-21 # lag time correction in seconds 
exo$exo_datetime_local<-format(exo$exo_datetime_utc,tz="America/Chicago")
exo$lake_time_local<-substr(exo$exo_datetime_local,12,19)
names(exo)[which(names(exo)=="temp")]<-"temp_c"
names(exo)[which(names(exo)=="chlor_rfu")]<-"chlorophyll_rfu"
names(exo)[which(names(exo)=="chlor_ug_l")]<-"chlorophyll_ugl"
names(exo)[which(names(exo)=="bg_apc_rfu")]<-"tal_pc_rfu"
names(exo)[which(names(exo)=="turb_fnu")]<-"turbidity_fnu"
names(exo)[which(names(exo)=="f_dom_rfu")]<-"f_dom_rfu"
names(exo)[which(names(exo)=="odo_mg_l")]<-"odo_mg_l"
names(exo)[which(names(exo)=="odo_percent")]<-"odo_percent"
names(exo)[which(names(exo)=="temp_c")]<-"temp_c"

exo$date <- substr(exo$exo_datetime_local, 1, 10)

# Garmin GPS puck data
table0<-fread("FLAMe_GPS.dat")#,skip=c(2:3))
table<-table0 %>% as.data.frame()
table<-table[-c(2:3),]
names(table)<-table[1,]
table<-table[-1,]
table$garmin_datetime_utc<-as.POSIXct(table$TIMESTAMP,tz="UTC")
table$garmin_datetime_local<-format(table$garmin_datetime_utc,tz="America/Chicago")
table_clean<-clean_names(table)
garmin<-table_clean
garmin$garmin_time_local<-substr(garmin$garmin_datetime_local,12,19)
garmin$xcoord<-as.numeric(garmin$longitude)
garmin$ycoord<-as.numeric(garmin$latitude)
garmin$date<-substr(garmin$garmin_datetime_local,1,10)
garmin$hour_dec<-hour(garmin$garmin_datetime_local)+
  minute(garmin$garmin_datetime_local)/60+
  second(garmin$garmin_datetime_local)/3600

garmin<-garmin %>% dplyr::select(garmin_datetime_utc,date,hour_dec,garmin_time_local,xcoord,ycoord)
#getUTMzone(mean(garmin$xcoord,na.rm=TRUE))


# 3. Merge data

data_merged<-merge(exo,garmin, 
                   by.x=c("lake_time_local","date"), 
                   by.y=c("garmin_time_local", "date"))


data_merged<-data_merged %>% dplyr::select(lat_dec=ycoord,lon_dec=xcoord,
                                           time_local=lake_time_local,date,datetime_utc=garmin_datetime_utc,
                                           wtemp=temp_c, chl=chlorophyll_rfu,chl_ugl=chlorophyll_ugl,pc=tal_pc_rfu,f_dom_rfu,
                                           turb=turbidity_fnu,do_mgl=odo_mg_l, do_sat=odo_percent, hour_dec
)

data_merged$date_time_c = paste(data_merged$date, data_merged$time_local)

# filtering garmin data by lake / date
bon1<-data_merged %>% filter(date=="2022-07-14", hour_dec>13.20, hour_dec<16.733)
bon2<-data_merged %>% filter(date=="2022-07-15", hour_dec>12, hour_dec<13.91)
waco<-data_merged %>% filter(date=="2022-07-23", hour_dec>10.767, hour_dec<13.25)
bw<-data_merged %>% filter(date=="2022-08-06", hour_dec>12.91, hour_dec<15.283)
iv<-data_merged %>% filter(date=="2022-08-07", hour_dec>10.67, hour_dec<13.68)
rb<-data_merged %>% filter(date=="2022-08-08", hour_dec>10.7, hour_dec<15.967)
ah<-data_merged %>% filter(date=="2022-08-12", hour_dec>9.8, hour_dec<12.8)

# combine exo data
#exo_list <- list(bon1, bw, iv,rb, ah)

#merge all data frames in list
#exo_df <- Reduce(function(x, y) merge(x, y, all=TRUE), exo_list)

#############################################################################################################################
# CRASR data - read in for all lakes
crasr_all <- read_csv("CRASR_summer2022.csv") %>% clean_names()

# separate by lake day
crasr_bon1 <- crasr_all %>% filter(datecollect == "7/14/2022") %>% add_column(lake = "bon") %>% 
 group_by(system,lat, lon, variable, site)%>% mutate(meanchl = mean(value))
cbon1 <- crasr_bon1 %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl)

crasr_bon2 <- crasr_all %>% filter(datecollect == "7/15/2022") %>% add_column(lake = "bon") %>% 
  dplyr::mutate(site = paste("BON", site, sep = "")) %>% group_by(system,lat, lon, variable, site) %>% mutate(meanchl = mean(value))
cbon2 <- crasr_bon2 %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl)

crasr_lw <- crasr_all %>% filter(datecollect == "7/23/2022")%>% add_column(lake = "lw") %>% 
  dplyr::mutate(site = paste("LW", site, sep = "")) %>% group_by(system,lat, lon, variable, site) %>% mutate(meanchl = mean(value))
clw <- crasr_lw %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl)

crasr_bw <- crasr_all %>% filter(datecollect == "8/6/2022")%>% add_column(lake = "bw") %>% 
  group_by(system,lat, lon, variable, site) %>% mutate(meanchl = mean(value))
cbw <- crasr_bw %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl)

crasr_iv <- crasr_all %>% filter(datecollect == "8/7/2022") %>% add_column(lake = "iv") %>% 
  group_by(system,lat, lon, variable, site) %>% mutate(meanchl = mean(value))
civ <- crasr_iv %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl)

crasr_rb <- crasr_all %>% filter(datecollect == "8/8/2022") %>% add_column(lake = "rb") %>% 
  group_by(system, lat, lon, variable, site) %>% mutate(meanchl = mean(value))
crb <- crasr_rb %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl)

crasr_ah <- crasr_all %>% filter(datecollect == "8/12/2022") %>% add_column(lake = "ah") %>% 
  group_by(system, lat, lon, variable, site) %>% mutate(meanchl = mean(value))
cah <- crasr_ah %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl)

# combine crasr data
crasr_list <- list(cbon2, clw, cbw, civ,crb, cah)

#merge all data frames in list
crasr_df <- Reduce(function(x, y) merge(x, y, all=TRUE), crasr_list)

##########################################################################################################
# merge sample data with sensor turbidity data to create validation data frames for all 4 lakes
#sampsens_merged <- merge(crasr_df, exo_df)

sampsensbon <- merge(cbon2, bon2)
sampsensbon$lat_diff<-sampsensbon$lat_dec-sampsensbon$lat
sampsensbon$lon_diff<-sampsensbon$lon_dec-sampsensbon$lon
samplesensorbon_merged<-sampsensbon %>% filter(abs(lat_diff)<0.0003,abs(lon_diff)<0.0003)

valdatabon<-samplesensorbon_merged %>% dplyr::select(system, site,lat,lon,turb_lab = turbidity,
                                                         turb_ysi=turb,chla_lab=`chl a`, chla_ysi = chl, chl_ugl) %>% 
  group_by(system, site,lat,lon,turb_lab,turb_ysi,chla_lab, chla_ysi) %>%
  dplyr::summarize(turb_ysi=mean(turb_ysi,na.rm=TRUE)) %>%
  dplyr::summarize(chla_ysi=mean(chla_ysi,na.rm=TRUE)) %>% 
  as.data.frame()

valbon2_means <- valdatabon %>% group_by(site) %>%
  mutate(turb_ysi_m=mean(turb_ysi,na.rm=TRUE)) %>% 
  mutate(chla_ysi_m=mean(chla_ysi,na.rm=TRUE))


# same but for chl_ugl
valdatabonugl<-samplesensorbon_merged %>% dplyr::select(system, site,lat,lon,turb_lab = turbidity,
                                                     turb_ysi=turb,chla_lab=`chl a`, chla_ysi = chl_ugl) %>% 
  group_by(system, site,lat,lon,turb_lab,turb_ysi,chla_lab, chla_ysi) %>%
  dplyr::summarize(turb_ysi=mean(turb_ysi,na.rm=TRUE)) %>%
  dplyr::summarize(chla_ysi=mean(chla_ysi,na.rm=TRUE)) %>% 
  as.data.frame()

valbon2_meansugl <- valdatabonugl %>% group_by(site) %>%
  mutate(turb_ysi_m=mean(turb_ysi,na.rm=TRUE)) %>% 
  mutate(chla_ysi_m=mean(chla_ysi,na.rm=TRUE))


lmt_all<-lm(turb_ysi_m~turb_lab,data=valbon2_means)
summary(lmt_all) #0.9772

lmc_all<-lm(chla_ysi_m~chla_lab,data=valbon2_means) 
summary(lmc_all) #0.2027

bnt<-valbon2_means %>%
  ggplot(aes(y=turb_ysi_m,x=turb_lab, label=site))+
  geom_smooth(method="lm",lty=1,col="dark gray",se=FALSE)+
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+
  ylab("Sensor Turbidity (NTU)")+
  xlab("Sample Turbidity (NTU)")+
  theme_bw()


bnc<-valbon2_means %>%
  ggplot(aes(y=chla_ysi_m,x=chla_lab, label=site))+
  geom_smooth(method="lm",lty=1,col="dark gray",se=FALSE)+
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+
  ylab("Sensor Chl-a (RFU)")+
  xlab("Sample Chl-a (mg/L)")+
  theme_bw()


##########################################################################################################

# merge sample data with sensor turbidity data to create validation data frames
# merging for bonham

samplesensorbon_merged<-merge(bon1,cbon1)
samplesensorbon_merged$lat_diff<-samplesensorbon_merged$lat_dec-samplesensorbon_merged$lat
samplesensorbon_merged$lon_diff<-samplesensorbon_merged$lon_dec-samplesensorbon_merged$lon
samplesensorbon_merged<-samplesensorbon_merged %>% filter(abs(lat_diff)<0.0003,abs(lon_diff)<0.0003)

validationdata<-samplesensorbon_merged %>% dplyr::select(system, site,lat,lon,turb_lab = turbidity,
                                                      turb_ysi=turb,chla_lab=`chl a`, chla_ysi = chl, chl_ugl) %>% 
  group_by(system, site,lat,lon,turb_lab,turb_ysi,chla_lab, chla_ysi) %>%
  dplyr::summarize(turb_ysi=mean(turb_ysi,na.rm=TRUE)) %>%
  dplyr::summarize(chla_ysi=mean(chla_ysi,na.rm=TRUE)) %>% 
  as.data.frame()

valbon1_means <- validationdata %>% group_by(site) %>%
  mutate(turb_ysi_m=mean(turb_ysi,na.rm=TRUE)) %>% 
  mutate(chla_ysi_m=mean(chla_ysi,na.rm=TRUE))

##
validationdataugl<-samplesensorbon_merged %>% dplyr::select(system, site,lat,lon,turb_lab = turbidity,
                                                         turb_ysi=turb,chla_lab=`chl a`, chla_ysi = chl_ugl) %>% 
  group_by(system, site,lat,lon,turb_lab,turb_ysi,chla_lab, chla_ysi) %>%
  dplyr::summarize(turb_ysi=mean(turb_ysi,na.rm=TRUE)) %>%
  dplyr::summarize(chla_ysi=mean(chla_ysi,na.rm=TRUE)) %>% 
  as.data.frame()

valbon1_meansugl <- validationdataugl %>% group_by(site) %>%
  mutate(turb_ysi_m=mean(turb_ysi,na.rm=TRUE)) %>% 
  mutate(chla_ysi_m=mean(chla_ysi,na.rm=TRUE))

lmt_bon<-lm(turb_ysi_m~turb_lab,data=valbon1_means) # r2 = 0.8768
summary(lmt_bon)

lmc_bon<-lm(chla_ysi_m~chla_lab,data=valbon1_means) # r2 = 0.9006
summary(lmc_bon)

bnt<-valbon1_means %>%
  ggplot(aes(y=turb_ysi_m,x=turb_lab, label=site))+
  geom_smooth(method="lm",lty=1,col="dark gray",se=FALSE)+
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+
  ylab("Sensor Turbidity (NTU)")+
  xlab("Sample Turbidity (NTU)")+
  theme_bw()

bnc<-valbon1_means %>%
  ggplot(aes(y=chla_ysi_m,x=chla_lab, label=site))+
  geom_smooth(method="lm",lty=1,col="dark gray",se=FALSE)+
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+
  ylab("Sensor Chl-a (RFU)")+
  xlab("Sample Chl-a (mg/L)")+
  theme_bw()


# merging for brownwood
samplesensorbw_merged<-merge(bw,cbw)
samplesensorbw_merged$lat_diff<-samplesensorbw_merged$lat_dec-samplesensorbw_merged$lat
samplesensorbw_merged$lon_diff<-samplesensorbw_merged$lon_dec-samplesensorbw_merged$lon
samplesensorbw_merged<-samplesensorbw_merged %>% filter(abs(lat_diff)<0.0002,abs(lon_diff)<0.0002)

validationdatabw<-samplesensorbw_merged %>% dplyr::select(system, site,lat,lon,turb_lab = turbidity,
                                                      turb_ysi=turb,chla_lab=`chl a`, chla_ysi = chl, chl_ugl) %>% 
  group_by(system, site,lat,lon,turb_lab,turb_ysi,chla_lab, chla_ysi) %>%
  dplyr::summarize(turb_ysi=mean(turb_ysi,na.rm=TRUE)) %>%
  dplyr::summarize(chla_ysi=mean(chla_ysi,na.rm=TRUE)) %>% 
  as.data.frame()

valbw_means <- validationdatabw %>% group_by(site) %>%
  mutate(turb_ysi_m=mean(turb_ysi,na.rm=TRUE)) %>% 
  mutate(chla_ysi_m=mean(chla_ysi,na.rm=TRUE))

##\
validationdatabwugl<-samplesensorbw_merged %>% dplyr::select(system, site,lat,lon,turb_lab = turbidity,
                                                          turb_ysi=turb,chla_lab=`chl a`, chla_ysi = chl, chl_ugl) %>% 
  group_by(system, site,lat,lon,turb_lab,turb_ysi,chla_lab, chla_ysi) %>%
  dplyr::summarize(turb_ysi=mean(turb_ysi,na.rm=TRUE)) %>%
  dplyr::summarize(chla_ysi=mean(chla_ysi,na.rm=TRUE)) %>% 
  as.data.frame()

valbw_meansugl <- validationdatabwugl %>% group_by(site) %>%
  mutate(turb_ysi_m=mean(turb_ysi,na.rm=TRUE)) %>% 
  mutate(chla_ysi_m=mean(chla_ysi,na.rm=TRUE))

lmt_bw<-lm(turb_ysi_m~turb_lab,data=valbw_means) # r2 = 0.4808
summary(lmt_bw)

lmc_bw<-lm(chla_ysi_m~chla_lab,data=valbw_means) # r2 = 0.5839
summary(lmc_bw)

bwt<-valbw_means %>%
  ggplot(aes(y=turb_ysi_m,x=turb_lab, label=site))+
  geom_smooth(method="lm",lty=1,col="dark gray",se=FALSE)+
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+
  ylab("Sensor Turbidity (NTU)")+
  xlab("Sample Turbidity (NTU)")+
  theme_bw()

bwc<-valbw_means %>%
  ggplot(aes(y=chla_ysi_m,x=chla_lab, label=site))+
  geom_smooth(method="lm",lty=1,col="dark gray",se=FALSE)+
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+
  ylab("Sensor Chl-a (RFU)")+
  xlab("Sample Chl-a (mg/L)")+
  theme_bw()

# merging for iv
samplesensoriv_merged<-merge(iv,civ)
samplesensoriv_merged$lat_diff<-samplesensoriv_merged$lat_dec-samplesensoriv_merged$lat
samplesensoriv_merged$lon_diff<-samplesensoriv_merged$lon_dec-samplesensoriv_merged$lon
samplesensoriv_merged<-samplesensoriv_merged %>% filter(abs(lat_diff)<0.0003,abs(lon_diff)<0.0003)

validationdataiv<-samplesensoriv_merged %>% dplyr::select(system, site,lat,lon,turb_lab = turbidity,
                                                          turb_ysi=turb,chla_lab=`chl a`, chla_ysi = chl, chl_ugl) %>% 
  group_by(system, site,lat,lon,turb_lab,turb_ysi,chla_lab, chla_ysi) %>%
  dplyr::summarize(turb_ysi=mean(turb_ysi,na.rm=TRUE)) %>%
  dplyr::summarize(chla_ysi=mean(chla_ysi,na.rm=TRUE)) %>% 
  as.data.frame()

valiv_means <- validationdataiv %>% group_by(site) %>%
  mutate(turb_ysi_m=mean(turb_ysi,na.rm=TRUE)) %>% 
  mutate(chla_ysi_m=mean(chla_ysi,na.rm=TRUE))

##
validationdataivugl<-samplesensoriv_merged %>% dplyr::select(system, site,lat,lon,turb_lab = turbidity,
                                                          turb_ysi=turb,chla_lab=`chl a`, chla_ysi = chl, chl_ugl) %>% 
  group_by(system, site,lat,lon,turb_lab,turb_ysi,chla_lab, chla_ysi) %>%
  dplyr::summarize(turb_ysi=mean(turb_ysi,na.rm=TRUE)) %>%
  dplyr::summarize(chla_ysi=mean(chla_ysi,na.rm=TRUE)) %>% 
  as.data.frame()

valiv_meansugl <- validationdataivugl %>% group_by(site) %>%
  mutate(turb_ysi_m=mean(turb_ysi,na.rm=TRUE)) %>% 
  mutate(chla_ysi_m=mean(chla_ysi,na.rm=TRUE))


lmt_iv<-lm(turb_ysi_m~turb_lab,data=valiv_means) # r2 = 0.9993
summary(lmt_iv)

lmc_iv<-lm(chla_ysi_m~chla_lab,data=valiv_means) # r2 = 0.8014
summary(lmc_iv)

ivt<-valiv_means %>%
  ggplot(aes(y=turb_ysi_m,x=turb_lab, label=site))+
  geom_smooth(method="lm",lty=1,col="dark gray",se=FALSE)+
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+
  ylab("Sensor Turbidity (NTU)")+
  xlab("Sample Turbidity (NTU)")+
  theme_bw()

ivc<-valiv_means %>%
  ggplot(aes(y=chla_ysi_m,x=chla_lab, label=site))+
  geom_smooth(method="lm",lty=1,col="dark gray",se=FALSE)+
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+
  ylab("Sensor Chl-a (RFU)")+
  xlab("Sample Chl-a (mg/L)")+
  theme_bw()

# merging for rb
samplesensorrb_merged<-merge(rb,crb)
samplesensorrb_merged$lat_diff<-samplesensorrb_merged$lat_dec-samplesensorrb_merged$lat
samplesensorrb_merged$lon_diff<-samplesensorrb_merged$lon_dec-samplesensorrb_merged$lon
samplesensorrb_merged<-samplesensorrb_merged %>% filter(abs(lat_diff)<0.0002,abs(lon_diff)<0.0002)

validationdatarb<-samplesensorrb_merged %>% dplyr::select(system, site,lat,lon,turb_lab = turbidity,
                                                          turb_ysi=turb,chla_lab=`chl a`, chla_ysi = chl, chl_ugl) %>% 
  group_by(system, site,lat,lon,turb_lab,turb_ysi,chla_lab, chla_ysi) %>%
  dplyr::summarize(turb_ysi=mean(turb_ysi,na.rm=TRUE)) %>%
  dplyr::summarize(chla_ysi=mean(chla_ysi,na.rm=TRUE)) %>% 
  as.data.frame()

valrb_means <- validationdatarb %>% group_by(site) %>%
  mutate(turb_ysi_m=mean(turb_ysi,na.rm=TRUE)) %>% 
  mutate(chla_ysi_m=mean(chla_ysi,na.rm=TRUE))

##
validationdatarbugl<-samplesensorrb_merged %>% dplyr::select(system, site,lat,lon,turb_lab = turbidity,
                                                          turb_ysi=turb,chla_lab=`chl a`, chla_ysi = chl, chl_ugl) %>% 
  group_by(system, site,lat,lon,turb_lab,turb_ysi,chla_lab, chla_ysi) %>%
  dplyr::summarize(turb_ysi=mean(turb_ysi,na.rm=TRUE)) %>%
  dplyr::summarize(chla_ysi=mean(chla_ysi,na.rm=TRUE)) %>% 
  as.data.frame()

valrb_meansugl <- validationdatarbugl %>% group_by(site) %>%
  mutate(turb_ysi_m=mean(turb_ysi,na.rm=TRUE)) %>% 
  mutate(chla_ysi_m=mean(chla_ysi,na.rm=TRUE))


lmt_rb<-lm(turb_ysi_m~turb_lab,data=valrb_means) # r2 = 0.999 (in nd_rb.R)
summary(lmt_rb)

lmc_rb<-lm(chla_ysi_m~chla_lab,data=valrb_means) # r2 = 0.9032
summary(lmc_rb)

rbt<- valrb_means %>%
  ggplot(aes(y=turb_ysi_m,x=turb_lab, label=site))+
  geom_smooth(method="lm",lty=1,col="dark gray",se=FALSE)+
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+
  ylab("Sensor Turbidity (NTU)")+
  xlab("Sample Turbidity (NTU)")+
  theme_bw()

rbc<-valrb_means %>%
  ggplot(aes(y=chla_ysi_m,x=chla_lab, label=site))+
  geom_smooth(method="lm",lty=1,col="dark gray",se=FALSE)+
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+
  ylab("Sensor Chl-a (RFU)")+
  xlab("Sample Chl-a (mg/L)")+
  theme_bw()

# merging for waco
samplesensorlw_merged<-merge(data_merged,clw)
samplesensorlw_merged$lat_diff<-samplesensorlw_merged$lat_dec-samplesensorlw_merged$lat
samplesensorlw_merged$lon_diff<-samplesensorlw_merged$lon_dec-samplesensorlw_merged$lon
samplesensorlw_merged<-samplesensorlw_merged %>% filter(abs(lat_diff)<0.0003,abs(lon_diff)<0.0003)

validationdatalw<-samplesensorlw_merged %>% dplyr::select(system, site,lat,lon,turb_lab = turbidity,
                                                          turb_ysi=turb,chla_lab=`chl a`, chla_ysi = chl, chl_ugl) %>% 
  group_by(system, site,lat,lon,turb_lab,turb_ysi,chla_lab, chla_ysi) %>%
  dplyr::summarize(turb_ysi=mean(turb_ysi,na.rm=TRUE)) %>%
  dplyr::summarize(chla_ysi=mean(chla_ysi,na.rm=TRUE)) %>% 
  as.data.frame()

vallw_means <- validationdatalw %>% group_by(site) %>%
  mutate(turb_ysi_m=mean(turb_ysi,na.rm=TRUE)) %>% 
  mutate(chla_ysi_m=mean(chla_ysi,na.rm=TRUE))

##
validationdatalwugl<-samplesensorlw_merged %>% dplyr::select(system, site,lat,lon,turb_lab = turbidity,
                                                          turb_ysi=turb,chla_lab=`chl a`, chla_ysi = chl, chl_ugl) %>% 
  group_by(system, site,lat,lon,turb_lab,turb_ysi,chla_lab, chla_ysi) %>%
  dplyr::summarize(turb_ysi=mean(turb_ysi,na.rm=TRUE)) %>%
  dplyr::summarize(chla_ysi=mean(chla_ysi,na.rm=TRUE)) %>% 
  as.data.frame()

vallw_meansugl <- validationdatalwugl %>% group_by(site) %>%
  mutate(turb_ysi_m=mean(turb_ysi,na.rm=TRUE)) %>% 
  mutate(chla_ysi_m=mean(chla_ysi,na.rm=TRUE))

lmt_lw<-lm(turb_ysi_m~turb_lab,data=vallw_means) # r2 = 0.9266
summary(lmt_lw)

lmc_lw<-lm(chla_ysi_m~chla_lab,data=vallw_means) # r2 = 0.8203
summary(lmc_lw)

vallw_means %>%
  ggplot(aes(y=turb_ysi_m,x=turb_lab, label=site))+
  geom_smooth(method="lm",lty=1,col="dark gray",se=FALSE)+
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+
  ylab("Sensor Turbidity (NTU)")+
  xlab("Sample Turbidity (NTU)")+
  theme_bw()



vallw_means %>%
  ggplot(aes(y=chla_ysi_m,x=chla_lab, label=site))+
  geom_smooth(method="lm",lty=1,col="dark gray",se=FALSE)+
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+
  ylab("Sensor Chl-a (RFU)")+
  xlab("Sample Chl-a (mg/L)")+
  theme_bw()

# merging for arrowhead
samplesensorah_merged<-merge(data_merged,cah)
samplesensorah_merged$lat_diff<-samplesensorah_merged$lat_dec-samplesensorah_merged$lat
samplesensorah_merged$lon_diff<-samplesensorah_merged$lon_dec-samplesensorah_merged$lon
samplesensorah_merged<-samplesensorah_merged %>% filter(abs(lat_diff)<0.0003,abs(lon_diff)<0.0003)

validationdataah<-samplesensorah_merged %>% dplyr::select(system, site,lat,lon,turb_lab = turbidity,
                                                          turb_ysi=turb,chla_lab=`chl a`, chla_ysi = chl, chl_ugl) %>% 
  group_by(system,site,lat,lon,turb_lab,turb_ysi,chla_lab, chla_ysi) %>%
  dplyr::summarize(turb_ysi=mean(turb_ysi,na.rm=TRUE)) %>%
  dplyr::summarize(chla_ysi=mean(chla_ysi,na.rm=TRUE)) %>% 
  as.data.frame()

valah_means <- validationdataah %>% group_by(site) %>%
  mutate(turb_ysi_m=mean(turb_ysi,na.rm=TRUE)) %>% 
  mutate(chla_ysi_m=mean(chla_ysi,na.rm=TRUE))

##
validationdataahugl<-samplesensorah_merged %>% dplyr::select(system, site,lat,lon,turb_lab = turbidity,
                                                          turb_ysi=turb,chla_lab=`chl a`, chla_ysi = chl, chl_ugl) %>% 
  group_by(system,site,lat,lon,turb_lab,turb_ysi,chla_lab, chla_ysi) %>%
  dplyr::summarize(turb_ysi=mean(turb_ysi,na.rm=TRUE)) %>%
  dplyr::summarize(chla_ysi=mean(chla_ysi,na.rm=TRUE)) %>% 
  as.data.frame()

valah_meansugl <- validationdataahugl %>% group_by(site) %>%
  mutate(turb_ysi_m=mean(turb_ysi,na.rm=TRUE)) %>% 
  mutate(chla_ysi_m=mean(chla_ysi,na.rm=TRUE))

lmt_ah<-lm(turb_ysi_m~turb_lab,data=valah_means) # r2 = 0.9983
summary(lmt_ah)

lmc_ah<-lm(chla_ysi_m~chla_lab,data=valah_means) # r2 = 0.2543
summary(lmc_ah)

valah_means %>%
  ggplot(aes(y=turb_ysi_m,x=turb_lab, label=site))+
  geom_smooth(method="lm",lty=1,col="dark gray",se=FALSE)+
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+
  ylab("Sensor Turbidity (NTU)")+
  xlab("Sample Turbidity (NTU)")+
  theme_bw()



valah_means %>%
  ggplot(aes(y=chla_ysi_m,x=chla_lab, label=site))+
  geom_smooth(method="lm",lty=1,col="dark gray",se=FALSE)+
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+
  ylab("Sensor Chl-a (RFU)")+
  xlab("Sample Chl-a (mg/L)")+
  theme_bw()

# combining all validation_means_dfs
valmeans_list <- list(valbon2_means,valbon1_means, vallw_means, valbw_means, valiv_means,valrb_means, valah_means)

valmeans_df <- Reduce(function(x, y) merge(x, y, all=TRUE), valmeans_list)

write_csv(valmeans_df, "sensor_sample_valmeans.csv")

# combining all valdfs for ugl
valmeansugl_list <- list(valbon2_meansugl,valbon1_meansugl, vallw_meansugl, valbw_meansugl, valiv_meansugl,valrb_meansugl, valah_meansugl)

valmeansugl_df <- Reduce(function(x, y) merge(x, y, all=TRUE), valmeansugl_list)

write_csv(valmeansugl_df, "senssamp_valmeansugl.csv")

sensamp<-read_csv("sensor_sample_valmeans.csv")
# running lm and plotting r2 for all lab~mean sensor values
lmt_valmeans<-lm(turb_ysi_m~turb_lab,data=sensamp) # r2 = 0.9893
summary(lmt_valmeans)

lm_ugl<-lm(chla_ysi_m~chla_lab,valmeansugl_df)

abb<-valmeansugl_df %>% filter(!system=="Lake Bonham")
abb_rfu<-valmeans_df%>% filter(!system=="Lake Bonham")

vmt_r2 <-ggplot(valmeans_df,aes(y=turb_ysi_m,x=turb_lab)) + 
  geom_point(aes(color = system)) + 
  geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(after_stat(rr.label))), 
               parse = TRUE, label.x.npc = "left", size = 5, rr.digits = 3) + ylab("Sensor Turbidity (NTU)")+
  xlab("Sample Turbidity (NTU)")+ theme_bw()
ggsave("vmt_r2t.png", vmt_r2)

lmc_valmeans<-lm(chla_ysi_m~chla_lab,data=valmeans_df) # r2 = 0.5398
summary(lmc_valmeans)

vmc_r2 <-ggplot(valmeans_df,aes(y=chla_ysi_m,x=chla_lab)) + 
  geom_point(aes(color = system)) + 
  geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(after_stat(rr.label))), 
               parse = TRUE, label.x.npc = "left", size = 5, rr.digits = 3) + ylab("Sensor Chl-a (RFU)")+
  xlab("Sample Chl-a (ug/L)")+ theme_bw()
ggsave("vm_r2c.png", vmc_r2)

# r2 plots
bnt_r2 <-ggplot(valbon_means,aes(y=turb_ysi_m,x=turb_lab, label=site)) + 
  geom_point(shape=21,size=8,fill="white") + 
  geom_text()+
  geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(after_stat(rr.label))), 
               parse = TRUE, label.x.npc = "left", size = 5, rr.digits = 3) + ylab("Sensor Turbidity (NTU)")+
  xlab("Sample Turbidity (NTU)")+ theme_bw()
ggsave("bon_r2t.png", bnt_r2)

bnc_r2<-ggplot(valbon_means,aes(y=chla_ysi_m,x=chla_lab, label=site)) + 
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+ 
  geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(after_stat(rr.label))), 
               parse = TRUE, label.x.npc = "left", size = 5, rr.digits = 3) + ylab("Sensor Chl-a (RFU)")+
  xlab("Sample Chl-a (mg/L)")+ theme_bw()
ggsave("bon_r2c.png", bnc_r2)

bwt_r2<-ggplot(valbw_means,aes(y=turb_ysi_m,x=turb_lab, label=site)) + 
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+ 
  geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(after_stat(rr.label))), 
               parse = TRUE, label.x.npc = "left", size = 5, rr.digits = 3) +ylab("Sensor Turbidity (NTU)")+
  xlab("Sample Turbidity (NTU)")+ theme_bw()
ggsave("bw_r2t.png", bwt_r2)

bwc_r2<-ggplot(valbw_means,aes(y=chla_ysi_m,x=chla_lab, label=site)) + 
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+ 
  geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(after_stat(rr.label))), 
               parse = TRUE, label.x.npc = "left", size = 5, rr.digits = 3) + ylab("Sensor Chl-a (RFU)")+
  xlab("Sample Chl-a (mg/L)")+ theme_bw()
ggsave("bw_r2c.png", bwc_r2)

ivt_r2<-ggplot(valiv_means,aes(y=turb_ysi_m,x=turb_lab, label=site)) + 
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+ 
  geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(after_stat(rr.label))), 
               parse = TRUE, label.x.npc = "left", size = 5, rr.digits = 3) + ylab("Sensor Turbidity (NTU)")+
  xlab("Sample Turbidity (NTU)")+theme_bw()
ggsave("iv_r2t.png", ivt_r2)

ivc_r2<-ggplot(valiv_means,aes(y=chla_ysi_m,x=chla_lab, label=site)) + 
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+
  geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(after_stat(rr.label))), 
               parse = TRUE, label.x.npc = "left", size = 5, rr.digits = 3) + ylab("Sensor Chl-a (RFU)")+
  xlab("Sample Chl-a (mg/L)")+ theme_bw()
ggsave("iv_r2c.png", ivc_r2)


rbt_r2<-ggplot(valrb_means,aes(y=turb_ysi_m,x=turb_lab, label=site)) + 
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+
  geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(after_stat(rr.label))), 
               parse = TRUE, label.x.npc = "left", size = 5, rr.digits = 3) + ylab("Sensor Turbidity (NTU)")+
  xlab("Sample Turbidity (NTU)")+theme_bw()
ggsave("rb_r2t.png", rbt_r2)

rbc_r2<-ggplot(valrb_means,aes(y=chla_ysi_m,x=chla_lab, label=site)) + 
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+ 
  geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(after_stat(rr.label))), 
               parse = TRUE, label.x.npc = "left", size = 5, rr.digits = 3) + ylab("Sensor Chl-a (RFU)")+
  xlab("Sample Chl-a (mg/L)")+ theme_bw()
ggsave("rb_r2c.png", rbc_r2)

lwt_r2<-ggplot(vallw_means,aes(y=turb_ysi_m,x=turb_lab, label=site)) + 
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+ 
  geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(after_stat(rr.label))), 
               parse = TRUE, label.x.npc = "left", size = 5, rr.digits = 3)  + ylab("Sensor Turbidity (NTU)")+
  xlab("Sample Turbidity (NTU)")+theme_bw()
ggsave("lw_r2t.png", lwt_r2)

lwc_r2<-ggplot(vallw_means,aes(y=chla_ysi_m,x=chla_lab, label=site)) + 
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+ 
  geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(after_stat(rr.label))), 
               parse = TRUE, label.x.npc = "left", size = 5, rr.digits = 3) + ylab("Sensor Chl-a (RFU)")+
  xlab("Sample Chl-a (mg/L)")+ theme_bw()
ggsave("lw_r2c.png", lwc_r2)

aht_r2<-ggplot(valah_means,aes(y=turb_ysi_m,x=turb_lab, label=site)) + 
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+ 
  geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(after_stat(rr.label))), 
               parse = TRUE, label.x.npc = "left", size = 5, rr.digits = 3) + 
  ylab("Sensor Turbidity (NTU)")+
  xlab("Sample Turbidity (NTU)")+theme_bw()
ggsave("ah_r2t.png", aht_r2)

ahc_r2<-ggplot(valah_means,aes(y=chla_ysi_m,x=chla_lab, label=site)) + 
  geom_point(shape=21,size=8,fill="white")+
  geom_text()+ 
  geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(after_stat(rr.label))), 
               parse = TRUE, label.x.npc = "left", size = 5, rr.digits = 3) + ylab("Sensor Chl-a (RFU)")+
  xlab("Sample Chl-a (mg/L)")+ theme_bw()
ggsave("ah_r2c.png", ahc_r2)


#5. Combine plots

combturb_plot<-grid.arrange(ncol=4,
                            bnt,bwt,
                            ivt,rbt,
                            padding = unit(0.0, "line"))


combchl_plot<-grid.arrange(ncol=4,
                            bnc,bwc,
                            ivc,rbc,
                            padding = unit(0.0, "line"))


turb_all6 <- grid.arrange(ncol=3,
                          bnt_r2,bwt_r2,
                          ivt_r2,rbt_r2,
                          lwt_r2, aht_r2,
                          padding = unit(0.0, "line"))
ggsave(filename = "combo_turb_lm.png", plot = turb_all6, width = 10, height = 5)

chla_all6<-grid.arrange(ncol=3,
                           bnc_r2,bwc_r2,
                           ivc_r2,rbc_r2,
                        lwc_r2, ahc_r2,
                           padding = unit(0.0, "line"))
ggsave(filename = "combo_chla_lm.png", plot = chla_all6, width = 10, height = 5)
