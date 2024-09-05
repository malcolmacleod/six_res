# A script to develop regression model between SDD and turb

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
#library(rgdal)
library(grid)
library(gridExtra)
library(lattice)
library(rioja)
library(zoo)
library(ggpubr)
library(OpenStreetMap)
library(basemaps)
library(RColorBrewer)
library(ggspatial)
library(ggpmisc)
library(cowplot)
library(gridExtra)
library(sf)

# 2. CRASR data - read in for all lakes -------------------------------------------
crasr_all <- read_csv("CRASR_summer2022.csv") %>% clean_names()

# separate by lake day
crasr_bon1 <- crasr_all %>% filter(datecollect == "7/14/2022") %>% add_column(lake = "bon") %>% group_by(system,lat, lon, variable, site) %>% mutate(meanchl = mean(value))
cbon1 <- crasr_bon1 %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl) %>% dplyr::select(lat, lon, system, site, chla = `chl a`, turbidity, tss, afdm, tp,tn, secchi)

crasr_bon2 <- crasr_all %>% filter(datecollect == "7/15/2022") %>% add_column(lake = "bon") %>% 
  dplyr::mutate(site = paste(lake, site, sep = "_")) %>% group_by(system,lat, lon, variable, site) %>% mutate(meanchl = mean(value))
cbon2 <- crasr_bon2 %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl)%>% dplyr::select(lat, lon, site, chla = `chl a`, turbidity, tss, afdm, tp,tn, secchi)

crasr_lw <- crasr_all %>% filter(datecollect == "7/23/2022")%>% add_column(lake = "lw") %>% 
  dplyr::mutate(site = paste("LW", site, sep = "")) %>% group_by(system,lat, lon, variable, site) %>% mutate(meanchl = mean(value))
clw <- crasr_lw %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl)%>% dplyr::select(lat, lon, system, site, chla = `chl a`, turbidity, tss, afdm, tp,tn)


crasr_bw <- crasr_all %>% filter(datecollect == "8/6/2022")%>% add_column(lake = "bw") %>% 
  dplyr::mutate(station = paste(lake, site, sep = "")) %>% group_by(system,lat, lon, variable, site) %>% mutate(meanchl = mean(value))
cbw <- crasr_bw %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl) %>% dplyr::select(lat, lon, system, site, chla = `chl a`, turbidity, tss, afdm, tp,tn, secchi)

crasr_iv <- crasr_all %>% filter(datecollect == "8/7/2022") %>% add_column(lake = "iv") %>% 
  dplyr::mutate(station = paste(lake, site, sep = "")) %>% group_by(system,lat, lon, variable, site) %>% mutate(meanchl = mean(value))
civ <- crasr_iv %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl) %>% dplyr::select(lat, lon, system, site, chla = `chl a`, turbidity, tp,tn, secchi)

crasr_rb <- crasr_all %>% filter(datecollect == "8/8/2022") %>% add_column(lake = "rb") %>% 
  dplyr::mutate(station = paste(lake, site, sep = "")) %>% group_by(system,lat, lon, variable, site) %>% mutate(meanchl = mean(value))
crb <- crasr_rb %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl) %>% dplyr::select(lat, lon, system, site, chla = `chl a`, turbidity, tss, afdm, tp,tn, secchi)

crasr_ah <- crasr_all %>% filter(datecollect == "8/12/2022") %>% add_column(lake = "ah") %>% 
  dplyr::mutate(station = paste(lake, site, sep = "")) %>% group_by(system, lat, lon, variable, site) %>% mutate(meanchl = mean(value))
cah <- crasr_ah %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl) %>% dplyr::select(lat, lon, system, site, chla = `chl a`, turbidity, tss, afdm, tp,tn, secchi)

# combine crasr data
crasr_list <- list(cbon1, cbon2, clw, civ, cbw, crb, cah)

#merge all data frames in list
crasr_df <- Reduce(function(x, y) merge(x, y, all=TRUE), crasr_list)
#write_csv(crasr_df, "crasr_df.csv")

# 3. trying lm on just the data from sample lakes -------------------------------
# reducing just stations, turb, and sdd for lm
turb_sdd <- crasr_df %>% dplyr::select(site, turbidity, secchi)

lmt_sdd<-lm(log(turbidity)~secchi,data=turb_sdd)  # log(turb) = -1.6149(x) + 3.4832 
summary(lmt_sdd) 
plot(lmt_sdd)

exp(-1.69149) # 0.1842448
exp(coef(lmt_sdd)["secchi"]) # 0.1989103

sdd<-turb_sdd$secchi
turb<-turb_sdd$turbidity

ggplot(turb_sdd,aes(y=log(turbidity),x=secchi)) + 
  geom_point()  +   geom_smooth(method = "lm", se=FALSE) +
    stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()


# 4. Bringing in NLA data -------------------------------------------------------
#nla_all <- read_csv("NLA_tsdd_all.csv") %>% rename(secchi = `secchi (m)`, turb = `turb (ntu)`)
#nla_tx <- read_csv("NLA_tsdd_tx.csv")%>% rename(secchi = `secchi (m)`, turb = `turb (ntu)`)
nla_res<-read_csv("NLA_tsdd_res.csv") %>% dplyr::rename(secchi = `secchi (m)`, turb = `turb (ntu)`)

# combining NLA data (all reservoirs) with our own
nla_crasr<-read_csv("nla_crasr.csv")

nla_crasr <- nla_crasr%>% dplyr::rename(secchi = `secchi (m)`, turb = `turb (ntu)`)

tsdlm<-lm(log(secchi)~log(turb), data = nla_crasr)
summary(tsdlm) # exp(y) = 1.13734 -0.68001(exp(x))

nla_crasr$residuals <- tsdlm$residuals
nla_crasr$predicted<-tsdlm$fitted.values
  
tsd_plot <- ggplot(nla_crasr,aes(log(turb), log(secchi))) + 
  geom_point() + 
  geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(after_stat(rr.label))), 
               parse = TRUE, label.x.npc = "left", size = 8, rr.digits = 3)+
  ggtitle(label = "Turbidity ~ Secchi for NLA Reservoirs 2012 and 2017")+
  xlab("log(Turbidity)") + ylab("log(Secchi Depth)") +
  theme_bw() +theme(plot.title = element_text(size = 20,hjust=0.5),
                    axis.title.x = element_text(size = 15),
                    axis.title.y = element_text(size = 15),
                    legend.title = element_text(size = 15))

# exp equation to predict SDD from turb based off lm is: secchi=3.118⋅(turb)^ −0.68001

#ggsave("nla_crasr_lm.png", tsd_plot)

# 5. a la arrowhead paper ---------------------------------------------------------
#x<-nla_crasr$turb
#y<-nla_crasr$secchi
#decay<-nls(y ~ 1/(b + x^c), start=list(b=1,c=1))
#xvec <- seq(min(nla_crasr$turb),max(nla_crasr$turb),length=1000)
#preds <- predict(decay,newdata=data.frame(x=xvec))
#pred<- data.frame(x = xvec, y = preds)

# nla_crasr %>% ggplot(aes(y=secchi,x=turb))+
# geom_line(data = pred, aes(x=x, y=y),lty=1,col="dark gray",lwd=1)+ geom_point()+
# ylab("Secchi Depth (m)") +xlab("Sample Turbidity (NTU)")+theme_bw()


# 6a Read in flame data files ------------------------------------------

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

garmin<-garmin %>% dplyr::select(garmin_datetime_utc,date,hour_dec,garmin_time_local,xcoord,ycoord, speed)


# 6b Merge data -----------------------------------------------------------

data_merged<-merge(exo,garmin, 
                   by.x=c("lake_time_local","date"), 
                   by.y=c("garmin_time_local", "date"))


data_merged<-data_merged %>% dplyr::select(lat_dec=ycoord,lon_dec=xcoord,speed,
                                           time_local=lake_time_local,date,datetime_utc=garmin_datetime_utc,
                                           wtemp=temp_c, chl=chlorophyll_rfu,pc=tal_pc_rfu,fdom=f_dom_rfu,
                                           turb=turbidity_fnu,do_mgl=odo_mg_l, do_sat=odo_percent, hour_dec
)

data_merged$date_time_c = paste(data_merged$date, data_merged$time_local)
# averaging and grouping by minute
#data_merged$hour <- hour(data_merged$date_time_c)
#data_merged$minute <- minute(data_merged$date_time_c)
#data_merged$hm <- paste(data_merged$hour, data_merged$minute, sep = "_")

# filtering garmin data by lake / date
bon1<-data_merged %>% filter(date=="2022-07-14", hour_dec>13.20, hour_dec<16.733) %>% filter(speed>0)
bon2<-data_merged %>% filter(date=="2022-07-15", hour_dec>12, hour_dec<13.91)%>% filter(speed>0)
waco<-data_merged %>% filter(date=="2022-07-23", hour_dec>10.767, hour_dec<13.25)%>% filter(speed>0)
bw<-data_merged %>% filter(date=="2022-08-06", hour_dec>12.91, hour_dec<15.283)%>% filter(speed>0)
iv<-data_merged %>% filter(date=="2022-08-07", hour_dec>10.67, hour_dec<13.68)%>% filter(speed>0)
rb<-data_merged %>% filter(date=="2022-08-08", hour_dec>10.7, hour_dec<15.967)%>% filter(speed>0)
ah<-data_merged %>% filter(date=="2022-08-12", hour_dec>9.8, hour_dec<12.8)%>% filter(speed>0)

# 6c predicting secchi from turb values ------------------------------------------------
ah$pred_sdd <- exp(predict(tsdlm, ah))
bon1$pred_sdd <- exp(predict(tsdlm, bon1))
bon2$pred_sdd <- exp(predict(tsdlm, bon2))
waco$pred_sdd <- exp(predict(tsdlm, waco))
bw$pred_sdd <- exp(predict(tsdlm, bw))
iv$pred_sdd <- exp(predict(tsdlm, iv))
rb$pred_sdd <- exp(predict(tsdlm, rb))

# calculating TSI from secchi
ah$tsi_sd <- 10 * (6-log(ah$pred_sdd)/log(2))
bon1$tsi_sd <- 10 * (6-log(bon1$pred_sdd)/log(2))
bon2$tsi_sd <- 10 * (6-log(bon2$pred_sdd)/log(2))
bw$tsi_sd <- 10 * (6-log(bw$pred_sdd)/log(2))
iv$tsi_sd <- 10 * (6-log(iv$pred_sdd)/log(2))
rb$tsi_sd <- 10 * (6-log(rb$pred_sdd)/log(2))
waco$tsi_sd <- 10 * (6-log(waco$pred_sdd)/log(2))

# also adding secchi from Sherjah paper
ah$sherjah<-(244.13 * ah$turb^-0.662) /100
bon2$sherjah<-(244.13 * bon2$turb^-0.662) /100
bw$sherjah<-(244.13 * bw$turb^-0.662) /100
iv$sherjah<-(244.13 * iv$turb^-0.662) /100
rb$sherjah<-(244.13 * rb$turb^-0.662) /100
waco$sherjah<-(244.13 * waco$turb^-0.662) /100
# combine crasr data
sdd_list <- list(bon2, waco, bw, iv,rb, ah)

#merge all data frames in list
sdd_df <- Reduce(function(x, y) merge(x, y, all=TRUE), sdd_list)

write_csv(sdd_df, "flm_sdd.csv")

# -------------------------------------------------------------------
# attempt to viz sdd to pred sdd
sampsensbon <- merge(bon2, cbon2)
sampsensbon$lat_diff<-sampsensbon$lat_dec-sampsensbon$lat
sampsensbon$lon_diff<-sampsensbon$lon_dec-sampsensbon$lon
ssbn_merged<-sampsensbon %>% filter(abs(lat_diff)<0.0004,abs(lon_diff)<0.0004)

valdatabon<-ssbn_merged %>% dplyr::select(system, site,lat,lon,
                                          turb_lab = turbidity,turb_ysi=turb,
                                          secchi,pred_sdd,afdm, tss,tn,tp, sherjah) 

valbon_means <- valdatabon %>% group_by(site) %>%
  mutate(turb_ysi_m=mean(turb_ysi,na.rm=TRUE)) %>% 
  mutate(pred_sdd_m=mean(pred_sdd,na.rm=TRUE)) %>% 
  mutate(sdd_s_m=mean(sherjah,na.rm=TRUE))

# with brownwood
sampsensbw <- merge(bw, cbw)
sampsensbw$lat_diff<-sampsensbw$lat_dec-sampsensbw$lat
sampsensbw$lon_diff<-sampsensbw$lon_dec-sampsensbw$lon
ssbw_merged<-sampsensbw %>% filter(abs(lat_diff)<0.0004,abs(lon_diff)<0.0004)

valdatabw<-ssbw_merged %>% dplyr::select(system, site,lat,lon,
                                         turb_lab = turbidity,turb_ysi=turb,
                                         secchi,pred_sdd,afdm, tss,tn,tp,sherjah) 

valbw_means <- valdatabw %>% group_by(site) %>%
  mutate(turb_ysi_m=mean(turb_ysi,na.rm=TRUE)) %>% 
  mutate(pred_sdd_m=mean(pred_sdd,na.rm=TRUE))%>% 
  mutate(sdd_s_m=mean(sherjah,na.rm=TRUE))

# now with iv
sampsensiv <- merge(iv, civ)
sampsensiv$lat_diff<-sampsensiv$lat_dec-sampsensiv$lat
sampsensiv$lon_diff<-sampsensiv$lon_dec-sampsensiv$lon
ssiv_merged<-sampsensiv %>% filter(abs(lat_diff)<0.0004,abs(lon_diff)<0.0004)

valdataiv<-ssiv_merged %>% dplyr::select(system, site,lat,lon,
                                         turb_lab = turbidity,turb_ysi=turb,
                                         secchi,pred_sdd,tn,tp,sherjah) 

valiv_means <- valdataiv %>% group_by(site) %>%
  mutate(turb_ysi_m=mean(turb_ysi,na.rm=TRUE)) %>% 
  mutate(pred_sdd_m=mean(pred_sdd,na.rm=TRUE))%>% 
  mutate(sdd_s_m=mean(sherjah,na.rm=TRUE))


# now with rb
sampsensrb <- merge(rb, crb)
sampsensrb$lat_diff<-sampsensrb$lat_dec-sampsensrb$lat
sampsensrb$lon_diff<-sampsensrb$lon_dec-sampsensrb$lon
ssrb_merged<-sampsensrb %>% filter(abs(lat_diff)<0.0004,abs(lon_diff)<0.0004)

valdatarb<-ssrb_merged %>% dplyr::select(system, site,lat,lon,
                                         turb_lab = turbidity,turb_ysi=turb,
                                         secchi,pred_sdd,afdm, tss,tn,tp,sherjah) 

valrb_means <- valdatarb %>% group_by(site) %>%
  mutate(turb_ysi_m=mean(turb_ysi,na.rm=TRUE)) %>% 
  mutate(pred_sdd_m=mean(pred_sdd,na.rm=TRUE))%>% 
  mutate(sdd_s_m=mean(sherjah,na.rm=TRUE))

# now with ah
sampsensah <- merge(ah, cah)
sampsensah$lat_diff<-sampsensah$lat_dec-sampsensah$lat
sampsensah$lon_diff<-sampsensah$lon_dec-sampsensah$lon
ssah_merged<-sampsensah %>% filter(abs(lat_diff)<0.0004,abs(lon_diff)<0.0004)

valdataah<-ssah_merged %>% dplyr::select(system, site,lat,lon,
                                         turb_lab = turbidity,turb_ysi=turb,
                                         secchi,pred_sdd,afdm, tss,tn,tp,sherjah) 

valah_means <- valdataah %>% group_by(site) %>%
  mutate(turb_ysi_m=mean(turb_ysi,na.rm=TRUE)) %>% 
  mutate(pred_sdd_m=mean(pred_sdd,na.rm=TRUE))%>% 
  mutate(sdd_s_m=mean(sherjah,na.rm=TRUE))

# now with waco
sampsenslw <- merge(waco, clw)
sampsenslw$lat_diff<-sampsenslw$lat_dec-sampsenslw$lat
sampsenslw$lon_diff<-sampsenslw$lon_dec-sampsenslw$lon
sslw_merged<-sampsenslw %>% filter(abs(lat_diff)<0.0004,abs(lon_diff)<0.0004)

valdatalw<-sslw_merged %>% dplyr::select(system, site,lat,lon,
                                         turb_lab = turbidity,turb_ysi=turb,
                                         pred_sdd,afdm, tss,tn,tp,sherjah) 

vallw_means <- valdatalw %>% group_by(site) %>%
  mutate(turb_ysi_m=mean(turb_ysi,na.rm=TRUE)) %>% 
  mutate(pred_sdd_m=mean(pred_sdd,na.rm=TRUE))%>% 
  mutate(sdd_s_m=mean(sherjah,na.rm=TRUE))

# combining all validation_means_dfs
valmeans_list <- list(valbon_means, valah_means, valbw_means, valiv_means,valrb_means, valah_means, vallw_means)

valmeans_df <- Reduce(function(x, y) merge(x, y, all=TRUE), valmeans_list)

# plotting out both turbidity and predicted SDD to stations

turbplot<-ggplot(valmeans_df,aes(y=turb_ysi_m,x=turb_lab)) + 
  geom_point(aes(color = system)) + 
  geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(after_stat(rr.label))), 
               parse = TRUE, label.x.npc = "left", size = 5, rr.digits = 3) + ylab("Sensor Turbidity (NTU)")+
  xlab("Sample Turbidity (NTU)")+ theme_bw()
ggsave("samp_ysi_turb.png", turbplot)

predsdd_plot<-valmeans_df %>% filter(!system=="Lake Waco") %>% 
  ggplot(aes(y=pred_sdd_m,x=secchi)) + 
  geom_point(aes(color = system)) + 
  geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(after_stat(rr.label))), 
               parse = TRUE, label.x.npc = "left", size = 5, rr.digits = 3) + ylab("Predicted SDD (m)")+
  xlab("Secchi (m)")+ theme_bw()
ggsave("sdd_predsdd.png", predsdd_plot)

#sdd_turb_ab<-plot_grid(turbplot, predsdd_plot, align = 'v', labels = c('A', 'B'))
#----------------------------------------------------------------------
# removing outlier value in iv and waco
iv <- iv %>% filter(iv$turb<25) %>% filter(date_time_c<"2022-08-07 13:28")
waco <- waco %>% filter(waco$turb<75)
# 7. visualizing secchi and TSI with maps-------------------------------------------------------
#dm_sf <- st_as_sf(iv, coords = c("lon_dec", "lat_dec"), crs = 4326)
#mapview(dm_sf, zcol= "tsi_sd")
# set some parameters to pass into ggplot 
ptsize<-0.25
ptalpha<-0.5

# setting gps boundaries for each lake-map
map_ah <- openmap(upperLeft = c(33.78, -98.462), 
                lowerRight = c(33.615, -98.312), 
                type = 'esri-imagery',zoom=14)

map_bon <- openmap(upperLeft = c(33.675,-96.175),
                  lowerRight = c(33.625,-96.125), 
                  type = 'esri-imagery',zoom=14)

map_bw <- openmap(upperLeft = c(31.93, -99.11),
                  lowerRight = c(31.74, -98.985), 
                  type = 'esri-imagery',zoom=14)

map_iv <- openmap(upperLeft = c(31.62, -99.84),
                  lowerRight = c(31.4, -99.6), 
                  type = 'esri-imagery',zoom=14)

map_rb <- openmap(upperLeft = c(32.026, -104.03),
                  lowerRight = c(31.84, -103.879), 
                  type = 'esri-imagery',zoom=14)

map_lw <- openmap(upperLeft = c(31.625,-97.3),
                lowerRight = c(31.47,-97.19),
                type = 'esri-imagery',zoom=14)


my_palette <- colorRampPalette(rev(brewer.pal(9, "YlOrBr")))

# plotting secchi
plot_sdah<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_ah)) +
  geom_jitter(data = ah, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon_dec, y = lat_dec,color=pred_sdd))+
  scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)")+ 
  ggtitle(label = "Arrowhead") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top"))+ 
  theme(plot.title = element_text(hjust = 0.5))

plot_sdbn1<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_bon)) +
  geom_jitter(data = bon1, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon_dec, y = lat_dec,color=pred_sdd))+
  scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)")+ 
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top"))

plot_sdbn2<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_bon)) +
  geom_jitter(data = bon2, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon_dec, y = lat_dec,color=pred_sdd))+
  ggtitle(label = "Bonham") +
  scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)")+ 
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top"))+ 
  theme(plot.title = element_text(hjust = 0.5))

plot_sdbw<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_bw)) +
  geom_jitter(data = bw, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon_dec, y = lat_dec,color=pred_sdd))+
  scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)")+ 
  ggtitle(label = "Brownwood") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top"))+ 
  theme(plot.title = element_text(hjust = 0.5))

plot_sdiv<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_iv)) +
  geom_jitter(data = iv, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon_dec, y = lat_dec,color=pred_sdd))+
  scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)")+ 
  ggtitle(label = "O.H. Ivie") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top"))+ 
  theme(plot.title = element_text(hjust = 0.5))

plot_sdrb<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_rb)) +
  geom_jitter(data = rb, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon_dec, y = lat_dec,color=pred_sdd))+
  scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)")+ 
  ggtitle(label = "Red Bluff") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top"))+ 
  theme(plot.title = element_text(hjust = 0.5))

plot_sdlw<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_lw)) +
  geom_jitter(data = waco, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon_dec, y = lat_dec,color=pred_sdd))+
  scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)")+ 
  ggtitle(label = "Waco") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top")) + 
  theme(plot.title = element_text(hjust = 0.5))

# Combine maps and save as external .png image

combined_plot<-grid.arrange(ncol=3,
                            plot_sdbn2,plot_sdlw,
                            plot_sdbw,
                            plot_sdiv,plot_sdrb,plot_sdah)+ annotation_scale(location = "bl", width_hint=0.4)

ggsave(filename = "combo_secchi_flm_imagery.png", plot = combined_plot, width = 10, height = 5)

################################################################################################################
# attempts to have a common scale bar and similar area covered in bounding box

# Define the bounding box size (latitude and longitude degrees)
bbox_width <- 0.25
bbox_height <- 0.25
# Define the center coordinates for each map
centers <- list(
  c(-96.15, 33.65), # Center for bon
  c(-97.22, 31.57), # Center for waco
  c(-99.03, 31.83),  # Center for bw
  c(-99.68, 31.55),  # Center for iv
  c(-103.94, 31.95),  # Center for rb
  c(-98.37,33.71)  # Center for ah
)


create_map <- function(center_lon, center_lat, bbox_width, bbox_height) {
  # Define the bounding box
  upperLeft <- c(center_lat + bbox_height / 2, center_lon - bbox_width / 2)
  lowerRight <- c(center_lat - bbox_height / 2, center_lon + bbox_width / 2)
  
  # Fetch the map
  map <- openmap(upperLeft, lowerRight, type = "esri-imagery")
  
  # Convert to ggplot object
  autoplot(map) + theme_minimal()
}

# create maps
map1 <- create_map(centers[[1]][1], centers[[1]][2], bbox_width, bbox_height) + ggtitle("Map 1")
map2 <- create_map(centers[[2]][1], centers[[2]][2], bbox_width, bbox_height) + ggtitle("Map 2")
map3 <- create_map(centers[[3]][1], centers[[3]][2], bbox_width, bbox_height) + ggtitle("Map 3")
map4 <- create_map(centers[[4]][1], centers[[4]][2], bbox_width, bbox_height) + ggtitle("Map 4")
map5 <- create_map(centers[[5]][1], centers[[5]][2], bbox_width, bbox_height) + ggtitle("Map 5")
map6 <- create_map(centers[[6]][1], centers[[6]][2], bbox_width, bbox_height) + ggtitle("Map 6")

# Combine plots into a single layout with grid.arrange
grid.arrange(map1, map2, map3, map4, map5, map6, ncol = 3)

# Optionally, add a common scale bar using cowplot
scale_bar <- ggplot() + 
  theme_void() + 
  theme(legend.position = "bottom") +
  guides(fill = guide_colourbar(barwidth = 10, barheight = 1))

combined_plot <- plot_grid(
  plot_grid(map1, map2, map3, map4, map5, map6, ncol = 3, align = 'v'),
  scale_bar,
  ncol = 1,
  rel_heights = c(1, 0.1)
)

print(combined_plot)

################################################################################################

flm_sdd<-read_csv("flm_sdd.csv")


# set some parameters to pass into ggplot 
ptsize<-0.25
ptalpha<-0.5


# attempts to have a common scale bar and similar area covered in bounding box

# Define the bounding box size (latitude and longitude degrees)
bbox_width <- 0.25
bbox_height <- 0.25
# Define the center coordinates for each map
centers <- list(
  c(-96.15, 33.65), # Center for bon
  c(-97.22, 31.57), # Center for waco
  c(-99.03, 31.83),  # Center for bw
  c(-99.68, 31.55),  # Center for iv
  c(-103.94, 31.95),  # Center for rb
  c(-98.37,33.71)  # Center for ah
)


df_bon <- data.frame(lon = c(-96.1, -96.2), lat = c(33.6, 33.7))
df_waco <- data.frame(lon = c(-97.2, -97.3), lat = c(31.5, 31.6))
df_bw <- data.frame(lon = c(-99.0, -99.1), lat = c(31.8, 31.9))
df_iv <- data.frame(lon = c(-99.6, -99.7), lat = c(31.5, 31.6))
df_rb <- data.frame(lon = c(-103.9, -104.0), lat = c(31.9, 32.0))
df_ah <- data.frame(lon = c(-98.3, -98.4), lat = c(33.7, 33.8))

create_plot <- function(center_lon, center_lat, bbox_width, bbox_height, df, title, add_scalebar = FALSE) {
  # Define the bounding box
  upperLeft <- c(center_lat + bbox_height / 2, center_lon - bbox_width / 2)
  lowerRight <- c(center_lat - bbox_height / 2, center_lon + bbox_width / 2)
  
  # Fetch the map
  map <- OpenStreetMap::openmap(upperLeft, lowerRight, type = "esri_imagery")
  map_projected <- OpenStreetMap::openproj(map)
  
  # Create the plot
  plot <- OpenStreetMap::autoplot.OpenStreetMap(map_projected) +
    geom_jitter(data = df, 
                size = psize,
                alpha = ptalpha,
                aes(x = lon, y = lat, color = pred_sdd)) +
    scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)") +
    ggtitle(label = title) +
    xlab("") + ylab("") +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.direction = "horizontal",
      legend.title.align = 0.5,
      legend.key.width = unit(0.75, "cm"),
      legend.position = c(0.5, 0.03)) +
    guides(col = guide_colorbar(title.position = "top")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  # Add scale bar if specified
  if (add_scalebar) {
    plot <- plot +
      scalebar(
        dist = 1,  # Length of the scale bar in kilometers
        dist_unit = "km",
        transform = TRUE,
        model = "WGS84",
        location = "bottomleft",
        units = "km",
        scale = 0.5,  # Adjust size of scale bar
        st.size = 3   # Adjust size of scale bar text
      )
  }
  
  return(plot)
}
# create maps
plot_bon <- create_plot(centers[[1]][1], centers[[1]][2], bbox_width, bbox_height, df_bon, "Bonham", add_scalebar = TRUE)
plot_waco <- create_plot(centers[[2]][1], centers[[2]][2], bbox_width, bbox_height, waco, "Waco")
plot_bw <- create_plot(centers[[3]][1], centers[[3]][2], bbox_width, bbox_height, bw, "Brownwood")
plot_iv <- create_plot(centers[[4]][1], centers[[4]][2], bbox_width, bbox_height, iv, "OH Ivie")
plot_rb <- create_plot(centers[[5]][1], centers[[5]][2], bbox_width, bbox_height, rb, "Red Bluff")
plot_ah <- create_plot(centers[[6]][1], centers[[6]][2], bbox_width, bbox_height, ah, "Arrowhead",add_scalebar = TRUE)
#####################################################################################################

# set some parameters to pass into ggplot 
ptsize<-0.25
ptalpha<-0.5

# attempts to have a common scale bar and similar area covered in bounding box

# Define the bounding box size (latitude and longitude degrees)
bbox_width_zoomed_in <- 0.05
bbox_height_zoomed_in <- 0.05
bbox_width <- 0.15
bbox_height <- 0.15
bbox_width_iv <- 0.25
bbox_height_iv <- 0.25

# trying to add scale by converting to sf
sf_bon2<-st_as_sf(bon2, coords = c("lon_dec", "lat_dec"), crs = 4326)
sf_waco<-st_as_sf(waco, coords = c("lon_dec", "lat_dec"), crs = 4326)
sf_bw<-st_as_sf(bw, coords = c("lon_dec", "lat_dec"), crs = 4326)
sf_iv<-st_as_sf(iv, coords = c("lon_dec", "lat_dec"), crs = 4326)
sf_rb<-st_as_sf(rb, coords = c("lon_dec", "lat_dec"), crs = 4326)
sf_ah<-st_as_sf(ah, coords = c("lon_dec", "lat_dec"), crs = 4326)
sf_bon2<-st_as_sf(bon2, coords = c("lon_dec", "lat_dec"), crs = 4326)

create_plot <- function(center_lon, center_lat, bbox_width, bbox_height, df, title, add_scalebar=FALSE) {
  # Define the bounding box
  upperLeft <- c(center_lat + bbox_height / 2, center_lon - bbox_width / 2)
  lowerRight <- c(center_lat - bbox_height / 2, center_lon + bbox_width / 2)
  
  # Fetch the map
  map <- openmap(upperLeft, lowerRight, type = "bing")
  map_projected <- OpenStreetMap::openproj(map)
  
  # Create the plot
  plot <- OpenStreetMap::autoplot.OpenStreetMap(map_projected) +
    geom_jitter(data = df, 
            size = ptsize,
            alpha = ptalpha,
            aes(x = lon_dec, y = lat_dec, color = pred_sdd)) +
    scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)") +
    ggtitle(label = title) +
    xlab("") + ylab("") +
    theme(
      plot.margin = unit(c(-1, 0, -1, 0), "lines"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.direction = "horizontal",
      legend.title.align = 0.5,
      legend.key.width = unit(0.75, "cm"),
      legend.position = c(0.5, 0.02),
      legend.text = element_text(size = 8)
    ) +
    annotation_scale(
      location = "tr",   # Start with top-right location
      width_hint = 0.5,
      unit = "km",       # Set unit to kilometers
      pad_x = unit(-3, "cm"),  # Adjust padding to move the scale towards the center
      pad_y = unit(-0.5, "cm") # Adjust padding to move the scale slightly down
    ) +
    guides(col = guide_colorbar(title.position = "top")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(plot)
}


# Define the center coordinates for each map
centers <- list(
  c(-96.15, 33.65), # Center for bon
  c(-97.25, 31.56), # Center for waco
  c(-99.07, 31.86),  # Center for bw
  c(-99.71, 31.55),  # Center for iv
  c(-103.94, 31.95),  # Center for rb
  c(-98.37,33.71)  # Center for ah
)

df_list <- list(
  df_bon = bon2,
  df_waco = waco,
  df_bw = bw,
  df_iv = iv,
  df_rb = rb,
  df_ah = ah
)


map1 <- create_plot(centers[[1]][1], centers[[1]][2], bbox_width_zoomed_in, bbox_height_zoomed_in, df_list$df_bon, "Bonham",add_scalebar=TRUE)
map2 <- create_plot(centers[[2]][1], centers[[2]][2], bbox_width, bbox_height, df_list$df_waco, "Waco") 
map3 <- create_plot(centers[[3]][1], centers[[3]][2], bbox_width, bbox_height, df_list$df_bw, "Brownwood") 
map4 <- create_plot(centers[[4]][1], centers[[4]][2], bbox_width_iv, bbox_height_iv,df_list$df_iv, "OH Ivie") 
map5 <- create_plot(centers[[5]][1], centers[[5]][2], bbox_width, bbox_height, df_list$df_rb, "Red Bluff") 
map6 <- create_plot(centers[[6]][1], centers[[6]][2], bbox_width, bbox_height, df_list$df_ah, "Arrowhead",add_scalebar=TRUE) 




# Combine plots into a single layout with grid.arrange
combo_sdd_lomarg<-grid.arrange(map1, map2, map3, map4, map5, map6, ncol = 3, padding = unit(c(-1, -1, -1, -1), "cm"))
ggsave("combo_sdd_improved.png", combo_sdd_lomarg)


# Add scale bar if specified
if (add_scalebar) {
  
  scale_length_km <- 10
  
  plot <- plot +
    annotation_scale(location = "tl", width_hint = 0.5,
                     text_col = "black",  # Color of the scale bar text
                     line_col = "black",  # Color of the scale bar line
                     unit = "km")  # Adjust width_hint as needed
}




# Fetch the map
map <- openmap(upperLeft, lowerRight, type = "bing")

# Convert to ggplot object
#autoplot(map) + theme_minimal()}

plot_sdah<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_ah)) +
  geom_jitter(data = ah, 
              size=ptsize,
              alpha=ptalpha,
              aes(x = lon_dec, y = lat_dec,color=pred_sdd))+
  scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)")+ 
  ggtitle(label = "Arrowhead") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top"))+ 
  theme(plot.title = element_text(hjust = 0.5))


plot_sdbn<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_bn)) +
  geom_jitter(data = bon2, 
              size=ptsize,
              alpha=ptalpha,
              aes(x = lon_dec, y = lat_dec,color=pred_sdd))+
  ggtitle(label = "Bonham") +
  scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)")+ 
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top"))+ 
  theme(plot.title = element_text(hjust = 0.5))

plot_sdbw<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_bw)) +
  geom_jitter(data = bw, 
              size=ptsize,
              alpha=ptalpha,
              aes(x = lon_dec, y = lat_dec,color=pred_sdd))+
  scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)")+ 
  ggtitle(label = "Brownwood") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top"))+ 
  theme(plot.title = element_text(hjust = 0.5))

plot_sdiv<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_iv)) +
  geom_jitter(data = iv, 
              size=ptsize,
              alpha=ptalpha,
              aes(x = lon_dec, y = lat_dec,color=pred_sdd))+
  scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)")+ 
  ggtitle(label = "O.H. Ivie") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top"))+ 
  theme(plot.title = element_text(hjust = 0.5))

plot_sdrb<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_rb)) +
  geom_jitter(data = rb, 
              size=ptsize,
              alpha=ptalpha,
              aes(x = lon_dec, y = lat_dec,color=pred_sdd))+
  scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)")+ 
  ggtitle(label = "Red Bluff") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top"))+ 
  theme(plot.title = element_text(hjust = 0.5))

plot_sdlw<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_lw)) +
  geom_jitter(data = waco, 
              size=ptsize,
              alpha=ptalpha,
              aes(x = lon_dec, y = lat_dec,color=pred_sdd))+
  scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)")+ 
  ggtitle(label = "Waco") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top")) + 
  theme(plot.title = element_text(hjust = 0.5))




create_map <- function(center_lon, center_lat, bbox_width, bbox_height) {
  # Define the bounding box
  upperLeft <- c(center_lat + bbox_height / 2, center_lon - bbox_width / 2)
  lowerRight <- c(center_lat - bbox_height / 2, center_lon + bbox_width / 2)
  
  # Fetch the map
  map <- openmap(upperLeft, lowerRight, type = "esri_imagery")

    #map <- openmap(upperLeft, lowerRight, type = "osm")

  
  # Convert to ggplot object
  autoplot(map) + theme_minimal()
}


  # create maps
  map1 <- create_map(centers[[1]][1], centers[[1]][2], bbox_width, bbox_height) + ggtitle("Map 1")
map2 <- create_map(centers[[2]][1], centers[[2]][2], bbox_width, bbox_height) + ggtitle("Map 2")
map3 <- create_map(centers[[3]][1], centers[[3]][2], bbox_width, bbox_height) + ggtitle("Map 3")
map4 <- create_map(centers[[4]][1], centers[[4]][2], bbox_width, bbox_height) + ggtitle("Map 4")
map5 <- create_map(centers[[5]][1], centers[[5]][2], bbox_width, bbox_height) + ggtitle("Map 5")
map6 <- create_map(centers[[6]][1], centers[[6]][2], bbox_width, bbox_height) + ggtitle("Map 6")


# Combine plots into a single layout with grid.arrange
grid.arrange(map1, map2, map3, map4, map5, map6, ncol = 3)

# Optionally, add a common scale bar using cowplot
scale_bar <- ggplot() + 
  theme_void() + 
  theme(legend.position = "bottom") +
  guides(fill = guide_colourbar(barwidth = 10, barheight = 1))

combined_plot <- plot_grid(
  plot_grid(map1, map2, map3, map4, map5, map6, ncol = 3, align = 'v'),
  scale_bar,
  ncol = 1,
  rel_heights = c(1, 0.1)
)

print(combined_plot)
##########################################################################################
flm_sdd<-read_csv("flm_sdd.csv")

# set some parameters to pass into ggplot 
ptsize<-0.25
ptalpha<-0.5


# attempts to have a common scale bar and similar area covered in bounding box

# Define the bounding box size (latitude and longitude degrees)
bbox_width <- 0.25
bbox_height <- 0.25
# Define the center coordinates for each map
centers <- list(
  c(-96.15, 33.65), # Center for bon
  c(-97.22, 31.57), # Center for waco
  c(-99.03, 31.83),  # Center for bw
  c(-99.68, 31.55),  # Center for iv
  c(-103.94, 31.95),  # Center for rb
  c(-98.37,33.71)  # Center for ah
)



df_bon <- data.frame(lon = c(-96.1, -96.2), lat = c(33.6, 33.7))
df_waco <- data.frame(lon = c(-97.2, -97.3), lat = c(31.5, 31.6))
df_bw <- data.frame(lon = c(-99.0, -99.1), lat = c(31.8, 31.9))
df_iv <- data.frame(lon = c(-99.6, -99.7), lat = c(31.5, 31.6))
df_rb <- data.frame(lon = c(-103.9, -104.0), lat = c(31.9, 32.0))
df_ah <- data.frame(lon = c(-98.3, -98.4), lat = c(33.7, 33.8))

create_plot <- function(center_lon, center_lat, bbox_width, bbox_height, df, title, add_scalebar = FALSE) {
  # Define the bounding box
  upperLeft <- c(center_lat + bbox_height / 2, center_lon - bbox_width / 2)
  lowerRight <- c(center_lat - bbox_height / 2, center_lon + bbox_width / 2)
  
  # Fetch the map
  map <- OpenStreetMap::openmap(upperLeft, lowerRight, type = "esri_imagery")
  map_projected <- OpenStreetMap::openproj(map)
  
  # Create the plot
  plot <- OpenStreetMap::autoplot.OpenStreetMap(map_projected) +
    geom_jitter(data = df, 
                size = psize,
                alpha = ptalpha,
                aes(x = lon, y = lat, color = pred_sdd)) +
    scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)") +
    ggtitle(label = title) +
    xlab("") + ylab("") +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.direction = "horizontal",
      legend.title.align = 0.5,
      legend.key.width = unit(0.75, "cm"),
      legend.position = c(0.5, 0.03)) +
    guides(col = guide_colorbar(title.position = "top")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  # Add scale bar if specified
  if (add_scalebar) {
    plot <- plot +
      scalebar(
        dist = 1,  # Length of the scale bar in kilometers
        dist_unit = "km",
        transform = TRUE,
        model = "WGS84",
        location = "bottomleft",
        units = "km",
        scale = 0.5,  # Adjust size of scale bar
        st.size = 3   # Adjust size of scale bar text
      )
  }
  
  return(plot)
}
# create maps
plot_bon <- create_plot(centers[[1]][1], centers[[1]][2], bbox_width, bbox_height, df_bon, "Bonham", add_scalebar = TRUE)
plot_waco <- create_plot(centers[[2]][1], centers[[2]][2], bbox_width, bbox_height, waco, "Waco")
plot_bw <- create_plot(centers[[3]][1], centers[[3]][2], bbox_width, bbox_height, bw, "Brownwood")
plot_iv <- create_plot(centers[[4]][1], centers[[4]][2], bbox_width, bbox_height, iv, "OH Ivie")
plot_rb <- create_plot(centers[[5]][1], centers[[5]][2], bbox_width, bbox_height, rb, "Red Bluff")
plot_ah <- create_plot(centers[[6]][1], centers[[6]][2], bbox_width, bbox_height, ah, "Arrowhead",add_scalebar = TRUE)
#####################################################################################################

# set some parameters to pass into ggplot 
ptsize<-0.25
ptalpha<-0.5

# attempts to have a common scale bar and similar area covered in bounding box

# Define the bounding box size (latitude and longitude degrees)
bbox_width_zoomed_in <- 0.05
bbox_height_zoomed_in <- 0.05
bbox_width <- 0.15
bbox_height <- 0.15
bbox_width_iv <- 0.25
bbox_height_iv <- 0.25


create_plot <- function(center_lon, center_lat, bbox_width, bbox_height, df, title, add_scalebar=FALSE) {
  # Define the bounding box
  upperLeft <- c(center_lat + bbox_height / 2, center_lon - bbox_width / 2)
  lowerRight <- c(center_lat - bbox_height / 2, center_lon + bbox_width / 2)
  
  # Fetch the map
  map <- openmap(upperLeft, lowerRight, type = "bing")
  map_projected <- OpenStreetMap::openproj(map)
  
  # Create the plot
  plot <- OpenStreetMap::autoplot.OpenStreetMap(map_projected) +
    geom_jitter(data = df, 
                size = ptsize,
                alpha = ptalpha,
                aes(x = lon_dec, y = lat_dec, color = pred_sdd)) +
    scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)") +
    ggtitle(label = title) +
    xlab("") + ylab("") +
    theme(
      plot.margin = unit(c(-1, 0, -1, 0), "lines"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.direction = "horizontal",
      legend.title.align = 0.5,
      legend.key.width = unit(0.75, "cm"),
      legend.position = c(0.5, 0.02)
    ) +
    guides(col = guide_colorbar(title.position = "top")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  return(plot)
}

# Define the center coordinates for each map
centers <- list(
  c(-96.15, 33.65), # Center for bon
  c(-97.25, 31.56), # Center for waco
  c(-99.07, 31.86),  # Center for bw
  c(-99.71, 31.55),  # Center for iv
  c(-103.94, 31.95),  # Center for rb
  c(-98.37,33.71)  # Center for ah
)

df_list <- list(
  df_bon = bon2,
  df_waco = waco,
  df_bw = bw,
  df_iv = iv,
  df_rb = rb,
  df_ah = ah
)


map1 <- create_plot(centers[[1]][1], centers[[1]][2], bbox_width_zoomed_in, bbox_height_zoomed_in, df_list$df_bon, "Bonham",add_scalebar=TRUE)
map2 <- create_plot(centers[[2]][1], centers[[2]][2], bbox_width, bbox_height, df_list$df_waco, "Waco") 
map3 <- create_plot(centers[[3]][1], centers[[3]][2], bbox_width, bbox_height, df_list$df_bw, "Brownwood") 
map4 <- create_plot(centers[[4]][1], centers[[4]][2], bbox_width_iv, bbox_height_iv,df_list$df_iv, "OH Ivie") 
map5 <- create_plot(centers[[5]][1], centers[[5]][2], bbox_width, bbox_height, df_list$df_rb, "Red Bluff") 
map6 <- create_plot(centers[[6]][1], centers[[6]][2], bbox_width, bbox_height, df_list$df_ah, "Arrowhead",add_scalebar=TRUE) 




# Combine plots into a single layout with grid.arrange
combo_sdd_lomarg<-grid.arrange(map1, map2, map3, map4, map5, map6, ncol = 3, padding = unit(c(-1, -1, -1, -1), "cm"))
ggsave("combo_sdd_improved.png", combo_sdd_lomarg)


# Add scale bar if specified
if (add_scalebar) {
  
  scale_length_km <- 10
  
  plot <- plot +
    annotation_scale(location = "tl", width_hint = 0.5,
                     text_col = "black",  # Color of the scale bar text
                     line_col = "black",  # Color of the scale bar line
                     unit = "km")  # Adjust width_hint as needed
}




# Fetch the map
map <- openmap(upperLeft, lowerRight, type = "bing")

# Convert to ggplot object
#autoplot(map) + theme_minimal()}

plot_sdah<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_ah)) +
  geom_jitter(data = ah, 
              size=ptsize,
              alpha=ptalpha,
              aes(x = lon_dec, y = lat_dec,color=pred_sdd))+
  scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)")+ 
  ggtitle(label = "Arrowhead") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top"))+ 
  theme(plot.title = element_text(hjust = 0.5))


plot_sdbn<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_bn)) +
  geom_jitter(data = bon2, 
              size=ptsize,
              alpha=ptalpha,
              aes(x = lon_dec, y = lat_dec,color=pred_sdd))+
  ggtitle(label = "Bonham") +
  scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)")+ 
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top"))+ 
  theme(plot.title = element_text(hjust = 0.5))

plot_sdbw<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_bw)) +
  geom_jitter(data = bw, 
              size=ptsize,
              alpha=ptalpha,
              aes(x = lon_dec, y = lat_dec,color=pred_sdd))+
  scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)")+ 
  ggtitle(label = "Brownwood") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top"))+ 
  theme(plot.title = element_text(hjust = 0.5))

plot_sdiv<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_iv)) +
  geom_jitter(data = iv, 
              size=ptsize,
              alpha=ptalpha,
              aes(x = lon_dec, y = lat_dec,color=pred_sdd))+
  scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)")+ 
  ggtitle(label = "O.H. Ivie") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top"))+ 
  theme(plot.title = element_text(hjust = 0.5))

plot_sdrb<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_rb)) +
  geom_jitter(data = rb, 
              size=ptsize,
              alpha=ptalpha,
              aes(x = lon_dec, y = lat_dec,color=pred_sdd))+
  scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)")+ 
  ggtitle(label = "Red Bluff") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top"))+ 
  theme(plot.title = element_text(hjust = 0.5))

plot_sdlw<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_lw)) +
  geom_jitter(data = waco, 
              size=ptsize,
              alpha=ptalpha,
              aes(x = lon_dec, y = lat_dec,color=pred_sdd))+
  scale_color_gradientn(colors = my_palette(100), name = "Secchi Disk Depth (m)")+ 
  ggtitle(label = "Waco") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top")) + 
  theme(plot.title = element_text(hjust = 0.5))




create_map <- function(center_lon, center_lat, bbox_width, bbox_height) {
  # Define the bounding box
  upperLeft <- c(center_lat + bbox_height / 2, center_lon - bbox_width / 2)
  lowerRight <- c(center_lat - bbox_height / 2, center_lon + bbox_width / 2)
  
  # Fetch the map

  map <- openmap(upperLeft, lowerRight, type = "esri_imagery")

    #map <- openmap(upperLeft, lowerRight, type = "osm")

  
  # Convert to ggplot object
  autoplot(map) + theme_minimal()
}


  # create maps
  map1 <- create_map(centers[[1]][1], centers[[1]][2], bbox_width, bbox_height) + ggtitle("Map 1")
map2 <- create_map(centers[[2]][1], centers[[2]][2], bbox_width, bbox_height) + ggtitle("Map 2")
map3 <- create_map(centers[[3]][1], centers[[3]][2], bbox_width, bbox_height) + ggtitle("Map 3")
map4 <- create_map(centers[[4]][1], centers[[4]][2], bbox_width, bbox_height) + ggtitle("Map 4")
map5 <- create_map(centers[[5]][1], centers[[5]][2], bbox_width, bbox_height) + ggtitle("Map 5")
map6 <- create_map(centers[[6]][1], centers[[6]][2], bbox_width, bbox_height) + ggtitle("Map 6")


# Combine plots into a single layout with grid.arrange
grid.arrange(map1, map2, map3, map4, map5, map6, ncol = 3)

# Optionally, add a common scale bar using cowplot
scale_bar <- ggplot() + 
  theme_void() + 
  theme(legend.position = "bottom") +
  guides(fill = guide_colourbar(barwidth = 10, barheight = 1))

combined_plot <- plot_grid(
  plot_grid(map1, map2, map3, map4, map5, map6, ncol = 3, align = 'v'),
  scale_bar,
  ncol = 1,
  rel_heights = c(1, 0.1)
)

print(combined_plot)


########################################################################################################


########################################################################################################
dm_sf <- st_as_sf(bwfd, coords = c("lon_dec", "lat_dec"), crs = 4326)
mapview(dm_sf, zcol= "fdom")
#filtering for fdom
wfd<- waco %>% filter(fdom>6)
rfd<- rb %>% filter(fdom>5.7)
ifd<-iv %>% filter(fdom>6)
bwfd<-bw %>% filter(fdom>7)
ahfd<-ah %>% filter(fdom>9)
# plotting fdom
plot_fdah<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_ah)) +
  geom_point(data = ahfd, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon_dec, y = lat_dec,color=fdom))+
  scale_color_viridis_c("fDOM (RFU)",direction=-1)+
  ggtitle(label = "Arrowhead, 08-12-2022") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top"))+ 
  theme(plot.title = element_text(hjust = 0.5))

plot_fdbn1<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_bon)) +
  geom_point(data = ah, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon_dec, y = lat_dec,color=fdom))+
  scale_color_viridis_c("fDOM (RFU)",direction=-1)+
  ggtitle(label = "Arrowhead, 08-12-2022") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top"))+ 
  theme(plot.title = element_text(hjust = 0.5))



plot_fdbn2<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_bon)) +
  geom_point(data = bon2, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon_dec, y = lat_dec,color=fdom))+
  ggtitle(label = "Bonham, 07-15-2022") +
  scale_color_viridis_c("fDOM (RFU)",direction=-1)+
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top"))+ 
  theme(plot.title = element_text(hjust = 0.5))

plot_fdbw<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_bw)) +
  geom_point(data = bwfd, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon_dec, y = lat_dec,color=fdom))+
  scale_color_viridis_c("fDOM (RFU)",direction=-1)+
  ggtitle(label = "Brownwood, 08-06-2022") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top"))+ 
  theme(plot.title = element_text(hjust = 0.5))

plot_fdiv<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_iv)) +
  geom_point(data = ifd, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon_dec, y = lat_dec,color=fdom))+
  scale_color_viridis_c("fDOM (RFU)",direction=-1)+
  ggtitle(label = "O.H. Ivie, 08-07-2022") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top"))+ 
  theme(plot.title = element_text(hjust = 0.5))

plot_fdrb<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_rb)) +
  geom_point(data = rfd, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon_dec, y = lat_dec,color=fdom))+
  scale_color_viridis_c("fDOM (RFU)",direction=-1)+
  ggtitle(label = "Red Bluff, 08-08-2022") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top"))+ 
  theme(plot.title = element_text(hjust = 0.5))

plot_fdlw<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_lw)) +
  geom_point(data = wfd, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon_dec, y = lat_dec,color=fdom))+
  scale_color_viridis_c("fDOM (RFU)",direction=-1)+
  ggtitle(label = "Waco, 07-23-2022") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top")) + 
  theme(plot.title = element_text(hjust = 0.5))

# Combine maps and save as external .png image

comb_fdom_plot<-grid.arrange(ncol=3,
                            plot_fdbn2,plot_fdlw,
                            plot_fdbw,
                            plot_fdiv,plot_fdrb,plot_fdah,
                            padding = unit(0.0, "line"))

ggsave(filename = "combo_fdom_viridis.png", plot = comb_fdom_plot, width = 10, height = 5)


###########################################################################################
# Load your GPX file
gpx_file <- "AllLogsThrough3May2023.gpx"
gpx <- xmlParse(gpx_file)

track <- htmlTreeParse(file = "AllLogsThrough3May2023.gpx", error = function(...)
{ }, useInternalNodes = T)

# Define namespaces
namespaces <- c(
  gpxtpx = "http://www.garmin.com/xmlschemas/TrackPointExtension/v1",
  gpxx = "http://www.garmin.com/xmlschemas/GpxExtensions/v3"
)
# Extract depths from gpxtpx:TrackPointExtension
depths_gpxtpx <- xpathSApply(gpx, "//gpxtpx:TrackPointExtension/gpxtpx:depth", xmlValue, namespaces = namespaces)
# Extract depths from gpxx:TrackPointExtension
depths_gpxx <- xpathSApply(gpx, "//gpxx:TrackPointExtension/gpxx:Depth", xmlValue, namespaces = namespaces)
# Combine both depth values
depths <- c(depths_gpxtpx, depths_gpxx)
# Convert to numeric
depths <- as.numeric(depths)
# Print depths to check
print(depths)

times <- xpathSApply(track, path = "//trkpt/time", xmlValue)
coords <- xpathSApply(track, path = "//trkpt", xmlAttrs)
## Extract lat and lon from coordinates
lats <- as.numeric(coords["lat",])
lons <- as.numeric(coords["lon",])

# Ensure all vectors have the same length by removing mismatched entries
min_length <- min(length(lats), length(lons), length(times), length(depths))
lats <- lats[1:min_length]
lons <- lons[1:min_length]
times <- times[1:min_length]
depths <- depths[1:min_length]


gpx <- data.frame(lat = lats, lon  = lons, time = times , depth = depths)
gpx$time <- as.POSIXct(gpx$time)

bon_gpx<-gpx %>% dplyr::filter(time>"2022-07-14" & time<"2022-07-16")
waco_gpx<-gpx %>% dplyr::filter(time>"2022-07-22" & time<"2022-07-24")
bw_gpx<-gpx%>% dplyr::filter(time>"2022-08-05" & time<"2022-08-07")
iv_gpx<-gpx%>% dplyr::filter(time>"2022-08-06" & time<"2022-08-08")
rb_gpx<-gpx%>% dplyr::filter(time>"2022-08-07" & time<"2022-08-09")
ah_gpx<-gpx%>% dplyr::filter(time>"2022-08-11" & time<"2022-08-13")

bon_depth<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_bon)) +
  geom_point(data = bon_gpx, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon, y = lat,color=depth))+
  scale_color_viridis_c("Depth (ft)",direction=-1)+
  ggtitle(label = "Bonham") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top")) + 
  theme(plot.title = element_text(hjust = 0.5))

lw_depth<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_lw)) +
  geom_point(data = waco_gpx, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon, y = lat,color=depth))+
  scale_color_viridis_c("Depth (ft)",direction=-1)+
  ggtitle(label = "Waco") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top")) + 
  theme(plot.title = element_text(hjust = 0.5))

bw_depth<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_bw)) +
  geom_point(data = bw_gpx, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon, y = lat,color=depth))+
  scale_color_viridis_c("Depth (ft)",direction=-1)+
  ggtitle(label = "Brownwood") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top")) + 
  theme(plot.title = element_text(hjust = 0.5))

iv_depth<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_iv)) +
  geom_point(data = iv_gpx, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon, y = lat,color=depth))+
  scale_color_viridis_c("Depth (ft)",direction=-1)+
  ggtitle(label = "OH Ivie") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top")) + 
  theme(plot.title = element_text(hjust = 0.5))

rb_depth<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_rb)) +
  geom_point(data = rb_gpx, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon, y = lat,color=depth))+
  scale_color_viridis_c("Depth (ft)",direction=-1)+
  ggtitle(label = "Red Bluff") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top")) + 
  theme(plot.title = element_text(hjust = 0.5))

ah_depth<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_ah)) +
  geom_point(data = ah_gpx, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon, y = lat,color=depth))+
  scale_color_viridis_c("Depth (ft)",direction=-1)+
  ggtitle(label = "Arrowhead") +
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.75,"cm"),
    legend.position = c(0.5, 0.03))+
  guides(col=guide_colorbar(title.position = "top")) + 
  theme(plot.title = element_text(hjust = 0.5))

# Combine maps and save as external .png image
comb_depth_plot<-grid.arrange(ncol=3,
                             bon_depth,lw_depth,
                             bw_depth,
                             iv_depth,rb_depth,ah_depth,
                             padding = unit(0.0, "line")) 
ggsave(filename = "combo_depth.png", plot = comb_depth_plot, width = 10, height = 5)

dm_sf <- st_as_sf(rb_gpx, coords = c("lon", "lat"), crs = 4326)
mapview(dm_sf, zcol= "depth")
