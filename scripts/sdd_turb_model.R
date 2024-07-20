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
devtools::install_github("16EAGLE/basemaps")

# 2. CRASR data - read in for all lakes -------------------------------------------
crasr_all <- read_csv("CRASR_summer2022.csv") %>% clean_names()

# separate by lake day
crasr_bon1 <- crasr_all %>% filter(datecollect == "7/14/2022") %>% add_column(lake = "bon") %>% group_by(system,lat, lon, variable, site) %>% mutate(meanchl = mean(value))
cbon1 <- crasr_bon1 %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl) %>% dplyr::select(lat, lon, system, site, chla = `chl a`, turbidity, tss, afdm, tp,tn, secchi)

crasr_bon2 <- crasr_all %>% filter(datecollect == "7/15/2022") %>% add_column(lake = "bon") %>% 
  dplyr::mutate(site = paste(lake, site, sep = "_")) %>% group_by(system,lat, lon, variable, site) %>% mutate(meanchl = mean(value))
cbon2 <- crasr_bon2 %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl)%>% dplyr::select(lat, lon, site, chla = `chl a`, turbidity, tss, afdm, tp,tn, secchi)

crasr_lw <- crasr_all %>% filter(datecollect == "7/23/2022")%>% add_column(lake = "lw") %>% 
  dplyr::mutate(station = paste(lake, site, sep = "")) %>% group_by(system,lat, lon, variable, station) %>% mutate(meanchl = mean(value))
clw <- crasr_lw %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl)%>% dplyr::select(lat, lon, station, chla = `chl a`, turbidity, tss, afdm, tp,tn)

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
  
#tsd_plot <- ggplot(nla_crasr,aes(log(turb), log(secchi))) + 
 # geom_point() + 
  #geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  #stat_poly_eq(formula = y ~ x, 
   #            aes(label = paste(after_stat(rr.label))), 
    #           parse = TRUE, label.x.npc = "left", size = 8, rr.digits = 3)+
  #ggtitle(label = "Turbidity ~ Secchi for NLA Reservoirs 2012 and 2017")+
  #xlab("log(Turbidity)") + ylab("log(Secchi Depth)") +
  #theme_bw() +theme(plot.title = element_text(size = 20,hjust=0.5),
   #                 axis.title.x = element_text(size = 15),
    #                axis.title.y = element_text(size = 15),
     #               legend.title = element_text(size = 15))

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
ah$secchi <- exp(predict(tsdlm, ah))
bon1$secchi <- exp(predict(tsdlm, bon1))
bon2$secchi <- exp(predict(tsdlm, bon2))
waco$secchi <- exp(predict(tsdlm, waco))
bw$secchi <- exp(predict(tsdlm, bw))
iv$secchi <- exp(predict(tsdlm, iv))
rb$secchi <- exp(predict(tsdlm, rb))

# calculating TSI from secchi
ah$tsi_sd <- 10 * (6-log(ah$secchi)/log(2))
bon1$tsi_sd <- 10 * (6-log(bon1$secchi)/log(2))
bon2$tsi_sd <- 10 * (6-log(bon2$secchi)/log(2))
bw$tsi_sd <- 10 * (6-log(bw$secchi)/log(2))
iv$tsi_sd <- 10 * (6-log(iv$secchi)/log(2))
rb$tsi_sd <- 10 * (6-log(rb$secchi)/log(2))
waco$tsi_sd <- 10 * (6-log(waco$secchi)/log(2))

# removing outlier value in iv and waco
iv <- iv %>% filter(iv$turb<25) %>% filter(date_time_c<"2022-08-07 13:28")
waco <- waco %>% filter(waco$turb<75)
# 7. visualizing secchi and TSI with maps-------------------------------------------------------
dm_sf <- st_as_sf(iv, coords = c("lon_dec", "lat_dec"), crs = 4326)
mapview(dm_sf, zcol= "tsi_sd")
# set some parameters to pass into ggplot 
ptsize<-0.9
ptalpha<-0.5

# setting gps boundaries for each lake-map
map_ah <- openmap(upperLeft = c(33.78, -98.462), 
                lowerRight = c(33.615, -98.312), 
                type = 'osm',zoom=12)

map_bon <- openmap(upperLeft = c(33.675,-96.175),
                  lowerRight = c(33.625,-96.125), 
                  type = 'osm',zoom=12)

map_bw <- openmap(upperLeft = c(31.93, -99.11),
                  lowerRight = c(31.74, -98.985), 
                  type = 'osm',zoom=12)

map_iv <- openmap(upperLeft = c(31.62, -99.84),
                  lowerRight = c(31.4, -99.6), 
                  type = 'osm',zoom=12)

map_rb <- openmap(upperLeft = c(32.026, -104.03),
                  lowerRight = c(31.84, -103.879), 
                  type = 'osm',zoom=12)

map_lw <- openmap(upperLeft = c(31.625,-97.3),
                lowerRight = c(31.47,-97.19),
                type = 'osm',zoom=12)


my_palette <- colorRampPalette(rev(brewer.pal(9, "YlOrBr")))

# plotting secchi
plot_sdah<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_ah)) +
  geom_point(data = ah, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon_dec, y = lat_dec,color=secchi))+
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
  geom_point(data = bon1, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon_dec, y = lat_dec,color=secchi))+
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
  geom_point(data = bon2, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon_dec, y = lat_dec,color=secchi))+
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
  geom_point(data = bw, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon_dec, y = lat_dec,color=secchi))+
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
  geom_point(data = iv, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon_dec, y = lat_dec,color=secchi))+
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
  geom_point(data = rb, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon_dec, y = lat_dec,color=secchi))+
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
  geom_point(data = waco, 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon_dec, y = lat_dec,color=secchi))+
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
                            plot_sdiv,plot_sdrb,plot_sdah)

ggsave(filename = "combo_secchi_flm.png", plot = combined_plot, width = 10, height = 5)

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
