# A script for calculating running linear models for satellite and in situ data
# Summer 2022 reservoirs

# clear workspace
rm(list=ls())

# load packages
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
library(sf)
library(OpenStreetMap)
library(leaflet)
library(leaflet.extras)
library(grid)
library(gridExtra)
#library(maptools)
library(maps)
#library(ggsn)
library(ggspatial)
library(factoextra)
library(automap)
library(sp)
library(lattice)
library(rioja)
library(zoo)
library(depmixS4)
library(ppclust)
library(factoextra)
library(cluster)
library(fclust)
library(mapview)
library(ggpubr)
library(plotly)
library(wesanderson)
library(ghibli)
install.packages("devtools")
devtools::install_github("johannesbjork/LaCroixColoR")
library(LaCroixColoR)
library(ggdist)
library(ggsignif)
library(XML)
library(gpx)
library(patchwork)
library(nlme)
library(anytime)



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

garmin<-garmin %>% dplyr::select(garmin_datetime_utc,date,hour_dec,garmin_time_local,xcoord,ycoord,speed)
#getUTMzone(mean(garmin$xcoord,na.rm=TRUE))


# 3. Merge data

data_merged<-merge(exo,garmin, 
                   by.x=c("lake_time_local","date"), 
                   by.y=c("garmin_time_local", "date"))


data_merged<-data_merged %>% dplyr::select(lat_dec=ycoord,lon_dec=xcoord, speed,
                                           time_local=lake_time_local,date,datetime_utc=garmin_datetime_utc,
                                           wtemp=temp_c, chl_rfu=chlorophyll_rfu,chl_ugl=chlorophyll_ugl,
                                           pc=tal_pc_rfu,fdom=f_dom_rfu,turb=turbidity_fnu,do_mgl=odo_mg_l, 
                                           do_sat=odo_percent, hour_dec
)

data_merged$date_time_c = paste(data_merged$date, data_merged$time_local)

dm_posc <- data_merged %>% filter(chl_ugl>0)




# 4. filtering garmin data by lake / date
ah<-data_merged %>% filter(date=="2022-08-12", hour_dec>9.8, hour_dec<12.8)

# assigning unique id label 
ah_id <- ah %>%
  mutate(label = row_number()) 
ah_id$system <- "arrowhead"

ah_unique<-ah_id %>% distinct(lat_dec,lon_dec, .keep_all=TRUE)

##########################################################################################################
# writing out lat/lon for df
ah_pts <- ah_id %>% dplyr::select(label,lat_dec,lon_dec) %>% na.omit()


# Split the full df into smaller subsets with <=5000 points each
ah_subsets <- split(ah_pts, rep(1:ceiling(nrow(ah_pts)/2000), each=2000, length.out=nrow(ah_pts)))
# Export each subset as a separate CSV file
for (i in seq_along(ah_subsets)) {
  subset <- ah_subsets[[i]]
  write.csv(subset, file = paste0("ahpts_", i, ".csv"), row.names = FALSE)
}
##########################################################################################################


# import csv results from eepoints
ah0 <- read_csv(c("s2bvs_ah1.csv","s2bvs_ah2.csv","s2bvs_ah3.csv","s2bvs_ah4.csv",
                  "s2bvs_ah5.csv","s2bvs_ah6.csv")) %>% dplyr::arrange(label)

# filtering just to bands that will be used and calculating normalized diffs

ah_nd <- ah0 %>% dplyr::select(label,lat_dec, lon_dec, B2, B3, B4, B5, SCL) %>% 
  mutate(ndti = (B4-B3)/(B4+B3),
         ndci = (B5-B4)/(B5+B4))


ah_s2flm <- merge(ah_id, ah_nd, by = "label") %>% dplyr::select(label, lat = lat_dec.x, lon = lon_dec.x,
                                                                system, speed, datetime_utc, date_time_c, fdom, wtemp,
                                                                chl_rfu, chl_ugl, turb,B2, B3, B4, B5, SCL, ndti, ndci) #%>% dplyr::filter(speed>0)
# trying sp's way of big merge
#ah_bm <- merge(ah_id, ah_nd)

#ah_bm$lat_absdev <- abs(ah_bm$lat_dec-ah_bm$lat)
#ah_bm$lon_absdev <- abs(ah_bm$lon_dec-ah_bm$lon)
#ah_bm$latlon_absdev <- sum(ah_bm$lat_absdev,ah_bm$lon_absdev)

#ah_merged <- ah_bm%>% 
 # slice_min(order_by = latlon_absdev, n = 1) %>%
  #arrange()

# adding flag for non-water pixels and filtering them out
ah_s2flm$flag<-0
ah_s2flm$flag[which(!ah_s2flm $SCL %in% c(4,5,6,11))]<-1
ah_s2flm$flag[which(ah_s2flm $B2<=500)]<-1

# neighbor flag
flagneighbors_lat<-ah_s2flm$lat[which(ah_s2flm$flag==1)]
flagneighbors_lon<-ah_s2flm$lon[which(ah_s2flm$flag==1)]
# changed to 0.0005 degree radius (45 to 55 m at this latitude)
flagneighbors_which<-c()
for(i in 1:length(flagneighbors_lat)){
  whichi<-which(ah_s2flm$lat_dec<(flagneighbors_lat[i]+0.0005) & 
                  
                  ah_s2flm$lat_dec>(flagneighbors_lat[i]-0.0005) & 
                  
                  ah_s2flm$lon_dec<(flagneighbors_lon[i]+0.0005) & 
                  
                  ah_s2flm$lon_dec>(flagneighbors_lon[i]-0.0005))
  
  flagneighbors_which<-c(flagneighbors_which,whichi)
}  

# add column for keeping track of neighbors of flagged pixels that weren't themselves flagged (yet)
ah_s2flm$flagneighbor<-0
ah_s2flm$flagneighbor[unique(flagneighbors_which)]<-1
ah_s2flm$flagneighbor[which(ah_s2flm$flag==1)]<-0

# now flag pixels that neighbor a flagged pixel   
ah_s2flm$flag[unique(flagneighbors_which)]<-1
ah_s2flm<-ah_s2flm %>% filter(flag==0) 

ah_unique<-ah_s2flm %>% distinct(lat,lon, .keep_all=TRUE)

write.csv(ah_s2flm, "ahpts_nf.csv")
write.csv(ah_unique, "ahpts_unq.csv")

# running linear model between ndti~turb
logEstimate <- lm(ndti~log(turb),data=data_use)
xvec <- seq(min(data_use$turb),max(data_use$turb),length=1000)
logpred <- predict(logEstimate,newdata=data.frame(turb=xvec))
pred <- data.frame(x = xvec, y = logpred)


turb_ndti_plot<-ggplot(data_use,aes(x=turb,y=ndti)) +
  geom_point(size=0.1,col=gray(0.3))+
  geom_smooth(method="lm",lty=5,col="black",se=FALSE,lwd=1)+
  geom_line(data = pred, aes(x=x, y=y),lty=1,col="black",lwd=1)+
  #xlim(c(0,99))+
  xlab("In-lake Turbidity (NTU)")+
  ylab("Satellite Turbidity (NDTI)")+
  theme_bw()
turb_ndti_plot

ggplotly(turb_ndti_plot)


ggplot(data_use,aes(log(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()



# removing negative chl-a
data_posc <- data_use %>% dplyr::filter(chl >0)
# visualizing linear model ndci~chla
logEst_c <- lm(ndci~chl,data=data_posc)
xvec <- seq(min(data_use$chl),max(data_use$chl),length=1000)
logpred <- predict(logEst_c,newdata=data.frame(chl=xvec))
pred <- data.frame(x = xvec, y = logpred)
summary(logEst_c)

chl_ndci_plot<-ggplot(data_posc,aes(x=chl,y=ndci)) +
  geom_point(size=0.1,col=gray(0.3))+
  geom_smooth(method="lm",lty=5,col="black",se=FALSE,lwd=1)+
  geom_line(data = pred, aes(x=x, y=y),lty=1,col="black",lwd=1)+
  #xlim(c(0,99))+
  xlab("In-lake chl-a (RFU)")+
  ylab("Satellite chl-a (NDCI)")+
  theme_bw()
chl_ndci_plot

ggplotly(chl_ndci_plot)




################################ doing the same for other lakes#######################
### Lake Bonham - 7/15/22
# first filtering to time on the lake
bon<-data_merged %>% filter(date=="2022-07-15", hour_dec>12, hour_dec<13.91)
# assigning unique id label 
bon_id <- bon %>%
  mutate(label = row_number()) 
bon_id$system <- "bonham"

bon_unique<-bon_id %>% distinct(lat_dec,lon_dec, .keep_all=TRUE)

## writing out lat/lon for all dfs
bon_pts <- bon_id %>% dplyr::select(label,lat_dec,lon_dec) %>% na.omit()
# Split the full df into smaller subsets with <=5000 points each
bon_subsets <- split(bon_pts, rep(1:ceiling(nrow(bon_pts)/5000), each=5000, length.out=nrow(bon_pts)))
# Export each subset as a separate CSV file
for (i in seq_along(bon_subsets)) {
  subset <- bon_subsets[[i]]
  write.csv(subset, file = paste0("bonpts_", i, ".csv"), row.names = FALSE)
}


# import csv results from eepoints
bon0 <- read_csv(c("s2bvs_bon1.csv","s2bvs_bon2.csv")) %>% dplyr::arrange(label)

# filtering just to bands that will be used and calculating normalized diffs

bon_nd <- bon0 %>% dplyr::select("label", "lat_dec", "lon_dec", "B2", "B3", "B4", "B5", "SCL") %>% 
  mutate(ndti = (B4-B3)/(B4+B3),
         ndci = (B5-B4)/(B5+B4))


bon_s2flm <- merge(bon_id, bon_nd, by = "label") %>% dplyr::select(label, system, lat = lat_dec.x, lon = lon_dec.x,
                                                                speed, datetime_utc, date_time_c, chl_rfu, chl_ugl, turb,fdom, wtemp,
                                                                B2, B3, B4, B5, SCL, ndti, ndci)#%>% dplyr::filter(speed>0)

# adding flag for non-water pixels and filtering them out
bon_s2flm$flag<-0
bon_s2flm$flag[which(!bon_s2flm $SCL %in% c(4,5,6,11))]<-1
#bon_s2flm$flag[which(bon_s2flm $B2<=500)]<-1

# neighbor flag
flagneighbors_lat<-bon_s2flm$lat[which(bon_s2flm$flag==1)]
flagneighbors_lon<-bon_s2flm$lon[which(bon_s2flm$flag==1)]
# changed to 0.0005 degree radius (45 to 55 m at this latitude)
flagneighbors_which<-c()
for(i in 1:length(flagneighbors_lat)){
  whichi<-which(bon_s2flm$lat_dec<(flagneighbors_lat[i]+0.0005) & 
                  
                  bon_s2flm$lat_dec>(flagneighbors_lat[i]-0.0005) & 
                  
                  bon_s2flm$lon_dec<(flagneighbors_lon[i]+0.0005) & 
                  
                  bon_s2flm$lon_dec>(flagneighbors_lon[i]-0.0005))
  
  flagneighbors_which<-c(flagneighbors_which,whichi)
}  

# add column for keeping track of neighbors of flagged pixels that weren't themselves flagged (yet)
bon_s2flm$flagneighbor<-0
bon_s2flm$flagneighbor[unique(flagneighbors_which)]<-1
bon_s2flm$flagneighbor[which(bon_s2flm$flag==1)]<-0

# now flag pixels that neighbor a flagged pixel   
bon_s2flm$flag[unique(flagneighbors_which)]<-1

data_bon<-bon_s2flm %>% filter(flag==0) %>% filter(SCL == 6) 
db_unq<-data_bon %>% distinct(lat,lon, .keep_all=TRUE)
#data_bon<-bon_s2flm %>% filter(SCL == 6) 

write.csv(data_bon, "bnpts_nf.csv")
write.csv(db_unq, "bnpts_unq.csv")


# running linear model between ndti~turb
logEstimate <- lm(ndti~log(turb),data=data_bon)
xvec <- seq(min(data_bon$turb),max(data_bon$turb),length=1000)
logpred <- predict(logEstimate,newdata=data.frame(turb=xvec))
pred <- data.frame(x = xvec, y = logpred)
summary(logEstimate)

turb_ndti_plot<-ggplot(data_bon,aes(x=turb,y=ndti)) +
  geom_point(size=0.1,col=gray(0.3))+
  geom_smooth(method="lm",lty=5,col="black",se=FALSE,lwd=1)+
  geom_line(data = pred, aes(x=x, y=y),lty=1,col="black",lwd=1)+
  #xlim(c(0,99))+
  xlab("In-lake Turbidity (NTU)")+
  ylab("Satellite Turbidity (NDTI)")+
  theme_bw()
turb_ndti_plot

ggplotly(turb_ndti_plot)

blm <- lm(ndti~log(turb), data_bon)

ggplot(data_use,aes(turb, ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()





### Lake Waco - 7/23/22
# first filtering to time on the lake
waco<-data_merged %>% filter(date=="2022-07-23", hour_dec>10.767, hour_dec<13.25)
# assigning unique id label 
waco_id <- waco %>%
  mutate(label = row_number()) 
waco_id$system <- "waco"
## writing out lat/lon for all dfs
waco_pts <- waco_id %>% dplyr::select(label,lat_dec,lon_dec) %>% na.omit()
# Split the full df into smaller subsets with <=5000 points each
waco_subsets <- split(waco_pts, rep(1:ceiling(nrow(waco_pts)/5000), each=5000, length.out=nrow(waco_pts)))
# Export each subset as a separate CSV file
for (i in seq_along(waco_subsets)) {
  subset <- waco_subsets[[i]]
  write.csv(subset, file = paste0("wacopts_", i, ".csv"), row.names = FALSE)
}

# import csv results from eepoints
lw0 <- read_csv(c("s2bvs_waco1.csv","s2bvs_waco2.csv")) %>% dplyr::arrange(label)

lw_nd <- lw0 %>% dplyr::select("label", "lat_dec", "lon_dec", "B2", "B3", "B4", "B5", "SCL") %>% 
  mutate(ndti = (B4-B3)/(B4+B3),
         ndci = (B5-B4)/(B5+B4))


waco_s2flm <- merge(waco_id, lw_nd, by = "label") %>% dplyr::select(label, system, lat = lat_dec.x, lon = lon_dec.x,
                                                                    datetime_utc, date_time_c, speed,chl_rfu, chl_ugl, turb,fdom, wtemp,
                                                                    B2, B3, B4, B5, SCL, ndti, ndci) %>% filter(turb<20) #%>% dplyr::filter(SCL==6)

#lwlm<-lm(ndti~log(turb), data_use)

# adding flag for non-water pixels and filtering them out
waco_s2flm$flag<-0
waco_s2flm$flag[which(!waco_s2flm $SCL %in% c(4,5,6,11))]<-1
#waco_s2flm$flag[which(waco_s2flm $B2<=500)]<-1

# neighbor flag
flagneighbors_lat<-waco_s2flm$lat[which(waco_s2flm$flag==1)]
flagneighbors_lon<-waco_s2flm$lon[which(waco_s2flm$flag==1)]
# changed to 0.0005 degree radius (45 to 55 m at this latitude)
flagneighbors_which<-c()
for(i in 1:length(flagneighbors_lat)){
  whichi<-which(waco_s2flm$lat_dec<(flagneighbors_lat[i]+0.0005) & 
                  
                  waco_s2flm$lat_dec>(flagneighbors_lat[i]-0.0005) & 
                  
                  waco_s2flm$lon_dec<(flagneighbors_lon[i]+0.0005) & 
                  
                  waco_s2flm$lon_dec>(flagneighbors_lon[i]-0.0005))
  
  flagneighbors_which<-c(flagneighbors_which,whichi)
}  

# add column for keeping track of neighbors of flagged pixels that weren't themselves flagged (yet)
waco_s2flm$flagneighbor<-0
waco_s2flm$flagneighbor[unique(flagneighbors_which)]<-1
waco_s2flm$flagneighbor[which(waco_s2flm$flag==1)]<-0

# now flag pixels that neighbor a flagged pixel   
waco_s2flm$flag[unique(flagneighbors_which)]<-1

waco_unique<-waco_s2flm %>% distinct(lat,lon, .keep_all=TRUE)

waco_data_use<-waco_unique %>% filter(flag==0) %>% filter(SCL==6)

write.csv(waco_s2flm, "wpts_all.csv")
write.csv(waco_unique, "wpts_unq.csv")
write.csv(waco_data_use, "wpts_nf_unq.csv")


waco_plot<-ggplot(data_use,aes(log(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()

ggplotly(waco_plot)


### Lake Brownwood - 8/6/22
# first filtering to time on the lake
bw<-data_merged %>% filter(date=="2022-08-06", hour_dec>12.91, hour_dec<15.283)
# assigning unique id label 
bw_id <- bw %>%
  mutate(label = row_number()) 
bw_id$system <- "brownwood"

## writing out lat/lon for all dfs
bw_pts <- bw_id %>% dplyr::select(label,lat_dec,lon_dec) %>% na.omit()
# Split the full df into smaller subsets with <=5000 points each
bw_subsets <- split(bw_pts, rep(1:ceiling(nrow(bw_pts)/5000), each=5000, length.out=nrow(bw_pts)))
# Export each subset as a separate CSV file
for (i in seq_along(bw_subsets)) {
  subset <- bw_subsets[[i]]
  write.csv(subset, file = paste0("bwptsgee_", i, ".csv"), row.names = FALSE)
}

# import csv results from eepoints
bw0 <- read_csv(c("s2bvs_bw1.csv","s2bvs_bw2.csv")) %>% dplyr::arrange(label)

#bwee0 <- read_csv(c("bw1_eepts.csv","bw2_eepts.csv","bw3_eepts.csv","bw4_eepts.csv","bw5_eepts.csv")) %>% dplyr::arrange(label)

# filtering just to bands that will be used and calculating normalized diffs

bw_nd <- bw0 %>% dplyr::select("label", "lat_dec", "lon_dec", "B2", "B3", "B4", "B5", "SCL") %>% 
  mutate(ndti = (B4-B3)/(B4+B3),
         ndci = (B5-B4)/(B5+B4))


bw_s2flm <- merge(bw_id, bw_nd, by = "label") %>% dplyr::select(label, system, lat = lat_dec.x, lon = lon_dec.x,
                                                                datetime_utc, date_time_c, speed,fdom, wtemp,
                                                                chl_rfu, chl_ugl, turb,B2, B3, B4, B5, SCL, ndti, ndci) %>% dplyr::filter(turb<40)
  #dplyr::filter(speed>0) #%>% dplyr::filter(SCL==6)

#bwee <- distinct(bwee0, label, lat, lon, distkm, variable, .keep_all = TRUE) %>% 
 # filter(variable == "B2" | variable == "B3" | variable == "B4" | variable == "B5" | variable == "SCL")
#bw_pw <- bwee %>% 
 # pivot_wider(id_cols = c(label,lat, lon, distkm), names_from = variable, values_from = value) %>% 
  #mutate(ndti = (B4-B3)/(B4+B3),
   #      ndci = (B5-B4)/(B5+B4)) %>% 
  #dplyr::filter(SCL == 6)


#bwee_s2flm <- merge(bw_id, bw_pw, by = "label") %>% dplyr::select(label, lat, lon,datetime_utc, date_time_c, chl, turb,B2, B3, B4, B5, SCL, ndti, ndci)


bwlm <- lm(ndti~log(turb), bw_s2flm)
summary(bwlm)

bweelm <- lm(ndti~turb, bwee_s2flm)
summary(bweelm)





# adding flag for non-water pixels and filtering them out
bw_s2flm$flag<-0
bw_s2flm$flag[which(!bw_s2flm $SCL %in% c(4,5,6,11))]<-1
#bw_s2flm$flag[which(bw_s2flm $B2<=500)]<-1

# neighbor flag
flagneighbors_lat<-bw_s2flm$lat[which(bw_s2flm$flag==1)]
flagneighbors_lon<-bw_s2flm$lon[which(bw_s2flm$flag==1)]
# changed to 0.0005 degree radius (45 to 55 m at this latitude)
flagneighbors_which<-c()
for(i in 1:length(flagneighbors_lat)){
  whichi<-which(bw_s2flm$lat_dec<(flagneighbors_lat[i]+0.0005) & 
                  
                  bw_s2flm$lat_dec>(flagneighbors_lat[i]-0.0005) & 
                  
                  bw_s2flm$lon_dec<(flagneighbors_lon[i]+0.0005) & 
                  
                  bw_s2flm$lon_dec>(flagneighbors_lon[i]-0.0005))
  
  flagneighbors_which<-c(flagneighbors_which,whichi)
}  

# add column for keeping track of neighbors of flagged pixels that weren't themselves flagged (yet)
bw_s2flm$flagneighbor<-0
bw_s2flm$flagneighbor[unique(flagneighbors_which)]<-1
bw_s2flm$flagneighbor[which(bw_s2flm$flag==1)]<-0

# now flag pixels that neighbor a flagged pixel   
bw_s2flm$flag[unique(flagneighbors_which)]<-1

bw_unique<-bw_s2flm %>% distinct(lat,lon, .keep_all=TRUE)

#bw_data_use<-bw_s2flm %>% filter(flag==0) 

data_bw<-bw_unique %>% filter(flag==0) %>% filter(turb < 40)
data_bw6 <- bw_s2flm %>% filter(SCL == 6) %>% filter(turb < 40)

write.csv(bw_s2flm, "bwpts_all.csv")
write.csv(bw_unique, "bwpts_unq.csv")
write.csv(data_bw, "bwpts_nf_unq.csv")

# running linear model between ndti~turb
logEstbw <- lm(ndti~log(turb),data=data_bw)
xvec <- seq(min(data_bw$turb),max(data_bw$turb),length=1000)
logpred <- predict(logEstbw,newdata=data.frame(turb=xvec))
pred <- data.frame(x = xvec, y = logpred)
summary(logEstbw)

turb_ndti_plot<-ggplot(data_bw,aes(x=turb,y=ndti)) +
  geom_point(size=0.1,col=gray(0.3))+
  geom_smooth(method="lm",lty=5,col="black",se=FALSE,lwd=1)+
  geom_line(data = pred, aes(x=x, y=y),lty=1,col="black",lwd=1)+
  xlim(c(0,25))+
  xlab("In-lake Turbidity (NTU)")+
  ylab("Satellite Turbidity (NDTI)")+
  theme_bw()
turb_ndti_plot

ggplot(bw_s2flm,aes(turb, ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()





### Lake O.H. Ivie - 8/7/22
# first filtering to time on the lake
iv<-data_merged %>% filter(date=="2022-08-07", hour_dec>10.67, hour_dec<13.68)
# assigning unique id label 
iv_id <- iv%>%
  mutate(label = row_number()) 
iv_id$system <- "ivie"
## writing out lat/lon for all dfs
iv_pts <- iv_id %>% dplyr::select(label,lat_dec,lon_dec) %>% na.omit()
# Split the full df into smaller subsets with <=5000 points each
iv_subsets <- split(iv_pts, rep(1:ceiling(nrow(iv_pts)/5000), each=5000, length.out=nrow(iv_pts)))
# Export each subset as a separate CSV file
for (i in seq_along(iv_subsets)) {
  subset <- iv_subsets[[i]]
  write.csv(subset, file = paste0("ivpts_", i, ".csv"), row.names = FALSE)
}

# import csv results from eepoints
iv0 <- read_csv(c("s2bvs_iv1.csv","s2bvs_iv2.csv","s2bvs_iv3.csv")) %>% dplyr::arrange(label)

iv_nd <- iv0 %>% dplyr::select("label", "lat_dec", "lon_dec", "B2", "B3", "B4", "B5", "SCL") %>% 
  mutate(ndti = (B4-B3)/(B4+B3), ndci = (B5-B4)/(B5+B4))


iv_s2flm <- merge(iv_id, iv_nd, by = "label") %>% dplyr::select(label, system, lat = lat_dec.x, lon = lon_dec.x,
                                                                datetime_utc, date_time_c, speed,chl_rfu, chl_ugl, turb,fdom, wtemp,
                                                                B2, B3, B4, B5, SCL, ndti, ndci) %>% 
  dplyr::filter(speed>0) #%>% dplyr::filter(SCL==6)

ivlm<-lm(ndti~log(turb), data_use)

# adding flag for non-water pixels and filtering them out
iv_s2flm$flag<-0
iv_s2flm$flag[which(!iv_s2flm $SCL %in% c(4,5,6,11))]<-1
#iv_s2flm$flag[which(iv_s2flm $B2<=500)]<-1

# neighbor flag
flagneighbors_lat<-iv_s2flm$lat[which(iv_s2flm$flag==1)]
flagneighbors_lon<-iv_s2flm$lon[which(iv_s2flm$flag==1)]
# changed to 0.0005 degree radius (45 to 55 m at this latitude)
flagneighbors_which<-c()
for(i in 1:length(flagneighbors_lat)){
  whichi<-which(iv_s2flm$lat_dec<(flagneighbors_lat[i]+0.0005) & 
                  
                  iv_s2flm$lat_dec>(flagneighbors_lat[i]-0.0005) & 
                  
                  iv_s2flm$lon_dec<(flagneighbors_lon[i]+0.0005) & 
                  
                  iv_s2flm$lon_dec>(flagneighbors_lon[i]-0.0005))
  
  flagneighbors_which<-c(flagneighbors_which,whichi)
}  

# add column for keeping track of neighbors of flagged pixels that weren't themselves flagged (yet)
iv_s2flm$flagneighbor<-0
iv_s2flm$flagneighbor[unique(flagneighbors_which)]<-1
iv_s2flm$flagneighbor[which(iv_s2flm$flag==1)]<-0

# now flag pixels that neighbor a flagged pixel   
iv_s2flm$flag[unique(flagneighbors_which)]<-1

iv_unique<-iv_s2flm %>% distinct(lat,lon, .keep_all=TRUE)
iv_nf<-iv_unique %>% filter(flag==0) 

write.csv(iv_s2flm, "ivpts_all.csv")
write.csv(iv_unique, "ivpts_unq.csv")
write.csv(iv_nf, "ivpts_nf_unq.csv")

### Red Bluff Reservoir - 8/8/22
# first filtering to time on the lake
rb<-data_merged %>% filter(date=="2022-08-08", hour_dec>10.7, hour_dec<15.967)
# assigning unique id label 
rb_id <- rb %>%
  mutate(label = row_number()) 
rb_id$system <- "redbluff"
## writing out lat/lon for all dfs
rb_pts <- rb_id %>% dplyr::select(label,lat_dec,lon_dec) %>% na.omit()
# Split the full df into smaller subsets with <=5000 points each
rb_subsets <- split(rb_pts, rep(1:ceiling(nrow(rb_pts)/5000), each=5000, length.out=nrow(rb_pts)))
# Export each subset as a separate CSV file
for (i in seq_along(rb_subsets)) {
  subset <- rb_subsets[[i]]
  write.csv(subset, file = paste0("rbpts_", i, ".csv"), row.names = FALSE)
}

rb0 <- read_csv(c("s2bvs_rb1.csv","s2bvs_rb2.csv")) %>% dplyr::arrange(label)
rb_nd <- rb0 %>% dplyr::select("label", "lat_dec", "lon_dec", "B2", "B3", "B4", "B5", "SCL") %>% 
  mutate(ndti = (B4-B3)/(B4+B3), ndci = (B5-B4)/(B5+B4))


rb_s2flm <- merge(rb_id, rb_nd, by = "label") %>% dplyr::select(label, system, lat = lat_dec.x, lon = lon_dec.x,
                                                                datetime_utc, date_time_c, speed,chl_rfu, chl_ugl, turb,fdom, wtemp,
                                                                B2, B3, B4, B5, SCL, ndti, ndci) %>% 
  dplyr::filter(speed>0)

# adding flag for non-water pixels and filtering them out
rb_s2flm$flag<-0
rb_s2flm$flag[which(!rb_s2flm $SCL %in% c(4,5,6,11))]<-1
#rb_s2flm$flag[which(rb_s2flm $B2<=500)]<-1

# neighbor flag
flagneighbors_lat<-rb_s2flm$lat[which(rb_s2flm$flag==1)]
flagneighbors_lon<-rb_s2flm$lon[which(rb_s2flm$flag==1)]
# changed to 0.0005 degree radius (45 to 55 m at this latitude)
flagneighbors_which<-c()
for(i in 1:length(flagneighbors_lat)){
  whichi<-which(rb_s2flm$lat_dec<(flagneighbors_lat[i]+0.0005) & 
                  
                  rb_s2flm$lat_dec>(flagneighbors_lat[i]-0.0005) & 
                  
                  rb_s2flm$lon_dec<(flagneighbors_lon[i]+0.0005) & 
                  
                  rb_s2flm$lon_dec>(flagneighbors_lon[i]-0.0005))
  
  flagneighbors_which<-c(flagneighbors_which,whichi)
}  

# add column for keeping track of neighbors of flagged pixels that weren't themselves flagged (yet)
rb_s2flm$flagneighbor<-0
rb_s2flm$flagneighbor[unique(flagneighbors_which)]<-1
rb_s2flm$flagneighbor[which(rb_s2flm$flag==1)]<-0

# now flag pixels that neighbor a flagged pixel   
rb_s2flm$flag[unique(flagneighbors_which)]<-1
data_use<-rb_s2flm %>% filter(flag==0) 

rb_unique<-rb_s2flm %>% distinct(lat,lon, .keep_all=TRUE)
rb_nf<-rb_unique %>% filter(flag==0) 

write.csv(rb_s2flm, "rbpts_all.csv")
write.csv(rb_unique, "rbpts_unq.csv")
write.csv(rb_nf, "rbpts_nf_unq.csv")


######## binding all s2flm  dfs to viz whole enchilada
# only ah_s2flm has been filtered for water only pixels
# all the rest have flag for non SCL 6 and neighbors but not for B2<500
## data_use<-s2flm %>% filter(flag==0) to filter flag
# filter(speed>0) has not been applied to any all _s2flm are just merged s2&flm for the lake

s2flm_list <- list(ah_s2flm, bon_s2flm, waco_s2flm, bw_s2flm, iv_s2flm, rb_s2flm)

#merge all data frames in list
s2flm_df <- Reduce(function(x, y) merge(x, y, all=TRUE), s2flm_list)

data_use<-s2flm_df %>% filter(flag==0) 
s2flm_nf_unique<-data_use%>% distinct(lat,lon, .keep_all=TRUE)

write_csv(s2flm_df, "s2flm6lakes.csv")
write_csv(data_use, "s2flm6lakes_noflag.csv")

nfuq_list <- list(ah_unique, db_unq, waco_data_use, data_bw, iv_nf, rb_nf)

#merge all data frames in list
nfuq_df <- Reduce(function(x, y) merge(x, y, all=TRUE), nfuq_list)

write_csv(nfuq_df, "s2flm6lakes_nf_unique.csv")



#################################################################################################
# calculating DWL, and TSI as dist from dam
# reading in points selected for a dist from dam transect on a clear day
waco_sarm<-read_csv("waco_sarm25july.csv")%>% distinct(label, lat, lon, variable, distkm, .keep_all = T)
waco_sa_pw<-waco_sarm %>% 
  pivot_wider(id_cols = c(label,lat, lon, distkm), names_from = variable, values_from = value) %>% 
  mutate(ndti = (B4-B3)/(B4+B3)) %>% filter(SCL == 6)
# normalizing the distance from 0-1
wsa_min<-min(waco_sa_pw$distkm)
wsa_max<-max(waco_sa_pw$distkm)
wsa<-waco_sa_pw %>% mutate(norm_dist = (distkm - wsa_min)/(wsa_max - wsa_min))
wsa$system<-"waco"

# arrowhead
ah_28july<-read_csv("ah_28july.csv")%>% distinct(label, lat, lon, variable, distkm, .keep_all = T)
ah_28july_pw<-ah_28july %>% 
  pivot_wider(id_cols = c(label,lat, lon, distkm), names_from = variable, values_from = value) %>% 
  mutate(ndti = (B4-B3)/(B4+B3)) %>% filter(SCL == 6)
# normalizing the distance from 0-1
ah_min<-min(ah_28july_pw$distkm)
ah_max<-max(ah_28july_pw$distkm)
ahpw<-ah_28july_pw %>% mutate(norm_dist = (distkm - ah_min)/(ah_max - ah_min))
ahpw$system<-"arrowhead"

# brownwood
bw_sarm<-read_csv("brownwood_sarm6aug.csv")%>% distinct(label, lat, lon, variable, distkm, .keep_all = T)
bw_sa_pw<-bw_sarm %>% 
  pivot_wider(id_cols = c(label,lat, lon, distkm), names_from = variable, values_from = value) %>% 
  mutate(ndti = (B4-B3)/(B4+B3)) %>% filter(SCL == 6)
# normalizing the distance from 0-1
bw_min<-min(bw_sa_pw$distkm)
bw_max<-max(bw_sa_pw$distkm)
bwsa<-bw_sa_pw %>% mutate(norm_dist = (distkm - bw_min)/(bw_max - bw_min))
bwsa$system<-"brownwood"

# ivie
iv_7aug<-read_csv("ohivie_7aug.csv")%>% distinct(label, lat, lon, variable, distkm, .keep_all = T)
iv_7aug_pw<-iv_7aug %>% 
  pivot_wider(id_cols = c(label,lat, lon, distkm), names_from = variable, values_from = value) %>% 
  mutate(ndti = (B4-B3)/(B4+B3)) %>% filter(SCL == 6) %>% filter(ndti<0.08) 
iv_7aug_pw<-iv_7aug_pw[-c(36,37), ]
# normalizing the distance from 0-1
iv_min<-min(iv_7aug_pw$distkm)
iv_max<-max(iv_7aug_pw$distkm)
ivpw<-iv_7aug_pw %>% mutate(norm_dist = (distkm - iv_min)/(iv_max - iv_min))
ivpw$system<-"ohivie"

#iv_28july<-read_csv("ivie_28july.csv")%>% distinct(label, lat, lon, variable, distkm, .keep_all = T)
#iv_28july_pw<-iv_28july %>% 
 # pivot_wider(id_cols = c(label,lat, lon, distkm), names_from = variable, values_from = value) %>% 
#  mutate(ndti = (B4-B3)/(B4+B3)) %>% filter(SCL == 6) #%>% filter(ndti<0.08)
# normalizing the distance from 0-1
#iv_min<-min(iv_28july_pw$distkm)
#iv_max<-max(iv_28july_pw$distkm)
#ivpw<-iv_28july_pw %>% mutate(norm_dist = (distkm - iv_min)/(iv_max - iv_min))
#ivpw$system<-"ohivie"

# bonham
bon_15jul<-read_csv("bonham_15jul.csv")%>% distinct(label, lat, lon, variable, distkm, .keep_all = T)
bon_15jul_pw<-bon_15jul%>% 
  pivot_wider(id_cols = c(label,lat, lon, distkm), names_from = variable, values_from = value) %>% 
  mutate(ndti = (B4-B3)/(B4+B3)) %>% filter(SCL == 6)
# normalizing the distance from 0-1
bn_min<-min(bon_15jul_pw$distkm)
bn_max<-max(bon_15jul_pw$distkm)
bnpw<-bon_15jul_pw %>% mutate(norm_dist = (distkm - bn_min)/(bn_max - bn_min))
bnpw$system<-"bonham"

# redbluff
rb_8aug<-read_csv("rb_8aug.csv")%>% distinct(label, lat, lon, variable, distkm, .keep_all = T)
rb_8aug_pw<-rb_8aug%>% 
  pivot_wider(id_cols = c(label,lat, lon, distkm), names_from = variable, values_from = value) %>% 
  mutate(ndti = (B4-B3)/(B4+B3)) %>% filter(SCL == 6)
# normalizing the distance from 0-1
rb_min<-min(rb_8aug_pw$distkm)
rb_max<-max(rb_8aug_pw$distkm)
rbpw<-rb_8aug_pw %>% mutate(norm_dist = (distkm - rb_min)/(rb_max - rb_min))
rbpw$system<-"redbluff"

#merge all data frames in list
s2t_list <- list(ahpw, bnpw, wsa, bwsa, ivpw, rbpw)

s2t_df <- Reduce(function(x, y) merge(x, y, all=TRUE), s2t_list)
dfd_transect<-write_csv(s2t_df, "dfd_transect.csv")
# setting factor of "system" to reflect date
s2t_df$system<- factor(s2t_df$system, levels = c("bonham", "waco", "brownwood","ohivie","redbluff", "arrowhead"))


line_ndti<- s2t_df %>% ggplot(aes(norm_dist, ndti, color = system)) + geom_line(size=1.25) + 
  xlab("Normalized Distance from Dam") + ylab("Normalized Difference Turbidity Index") +  
  theme_minimal()+scale_color_manual(values = rev(lacroix_palette("PeachPear", type = "discrete")),
                                     labels=c("bonham" = "Bonham", "waco"="Waco","brownwood"="Brownwood",
                                              "ohivie"="O.H. Ivie","redbluff"="Red Bluff","arrowhead"="Arrowhead"),
                                     name = "System")+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 14))

ggsave("line_ndti.png",line_ndti)

# doing the same but in viridis - plasma
line_ndti_p<- s2t_df %>% ggplot(aes(norm_dist, ndti, color = system)) + geom_line(size=1.25) + 
  xlab("Normalized Distance from Dam") + ylab("Normalized Difference Turbidity Index") +  
  theme_bw()+scale_color_manual(values = viridis::viridis(6, option = "C"),
                                     labels=c("bonham" = "Bonham", "waco"="Waco","brownwood"="Brownwood",
                                              "ohivie"="O.H. Ivie","redbluff"="Red Bluff","arrowhead"="Arrowhead"),
                                     name = "System")+
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))

ggsave("line_ndti_plasma.png",line_ndti_p,width = 7.5, height = 5)

# adding chromaticity DWL values to each point
#Convert R,G, and B spectral reflectance to dwl based on CIE chromaticity color space, see Wang et al 2015.
fui.hue <- function(R, G, B) {
  
  require(colorscience)
  # chromaticity.diagram.color.fill()
  Xi <- 2.7689*R + 1.7517*G + 1.1302*B
  Yi <- 1.0000*R + 4.5907*G + 0.0601*B
  Zi <- 0.0565*G + 5.5943*B
  
  # calculate coordinates on chromaticity diagram
  x <-  Xi / (Xi + Yi +  Zi)
  y <-  Yi / (Xi + Yi +  Zi)
  z <-  Zi / (Xi + Yi +  Zi)
  
  # calculate hue angle
  alpha <- atan2( (x - 0.33), (y - 0.33)) * 180/pi
  
  # make look up table for hue angle to wavelength conversion
  cie <- cccie31 %>%
    mutate(a = atan2( (x - 0.33), (y - 0.33)) * 180/pi) %>%
    dplyr::filter(wlnm <= 700) %>%
    dplyr::filter(wlnm >=380)
  
  # find nearest dominant wavelength to hue angle
  wl <- cie[as.vector(sapply(alpha,function(x) which.min(abs(x - cie$a)))), 'wlnm']
  #out <- cbind(as.data.frame(alpha), as.data.frame(wl))
  
  return(wl)
}


# assigning bands to R,G,B
R <- s2t_df$B4
G <- s2t_df$B3
B <- s2t_df$B2

s2t_df <- s2t_df %>% mutate(dwl = fui.hue(R, G, B))


# now seeing DWL as normDFD
line_dwl<- s2t_df %>% ggplot(aes(norm_dist, dwl, color = system)) + geom_line(size=1.25) + 
  xlab("Normalized Distance from Dam") + ylab("Dominant Wavelength") + theme_minimal() +
  scale_color_manual(values = rev(lacroix_palette("PeachPear", type = "discrete")),
                     labels=c("bonham" = "Bonham", "waco"="Waco","brownwood"="Brownwood","ohivie"="O.H. Ivie",
                              "redbluff"="Red Bluff","arrowhead"="Arrowhead"),name = "System")+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 14))

ggsave("line_dwl.png", line_dwl)

# doing the same in viridis - plasma
line_dwl_p<- s2t_df %>% ggplot(aes(norm_dist, dwl, color = system)) + geom_line(size=1.25) + 
  xlab("Normalized Distance from Dam") + ylab("Dominant Wavelength") + theme_bw() +
  scale_color_manual(values = viridis::viridis(6, option = "C"),
                     labels=c("bonham" = "Bonham", "waco"="Waco","brownwood"="Brownwood",
                              "ohivie"="O.H. Ivie","redbluff"="Red Bluff","arrowhead"="Arrowhead"),
                     name = "System")+
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))

ggsave("line_dwl_plasma.png",line_dwl_p,width = 7.5, height = 5)

comb_linedfd_plot<- line_ndti_p + line_dwl_p + plot_layout(guides = "collect") & theme(legend.position = "right")
ggsave("combo_line_plasma.png",comb_linedfd_plot,width = 12, height = 5)
# Display the combined plot
print(combined_plot)
# now adding AVW with 10 bands
s2t_df$avw10band <- ((s2t_df$B2 + s2t_df$B3 + s2t_df$B4 + s2t_df$B5 + s2t_df$B6 + s2t_df$B7 + s2t_df$B8 + s2t_df$B8A + s2t_df$B11 + s2t_df$B12) / 
                       ((s2t_df$B2/496.6) + (s2t_df$B3/560) + (s2t_df$B4/664.5) + (s2t_df$B5/703.9)+ (s2t_df$B6/740.2)+ (s2t_df$B7/782.5)+ (s2t_df$B8/835.1)+ 
                          (s2t_df$B8A/864.8)+ (s2t_df$B11/1613.7)+ (s2t_df$B12/2202.4)))

line_avw<- s2t_df %>% ggplot(aes(norm_dist, avw4band, color = system)) + geom_line(size=1.25) + 
  xlab("Normalized Distance from Dam") + ylab("Dominant Wavelength") + theme_minimal()
line_avw + scale_color_manual(values = lacroix_palette("PeachPear", type = "discrete"))

## trying SP's AVW way
sums_avw<-s2t_df$B2+s2t_df$B3+s2t_df$B4
s2t_df$spc_nm<-(490*s2t_df$B2/sums_avw)+(560*s2t_df$B3/sums_avw)+(665*s2t_df$B4/sums_avw)

## AVW RBG
s2t_df$avw_rbg <- ((s2t_df$B2 + s2t_df$B3 + s2t_df$B4) / 
                       ((s2t_df$B2/496.6) + (s2t_df$B3/560) + (s2t_df$B4/664.5)))

########################################################################################################################
## testing out with all 6 lakes
all6lakes<-read_csv("s2flm6lakes.csv")
noflaglakes <- read_csv("s2flm6lakes_noflag.csv")

s2_df<-read_csv("sixres_noflag.csv")

nfunq<-read_csv("s2flm6lakes_nf_unique.csv")
# adding dwl to all6lakes
# assigning bands to R,G,B
R <- all6lakes$B4
G <- all6lakes$B3
B <- all6lakes$B2

s2t_df <- all6lakes %>% mutate(dwl = fui.hue(R, G, B)) 

R <- nfunq$B4
G <- nfunq$B3
B <- nfunq$B2

s2nfu<-nfunq%>% mutate(dwl = fui.hue(R, G, B)) 
#waco723<-all6lakes %>% filter(system=="waco") %>% filter(flag==0) %>% filter(turb<20) %>% filter(speed>0)

#ggplot(waco723,aes(log(turb), ndti)) + 
 # geom_point() + 
  #geom_smooth(method = "lm", se=FALSE) +
  #stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()

s2flm_noflag <- s2t_df%>% filter(speed>0) %>% filter(flag == 0)%>% filter(turb<100 & ndti < 0.3) %>% filter(SCL == 6)
s2flm_unique<-s2flm_noflag %>% distinct(lat,lon, .keep_all=TRUE)
write_csv(s2flm_unique, "s2flm_6lakes_datause.csv")

s2flm_5lake <- s2flm_noflag %>% filter(!system == "waco") %>% filter(turb<100 & ndti < 0.3) %>% filter(speed>0) #%>% filter(chl_rfu>0)

s2flm_nobon <- s2flm_noflag %>% filter(!system == "bonham")

ggplot(s2flm_unique,aes(log(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()


lm5lake<-lm(ndti~turb,data=all6lakes)

# curating full transects to write out
s2flm_select<-s2flm_noflag %>% dplyr::select(lat, lon, system, datetime_cst=date_time_c,turb,ndti,dwl) %>% dplyr::arrange(datetime_cst)
write_csv(s2flm_select, "sixres_noflag.csv")

s2flm_du_select<-s2flm_unique %>% dplyr::select(lat, lon, system, datetime_cst=date_time_c,turb,ndti,dwl) %>% dplyr::arrange(datetime_cst)
write_csv(s2flm_du_select, "sixres_select.csv")


#####################################################################################################################################
# separating lake data between river arm and main body  
s2_df<-read_csv("s2flm_6lakes_datause.csv")

#dm_sf <- st_as_sf(bnras, coords = c("lon", "lat"), crs = 4326)
#mapview(dm_sf)

s2_w<-s2_df %>% filter(system == "waco")
s2_a<-s2_df %>% filter(system == "arrowhead")
s2_bw<-s2_df %>% filter(system == "brownwood")
s2_iv<-s2_df %>% filter(system == "ivie")
s2_bn<-s2_df %>% filter(system == "bonham")
s2_rb<-s2_df %>% filter(system == "redbluff")

## defining bboxes for waco ##
#wmbbox <- st_bbox(c(xmin = -97.25316, ymin = 31.539846, xmax = -97.1882756, ymax = 31.589283782), crs = st_crs(s2_w))

wmb <- s2_w[s2_w$lat >= 31.530244 & s2_w$lat <= 31.6095146 & s2_w$lon >= -97.2520428 & s2_w$lon <= -97.1809749,]
wna<-s2_w[s2_w$lat >= 31.5821285 & s2_w$lat <= 31.6157569 & s2_w$lon >= -97.308004 & s2_w$lon <= -97.2523861,]
wsa<-s2_w[s2_w$lat >= 31.4886804 & s2_w$lat <= 31.5311218 & s2_w$lon >= -97.265776 & s2_w$lon <= -97.2266369,]

wmb$zone="body"
wna$zone="arm"
wsa$zone="arm"

waco_bind<-rbind(wmb,wna,wsa)
#write_csv(waco_bind, "waco_zone.csv")

## defining bboxes for ah ##
amb<-s2_a[s2_a$lat >= 33.6933513 & s2_a$lat <= 33.7744378 & s2_a$lon >= -98.41832729 & s2_a$lon <= -98.3015976,]
ara<-s2_a[s2_a$lat >= 33.6339153 & s2_a$lat <= 33.6936369 & s2_a$lon >= -98.4667358 & s2_a$lon <= -98.351036,]

amb$zone="body"
ara$zone="arm"

ah_bind<-rbind(amb,ara)

## defining bboxes for bw ##
bwmb<-s2_bw[s2_bw$lat >= 31.7955865 & s2_bw$lat <= 31.865884 & s2_bw$lon >= -99.0787206 & s2_bw$lon <= -98.9753804,]
bwna<-s2_bw[s2_bw$lat >= 31.8661754 & s2_bw$lat <= 31.928552 & s2_bw$lon >= -99.0583197 & s2_bw$lon <= -99.0109412,]
bwsa<-s2_bw[s2_bw$lat >= 31.8200945 & s2_bw$lat <= 31.857136 & s2_bw$lon >= -99.1343388 & s2_bw$lon <= -99.0787206,]

bwmb$zone="body"
bwna$zone="arm"
bwsa$zone="arm"

bw_bind<-rbind(bwmb,bwna,bwsa)

## defining bboxes for iv ##
imb<-s2_iv[s2_iv$lat >= 31.479319 & s2_iv$lat <= 31.5557074 & s2_iv$lon >= -99.7077124 & s2_iv$lon <= -99.574503,]
ira<-s2_iv[s2_iv$lat >= 31.5559999 & s2_iv$lat <= 31.6256021 & s2_iv$lon >= -99.8392049 & s2_iv$lon <= -99.659991,]

imb$zone="body"
ira$zone="arm"

iv_bind<-rbind(imb,ira)

## defining bboxes for bon ##
bnmb1<-s2_bn[s2_bn$lat >= 33.6512746 & s2_bn$lat <= 33.6556328 & s2_bn$lon >= -96.1523812 & s2_bn$lon <= -96.13573,]
bnmb2<-s2_bn[s2_bn$lat >= 33.6442593 & s2_bn$lat <= 33.6515473 & s2_bn$lon >= -96.1636449 & s2_bn$lon <= -96.1349774,]
bnmb1$zone="body"
bnmb2$zone="body"


bnrane<-s2_bn[s2_bn$lat >= 33.6556197 & s2_bn$lat <= 33.67240739 & s2_bn$lon >= -96.15240104 & s2_bn$lon <= -96.13789,]
bnranw<-s2_bn[s2_bn$lat >= 33.6516187 & s2_bn$lat <= 33.66576413 & s2_bn$lon >= -96.1694813 & s2_bn$lon <= -96.152401,]
bnraw<-s2_bn[s2_bn$lat >= 33.63775681 & s2_bn$lat <= 33.65047555 & s2_bn$lon >= -96.17806441 & s2_bn$lon <= -96.163731,] 
bnras<-s2_bn[s2_bn$lat >= 33.6397708 & s2_bn$lat <= 33.644201 & s2_bn$lon >= -96.1537165 & s2_bn$lon <= -96.13637867,]
bnrane$zone="arm"
bnranw$zone="arm"
bnraw$zone="arm"
bnras$zone="arm"

bn_bind<-rbind(bnmb1,bnmb2,bnrane,bnranw,bnraw,bnras)

## defining bboxes for rb ##
rmb<-s2_rb[s2_rb$lat >= 31.895208 & s2_rb$lat <= 31.9744585 & s2_rb$lon >= -103.9684497 & s2_rb$lon <= -103.8946353,]
rra<-s2_rb[s2_rb$lat >= 31.974749 & s2_rb$lat <= 32.0443283 & s2_rb$lon >= -104.01239497 & s2_rb$lon <= -103.9324008,]

rmb$zone="body"
rra$zone="arm"

rb_bind<-rbind(rmb,rra)


all_bind <- rbind(bn_bind,waco_bind,ah_bind,bw_bind,iv_bind,rb_bind)
#write_csv(all_bind,"sixres_zone.csv")

ndti_aov<-aov(ndti ~ zone + system, all_bind)
summary(ndti_aov)

plot(ndti_aov)

TukeyHSD(ndti_aov)


##############################################################################################
# adding NLA derived TSI to each flame point
# reading in data (all NLA reservoirs for 2012 + 2017 and CRASR sample data)
nla_crasr<-read_csv("nla_crasr.csv")
nla_crasr <- nla_crasr%>% dplyr::rename(secchi = `secchi (m)`, turb = `turb (ntu)`)
# lm for log(secchi)~log*turb
tsdlm<-lm(log(secchi)~log(turb), data = nla_crasr)
# predicting secchi values from EXO turbidity data
all_bind$secchi <- exp(predict(tsdlm, all_bind))
# calculating TSI from newly derived secchi
all_bind$tsi_sd <- 10 * (6-log(all_bind$secchi)/log(2))

# writing zoned csv with all variables
write_csv(all_bind,"six4m_zone.csv")

sd_aov<-aov(secchi ~ zone + system + zone:system, all_bind)
summary(sd_aov)
plot(sd_aov)

TukeyHSD(sd_aov)

# aov for turb
turb_aov<-aov(turb ~ zone + system+ zone:system, all_bind)
summary(turb_aov)

plot(turb_aov)

tukey_turb<-TukeyHSD(turb_aov)



# aov for dwl
dwl_aov<-aov(dwl ~ zone + system + zone:system, all_bind)
dwl_summary<-summary(dwl_aov)
plot(dwl_aov)

TukeyHSD(dwl_aov)


F_value <- dwl_summary[[1]]["zone", "F value"]
df1 <- dwl_summary[[1]]["zone", "Df"]
df2 <- dwl_summary[[1]]["Residuals", "Df"]

# Calculate the exact p-value
p_value <- pf(F_value, df1, df2, lower.tail = FALSE)
p_value

###################################################################################################
# making plots with "six4m_zone.csv"
zn6res_df<-read_csv("six4m_zone.csv")

ggplot(zn6res_df,aes(log(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()

# ndti as normDFD
line_ndti<-zn6res_df %>% ggplot(aes(norm_dist, ndti, color = system)) + geom_line(size=1.25) + 
  xlab("Normalized Distance from Dam") + ylab("Normalized Difference Turbidity Index") +  
  theme_minimal()+scale_color_viridis_d() +
                                     labels=c("bonham" = "Bonham", "waco"="Waco","brownwood"="Brownwood",
                                              "ohivie"="O.H. Ivie","redbluff"="Red Bluff","arrowhead"="Arrowhead"),name = "System"+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 14))

ggsave("line_ndti.png",line_ndti)


# now seeing DWL as normDFD
line_dwl<- zn6res_df %>% ggplot(aes(norm_dist, dwl, color = system)) + geom_line(size=1.25) + 
  xlab("Normalized Distance from Dam") + ylab("Dominant Wavelength") + theme_minimal() +
  scale_color_manual(values = rev(lacroix_palette("PeachPear", type = "discrete")),
                     labels=c("bonham" = "Bonham", "waco"="Waco","brownwood"="Brownwood","ohivie"="O.H. Ivie",
                              "redbluff"="Red Bluff","arrowhead"="Arrowhead"),name = "System")+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.title = element_text(size = 14))

ggsave("line_dwl.png", line_dwl)


######## viz box plots with zones ##############################################
ggboxplot(all_bind, x = "system", y = "secchi", fill = "zone",
          palette = c("#00AFBB", "#E7B800"))

p<- ggboxplot(all_bind, x = "system", y = "secchi", fill = "zone",
          palette = c("#E7B800","#00AFBB"), order = c("bonham", "waco","brownwood",
                                                       "ivie", "redbluff","arrowhead")) + 
  theme(legend.position = "right")

all_bind %>% ggplot(aes(x=system, y=secchi, fill =zone)) + 
  labs(fill = "Zone", labels = c("Arm", "Body")) +
  geom_boxplot() +theme_classic() + scale_fill_manual(values = c("#E7B800","#00AFBB"))

#p<-all_bind %>% ggplot(aes(system, secchi, fill = zone)) + geom_boxplot()


anova_results <- compare_means(secchi ~ zone, data = all_bind, group.by = "system", method = "anova")

# Add the significance annotations to the plot
p <- p + stat_pvalue_manual(anova_results, label = "p.signif", tip.length = 0.01)

# Move the legend position
p <- p + theme(legend.position = "bottom")

## visualizing turbidity to dwl relationship ##
gtz<-s2t_df %>% filter(flag == 0)%>% filter(turb<100 & ndti < 0.3) %>% filter(speed>0)

twl<-lm(log(turb)~dwl, data=all_bind)
summary(twl)

ggplot(all_bind,aes(log(turb), dwl)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()

# summary stats for the following: all6lakes, unique_nf, byzone
all6lakes %>% dplyr::group_by(system) %>% 
  select(speed, turb) %>%  summary()

s2t_df_select<-s2t_df %>% select(system, speed, turb, ndti, dwl)
s2df_by_summary<-by(s2t_df_select, s2t_df$system, summary)

zn_df_select<-zn6res_df %>% select(system, zone,speed, turb, ndti, dwl, secchi, tsi_sd)
zn_by_summary<-by(zn_df_select, zn6res_df$system, summary)

summary_df <- zn_df_select %>%
  group_by(system) %>%
  summarise(
    mean_secchi = mean(secchi),
    sd_secchi = sd(secchi),
    mean_turb = mean(turb),
    sd_turb = sd(turb),
    mean_dwl = mean(dwl),
    sd_dwl = sd(dwl))


# means, CI, and range for variables by zone
# defining confidence level
confidence_level <- 0.95
alpha <- 1-confidence_level

# summarizing
systemstats_turb<- zn_df_select %>% group_by(system) %>% 
  summarise(n=n(),mean=mean(turb),sd=sd(turb),se=sd/sqrt(n),
            lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
            upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
            min=min(turb), max=max(turb),
            margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))

zonestats_turb<- zn_df_select %>% group_by(zone) %>% 
  summarise(n=n(),mean=mean(turb),sd=sd(turb),se=sd/sqrt(n),
            lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
            upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
            min=min(turb), max=max(turb),
            margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))

systemstats_secchi<- zn_df_select %>% group_by(system) %>% 
  summarise(n=n(),mean=mean(secchi),sd=sd(secchi),se=sd/sqrt(n),
            lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
            upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
            min=min(secchi), max=max(secchi),
            margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))

zonestats_secchi<- zn_df_select %>% group_by(zone) %>% 
  summarise(n=n(),mean=mean(secchi),sd=sd(secchi),se=sd/sqrt(n),
            lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
            upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
            min=min(secchi), max=max(secchi),
            margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))

systemstats_dwl<- zn_df_select %>% group_by(system) %>% 
  summarise(n=n(),mean=mean(dwl),sd=sd(dwl),se=sd/sqrt(n),
            lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
            upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
            min=min(dwl), max=max(dwl),
            margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))

zonestats_dwl<- zn_df_select %>% group_by(zone) %>% 
  summarise(n=n(),mean=mean(dwl),sd=sd(dwl),se=sd/sqrt(n),
            lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
            upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
            min=min(dwl), max=max(dwl),
            margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))

zonestats_ndti<- zn_df_select %>% group_by(zone) %>% 
  summarise(n=n(),mean=mean(ndti),sd=sd(ndti),se=sd/sqrt(n),
            lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
            upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
            min=min(dwl), max=max(dwl),
            margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))
###############################################################################################################
# comparing predicted SD values with measurements

# CRASR data - read in for all lakes
crasr_all <- read_csv("CRASR_summer22all.csv") %>% clean_names()

# separate by lake day
#crasr_bon1 <- crasr_all %>% filter(datecollect == "7/14/2022") %>% add_column(lake = "bon") %>% 
 # group_by(system,lat, lon, variable, station=site)%>% mutate(meanchl = mean(value))
#cbon1 <- crasr_bon1 %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl)

crasr_bon2 <- crasr_all %>% filter(datecollect == "7/15/2022") %>% add_column(lake = "bon") %>% 
  dplyr::mutate(station = paste("B", site, sep = "")) %>% group_by(system,lat, lon, variable, site) %>% mutate(meanchl = mean(value))
cbon2 <- crasr_bon2 %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl)

crasr_lw <- crasr_all %>% filter(datecollect == "7/23/2022")%>% add_column(lake = "lw") %>% 
  dplyr::mutate(station = paste("LW", site, sep = "")) %>% group_by(system,lat, lon, variable, site) %>% mutate(meanchl = mean(value))
clw <- crasr_lw %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl)

crasr_bw <- crasr_all %>% filter(datecollect == "8/6/2022")%>% add_column(lake = "bw") %>% 
  group_by(system,lat, lon, variable, station=site) %>% mutate(meanchl = mean(value))
cbw <- crasr_bw %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl)

crasr_iv <- crasr_all %>% filter(datecollect == "8/7/2022") %>% add_column(lake = "iv") %>% 
  group_by(system,lat, lon, variable, station=site) %>% mutate(meanchl = mean(value))
civ <- crasr_iv %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl)

crasr_rb <- crasr_all %>% filter(datecollect == "8/8/2022") %>% add_column(lake = "rb") %>% 
  group_by(system, lat, lon, variable, station=site) %>% mutate(meanchl = mean(value))
crb <- crasr_rb %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl)

crasr_ah <- crasr_all %>% filter(datecollect == "8/12/2022") %>% add_column(lake = "ah") %>% 
  group_by(system, lat, lon, variable, station=site) %>% mutate(meanchl = mean(value))
cah <- crasr_ah %>% distinct(meanchl) %>% pivot_wider(names_from = variable, values_from = meanchl)

# combine crasr data
crasr_list <- list(cbon2, clw, cbw, civ,crb, cah)

#merge all data frames in list
crasr_df <- Reduce(function(x, y) merge(x, y, all=TRUE), crasr_list)

# adding a predicted sdd and tsi column
# predicting secchi values from EXO turbidity data
bw$sdd <- exp(predict(tsdlm, bw))
# calculating TSI from newly derived secchi
bw$tsi_sd <- 10 * (6-log(bw$sdd)/log(2))

samplesensorbw_merged<-merge(bw, cbw) 


samplesensorbw_merged$lat_diff<-samplesensorbw_merged$lat_dec-samplesensorbw_merged$lat
samplesensorbw_merged$lon_diff<-samplesensorbw_merged$lon_dec-samplesensorbw_merged$lon
samplesensorbw_merged<-samplesensorbw_merged %>% filter(abs(lat_diff)<0.0001,abs(lon_diff)<0.0001)

validationdatabw<-samplesensorbw_merged %>% dplyr::select(system,station,lat,lon,turb_lab = turbidity,
                                                          turb_ysi=turb,sdd_pred = sdd, secchi) %>% 
  group_by(station,lat,lon,turb_lab,turb_ysi,sdd_pred,secchi) %>%
  dplyr::summarize(turb_ysi=mean(turb_ysi,na.rm=TRUE)) %>%
  dplyr::summarize(sdd_pred=mean(sdd_pred,na.rm=TRUE)) %>% 
  as.data.frame()

valbw_means <- validationdatabw %>% group_by(station) %>%
  mutate(turb_ysi_m=mean(turb_ysi,na.rm=TRUE)) %>% 
  mutate(sdd_pred_m=mean(sdd_pred,na.rm=TRUE))

lmc_bw<-lm(sdd_pred_m~secchi,data=valbw_means) # r2 = 0.489
summary(lmc_bw)


# calculating DON + DOP per site


##################################################################################################################
# determining zone using a depth threshold
#echomap_depth<- read_csv("AllLogsThrough3May2023.gpx")

## function to shift vectors
shift.vec <- function (vec, shift) {
  if(length(vec) <= abs(shift)) {
    rep(NA ,length(vec))
  }else{
    if (shift >= 0) {
      c(rep(NA, shift), vec[1:(length(vec)-shift)]) }
    else {
      c(vec[(abs(shift)+1):length(vec)], rep(NA, abs(shift))) } } }

options(digits=10)
## reading in gpx track
track <- htmlTreeParse(file = "AllLogsThrough3May2023.gpx", error = function(...)
{ }, useInternalNodes = T)

## Getting elevations, times, and coordinates via respective xpath
depth <- as.numeric(xpathSApply(track, path = "//trkpt/Depth", xmlValue))
times <- xpathSApply(track, path = "//trkpt/time", xmlValue)
coords <- xpathSApply(track, path = "//trkpt", xmlAttrs)

## Extract lat and lon from coordinates
lats <- as.numeric(coords["lat",])
lons <- as.numeric(coords["lon",])

#gpx <- data.frame(lat = lats, lon  = lons, , time = times)
#rm(list = c("elev", "lats", "lons", "times", "track", "coords"))
#head(gpx)

# trying other way
# Load your GPX file
gpx_file <- "AllLogsThrough3May2023.gpx"
gpx <- xmlParse(gpx_file)
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





##########################################################################################################
# lake waco 6/23 JUNK DATA :'(
# Exo 2 data
getUTMzone <- function(lon) {
  (floor((lon + 180)/6) %% 60) + 1
}

table0<-read_xlsx(skip=8,"KorExo Measurement File Export 062322 - Lake Waco NLA 2022 flame 165134.xlsx")
table<-table0 %>% as.data.frame()
table_clean<-clean_names(table)
exo<-table_clean
exo<-exo %>% rename(exo_time_hh_mm_ss=time_hh_mm_ss,exo_date_mm_dd_yyyy=date_mm_dd_yyyy)

# Correct for lag time between intake and sensor using an approximation
# the Exo 2 flow cell has a volume of 1L
exo$lake_time_local<-substr(exo$exo_time_hh_mm_ss-15,12,19) # lag time correction in seconds 
#exo1$time_central<-substr(exo1$time_hh_mm_ss,12,19)

#exo1$laketime_central<-exo1$time_hh_mm_ss-10
names(exo)[which(names(exo)=="p_h")]<-"pH"

#exo<-exo %>% filter(lake_time_central>12)

# Garmin GPS puck data
#table0<-fread("data/track_points_20210903_LakeWaco.csv")
table0<-fread("FLAMe_GPS.dat")#,skip=c(2:3))
table<-table0 %>% as.data.frame()
table<-table[-c(2:3),]
names(table)<-table[1,]
table<-table[-1,]
#table$date<-rownames(table)
table$garmin_datetime_utc<-as.POSIXct(table$TIMESTAMP,tz="UTC")
table$garmin_datetime_local<-format(table$garmin_datetime_utc,tz="America/Chicago")
table_clean<-clean_names(table)
garmin<-table_clean
#garmin$garmin_datetime_utc <- ymd_hms(garmin$time_local_posi_xct)
#garmin$garmin_datetime_central<-format(garmin$garmin_datetime_utc, tz="UTC+5")
garmin$garmin_time_local<-substr(garmin$garmin_datetime_local,12,19)
garmin$xcoord<-as.numeric(garmin$longitude)
garmin$ycoord<-as.numeric(garmin$latitude)
garmin$date<-substr(garmin$garmin_datetime_local,1,10)
garmin$hour_dec<-hour(garmin$garmin_datetime_local)+
  minute(garmin$garmin_datetime_local)/60+
  second(garmin$garmin_datetime_local)/3600

garmin<-garmin %>% dplyr::select(garmin_datetime_utc,date,hour_dec,garmin_time_local,xcoord,ycoord,speed)
garmin<-garmin %>% filter(date=="2022-06-23", hour_dec>11.36, hour_dec<15.38)

getUTMzone(mean(garmin$xcoord,na.rm=TRUE))
# Zone 14, which is needed for UTM projection later



# 3. Merge data

data_merged<-merge(exo,garmin,by.x="lake_time_local",by.y="garmin_time_local")
data_merged<-data_merged %>% dplyr::select(lat_dec=ycoord,lon_dec=xcoord, speed,
                                           time_central=lake_time_local,date=exo_date_mm_dd_yyyy,datetime_utc=garmin_datetime_utc,
                                           wtemp_c=temp_c,sp_cond_m_s_cm,odo_mg_l,odo_percent_sat,
                                           chl=chlorophyll_rfu,tal_pc_rfu,pH,
                                           cond_m_s_cm,f_dom_rfu,#odo_percent_local,
                                           sal_psu,turb=turbidity_fnu, 
)

data_merged$date_time_c = paste(data_merged$date, data_merged$time_local)

lw623<- data_merged

# assigning unique id label 
lw623_id <- lw623 %>%
  mutate(label = row_number()) 
lw623_id$system <- "waco"

##########################################################################################################
# writing out lat/lon for df
lw623_pts <- lw623_id %>% dplyr::select(label,lat_dec,lon_dec) %>% na.omit()

# Split the full df into smaller subsets with <=5000 points each
lw623_subsets <- split(lw623_pts, rep(1:ceiling(nrow(lw623_pts)/5000), each=5000, length.out=nrow(lw623_pts)))
# Export each subset as a separate CSV file
for (i in seq_along(lw623_subsets)) {
  subset <- lw623_subsets[[i]]
  write.csv(subset, file = paste0("lw623pts_", i, ".csv"), row.names = FALSE)
}
##########################################################################################################


# import csv results from eepoints
lw6230 <- read_csv(c("s2bvs_lw623_1.csv","s2bvs_lw623_2.csv","s2bvs_lw623_3.csv")) %>% dplyr::arrange(label)

# filtering just to bands that will be used and calculating normalized diffs

lw623_nd <- lw6230 %>% dplyr::select(label,lat_dec, lon_dec, B2, B3, B4, B5, SCL) %>% 
  mutate(ndti = (B4-B3)/(B4+B3),
         ndci = (B5-B4)/(B5+B4))


lw623_s2flm <- merge(lw623_id, lw623_nd, by = "label") %>% dplyr::select(label, lat = lat_dec.x, lon = lon_dec.x,
                                                                system, speed, datetime_utc, date_time_c, 
                                                                chl, turb,B2, B3, B4, B5, SCL, ndti, ndci)

lw623s2flm <- lw623_s2flm %>% filter(SCL == 6) %>% filter(speed>0)

ggplot(lw62322,aes(turb, ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()

# assigning bands to R,G,B
R <- lw623s2flm$B4
G <- lw623s2flm$B3
B <-lw623s2flm$B2

s2lw623_df <- lw623s2flm %>% mutate(dwl = fui.hue(R, G, B))
write_csv(s2lw623_df, "lw623_s2flm.csv")

# combining with all other lakes just to see
lw62322<-read_csv("lw623_s2flm.csv")

s2flm_list623 <- list(ah_s2flm, bon_s2flm, lw62322, bw_s2flm, iv_s2flm, rb_s2flm)

#merge all data frames in list
s2flm623_df <- Reduce(function(x, y) merge(x, y, all=TRUE), s2flm_list623)







# lake waco 6/23 data is junk
###################################################################################################################
###############################################################################################################


# first viz transect points for each lake 
# ah (dam = 4800, S.A. = 2000) other option is (main body = 300, S.A.= 1780)
# bon (dam = 792, arm = 223)
# bw (dam = 5033, S.A. = 2393)
# iv (dam = 5078, arm1 = 7824)
# waco (dam = 6250, N.A. = 5295)
# rb (dam = 152, arm = 2350)
dm_sf <- st_as_sf(waco, coords = c("lon_dec", "lat_dec"), crs = 4326)
mapview(dm_sf)


# filtering to a one way transect from dam to arm
ah_dfd <- ah %>% filter(time_local>"10:20:00" & time_local<"11:09:00")%>%
  mutate(label = row_number()) %>% dplyr::select(label, lat_dec,lon_dec)
write_csv(ah_dfd, "ah_dfd.csv")

ah_dfd_id <- ah %>% filter(time_local>"10:20:00" & time_local<"11:09:00")%>%
  mutate(label = row_number())#%>% rename(lat=lat_dec, lon = lon_dec)
ah_dfd_id$system="arrowhead"

###
bon_dfd<-bon %>% filter(time_local>"12:03:44" & time_local<"12:13:15")%>%
  mutate(label = row_number()) %>% dplyr::select(label, lat_dec,lon_dec)
write_csv(bon_dfd, "bon_dfd.csv")

bon_dfd_id <- bon %>% filter(time_local>"12:03:44" & time_local<"12:13:15")%>%
  mutate(label = row_number())#%>% rename(lat=lat_dec, lon = lon_dec)
bon_dfd_id$system="bonham"
###
bw_dfd<-bw %>% filter(time_local>"13:34:35" & time_local<"14:19:06")%>%
  mutate(label = row_number()) %>% dplyr::select(label, lat_dec,lon_dec)
write_csv(bw_dfd, "bw_dfd.csv")

bw_dfd_id<-bw %>% filter(time_local>"13:34:35" & time_local<"14:19:06")%>%
  mutate(label = row_number())# %>% rename(lat=lat_dec, lon = lon_dec)
bw_dfd_id$system="brownwood"
###
iv_dfd<-iv %>% filter(time_local>"12:16:17" & time_local<"13:02:06")%>%
  mutate(label = row_number()) %>% dplyr::select(label, lat_dec,lon_dec)
write_csv(iv_dfd, "iv_dfd.csv")

iv_dfd_id<-iv %>% filter(time_local>"12:16:17" & time_local<"13:02:06")%>%
  mutate(label = row_number())
iv_dfd_id$system="ivie"
###
waco_dfd<-waco %>% filter(time_local>"12:14:44" & time_local<"12:30:56")%>%
  mutate(label = row_number()) %>% dplyr::select(label, lat_dec,lon_dec)
write_csv(waco_dfd, "waco_dfd.csv")

waco_dfd_id<-waco %>% filter(time_local>"12:14:44" & time_local<"12:30:56")%>%
  mutate(label = row_number())
waco_dfd_id$system="waco"
###
rb_dfd<-rb%>% filter(time_local>"10:44:29" & time_local<"11:21:36")%>%
  mutate(label = row_number()) %>% dplyr::select(label, lat_dec,lon_dec)
write_csv(rb_dfd, "rb_dfd.csv")

rb_dfd_id<-rb%>% filter(time_local>"10:44:29" & time_local<"11:21:36")%>%
  mutate(label = row_number())
rb_dfd_id$system="redbluff"

# reading in transect
# first arrowhead
ah_dfd_s2<-read_csv("ah_dfd_s2.csv")  %>% distinct(label, lat, lon, variable, distkm, .keep_all = T)
ah_pw<-ah_dfd_s2 %>% 
  pivot_wider(id_cols = c(label,lat, lon, distkm), names_from = variable, values_from = value) %>% 
  mutate(ndti = (B4-B3)/(B4+B3)) 

ah_merge<-merge(ah_dfd_id, ah_pw, by = "label") %>% dplyr::select(label, system, distkm, lat, lon,
                                                                  datetime_utc, date_time_c, speed,turb,
                                                                  B2, B3, B4, B5, SCL, ndti)

# bonham
bon_dfd_s2<-read_csv("bon_dfd_s2.csv")%>% distinct(label, lat, lon, variable, distkm, .keep_all = T)
bon_pw<-bon_dfd_s2%>% 
  pivot_wider(id_cols = c(label,lat, lon, distkm), names_from = variable, values_from = value) %>% 
  mutate(ndti = (B4-B3)/(B4+B3))   

bon_merge<-merge(bon_dfd_id, bon_pw, by = "label") %>% dplyr::select(label, system, distkm, lat, lon,
                                                                     datetime_utc, date_time_c, speed,turb,
                                                                     B2, B3, B4, B5, SCL, ndti)


# brownwood
bw_dfd_s2<-read_csv("bw_dfd_s2.csv") %>% distinct(label, lat, lon, variable, distkm, .keep_all = T)
bw_pw<-bw_dfd_s2 %>% 
  pivot_wider(id_cols = c(label,lat, lon, distkm), names_from = variable, values_from = value) %>% 
  mutate(ndti = (B4-B3)/(B4+B3))

bw_merge<-merge(bw_dfd_id, bw_pw, by = "label") %>% dplyr::select(label, system, distkm, lat, lon,
                                                                  datetime_utc, date_time_c, speed,turb,
                                                                  B2, B3, B4, B5, SCL, ndti)

# ivie
iv_dfd_s2<-read_csv("iv_dfd_s2.csv") %>% distinct(label, lat, lon, variable, distkm, .keep_all = T)
iv_pw<-iv_dfd_s2 %>% 
  pivot_wider(id_cols = c(label,lat, lon, distkm), names_from = variable, values_from = value) %>% 
  mutate(ndti = (B4-B3)/(B4+B3))

iv_merge<-merge(iv_dfd_id, iv_pw, by = "label") %>% dplyr::select(label, system, distkm, lat, lon,
                                                                  datetime_utc, date_time_c, speed,turb,
                                                                  B2, B3, B4, B5, SCL, ndti)
# waco
waco_dfd_s2<-read_csv("waco_dfd_s2.csv")%>% distinct(label, lat, lon, variable, distkm, .keep_all = T)
waco_pw<-waco_dfd_s2 %>% 
  pivot_wider(id_cols = c(label,lat, lon, distkm), names_from = variable, values_from = value) %>% 
  mutate(ndti = (B4-B3)/(B4+B3))

waco_merge<-merge(waco_dfd_id, waco_pw, by = "label") %>% dplyr::select(label, system, distkm, lat, lon,
                                                                        datetime_utc, date_time_c, speed,turb,
                                                                        B2, B3, B4, B5, SCL, ndti)
# redbluff
rb_dfd_s2<-read_csv("rb_dfd_s2.csv") %>% distinct(label, lat, lon, variable, distkm, .keep_all = T)
rb_pw<-rb_dfd_s2 %>% 
  pivot_wider(id_cols = c(label,lat, lon, distkm), names_from = variable, values_from = value) %>% 
  mutate(ndti = (B4-B3)/(B4+B3))

rb_merge<-merge(rb_dfd_id, rb_pw, by = "label") %>% dplyr::select(label, system, distkm, lat, lon,
                                                                  datetime_utc, date_time_c, speed,turb,
                                                                  B2, B3, B4, B5, SCL, ndti)



s2dfd_list <- list(ah_merge, bon_merge, waco_merge, bw_merge, iv_merge, rb_merge)

#merge all data frames in list
s2dfd_df <- Reduce(function(x, y) merge(x, y, all=TRUE), s2dfd_list)

s2transects<-s2dfd_df %>% filter(SCL == 6)

s2transects %>% ggplot(aes(distkm, turb, color = system)) + geom_point()

# this attempt was junk, does not show much. Trying now with points 