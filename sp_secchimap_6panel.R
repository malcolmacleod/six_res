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
library(ggmap)
library(sf)
library(viridis)
library(raster)
library(patchwork)


# 2. CRASR data - read in for all lakes --------------------------------------------
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

# combining NLA data (all reservoirs) with our own
nla_crasr<-read_csv("nla_crasr.csv")

nla_crasr <- nla_crasr%>% dplyr::rename(secchi = `secchi (m)`, turb = `turb (ntu)`)

# fit natural log log model with secchi and turbidity
tsdlm<-lm(log(secchi)~log(turb), data = nla_crasr)
summary(tsdlm) # exp(y) = 1.13734 -0.68001(exp(x))

nla_crasr$residuals <- tsdlm$residuals
nla_crasr$predicted<-tsdlm$fitted.values

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

# filtering garmin data by lake / date
bon1<-data_merged %>% filter(date=="2022-07-14", hour_dec>13.20, hour_dec<16.733) %>% filter(speed>0)
bon2<-data_merged %>% filter(date=="2022-07-15", hour_dec>12, hour_dec<13.91)%>% filter(speed>0)
waco<-data_merged %>% filter(date=="2022-07-23", hour_dec>10.767, hour_dec<13.25)%>% filter(speed>0)
bw<-data_merged %>% filter(date=="2022-08-06", hour_dec>12.91, hour_dec<15.283)%>% filter(speed>0)
iv<-data_merged %>% filter(date=="2022-08-07", hour_dec>10.67, hour_dec<13.68)%>% filter(speed>0)
rb<-data_merged %>% filter(date=="2022-08-08", hour_dec>10.7, hour_dec<15.967)%>% filter(speed>0)
ah<-data_merged %>% filter(date=="2022-08-12", hour_dec>9.8, hour_dec<12.8)%>% filter(speed>0)

# removing anamolous boat times for AH, IV, and LW
ah<- ah %>% filter(hour_dec<12.56 | hour_dec>12.59)
iv<-iv %>% filter(hour_dec<10.8 | hour_dec>10.83) %>% 
  filter(hour_dec<10.68 | hour_dec>10.69) %>% 
  filter(hour_dec<13.6 | hour_dec>13.65)
waco<-waco %>% filter(hour_dec<11.16 | hour_dec>11.18) %>% 
  filter(hour_dec<11.25 | hour_dec>11.27)

# 6c predicting secchi from turb values ------------------------------------------------
# ah$pred_sdd <- exp(predict(tsdlm, ah))
# bon1$pred_sdd <- exp(predict(tsdlm, bon1))
# bon2$pred_sdd <- exp(predict(tsdlm, bon2))
# waco$pred_sdd <- exp(predict(tsdlm, waco))
# bw$pred_sdd <- exp(predict(tsdlm, bw))
# iv$pred_sdd <- exp(predict(tsdlm, iv))
# rb$pred_sdd <- exp(predict(tsdlm, rb))
# 
# # calculating TSI from secchi
# ah$tsi_sd <- 10 * (6-log(ah$pred_sdd)/log(2))
# bon1$tsi_sd <- 10 * (6-log(bon1$pred_sdd)/log(2))
# bon2$tsi_sd <- 10 * (6-log(bon2$pred_sdd)/log(2))
# bw$tsi_sd <- 10 * (6-log(bw$pred_sdd)/log(2))
# iv$tsi_sd <- 10 * (6-log(iv$pred_sdd)/log(2))
# rb$tsi_sd <- 10 * (6-log(rb$pred_sdd)/log(2))
# waco$tsi_sd <- 10 * (6-log(waco$pred_sdd)/log(2))
# 
# # also adding secchi from Sherjah paper
# ah$sherjah<-(244.13 * ah$turb^-0.662) /100
# bon2$sherjah<-(244.13 * bon2$turb^-0.662) /100
# bw$sherjah<-(244.13 * bw$turb^-0.662) /100
# iv$sherjah<-(244.13 * iv$turb^-0.662) /100
# rb$sherjah<-(244.13 * rb$turb^-0.662) /100
# waco$sherjah<-(244.13 * waco$turb^-0.662) /100

flame_dwl<-read_csv("s2cp_boatpaths_zone.csv") 
flame_dwl$dwlgroup<-as.factor(flame_dwl$dwlgroup)

bon1<-flame_dwl %>% filter(system=="bonham")
bon2<-flame_dwl %>% filter(system=="bonham")
waco<-flame_dwl %>% filter(system=="waco")
bw<-flame_dwl %>% filter(system=="brownwood")
iv<-flame_dwl %>% filter(system=="ivie")
rb<-flame_dwl %>% filter(system=="redbluff")
ah<-flame_dwl %>%filter(system=="arrowhead")

###################################################################

my_palette <- colorRampPalette(rev(brewer.pal(9, "YlOrBr")))

# set some parameters to pass into ggplot 
ptsize<-0.4 #0.25
ptalpha<-0.5
ptalpha<-1

# attempts to have a common scale bar and similar area covered in bounding box

# Define the bounding box size (latitude and longitude degrees)
bbox_width_zoomed_in <- 0.05
bbox_height_zoomed_in <- 0.05
bbox_width <- 0.15
bbox_height <- 0.15
bbox_width_iv <- 0.25
bbox_height_iv <- 0.25
#bbox_width_bw <- 0.2

# create_plot <- function(center_lon, center_lat, bbox_width, bbox_height, df, title, add_scalebar = FALSE) {
#   require(OpenStreetMap)
#   require(ggplot2)
#   require(viridis)
#   
#   # Define bounding box
#   upperLeft <- c(center_lat + bbox_height / 2, center_lon - bbox_width / 2)
#   lowerRight <- c(center_lat - bbox_height / 2, center_lon + bbox_width / 2)
#   
#   message("Attempting to download map for: ", title)
#   message("UpperLeft: ", paste(upperLeft, collapse = ", "), " | LowerRight: ", paste(lowerRight, collapse = ", "))
#   
#   # Safe tile download
#   map <- tryCatch({
#     openmap(upperLeft, lowerRight, type = "apple-iphoto")
#   }, error = function(e) {
#     message("Tile download failed for ", title, ": ", e$message)
#     return(NULL)
#   })
#   
#   if (is.null(map)) {
#     return(ggplot() + ggtitle(paste("Map unavailable for", title)))
#   }
#   
#   map_projected <- openproj(map)
#   
#   plot <- autoplot.OpenStreetMap(map_projected) +
#     geom_jitter(
#       data = df,
#       size = ptsize,
#       alpha = ptalpha,
#       aes(x = lon_dec, y = lat_dec, color = turb),
#       inherit.aes = FALSE
#     ) +
#     scale_color_gradientn(
#       colors = viridis(100, option = "plasma"),
#       name = "Turbidity (NTU)",
#       limits = c(0, max(df$turb, na.rm = TRUE)),
#       breaks = pretty(range(df$turb, na.rm = TRUE), n = 5)
#     ) +
#     ggtitle(label = title) +
#     xlab("") + ylab("") +
#     theme(
#       axis.text = element_blank(),
#       axis.ticks = element_blank(),
#       legend.direction = "horizontal",
#       legend.text = element_text(size = 13),
#       legend.title = element_text(size = 13, hjust = 0.5),
#       legend.key.width = unit(1.5, "cm"),
#       legend.position = c(0.5, -0.17),
#       plot.title = element_text(hjust = 0.5),
#       plot.margin = ggplot2::margin(0, 0, 0, 0)
#     ) +
#     guides(col = guide_colorbar(title.position = "top"))
#   
#   return(plot)
# }

create_plot <- function(center_lon, center_lat, bbox_width, bbox_height, df, title, add_scalebar=FALSE) {
  # Define the bounding box
  upperLeft <- c(center_lat + bbox_height / 2, center_lon - bbox_width / 2)
  lowerRight <- c(center_lat - bbox_height / 2, center_lon + bbox_width / 2)

  # Fetch the map
  map <- openmap(upperLeft, lowerRight, type = "apple-iphoto")
  map_projected <- OpenStreetMap::openproj(map)
  #  map_projected <- OpenStreetMap::openproj(map, projection = "+init=epsg:4326")
  #  map_projected <- OpenStreetMap::openproj(map, projection = "EPSG:4326")

  # Create the plot
  plot <- OpenStreetMap::autoplot.OpenStreetMap(map_projected) +
    geom_jitter(data = df,
                size = ptsize,
                alpha = ptalpha,
                aes(x = lon_dec, y = lat_dec, color = turb)) +
    # set the same color limits for every map
    #    scale_color_gradientn(colors = my_palette(100), name = "Secchi Depth (m)",limits=c(0,1.5)) +
    scale_color_gradientn(colors = rev(viridis(100,option="plasma")), name = "Turbidity (NTU)") + #,limits=c(0.1,1.5),breaks=c(0.2,0.4,0.6,0.8,1,1.2,1.4)) +
    #    scale_color_gradientn(colors = rev(viridis(256)), name = "Secchi Depth (m)",limits=c(0.0,1.5)) +

    ggtitle(label = title) +
    xlab("") + ylab("") +
    theme(
      #      plot.margin = unit(c(-1, 0, -1, 0), "lines"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.direction = "horizontal",
      legend.title.align = 0.5,
      legend.text = element_text(size = 13),     # Set the legend text size
      legend.title = element_text(size = 13),     # Set the legend title size
      #      legend.key.width = unit(0.75, "cm"),
      legend.key.width = unit(1.8, "cm"),
      #      legend.position = c(0.5, 0.02),
      legend.position = c(0.5, -0.17),
      #        legend.position = "bottom",              # Position the legend at the bottom
      #        legend.justification = "center"
    ) +
    guides(col = guide_colorbar(title.position = "top")) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = margin(0, 0, 0, 0)) # reduce margins of individual plots

  return(plot)
}

# Define the center coordinates for each map
centers <- list(
  c(-96.15, 33.65), # Center for bon
  c(-97.25, 31.56), # Center for waco
  c(-99.07, 31.86),  # Center for bw
#  c(-99.06, 31.86),  # Center for bw
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

# prepare dimensions for scale bars - doing it the manual way, a bit verbose
barheight<-0.92
baredge<-0.98
textheight<-0.95
map1_scalexcoords<-c(map1$layers[[1]]$data$x[1] + ((map1$layers[[1]]$data$x[2] - map1$layers[[1]]$data$x[1])*(baredge-0.898)), # left edge of scale bar
                     map1$layers[[1]]$data$x[1] + ((map1$layers[[1]]$data$x[2] - map1$layers[[1]]$data$x[1])*baredge))  # right edge of scale bat
map1_scaleycoords<-c(map1$layers[[1]]$data$y[1] + ((map1$layers[[1]]$data$y[2] - map1$layers[[1]]$data$y[1])*barheight), # top of scale bar (left)
                     map1$layers[[1]]$data$y[1] + ((map1$layers[[1]]$data$y[2] - map1$layers[[1]]$data$y[1])*barheight), # top of scale bar (right)
                     map1$layers[[1]]$data$y[1] + ((map1$layers[[1]]$data$y[2] - map1$layers[[1]]$data$y[1])*textheight)) # label height

map2_scalexcoords<-c(map2$layers[[1]]$data$x[1] + ((map2$layers[[1]]$data$x[2] - map2$layers[[1]]$data$x[1])*(baredge-0.3746)),
                     map2$layers[[1]]$data$x[1] + ((map2$layers[[1]]$data$x[2] - map2$layers[[1]]$data$x[1])*baredge))
map2_scaleycoords<-c(map2$layers[[1]]$data$y[1] + ((map2$layers[[1]]$data$y[2] - map2$layers[[1]]$data$y[1])*barheight),
                     map2$layers[[1]]$data$y[1] + ((map2$layers[[1]]$data$y[2] - map2$layers[[1]]$data$y[1])*barheight),
                     map2$layers[[1]]$data$y[1] + ((map2$layers[[1]]$data$y[2] - map2$layers[[1]]$data$y[1])*textheight))
map3_scalexcoords<-c(map3$layers[[1]]$data$x[1] + ((map3$layers[[1]]$data$x[2] - map3$layers[[1]]$data$x[1])*(baredge-0.3757)),
                     map3$layers[[1]]$data$x[1] + ((map3$layers[[1]]$data$x[2] - map3$layers[[1]]$data$x[1])*baredge))
map3_scaleycoords<-c(map3$layers[[1]]$data$y[1] + ((map3$layers[[1]]$data$y[2] - map3$layers[[1]]$data$y[1])*barheight),
                     map3$layers[[1]]$data$y[1] + ((map3$layers[[1]]$data$y[2] - map3$layers[[1]]$data$y[1])*barheight),
                     map3$layers[[1]]$data$y[1] + ((map3$layers[[1]]$data$y[2] - map3$layers[[1]]$data$y[1])*textheight))
map4_scalexcoords<-c(map4$layers[[1]]$data$x[1] + ((map4$layers[[1]]$data$x[2] - map4$layers[[1]]$data$x[1])*(baredge-0.2663)),
                     map4$layers[[1]]$data$x[1] + ((map4$layers[[1]]$data$x[2] - map4$layers[[1]]$data$x[1])*baredge))
map4_scaleycoords<-c(map4$layers[[1]]$data$y[1] + ((map4$layers[[1]]$data$y[2] - map4$layers[[1]]$data$y[1])*barheight),
                     map4$layers[[1]]$data$y[1] + ((map4$layers[[1]]$data$y[2] - map4$layers[[1]]$data$y[1])*barheight),
                     map4$layers[[1]]$data$y[1] + ((map4$layers[[1]]$data$y[2] - map4$layers[[1]]$data$y[1])*textheight))
map5_scalexcoords<-c(map5$layers[[1]]$data$x[1] + ((map5$layers[[1]]$data$x[2] - map5$layers[[1]]$data$x[1])*(baredge-0.3757)),
                     map5$layers[[1]]$data$x[1] + ((map5$layers[[1]]$data$x[2] - map5$layers[[1]]$data$x[1])*baredge))
map5_scaleycoords<-c(map5$layers[[1]]$data$y[1] + ((map5$layers[[1]]$data$y[2] - map5$layers[[1]]$data$y[1])*barheight),
                     map5$layers[[1]]$data$y[1] + ((map5$layers[[1]]$data$y[2] - map5$layers[[1]]$data$y[1])*barheight),
                     map5$layers[[1]]$data$y[1] + ((map5$layers[[1]]$data$y[2] - map5$layers[[1]]$data$y[1])*textheight))
map6_scalexcoords<-c(map6$layers[[1]]$data$x[1] + ((map6$layers[[1]]$data$x[2] - map6$layers[[1]]$data$x[1])*(baredge-0.3812)),
                     map6$layers[[1]]$data$x[1] + ((map6$layers[[1]]$data$x[2] - map6$layers[[1]]$data$x[1])*baredge))
map6_scaleycoords<-c(map6$layers[[1]]$data$y[1] + ((map6$layers[[1]]$data$y[2] - map6$layers[[1]]$data$y[1])*barheight),
                     map6$layers[[1]]$data$y[1] + ((map6$layers[[1]]$data$y[2] - map6$layers[[1]]$data$y[1])*barheight),
                     map6$layers[[1]]$data$y[1] + ((map6$layers[[1]]$data$y[2] - map6$layers[[1]]$data$y[1])*textheight))

# check lengths of each scale bar before adding to map
pointDistance(p1=c(map1_scalexcoords[1],map1_scaleycoords[1]),
              p2=c(map1_scalexcoords[2],map1_scaleycoords[2]),
              lonlat=TRUE)
pointDistance(p1=c(map2_scalexcoords[1],map2_scaleycoords[1]),
              p2=c(map2_scalexcoords[2],map2_scaleycoords[2]),
              lonlat=TRUE)
pointDistance(p1=c(map3_scalexcoords[1],map3_scaleycoords[1]),
              p2=c(map3_scalexcoords[2],map3_scaleycoords[2]),
              lonlat=TRUE)
pointDistance(p1=c(map4_scalexcoords[1],map4_scaleycoords[1]),
              p2=c(map4_scalexcoords[2],map4_scaleycoords[2]),
              lonlat=TRUE)
pointDistance(p1=c(map5_scalexcoords[1],map5_scaleycoords[1]),
              p2=c(map5_scalexcoords[2],map5_scaleycoords[2]),
              lonlat=TRUE)
pointDistance(p1=c(map6_scalexcoords[1],map6_scaleycoords[1]),
              p2=c(map6_scalexcoords[2],map6_scaleycoords[2]),
              lonlat=TRUE)

# add scale bar and dam location to each map, use one legend
map1<- map1 + theme(legend.position = "none") +
  geom_segment(x = -96.135058, y = 33.646397, xend = -96.135981, yend = 33.655287,lwd=2,color=gray(0.6))+ # dam location
  geom_segment(x = map1_scalexcoords[1],y=map1_scaleycoords[2],xend=map1_scalexcoords[2],yend=map1_scaleycoords[2],lwd=1.5,color=gray(1)) #+ # scale bar part
#  annotate("text",x = mean(map1_scalexcoords), y = map1_scaleycoords[3],color=gray(1),label = "4 km") # label for scale bar
map2<- map2 + theme(legend.position = "none") +
  geom_segment(x = -97.1910, y = 31.576152, xend = -97.2136, yend = 31.594505,lwd=2,color=gray(0.6))+
  geom_segment(x = map2_scalexcoords[1],y=map2_scaleycoords[2],xend=map2_scalexcoords[2],yend=map2_scaleycoords[2],lwd=1.5,color=gray(1)) #+
#  annotate("text",x = mean(map2_scalexcoords), y = map2_scaleycoords[3],color=gray(1),label = "4 km")
map3<- map3 + theme(legend.position = "none") +
  geom_segment(x = -99.0031, y = 31.8352, xend = -99.0002, yend = 31.8421,lwd=2,color=gray(0.6))+
  geom_segment(x = map3_scalexcoords[1],y=map3_scaleycoords[2],xend=map3_scalexcoords[2],yend=map3_scaleycoords[2],lwd=1.5,color=gray(1)) #+
#  annotate("text",x = mean(map3_scalexcoords), y = map3_scaleycoords[3],color=gray(1),label = "4 km")
map4<- map4 + theme(legend.position = "none") +
  geom_segment(x = -99.6684, y = 31.4970, xend = -99.6640, yend = 31.5025,lwd=2,color=gray(0.6))+
  geom_segment(x = map4_scalexcoords[1],y=map4_scaleycoords[2],xend=map4_scalexcoords[2],yend=map4_scaleycoords[2],lwd=1.5,color=gray(1)) #+
#  annotate("text",x = mean(map4_scalexcoords), y = map4_scaleycoords[3],color=gray(1),label = "4 km")
map5<- map5 + theme(legend.position = "none") +
  geom_segment(x = -103.9121, y = 31.8962, xend = -103.9080, yend = 31.9069,lwd=2,color=gray(0.6))+
  geom_segment(x = map5_scalexcoords[1],y=map5_scaleycoords[2],xend=map5_scalexcoords[2],yend=map5_scaleycoords[2],lwd=1.5,color=gray(1)) #+
#  annotate("text",x = mean(map5_scalexcoords), y = map5_scaleycoords[3],color=gray(1),label = "4 km")
map6<- map6 + theme(legend.position = "none") +
  geom_segment(x = -98.3539, y = 33.7653, xend = -98.3753, yend = 33.7657,lwd=2,color=gray(0.6))+
  geom_segment(x = map6_scalexcoords[1],y=map6_scaleycoords[2],xend=map6_scalexcoords[2],yend=map6_scaleycoords[2],lwd=1.5,color=gray(1)) #+
#  annotate("text",x = mean(map6_scalexcoords), y = map6_scaleycoords[3],color=gray(1),label = "4 km")


# Combine plots into a single layout with grid.arrange
# Does this panel order make sense? 
combo_sdd_lomarg<-grid.arrange(map5, map4, map3, map6, map2, map1, ncol = 3)#, padding = unit(c(-1, -1, -1, -1), "cm"))
ggsave("combo_turb.png", combo_sdd_lomarg,width=9,height=8)



ggplot(ah, aes(x = lon_dec, y = lat_dec, color = turb)) +
  geom_jitter(size = 3) +
  scale_color_gradientn(
    colors = viridis::viridis(100, option = "plasma"),
    name = "Turbidity (NTU)",
    limits = c(0, max(ah$turb, na.rm = TRUE)),
    breaks = pretty(range(ah$turb, na.rm = TRUE), n = 5)
  ) +
  theme_minimal()

##################################################################################################
# set some parameters to pass into ggplot 
ptsize<-0.4 #0.25
ptalpha<-0.5
ptalpha<-1

# attempts to have a common scale bar and similar area covered in bounding box

# Define the bounding box size (latitude and longitude degrees)
bbox_width_zoomed_in <- 0.05
bbox_height_zoomed_in <- 0.05
bbox_width <- 0.15
bbox_height <- 0.15
bbox_width_iv <- 0.25
bbox_height_iv <- 0.25
#bbox_width_bw <- 0.2

create_plot <- function(center_lon, center_lat, bbox_width, bbox_height, df, title, add_scalebar=FALSE) {
  # Define the bounding box
  upperLeft <- c(center_lat + bbox_height / 2, center_lon - bbox_width / 2)
  lowerRight <- c(center_lat - bbox_height / 2, center_lon + bbox_width / 2)
  
  # Fetch the map
  map <- openmap(upperLeft, lowerRight, type = "esri-imagery")
  map_projected <- OpenStreetMap::openproj(map)
  #  map_projected <- OpenStreetMap::openproj(map, projection = "+init=epsg:4326")
  #  map_projected <- OpenStreetMap::openproj(map, projection = "EPSG:4326")
  
  # Create the plot
  plot <- OpenStreetMap::autoplot.OpenStreetMap(map_projected) +
    geom_jitter(data = df, 
                size = ptsize,
                alpha = ptalpha,
                aes(x = lon_dec, y = lat_dec, color = turb)) +
    # set the same color limits for every map
    #    scale_color_gradientn(colors = my_palette(100), name = "Secchi Depth (m)",limits=c(0,1.5)) +
    scale_color_gradientn(colors = rev(viridis(100,option="D")), name = "Turbidity (NTU)") +
    #    scale_color_gradientn(colors = rev(viridis(256)), name = "Secchi Depth (m)",limits=c(0.0,1.5)) +
    
    ggtitle(label = title) +
    xlab("") + ylab("") +
    theme(
      plot.margin = unit(c(-1, 0, -1, 0), "lines"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.direction = "horizontal",
      legend.title.align = 0.5,
      legend.text = element_text(size = 13),     # Set the legend text size
      legend.title = element_text(size = 13),     # Set the legend title size
      # legend.key.width = unit(0.75, "cm"),
      legend.key.width = unit(1.5, "cm"),
      legend.position = c(0.5, 1.3),
      #legend.position = c(0.5, -0.17),
      #        legend.position = "bottom",              # Position the legend at the bottom
      #  legend.justification = "center"
    ) +
    guides(col = guide_colorbar(title.position = "top")) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = margin(0, 0, 0, 0)) # reduce margins of individual plots
  
  return(plot)
}

# Define the center coordinates for each map
centers <- list(
  c(-96.15, 33.65), # Center for bon
  c(-97.25, 31.56), # Center for waco
  c(-99.07, 31.86),  # Center for bw
  #  c(-99.06, 31.86),  # Center for bw
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

# prepare dimensions for scale bars - doing it the manual way, a bit verbose
barheight<-0.92
baredge<-0.98
textheight<-0.95
map1_scalexcoords<-c(map1$layers[[1]]$data$x[1] + ((map1$layers[[1]]$data$x[2] - map1$layers[[1]]$data$x[1])*(baredge-0.898)), # left edge of scale bar
                     map1$layers[[1]]$data$x[1] + ((map1$layers[[1]]$data$x[2] - map1$layers[[1]]$data$x[1])*baredge))  # right edge of scale bat
map1_scaleycoords<-c(map1$layers[[1]]$data$y[1] + ((map1$layers[[1]]$data$y[2] - map1$layers[[1]]$data$y[1])*barheight), # top of scale bar (left)
                     map1$layers[[1]]$data$y[1] + ((map1$layers[[1]]$data$y[2] - map1$layers[[1]]$data$y[1])*barheight), # top of scale bar (right)
                     map1$layers[[1]]$data$y[1] + ((map1$layers[[1]]$data$y[2] - map1$layers[[1]]$data$y[1])*textheight)) # label height

map2_scalexcoords<-c(map2$layers[[1]]$data$x[1] + ((map2$layers[[1]]$data$x[2] - map2$layers[[1]]$data$x[1])*(baredge-0.3746)),
                     map2$layers[[1]]$data$x[1] + ((map2$layers[[1]]$data$x[2] - map2$layers[[1]]$data$x[1])*baredge))
map2_scaleycoords<-c(map2$layers[[1]]$data$y[1] + ((map2$layers[[1]]$data$y[2] - map2$layers[[1]]$data$y[1])*barheight),
                     map2$layers[[1]]$data$y[1] + ((map2$layers[[1]]$data$y[2] - map2$layers[[1]]$data$y[1])*barheight),
                     map2$layers[[1]]$data$y[1] + ((map2$layers[[1]]$data$y[2] - map2$layers[[1]]$data$y[1])*textheight))
map3_scalexcoords<-c(map3$layers[[1]]$data$x[1] + ((map3$layers[[1]]$data$x[2] - map3$layers[[1]]$data$x[1])*(baredge-0.3757)),
                     map3$layers[[1]]$data$x[1] + ((map3$layers[[1]]$data$x[2] - map3$layers[[1]]$data$x[1])*baredge))
map3_scaleycoords<-c(map3$layers[[1]]$data$y[1] + ((map3$layers[[1]]$data$y[2] - map3$layers[[1]]$data$y[1])*barheight),
                     map3$layers[[1]]$data$y[1] + ((map3$layers[[1]]$data$y[2] - map3$layers[[1]]$data$y[1])*barheight),
                     map3$layers[[1]]$data$y[1] + ((map3$layers[[1]]$data$y[2] - map3$layers[[1]]$data$y[1])*textheight))
map4_scalexcoords<-c(map4$layers[[1]]$data$x[1] + ((map4$layers[[1]]$data$x[2] - map4$layers[[1]]$data$x[1])*(baredge-0.2663)),
                     map4$layers[[1]]$data$x[1] + ((map4$layers[[1]]$data$x[2] - map4$layers[[1]]$data$x[1])*baredge))
map4_scaleycoords<-c(map4$layers[[1]]$data$y[1] + ((map4$layers[[1]]$data$y[2] - map4$layers[[1]]$data$y[1])*barheight),
                     map4$layers[[1]]$data$y[1] + ((map4$layers[[1]]$data$y[2] - map4$layers[[1]]$data$y[1])*barheight),
                     map4$layers[[1]]$data$y[1] + ((map4$layers[[1]]$data$y[2] - map4$layers[[1]]$data$y[1])*textheight))
map5_scalexcoords<-c(map5$layers[[1]]$data$x[1] + ((map5$layers[[1]]$data$x[2] - map5$layers[[1]]$data$x[1])*(baredge-0.3757)),
                     map5$layers[[1]]$data$x[1] + ((map5$layers[[1]]$data$x[2] - map5$layers[[1]]$data$x[1])*baredge))
map5_scaleycoords<-c(map5$layers[[1]]$data$y[1] + ((map5$layers[[1]]$data$y[2] - map5$layers[[1]]$data$y[1])*barheight),
                     map5$layers[[1]]$data$y[1] + ((map5$layers[[1]]$data$y[2] - map5$layers[[1]]$data$y[1])*barheight),
                     map5$layers[[1]]$data$y[1] + ((map5$layers[[1]]$data$y[2] - map5$layers[[1]]$data$y[1])*textheight))
map6_scalexcoords<-c(map6$layers[[1]]$data$x[1] + ((map6$layers[[1]]$data$x[2] - map6$layers[[1]]$data$x[1])*(baredge-0.3812)),
                     map6$layers[[1]]$data$x[1] + ((map6$layers[[1]]$data$x[2] - map6$layers[[1]]$data$x[1])*baredge))
map6_scaleycoords<-c(map6$layers[[1]]$data$y[1] + ((map6$layers[[1]]$data$y[2] - map6$layers[[1]]$data$y[1])*barheight),
                     map6$layers[[1]]$data$y[1] + ((map6$layers[[1]]$data$y[2] - map6$layers[[1]]$data$y[1])*barheight),
                     map6$layers[[1]]$data$y[1] + ((map6$layers[[1]]$data$y[2] - map6$layers[[1]]$data$y[1])*textheight))

# check lengths of each scale bar before adding to map
pointDistance(p1=c(map1_scalexcoords[1],map1_scaleycoords[1]),
              p2=c(map1_scalexcoords[2],map1_scaleycoords[2]),
              lonlat=TRUE)
pointDistance(p1=c(map2_scalexcoords[1],map2_scaleycoords[1]),
              p2=c(map2_scalexcoords[2],map2_scaleycoords[2]),
              lonlat=TRUE)
pointDistance(p1=c(map3_scalexcoords[1],map3_scaleycoords[1]),
              p2=c(map3_scalexcoords[2],map3_scaleycoords[2]),
              lonlat=TRUE)
pointDistance(p1=c(map4_scalexcoords[1],map4_scaleycoords[1]),
              p2=c(map4_scalexcoords[2],map4_scaleycoords[2]),
              lonlat=TRUE)
pointDistance(p1=c(map5_scalexcoords[1],map5_scaleycoords[1]),
              p2=c(map5_scalexcoords[2],map5_scaleycoords[2]),
              lonlat=TRUE)
pointDistance(p1=c(map6_scalexcoords[1],map6_scaleycoords[1]),
              p2=c(map6_scalexcoords[2],map6_scaleycoords[2]),
              lonlat=TRUE)

# add scale bar and dam location to each map, use one legend
map1<- map1 + theme(legend.position = "none") +
  geom_segment(x = -96.135058, y = 33.646397, xend = -96.135981, yend = 33.655287,lwd=2,color=gray(0.6))+ # dam location
  geom_segment(x = map1_scalexcoords[1],y=map1_scaleycoords[2],xend=map1_scalexcoords[2],yend=map1_scaleycoords[2],lwd=1.5,color=gray(1)) #+ # scale bar part
#  annotate("text",x = mean(map1_scalexcoords), y = map1_scaleycoords[3],color=gray(1),label = "4 km") # label for scale bar
map2<- map2 + #theme(legend.position = "none") +
  geom_segment(x = -97.1910, y = 31.576152, xend = -97.2136, yend = 31.594505,lwd=2,color=gray(0.6))+
  geom_segment(x = map2_scalexcoords[1],y=map2_scaleycoords[2],xend=map2_scalexcoords[2],yend=map2_scaleycoords[2],lwd=1.5,color=gray(1)) #+
#  annotate("text",x = mean(map2_scalexcoords), y = map2_scaleycoords[3],color=gray(1),label = "4 km")
map3<- map3 + theme(legend.position = "none") +
  geom_segment(x = -99.0031, y = 31.8352, xend = -99.0002, yend = 31.8421,lwd=2,color=gray(0.6))+
  geom_segment(x = map3_scalexcoords[1],y=map3_scaleycoords[2],xend=map3_scalexcoords[2],yend=map3_scaleycoords[2],lwd=1.5,color=gray(1)) #+
#  annotate("text",x = mean(map3_scalexcoords), y = map3_scaleycoords[3],color=gray(1),label = "4 km")
map4<- map4 + theme(legend.position = "none") +
  geom_segment(x = -99.6684, y = 31.4970, xend = -99.6640, yend = 31.5025,lwd=2,color=gray(0.6))+
  geom_segment(x = map4_scalexcoords[1],y=map4_scaleycoords[2],xend=map4_scalexcoords[2],yend=map4_scaleycoords[2],lwd=1.5,color=gray(1)) #+
#  annotate("text",x = mean(map4_scalexcoords), y = map4_scaleycoords[3],color=gray(1),label = "4 km")
map5<- map5 + theme(legend.position = "none") +
  geom_segment(x = -103.9121, y = 31.8962, xend = -103.9080, yend = 31.9069,lwd=2,color=gray(0.6))+
  geom_segment(x = map5_scalexcoords[1],y=map5_scaleycoords[2],xend=map5_scalexcoords[2],yend=map5_scaleycoords[2],lwd=1.5,color=gray(1)) #+
#  annotate("text",x = mean(map5_scalexcoords), y = map5_scaleycoords[3],color=gray(1),label = "4 km")
map6<- map6 + theme(legend.position = "none") +
  geom_segment(x = -98.3539, y = 33.7653, xend = -98.3753, yend = 33.7657,lwd=2,color=gray(0.6))+
  geom_segment(x = map6_scalexcoords[1],y=map6_scaleycoords[2],xend=map6_scalexcoords[2],yend=map6_scaleycoords[2],lwd=1.5,color=gray(1)) #+
#  annotate("text",x = mean(map6_scalexcoords), y = map6_scaleycoords[3],color=gray(1),label = "4 km")


# Combine plots into a single layout with grid.arrange
# Does this panel order make sense? 
combo_turb_lomarg<-grid.arrange(map5, map4, map3, map6, map2, map1, ncol = 3)#, padding = unit(c(-1, -1, -1, -1), "cm"))
ggsave("combo_turb_update.png", combo_turb_lomarg,width=9,height=8)

###############################################################################

fui_palette<-c("1" = "#2158bc","2" = "#316dc5","3" = "#327cbb","4" = "#4b80a0",
               "5" = "#568f96","6" = "#6d9298", "7" = "#698c86","8" = "#759e72",
               "9" = "#7ba654","10" = "#7dae38","11" = "#94b660","12" = "#94b660", 
               "13" = "#a5bc76", "14" = "#aab86d","15" = "#adb55f","16" = "#a8a965",
               "17" = "#ae9f5c","18" = "#b3a053"
               ,"19" = "#af8a44","20" = "#a46905","21" = "#9f4d04")





# set some parameters to pass into ggplot 
ptsize<-0.4 #0.25
ptalpha<-0.5
ptalpha<-1

# attempts to have a common scale bar and similar area covered in bounding box

# Define the bounding box size (latitude and longitude degrees)
bbox_width_zoomed_in <- 0.05
bbox_height_zoomed_in <- 0.05
bbox_width <- 0.15
bbox_height <- 0.15
bbox_width_iv <- 0.25
bbox_height_iv <- 0.25
#bbox_width_bw <- 0.2

create_plot <- function(center_lon, center_lat, bbox_width, bbox_height, df, title, add_scalebar=FALSE) {
  # Define the bounding box
  upperLeft <- c(center_lat + bbox_height / 2, center_lon - bbox_width / 2)
  lowerRight <- c(center_lat - bbox_height / 2, center_lon + bbox_width / 2)
  
  # Fetch the map
  map <- openmap(upperLeft, lowerRight, type = "osm")
  map_projected <- OpenStreetMap::openproj(map)
  #  map_projected <- OpenStreetMap::openproj(map, projection = "+init=epsg:4326")
  #  map_projected <- OpenStreetMap::openproj(map, projection = "EPSG:4326")
  
  # Create the plot
  plot <- OpenStreetMap::autoplot.OpenStreetMap(map_projected) +
    geom_jitter(data = df, 
                size = ptsize,
                alpha = ptalpha,
                aes(x = lon_dec, y = lat_dec, color = dwlgroup)) +
    # set the same color limits for every map
    #    scale_color_gradientn(colors = my_palette(100), name = "Secchi Depth (m)",limits=c(0,1.5)) +
    scale_color_manual(values = fui_palette, name = "Dominant Wavelength (nm)") +
    #    scale_color_gradientn(colors = rev(viridis(256)), name = "Secchi Depth (m)",limits=c(0.0,1.5)) +
    
    ggtitle(label = title) +
    xlab("") + ylab("") +
    theme(
      plot.margin = unit(c(-1, 0, -1, 0), "lines"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.direction = "horizontal",
      legend.title.align = 0.5,
      legend.text = element_text(size = 13),     # Set the legend text size
      legend.title = element_text(size = 13),     # Set the legend title size
      # legend.key.width = unit(0.75, "cm"),
      legend.key.width = unit(1.5, "cm"),
      legend.position = c(0.5, 1.3),
      #legend.position = c(0.5, -0.17),
      #        legend.position = "bottom",              # Position the legend at the bottom
      #  legend.justification = "center"
    ) +
    guides(color = guide_legend(title.position = "top")) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = margin(0, 0, 0, 0)) # reduce margins of individual plots
  
  return(plot)
}

# Define the center coordinates for each map
centers <- list(
  c(-96.15, 33.65), # Center for bon
  c(-97.25, 31.56), # Center for waco
  c(-99.07, 31.86),  # Center for bw
  #  c(-99.06, 31.86),  # Center for bw
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

# prepare dimensions for scale bars - doing it the manual way, a bit verbose
barheight<-0.92
baredge<-0.98
textheight<-0.95
map1_scalexcoords<-c(map1$layers[[1]]$data$x[1] + ((map1$layers[[1]]$data$x[2] - map1$layers[[1]]$data$x[1])*(baredge-0.898)), # left edge of scale bar
                     map1$layers[[1]]$data$x[1] + ((map1$layers[[1]]$data$x[2] - map1$layers[[1]]$data$x[1])*baredge))  # right edge of scale bat
map1_scaleycoords<-c(map1$layers[[1]]$data$y[1] + ((map1$layers[[1]]$data$y[2] - map1$layers[[1]]$data$y[1])*barheight), # top of scale bar (left)
                     map1$layers[[1]]$data$y[1] + ((map1$layers[[1]]$data$y[2] - map1$layers[[1]]$data$y[1])*barheight), # top of scale bar (right)
                     map1$layers[[1]]$data$y[1] + ((map1$layers[[1]]$data$y[2] - map1$layers[[1]]$data$y[1])*textheight)) # label height

map2_scalexcoords<-c(map2$layers[[1]]$data$x[1] + ((map2$layers[[1]]$data$x[2] - map2$layers[[1]]$data$x[1])*(baredge-0.3746)),
                     map2$layers[[1]]$data$x[1] + ((map2$layers[[1]]$data$x[2] - map2$layers[[1]]$data$x[1])*baredge))
map2_scaleycoords<-c(map2$layers[[1]]$data$y[1] + ((map2$layers[[1]]$data$y[2] - map2$layers[[1]]$data$y[1])*barheight),
                     map2$layers[[1]]$data$y[1] + ((map2$layers[[1]]$data$y[2] - map2$layers[[1]]$data$y[1])*barheight),
                     map2$layers[[1]]$data$y[1] + ((map2$layers[[1]]$data$y[2] - map2$layers[[1]]$data$y[1])*textheight))
map3_scalexcoords<-c(map3$layers[[1]]$data$x[1] + ((map3$layers[[1]]$data$x[2] - map3$layers[[1]]$data$x[1])*(baredge-0.3757)),
                     map3$layers[[1]]$data$x[1] + ((map3$layers[[1]]$data$x[2] - map3$layers[[1]]$data$x[1])*baredge))
map3_scaleycoords<-c(map3$layers[[1]]$data$y[1] + ((map3$layers[[1]]$data$y[2] - map3$layers[[1]]$data$y[1])*barheight),
                     map3$layers[[1]]$data$y[1] + ((map3$layers[[1]]$data$y[2] - map3$layers[[1]]$data$y[1])*barheight),
                     map3$layers[[1]]$data$y[1] + ((map3$layers[[1]]$data$y[2] - map3$layers[[1]]$data$y[1])*textheight))
map4_scalexcoords<-c(map4$layers[[1]]$data$x[1] + ((map4$layers[[1]]$data$x[2] - map4$layers[[1]]$data$x[1])*(baredge-0.2663)),
                     map4$layers[[1]]$data$x[1] + ((map4$layers[[1]]$data$x[2] - map4$layers[[1]]$data$x[1])*baredge))
map4_scaleycoords<-c(map4$layers[[1]]$data$y[1] + ((map4$layers[[1]]$data$y[2] - map4$layers[[1]]$data$y[1])*barheight),
                     map4$layers[[1]]$data$y[1] + ((map4$layers[[1]]$data$y[2] - map4$layers[[1]]$data$y[1])*barheight),
                     map4$layers[[1]]$data$y[1] + ((map4$layers[[1]]$data$y[2] - map4$layers[[1]]$data$y[1])*textheight))
map5_scalexcoords<-c(map5$layers[[1]]$data$x[1] + ((map5$layers[[1]]$data$x[2] - map5$layers[[1]]$data$x[1])*(baredge-0.3757)),
                     map5$layers[[1]]$data$x[1] + ((map5$layers[[1]]$data$x[2] - map5$layers[[1]]$data$x[1])*baredge))
map5_scaleycoords<-c(map5$layers[[1]]$data$y[1] + ((map5$layers[[1]]$data$y[2] - map5$layers[[1]]$data$y[1])*barheight),
                     map5$layers[[1]]$data$y[1] + ((map5$layers[[1]]$data$y[2] - map5$layers[[1]]$data$y[1])*barheight),
                     map5$layers[[1]]$data$y[1] + ((map5$layers[[1]]$data$y[2] - map5$layers[[1]]$data$y[1])*textheight))
map6_scalexcoords<-c(map6$layers[[1]]$data$x[1] + ((map6$layers[[1]]$data$x[2] - map6$layers[[1]]$data$x[1])*(baredge-0.3812)),
                     map6$layers[[1]]$data$x[1] + ((map6$layers[[1]]$data$x[2] - map6$layers[[1]]$data$x[1])*baredge))
map6_scaleycoords<-c(map6$layers[[1]]$data$y[1] + ((map6$layers[[1]]$data$y[2] - map6$layers[[1]]$data$y[1])*barheight),
                     map6$layers[[1]]$data$y[1] + ((map6$layers[[1]]$data$y[2] - map6$layers[[1]]$data$y[1])*barheight),
                     map6$layers[[1]]$data$y[1] + ((map6$layers[[1]]$data$y[2] - map6$layers[[1]]$data$y[1])*textheight))

# check lengths of each scale bar before adding to map
pointDistance(p1=c(map1_scalexcoords[1],map1_scaleycoords[1]),
              p2=c(map1_scalexcoords[2],map1_scaleycoords[2]),
              lonlat=TRUE)
pointDistance(p1=c(map2_scalexcoords[1],map2_scaleycoords[1]),
              p2=c(map2_scalexcoords[2],map2_scaleycoords[2]),
              lonlat=TRUE)
pointDistance(p1=c(map3_scalexcoords[1],map3_scaleycoords[1]),
              p2=c(map3_scalexcoords[2],map3_scaleycoords[2]),
              lonlat=TRUE)
pointDistance(p1=c(map4_scalexcoords[1],map4_scaleycoords[1]),
              p2=c(map4_scalexcoords[2],map4_scaleycoords[2]),
              lonlat=TRUE)
pointDistance(p1=c(map5_scalexcoords[1],map5_scaleycoords[1]),
              p2=c(map5_scalexcoords[2],map5_scaleycoords[2]),
              lonlat=TRUE)
pointDistance(p1=c(map6_scalexcoords[1],map6_scaleycoords[1]),
              p2=c(map6_scalexcoords[2],map6_scaleycoords[2]),
              lonlat=TRUE)

# add scale bar and dam location to each map, use one legend
map1<- map1 + theme(legend.position = "none") +
  geom_segment(x = -96.135058, y = 33.646397, xend = -96.135981, yend = 33.655287,lwd=2,color=gray(0.6))+ # dam location
  geom_segment(x = map1_scalexcoords[1],y=map1_scaleycoords[2],xend=map1_scalexcoords[2],yend=map1_scaleycoords[2],lwd=1.5,color=gray(1)) #+ # scale bar part
#  annotate("text",x = mean(map1_scalexcoords), y = map1_scaleycoords[3],color=gray(1),label = "4 km") # label for scale bar
map2<- map2 + #theme(legend.position = "none") +
  geom_segment(x = -97.1910, y = 31.576152, xend = -97.2136, yend = 31.594505,lwd=2,color=gray(0.6))+
  geom_segment(x = map2_scalexcoords[1],y=map2_scaleycoords[2],xend=map2_scalexcoords[2],yend=map2_scaleycoords[2],lwd=1.5,color=gray(1)) #+
#  annotate("text",x = mean(map2_scalexcoords), y = map2_scaleycoords[3],color=gray(1),label = "4 km")
map3<- map3 + theme(legend.position = "none") +
  geom_segment(x = -99.0031, y = 31.8352, xend = -99.0002, yend = 31.8421,lwd=2,color=gray(0.6))+
  geom_segment(x = map3_scalexcoords[1],y=map3_scaleycoords[2],xend=map3_scalexcoords[2],yend=map3_scaleycoords[2],lwd=1.5,color=gray(1)) #+
#  annotate("text",x = mean(map3_scalexcoords), y = map3_scaleycoords[3],color=gray(1),label = "4 km")
map4<- map4 + theme(legend.position = "none") +
  geom_segment(x = -99.6684, y = 31.4970, xend = -99.6640, yend = 31.5025,lwd=2,color=gray(0.6))+
  geom_segment(x = map4_scalexcoords[1],y=map4_scaleycoords[2],xend=map4_scalexcoords[2],yend=map4_scaleycoords[2],lwd=1.5,color=gray(1)) #+
#  annotate("text",x = mean(map4_scalexcoords), y = map4_scaleycoords[3],color=gray(1),label = "4 km")
map5<- map5 + theme(legend.position = "none") +
  geom_segment(x = -103.9121, y = 31.8962, xend = -103.9080, yend = 31.9069,lwd=2,color=gray(0.6))+
  geom_segment(x = map5_scalexcoords[1],y=map5_scaleycoords[2],xend=map5_scalexcoords[2],yend=map5_scaleycoords[2],lwd=1.5,color=gray(1)) #+
#  annotate("text",x = mean(map5_scalexcoords), y = map5_scaleycoords[3],color=gray(1),label = "4 km")
map6<- map6 + theme(legend.position = "none") +
  geom_segment(x = -98.3539, y = 33.7653, xend = -98.3753, yend = 33.7657,lwd=2,color=gray(0.6))+
  geom_segment(x = map6_scalexcoords[1],y=map6_scaleycoords[2],xend=map6_scalexcoords[2],yend=map6_scaleycoords[2],lwd=1.5,color=gray(1)) #+
#  annotate("text",x = mean(map6_scalexcoords), y = map6_scaleycoords[3],color=gray(1),label = "4 km")


# Combine plots into a single layout with grid.arrange
# Does this panel order make sense? 
combo_turb_lomarg<-grid.arrange(map5, map4, map3, map6, map2, map1, ncol = 3)#, padding = unit(c(-1, -1, -1, -1), "cm"))
ggsave("combo_turb_update.png", combo_turb_lomarg,width=9,height=8)



