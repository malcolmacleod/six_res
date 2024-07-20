rm(list = ls())

# working with shpfiles clipped in qgis
require(data.table)
require(dplyr)
require(rstac)
require(terra)
require(mapview)
require(httr)
require(Metrics)
require(geodrawr)
require(svDialogs)
require(rstac)
require(randomForest)
require(rasterVis)
require(RColorBrewer)
require(terrainr)
require(sf)
library(raster)
library(rgdal)
library(plotly)
library(tidyverse)

# reading in GEE output
# first with Lake Brownwood
bwdwl_bands <- rast("bw_dwl.tif")
#plot(bwdwl_bands$dwLehmann)

bw_df <-as.data.frame(bwdwl_bands, xy=TRUE)
bw_df_fui <- bw_df %>% mutate(dwlgroup = cut(dwLehmann, breaks =c(470, 475, 480, 485, 489,495, 509,530,549,559,564,567,568,569,570,571,573,575,577,579,581,583), right = T, labels = F))
bw_df_fui$dwlgroup=as.factor(bw_df_fui$dwlgroup)
bw_df_fui$system="brownwood"

bwp<-bw_df_fui %>%  ggplot(aes(dwLehmann, fill = dwlgroup)) + geom_histogram()+
  scale_fill_manual(values = c("1" = "#2158bc","2" = "#316dc5","3" = "#327cbb","4" = "#4b80a0","5" = "#568f96","6" = "#6d9298", "7" = "#698c86","8" = "#759e72",
                               "9" = "#7ba654","10" = "#7dae38","11" = "#94b660","12" = "#94b660", "13" = "#a5bc76", "14" = "#aab86d","15" = "#adb55f","16" = "#a8a965",
                               "17" = "#ae9f5c","18" = "#b3a053","19" = "#af8a44","20" = "#a46905","21" = "#9f4d04")) + theme_minimal()+ 
  xlim(540, 583) +
  labs(title = "Brownwood",
       x = "Dominant Wavelength (nm)",y = "Count",fill = "Forel-Ule Index")+
  theme(plot.title = element_text(size = 15,hjust=0.5),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.title = element_text(size = 13))
ggsave("bw_dwl_histo.png", bwp)
ggplotly(bwp)

# 1& = 473 (class 8-20), .05%=236 (class6-20)

# next with oh ivie
ivdwl_bands <- rast("iv_dwl.tif")
plot(ivdwl_bands$dwLehmann)

iv_df <-as.data.frame(ivdwl_bands, xy=TRUE)
iv_df_fui <- iv_df %>% mutate(dwlgroup = cut(dwLehmann, breaks =c(470, 475, 480, 485, 489,495, 509,530,549,559,564,567,568,569,570,571,573,575,577,579,581,583), right = T, labels = F))
iv_df_fui$dwlgroup=as.factor(iv_df_fui$dwlgroup)
iv_df_fui$system="ivie"

ivp<-iv_df_fui %>%  ggplot(aes(dwLehmann, fill = dwlgroup)) + geom_histogram()+
  scale_fill_manual(values = c("1" = "#2158bc","2" = "#316dc5","3" = "#327cbb","4" = "#4b80a0","5" = "#568f96","6" = "#6d9298", "7" = "#698c86","8" = "#759e72",
                               "9" = "#7ba654","10" = "#7dae38","11" = "#94b660","12" = "#94b660", "13" = "#a5bc76", "14" = "#aab86d","15" = "#adb55f","16" = "#a8a965",
                               "17" = "#ae9f5c","18" = "#b3a053","19" = "#af8a44","20" = "#a46905","21" = "#9f4d04")) + theme_minimal()+ 
  labs(title = "O.H. Ivie",
       x = "Dominant Wavelength (nm)",y = "Count",fill = "Forel-Ule Index")+
  xlim(540, 583) +
  theme(plot.title = element_text(size = 15,hjust=0.5),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.title = element_text(size = 13))
ggsave("iv_dwl_histo.png", ivp)
ggplotly(ivp)



# 1% = 905 (9-19); 0.05% = 452 (6-20)

# now with redbluff
rbdwl_bands <- rast("rb_dwl.tif")
plot(rbdwl_bands$dwLehmann)

rb_df <-as.data.frame(rbdwl_bands, xy=TRUE)
rb_df_fui <- rb_df %>% mutate(dwlgroup = cut(dwLehmann, breaks =c(470, 475, 480, 485, 489,495, 509,530,549,559,564,567,568,569,570,571,573,575,577,579,581,583), right = T, labels = F))
rb_df_fui$dwlgroup=as.factor(rb_df_fui$dwlgroup)
rb_df_fui$system="redbluff"

rbp<-rb_df_fui %>%  ggplot(aes(dwLehmann, fill = dwlgroup)) + geom_histogram()+
  scale_fill_manual(values = c("1" = "#2158bc","2" = "#316dc5","3" = "#327cbb","4" = "#4b80a0","5" = "#568f96","6" = "#6d9298", "7" = "#698c86","8" = "#759e72",
                               "9" = "#7ba654","10" = "#7dae38","11" = "#94b660","12" = "#94b660", "13" = "#a5bc76", "14" = "#aab86d","15" = "#adb55f","16" = "#a8a965",
                               "17" = "#ae9f5c","18" = "#b3a053","19" = "#af8a44","20" = "#a46905","21" = "#9f4d04")) + theme_minimal()+ 
  labs(title = "Red Bluff",
       x = "Dominant Wavelength (nm)",y = "Count",fill = "Forel-Ule Index")+
  xlim(540, 583) +
  theme(plot.title = element_text(size = 15,hjust=0.5),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.title = element_text(size = 13))
ggsave("rb_dwl_histo.png", rbp)
ggplotly(rbp)

# 1% = 494.98 (9-19); .05% = 247.5 (9-20)

# now with bonham
bndwl_bands <- rast("bon_dwl.tif")
plot(bndwl_bands$dwLehmann)

bn_df <-as.data.frame(bndwl_bands, xy=TRUE)
bn_df_fui <- bn_df %>% mutate(dwlgroup = cut(dwLehmann, breaks =c(470, 475, 480, 485, 489,495, 509,530,549,559,564,567,568,569,570,571,573,575,577,579,581,583), right = T, labels = F))
bn_df_fui$dwlgroup=as.factor(bn_df_fui$dwlgroup)
bn_df_fui$system="bonham"

bnp<-bn_df_fui %>%  ggplot(aes(dwLehmann, fill = dwlgroup)) + geom_histogram()+
  scale_fill_manual(values = c("1" = "#2158bc","2" = "#316dc5","3" = "#327cbb","4" = "#4b80a0","5" = "#568f96","6" = "#6d9298", "7" = "#698c86","8" = "#759e72",
                               "9" = "#7ba654","10" = "#7dae38","11" = "#94b660","12" = "#94b660", "13" = "#a5bc76", "14" = "#aab86d","15" = "#adb55f","16" = "#a8a965",
                               "17" = "#ae9f5c","18" = "#b3a053","19" = "#af8a44","20" = "#a46905","21" = "#9f4d04")) + theme_minimal()+ 
  xlim(470, 583) +
  labs(title = "Bonham",
       x = "Dominant Wavelength (nm)",y = "Count",fill = "Forel-Ule Index")+
  theme(plot.title = element_text(size = 15,hjust=0.5),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.title = element_text(size = 13))
ggsave("bon_dwl_histo.png", bnp)
ggplotly(bnp)

#1% = 73 (4-11); .05% = (4-15)

# now with waco
lwdwl_bands <- rast("lw_dwl.tif")
plot(lwdwl_bands$dwLehmann)

lw_df <-as.data.frame(lwdwl_bands, xy=TRUE)
lw_df_fui <- lw_df %>% mutate(dwlgroup = cut(dwLehmann, breaks =c(470, 475, 480, 485, 489,495, 509,530,549,559,564,567,568,569,570,571,573,575,577,579,581,583), right = T, labels = F))
lw_df_fui$dwlgroup=as.factor(lw_df_fui$dwlgroup)
lw_df_fui$system="waco"

lwp<-lw_df_fui %>%  ggplot(aes(dwLehmann, fill = dwlgroup)) + geom_histogram()+
  scale_fill_manual(values = c("1" = "#2158bc","2" = "#316dc5","3" = "#327cbb","4" = "#4b80a0","5" = "#568f96","6" = "#6d9298", "7" = "#698c86","8" = "#759e72",
                               "9" = "#7ba654","10" = "#7dae38","11" = "#94b660","12" = "#94b660", "13" = "#a5bc76", "14" = "#aab86d","15" = "#adb55f","16" = "#a8a965",
                               "17" = "#ae9f5c","18" = "#b3a053","19" = "#af8a44","20" = "#a46905","21" = "#9f4d04")) + theme_minimal()+ 
  xlim(555, 577) +
  labs(title = "Waco",
       x = "Dominant Wavelength (nm)",y = "Count",fill = "Forel-Ule Index")+
  theme(plot.title = element_text(size = 15,hjust=0.5),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.title = element_text(size = 13))
ggsave("lw_dwl_histo.png", lwp)
ggplotly(lwp)

# 1% = 608 (9-17); .05% = 304 (6-18)

# lastly with arrowhead
ahdwl_bands <- rast("ah_dwl.tif")
plot(ahdwl_bands$dwLehmann)

ah_df <-as.data.frame(ahdwl_bands, xy=TRUE)
ah_df_fui <- ah_df %>% mutate(dwlgroup = cut(dwLehmann, breaks =c(470, 475, 480, 485, 489,495, 509,530,549,559,564,567,568,569,570,571,573,575,577,579,581,583), right = T, labels = F))
ah_df_fui$dwlgroup=as.factor(ah_df_fui$dwlgroup)
ah_df_fui$system="arrowhead"

ahp<-ah_df_fui %>%  ggplot(aes(dwLehmann, fill = dwlgroup)) + geom_histogram()+
  scale_fill_manual(values = c("1" = "#2158bc","2" = "#316dc5","3" = "#327cbb","4" = "#4b80a0","5" = "#568f96","6" = "#6d9298", "7" = "#698c86","8" = "#759e72",
                               "9" = "#7ba654","10" = "#7dae38","11" = "#94b660","12" = "#94b660", "13" = "#a5bc76", "14" = "#aab86d","15" = "#adb55f","16" = "#a8a965",
                               "17" = "#ae9f5c","18" = "#b3a053","19" = "#af8a44","20" = "#a46905","21" = "#9f4d04")) + theme_minimal()+ 
  xlim(470, 583) +
  labs(title = "Arrowhead",
       x = "Dominant Wavelength (nm)",y = "Count",fill = "Forel-Ule Index")+
  theme(plot.title = element_text(size = 15,hjust=0.5),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.title = element_text(size = 13))
ggsave("ah_dwl_histo.png", ahp)
ggplotly(ahp)

# 1% = 1150 (3-21); .05% = 575.11 (2-21)

dwl_all6<-grid.arrange(ncol=3,
                        bnp,lwp,
                       bwp,ivp,
                       rbp,ahp,
                        padding = unit(8.0, "line"))
ggsave(filename = "combo_dwl.png", plot = dwl_all6, , width = 18, height = 11) 

######################################################################################
df_bind <- rbind(bn_df_fui,lw_df_fui,ah_df_fui,bw_df_fui,iv_df_fui)

# Define the UTM zone and datum (replace with your specific UTM EPSG code)
# UTM zone 14N, datum WGS84, RedBluff is the only one that is zone 13N
utm_crs <- st_crs(32614)  # EPSG code for UTM zone 33N# Convert data frame to sf object
sf_df <- st_as_sf(df_bind, coords = c("x", "y"), crs = utm_crs)

latlon_crs <- st_crs(4326)  # EPSG code for WGS84
sf_df_latlon <- st_transform(sf_df, crs = latlon_crs)
# Convert back to a regular data frame
df_latlon <- as.data.frame(st_coordinates(sf_df_latlon))
df_latlon$dwl <- sf_df_latlon$dwLehmann  
df_latlon$fui<-sf_df_latlon$dwlgroup
df_latlon$system<-sf_df_latlon$system

utm_crs_rb<-st_crs(32613)
rb_df_sf<-st_as_sf(rb_df_fui, coords = c("x", "y"), crs = utm_crs_rb)
rblatlon_crs <- st_crs(4326)  # EPSG code for WGS84
rbsf_df_latlon <- st_transform(rb_df_sf, crs = rblatlon_crs)
# Convert back to a regular data frame
rb_sf <- as.data.frame(st_coordinates(rbsf_df_latlon))
rb_sf$dwl <- rbsf_df_latlon$dwLehmann  
rb_sf$fui<-rbsf_df_latlon$dwlgroup
rb_sf$system<-rbsf_df_latlon$system

#####################################################################################################################################
# separating lake data between river arm and main body  

s2_w<-df_latlon %>% filter(system=="waco")
s2_a<-df_latlon %>% filter(system=="arrowhead")
s2_bw<-df_latlon %>% filter(system=="brownwood")
s2_iv<-df_latlon %>% filter(system=="ivie")
s2_bn<-df_latlon %>% filter(system=="bonham")
s2_rb<-df_latlon %>% filter(system=="redbluff")

## defining bboxes for waco ##
#wmbbox <- st_bbox(c(xmin = -97.25316, ymin = 31.539846, xmax = -97.1882756, ymax = 31.589283782), crs = st_crs(s2_w))

wmb <- s2_w[s2_w$Y >= 31.530244 & s2_w$Y <= 31.6095146 & s2_w$X >= -97.2520428 & s2_w$X <= -97.1809749,]
wna<-s2_w[s2_w$Y >= 31.5821285 & s2_w$Y <= 31.6157569 & s2_w$X >= -97.308004 & s2_w$X <= -97.2523861,]
wsa<-s2_w[s2_w$Y >= 31.4886804 & s2_w$Y <= 31.5311218 & s2_w$X >= -97.265776 & s2_w$X <= -97.2266369,]

wmb$zone="body"
wna$zone="arm"
wsa$zone="arm"

waco_bind<-rbind(wmb,wna,wsa)
#write_csv(waco_bind, "waco_zone.csv")

## defining bboxes for ah ##
amb<-s2_a[s2_a$Y >= 33.6933513 & s2_a$Y <= 33.7744378 & s2_a$X >= -98.41832729 & s2_a$X <= -98.3015976,]
ara<-s2_a[s2_a$Y >= 33.6339153 & s2_a$Y <= 33.6936369 & s2_a$X >= -98.4667358 & s2_a$X <= -98.351036,]

amb$zone="body"
ara$zone="arm"

ah_bind<-rbind(amb,ara)

## defining bboxes for bw ##
bwmb<-s2_bw[s2_bw$Y >= 31.7955865 & s2_bw$Y <= 31.865884 & s2_bw$X >= -99.0787206 & s2_bw$X <= -98.9753804,]
bwna<-s2_bw[s2_bw$Y >= 31.8661754 & s2_bw$Y <= 31.928552 & s2_bw$X >= -99.0583197 & s2_bw$X <= -99.0109412,]
bwsa<-s2_bw[s2_bw$Y >= 31.8200945 & s2_bw$Y <= 31.857136 & s2_bw$X >= -99.1343388 & s2_bw$X <= -99.0787206,]

bwmb$zone="body"
bwna$zone="arm"
bwsa$zone="arm"

bw_bind<-rbind(bwmb,bwna,bwsa)

## defining bboxes for iv ##
imb<-s2_iv[s2_iv$Y >= 31.479319 & s2_iv$Y <= 31.5557074 & s2_iv$X >= -99.7077124 & s2_iv$X <= -99.574503,]
ira<-s2_iv[s2_iv$Y >= 31.5559999 & s2_iv$Y <= 31.6256021 & s2_iv$X >= -99.8392049 & s2_iv$X <= -99.659991,]

imb$zone="body"
ira$zone="arm"

iv_bind<-rbind(imb,ira)

## defining bboxes for bon ##
bnmb1<-s2_bn[s2_bn$Y >= 33.6512746 & s2_bn$Y <= 33.6556328 & s2_bn$X >= -96.1523812 & s2_bn$X <= -96.13573,]
bnmb2<-s2_bn[s2_bn$Y >= 33.6442593 & s2_bn$Y <= 33.6515473 & s2_bn$X >= -96.1636449 & s2_bn$X <= -96.1349774,]
bnmb1$zone="body"
bnmb2$zone="body"


bnrane<-s2_bn[s2_bn$Y >= 33.6556197 & s2_bn$Y <= 33.67240739 & s2_bn$X >= -96.15240104 & s2_bn$X <= -96.13789,]
bnranw<-s2_bn[s2_bn$Y >= 33.6516187 & s2_bn$Y <= 33.66576413 & s2_bn$X >= -96.1694813 & s2_bn$X <= -96.152401,]
bnraw<-s2_bn[s2_bn$Y >= 33.63775681 & s2_bn$Y <= 33.65047555 & s2_bn$X >= -96.17806441 & s2_bn$X <= -96.163731,] 
bnras<-s2_bn[s2_bn$Y >= 33.6397708 & s2_bn$Y <= 33.644201 & s2_bn$X >= -96.1537165 & s2_bn$X <= -96.13637867,]
bnrane$zone="arm"
bnranw$zone="arm"
bnraw$zone="arm"
bnras$zone="arm"

bn_bind<-rbind(bnmb1,bnmb2,bnrane,bnranw,bnraw,bnras)

## defining bboxes for rb ##
rmb<-rb_sf[rb_sf$Y >= 31.895208 & rb_sf$Y <= 31.9744585 & rb_sf$X >= -103.9684497 & rb_sf$X <= -103.8946353,]
rra<-rb_sf[rb_sf$Y >= 31.974749 & rb_sf$Y <= 32.0443283 & rb_sf$X >= -104.01239497 & rb_sf$X <= -103.9324008,]

rmb$zone="body"
rra$zone="arm"

rb_bind<-rbind(rmb,rra)

# this bbind product has the correct lat/lon and dwl and zone
bb_bind <- rbind(bn_bind,waco_bind,ah_bind,bw_bind,iv_bind,rb_bind)
bb_bind$dwl<-as.numeric(bb_bind$dwl)

bb_nona<-bb_bind %>% drop_na()

# means, CI, and range for variables by zone
# defining confidence level
confidence_level <- 0.95
alpha <- 1-confidence_level

# summary stats for them all
systemstats_dwl<- bb_nona %>% group_by(system) %>% 
  summarise(n=n(),mean=mean(dwl),sd=sd(dwl),se=sd/sqrt(n),
            lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
            upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
            min=min(dwl), max=max(dwl),
            margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))

zonestats_dwl<- bb_nona %>% group_by(zone) %>% 
  summarise(n=n(),mean=mean(dwl),sd=sd(dwl),se=sd/sqrt(n),
            lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
            upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
            min=min(dwl), max=max(dwl),
            margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))













# separating arm from body by bounding box
# creating bounding boxes
wmb <- extent(c(31.5477465,31.589283782,-97.2476704,-97.1882756))
wna<-extent(c(31.5886988,31.6123845,-97.2926457,-97.236684))
wsa<-extent(c(31.50619073,31.5351646,-97.2541936,-97.219518))

amb<-extent(c(33.7094465,33.77168417,-98.409005,-98.306695))
ara<-extent(c(33.63058548,33.6991644,-98.4687434,-98.328667))

bwmb<-extent(c(31.81675433,31.8587545,-99.0424825,-98.982401))
bwna<-extent(c(31.86079574,31.90889642,-99.0569021,-98.9964772))
bwsa<-extent(c(31.79749784,31.84825629,-99.1386129,-99.043169))

imb<-extent(c(31.47982501,31.545356,-99.7080943,-99.633277))
ira<-extent(c(31.5582608,31.6068112,-99.771952,-99.660029))

rmb<-extent(c(31.47982501,31.9580921,-103.962899,-103.89335))
rra<-extent(c(31.96781824,32.04118549,-104.001352,-103.932687))

# Check if CRS transformation is needed
raster_crs <- crs(lwdwl_bands)
wmb <- CRS("+proj=longlat +datum=WGS84")

# cropping raster to bounding box
cwmb<-crop(lwdwl_bands, wmb)
