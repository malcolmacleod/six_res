# A script for with working code for chapter 1 data analysis
# This combines from scripts to include working with FLAMe data, CRASR samples and ee.points transect queries 

# Note that everything with GEE TIF output is in its own script "dwl6lakes.R"
# Processing ACOLITE output is in "acolite6lakes.R"

# clear workspace
rm(list=ls())


# 1. load essential packages

library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(OpenStreetMap)
library(colorscience)
library(sf)
library(raster)
require(rasterVis)
require(RColorBrewer)
require(terrainr)
require(rstac)
require(terra)
library(nlme)
library(anytime)
library(patchwork)
library(ggforce)

# 2. Reading in data

# Reading the combined FLAMe and Sentinel-2 data for pts along boat path with essential data variables:
# lat, lon, DWL, predicted SD, NDTI, turbidity
# this file for "s2flame_znpts" is all the points filtered to "no flag" and unique lat/lons (r2=0.59, n = 28328)
# applying SCL = 6 (r2 = 0.6, n = 28227), after removing Waco (which should be done) r2=0.68, n =26639
# filtering to within visible range (r2=0.7, n=26593 or r2 = 0.72 not log10 transformed)
# applied the band2 threshold > 300 (r2 = 0.75, n=18677)

# first reading in no flag unique data with zones specified (narm etc)
s2flame_znpts<-read_csv("znspcf_s2flm.csv")%>% 
  dplyr::rename(group=zone)%>%
  mutate(zone=case_when(group %in% c("north_arm","south_arm","arm") ~ 
                          "arm",group == "body" ~ "body"))

s2f_znspc_filter<- s2flame_znpts%>%
  filter(SCL == "6")%>% 
  filter(!system=="waco")%>% filter(dwl<=583 & dwl>=475)

s2f_avg<-s2flame_znpts%>% group_by(system, group) %>% dplyr::summarise(mean_turb = mean(turb),
                                                                   mean_secchi = mean(secchi),
                                                                   mean_dwl = mean(dwl),
                                                                   mean_ndti = mean(ndti),
                                                                   mean_tsi=mean(tsi_sd),
                                                                   mean_fdom = mean(fdom),
                                                                   mean_chla = mean(chl_rfu))

s2f_avg_zn<- s2f_avg %>% mutate(zone=case_when(group %in% c("north_arm","south_arm","arm") ~ "arm",
                                               group == "body" ~ "body"))

ggplot(s2f_znspc_filter,aes(log10(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()


# # reading in boat path points with only 2 zones - mb and ra
# s2flame_znpts<-read_csv("s2flm6lakes_nfu_4mzone.csv")  
# 
# s2flm_s2c_filter<-read_csv("six4m_zone_b2thresh.csv") %>%
#   filter(SCL == "6")%>% 
#   filter(!system=="waco")
# 
# s2flm_allfilter<-s2flame_znpts  %>% filter(SCL == "6")%>% 
#   filter(!system=="waco") %>% filter(dwl<=583 & dwl>=475)
#   #filter(B2>300) 
# 
# s2flm_iv<-s2flm_allfilter %>% filter(system=="ivie") %>% filter(dwl<584)
# 
# ggplot(s2flm_iv,aes(turb, ndti)) + 
#   geom_point()
# 
# 
# 
# ggplot(s2flm_s2c_filter,aes(log10(turb), ndti)) + 
#   geom_point() + 
#   geom_smooth(method = "lm", se=FALSE) +
#   stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()
# 
# ggplot(s2flm_allfilter,aes(turb, ndti)) + 
#   geom_point() + 
#   geom_smooth(method = "lm", se=FALSE) +
#   stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()
# 
# s2flame_znpts %>% dplyr::group_by(system) %>% dplyr::summarise(count = n())
# s2flm_allfilter %>% dplyr::group_by(system) %>% dplyr::summarise(count = n())

# reading in the Sentinel-2 Cloud Probability masked sen2cor data along boat paths. n = 27620
s2flm_nc<-read_csv("s2cloudprob_flame_pts.csv") 

s2flm_ncfilter<-s2flm_nc  %>%
  filter(!system=="waco")

ggplot(s2flm_ncfilter,aes(turb, ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()

s2flm_nc %>% dplyr::group_by(system) %>% dplyr::summarise(count = n())

# Reading in the ACOLITE L2 Water Product data along the boat paths. n=25374
s2flm_l2w<-read_csv("l2w6lakes.csv")

s2flm_l2_scale<-s2flm_l2w %>% mutate(B2=rrs_blue*100000,
                                     B3=rrs_green*100000,
                                     B4=rrs_red*100000)

s2flm_l2w_nona<-s2flm_l2w %>% drop_na() # n = 14727

excluded_na <- s2flm_l2w %>%
  filter(is.na(ndti))

ggplot(s2flm_l2w_nona,aes(turb, ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()

s2flm_l2w %>% dplyr::group_by(system) %>% dplyr::summarise(count = n())
s2flm_l2w_nona %>% dplyr::group_by(system) %>% dplyr::summarise(count = n())

# excluding the 10647 pts from 
small_sf <- st_as_sf(small_data, coords = c("lon", "lat"), crs = 4326)

# Reading in data for distance from dam transects queried in ee.points Shiny
# inlcudes all output data plus NDTI and normalized DFD
dfd_transect<- read_csv("dfd_transect.csv")

# Reading in sample to sensor validation data (using method from AH paper) 
# refer to "station_YSI.R" for making of sausage
sample_ysi<-read_csv("valmeans_secchi.csv")

# Reading in sentinel-2 data for the station lat lons
# note that this calculates dwl/ndti for hi-qual images for AH (7/28) and Waco(7/25)
s2crasr<-read_csv("s2crasr_merge.csv")

# Reading full system Sen2Cor data for DWL and NDTI
dwl6lake_all_nozone<-read_csv("dwl6lake.csv")

dwl6lake_zone<-read_csv("dwl6lake_zn_nona.csv")

# reading in full system sen2cor using the ee.ImageCollection("COPERNICUS/S2_CLOUD_PROBABILITY") approach
# this is MUCH more effective cloud mask than the previous approach "masks2clouds". dwlgroup is factored to viz
dwl_cloudmask<-read_csv("dwl6lake_truecloudmask.csv")
dwl_cloudmask$dwlgroup<-factor(dwl_cloudmask$dwlgroup)

# reading in 2022 NLA water chem data
nla2022<-read_csv("nla22_waterchem_wide.csv") %>% clean_names()
nla_ourlakes<-nla2022 %>% filter(site_id %in% c("NLA22_TX-10002","NLA22_TX-10038",
                                             "NLA22_TX-10030", "NLA22_TX-10001","NLA22_TX-10054"))

# Arrowhead: site_id = NLA22_TX-10002, unique_id =NLA_TX-10216, gnis_name = `Lake Arrowhead`
# Bonham: site_id = NLA22_TX-10038, unique_id = NLA_TX-10332, gnis_name = 'Lake Bonham'
# Waco: site_id = NLA22_TX-10030, unique_id = NLA_TX-10324, gnis_name = 'Waco Lake'
# Red Bluff: site_id = NLA22_TX-10001, unique_id = NLA_TX-10215, gnis_name = 'Red Bluff Reservoir'
# O.H. Ivie: not listed by name but by process of elimination from county I think it is site_id = NLA22_TX-10054, unique_id=NLA_TX-10243
# Brownwood: not found

# 4. Plotting data

ggplot(s2flame_znpts,aes(log10(turb), ndti,color=system)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()
s2flame_znpts %>% dplyr::group_by(system) %>% dplyr::summarise(count = n())

ggplot(s2flm_l2w,aes(log10(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()
s2flm_l2w %>% dplyr::group_by(system) %>% dplyr::summarise(count = n())

s2flm_nfunq %>% 
  ggplot(aes(turb, ndti, color=system)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + 
  coord_cartesian(ylim = range(s2flm_nfunq$ndti))+theme_bw()
s2flm_nfunq %>% dplyr::group_by(system) %>% dplyr::summarise(count = n())

sysndti<-ggplot(s2flm_nfunq,aes(turb, ndti, color=system)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + 
  coord_cartesian(ylim = range(s2flm_nfunq$ndti))+ theme_bw()
ggplotly(sysndti)

s2flm_nfunq %>% dplyr::group_by(system) %>% dplyr::summarise(count = n())


# First we have log10(turb) ~ NDTI plotted along the boat path for all lakes using sen2cor
s2c_turb_ndti_plot<-ggplot(s2f_znspc_filter,aes(log10(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  xlab("log10(Turbidity)") + ylab("Normalized Difference Turbidity Index") + theme_bw()

# plotting log10(turb) ~ NDTI along boat path for all lakes using ACOLITE L2W
l2w_turb_ndti_plot<-s2flm_l2w %>%  ggplot(aes(log10(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  xlab("log10(Turbidity)") + ylab("Normalized Difference Turbidity Index") + theme_bw()

boatpath_turbndti_combo<- s2c_turb_ndti_plot + l2w_turb_ndti_plot +
  plot_annotation(tag_levels = 'A') 
ggsave("boatpath_ndt_turb_combo.png", boatpath_turbndti_combo, width=14, height = 8)

# plotting NDTI x norm DFD
line_ndti<- dfd_transect %>% ggplot(aes(norm_dist, ndti, color = system)) + geom_line(size=1.25) + 
  xlab("Normalized Distance from Dam") + ylab("Normalized Difference Turbidity Index") +  
  theme_bw()+scale_color_manual(values = viridis::viridis(6, option = "C"),
                                labels=c("bonham" = "Bonham", "waco"="Waco","brownwood"="Brownwood",
                                         "ohivie"="O.H. Ivie","redbluff"="Red Bluff","arrowhead"="Arrowhead"),
                                name = "System")+
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))

# plotting DWL x norm DFD
line_dwl_p<- dfd_transect %>% ggplot(aes(norm_dist, dwl, color = system)) + geom_line(size=1.25) + 
  xlab("Normalized Distance from Dam") + ylab("Dominant Wavelength") + theme_bw() +
  scale_color_manual(values = viridis::viridis(6, option = "C"),
                     labels=c("bonham" = "Bonham", "waco"="Waco","brownwood"="Brownwood",
                              "ohivie"="O.H. Ivie","redbluff"="Red Bluff","arrowhead"="Arrowhead"),
                     name = "System")+
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))

# plotting relationship between grab sample data and mean sensor values
lmt_valmeans<-lm(turb_ysi_m~turb_lab,data=sample_ysi) 
summary(lmt_valmeans) # r2 = 0.9893

# plotting sensor~sample for all 6 lakes - same lm as above
all6t_r2 <-ggplot(sample_ysi,aes(y=turb_ysi_m,x=turb_lab)) + 
  geom_point(aes(color = system)) + 
  geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(after_stat(rr.label))), 
               parse = TRUE, label.x.npc = "left", size = 5, rr.digits = 3) + ylab("Sensor Turbidity (NTU)")+
  xlab("Sample Turbidity (NTU)")+ theme_bw()
# code for plots of individual lakes and their stations + all combined are found in "station_ysi.R" 

# plotting DWL histograms and maps
dwl_cloudmask$system <- factor(dwl_cloudmask$system, levels = 
                           c("Red Bluff","OH Ivie","Brownwood","Arrowhead","Waco","Bonham"))
# making a FUI palette
fui_palette<-c("1" = "#2158bc","2" = "#316dc5","3" = "#327cbb","4" = "#4b80a0",
               "5" = "#568f96","6" = "#6d9298", "7" = "#698c86","8" = "#759e72",
               "9" = "#7ba654","10" = "#7dae38","11" = "#94b660","12" = "#94b660", 
               "13" = "#a5bc76", "14" = "#aab86d","15" = "#adb55f","16" = "#a8a965",
               "17" = "#ae9f5c","18" = "#b3a053","19" = "#af8a44","20" = "#a46905","21" = "#9f4d04")

fui_palette_cont<-c("#2158bc","#316dc5","#327cbb", "#4b80a0",
                    "#568f96", "#6d9298", "#698c86","#759e72",
                    "#7ba654", "#7dae38", "#94b660", "#94b660", 
                    "#a5bc76", "#aab86d","#adb55f", "#a8a965",
                    "#ae9f5c","#b3a053", "#af8a44","#a46905", "#9f4d04")

dwl_hist<-ggplot(dwl_cloudmask, aes(x = dwLehmann, fill = dwlgroup)) +
  geom_histogram(aes(y = after_stat(count / tapply(count, PANEL, sum)[PANEL]), fill = ..x..), 
                 color = "black", bins = 21) +  # Normalize bin heights within each facet
  scale_fill_gradientn(colors = fui_palette,
                       values = scales::rescale(c(475, 480, 485,489,495,509,530,549,559,564,567,568,569,570,571,573,575,577,579,581,583))) +  # Apply the palette to the fill aesthetic
  scale_x_continuous(breaks = c(475, 500, 525, 550, 575)) + # Custom x-axis tick labels
  coord_cartesian(xlim = c(520, 583))+
  facet_wrap(~system) + 
  ylab("Proportion of surface area") + 
  xlab("Dominant wavelength (nm)") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("dwl_histo_xlim.png", dwl_hist)

dwl_hist_att2<-ggplot(dwl_cloudmask, aes(x = dwLehmann, fill = dwlgroup)) +
  geom_histogram(aes(y = after_stat(count / tapply(count, PANEL, sum)[PANEL]), fill = ..x..)) +  # Normalize bin heights within each facet                  color = "black", bins = 21 
  scale_fill_gradientn(colors = fui_palette,
                       values = scales::rescale(c(475, 480, 485,489,495,509,530,549,559,564,567,568,569,570,571,573,575,577,579,581,583)))+ # Custom x-axis tick labels
  coord_cartesian(xlim = c(520, 583))+
  facet_wrap(~system) + 
  ylab("Proportion of surface area") + 
  xlab("Dominant wavelength (nm)") +
  theme_bw() +
  theme(legend.position = "none")
ggsave("dwl_histo_nobin.png",dwl_hist_att2)

ggplotly(dwl_hist)

# trying to make histo without the bars


all_map_facet<-dwl_cloudmask %>% ggplot(aes(x,y,color=dwlgroup)) + geom_point(shape=15, size=0.35)+
  scale_color_manual(values = c(fui_palette)) +guides(color = guide_legend(override.aes = list(size = 6))) +
  ylab("Latitude") + 
  xlab("Longitude") +
  labs(color = "Forel-Ule Scale") +
  theme_classic() +facet_wrap(~system, scales = "free")


# verifying ndti pts
map0 <- openmap(upperLeft = c(33.672, -96.18538), 
                lowerRight = c(33.63142, -96.12549),
                #                type = 'stamen-terrain',zoom=12) # stamen-terrain is deprecated
                type = 'osm',zoom=14)

ndti_map_bonham<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map0)) +
  geom_point(data = s2flame_znpts %>% filter(system=="bonham"), 
             #             size=ptsize,
             #             alpha=ptalpha,
             aes(x = lon, y = lat,color=ndti))+
  scale_color_viridis_c("ndti",direction=)+
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.6,"cm"),
    legend.position = c(0.28, 0.8))+
  guides(col=guide_colorbar(title.position = "top"))
ndti_map_bonham

map0 <- openmap(upperLeft = c(31.61010, -99.8149), 
                lowerRight = c(31.45884, -99.58825),
                #                type = 'stamen-terrain',zoom=12) # stamen-terrain is deprecated
                type = 'osm',zoom=12)

ndti_map_ivie<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map0)) +
  geom_point(data = s2flame_znpts %>% filter(system=="ivie"), 
             #             size=ptsize,
             #             alpha=ptalpha,
             aes(x = lon, y = lat,color=ndti))+
  scale_color_viridis_c("ndti",direction=)+
  xlab("")+ylab("")+
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.direction="horizontal",
    legend.title.align=0.5,
    legend.key.width = unit(0.6,"cm"),
    legend.position = c(0.28, 0.1))+
  guides(col=guide_colorbar(title.position = "top"))
ndti_map_ivie

# NDTI~turbidity from boat path points
logEstimate <- lm(ndti~log(turb),data=s2flame_znpts)
xvec <- seq(min(s2flame_znpts$turb),max(s2flame_znpts$turb),length=1000)
logpred <- predict(logEstimate,newdata=data.frame(turb=xvec))
pred <- data.frame(x = xvec, y = logpred)

# simple plot showing R2 = 0.59 for all lakes
turbndti_p<-ggplot(s2flame_znpts,aes(log(turb), ndti)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) +
  stat_regline_equation(aes(label = ..rr.label..)) + theme_bw()

# plot done in the style of AH paper
turb_ndti_plot<-ggplot(s2flame_znpts,aes(x=turb,y=ndti)) +
  geom_point(size=0.1,col=gray(0.3))+
  geom_smooth(method="lm",lty=5,col="black",se=FALSE,lwd=1)+
  geom_line(data = pred, aes(x=x, y=y),lty=1,col="black",lwd=1)+
  #xlim(c(0,99))+
  xlab("In-lake Turbidity (NTU)")+
  ylab("Satellite Turbidity (NDTI)")+
  theme_bw()

# plotting out YSI to sample
turb_stations <-ggplot(sample_ysi,aes(y=turb_ysi_m,x=turb_lab)) + 
  geom_point(aes(color = system)) + 
  geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(after_stat(rr.label))), 
               parse = TRUE, label.x.npc = "left", size = 5, rr.digits = 3) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + 
  ylab("Sensor Turbidity (NTU)")+
  xlab("Sample Turbidity (NTU)")+
  labs(color = "System") + theme_bw()


# Now it's predicted Secchi to in situ Secchi
sdd_stations<-sample_ysi %>% #filter(!system =="Lake Waco") %>% 
  ggplot(aes(y=pred_sdd_m,x=secchi)) + 
  geom_point(aes(color = system)) + 
  geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(after_stat(rr.label))), 
               parse = TRUE, label.x.npc = "left", size = 5, rr.digits = 3) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  ylab("Predicted Secchi (m)")+
  xlab("Secchi (m)")+
  labs(color = "System")+ theme_bw()

comb_sampsens_plot<- turb_stations + sdd_stations + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & theme(legend.position = "right")
ggsave("sdd_turb_stations.png", comb_sampsens_plot, width = 15, height = 7)


# 4a. Summary statistics

# defining confidence level
confidence_level <- 0.95
alpha <- 1-confidence_level

# calculating means, CI, and range for variables by zone, can also specify group_by(system)

zonestats_turb<- s2flame_znpts %>% group_by(zone) %>% 
  dplyr::summarise(n=n(),mean=mean(turb),sd=sd(turb),se=sd/sqrt(n),
            lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
            upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
            min=min(turb), max=max(turb),
            margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))


zonestats_secchi<- s2flame_znpts %>% #group_by(zone) %>% 
  dplyr::summarise(n=n(),mean=mean(secchi),sd=sd(secchi),se=sd/sqrt(n),
            lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
            upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
            min=min(secchi), max=max(secchi),
            margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))

zonestats_ndti<- s2flame_znpts %>% group_by(zone) %>% 
  dplyr::summarise(n=n(),mean=mean(ndti),sd=sd(ndti),se=sd/sqrt(n),
            lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
            upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
            min=min(ndti), max=max(ndti),
            margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))

zonestats_dwl<- s2flame_znpts %>% group_by(zone) %>% 
  dplyr::summarise(n=n(),mean=mean(dwl),sd=sd(dwl),se=sd/sqrt(n),
            lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
            upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
            min=min(dwl), max=max(dwl),
            margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))

zonestats_tsi<- s2flame_znpts %>% group_by(zone) %>% 
  dplyr::summarise(n=n(),mean=mean(tsi_sd),sd=sd(tsi_sd),se=sd/sqrt(n),
            lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
            upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
            min=min(tsi_sd), max=max(tsi_sd),
            margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))

s2flm_river<- s2flame_znpts %>% filter(zone=="arm")
s2flm_lacus<- s2flame_znpts %>% filter(zone=="body")

tsi_r<-s2flm_river$tsi_sd
tsi_l<-s2flm_lacus$tsi_sd
quantile(tsi_r, c(.05, 0.5, .95))
quantile(tsi_l, c(.05, 0.5, .95))

confidence_level <- 0.95
alpha <- 1-confidence_level

# calculating means, CI, and range for variables by zone, can also specify group_by(system)

syststats_turb<- s2flame_znpts %>% group_by(system) %>% 
  dplyr::summarise(n=n(),mean=mean(turb),sd=sd(turb),se=sd/sqrt(n),
                   lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
                   upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
                   min=min(turb), max=max(turb),
                   margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))


syststats_secchi<- s2flame_znpts %>% group_by(system) %>% 
  dplyr::summarise(n=n(),mean=mean(secchi),sd=sd(secchi),se=sd/sqrt(n),
                   lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
                   upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
                   min=min(secchi), max=max(secchi),
                   margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))

syststats_ndti<- s2flame_znpts %>% group_by(system) %>% 
  dplyr::summarise(n=n(),mean=mean(ndti),sd=sd(ndti),se=sd/sqrt(n),
                   lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
                   upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
                   min=min(ndti), max=max(ndti),
                   margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))

syststats_dwl<- s2flame_znpts %>% group_by(system) %>% 
  dplyr::summarise(n=n(),mean=mean(dwl),sd=sd(dwl),se=sd/sqrt(n),
                   lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
                   upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
                   min=min(dwl), max=max(dwl),
                   margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))

syststats_tsi<- s2flame_znpts %>% group_by(system) %>% 
  dplyr::summarise(n=n(),mean=mean(tsi_sd),sd=sd(tsi_sd),se=sd/sqrt(n),
                   lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
                   upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
                   min=min(tsi_sd), max=max(tsi_sd),
                   margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))

# note that results for DWL along boat path are not consistent with full system analysis
# same summary stats were calculated for the zone designated water pixel data from GEE output
## see "dwl5lakes.R" for processing of GEE output - histogram creation, zone designation, and stats ##
s2flame_znpts$system<- str_to_title(s2flame_znpts$system)
s2flm_allfilter$system<- str_to_title(s2flm_allfilter$system)

dwl_system_boxplot<-s2flame_znpts %>% filter(dwl>469 & dwl<584) %>% 
  ggplot(aes(x=system, y=dwl, fill =zone)) + 
  labs(fill = "Zone", labels = c("Arm", "Body")) + 
  xlab("System") + ylab("Dominant Wavelength (nm)") +
  geom_boxplot() +
  scale_fill_manual(values = c("#E7B800","#00AFBB"),
                    labels=c("arm" = "Arm", "body"="Body"))+
  theme_classic()
#ggsave("dwl_system_boxplot_filter.png",dwl_system_boxplot, width = 10, height = 7)

dwl_system_boxplot_nc<-s2flm_nc %>% filter(dwl>469 & dwl<584) %>% 
  ggplot(aes(x=system, y=dwl, fill =zone)) + 
  labs(fill = "Zone", labels = c("Arm", "Body")) + 
  xlab("System") + ylab("Dominant Wavelength (nm)") +
  geom_boxplot() +
  scale_fill_manual(values = c("#E7B800","#00AFBB"),
                    labels=c("arm" = "Arm", "body"="Body"))+
  theme_classic()

s2flm_nc %>% filter(dwl>469 & dwl<584) %>% 
  ggplot(aes(x=system, y=dwl, fill =zone)) + 
  labs(fill = "Zone", labels = c("Arm", "Body")) + 
  xlab("System") + ylab("Dominant Wavelength (nm)") +
  geom_boxplot() + geom_boxen() +
  scale_fill_manual(values = c("#E7B800","#00AFBB"),
                    labels=c("arm" = "Arm", "body"="Body"))+
  theme_classic()

ggboxplot(s2flm_nc , x = "system", y = "dwl", fill = "zone",
          palette = c("#E7B800", "#00AFBB"))

s2flame_znpts %>% filter(dwl>469 & dwl<584) %>% 
  ggplot(aes(system, dwl, color = zone)) + geom_boxplot()

# 4... run AOV on avg data
avg_aov_turb <- aov(mean_turb~zone*system,s2f_avg_zn)
summary(avg_aov_turb)
aovt_posthoc <- TukeyHSD(avg_aov_turb)
aovt_posthoc

avg_aov_secchi <- aov(mean_secchi~zone*system,s2f_avg_zn)
summary(avg_aov_secchi)
aovsecchi_posthoc <- TukeyHSD(avg_aov_secchi)
aovsecchi_posthoc

avg_aov_dwl <- aov(mean_dwl~zone*system,s2f_avg_zn)
summary(avg_aov_dwl)
aovdwl_posthoc <- TukeyHSD(avg_aov_dwl)
aovdwl_posthoc

avg_aov_ndti <- aov(mean_ndti~zone*system,s2f_avg_zn)
summary(avg_aov_ndti)
aovt_posthoc <- TukeyHSD(avg_aov_ndti)
aovt_posthoc


# 4b. AOV and Tukeys HSD post-hoc

# using the subset method (10% of data) applying SP code - full can be found in "s2flm_pts.R"

# concatenate lat/lon and remove duplicated lat/lon pairs
s2flame_znpts$latlon<-paste(s2flame_znpts$lat,s2flame_znpts$lon)
s2flame_znpts<-s2flame_znpts[-which(duplicated(s2flame_znpts$latlon)==TRUE),]
s2flame_znpts<-data.frame(s2flame_znpts)
s2flame_znpts$index<-1:length(s2flame_znpts[,1])

# use a random subset of data to avoid autocorrelation and speed up processing
# see Loken et al. 2019 https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JG005186
set.seed(1) # makes the subsetting reproducible
sampled <- sample(1:length(s2flame_znpts[,1]),size=2537) #2528) #200) # adjust subset size as needed #total is 28728
data_sampled<-s2flame_znpts[sampled,]

# repeating for s2flmnc
# concatenate lat/lon and remove duplicated lat/lon pairs
s2flm_nc$latlon<-paste(s2flm_nc$latitude,s2flm_nc$longitude)
s2flm_nc<-s2flm_nc[-which(duplicated(s2flm_nc$latlon)==TRUE),]
s2flm_nc<-data.frame(s2flm_nc)
s2flm_nc$index<-1:length(s2flm_nc[,1])

# use a random subset of data to avoid autocorrelation and speed up processing
# see Loken et al. 2019 https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JG005186
set.seed(1) # makes the subsetting reproducible
sampled <- sample(1:length(s2flm_nc[,1]),size=2537) #2528) #200) # adjust subset size as needed #total is 28728
data_sampled_nc<-s2flm_nc[sampled,]


noivorbon_aov<-s2flame_znpts %>% filter(!system=="bonham") 
noivorbon_aov<-noivorbon_aov %>% filter(!system=="ivie") 
baov<-aov(turb~zone*system, noivorbon_aov)
summary(baov)
bonaov_posthoc <- TukeyHSD(baov)

# fit new anova on random subset for turb
aov_turb <- aov(turb~zone*system,data_sampled_nc)
summary(aov_turb)
aovt_posthoc <- TukeyHSD(aov_turb)
aovt_posthoc

interaction.plot(
  x.factor = data_sampled$zone,
  trace.factor = data_sampled$system,
  response = data_sampled$turb,
  fun = median,
  ylab = "turbidity",
  xlab = "zone",
  trace.label = "system",
  col = c("#0198f9", "#f95801", "forestgreen", "yellow", "violet", "black"),
  lyt = 1,
  lwd = 3
)

interaction.plot(
  x.factor = data_sampled$system,
  trace.factor = data_sampled$zone,
  response = data_sampled$turb,
  fun = median,
  ylab = "turbidity",
  xlab = "system",
  trace.label = "zone",
  col = c("#0198f9", "#f95801"),
  lyt = 1,
  lwd = 3
)

# do diagnostics
plot(aov_turb,which=1)
plot(aov_turb,which=2)
plot(aov_turb,which=3)

# buiid data frame with residuals and fitted, to check model diagnostics
resid <- as.numeric(aov_turb$residuals)
resid_st <- as.numeric(rstudent(aov_turb)) #
fitted <-as.numeric(aov_turb$fitted.values)
value <-as.numeric(aov_turb$model$turb)
aov1_df<- data.frame(zone=aov_turb$model$zone,
                     system=aov_turb$model$system,
                     value,
                     resid,resid_st,fitted)

# this is a ggplot version of the above residual vs fitted, plot(aov1,which=1)
resid_vs_fitted<-aov1_df %>%
  ggplot(aes(x=fitted,y=resid))+
  geom_point()+
  theme_bw()
resid_vs_fitted

# do more diagnostics - acf and pacf of the fitted and residuals
png("aov1_acf_fitted.png")
acf_fitted<-acf(aov1_df$fitted)
dev.off()
png("aov1_pacf_fitted.png")
pacf_fitted<-pacf(aov1_df$fitted)
dev.off()
png("aov1_acf_resid.png")
acf_resid<-acf(aov1_df$resid)
dev.off()
png("aov1_pacf_resid.png")
pacf_resid<-pacf(aov1_df$resid)
dev.off()

# fit anova on random subset for predicted Secchi
aov_secchi <- aov(secchi~zone*system,data_sampled)
summary(aov_secchi)
secchi_posthoc <- TukeyHSD(aov_secchi)
secchi_posthoc
# buiid data frame with residuals and fitted, to check model diagnostics
resid <- as.numeric(aov_secchi$residuals)
resid_st <- as.numeric(rstudent(aov_secchi)) #
fitted <-as.numeric(aov_secchi$fitted.values)
value <-as.numeric(aov_secchi$model$secchi)
aov_secchi_df<- data.frame(zone=aov_secchi$model$zone,
                           system=aov_secchi$model$system,
                           value,
                           resid,resid_st,fitted)

# this is a ggplot version of the above residual vs fitted, plot(aov1,which=1)
resid_vs_fitted<-aov_secchi_df %>%
  ggplot(aes(x=fitted,y=resid))+
  geom_point()+
  theme_bw()
resid_vs_fitted

# do more diagnostics - acf and pacf of the fitted and residuals

acf_fitted<-acf(aov_secchi_df$fitted)
pacf_fitted<-pacf(aov_secchi_df$fitted)
acf_resid<-acf(aov_secchi_df$resid)
pacf_resid<-pacf(aov_secchi_df$resid)

## doing the same as above but with sat variables including 7/25 waco
# concatenate lat/lon and remove duplicated lat/lon pairs
s2flm_withwaco725$latlon<-paste(s2flm_withwaco725$lat,s2flm_withwaco725$lon)
s2flm_withwaco725 <- s2flm_withwaco725[!duplicated(s2flm_withwaco725$latlon),]

total_rows <- nrow(s2flm_withwaco725)
# use a random subset of data to avoid autocorrelation and speed up processing
# see Loken et al. 2019 https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JG005186
set.seed(1) # makes the subsetting reproducible
sampled <- sample(1:total_rows,size=3055) #2528) #200) # adjust subset size as needed #total is 28728
data_sampled725<-s2flm_withwaco725[sampled,]


# fit anova on random subset for predicted DWL
aov_dwl <- aov(dwl~zone*system,data_sampled)
summary(aov_dwl)
dwl_posthoc <- TukeyHSD(aov_dwl)
dwl_posthoc

# buiid data frame with residuals and fitted, to check model diagnostics
resid <- as.numeric(aov_dwl$residuals)
resid_st <- as.numeric(rstudent(aov_dwl)) #
fitted <-as.numeric(aov_dwl$fitted.values)
value <-as.numeric(aov_dwl$model$dwl)
aov_dwl_df<- data.frame(zone=aov_dwl$model$zone,
                        system=aov_dwl$model$system,
                        value,
                        resid,resid_st,fitted)

# this is a ggplot version of the above residual vs fitted, plot(aov1,which=1)
dwl_resid_vs_fitted<-aov_dwl_df %>%
  ggplot(aes(x=fitted,y=resid))+
  geom_point()+
  theme_bw()
dwl_resid_vs_fitted

# do more diagnostics - acf and pacf of the fitted and residuals

acf_fitted<-acf(aov_dwl_df$fitted)
pacf_fitted<-pacf(aov_dwl_df$fitted)
acf_resid<-acf(aov_dwl_df$resid)
pacf_resid<-pacf(aov_dwl_df$resid)


aov_ndti <- aov(ndti~zone*system,data_sampled)
summary(aov_ndti)
ndti_posthoc <- TukeyHSD(aov_ndti)
ndti_posthoc

# buiid data frame with residuals and fitted, to check model diagnostics
resid <- as.numeric(aov_ndti$residuals)
resid_st <- as.numeric(rstudent(aov_ndti)) #
fitted <-as.numeric(aov_ndti$fitted.values)
value <-as.numeric(aov_ndti$model$ndti )
aov_ndti_df<- data.frame(zone=aov_ndti$model$zone,
                         system=aov_ndti$model$system,
                         value,
                         resid,resid_st,fitted)

# this is a ggplot version of the above residual vs fitted, plot(aov1,which=1)
ndti_resid_vs_fitted<-aov_ndti_df %>%
  ggplot(aes(x=fitted,y=resid))+
  geom_point()+
  theme_bw()
ndti_resid_vs_fitted

# do more diagnostics - acf and pacf of the fitted and residuals

acf_fitted<-acf(aov_ndti_df$fitted)
pacf_fitted<-pacf(aov_ndti_df$fitted)
acf_resid<-acf(aov_ndti_df$resid)
pacf_resid<-pacf(aov_ndti_df$resid)

################################################################################
# testing interaction for whole data
aov_turb_all <- aov(turb~zone*system,s2flame_znpts)
summary(aov_turb_all)
aov_turb_all_posthoc <- TukeyHSD(aov_turb_all)
aov_turb_all_posthoc

aov_secchi_all <- aov(secchi~zone*system,s2flame_znpts)
summary(aov_secchi_all)
aov_secchi_all_posthoc <- TukeyHSD(aov_secchi_all)
aov_secchi_all_posthoc
################################################################################
# summary stats for full system 
# summary stats for them all
systemstats_dwl<- dwl6lake_zone %>% group_by(system) %>% 
  dplyr::summarise(n=n(),mean=mean(dwl),sd=sd(dwl),se=sd/sqrt(n),
                   lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
                   upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
                   min=min(dwl), max=max(dwl),
                   margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))

zonestats_dwl<- dwl6lake_zone %>% group_by(zone) %>% 
  dplyr::summarise(n=n(),mean=mean(dwl),sd=sd(dwl),se=sd/sqrt(n),
                   lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
                   upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
                   min=min(dwl), max=max(dwl),
                   margin_of_error = qt(1 - alpha / 2, df = n - 1) * se,
                   moe= z*se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))

systemstats_ndti<- dwl6lake_zone %>% group_by(system) %>% 
  dplyr::summarise(n=n(),mean=mean(ndti),sd=sd(ndti),se=sd/sqrt(n),
                   lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
                   upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
                   min=min(ndti), max=max(ndti),
                   margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))

zonestats_ndti<- dwl6lake_zone %>% group_by(zone) %>% 
  dplyr::summarise(n=n(),mean=mean(ndti),sd=sd(ndti),se=sd/sqrt(n),
                   lower_ci=mean-qt(1-alpha/2, df = n-1)* se,
                   upper_ci=mean+qt(1-alpha/2, df=n-1)*se,
                   min=min(ndti), max=max(ndti),
                   margin_of_error = qt(1 - alpha / 2, df = n - 1) * se)%>%
  mutate(mean_ci = paste0(round(mean, 2), " ± ", round(margin_of_error, 2)))


aov_dwl_zn<-aov(dwl~zone*system,bb_nona)

summary(aov_dwl_zn)
TukeyHSD(aov_dwl_zn)

# # checking quantiles for each lake
# quantile(ah_df_fui$dwLehmann, c(.05, .1, 0.2, .25, 0.3, 0.5, .75, .95))
# quantile(bn_df_fui$dwLehmann, c(.05, .1, .25, 0.5, .75, .98))
# quantile(bw_df_fui$dwLehmann, c(.05, .1, .25, 0.5, .8, 0.85, 0.9, .95), na.rm=TRUE)
# quantile(iv_df_fui$dwLehmann, c(.05, .1, .25, 0.29,0.5, .78, 0.89, .95), na.rm=TRUE)
# quantile(lw_df_fui$dwLehmann, c(.02, .1, .25, 0.5, .75, 0.9, .975), na.rm=TRUE)
# quantile(rb_df_fui$dwLehmann, c(.05, .1, .25, 0.5, .8, 0.84, 0.9, .95), na.rm=TRUE)



################################################################################
# Just to show how some things were determined

# 1. predicted Secchi
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

# 2. Chromaticity calculation for Dominant Wavelength 
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
R <-df$B4
G <-df$B3
B <-df$B2

# assigning bands to R,G,B
R <-s2flm_l2_scale$B4
G <-s2flm_l2_scale$B3
B <-s2flm_l2_scale$B2

s2flm_l2w_dwl <- s2flm_l2_scale %>% 
  dplyr::mutate(dwl = fui.hue(unlist(R), unlist(G), unlist(B)))

s2flm_l2w_dwl<- s2flm_l2_scale %>% dplyr::mutate(dwl = fui.hue(R, G, B))

# 3. defining zones using bounding boxes, these include entire system
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





#####
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


####################
plot_ndtibn<-OpenStreetMap::autoplot.OpenStreetMap(OpenStreetMap::openproj(map_bon)) +
  geom_point(data = s2flame_znpts %>% filter(system=="bonham"), 
             size=ptsize,
             alpha=ptalpha,
             aes(x = lon, y = lat,color=ndti))+
  ggtitle(label = "Bonham") +
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

